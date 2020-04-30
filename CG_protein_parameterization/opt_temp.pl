#!/usr/bin/perl
#PBS -r n
#PBS -m b
#PBS -m e
#PBS -o $PBS_JOBID.out
#PBS -e $PBS_JOBID.err

chdir($ENV{PBS_O_WORKDIR});

my $n_prod = 5000;
my $skip_step = 2000;
my $nsteps_prod = 2000;

open(LOG, ">opt_temp.log");

my $temp_set = "";
open(TEMP, "<mubrex.cntrl") || die("Error: Cannot find mubrex.cntrl\n\n");
while(my $line = <TEMP>)
{
  if($line =~ /^temps =/)
  {
    print INP $line;
    chomp($line);
    my @str = split(/=/, $line);
    $temp_set = $str[1];
    $temp_set =~ s/^\s+|\s+$//g;
    last;
  }
}
close(TEMP);

my @temp_array = split(/\s+/, $temp_set);
my $n_temp = @temp_array;

my $n_step = 0;

while(@temp_array <= 16)
{
  for(my $i = 1; $i <= $n_temp; $i++)
  {
  	system("rm -rf aa$i");
  }
  system("rm -rf analysis/");
  system("rm -f info.log");
  system("rm -f logs/stats-tswap-1-1.log");

  $n_temp = @temp_array;

  open(TEMP, "<mubrex.cntrl") || die("Error: Cannot find mubrex.cntrl\n\n");
  open(INP, ">mubrex_opt_temp.cntrl") || die("Error: Cannot create mubrex_opt_temp.cntrl\n\n");
  while(my $line = <TEMP>)
  {
    if($line =~ /^nexch_prod =/)
    {
  	  print INP "nexch_prod = $n_prod\n";
    }
    elsif($line =~ /^nwindows =/)
    {
      print INP "nwindows = $n_temp\n";
    }
    elsif($line =~ /^temps =/)
    {
      $temp_set = join(" ", @temp_array);
  	  print INP "temps = $temp_set\n";
    }
    elsif($line =~ /^tpn =/)
    {
      print INP "tpn = $n_temp\n";
    }
    elsif($line =~ /^nsteps_prod =/)
    {
      print INP "nsteps_prod = $nsteps_prod\n";
    }
    elsif($line =~ /^starting_strucs/)
    {
      chomp($line);
      my @str = split(/ = /, $line);
      my $cor = $str[1];
      $cor =~ s/^\s+|\s+$//g;
      for(my $i = 1; $i <= $n_temp; $i++)
      {
        print INP "starting_strucs_t$i = $cor\n";
      }
      last;
    }
    else
    {
      print INP $line;
    }
  }
  close(TEMP);
  close(INP);

  print LOG "TEMP $n_step: $temp_set\n";
  system("parallel_temperature_REX.py -f mubrex_opt_temp.cntrl");

  my @new_temp_array = ();
  
  my $melting_T = check_melting_temp($skip_step+1);
  print LOG "MELT $n_step: $melting_T\n";

  if($melting_T <= $temp_array[2])
  {
    my $max_T = $temp_array[$#temp_array];
    my $min_T = $max_T - 2 * ($max_T - int($melting_T));
    @new_temp_array = regen_temp($min_T, $max_T, $n_temp);
  }
  elsif($melting_T >= $temp_array[$#temp_array-1])
  {
    my $min_T = $temp_array[0];
    my $max_T = $min_T + 2 * (int($melting_T) - $min_T);
    @new_temp_array = regen_temp($min_T, $max_T, $n_temp);
  }
  else
  {
    my @exchange_rate = check_exc_rate($n_temp, $n_prod);
    if(!@exchange_rate)
    {
      die("Error: Cannot find T-EXCHNG $n_prod in logs/stats-tswap-1-1.log\n\n");
    }
    else
    {
      my $str = join(" ", @exchange_rate);
      print LOG "EXCR $n_step: $str\n";
    }

    my @insert_points = ();
    my @delete_points = ();

    my @ex_int = ();
    foreach my $ex (@exchange_rate)
    {
      push@ex_int, int($ex * 10);
    }

    my @freq = ();
    for(my $i = 0; $i < 10; $i++)
    {
      my $n = 0;
      foreach my $ex (@ex_int)
      {
        if($ex eq $i)
        {
          $n++;
        }
      }
      push@freq, $n;
    }

    my $max_freq = 0;
    my @max_i = ();
    for(my $i = 0; $i < 10; $i++)
    {
      if($freq[$i] > $max_freq)
      {
        $max_freq = $freq[$i];
        @max_i = ($i);
      }
      elsif($freq[$i] eq $max_freq)
      {
        push@max_i, $i;
      }
    }

    my $mean = 0;
    my $n = 0;
    my $i = 0;
    foreach my $ex (@ex_int)
    {
      if($ex ~~ @max_i)
      {
        $mean += $exchange_rate[$i];
        $n++;
      }
      $i++;
    }
    $mean /= $n;

    for(my $i = 1; $i <= @exchange_rate; $i++)
    {
      if(($exchange_rate[$i-1] - $mean >= 0.2)&&($exchange_rate[$i-1] >= 0.8)&&($mean >= 0.6))
      {
        if($i eq 1)
        {
          push@delete_points, 2;
        }
        elsif($i eq @exchange_rate)
        {
          push@delete_points, $i-1;
        }
        else
        {
          push@delete_points, $i;
        }
      }
      elsif(($mean - $exchange_rate[$i-1] >= 0.2)||($exchange_rate[$i-1] <= 0.4))
      {
        if($i eq @exchange_rate)
        {
          push@insert_points, $i-1;
        }
        elsif($i eq 1)
        {
          push@insert_points, $i;
        }
        else
        {
          push@insert_points, $i-1;
          push@insert_points, $i;
        }
      }
    }
    @insert_points = do { my %seen; grep {!$seen{$_}++} @insert_points }; # Get unique elements in array
    @delete_points = do { my %seen; grep {!$seen{$_}++} @delete_points };

    my $insert_str = join(" ", @insert_points);
    my $delete_str = join(" ", @delete_points);
    print LOG "INTW $n_step: $insert_str\n";
    print LOG "DELW $n_step: $delete_str\n";

    if(@insert_points eq 0 && @delete_points eq 0)
    {
      print LOG "Optimization completed.\n";
      last;
    }

    for(my $i = 1; $i <= $n_temp; $i++)
    {
      push@new_temp_array, $temp_array[$i-1];
      if($i ~~ @insert_points)
      {
        my $new_t = int(($temp_array[$i-1]+$temp_array[$i])/2);
        push@new_temp_array, $new_t;
      }

      if($i ~~ @delete_points)
      {
        pop @new_temp_array;
      }
    }
  }

  #my @results = check_sampling("info.log", 0, 100000);
  #my @T_access = @{$results[0]};
  #my @insert_points = @{$results[1]};

  #my $T_access_str = join(" ", @T_access);
  #my $insert_str = join(" ", @insert_points);
  #print LOG "TACT $n_step: $T_access_str\n";
  #print LOG "INTW $n_step: $insert_str\n";

  #if(@insert_points eq 0)
  #{
  #  last;
  #}

  #for(my $i = 1; $i <= $n_temp; $i++)
  #{
  #  push@new_temp_array, $temp_array[$i-1];
  #  if($i ~~ @insert_points)
  #  {
  #    my $new_t = int(($temp_array[$i-1]+$temp_array[$i])/2);
  #    push@new_temp_array, $new_t;
  #  }
  #}

  @temp_array = @new_temp_array;
  
  $n_step++;
}

close(LOG);

######################################################################
sub check_sampling
{
  my $log_file = $_[0];
  my $skip_step = $_[1];
  my $end_step = $_[2];

  my $step = 0;
  my $n_windows;
  my @counter = ();
  my @tag = ();
  my @max_window = ();
  my @min_window = ();
  open(DAT, "<$log_file") || die("Error: Cannot find $log_file\n\n");
  while(my $line = <DAT>)
  {
    if($line =~ /^nwindows =/)
    {
      chomp($line);
      my @str = split(/=/, $line);
      $n_windows = $str[1];
      $n_windows =~ s/^\s+|\s+$//g;
      for(my $i = 0; $i < $n_windows; $i++)
      {
        push@counter, 0;
        push@tag, -1;
        push@max_window, 0;
        push@min_window, 1000;
      }
    }
    elsif($line =~ /^PROD [0-9]+:/)
    {
      $step++;
      if($step > $skip_step && $step <= $end_step)
      {
        my @str = split(/: /, $line);
        my $window = $str[1];
        $window =~ s/\|//g;
        @str = split(/\s+/, $window);
        my @windows = @str[0..$n_windows-1];

        for(my $i = 1; $i <= $n_windows; $i++)
        {
          for(my $j = 1; $j <= $n_windows; $j++)
          {
            if($windows[$j] eq $i)
            {
              if($j > $max_window[$i-1])
              {
                $max_window[$i-1] = $j;
              }
              if($j < $min_window[$i-1])
              {
                $min_window[$i-1] = $j;
              }
              last;
            }
          }

          if($line =~ /: \Q$i\E\|/)
          {
            if($tag[$i-1] eq -1 && $counter[$i-1] eq 0)
            {
              $tag[$i-1] = 0;
            }
            if($tag[$i-1] eq 0 && $counter[$i-1] % 2 ne 0)
            {
              $counter[$i-1]++;
            }
            elsif($tag[$i-1] eq 1 && $counter[$i-1] % 2 eq 0)
            {
              $counter[$i-1]++;
            }
          }
          elsif($line =~ / \Q$i\E\| /)
          {
            if($tag[$i-1] eq -1 && $counter[$i-1] eq 0)
            {
              $tag[$i-1] = 1;
            }
            if($tag[$i-1] eq 0 && $counter[$i-1] % 2 eq 0)
            {
              $counter[$i-1]++;
            }
            elsif($tag[$i-1] eq 1 && $counter[$i-1] % 2 ne 0)
            {
              $counter[$i-1]++;
            }
          }
        }
      }
    }
  }
  close(DAT);
  
  my @insert_points = ();
  for(my $i = 0; $i < $n_windows; $i++)
  {
    $counter[$i] = $counter[$i] / 2;
    if($counter[$i] eq 0)
    {
      my $max = $max_window[$i];
      my $min = $min_window[$i];
      if($min ne 1)
      {
        if(!($min-1 ~~ @insert_points))
        {
          push@insert_points, $min-1;
        }
      }
      if($max ne $n_windows)
      {
        if(!($max ~~ @insert_points))
        {
          push@insert_points, $max;
        }
      }
    }
  }
  
  @insert_points = sort {$a <=> $b} @insert_points;
  return (\@counter, \@insert_points);
}

###############################################
sub check_exc_rate
{
  my $n_temp = $_[0];
  my $n_prod = $_[1];
  
  my @exchange_rate = ();
  open(DAT, "<logs/stats-tswap-1-1.log") || die("Error: Cannot find logs/stats-tswap-1-1.log");
  while(my $line = <DAT>)
  {
    if($line =~ /^T-EXCHNG $n_prod:/)
    {
      chomp($line);
      $line =~ s/^\s+|\s+$//g;
      my @str = split(/\s+/, $line);
      my $start_index = $#str-$n_temp+1;
      @exchange_rate = @str[$start_index..$#str];
    }
  }
  close(DAT);

  return @exchange_rate;
}

###############################################
sub check_melting_temp
{
  my $skip_step = $_[0];
  my $melting_T;

  my $screen = `analysis_folding_stability.pl -i info.log -m ca -t 0 -s $skip_step -e "pywham.py" 2>&1`;
  if($screen =~ /error/i)
  {
    die($screen);
  }
  else
  {
    my @str = split(/\n/, $screen);
    foreach my $line (@str)
    {
      if($line =~ /Done. Melting temperature is /)
      {
        my @s = split(/ is /, $line);
        $ss = $s[1];
        @s = split(/ K/, $ss);
        $melting_T = $s[0];
      }
    }
  }

  return $melting_T;
}

##################################################
sub regen_temp
{
  my $min_T = $_[0];
  my $max_T = $_[1];
  my $n_windows = $_[2];

  my $c = 1/$n_windows;
  my $b = ($min_T * exp($c * ($n_windows-1)) - $max_T) / (exp($c * ($n_windows-1)) - 1);
  my $a = ($min_T - $b) / exp($c);
  
  my @temp_array = ();
  for(my $i = 1; $i <= $n_windows; $i++)
  {
    if($i eq 1)
    {
      push@temp_array, $min_T;
    }
    elsif($i eq $n_windows)
    {
      push@temp_array, $max_T;
    }
    else
    {
      my $t = int($a * exp($c * $i) + $b);
      push@temp_array, $t;
    }
  }

  return @temp_array;
}
