#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;

my ($help, $input_cor, $sec_def, $dom_def, $traj, $start, $end, $if_meaningful, $mask, $output_dir, $restart);
GetOptions(
 'help|h!' => \$help,
 'input|i=s' => \$input_cor,
 'secs|s=s' => \$sec_def,
 'domain|d=s' => \$dom_def,
 'traj|t=s' => \$traj,
 'begin|b=s' => \$start,
 'end|e=s' => \$end,
 'meaningful|m=s' => \$if_meaningful,
 'mask|k=s' => \$mask,
 'outdir|o=s' => \$output_dir,
 'restart|r=s' => \$restart,
);

my $usage = "
  Usage: perl calc_native_contact_fraction.pl 
              --input | -i <INPUT.COR> for identify native contacts
              --domain | -d <DOMAIN.DAT> for domain defination
              [--secs | -s] <SECONDARY STRUCTURE> for secondary 
                          structure defination. If no file specified,
                          it will calculate all native contacts regardless
                          of the secondary structure.
              --traj | -t <TRAJ.DCD> for simulation trajectory
              [--begin | -b] <START STEP> to split trajectory. 
                                  Default is 1
              [--end | -e] <END STEP> to split trajectory. 
                                Default is the end frame.
              [--meaningful | -m] <1 or 0> 1 to only calculate contacts
                                  in the interface has native contacts;
                                  0 to calculate all possible interfaces.
                                  The interface has no native contacts will
                                  always has Q = 1.
                                  Default is 0. 
              [--mask | -k] <MASK> for a subset of residues to calculate their
                                   native contacts with all other residues. 
                                   Default is 'all'.
                                   Example MASK: 
                                   '1-45' for selection from resid 1 to 45;
                                   '1,3,5-9' for selection of resid 1, 3 
                                   and from 5 to 9.
              [--outdir | -o] <DIRECTORY> for the outputs. Default is the directory
                                          of your trajectories
              [--restart | -r] <0 or 1> 0: Do not restart calculation; 1: Restart 
                                        calculation. Default is 0.
              [--help | -h]\n
  Example for DOMAIN.DAT:
  1:96 302:350 a #Domain 1 is from resid 1 to 96 and 302 to 350, and in
                 #alpha-helix class
  97:155 b #Domain 2 is from resid 97 to 155 and in beta-sheet class
  156:301 c #Domain 3 is from resid 156 to 301 and in alpha-beta class\n\n";

if($help)
{
  die($usage);
}
elsif(!(defined($input_cor) && defined($dom_def) && defined($traj)))
{
  die($usage);
}

if(!defined($sec_def))
{
  $sec_def = "";
}
if(!defined($start))
{
  $start = 1;
}
if(!defined($end))
{
  $end = 9999999999999999999;
}
if(!defined($if_meaningful))
{
  $if_meaningful = 0;
}
if(!defined($mask))
{
  $mask = "all";
}

my @str = ();
if($traj =~ /\.cor$/)
{
  @str = split(/\.cor/, $traj);
}
else
{
  @str = split(/\.dcd/, $traj);
}
@str = split(/\//, $str[0]);
my $name = $str[$#str];
my $dir = join("/", @str[0...$#str-1]);
if(@str eq 1)
{
  $dir = './';
}

if(!defined($output_dir))
{
  $output_dir = $dir;
}

if(!defined($restart))
{
  $restart = 0;
}
elsif($restart != 0)
{
  $restart = 1;
}

my $cutoff = 8;

my @domain = parse_domain($dom_def);
my @sel_idx = parse_mask($mask);
if($sel_idx[0] eq "all")
{
  print "-> Select all residues\n";
}
else
{
  print "-> Select ".@sel_idx." residues\n";
}

my @native_cor = parse_cor($input_cor);
my $native_natom = @native_cor;

my @results = calc_contact_dist_map(\@native_cor, $cutoff);
my @native_contact_map = @{$results[0]};
my @native_distance_map = @{$results[1]};

my @sec_struc = parse_secondary_structure($sec_def);

my @native_contact_number = calc_contact_number(\@native_contact_map, \@domain, \@sec_struc, \@sel_idx);
my @meaningful_domain_idx = ();
print "-> native contacts found in $input_cor:\n";
for(my $i = 1; $i <= @domain; $i++)
{
  my $dom = $domain[$i-1];
  if($native_contact_number[$i-1] ne 0)
  {
    push@meaningful_domain_idx, $i;
  }
  
  if($dom->{"class"} eq "i")
  {
    print "   Interface ".$dom->{"range"}->[0]."|".$dom->{"range"}->[1].": ".$native_contact_number[$i-1]."\n";
  }
  else
  {
    print "   Domain ".$i.": ".$native_contact_number[$i-1]."\n";
  }
}

if(!$if_meaningful)
{
  @meaningful_domain_idx = (1..@domain);
}

if($traj =~ /\.cor$/)
{
  open(DAT, ">${output_dir}/qbb_${name}.dat")||die("Error: cannot create ${output_dir}/qbb_${name}.dat\n\n");
  foreach my $i (@meaningful_domain_idx)
  {
    my $dom = $domain[$i-1];
    if($dom->{"class"} eq "i")
    {
      printf DAT ("%10s ", $dom->{"range"}->[0]."|".$dom->{"range"}->[1]);
    }
    else
    {
      printf DAT ("%10s ", "D_$i");
    }
  }
  printf DAT ("%10s\n", "total");
  my @traj_cor = parse_cor($traj);
  my $natom = @traj_cor;
  if($natom ne $native_natom)
  {
    print("Warning: atom numbers mismatch in $input_cor ($native_natom) and $traj ($natom)\n")
  }
  my @traj_contact_number = calc_traj_contact_number(\@traj_cor, \@domain, \@sec_struc, \@native_contact_map, \@native_distance_map, \@sel_idx);
  my $tot_tcn = 0;
  my $tot_ncn = 0;
  foreach my $j (@meaningful_domain_idx)
  {
    my $fraction = 1;
    if($native_contact_number[$j-1] != 0)
    {
      $fraction = $traj_contact_number[$j-1] / $native_contact_number[$j-1];
    }
    else
    {
      $fraction = -1;
    }
    $fraction = sprintf("%.4f", $fraction);
    printf DAT ("%10s ", $fraction);
    $tot_tcn += $traj_contact_number[$j-1];
    $tot_ncn += $native_contact_number[$j-1];
  }
  my $fraction = 1;
  if($tot_ncn != 0)
  {
    $fraction = $tot_tcn / $tot_ncn;
  }
  else
  {
    $fraction = -1;
  }
  $fraction = sprintf("%.4f", $fraction);
  printf DAT ("%10s\n", $fraction);
  exit;
}

open(my $dcdfile, "<$traj") || die("Error: Cannot find $traj\n\n");
binmode($dcdfile);
my ($natom, $tstep, $nframe, $first, $delta, $crystal, $fixed) = read_DCD_head($dcdfile);

if($natom ne $native_natom)
{
  print("Warning: atom numbers mismatch in $input_cor ($native_natom) and $traj ($natom)\n")
}

if($restart && -s "${output_dir}/qbb_${name}.dat")
{
  open(DAT, "<${output_dir}/qbb_${name}.dat")||die("Error: cannot find ${output_dir}/qbb_${name}.dat\n\n");
  open(DAT_COPY, ">${output_dir}/qbb_${name}_copy.dat")||die("Error: cannot create ${output_dir}/qbb_${name}_copy.dat\n\n");
  
  foreach my $i (@meaningful_domain_idx)
  {
    my $dom = $domain[$i-1];
    if($dom->{"class"} eq "i")
    {
      printf DAT_COPY ("%10s ", $dom->{"range"}->[0]."|".$dom->{"range"}->[1]);
    }
    else
    {
      printf DAT_COPY ("%10s ", "D_$i");
    }
  }
  printf DAT_COPY ("%10s\n", "total");
  
  my $num_lines = 0;
  while(my $line = <DAT>)
  {
    if($line =~ /^\s*\-?[0-9].+\n$/)
    {
      $num_lines++;
      print DAT_COPY "$line";
    }
  }
  close(DAT);
  close(DATA_COPY);
  `rm -f ${output_dir}/qbb_${name}.dat`;
  `mv ${output_dir}/qbb_${name}_copy.dat ${output_dir}/qbb_${name}.dat`;
  
  $start += $num_lines;
  open(DAT, ">>${output_dir}/qbb_${name}.dat")||die("Error: cannot find ${output_dir}/qbb_${name}.dat\n\n");
}
else
{
  open(DAT, ">${output_dir}/qbb_${name}.dat")||die("Error: cannot create ${output_dir}/qbb_${name}.dat\n\n");
  foreach my $i (@meaningful_domain_idx)
  {
    my $dom = $domain[$i-1];
    if($dom->{"class"} eq "i")
    {
      printf DAT ("%10s ", $dom->{"range"}->[0]."|".$dom->{"range"}->[1]);
    }
    else
    {
      printf DAT ("%10s ", "D_$i");
    }
  }
  printf DAT ("%10s\n", "total");
}

my $deltat = $delta * $tstep;
my $firstframe = $first / $delta;
for(my $i = 1; $i <= $nframe; $i++)
{
  my @traj_cor = read_DCD_frame($dcdfile, $crystal);
  if($i >= $start && $i <= $end)
  {
    my @traj_contact_number = calc_traj_contact_number(\@traj_cor, \@domain, \@sec_struc, \@native_contact_map, \@native_distance_map, \@sel_idx);
    my $tot_tcn = 0;
    my $tot_ncn = 0;
    foreach my $j (@meaningful_domain_idx)
    {
      my $fraction = 1;
      if($native_contact_number[$j-1] != 0)
      {
        $fraction = $traj_contact_number[$j-1] / $native_contact_number[$j-1];
      }
      else
      {
        $fraction = -1;
      }
      $fraction = sprintf("%.4f", $fraction);
      printf DAT ("%10s ", $fraction);
      $tot_tcn += $traj_contact_number[$j-1];
      $tot_ncn += $native_contact_number[$j-1];
    }
    my $fraction = 1;
    if($tot_ncn != 0)
    {
      $fraction = $tot_tcn / $tot_ncn;
    }
    else
    {
      $fraction = -1;
    }
    $fraction = sprintf("%.4f", $fraction);
    printf DAT ("%10s\n", $fraction);
  }
}
close(DAT);

######################################################
sub parse_domain
{
  my $dom_def = $_[0];

  my @domain = ();

  my %class = ("a" => "Alpha-helix",
               "b" => "Beta-sheet",
               "c" => "Alpha-Beta");

  print "-> Domain defination:\n";

  open(DAT, "<$dom_def") || die("Error: Cannot find $dom_def\n\n");
  my $n = 1;
  my $start_min = 10000;
  while(my $line = <DAT>)
  {
    chomp($line);
    $line =~ s/^\s+|\s+$//g;
    if($line !~ /^#/ && $line ne '')
    {
      my @str = split(/\s*#/, $line);
      my @data = split(/\s+/, $str[0]);
      my @res_range = @data[0..$#data-1];
      my $sec_class = $data[$#data];
      print "   Domain $n: ";
      if($sec_class ne "a" && $sec_class ne "b" && $sec_class ne "c")
      {
        die("Error: Forget to specify class as a, b or c?\n\n");
      }
      my %dom = ("range" => [],
                 "class" => $sec_class);
      foreach my $d (@res_range)
      {
        my @data = split(/\s*:\s*/, $d);
        my $start = $data[0];
        my $end = $data[1];
        print "$start ~ $end; ";
        if($start < $start_min)
        {
          $start_min = $start;
        }
        push@{$dom{"range"}}, [$start, $end];
      }
      print "Class: ".$class{$sec_class}."\n";
      push@domain, \%dom;
      $n++;
    }
  }
  close(DAT);

  my $offset = $start_min - 1;
  foreach my $d (@domain)
  {
    foreach my $r (@{$d->{"range"}})
    {
      $r->[0] -= $offset;
      $r->[1] -= $offset;
    }
  }

  for(my $i = 1; $i < $n-1; $i++)
  {
    for(my $j = $i+1; $j < $n; $j++)
    {
      my %hash = ("range" => [$i, $j],
                  "class" => "i",);
      push@domain, \%hash;
    }
  }
  return @domain;
}

######################################################
sub parse_mask
{
  my $mask = $_[0];
  my @sel_idx = ();
  $mask =~ s/^\s+|\s+$//g;
  if($mask eq "all")
  {
    @sel_idx = ("all");
  }
  else
  {
    my @str = split(/,\s*/, $mask);
    foreach my $s (@str)
    {
      if($s =~ /\-/)
      {
        my @str_1 = split(/\s*\-\s*/, $s);
        for(my $i = $str_1[0]; $i <= $str_1[1]; $i++)
        {
          push@sel_idx, $i;
        }
      }
      else
      {
        push@sel_idx, $s;
      }
    }
  }
  return @sel_idx;
}

######################################################
sub parse_cor
{
  my $cor = $_[0];
  
  my @native_cor = ();
  my $num_atom = 0;
  my $if_ext = 0;
  open(COR, "<$cor") || die("Error: cannot find $cor\n\n");
  while(my $line = <COR>)
  {
  	chomp($line);
  	if($line !~ /^\*/)
  	{
  	  if($line =~ /^\s*[0-9]+$/)
  	  {
  	  	$line =~ s/^\s+|\s+$//g;
  	  	$num_atom = $line;
  	  }
  	  elsif($line =~ /^\s*[0-9]+\s+EXT$/)
  	  {
  	    $line =~ s/^\s+|\s+$//g;
  	    $if_ext = 1;
  	    my @str = split(/\s+/, $line);
  	    $num_atom = $str[0];
  	  }
  	  else
  	  {
  	    my ($x, $y, $z) = (0, 0, 0);
  	    if($if_ext eq 0)
  	    {
  	      ($x, $y, $z) = unpack("x5 x5 x6 x1 x3 a10 a10 a10 x20", $line);
  	    }
  	  	else
  	  	{
  	  	  ($x, $y, $z) = unpack("x10 x10 x10 x10 a20 a20 a20 x10 x10 x20", $line);
  	  	}
  	  	$x =~ s/^\s+|\s+$//g;
  	  	$y =~ s/^\s+|\s+$//g;
  	  	$z =~ s/^\s+|\s+$//g;
  	  	push@native_cor, [$x, $y, $z];
  	  }
  	}
  }
  close(COR);

  return @native_cor;
}

######################################################
sub parse_secondary_structure
{
  my $sec_def = $_[0];
  my @sec_struc = ();
  print "-> Secondary structure defination:\n";
  print "   ";
  if($sec_def eq "")
  {
    print "None\n";
    return ([1, 99999999999999]);
  }
  open(DAT, "<$sec_def") || die("Error: cannot find $sec_def\n\n");
  while(my $line = <DAT>)
  {
  	chomp($line);
  	my @str = split(/\s+/, $line);
  	push@sec_struc, [$str[1], $str[2]];
    print $str[1]." ~ ".$str[2]."; ";
  }
  close(DAT);
  print "\n";
  return @sec_struc;
}

#########################################################
sub calc_contact_dist_map
{
  my @coor = @{$_[0]};
  my $cutoff = $_[1];

  my @c_map = ();
  my @d_map = ();
  for(my $i = 0; $i < @coor; $i++)
  {
    $c_map[$i] = [];
    $d_map[$i] = [];
  }

  for(my $i = 0; $i < @coor-4; $i++)
  {
    for(my $j = $i + 4; $j < @coor; $j++)
    {
      my $distance = ($coor[$i]->[0] - $coor[$j]->[0]) ** 2 + ($coor[$i]->[1] - $coor[$j]->[1]) ** 2 + ($coor[$i]->[2] - $coor[$j]->[2]) ** 2;
      $distance = $distance ** 0.5;

      if($distance <= $cutoff)
      {
        $c_map[$i]->[$j] = 1;
        $c_map[$j]->[$i] = 1;
        $d_map[$i]->[$j] = $distance;
        $d_map[$j]->[$i] = $distance;
      }
      else
      {
        $c_map[$i]->[$j] = 0;
        $c_map[$j]->[$i] = 0;
        $d_map[$i]->[$j] = 0;
        $d_map[$j]->[$i] = 0;
      }
    }
  }

  my @results = (\@c_map, \@d_map);
  return @results;
}

#########################################################
sub calc_contact_number
{
  my @map = @{$_[0]};
  my @domain = @{$_[1]};
  my @sec_struc = @{$_[2]};
  my @sel_idx = @{$_[3]};
  
  my @result = ();
  foreach my $dom (@domain)
  {
    my $contact_num = 0;
    my @range = @{$dom->{"range"}};
    if($dom->{"class"} eq "i")
    {
      my @range_1 = @{$domain[$range[0]-1]->{"range"}};
      my @range_2 = @{$domain[$range[1]-1]->{"range"}};
      foreach my $rd_1 (@range_1)
      {
        for(my $i = $rd_1->[0]; $i <= $rd_1->[1]; $i++)
        {
          my $tag_i = 0;
          foreach my $rs (@sec_struc)
          {
            if($i >= $rs->[0] && $i <= $rs->[1])
            {
              $tag_i = 1;
              last;
            }
          }

          if($tag_i eq 0)
          {
            next;
          }

          foreach my $rd_2 (@range_2)
          {
            for(my $j = $rd_2->[0]; $j <= $rd_2->[1]; $j++)
            {
              my $tag_j = 0;
              foreach my $rs (@sec_struc)
              {
                if($j >= $rs->[0] && $j <= $rs->[1])
                {
                  $tag_j = 1;
                  last;
                }
              }

              if($tag_j eq 0)
              {
                next;
              }

              if($map[$i-1]->[$j-1] == 1)
              {
                if($sel_idx[0] eq "all")
                {
                  $contact_num++;
                }
                elsif(($i ~~ @sel_idx) || ($j ~~ @sel_idx))
                {
                  $contact_num++;
                }
              }
            }
          }
        }
      }
    }
    else
    {
      foreach my $rd (@range)
      {
        for(my $i = $rd->[0]; $i <= $rd->[1]-1; $i++)
        {
          my $tag_i = 0;
          foreach my $rs (@sec_struc)
          {
            if($i >= $rs->[0] && $i <= $rs->[1])
            {
              $tag_i = 1;
              last;
            }
          }

          if($tag_i eq 0)
          {
            next;
          }

          for(my $j = $i+1; $j <= $rd->[1]; $j++)
          {
            my $tag_j = 0;
            foreach my $rs (@sec_struc)
            {
              if($j >= $rs->[0] && $j <= $rs->[1])
              {
                $tag_j = 1;
                last;
              }
            }

            if($tag_j eq 0)
            {
              next;
            }

            if($map[$i-1]->[$j-1] == 1)
            {
              if($sel_idx[0] eq "all")
              {
                $contact_num++;
              }
              elsif(($i ~~ @sel_idx) || ($j ~~ @sel_idx))
              {
                $contact_num++;
              }
            }
          }
        }
      }
    }
    push@result, $contact_num;
  }
  return @result;
}

#########################################################
sub calc_traj_contact_number
{
  my @coor = @{$_[0]};
  my @domain = @{$_[1]};
  my @sec_struc = @{$_[2]};
  my @c_map = @{$_[3]};
  my @d_map = @{$_[4]};
  my @sel_idx = @{$_[5]};
  my $sdist = 1.2;
  
  my @result = ();
  foreach my $dom (@domain)
  {
    my $contact_num = 0;
    my @range = @{$dom->{"range"}};
    if($dom->{"class"} eq "i")
    {
      my @range_1 = @{$domain[$range[0]-1]->{"range"}};
      my @range_2 = @{$domain[$range[1]-1]->{"range"}};
      foreach my $rd_1 (@range_1)
      {
        for(my $i = $rd_1->[0]; $i <= $rd_1->[1]; $i++)
        {
          my $tag_i = 0;
          foreach my $rs (@sec_struc)
          {
            if($i >= $rs->[0] && $i <= $rs->[1])
            {
              $tag_i = 1;
              last;
            }
          }

          if($tag_i eq 0)
          {
            next;
          }

          foreach my $rd_2 (@range_2)
          {
            for(my $j = $rd_2->[0]; $j <= $rd_2->[1]; $j++)
            {
              my $tag_j = 0;
              foreach my $rs (@sec_struc)
              {
                if($j >= $rs->[0] && $j <= $rs->[1])
                {
                  $tag_j = 1;
                  last;
                }
              }

              if($tag_j eq 0)
              {
                next;
              }

              if($c_map[$i-1]->[$j-1] == 1)
              {
                my $distance = 0;
                if($coor[$i-1] ne "" && $coor[$j-1] ne "")
                {
                  $distance = ($coor[$i-1]->[0] - $coor[$j-1]->[0]) ** 2 + ($coor[$i-1]->[1] - $coor[$j-1]->[1]) ** 2 + ($coor[$i-1]->[2] - $coor[$j-1]->[2]) ** 2;
                  $distance = $distance ** 0.5;
                }
                if($distance <= $sdist * $d_map[$i-1]->[$j-1] && $distance > 0)
                {
                  if($sel_idx[0] eq "all")
                  {
                    $contact_num++;
                  }
                  elsif(($i ~~ @sel_idx) || ($j ~~ @sel_idx))
                  {
                    $contact_num++;
                  }
                }
              }
            }
          }
        }
      }
    }
    else
    {
      foreach my $rd (@range)
      {
        for(my $i = $rd->[0]; $i <= $rd->[1]-1; $i++)
        {
          my $tag_i = 0;
          foreach my $rs (@sec_struc)
          {
            if($i >= $rs->[0] && $i <= $rs->[1])
            {
              $tag_i = 1;
              last;
            }
          }

          if($tag_i eq 0)
          {
            next;
          }

          for(my $j = $i+1; $j <= $rd->[1]; $j++)
          {
            my $tag_j = 0;
            foreach my $rs (@sec_struc)
            {
              if($j >= $rs->[0] && $j <= $rs->[1])
              {
                $tag_j = 1;
                last;
              }
            }

            if($tag_j eq 0)
            {
              next;
            }

            if($c_map[$i-1]->[$j-1] == 1)
            {
              my $distance = 0;
              if($coor[$i-1] ne "" && $coor[$j-1] ne "")
              {
                $distance = ($coor[$i-1]->[0] - $coor[$j-1]->[0]) ** 2 + ($coor[$i-1]->[1] - $coor[$j-1]->[1]) ** 2 + ($coor[$i-1]->[2] - $coor[$j-1]->[2]) ** 2;
                $distance = $distance ** 0.5;
              }
              if($distance <= $sdist * $d_map[$i-1]->[$j-1] && $distance > 0)
              {
                if($sel_idx[0] eq "all")
                {
                  $contact_num++;
                }
                elsif(($i ~~ @sel_idx) || ($j ~~ @sel_idx))
                {
                  $contact_num++;
                }
              }
            }
          }
        }
      }
    }
    push@result, $contact_num;
  }
  return @result;
}

#########################################################
sub readFortran 
{
  my $handle=shift;

  my $dat;
  my $tdat;

  read($handle,$tdat,4) || die "cannot read data";
  my $len=unpack("L",$tdat);
  read($handle,$dat,$len) || die "cannot read data";
  read($handle,$tdat,4) || die "cannot read data";

#  printf STDERR "Fread %d bytes\n",$len;

  return ($dat,$len);
}

################################################################
sub read_DCD
{
  my $traj = $_[0];
  my $boxsize = $_[1];
  my $from = 1;
  my $to = 999999999999999;
  my $step = 1;

  open(my $dcdfile, "<$traj") || die("Error: Cannot find $traj\n\n");
  binmode($dcdfile);
  my $buffer;
  my $len;
  ($buffer,$len)=readFortran($dcdfile);
  my ($tag,@icontrol)=unpack("A4L*",$buffer);
  
  ($buffer,$len)=readFortran($dcdfile);
  ($buffer,$len)=readFortran($dcdfile);
  my $natom=unpack("L",$buffer);
  
  my $tstep=unpack("f",pack("L",$icontrol[9]))*4.88882129E-02;
  my $nfiles=$icontrol[0];
  my $first=$icontrol[1];
  my $delta=$icontrol[2];
  my $deltat=$icontrol[2]*$tstep;
  my $crystal=$icontrol[10];
  my $fixed=$icontrol[8];

  my $firstframe=$first/$delta;

  my ($xbuf,$ybuf,$zbuf);
  
  my $itot = 0;
  for (my $i=1; $itot<=$to && $i<=$nfiles; $i++) 
  {
    $itot++;
    if ($crystal) 
    {
      my ($tbuf,$tlen)=readFortran($dcdfile);
      if ($boxsize)
      {
        my @cdat=unpack("d*",$tbuf);
        printf "%f %f %f\n",$cdat[0],$cdat[2],$cdat[5];
      }
    }

    ($xbuf,$len)=readFortran($dcdfile); # printf STDERR "%d ",$len;
    ($ybuf,$len)=readFortran($dcdfile); # printf STDERR "%d ",$len;
    ($zbuf,$len)=readFortran($dcdfile); # printf STDERR "%d \n",$len;
      
    if ($itot>=$from && $itot<=$to && ($itot%$step)==0 && !$boxsize)
    {
      my @xcoor=unpack("f*",$xbuf);
      my @ycoor=unpack("f*",$ybuf);
      my @zcoor=unpack("f*",$zbuf);
      for(my $j = 0; $j < @xcoor; $j++)
      {
        print "$xcoor[$j] $ycoor[$j] $zcoor[$j]\n";
      }
      print "\n";
    }
  }
  close($dcdfile);
}

################################################################
sub read_DCD_head
{
  my $dcdfile = $_[0];
  my $buffer;
  my $len;
  ($buffer,$len)=readFortran($dcdfile);
  my ($tag,@icontrol)=unpack("A4L*",$buffer);
  
  ($buffer,$len)=readFortran($dcdfile);
  ($buffer,$len)=readFortran($dcdfile);
  my $natom=unpack("L",$buffer);
  
  my $tstep=unpack("f",pack("L",$icontrol[9]))*4.88882129E-02;
  my $nframe=$icontrol[0];
  my $first=$icontrol[1];
  my $delta=$icontrol[2];
  my $crystal=$icontrol[10];
  my $fixed=$icontrol[8];

  my @results = ($natom, $tstep, $nframe, $first, $delta, $crystal, $fixed);
  return @results;
}

################################################################
sub read_DCD_frame
{
  my $dcdfile = $_[0];
  my $crystal = $_[1];

  my @coor = ();

  my ($xbuf,$ybuf,$zbuf);
  
  if ($crystal) 
  {
    my ($tbuf,$tlen)=readFortran($dcdfile);
  }

  ($xbuf,$len)=readFortran($dcdfile);
  ($ybuf,$len)=readFortran($dcdfile);
  ($zbuf,$len)=readFortran($dcdfile);
      
  my @xcoor=unpack("f*",$xbuf);
  my @ycoor=unpack("f*",$ybuf);
  my @zcoor=unpack("f*",$zbuf);
  for(my $j = 0; $j < @xcoor; $j++)
  {
    push@coor, [$xcoor[$j], $ycoor[$j], $zcoor[$j]];
  }
  
  return @coor;
}
