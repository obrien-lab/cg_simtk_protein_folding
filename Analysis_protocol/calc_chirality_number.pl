#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;

my ($help, $input_cor, $sec_def, $dom_def, $traj, $start, $end, $output_dir);
GetOptions(
 'help|h!' => \$help,
 'input|i=s' => \$input_cor,
 'secs|s=s' => \$sec_def,
 'domain|d=s' => \$dom_def,
 'traj|t=s' => \$traj,
 'begin|b=s' => \$start,
 'end|e=s' => \$end,
 'outdir|o=s' => \$output_dir,
);

my $usage = "
  Usage: perl calc_chirality_number.pl 
              --input | -i <INPUT.COR> for identify native contacts
              --domain | -d <DOMAIN.DAT> for domain defination
              [--secs | -s] <SECONDARY STRUCTURE> for secondary 
                          structure defination. If no file specified,
                          it will calculate all residues regardless
                          of the secondary structure.
              --traj | -t <TRAJ.DCD> for simulation trajectory
              [--begin | -b] <START STEP> to split trajectory. 
                                  Default is 1
              [--end | -e] <END STEP> to split trajectory. 
                                Default is the end frame.
              [--outdir | -o] <DIRECTORY> for the outputs. Default is the directory
                                          of your trajectories
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

my $cutoff = 0;

my @domain = parse_domain($dom_def);

my @native_cor = parse_cor($input_cor);
my $native_natom = @native_cor;

my @sec_struc = parse_secondary_structure($sec_def);

my @native_c_list = calc_traj_c_list(\@native_cor, \@domain, \@sec_struc);
#print "@{$native_c_list[$#native_c_list]} \n";

if($traj =~ /\.cor$/)
{
  open(DAT, ">${output_dir}/K_${name}.dat")||die("Error: cannot create ${output_dir}/K_${name}.dat\n\n");
  for(my $i = 1; $i <= @domain; $i++)
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
  my @traj_c_list = calc_traj_c_list(\@traj_cor, \@domain, \@sec_struc);

  for(my $i = 0; $i < @traj_c_list; $i++)
  {
    my $K = calc_K($native_c_list[$i], $traj_c_list[$i], $cutoff);
    my $fraction = sprintf("%.4f", $K);
    printf DAT ("%10s ", $fraction);
  }
  printf DAT ("\n");
  exit;
}

open(my $dcdfile, "<$traj") || die("Error: Cannot find $traj\n\n");
binmode($dcdfile);
my ($natom, $tstep, $nframe, $first, $delta, $crystal, $fixed) = read_DCD_head($dcdfile);

if($natom ne $native_natom)
{
  print("Warning: atom numbers mismatch in $input_cor ($native_natom) and $traj ($natom)\n")
}

open(DAT, ">${output_dir}/K_${name}.dat")||die("Error: cannot create ${output_dir}/K_${name}.dat\n\n");
for(my $i = 1; $i <= @domain; $i++)
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

my $deltat = $delta * $tstep;
my $firstframe = $first / $delta;
for(my $i_frame = 1; $i_frame <= $nframe; $i_frame++)
{
  my @traj_cor = read_DCD_frame($dcdfile, $crystal);
  if($i_frame >= $start && $i_frame <= $end)
  {
    my @traj_c_list = calc_traj_c_list(\@traj_cor, \@domain, \@sec_struc);
    #print "@{$traj_c_list[$#traj_c_list]} \n";

    for(my $i = 0; $i < @traj_c_list; $i++)
    {
      my $K = calc_K($native_c_list[$i], $traj_c_list[$i], $cutoff);
      my $fraction = sprintf("%.4f", $K);
      printf DAT ("%10s ", $fraction);
    }
    printf DAT ("\n");
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
    if($line !~ /^#/)
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
sub V_dot
{
  my @V1 = @{$_[0]};
  my @V2 = @{$_[1]};
  
  my $result = 0;
  for(my $i = 0; $i < @V1; $i++)
  {
    $result += $V1[$i] * $V2[$i];
  }
  
  return $result;
}

#########################################################
sub V_cross
{
  my @V1 = @{$_[0]};
  my @V2 = @{$_[1]};
  
  my @result = ();
  $result[0] = $V1[1]*$V2[2]-$V1[2]*$V2[1];
  $result[1] = $V1[2]*$V2[0]-$V1[0]*$V2[2];
  $result[2] = $V1[0]*$V2[1]-$V1[1]*$V2[0];
  
  return @result;
}

#########################################################
sub V_norm
{
  my @V1 = @{$_[0]};
  
  my $result = sqrt(V_dot(\@V1, \@V1));
  
  return $result;
}

#########################################################
sub V_add
{
  my @V1 = @{$_[0]};
  my @V2 = @{$_[1]};
  
  my @result = ();
  for(my $i = 0; $i < @V1; $i++)
  {
    $result[$i] = $V1[$i] + $V2[$i];
  }
  
  return @result;
}

#########################################################
sub V_minus
{
  my @V1 = @{$_[0]};
  my @V2 = @{$_[1]};
  
  my @result = ();
  for(my $i = 0; $i < @V1; $i++)
  {
    $result[$i] = $V1[$i] - $V2[$i];
  }
  
  return @result;
}

#########################################################
sub calc_c_list
{
  my @coor = @{$_[0]};
  my @sel = @{$_[1]};

  my @c_list = ();

  for(my $i = 0; $i < @sel-3; $i++)
  {
    my @v1 = V_minus($coor[$sel[$i+1]], $coor[$sel[$i]]);
    my @v2 = V_minus($coor[$sel[$i+2]], $coor[$sel[$i+1]]);
    my @v3 = V_minus($coor[$sel[$i+3]], $coor[$sel[$i+2]]);
    my @a = V_cross(\@v1, \@v2);
    my $c = V_dot(\@a, \@v3);
    my $cn = V_norm(\@v1) * V_norm(\@v2) * V_norm(\@v3);
    if($cn == 0)
    {
      $c = 'nan';
    }
    else
    {
      # $c = $c / (V_norm(\@v1) * V_norm(\@v2) * V_norm(\@v3));
      $c = $c / $cn;
    }
    push@c_list, $c;
  }

  return @c_list;
}

#########################################################
sub calc_traj_c_list
{
  my @coor = @{$_[0]};
  my @domain = @{$_[1]};
  my @sec_struc = @{$_[2]};
  
  my @result = ();
  foreach my $dom (@domain)
  {
    my @range = @{$dom->{"range"}};
    if($dom->{"class"} eq "i")
    {
      my @range_1 = @{$domain[$range[0]-1]->{"range"}};
      my @range_2 = @{$domain[$range[1]-1]->{"range"}};
      @range = (@range_1, @range_2);
    }
    
    my @sel = ();
    foreach my $rd (@range)
    {
      for(my $i = $rd->[0]; $i <= $rd->[1]; $i++)
      {
        my $tag_i = 0;
        foreach my $rs (@sec_struc)
        {
          #if($i >= $rs->[0] && $i <= $rs->[1])
          if($i == $rs->[0] || $i == $rs->[1])
          {
            $tag_i = 1;
            last;
          }
        }

        if($tag_i eq 1)
        {
          push@sel, $i-1;
        }
      }
    }
    my @c_list = calc_c_list(\@coor, \@sel);
    push@result, \@c_list;
  }
  my @sel = ();
  for(my $i = 1; $i <= $sec_struc[$#sec_struc]->[1]; $i++)
  {
    my $tag_i = 0;
    foreach my $rs (@sec_struc)
    {
      #if($i >= $rs->[0] && $i <= $rs->[1])
      if($i == $rs->[0] || $i == $rs->[1])
      {
        $tag_i = 1;
        last;
      }
    }

    if($tag_i eq 1)
    {
      push@sel, $i-1;
    }
  }
  my @c_list = calc_c_list(\@coor, \@sel);
  push@result, \@c_list;
  return @result;
}

#########################################################
sub calc_K
{
  my @native_c_list = @{$_[0]};
  my @traj_c_list = @{$_[1]};
  my $cutoff = $_[2];
  
  my $K = 0;
  my $not_nan = 0;
  for(my $i = 0; $i < @native_c_list; $i++)
  {
    if($traj_c_list[$i] ne 'nan' && $traj_c_list[$i]*$native_c_list[$i] >= 0)
    {
      $K += 1;
    }
    if($native_c_list[$i] ne 'nan')
    {
      $not_nan += 1;
    }
  }
  if($not_nan == 0)
  {
    $K = -1;
  }
  else
  {
    $K /= $not_nan;
  }
  return $K;
}

#########################################################
sub calc_K_2
{
  my @native_c_list = @{$_[0]};
  my @traj_c_list = @{$_[1]};
  my $cutoff = $_[2];
  
  my $K_0 = 0;
  my $K_1 = 0;
  for(my $i = 0; $i < @native_c_list; $i++)
  {
    $K_0 += $native_c_list[$i];
    $K_1 += $traj_c_list[$i];
  }
  if($K_0 * $K_1 < 0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
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
