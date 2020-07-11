#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;

my ($help, $input_cor, $traj, $start, $end, $mask, $output_dir);
GetOptions(
 'help|h!' => \$help,
 'input|i=s' => \$input_cor,
 'traj|t=s' => \$traj,
 'begin|b=s' => \$start,
 'end|e=s' => \$end,
 'mask|k=s' => \$mask,
 'outdir|o=s' => \$output_dir,
);

my $usage = "
  Usage: perl calc_entanglement_number.pl
              --input | -i <INPUT.COR> for identify native contacts
              --traj | -t <TRAJ.DCD> for simulation trajectory
              [--begin | -b] <START STEP> to split trajectory. 
                                  Default is 1
              [--end | -e] <END STEP> to split trajectory. 
                                Default is the end frame.
              [--mask | -k] <MASK> for a subset of residues to calculate their
                                   native contacts with all other residues. 
                                   Default is 'all'.
                                   Example MASK: 
                                   '1-45' for selection from resid 1 to 45;
                                   '1,3,5-9' for selection of resid 1, 3 
                                   and from 5 to 9.
              [--outdir | -o] <DIRECTORY> for the outputs. Default is the directory
                                          of your trajectories
              [--help | -h]\n\n";

if($help)
{
  die($usage);
}
elsif(!(defined($input_cor) && defined($traj)))
{
  die($usage);
}

if(!defined($start))
{
  $start = 1;
}
if(!defined($end))
{
  $end = 9999999999999999999;
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

my $cutoff = 8;

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

my @native_contact_list = calc_contact_list(\@native_cor, $cutoff);

my @G_native_list = calc_G_list(\@native_cor, \@native_contact_list, $cutoff);
#print "@G_native_list \n";

if($traj =~ /\.cor$/)
{
  open(DAT, ">${output_dir}/G_${name}.dat")||die("Error: cannot create ${output_dir}/G_${name}.dat\n\n");
  printf DAT ("# G1: Gain entanglement and not switch chirality (e.g. 0 to 1)\n");
  printf DAT ("# G2: Gain entanglement and switch chirality (e.g. -1 to 2)\n");
  printf DAT ("# G3: Lose entanglement and not switch chirality (e.g. -2 to -1)\n");
  printf DAT ("# G4: Lose entanglement and switch chirality (e.g. -2 to 1)\n");
  printf DAT ("# G5: Only switch chirality (e.g. -1 to 1)\n");
  printf DAT ("# fGtot: Total number of change of entanglement / Total number of native contacts in this structure\n");
  printf DAT ("# Total native contact: %d\n", $#native_contact_list+1);
  printf DAT ("%6s %6s %6s %6s %6s %6s %6s\n", "max_g", "G1", "G2", "G3", "G4", "G5", "fGtot");
  close(DAT);
  my @traj_cor = parse_cor($traj);
  my $natom = @traj_cor;
  if($natom != $native_natom)
  {
    print("Warning: atom numbers mismatch in $input_cor ($native_natom) and $traj ($natom)\n")
  }
  my @G_traj_list = calc_G_list(\@traj_cor, \@native_contact_list, $cutoff);
  my $G1 = 0; # Gain entanglement and not switch chirality (e.g. 0 to 1)
  my $G2 = 0; # Gain entanglement and switch chirality (e.g. -1 to 2)
  my $G3 = 0; # Lose entanglement and not switch chirality (e.g. -2 to -1)
  my $G4 = 0; # Lose entanglement and switch chirality (e.g. -2 to 1)
  my $G5 = 0; # Only switch chirality (e.g. -1 to 1)
  my $max_g = -999;
  my $min_g = 999;
  my $num_nc_traj = 0;
  for(my $i = 0; $i < @native_contact_list; $i++)
  {
    if(($G_traj_list[$i] ne 'nan'))
    {
      $num_nc_traj += 1;
      if(($G_native_list[$i] * $G_traj_list[$i]) >= 0)
      {
        if(abs($G_native_list[$i]) < abs($G_traj_list[$i]))
        {
          $G1++;
        }
        elsif(abs($G_native_list[$i]) > abs($G_traj_list[$i]))
        {
          $G3++;
        }
      }
      else
      {
        if(abs($G_native_list[$i]) < abs($G_traj_list[$i]))
        {
          $G2++;
        }
        elsif(abs($G_native_list[$i]) > abs($G_traj_list[$i]))
        {
          $G4++;
        }
        elsif(abs($G_native_list[$i]) == abs($G_traj_list[$i]))
        {
          $G5++;
        }
      }
    }
    
    if($max_g < $G_traj_list[$j])
    {
      $max_g = $G_traj_list[$j];
    }
    if($min_g > $G_traj_list[$j])
    {
      $min_g = $G_traj_list[$j];
    }
  }
  my $fGtot = 0;
  if($num_nc_traj != 0)
  {
    $fGtot = ($G1+$G2+$G3+$G4+$G5)/$num_nc_traj;
  }
  if($max_g == -999)
  {
    $max_g = 0;
  }
  if($min_g == 999)
  {
    $min_g = 0;
  }
  open(DAT, ">>${output_dir}/G_${name}.dat");
  printf DAT ("%6s %6d %6d %6d %6d %6d %6.4f\n", "$min_g,$max_g", $G1, $G2, $G3, $G4, $G5, $fGtot);
  close(DAT);
  exit;
}

open(my $dcdfile, "<$traj") || die("Error: Cannot find $traj\n\n");
binmode($dcdfile);
my ($natom, $tstep, $nframe, $first, $delta, $crystal, $fixed) = read_DCD_head($dcdfile);

if($natom ne $native_natom)
{
  print("Warning: atom numbers mismatch in $input_cor ($native_natom) and $traj ($natom)\n")
}

open(DAT, ">${output_dir}/G_${name}.dat")||die("Error: cannot create ${output_dir}/G_${name}.dat\n\n");
printf DAT ("# G1: Gain entanglement and not switch chirality (e.g. 0 to 1)\n");
printf DAT ("# G2: Gain entanglement and switch chirality (e.g. -1 to 2)\n");
printf DAT ("# G3: Lose entanglement and not switch chirality (e.g. -2 to -1)\n");
printf DAT ("# G4: Lose entanglement and switch chirality (e.g. -2 to 1)\n");
printf DAT ("# G5: Only switch chirality (e.g. -1 to 1)\n");
printf DAT ("# fGtot: Total number of change of entanglement / Total number of native contacts in this structure\n");
printf DAT ("# Total native contact: %d\n", $#native_contact_list+1);
printf DAT ("%6s %6s %6s %6s %6s %6s %6s\n", "max_g", "G1", "G2", "G3", "G4", "G5", "fGtot");
close(DAT);

my $deltat = $delta * $tstep;
my $firstframe = $first / $delta;
for(my $i = 1; $i <= $nframe; $i++)
{
  my @traj_cor = read_DCD_frame($dcdfile, $crystal);
  if($i >= $start && $i <= $end)
  {
    my @G_traj_list = calc_G_list(\@traj_cor, \@native_contact_list, $cutoff);
    #print "@G_traj_list \n";
    my $G1 = 0; # Gain entanglement and not switch chirality (e.g. 0 to 1)
    my $G2 = 0; # Gain entanglement and switch chirality (e.g. -1 to 2)
    my $G3 = 0; # Lose entanglement and not switch chirality (e.g. -2 to -1)
    my $G4 = 0; # Lose entanglement and switch chirality (e.g. -2 to 1)
    my $G5 = 0; # Only switch chirality (e.g. -1 to 1)
    my $max_g = -999;
    my $min_g = 999;
    my $num_nc_traj = 0;
    for(my $j = 0; $j < @native_contact_list; $j++)
    {
      if(($G_traj_list[$j] ne 'nan'))
      {
        $num_nc_traj += 1;
        if(($G_native_list[$j] * $G_traj_list[$j]) >= 0)
        {
          if(abs($G_native_list[$j]) < abs($G_traj_list[$j]))
          {
            $G1++;
          }
          elsif(abs($G_native_list[$j]) > abs($G_traj_list[$j]))
          {
            $G3++;
          }
        }
        else
        {
          if(abs($G_native_list[$j]) < abs($G_traj_list[$j]))
          {
            $G2++;
          }
          elsif(abs($G_native_list[$j]) > abs($G_traj_list[$j]))
          {
            $G4++;
          }
          elsif(abs($G_native_list[$j]) == abs($G_traj_list[$j]))
          {
            $G5++;
          }
        }
      }
      
      if($max_g < $G_traj_list[$j])
      {
        $max_g = $G_traj_list[$j];
      }
      if($min_g > $G_traj_list[$j])
      {
        $min_g = $G_traj_list[$j];
      }
    }
    my $fGtot = 0;
    if($num_nc_traj ne 0)
    {
      $fGtot = ($G1+$G2+$G3+$G4+$G5)/$num_nc_traj;
    }
    if($max_g == -999)
    {
      $max_g = 0;
    }
    if($min_g == 999)
    {
      $min_g = 0;
    }
    open(DAT, ">>${output_dir}/G_${name}.dat");
    printf DAT ("%6s %6d %6d %6d %6d %6d %6.4f\n", "$min_g,$max_g", $G1, $G2, $G3, $G4, $G5, $fGtot);
    close(DAT);
  }
}
close(DAT);

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

#########################################################
sub calc_contact_list
{
  my @coor = @{$_[0]};
  my $cutoff = $_[1];

  my @c_map = ();

  for(my $i = 0; $i < @coor-4; $i++)
  {
    for(my $j = $i + 4; $j < @coor; $j++)
    {
      my $distance = ($coor[$i]->[0] - $coor[$j]->[0]) ** 2 + ($coor[$i]->[1] - $coor[$j]->[1]) ** 2 + ($coor[$i]->[2] - $coor[$j]->[2]) ** 2;
      $distance = $distance ** 0.5;

      if($distance <= $cutoff)
      {
        push@c_map, [$i, $j];
      }
    }
  }

  return @c_map;
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
sub calc_G_list
{
  my @coor = @{$_[0]};
  my @sel = @{$_[1]};
  my $cutoff = $_[2];
  my $terminal_cutoff = 5;
  
  my $n_atom = @coor;
  # Generate contact matrix
  my @R = ();
  my @dR = ();
  for(my $i = 0; $i < $n_atom-1; $i++)
  {
    push@R, [($coor[$i]->[0]+$coor[$i+1]->[0])/2, ($coor[$i]->[1]+$coor[$i+1]->[1])/2, ($coor[$i]->[2]+$coor[$i+1]->[2])/2];
    push@dR, [$coor[$i+1]->[0]-$coor[$i]->[0], $coor[$i+1]->[1]-$coor[$i]->[1], $coor[$i+1]->[2]-$coor[$i]->[2]];
  }
  my @M = ();
  for(my $i = 0; $i < $n_atom-1; $i++)
  {
    $M[$i] = [];
    for(my $j = 0; $j < $n_atom-1; $j++)
    {
      $M[$i]->[$j] = 0;
    }
  }
  
  for(my $i = 0; $i < $n_atom-2; $i++)
  {
    for(my $j = $i+1; $j < $n_atom-1; $j++)
    {
      my @v1 = ($R[$i]->[0]-$R[$j]->[0], $R[$i]->[1]-$R[$j]->[1], $R[$i]->[2]-$R[$j]->[2]);
      my $v1_n = ($v1[0]**2+$v1[1]**2+$v1[2]**2)**(3/2);
      @v1 = ($v1[0]/$v1_n, $v1[1]/$v1_n, $v1[2]/$v1_n);
      my @v2 = V_cross($dR[$i], $dR[$j]);
      $M[$i]->[$j] = V_dot(\@v1, \@v2);
      $M[$j]->[$i] = $M[$i]->[$j];
    }
  }

  # Calculate G
  my @G_list = ();
  for(my $i = 0; $i < @sel; $i++)
  {
    my $r1_i = $sel[$i]->[0];
    my $r1_j = $sel[$i]->[1];
    my @coor_1 = @{$coor[$r1_i]};
    my @coor_2 = @{$coor[$r1_j]};
    my $dist = (($coor_1[0] - $coor_2[0]) ** 2 + ($coor_1[1] - $coor_2[1]) ** 2 + ($coor_1[2] - $coor_2[2]) ** 2)**(1/2);
    if($dist <= $cutoff)
    {
      my @G0 = (0, 0);
      my @r2_range = ([$terminal_cutoff, $r1_i-4], [$r1_j+4, $n_atom-$terminal_cutoff-1]);
      for(my $j = 0; $j < @r2_range; $j++)
      {
        my $r2_i = $r2_range[$j]->[0];
        my $r2_j = $r2_range[$j]->[1];
        for(my $r1 = $r1_i; $r1 < $r1_j; $r1++)
        {
          for(my $r2 = $r2_i; $r2 < $r2_j; $r2++)
          {
            $G0[$j] += $M[$r1]->[$r2];
          }
        }
      }
      ####
      for(my $k = 0; $k < 2; $k++)
      {
        $G0[$k] = sprintf('%.0f',$G0[$k]/4/3.14);
      }
      #my $G0 = ($G0[0]+$G0[1])/4/3.14;
      push@G_list, ($G0[0]+$G0[1]);
    }
    else
    {
      push@G_list, 'nan';
    }
  }

  return @G_list
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
