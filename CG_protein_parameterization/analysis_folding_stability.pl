#!/usr/bin/perl

use Getopt::Long;
use threads;

my ($help, $input_file, $is_ca, $target_temp, $wham_exec, $start);
GetOptions(
 'help|h!' => \$help,
 'input|i=s' => \$input_file,
 'model|m=s' => \$is_ca,
 'temperature|t=s' => \$target_temp,
 'exec|e=s' => \$wham_exec,
 'start|s=s' => \$start,
);

if($help || !(defined($input_file) && defined($is_ca) && defined($target_temp)))
{
  die("
  Usage: perl analysis_folding_stability.pl 
              --input | -i <INPUT.log> from REX simulation
              --model | -m <ca or cacb> ca for CA cg model; cacb for CA-CB cg model
              --temperature | -t <TEMPERATURE> for stability calculation;
                                 0 means skip stability calculation and only
                                 calculate CV.
              [--exec | -e <EXECUTABLE>] for running wham. 
                                         Can be 'wham_general_v1.29.pl' or 'pywham.py'.
                                         Default is 'wham_general_v1.29.pl'.
              [--start | -s] <STEP> start for trajectories selection. Default is 1.
              [--help | -h]\n\n");
}

if(!defined($start))
{
  $start = 1;
}
my $num_restart = 1;

if(!defined($wham_exec))
{
  $wham_exec = "wham_general_v1.29.pl";
}

my $start_time = time();

mkdir("analysis");
system("echo '' > analysis.log");

print "-> Parsing information from $input_file...\n";
my ($nwindows, $temps, $prefix, $param_file, $tot_prod) = parse_rex_info($input_file);
print "   Done. Number of windows: $nwindows\n";

if($start >= $tot_prod)
{
  die("Error: Start step $start >= total step number $tot_prod");
}

merge_restart($nwindows);

print "-> get Ep from dcd...\n";
Get_Ep($prefix, $param_file, $nwindows);
print "   All Done.";

print "-> Printing REX energy...\n";
my $max_ene = -9999999999999;
my $min_ene = 9999999999999;
for(my $i = 1; $i <= $nwindows; $i++)
{
  open(D1, "<aa$i/ene_1.log");
  open(D2, ">aa$i/ene.dat");
  print D2 "nsteps   energy\n";
  my $n = 1;
  while(my $line = <D1>)
  {
    if($n >= $start)
    {
      print D2 $n." ".$line;
      chomp($line);
      my @str = split(/\s+/, $line);
      if($str[$#str] > $max_ene)
      {
        $max_ene = $str[$#str];
      }
      elsif($str[$#str] < $min_ene)
      {
        $min_ene = $str[$#str];
      }
    }
    $n++;
  }
  close(D1);
  close(D2);
}
$avg_ene = ($max_ene + $min_ene) / 2;
print "   Done. minE is $min_ene, maxE is $max_ene, average is $avg_ene\n";

my $wham_folder = "wham_analysis";
if($wham_exec eq "pywham.py")
{
  $wham_folder = "pywham_analysis";
}

chdir("analysis");
mkdir($wham_folder);
if($wham_exec eq "pywham.py" || $wham_exec eq "mix")
{
  open(OUT,">${wham_folder}/myreaders.py");
  print OUT "def ReadLastColumn(filename):
    QList = []
    results = []
    n = 0
    f = open(filename, \"r\")
    for line in f:
        if n == 0:
            n = n + 1
            continue
        sline = line.strip().split()
        QList.append(float(sline[-1]))
    f.close()
    results.append(QList)
    return results";
  close(OUT);
}

if($target_temp eq 0)
{
  chdir($wham_folder);
  
  my $melting_t = 0;

  if($wham_exec eq "mix")
  {
  	print "-> Pre-calculate WHAM equation via pywham...\n";
  	$melting_t = Call_wham_cv("pywham.py", $temps);
  }
  else
  {
  	print "-> Calculating CV via WHAM...\n";
    $melting_t = Call_wham_cv($wham_exec, $temps);
  }
  print "   Done. Melting temperature is $melting_t K\n";
}
else
{
  mkdir("final_ts");

  print "-> Calculating properties...\n";
  Calc_props($prefix, $nwindows, $start);
  print "   All Done\n";

  chdir($wham_folder);
  
  my $melting_t = 0;

  if($wham_exec eq "mix")
  {
  	print "-> Calculating CV via pywham...\n";
  	$melting_t = Call_wham_cv("pywham.py", $temps);
  	print "-> Re-calculate free energy via WHAM...\n";
    Call_wham_cv($wham_exec, $temps);
  }
  else
  {
  	print "-> Calculating CV via WHAM...\n";
    $melting_t = Call_wham_cv($wham_exec, $temps);
  }
  print "   Done. Melting temperature is $melting_t K\n";
  print "-> Calculating Q threshold where the cumulative probability equals 0.5 at the $melting_t K...\n";
  my $Q_threshold = Find_Q_threshold($wham_exec, $temps, $melting_t, 0);
  print "   Done. Q threshold is $Q_threshold\n";
  print "-> Calculating Q310 where the probability of Q reaches maximun at the 310 K...\n";
  my $Q_310 = Find_Q_310($wham_exec, $temps, '310.0', 0);
  print "   Done. Q310 is $Q_310\n";
  print "-> Calculating the probability of being in the native state at the $target_temp K...\n";
  my $probability = Find_probability($wham_exec, $temps, $target_temp, $Q_threshold);
  print "   Done. Probability is $probability\n";

  my $stability = `dGns.pl $target_temp $probability 2>&1`;
  if($stability !~ /kcal\/mol at/)
  {
    die("Error: Failed to calculate stability using probability $probability\n\n");
  }
  else
  {
    print "-> Protein stability is $stability";
  }
}

my $end_time = time();
my $cost_time = $end_time - $start_time;
print "-> Total time usage: " . convert_time($cost_time) . "\n";

#####################################################################
sub parse_rex_info
{
  my $info_file = $_[0];

  my $nwindows;
  my @temps;
  my $tot_prod;
  my $prefix;
  my $param_file;

  open(INFO, "<$info_file") || die("Error: Cannot find $info_file\n\n");
  while(my $line = <INFO>)
  {
  	chomp($line);
  	if($line =~ /^Number of windows:/)
  	{
  	  my @str = split(/\:/, $line);
  	  $nwindows = $str[1];
  	  $nwindows =~ s/^\s+|\s+$//g;
  	}
  	elsif($line =~ /^Temperatures:/)
  	{
  	  my @str = split(/:/, $line);
  	  my $t = $str[1];
  	  $t =~ s/^\s+|\s+$//g;
  	  @temps = split(/\s+/, $t);
  	}
  	elsif($line =~ /^Number of exchanges in production:/)
  	{
  	  my @str = split(/:/, $line);
  	  $tot_prod = $str[1];
  	  $tot_prod =~ s/^\s+|\s+$//g;
  	}
  	elsif($line =~ /^Charmm psf file:/)
  	{
  	  my @str = split(/:/, $line);
  	  my $psf = $str[1];
  	  $psf =~ s/^\s+|\s+$//g;
  	  @str = split(/.psf/, $psf);
  	  $prefix = $str[0];
  	}
  	elsif($line =~ /^Charmm prm file:/)
  	{
  	  my @str = split(/:/, $line);
  	  $param_file = $str[1];
  	  $param_file =~ s/^\s+|\s+$//g;
  	}
  }

  my @result = ($nwindows, \@temps, $prefix, $param_file, $tot_prod);

  return @result;
}

#####################################################################
sub merge_restart
{
  my $nwindows = $_[0];
  for(my $i = 1; $i <= $nwindows; $i++)
  {
    my $list = `ls aa$i/ene_r* 2>&1`;
    if($list =~ /No such file or directory/)
    {
      next;
    }
    my @ene_list = split(/\n/, $list);
    foreach my $l (@ene_list)
    {
      system("cat $l >> aa$i/ene_1.log");
      system("rm -f $l");
    }
  }
}

#####################################################################
sub Get_Ep
{
  my $prefix = $_[0];
  my $param_file = $_[1];
  my $nwindows = $_[2];
  my @str = split(/\.prm/, $param_file);
  my $xml_file = $str[0].".xml";
  my @threads = ();
  for(my $i = 1; $i <= $nwindows; $i++)
  {
    my $code = "cd aa${i}/; get_Ep_from_dcd.py -p ../$prefix.psf -c mc1.dcd -x ../$xml_file -o ene_1.log -m ':*' >/dev/null";
    my $t = threads->create('run_exe', $code);
    push@threads, $t;
    $t->detach();
  }
  
  while()
  {
    sleep(10);
    my $tag = 0;
    foreach my $t (@threads)
    {
      if($t->is_running())
      {
        $tag = 1;
        last;
      }
    }
    if($tag eq 0)
    {
      last;
    }
  }
}

#####################################################################
sub Calc_props
{
  my $prefix = $_[0];
  my $nwindows = $_[1];
  my $start = $_[2];

  my $natom = 0;
  my $nresidue = 0;
  
  open(PSF, "<../$prefix.psf") || die("Error: Cannot find $prefix.psf\n\n");
  while(my $line = <PSF>)
  {
  	chomp($line);
    if($line =~ /!NATOM/)
    {
      my @str = split(/!/, $line);
  	  $natom = $str[0];
  	  $natom =~ s/^\s+|\s+$//g;
      
      my $n = 0;
  	  while($n < $natom)
  	  {
  	  	$line = <PSF>;
  	  	$n++;
  	  }
  	  $nresidue = substr($line, 14, 4);
  	  $nresidue =~ s/^\s+|\s+$//g;
    }
  }
  close(PSF);
  
  open(DOM, ">domain_def.dat") || die("Error: Cannot create domain_def.dat\n\n");
  print DOM "1:$nresidue a\n";
  close(DOM);
  
  my @threads = ();
  for(my $i = 1; $i <= $nwindows; $i++)
  {
    my $t = threads->create('run_exe', "calc_native_contact_fraction.pl -i ../$prefix.cor -d domain_def.dat -s ../setup/secondary_struc_defs.txt -t ../aa${i}/mc1.dcd -b $start >> ../analysis.log 2>&1");
    push@threads, $t;
    $t->detach();
  }
  
  while()
  {
    sleep(10);
    my $tag = 0;
    foreach my $t (@threads)
    {
      if($t->is_running())
      {
        $tag = 1;
        last;
      }
    }
    if($tag eq 0)
    {
      last;
    }
  }

  for(my $i = 1; $i <= $nwindows; $i++)
  {
    chdir("../aa$i/");
    system("mv qbb_mc1.dat ../analysis/final_ts/aa${i}_props_vs_time.dat");
    chdir("../analysis/");
  }
}

#####################################################################
sub run_exe
{
  my $cmd = $_[0];
  system($cmd);
}

#####################################################################
sub Call_wham_cv
{
  my $wham_exec = $_[0];
  my @temps = @{$_[1]};
  
  my $melting_t = 0;

  my $max_temp = $temps[$#temps];
  $max_temp += 20;
  my $min_temp = $temps[0];
  $min_temp -= 20;
  my $nsim_temps = @temps;
  
  if($wham_exec eq "pywham.py")
  {
    open(CNTRL, ">wham_cv.xml") || die("Error: Cannot create analysis/pywham_analysis/wham_cv.xml\n\n");
    print CNTRL "<?xml version=\"1.0\" ?>
<WhamSpec>
    <General>
        <Coordinates>
            <Coordinate name=\"V\" />
        </Coordinates>
        <CoordinateFileReaders pythonModule=\"myreaders\" returnsTime=\"false\">
            <Reader name=\"ReadLastColumn\">
                <ReturnList name=\"V\" />
            </Reader>
        </CoordinateFileReaders>
        <Binnings>
            <Binning name=\"V\">
                <Interval>1</Interval>
            </Binning>
        </Binnings>
        <Parameters>
            <Parameter name =\"reference_energy\">$avg_ene</Parameter>
            <Parameter name=\"kB\">0.0019872041</Parameter>
        </Parameters>
    </General>
    <Trajectories>\n";
    for(my $j=1;$j<=$nsim_temps;$j++)
    {
      print CNTRL " " x 8 . "<Trajectory T=\"" . $temps->[$j-1] . "\">
            <EnergyFunction>V</EnergyFunction>
            <CoordinateFiles>
                <CoordinateFile>../../aa${j}/ene.dat</CoordinateFile>
            </CoordinateFiles>
        </Trajectory>\n";
    }
    print CNTRL "    </Trajectories>
    <Jobs>
        <HeatCapacity outFile=\"cv/cv.dat\">
            <EnergyFunction>V</EnergyFunction>
            <Temperatures>${min_temp}:0.1:${max_temp}</Temperatures>
        </HeatCapacity>
    </Jobs>
</WhamSpec>";
    close(CNTRL);
    
    system("$wham_exec wham_cv.xml  >> ../../analysis.log 2>&1");
    
    open(DAT, "<cv/cv.dat") || die("Error: Cannot find analysis/pywham_analysis/cv/cv.dat\n\n");
    my $cv_max = 0;
    while(my $line = <DAT>)
    {
      chomp($line);
      $line =~ s/^\s+|\s+$//g;
      my ($t, $cv) = split(/\s+/, $line);
      if($cv > $cv_max)
      {
  	    $cv_max = $cv;
        $melting_t = $t;
      }
    }
    close(DAT);
  }
  else
  {
    mkdir("cv");

    open(OUT,">wham_ene_files.txt");
    for(my $j=1;$j<=$nsim_temps;$j++) 
    {
      print OUT "../../aa$j/ene.dat " . $temps->[$j-1] . "\n";
    }
    close OUT;
    
    open(TEMP, "<../../setup/wham_cv_template.cntrl") || die("Error: Cannot find setup/wham_cv_template.cntrl\n\n");
    open(CNTRL, ">wham_cv.cntrl") || die("Error: Cannot create analysis/wham_analysis/wham_cv.cntrl\n\n");
    while(my $line = <TEMP>)
    {
      if($line =~ /^nsim_temps /)
      {
  		  print CNTRL "nsim_temps = $nsim_temps\n";
      }
      elsif($line =~ /^tmin /)
      {
  		  print CNTRL "tmin = $min_temp\n";
      }
      elsif($line =~ /^tmax /)
      {
  		  print CNTRL "tmax = $max_temp\n";
      }
      elsif($line =~ /^ppn /)
      {
  		  print CNTRL "ppn = $nsim_temps\n";
      }
      elsif($line =~ /^initial_guess / && $wham_exec eq "mix")
      {
      	print CNTRL "initial_guess = yes\n";
      	print CNTRL "freenrg_file = free_energy.dat\n";
      }
      elsif($line =~ /^heat_capacity / && $wham_exec eq "mix")
      {
      	print CNTRL "heat_capacity = off\n";
      }
      elsif($line =~ /^dt /)
      {
  		  print CNTRL "dt = 0.1\n";
      }
      elsif($line =~ /^skip_first_nlines /)
      {
        # do nothing
      }
      else
      {
  		  print CNTRL $line;
      }
    }
    print CNTRL "skip_first_nlines = 1\n";
    close(TEMP);
    close(CNTRL);

    if($wham_exec eq "mix")
    {
      system("wham_general_v1.29.pl wham_cv.cntrl  >> ../../analysis.log 2>&1");
      return 0;
    }
    else
    {
      system("$wham_exec wham_cv.cntrl  >> ../../analysis.log 2>&1");
    }

    open(DAT, "<cv/prot_cv_vs_temp.dat") || die("Error: Cannot find analysis/wham_analysis/cv/prot_cv_vs_temp.dat\n\n");
    my $cv_max = 0;
    while(my $line = <DAT>)
    {
      chomp($line);
      my ($t, $cv) = split(/\s+/, $line);
      if($cv > $cv_max)
      {
  		$cv_max = $cv;
        $melting_t = $t;
      }
    }
    close(DAT);
  }
  
  $melting_t *= 1;
  if($melting_t !~ /\./)
  {
  	$melting_t .= ".0";
  }
  return $melting_t;
}

#####################################################################
sub Call_wham_qbb
{
  my $wham_exec = $_[0];
  my @temps = @{$_[1]};
  my $melting_t = $_[2];
  my $Q_threshold = $_[3];
  
  my $nsim_temps = @temps;

  if($wham_exec eq "pywham.py")
  {
    open(CNTRL, ">wham_qbb.xml") || die("Error: Cannot create analysis/pywham_analysis/wham_qbb.xml\n\n");
    print CNTRL "<?xml version=\"1.0\" ?>
<WhamSpec>
    <General>
        <Coordinates>
            <Coordinate name=\"V\" />
            <Coordinate name=\"Q\" />
        </Coordinates>
        <CoordinateFileReaders pythonModule=\"myreaders\" returnsTime=\"false\">
            <Reader name=\"ReadLastColumn\">
                <ReturnList name=\"V\" />
            </Reader>
            <Reader name=\"ReadLastColumn\">
                <ReturnList name=\"Q\" />
            </Reader>
        </CoordinateFileReaders>
        <Binnings>
            <Binning name=\"V\">
                <Interval>0.0001</Interval>
            </Binning>
            <Binning name=\"Q\">
                <Begin>0</Begin>
                <End>1.02</End>
                <Interval>";
    if($Q_threshold)
    {
      if($Q_threshold >= 0.5)
      {
      	print CNTRL $Q_threshold/8;
      }
      else
      {
      	print CNTRL $Q_threshold/4;
      }
    }
    else
    {
      print CNTRL "0.04";
    }
    print CNTRL "</Interval>
            </Binning>
        </Binnings>
        <Parameters>
            <Parameter name =\"reference_energy\">$avg_ene</Parameter>
            <Parameter name=\"kB\">0.0019872041</Parameter>
        </Parameters>
    </General>
    <Trajectories>\n";
    for(my $j=1;$j<=$nsim_temps;$j++)
    {
      print CNTRL " " x 8 . "<Trajectory T=\"" . $temps->[$j-1] . "\">
            <EnergyFunction>V</EnergyFunction>
            <CoordinateFiles>
                <CoordinateFile>../../aa${j}/ene.dat</CoordinateFile>
                <CoordinateFile>../final_ts/aa${j}_props_vs_time.dat</CoordinateFile>
            </CoordinateFiles>
        </Trajectory>\n";
    }
    print CNTRL "    </Trajectories>
    <Jobs>
        <FreeEnergy outFilePrefix=\"qbb/fe_${Q_threshold}_\">
            <Coordinates>
                <Coordinate name=\"Q\" />
            </Coordinates>
            <EnergyFunction>V</EnergyFunction>
            <Temperatures>$melting_t</Temperatures>
            <Parameters>
                <Parameter name=\"in_kT\">false</Parameter>
            </Parameters>
        </FreeEnergy>
    </Jobs>
</WhamSpec>";
    close(CNTRL);
    
    system("$wham_exec wham_qbb.xml  >> ../../analysis.log 2>&1");
  }
  else
  {
    open(TEMP, "<../../setup/wham_qbb_template.cntrl") || die("Error: Cannot find setup/wham_qbb_template.cntrl\n\n");
    open(CNTRL, ">wham_qbb.cntrl") || die("Error: Cannot create analysis/wham_analysis/wham_qbb.cntrl\n\n");
    while(my $line = <TEMP>)
    {
      if($line =~ /^nsim_temps /)
      {
        print CNTRL "nsim_temps = $nsim_temps\n";
      }
      elsif($line =~ /^tmin /)
      {
        print CNTRL "tmin = $melting_t\n";
      }
      elsif($line =~ /^tmax /)
      {
        print CNTRL "tmax = $melting_t\n";
      }
      elsif($line =~ /^skip_first_nlines /)
      {
        # do nothing
      }
      elsif($line =~ /^cond_file =/)
      {
        if($Q_threshold eq 0)
        {
          last;
        }
        else
        {
          print CNTRL $line;
        }
      }
      elsif($line =~ /^cond_val1 =/)
      {
        print CNTRL "cond_val1 = $Q_threshold\n";
      }
      else
      {
        print CNTRL $line;
      }
    }
    print CNTRL "skip_first_nlines = 1\n";
    close(TEMP);
    close(CNTRL);
    
    if($wham_exec eq "mix")
    {
      system("wham_general_v1.29.pl wham_qbb.cntrl >> ../../analysis.log 2>&1");
    }
    else
    {
      system("$wham_exec wham_qbb.cntrl >> ../../analysis.log 2>&1");
    }
  }
}

#####################################################################
sub Find_Q_threshold
{
  my ($wham_exec, $temps, $target_temp, $Q_threshold) = @_;
  
  if($wham_exec ne "pywham.py")
  {
    mkdir("qbb");    

    open(OUT,">wham_props_files.txt");
    for(my $j=1;$j<=@{$temps};$j++) 
    {
      print OUT "../final_ts/aa${j}_props_vs_time.dat " . $temps->[$j-1] . "\n";
    }
    close OUT;
  }

  Call_wham_qbb($wham_exec, $temps, $target_temp, $Q_threshold);
  
  if($wham_exec eq "pywham.py")
  {
    open(DAT, "<qbb/fe_${Q_threshold}_T${target_temp}") || die("Error: Cannot find analysis/pywham_analysis/qbb/fe_${Q_threshold}_T${target_temp}\n\n");
  }
  else
  {
    $target_temp *= 1;
    open(DAT, "<qbb/qbb_prob_vs_qbb_${target_temp}k.dat") || die("Error: Cannot find analysis/wham_analysis/qbb/qbb_prob_vs_qbb_${target_temp}k.dat\n\n");
  }

  my @P = ();
  my @Q = ();
  while(my $line = <DAT>)
  {
    chomp($line);
    $line =~ s/^\s+|\s+$//g;
    my @str = split(/\s+/, $line);
    my $q = $str[0];
    my $p = $str[1];
    push@P, $p;
    push@Q, $q;
  }
  close(DAT);

  my @c_P = ();
  my $sum = 0;
  foreach my $p (@P)
  {
    $sum += $p;
    push@c_P, $sum;
  }

  my $index = 0;

  for(my $i = 0; $i < @Q; $i++)
  {
    if($c_P[$i] > 0.5)
    {
      $index = $i;
      last;
    }
    elsif($c_P[$i] eq 0.5)
    {
      return $Q[$i];
    }
  }

  my $q1 = $Q[$index-1];
  my $q2 = $Q[$index];
  my $p1 = $c_P[$index-1];
  my $p2 = $c_P[$index];

  my $k = ($p2-$p1) / ($q2-$q1);
  my $b = $p1 - $k * $q1;

  my $Q_eq = (0.5 - $b) / $k;

  return $Q_eq;
}

#####################################################################
sub Find_Q_310
{
  my ($wham_exec, $temps, $target_temp, $Q_threshold) = @_;
  Call_wham_qbb($wham_exec, $temps, $target_temp, $Q_threshold);
  
  if($wham_exec eq "pywham.py")
  {
    open(DAT, "<qbb/fe_${Q_threshold}_T${target_temp}") || die("Error: Cannot find analysis/pywham_analysis/qbb/fe_${Q_threshold}_T${target_temp}\n\n");
  }
  else
  {
    $target_temp *= 1;
    open(DAT, "<qbb/qbb_prob_vs_qbb_${target_temp}k.dat") || die("Error: Cannot find analysis/wham_analysis/qbb/qbb_prob_vs_qbb_${target_temp}k.dat\n\n");
  }

  my $probability = 0;
  my $Q_310 = 0;
  while(my $line = <DAT>)
  {
    chomp($line);
    $line =~ s/^\s+|\s+$//g;
    my @str = split(/\s+/, $line);
    my $q = $str[0];
    my $p = $str[1];
    if($p > $probability)
    {
      $probability = $p;
      $Q_310 = $q * 1;
    }
  }
  close(DAT);

  return $Q_310;
}

#####################################################################
sub Find_probability
{
  my ($wham_exec, $temps, $target_temp, $Q_threshold) = @_;
  Call_wham_qbb($wham_exec, $temps, $target_temp, $Q_threshold);

  my $probability = 0;
  if($wham_exec eq "pywham.py")
  {
  	$target_temp *= 1;
    if($target_temp !~ /\./)
    {
  	  $target_temp .= ".0";
    }
    open(DAT, "<qbb/fe_${Q_threshold}_T${target_temp}") || die("Error: Cannot find analysis/pywham_analysis/qbb/fe_${Q_threshold}_T${target_temp}\n\n");
    while(my $line = <DAT>)
    {
      chomp($line);
      $line =~ s/^\s+|\s+$//g;
      my @str = split(/\s+/, $line);
      my $q = $str[0];
      my $p = $str[1];
      if($q > $Q_threshold)
      {
        $probability += $p;
      }
    }
    close(DAT);
  }
  else
  {
    open(DAT, "<qbb/qbb_conditional_probability_vs_temp.dat") || die("Error: Cannot find analysis/wham_analysis/qbb/qbb_conditional_probability_vs_temp.dat\n\n");
    while(my $line = <DAT>)
    {
      chomp($line);
      my @str = split(/\s+/, $line);
      my $t = $str[0];
      my $p = $str[1];
      if($t eq $target_temp)
      {
        $probability = $p * 1;
      }
    }
    close(DAT);
  }

  return $probability;
}

#####################################################################
sub convert_time
{
  my $time = $_[0];
  my $hours = int($time / 3600);
  my $minutes = int(($time - $hours * 3600) / 60);
  my $seconds = sprintf("%.2f", $time - $hours * 3600 - $minutes * 60);
  return "$hours:$minutes:$seconds";
}
