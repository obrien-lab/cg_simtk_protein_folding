#!/usr/bin/perl

use Getopt::Long;

my ($help, $input_file, $is_ca, $nscal_set, $temp_set, $restart, $opt_temp, $if_analyze, $allocation, $walltime);
GetOptions(
 'help|h!' => \$help,
 'input|i=s' => \$input_pdb,
 'model|m=s' => \$is_ca,
 'nscalset|n=s' => \$nscal_set,
 'tempset|t=s' => \$temp_set,
 'restart|r=s' => \$restart,
 'opttemp|o=s' => \$opt_temp,
 'unfolded|u=s' => \$if_unfolded,
 'analysis|x=s' => \$if_analyze,
 'allocation|a=s' => \$allocation,
 'walltime|w=s' => \$walltime,
);

my $usage = "
  Usage: perl scan_nscal.pl 
              --input | -i <INPUT.PDB> for CG model creation
              --model | -m <ca or cacb> ca for CA cg model;
                           cacb for CA-CB cg model
              --nscalset | -n <\"NUM1 NUM2 ...\"> for nscal scan
              --tempset | -t <\"NUM1 NUM2 ...\"> for temperature windows
              [--opttemp | -o] <1 or 0> for temperature optimization starting 
                                        from native state.
                               Default is 0: not optimize temperature distribution.
                               Cannot be used in restart calculation.
              [--unfolded | -u] <1 or 0> for REX starting from unfolded state.
                                Default is 0: starts from native state.
              [--restart | -r] <RESTART_STEPS> for restart. Must be an integer.
                               Must specify -n. Default is not restart.
              [--analysis | -x] <TARGET_T> for analysis. Defalt is 0, not run analysis.
              [--allocation | -a] <\"ACCOUNT\"> for ACI cluster. 
                                  Default is \"epo2_a_g_bc_default\"
              [--walltime | -w] <\"TIME\"> for ACI cluster. 
                                Default is \"10:00:00:00\"
              [--help | -h]\n
  Example for T windows optimization: perl scan_nscal.pl -i ***.pdb -m ca -n \"1.2 1.3 1.4\" -o 1 -a \"open\" -w \"2:00:00:00\"\n
  Example for running nscal scan: perl scan_nscal.pl -i ***.pdb -m ca -n \"1.2 1.3 1.4\" -t \"290 300 310 320 330 340 350 360\" -a \"open\" -w \"2:00:00:00\"\n
  Example for restarting nscal scan at step 12001: perl scan_nscal.pl -m ca -n \"1.2\" -r 12001 -a \"open\" -w \"2:00:00:00\"\n
  Example for analysis at 298K: perl scan_nscal.pl -m ca -n \"1.2 1.3 1.4\" -x 298 -a \"open\" -w \"2:00:00:00\"\n\n";

if($help)
{
  die($usage);
}
elsif(defined($restart))
{
  if($restart eq 0)
  {
    if(!(defined($input_pdb) && defined($is_ca) && defined($nscal_set) && defined($temp_set)))
    {
      die($usage);
    }
  }
  elsif($restart > 0)
  {
    if(!defined($nscal_set) || $opt_temp eq 1)
    {
      die($usage);
    }
  }
  else
  {
    die($usage);
  }
}
elsif($opt_temp eq 1)
{
  if(!(defined($input_pdb) && defined($is_ca) && defined($nscal_set)))
  {
    die($usage);
  }
}
elsif($if_analyze)
{
  if(!defined($is_ca) || !defined($nscal_set))
  {
    die($usage);
  }
}
elsif(!(defined($input_pdb) && defined($is_ca) && defined($nscal_set) && defined($temp_set)))
{
  die($usage);
}


if(!defined($opt_temp))
{
  $opt_temp = 0;
}
elsif($opt_temp ne 0 && $opt_temp ne 1)
{
  die($usage);
}

if(!defined($if_analyze))
{
  $if_analyze = 0;
}

if(!defined($allocation))
{
  $allocation = "cyberlamp";
}
if(!defined($walltime))
{
  $walltime = "10:00:00:00";
}

if(!defined($if_unfolded))
{
  $if_unfolded = 0;
}

$nscal_set =~ s/^\s+|\s+$//g;
my @nscal_array = split(/\s+/, $nscal_set);

if($restart > 0)
{
  my $opt_str = "";
  if($if_unfolded eq 1)
  {
    $opt_str = "unfolded_";
  }
  foreach my $nscal (@nscal_array)
  {
    my $pdb_name = "";
    chdir("${opt_str}data/nscal_$nscal");
    open(TEMP, "<mubrex.cntrl") || die("Error: Cannot find mubrex.cntrl\n\n");
    open(INP, ">mubrex_restart_$restart.cntrl") || die("Error: Cannot create mubrex_restart_$restart.cntrl\n\n");
    my $n = 0;
    while(my $line = <TEMP>)
    {
      if($line =~ /^nwindows =/)
      {
        print INP $line;
        chomp($line);
        my @str = split(/=/, $line);
        $n = $str[1];
        $n =~ s/^\s+|\s+$//g;
      }
      elsif($line =~ /^psf =/)
      {
        print INP $line;
        chomp($line);
        my @str = split(/=\s+/, $line);
        @str = split(/\//, $str[1]);
        @str = split(/\.psf/, $str[$#str]);
        @str = split(/_/, $str[0]);
        $pdb_name = $str[0];
      }
      elsif($line =~ /^nexch_prod =/)
      {
        chomp($line);
        my @str = split(/=\s+/, $line);
        my $old_prod_n = $str[1];
        $old_prod_n =~ s/^\s+|\s+$//g;
        my $new_prod_n = $old_prod_n;
        if($restart > $old_prod_n)
        {
          $new_prod_n = (int($restart / $old_prod_n) + 1) * $old_prod_n;
        }
        print INP "nexch_prod = $new_prod_n\n";
        print "-> Set nexch_prod = $new_prod_n\n";
      }
      elsif($line =~ /^starting_strucs/)
      {
        print INP "restart = 1\n";
        print INP "nsteps_start = $restart\n";
        last;
      }
      else
      {
        print INP $line;
      }
    }
    close(TEMP);
    close(INP);

    open(TEMP, "<rex.pbs") || die("Error: Cannot find rex.pbs\n\n");
    open(INP, ">rex_restart_$restart.pbs") || die("Error: Cannot create rex_restart_$restart.pbs\n\n");
    while(my $line = <TEMP>)
    {
      if($line =~ /mubrex\.cntrl/)
      {
        print INP "parallel_temperature_REX.py -f mubrex_restart_$restart.cntrl\n";
      }
      elsif($line !~ /run_REX_LD\.py/)
      {
        print INP $line;
      }
    }
    close(TEMP);
    close(INP);

    print "-> Going to restart REX calculations at step $restart for nscal=$nscal\n";
    system("qsub -N n_${pdb_name}_${nscal}_r -l nodes=1:ppn=".($n+1)." -l walltime=$walltime -f rex_restart_$restart.pbs -A $allocation");
    chdir("../../");
  }
}
elsif($if_analyze eq 0)
{
  if($opt_temp eq 1 && !defined($temp_set))
  {
    $temp_set = gen_temp(280, 400, 8);
  }
  my $opt_str = "";
  if($opt_temp eq 1)
  {
    $opt_str = "opt_";
  }
  elsif($if_unfolded eq 1)
  {
    $opt_str = "unfolded_";
  }

  $temp_set =~ s/^\s+|\s+$//g;
  my @temp_array = split(/\s+/, $temp_set);
  if(@temp_array >= 20)
  {
    die("Error: Number of temperature windows cannot be excess 20\n\n");
  }

  $is_ca = lc $is_ca;
  if($is_ca ne "ca" && $is_ca ne "cacb")
  {
    die("Error: -m must be ca or cacb\n\n");
  }

  my $clean_pdb = clean_pdb($input_pdb);
  get_secondary_structure($clean_pdb);
  
  my @str = split(/\//, $input_pdb);
  my $pdb_name = $str[$#str];
  @str = split(/\.pdb/, $pdb_name);
  $pdb_name = $str[0];

  mkdir("${opt_str}data");

  foreach my $nscal (@nscal_array)
  {
    mkdir("${opt_str}data/nscal_$nscal");
    chdir("${opt_str}data/nscal_$nscal");
    mkdir("setup");
    my ($prefix, $prm_name);
    if($is_ca eq "ca")
    {
      ($prefix, $prm_name) = creat_CG_model($clean_pdb, 0, $nscal);
    }
    else
    {
      ($prefix, $prm_name) = creat_CG_model($clean_pdb, 1, $nscal);
    }
    #`cp ../../setup/rmsd_vs_time_template.inp ./setup`;
    `cp ../../setup/secondary_struc_defs.txt ./setup`;
    `cp ../../setup/wham_cv_template.cntrl ./setup`;
    `cp ../../setup/wham_qbb_template.cntrl ./setup`;
    #`cp ../../setup/wham_rmsd_template.cntrl ./setup`;
    #`cp ../../setup/analysis_folding_stability.pl ./setup`;

    open(TEMP, "<../../setup/mubrex_template.cntrl") || die("Error: Cannot find setup/mubrex_template.cntrl\n\n");
    open(INP, ">mubrex.cntrl") || die("Error: Cannot create mubrex.cntrl\n\n");
    my $n = $#temp_array + 1;
    while(my $line = <TEMP>)
    {
      if($line =~ /^nwindows =/)
      {
        print INP "nwindows = $n\n";
      }
      elsif($line =~ /^temps = /)
      {
        print INP "temps = $temp_set\n";
      }
      elsif($line =~ /^tpn = /)
      {
        print INP "tpn = $n\n";
      }
      elsif($line =~ /^psf = /)
      {
        print INP "psf = setup/$prefix.psf\n";
      }
      elsif($line =~ /^top = /)
      {
        print INP "top = setup/$prefix.top\n";
      }
      elsif($line =~ /^param = /)
      {
        print INP "param = setup/$prm_name\n";
      }
      elsif($line =~ /^starting_strucs/)
      {
        for(my $i = 1; $i <= $n; $i++)
        {
          if($if_unfolded eq 0 || $opt_temp eq 1)
          {
            print INP "starting_strucs_t$i = setup/$prefix.cor\n";
          }
          elsif($if_unfolded eq 1)
          {
            print INP "starting_strucs_t$i = setup/unfolded.cor\n";
          }
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
    
    if($opt_temp eq 1)
    {
      print "-> Going to optmize temperature distribution for nscal=$nscal\n";
      #`cp \$CG_MODEL_HOME/perl/Yang/opt_temp.pl ./`;
      system("qsub -N o_${pdb_name}_$nscal -l nodes=1:ppn=17 -l mem=8gb -l walltime=$walltime -f opt_temp.pl -A $allocation");
    }
    else
    {
      print "-> Going to run REX calculations for nscal=$nscal\n";
      open(PBS, ">rex.pbs") || die("Error: Cannot create rex.pbs\n\n");
      print PBS "#!/bin/csh
#PBS -r n
#PBS -m b
#PBS -m e
#PBS -o \$PBS_JOBID.out
#PBS -e \$PBS_JOBID.err
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat \$PBS_NODEFILE`
cp \$PBS_NODEFILE pbs_nodefile
set NPROCS = `wc -l < \$PBS_NODEFILE`
echo This job has allocated \$NPROCS nodes

cd \$PBS_O_WORKDIR\n";
      if($if_unfolded eq 1)
      {
        system("parse_cg_prm.py -t setup/$prefix.top -p setup/$prm_name");
        my @str = split(/\.prm/, $prm_name);
        my $xml_name = $str[0].".xml";
        print PBS "run_REX_LD.py -p setup/$prefix.psf -c setup/$prefix.cor -x setup/$xml_name -t 1000 -b 1000 -s 100000 -n $n -o setup/unfolded\n";
      }
      print PBS "parallel_temperature_REX.py -f mubrex.cntrl\n";
      close(PBS);
      system("qsub -N n_${pdb_name}_$nscal -l nodes=1:ppn=".($n+1)." -l walltime=$walltime -f rex.pbs -A $allocation");
    }
    chdir("../../");
  }
}
elsif($if_analyze > 0)
{
  my $opt_str = "";
  if($if_unfolded eq 1)
  {
    $opt_str = "unfolded_";
  }
  foreach my $nscal (@nscal_array)
  {
    chdir("${opt_str}data/nscal_$nscal");
    `cp ../../setup/secondary_struc_defs.txt ./setup`;
    `cp ../../setup/wham_cv_template.cntrl ./setup`;
    `cp ../../setup/wham_qbb_template.cntrl ./setup`;
    #`cp ../../setup/analysis_folding_stability.pl ./setup`;
    analysis_folding_stability($is_ca, $if_analyze, 20001);
    open(TEMP, "<mubrex.cntrl") || die("Error: Cannot find mubrex.cntrl\n\n");
    my $n = 0;
    while(my $line = <TEMP>)
    {
      if($line =~ /^nwindows =/)
      {
        print INP $line;
        chomp($line);
        my @str = split(/=/, $line);
        $n = $str[1];
        $n =~ s/^\s+|\s+$//g;
        last;
      }
    }
    close(TEMP);
    print "-> Going to analyze REX results at $if_analyze K for nscal=$nscal\n";
    system("qsub -N analyze_$nscal -l nodes=1:ppn=$n -l mem=8gb -l walltime=$walltime -f analysis.pbs -A $allocation");
    chdir("../../");
  }
}



########################################################################
sub clean_pdb
{
  my $pdb = $_[0];
  
  my $clean_data = "";

  print "-> Cleaning PDB file $pdb\n";

  open (PDB, "<$pdb") || die("Error: cannot open $pdb\n\n");
  while(my $line = <PDB>)
  {
  	if($line =~ /^ATOM /)
  	{
  	  my $alt_loc = substr($line, 16, 1);
      if ($alt_loc eq " " || $alt_loc eq "A" || $alt_loc eq "1")
      {
      	$clean_data .= substr($line, 0, 16) . " " . substr($line, 17);
      }
  	}

  	if($line =~ /^ENDMDL /)
  	{
  	  last;
  	}
  	if($line =~ /^TER /)
  	{
  	  last;
  	}
  }
  close(PDB);

  my @str = split(/\//, $pdb);
  my $pdb_name = $str[$#str];
  @str = split(/\./, $pdb_name);
  my $pdb_code = $str[0];

  open (PDB, ">${pdb_code}_clean.pdb") || die("Error: cannot create ${pdb_code}_clean.pdb\n\n");
  print PDB $clean_data;
  close(PDB);
  
  print "   PDB file cleaned\n";

  return "${pdb_code}_clean.pdb";
}

########################################################################
sub creat_CG_model
{
  my $pdb = $_[0];
  my $is_ca = $_[1];
  my $nscal = $_[2];
  
  my $potential_name = '';
  if($is_ca eq 0)
  {
  	print "-> Creating Ca model\n";
  	$potential_name = "bt";
  }
  else
  {
  	print "-> Creating Ca-SCM model\n";
  	$potential_name = "mj";
  }
  
  mkdir("create_model");
  chdir("create_model");

  open (IN, ">go_model.cntrl") || die("Error: cannot create Go model input file\n\n");
  print IN "pdbfile = ../$pdb
nscal = $nscal
casm = $is_ca 
potential_name = $potential_name
fnn = 1\n";
  close(IN);
  
  system("create_cg_protein_model.py -f go_model.cntrl > go_model.log 2>&1");

  my @str = split(/\./, $pdb);
  my $pdb_code = $str[0];
  my $prefix = lc $pdb_code;

  if($is_ca eq 0)
  {
    $prefix .= "_ca";
  }
  else
  {
  	$prefix .= "_ca-cb";
  }
  my $prm_name = lc $pdb_code;
  $prm_name .= "_nscal" . $nscal . "_fnn1_go_" . $potential_name . ".prm";

  if(-e "$prefix.psf")
  {
  	`cp $prefix.psf ../setup/$prefix.psf`;
  	`cp $prefix.top ../setup/$prefix.top`;
  	`cp $prefix.cor ../setup/$prefix.cor`;
  	`cp $prm_name ../setup/$prm_name`;
  	chdir("../");
  	`rm -rf create_model`;
  	print "   CG model created\n";
  }
  else
  {
  	die("Error: failed to create CG model from $pdb\n\n");
  }
  
  return ($prefix, $prm_name);
}

######################################################################
sub gen_temp
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

  my $result = join(" ", @temp_array);
  return $result;
}

######################################################################
sub get_secondary_structure
{
  my $pdb = $_[0];

  print "-> Getting secondary structure information\n";

  my $screen = `stride $pdb`;

  if($screen !~ /^\QREM  -------------------- Secondary structure summary\E/)
  {
    die($screen);
  }

  my @str = split(/\n/, $screen);
  my $sec_str = "";
  foreach my $line (@str)
  {
    if($line =~ /^STR /)
    {
      my $ss = unpack "x10 a50", $line;
      $sec_str .= $ss;
    }
  }

  my @sec = ();
  while($sec_str =~ m/(E{4,})/g)
  {
    my $len = length($1);
    my $end = pos($sec_str);
    my $start = $end - $len + 1;
    push@sec, [$start, $end];
  }

  while($sec_str =~ m/(H{4,})/g)
  {
    my $len = length($1);
    my $end = pos($sec_str);
    my $start = $end - $len + 1;
    push@sec, [$start, $end];
  }

  while($sec_str =~ m/(G{4,})/g)
  {
    my $len = length($1);
    my $end = pos($sec_str);
    my $start = $end - $len + 1;
    push@sec, [$start, $end];
  }

  while($sec_str =~ m/(I{4,})/g)
  {
    my $len = length($1);
    my $end = pos($sec_str);
    my $start = $end - $len + 1;
    push@sec, [$start, $end];
  }

  @sec = sort {$a->[0] <=> $b->[0]} @sec;

  open(SEC, ">setup/secondary_struc_defs.txt") || die("Error: Cannot create setup/secondary_struc_defs.txt\n\n");
  my $n = 1;
  foreach my $ss (@sec)
  {
    my $start = $ss->[0];
    my $end = $ss->[1];
    print SEC "$n $start $end\n";
    $n++;
  }
  close(SEC);
  print "   Done.\n";
}

######################################################################
sub analysis_folding_stability
{
  my $is_ca = $_[0];
  my $target_temp = $_[1];
  my $start = $_[2];

  my $info_file = "";

  my $info_list = `ls ./info*.log`;
  chomp($info_list);
  my @info_array = split(/\n/, $info_list);
  if(@info_array eq 1)
  {
    $info_file = $info_array[0];
  }
  elsif(@info_array > 1)
  {
    my $restart_num = 0;
    foreach my $in (@info_array)
    {
      if($in eq "info.log")
      {
        next;
      }
      my @str = split(/\.log/, $in);
      my $num = $str[0];
      @str = split(/\Qinfo_r_\E/, $num);
      $num = $str[1];
      if($num > $restart_num)
      {
        $restart_num = $num;
      }
    }
    $info_file = "info_r_${restart_num}.log";
  }

  open(SCR, ">analysis.pbs") || die("Error: Cannot create analysis.pbs\n\n");
  print SCR "#PBS -r n
#PBS -m b
#PBS -m e
#PBS -o \$PBS_JOBID.out
#PBS -e \$PBS_JOBID.err
cd \$PBS_O_WORKDIR
analysis_folding_stability.pl -i $info_file -m $is_ca -t $target_temp -s $start -e \"mix\" > analysis.out 2>&1\n";
  close(SCR);
}
