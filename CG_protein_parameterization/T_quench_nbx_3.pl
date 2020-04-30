#!/usr/bin/perl

use Getopt::Long;

my ($help, $input_pdb, $is_ca, $nscal, $num_traj, $Q_k, $tpn, $allocation, $walltime);
GetOptions(
 'help|h!' => \$help,
 'input|i=s' => \$input_pdb,
 'model|m=s' => \$is_ca,
 'nscal|n=s' => \$nscal,
 'ntraj|t=s' => \$num_traj,
 'qthreshold|q=s' => \$Q_k,
 'totalprocnum|p=s' => \$tpn,
 'allocation|a=s' => \$allocation,
 'walltime|w=s' => \$walltime,
);

my $usage = "
  Usage: perl T_quench.pl 
              --input | -i <INPUT.PDB> for CG model creation
              --model | -m <ca or cacb> ca for CA cg model;
                           cacb for CA-CB cg model
              --nscal | -n <NUM> for nscal
              --ntraj | -t <NUM> for number of trajectories
              --qthreshold | -q <NUM> for determining folded state
              [--totalprocnum | -p] <NUM> for total cpu number required.
              [--allocation | -a] <\"ACCOUNT\"> for ACI cluster. 
                                  Default is \"epo2_a_g_bc_default\"
              [--walltime | -w] <\"TIME\"> for ACI cluster. 
                                Default is \"10:00:00:00\"
              [--help | -h]\n\n";

my $ppn = 1;

if($help)
{
  die($usage);
}
elsif(!(defined($input_pdb) && defined($is_ca) && defined($nscal) && defined($num_traj) && defined($Q_k)))
{
  die($usage);
}

if(!defined($allocation))
{
  $allocation = "cyberlamp";
}
if(!defined($walltime))
{
  $walltime = "10:00:00:00";
}
if(!defined($tpn))
{
  $tpn = 20;
}

$is_ca = lc $is_ca;
if($is_ca ne "ca" && $is_ca ne "cacb")
{
  die("Error: -m must be ca or cacb\n\n");
}

my $clean_pdb = clean_pdb($input_pdb);

mkdir("setup");

get_secondary_structure($clean_pdb);

my ($prefix, $prm_name);
if($is_ca eq "ca")
{
  ($prefix, $prm_name) = creat_CG_model($clean_pdb, 0, 0, 1, 0, 0, "bt", 1, $nscal);
}
else
{
  ($prefix, $prm_name) = creat_CG_model($clean_pdb, 1, 0, 0, 1, 1, "mj", 1, $nscal);
}

my $psf = "setup/$prefix.psf";
my $strt = "setup/$prefix.cor";
my $top = "setup/$prefix.top";
my $prm = "setup/$prm_name";

open(IN, ">setup/T_quench.cntrl");
print IN "tpn = $tpn
ppn = $ppn
num_traj = $num_traj
psf = $psf
top = $top
param = $prm
starting_strucs = $strt
temp_equil = 1000
temp_prod = 310
Q_threshold = $Q_k
secondary_structure_def = setup/secondary_struc_defs.txt
log_file = info.log
restart = 1\n";
close(IN);

open(PBS, ">Tq.pbs") || die("Error: Cannot create Tq.pbs\n\n");
print PBS "#PBS -r n
#PBS -m b
#PBS -m e
#PBS -o \$PBS_JOBID.out
#PBS -e \$PBS_JOBID.err
cd \$PBS_O_WORKDIR
temperature_quenching.py -f setup/T_quench.cntrl\n";
close(PBS);
  
print "-> Going to run Tq using $tpn CPUs\n";
system("qsub -N Tq_$input_pdb -l nodes=1:ppn=$tpn -l walltime=$walltime -f Tq.pbs -A $allocation");

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
  my $bondlength_go = $_[2];
  my $angle_dw = $_[3];
  my $dihedral_go = $_[4];
  my $improper_go = $_[5];
  my $native_energy = $_[6];
  my $charges = $_[7];
  my $nscal = $_[8];

  if($is_ca eq 0)
  {
  	print "-> Creating Ca model ";
  }
  else
  {
  	print "-> Creating Ca-SCM model ";
  }
  
  my $combined_name = "";
  if($bondlength_go eq 0)
  {
  	$combined_name .= "b_ng_";
  }
  else
  {
  	$combined_name .= "b_go_";
  }

  if($angle_dw eq 0)
  {
  	$combined_name .= "a_go_";
  }
  else
  {
  	$combined_name .= "a_dw_";
  }

  if($dihedral_go eq 0)
  {
  	$combined_name .= "d_ng_";
  }
  else
  {
  	$combined_name .= "d_go_";
  }

  if($improper_go eq 0)
  {
  	$combined_name .= "i_ng_";
  }
  else
  {
  	$combined_name .= "i_go_";
  }

  if($native_energy eq "mj")
  {
  	$combined_name .= "n_mj_";
  }
  elsif($native_energy eq "bt")
  {
  	$combined_name .= "n_bt_";
  }

  if($charges eq 0)
  {
  	$combined_name .= "c_0_";
  }
  else
  {
  	$combined_name .= "c_1_";
  }

  $combined_name .= "s_$nscal";

  print "$combined_name\n";
  
  mkdir("create_model");
  chdir("create_model");

  open (IN, ">go_model.cntrl") || die("Error: cannot create Go model input file\n\n");
  print IN "charmm = \$c35b5_dhdwp
pdb=../$pdb
nscal=$nscal
pot=$native_energy
bondlength_go = $bondlength_go
dihedral_go = $dihedral_go
improperdihed_go = $improper_go
casm = $is_ca
charges = $charges
angle_dw = $angle_dw
fnn = 1\n";
  close(IN);
  
  system("create_cg_protein_model_v34_0.37_nbx3.pl go_model.cntrl > go_model.log 2>&1");

  if(-e "create_psf.inp")
  {
  	system("\$c35b5_dhdwp < create_psf.inp > create_psf.log 2>&1");
  }
  else
  {
  	die("Error: failed to create CG model from $pdb\n\n");
  }

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
  $prm_name .= "_nscal" . $nscal . "_fnn1_go_" . $native_energy . ".prm";

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
