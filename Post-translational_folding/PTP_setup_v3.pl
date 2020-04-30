#!/usr/bin/perl

use POSIX;
use Getopt::Long;
use Data::Dumper;

my ($help, $CSP_dir, $pdb_name, $mutant_type, $ppn, $exp_time, $temp, $Q_threshold, $allocation, $rep_per_traj, $start_rep_idx);
GetOptions(
 'help|h!' => \$help,
 'cspdir|d=s' => \$CSP_dir,
 'pdbname|p=s' => \$pdb_name,
 'type|m=s' => \$mutant_type,
 'temp|t=s' => \$temp,
 'exp_time|e=s' => \$exp_time,
 'ppn|n=s' => \$ppn,
 'Q_th|q=s' => \$Q_threshold,
 'allocation|A=s' => \$allocation,
 'rep_per_traj|j=s' => \$rep_per_traj,
 'start_rep_idx|s=s' => \$start_rep_idx,
);

my $usage = "
  Usage: perl PTP_setup_v2.pl 
              --cspdir | -d <DIR> for continuous synthesis 
              --pdbname | -p <XXX> used in CG model creation
              --temp | -t <temperature> used for MD simulation
              --exp_time | e <seconds> for simulation
              --Q_th | -q <Threshold of Q> to determine folding status
              [--type | -m <mutant type>] Default 'wild'
              [--allocation | -A <Allocation> for job to run] Default 'cyberlamp'
              [--ppn | -n <number of CPUs> for each trajectory] Default 1.
              [--rep_per_traj | -j <Num of replicas for each traj>] Default 1.
              [--start_rep_idx | -s <Start replica index for each traj>] Default 1.
              [--help | -h]\n\n";

if($help)
{
  die($usage);
}
elsif(!(defined($CSP_dir) && defined($pdb_name) && defined($exp_time) && defined($temp) && defined($Q_threshold)))
{
  die($usage);
}

if(!defined($ppn))
{
  $ppn = 1;
}
if(!defined($mutant_type))
{
  $mutant_type = 'wild';
}
if(!defined($allocation))
{
  $allocation = 'cyberlamp';
}
if(!defined($rep_per_traj))
{
  $rep_per_traj = 1;
}
if(!defined($start_rep_idx))
{
  $start_rep_idx = 1;
}

my $alpha = 4331293;
my $dt = 0.015/1000;
my $nsave = 5000;
$dt = $dt*$nsave*$alpha/1e9;
my $sim_step = POSIX::ceil($exp_time / $dt) * $nsave;

my $walltime = "14:00:00:00";
if($allocation != 'cyberlamp')
{
  $walltime = "30:00:00:00";
}

my $job_name = "P$pdb_name";

my $prefix = "${pdb_name}_model_clean";
my $screen = `ls $CSP_dir/setup/$prefix* 2>&1`;
if($screen =~ /No such file or directory/)
{
  $prefix = "${pdb_name}_m_clean";
  $screen = `ls $CSP_dir/setup/$prefix* 2>&1`;
  if($screen =~ /No such file or directory/)
  {
    die("Cannot find $CSP_dir/setup/$prefix*");
  }
}

mkdir("setup");
system("cp $CSP_dir/setup/$prefix* ./setup/");
system("cp $CSP_dir/setup/domain_def.dat ./setup/");
system("cp $CSP_dir/setup/secondary_struc_defs.txt ./setup/");

my $psf = "../setup/${prefix}_ca.psf";
my $top = "setup/${prefix}_ca.top";
my $prm = "setup/${prefix}_nscal1_fnn1_go_bt.prm";
system("parse_cg_prm.py -t $top -p $prm");
my $xml = "../setup/${prefix}_nscal1_fnn1_go_bt.xml";
my $native_cor = "../setup/${prefix}_ca.cor";
my $sec_struc_def = "../setup/secondary_struc_defs.txt";

my $n_traj = `ls $CSP_dir/traj/ | wc -l`;

if($mutant_type =~ /^silent/)
{
  $job_name .= "m";
}
else
{
  $job_name .= substr($mutant_type, 0, 1);
}

for(my $i = 1; $i <= $n_traj; $i++)
{
  mkdir($i);
  chdir($i);
  system("cp  ../$CSP_dir/traj/$i/prot_* ./");
  my $cor = `ls | grep .cor`;
  chomp($cor);
  my $vec = `ls | grep .vel`;
  chomp($vec);
  for(my $j = $start_rep_idx; $j <= $rep_per_traj; $j++)
  {
    # check output
    my $tag_finish = 0;
    if(-e "traj_${i}_${j}.out")
    {
      open(OUTPUT, "<traj_${i}_${j}.out");
      my @lines = <OUTPUT>;
      my $last_line = $lines[$#lines];
      chomp($last_line);
      $last_line =~ s/^\s+|\s+$//g;
      close(OUTPUT);
      if($last_line !~ /^Done/)
      {
        my @words = split(/\s+/, $last_line);
        my $current_step = $words[1];
        #print "$current_step $sim_step\n";
        if($current_step >= $sim_step)
        {
          $tag_finish = 1;
        }
      }
      else
      {
        $tag_finish = 1;
      }
    }
    
    if($tag_finish eq 1)
    {
      print("traj ${i}_${j} has already finished\n");
      next;
    }
    #
    my $rand_num = int(rand()*10000000);
    open(JOB, ">job_${j}.pbs");
    print JOB "#PBS -N ${job_name}_${i}_${j}
#PBS -r n
#PBS -o \$PBS_JOBID.out
#PBS -e \$PBS_JOBID.err
#PBS -l nodes=1:ppn=$ppn 
#PBS -l walltime=$walltime
#PBS -A $allocation

cd \$PBS_O_WORKDIR
post_trans_single_run_v2.py $psf $cor $xml $temp $ppn traj_${i}_${j} $rand_num $vec $sim_step $sec_struc_def $Q_threshold $native_cor
";
    close(JOB);
    print("submit traj ${i}_${j}\n");
    system("qsub job_${j}.pbs");
  }
  
  chdir("../");
}
