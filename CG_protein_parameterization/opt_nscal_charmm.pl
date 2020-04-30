#!/usr/bin/perl

use threads;
use Getopt::Long;
use Data::Dumper;

my ($help, $input_pdb, $dom_def, $restart);
GetOptions(
 'help|h!' => \$help,
 'input|i=s' => \$input_pdb,
 'domain|d=s' => \$dom_def,
 'restart|r=s' => \$restart,
);

my $usage = "
  Usage: perl opt_nscal.pl 
              --input | -i <INPUT.PDB> for CG model creation
              --domain | -d <DOMAIN.DAT> for domain defination
              [--restart | -r] <0 or 1> restart optimization. Default 0, not restart. 
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
elsif(!(defined($input_pdb) && defined($dom_def)))
{
  die($usage);
}

if(!defined($restart))
{
  $restart = 0;
}

my $is_ca = "ca";

my %nscal_set = ("a" => [1.17, 1.427, 2.012, 2.48, 1.17],
                 "b" => [1.442, 1.759, 2.48, 2.48, 1.17],
                 "c" => [1.114, 1.359, 1.916, 2.48, 1.17],
                 "i" => [1.235, 1.507, 2.124, 2.48, 1.507]);
my $sim_step = 66666667; #for 1000 ns
my $num_traj = 10;
my $Q_threshold = 0.6908;
my $frame_threshold = 0.98;

mkdir("setup");

my @domain = parse_domain($dom_def);
my @nscal_index = ();
my $start_n = 1;
if($restart eq 0)
{
  for(my $i = 0; $i < @domain; $i++)
  {
    $nscal_index[$i] = 0;
  }
  open(INFO, ">opt_nscal.log")||die("Error: cannot create opt_nscal.log\n\n");
  print INFO "## Start at ".localtime()."\n";
}
else
{
  open(INFO, "<opt_nscal.log")||die("Error: cannot find opt_nscal.log\n\n");
  while(my $line = <INFO>)
  {
    if($line =~ /^## Round /)
    {
      chomp($line);
      my @str = split(/:/, $line);
      @str = split(/\s+/, $str[0]);
      $start_n = $str[$#str];
    }
    elsif($line =~ /^## Final nscal:/)
    {
      die("Finished. No need to restart.\n");
    }
    elsif($line =~ /^\-\> Set nscal as:/)
    {
      my $nn = 0;
      my @nscal_index_1 = ();
      while($line = <INFO>)
      {
        if($line =~ /^\-\>/)
        {
          last;
        }
        chomp($line);
        my @str = split(/ = /, $line);
        my $ns = $str[1];
        print "$ns\n";
        my @array = @{$nscal_set{$domain[$nn]->{"class"}}};
        for(my $i = 0; $i < @array; $i++)
        {
          if($array[$i] eq $ns)
          {
            push@nscal_index_1, $i;
            last;
          }
        }
        $nn += 1;  
      }
      @nscal_index = @nscal_index_1;
    }
  }
  close(INFO);
  open(INFO, "<opt_nscal.log")||die("Error: cannot find opt_nscal.log\n\n");
  open(INFO_2, ">opt_nscal_new.log")||die("Error: cannot create opt_nscal_new.log\n\n");
  while(my $line = <INFO>)
  {
    if($line =~ /^## Round \E$start_n\Q/)
    {
      last;
    }
    else
    {
      print INFO_2 $line;
    }
  }
  close(INFO);
  close(INFO_2);
  system("rm -f opt_nscal.log");
  system("mv opt_nscal_new.log opt_nscal.log");
  open(INFO, ">>opt_nscal.log")||die("Error: cannot create opt_nscal.log\n\n");
}
my $clean_pdb = clean_pdb($input_pdb);
get_secondary_structure($clean_pdb);

print "@nscal_index \n";
for(my $interation = $start_n; $interation <= @{$nscal_set{"a"}}; $interation++)
{
  mkdir("round_$interation");
  chdir("round_$interation");
  mkdir("setup");
  
  print INFO "## Round $interation:\n";
  print INFO "-> Set nscal as:\n";

  for(my $i = 0; $i < @domain; $i++)
  {
    my $dom = $domain[$i];
    $dom->{"nscal"} = $nscal_set{$dom->{"class"}}->[$nscal_index[$i]];
    if($dom->{"class"} eq "i")
    {
      print INFO "   Interface ".$dom->{"range"}->[0]."|".$dom->{"range"}->[1].": nscal = ".$dom->{"nscal"}."\n";
    }
    else
    {
      print INFO "   Domain ".($i+1).": nscal = ".$dom->{"nscal"}."\n";
    }
  }

  my ($prefix, $prm_name) = creat_CG_model("../$clean_pdb", 0, 0, 1, 0, 0, "bt", 1, \@domain);

  my $psf = "setup/$prefix.psf";
  my $strt = "setup/$prefix.cor";
  my $top = "setup/$prefix.top";
  my $prm = "setup/$prm_name";

  open(IN, ">setup/simulation.inp");
  print IN "*
!======== BEGIN: Set parameters ===========
prnlev 3
bomlev 0

!! Set up temperature run
set timestp 0.015
set fbsolu 0.05
set ftemp 310 

 rand clcg

 ! read parameter and topology files
 open unit 10 read form name \@top
 read rtf unit 10 card
 close unit 10

 open unit 10 read form name \@prm
 read param unit 10 card flex
 close unit 10

 ! Get psf
 read psf card name \@psf

 ! Get the coordinates
 open read unit 3 card name \@strt
 read coor unit 3 card
 close unit 3

 SHAKE BOND PARA

 ! Set friction forces
 scalar fbeta set \@fbsolu sele segid A end
 SCAL FBETA SHOW sele segid A end

 ! speeds up the non-bonded list generation
 nbond bygr
 update DEBY SCRE 10.0
 eten on 

 open unit 81 write file name traj_\@istep.dcd
 dyna leap langevin start nstep $sim_step timestep \@timestp -
    iprfrq 0 ieqfrq 0 ntrfrq 0  -
    iunrea -1 iunwri 31 iuncrd 81 iunvel -1 kunit -1 -
    ihtfrq 0 teminc 0 nprint 50000 nsavc 5000 nsavv 0 ihbfrq 0 -
    inbfrq -1 imgfrq 0 -
    ilbfrq 0 - ! and step frequency for writing restart file
    rbuffer 0.0 tbath \@ftemp tstruc \@ftemp -
    firstt \@ftemp finalt \@ftemp -
    iasors 0 iasvel 1 iscvel 0 ichecw 0 twindh 0.0 twindl 0.0 iseed \@rand -
    echeck 5000

  ! write crd file
  open write unit 10 card name traj_\@istep.cor
  write coor card unit 10
  * temp = \@ftemp and potential energy = ?ener
  *
  close unit 10

STOP\n";
  close(IN);
  
  my @threads = ();
  for(my $i = 1; $i <= $num_traj; $i++)
  {
    my $rand_num = int(rand()*10000000000);
    my $t = threads->create('run_exe', "\$c35b5_dhdwp prm=$prm top=$top psf=$psf strt=$strt rand=$rand_num istep=$i < setup/simulation.inp > $i.out 2>&1");
    push@threads, $t;
    $t->detach();
  }

  print INFO "-> Running simulations...\n";

  my $start_time = time();
  while()
  {
    sleep(30);
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

    open(LOG, ">simutaion.log") || die("Error: Cannot create simulation.log\n\n");
    printf LOG ("%-20s %-20s %-20s %-20s %-20s\n", "#Num", "Step", "Speed(Step/s)", "Used_Time", "Rest_Time");
    my $used_time = time() - $start_time;
    for(my $i = 1; $i <= $num_traj; $i++)
    {
      my $stage = `grep "DYNA>" $i.out | grep " 0.00000" |wc -l`;
      chomp($stage);
      my $info = `grep "DYNA>" $i.out | tail -n 1`;
      chomp($info);
      my @str = split(/\s+/, $info);
      my $step = $str[1];

      if($stage eq 0)
      {
        printf LOG ("%-20s %-20s %-20s %-20s %-20s\n", $i, "NaN", "NaN", convert_time($used_time), "NaN");
        next;
      }

      if($step eq 0)
      {
        printf LOG ("%-20s %-20s %-20s %-20s %-20s\n", $i, $step, "NaN", convert_time($used_time), "NaN");
        next;
      }

      my $speed = $step / $used_time;
      my $rest_time = ($sim_step - $step) / $speed;
      printf LOG ("%-20s %-20s %-20s %-20s %-20s\n", $i, $step, int($speed), convert_time($used_time), convert_time($rest_time));
    }
    close(LOG);
  }
  
  print INFO "-> Probability of domain stability:\n";

  my @if_stable = ();
  my @stable_frac = ();
  for(my $j = 0; $j < @domain; $j++)
  {
    $if_stable[$j] = [];
    $stable_frac[$j] = [];
  }
  
  my @threads = ();
  for(my $i = 1; $i <= $num_traj; $i++)
  {
    my $t = threads->create('run_exe', "calc_native_contact_fraction.pl -i $strt -d ../$dom_def -s ../setup/secondary_struc_defs.txt -t traj_$i.dcd");
    push@threads, $t;
    $t->detach();
  }
  
  while()
  {
    sleep(30);
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
  
  for(my $i = 1; $i <= $num_traj; $i++)
  {
    open(DAT, "<qbb_traj_$i.dat")||die("Error: cannot find qbb_traj_$i.dat\n\n");
    my $line = <DAT>;
    my $n_tot = 0;
    my @n_domain = ();
    foreach my $dom (@domain)
    {
      push@n_domain, 0;
    }
    while($line = <DAT>)
    {
      chomp($line);
      $line =~ s/^\s+|\s+$//g;
      my @qbb = split(/\s+/, $line);
      $n_tot++;
      for(my $j = 0; $j < @qbb-1; $j++)
      {
        if($qbb[$j] > $Q_threshold)
        {
          $n_domain[$j]++;
        }
      }
    }
    close(DAT);

    for(my $j = 0; $j < @domain; $j++)
    {
      my $frac = $n_domain[$j] / $n_tot;
      $stable_frac[$j]->[$i-1] = $frac;
      if($frac >= $frame_threshold)
      {
        $if_stable[$j]->[$i-1] = 1;
      }
      else
      {
        $if_stable[$j]->[$i-1] = 0;
      }
    }
  }
  
  my $break_tag = 1;
  for(my $j = 0; $j < @domain; $j++)
  {
    my $dom = $domain[$j];
    if($dom->{"class"} eq "i")
    {
      print INFO "   Interface ".$dom->{"range"}->[0]."|".$dom->{"range"}->[1].": ";
    }
    else
    {
      print INFO "   Domain ".($j+1).": ";
    }
    my $tag = 1;
    for(my $i = 0; $i < $num_traj; $i++)
    {
      $tag *= $if_stable[$j]->[$i];

      printf INFO ("%.3f ", $stable_frac[$j]->[$i]);
    }
    if($tag eq 0)
    {
      print INFO "instable\n";
      $nscal_index[$j]++;
      $break_tag = 0;
    }
    else
    {
      print INFO "stable\n";
    }
  }

  if($break_tag)
  {
    last;
  }

  chdir("../");
}

print INFO "## Final nscal:\n";
for(my $i = 0; $i < @domain; $i++)
{
  my $dom = $domain[$i];
  if($dom->{"class"} eq "i")
  {
    print INFO "   Interface ".$dom->{"range"}->[0]."|".$dom->{"range"}->[1].": nscal = ".$dom->{"nscal"}."\n";
  }
  else
  {
    print INFO "   Domain ".($i+1).": nscal = ".$dom->{"nscal"}."\n";
  }
}
close(INFO);

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

########################################################################
sub clean_pdb
{
  my $pdb = $_[0];
  
  my $clean_data = "";

  my @AA = ("GLY","ALA","VAL","LEU","ILE","MET","PHE","PRO","SER","THR","CYS","ASN","GLN","TYR","TRP","ASP","GLU","HIS","LYS","ARG");

  print "-> Cleaning PDB file $pdb\n";

  open (PDB, "<$pdb") || die("Error: cannot open $pdb\n\n");
  while(my $line = <PDB>)
  {
    if($line =~ /^ATOM /)
    {
      my $alt_loc = substr($line, 16, 1);
      my $resname = substr($line, 17, 3);
      if (($alt_loc eq " " || $alt_loc eq "A" || $alt_loc eq "1") && $resname ~~ @AA)
      {
        $clean_data .= substr($line, 0, 16) . " " . substr($line, 17);
      }
    }

    if($line =~ /^ENDMDL /)
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
  my @domain = @{$_[8]};

  my $nscal = 1;

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

  open (IN, ">domain_def.dat") || die("Error: cannot create domain_def.dat\n\n");
  my $ndomain = 0;
  my $n = 1;
  my @index_1 = ();
  my @n_array = ();
  foreach my $d (@domain)
  {
    if($d->{"class"} eq "i")
    {
      $n_array[$d->{"range"}->[0]]->[$d->{"range"}->[1]] = $d->{"nscal"};
    }
    else
    {
      foreach my $a (@{$d->{"range"}})
      {
        print IN "domain = ".$a->[0]."-".$a->[1]."\n";
        push@index_1, $n;
        $ndomain++;
        print IN "scale factor = ".$d->{"nscal"}."\n";
        $n_array[$n] = [];
        $n_array[$n]->[$n] = $d->{"nscal"};
      }
      $n++;
    }
  }

  for(my $i = 1; $i <= $ndomain-1; $i++)
  {
    for(my $j = $i + 1; $j <= $ndomain; $j++)
    {
      print IN "scale factor = ".$n_array[$index_1[$i-1]]->[$index_1[$j-1]]."\n";
    }
  }
  close(IN);
  
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
fnn = 1
ndomain = $ndomain
int_file = domain_def.dat\n";
  close(IN);
  
  system("create_cg_protein_model_v34_nbx3_multidomain.pl go_model.cntrl > go_model.log 2>&1");

  if(-e "create_psf.inp")
  {
    system("\$c35b5_dhdwp < create_psf.inp > create_psf.log 2>&1");
  }
  else
  {
    die("Error: failed to create CG model from $pdb\n\n");
  }

  my @str = split(/\.pdb/, $pdb);
  my $pdb_code = $str[0];
  @str = split(/\//, $pdb_code);
  $pdb_code = $str[$#str];
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

###########################################################
sub run_exe
{
  my $cmd = $_[0];
  system($cmd);
}

###########################################################
sub convert_time
{
  my $time = $_[0];
  my $hours = int($time / 3600);
  my $minutes = int(($time - $hours * 3600) / 60);
  my $seconds = int($time - $hours * 3600 - $minutes * 60);
  return "$hours:$minutes:$seconds";
}

###########################################################
sub analysis_traj
{
  my $num_traj = $_[0];
  my $Q_threshold = $_[1];
  my $frame_threshold = $_[2];


}
