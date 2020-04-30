#!/usr/bin/perl
use File::Basename;
use Math::Trig;
use Cwd;

#
$version = 1.34;
#
# Take in a PDB file, create Calpha-Sidechain GO-type
# pdb, psf, prm, and top files, with
# Miyazawa-Jernigan parameters for non-bonds.
#
# Purpose: Create a Ca-Cb model to be used in
#   CHARMM.
#
# Assumptions: (1) cutoff to define native is
#    dependent on residue pair types.
#    (2) MMTSB package is installed
#    (3) Center-of-mass of sidechain computed
#    from heavy atoms only
#
# Inputs: (1) regular PDB file
#
# Author:
# Edward P. O'Brien Jr., D. Thirumalai, B.R. Brooks groups, NIH/UMD
# 2/7/07
#
# Please cite the following article when referencing this model:
# Title: Effects of denaturants and osmolytes on proteins are accurately predicted by the molecular transfer model
# Author(s): O'Brien EP, Ziv G, Haran G, et al.
# Source: PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES OF THE UNITED STATES OF AMERICA   Volume: 105   Issue: 36   Pages: 13403-13408
#
# Modifications:
#   (3/7) I now define 'native' side-chain contacts
#     as two sc's having heavy atoms within 4.5 Angstroms of
#     each other.
#   (3/8) Compute dihedral prefactors independent of <s.c. interaction energy>,
#     instead have the energy difference between wells be equal to ~kbT at
#     350 K, ~0.7 kcal/mol
#   (3/15) The vdw radius of the side-chains no longer scales with
#     <sc interaction energy>. Instead we can now independently scale
#     it with 'fnn' arguement.
#   (3/15) Hbond energies are set to 1, instead of the average
#     energy of the side-chain interaction matrix.
#   (3/26) Introduce attracive backbone-sidechain interactions
#   (3/26) Weaken side-chain interactions and strengthen backbone
#     interactions.
#   (3/26) Modified definition of hbond energies. Now, if one residue
#     takes part in two hydrogen bonds with another residue (as in the
#     case of resiudes in anti-parallel beta-sheets), I double the
#     strength of the Hydrogen bond for those residues.
#   (3/30) Made angle and bond force constants independent of non-bonded 
#     potential.
#   (4/11) Increased angle force constant because they were to flexible.
#     This caused huge forces when computing the dihedral forces, because
#     some angles would go to ~180 degress! 
#   (2/1/08) Changed dihedral and hydrogen bond energy so that a helix
#     is not stable, by itself, at 300 K. CHANGED CUTNB variable from
#     399 to 24 angstroms!!!
#   (2/26/08) Changed atom nomenclature from CA -> A and CB -> B to
#     allow larger proteins to be used. Currently charmm reads in
#     only 4 character fields per atom, so CA999 will be read in
#     as CA99.
#   (3/6/08) Only use the flexible parameter reader in charmm. This
#     requires that the parameter file now contain the masses and
#     codes of each atom (previously in TOPOLOGY file).
#   (3/17/08) Allow multiple proteins, non-covalently linked, to
#     be created.
#   (3/24/08) Allow interfaces to be defined that have a different scaling
#     for non-bonded interactions.
#   (2/21/09) Changed dihedral potential to {0.55,0.275} and {0.3,0.15} for helices.
#   Set H-bonds to 0.5, and 0.25 for helices. I now use STRIDE to identify secondary
#   structure elements and with this differentiate helices from everything else.
#   (3/7/09) Allow user to use Non-Go (generic) dihedral angles by changing logical flag
#   'dihedral_go' to 0. I use the CL Brooks dihedral poential derived from the PDB.
#   (5/10/2010) Allow user to define a homogeneous side-chain side-chain interaction energy
#   using the 'CT' command line keyword for the 'Potential' argument. In this case the 'scale factor'
#   argument becomes the sc-sc energy.
#   (5/10/2010) Allow user to create a C-alpha model only (i.e. without side chain interaction site)
#   by modifying the hard-coded parameter '$casm=0'.
#   (5/13/2010) Allow the user to define a sequence, with out a PDB file, for which this program
#   will create a generic set of parameter and topology files using the 'generic' keyword on the command line.
#   (5/14/2010) Changed dihedral potential to {0.3,0.15} for all secondary structures WHEN 'casm=1'.
#   (10/26/2010) Major changes: (1) Read in a control file to set all options or use command line. This is necessary
#   as the number of options is very large now. (2) Allow charges to be added to protein with  
#   keyword 'charges = 1'. (3) Allow the double well backbone angle potential of Best and Hummer
#   to be used with 'angle_dw = 1'. (4) Allow non-go Calpha-Calpha bond length to be used (3.81 angstroms)
#   with keyword 'bondlength_go = 1'. (5) Create PSF charmm script and file automatically if 'charmm' path defined.
#   (10/27/2010) The logic in writing out NBFIX values was messed up when the c-alpha only model
#   was invoked. I wrote out energies for the same residue pairs twice some times. Now, when casm=1,
#   add the energies in an NxN matrix and write out non-zero values in NBFIX section; where N is
#   the number of residues in the protein.
#   (10/29/2010) When using the C_alpha only model use the Karanacolis-Brooks definition for
#   the collision diameter of non-native interactions.
#   (11/15/2010) Allow the user to define multiple PDB files that contain both unstructured protein
#   fragments and structured protein sections that this code stiches together. There are only
#   two main assumptions: (1) that the order of the PDB files starts from the N-terminus of
#   of the final stiched together protein, and (2) that you define the non-bonded potential type
#   for each PDB. NOTE WELL, that for the unstructured PDB's you only need the C-alpha positions
#   listed in the PDB. You can automatically create the input PDB file for a given amino acid sequence
#   using the script "create_unstructured_pdb.pl".
#   (11/16/2010) When the 'GENERIC' option is invoked use the full non-bonded hamiltonian
#   in NBFIX - that is allow unstructured portions of the protein to have attractive interactions
#   with as strength of 0.3*nscal. You must no specify in the control file the generic
#   non-bonded potential as either 'generic-bt', 'generic-mj', or 'generic-kgs'.
#   (11/21/2010) Allow the user to request a different Atom name prefix in the cntrl file using
#   'ca_name' or 'sc_name' for C-alpha or side chain interaction sites, respectively.
#   (5/2/2011) Allow the user to define multiple 'ca_name's in the cntrl file so that
#   they can name the residues with in different PDB files with different names. Allow
#   the user to define whether interaction parameters between two PDB structures
#   are printed out to the parameter file by defining 'interact_scale' keyword in
#   the control file. Without this keyword these parameters will not be written out
#   even if there are two PDB files listed in the cntrl file.
#   (5/19/2014) We only calculate the side chain center of mass if we are using the calpha side chain model
#   Ed, will this screw up quality control if there are missing side chain atoms that we use 
#   to identify native contacts?
#   (5/20/2014) Have this code automatically determine the format of the PDB file to accurately extract
#   the PDB coordinates. We do this by checking if the first time we read in the x,y,z coordinates are
#   all floating points.
#   (10/6/2014) Modify Determine_CM subroutine to avoid bug with the non-standard arrangment of 
#   side chain and backbone atoms in Proline.
#   (10/8/2014) Allow the user to define mixed potentials containing 'go' and 'non-go' segments.
#   Also define a keyword 'improperdihed_go' to allow go and non-go potentials to be used
#   for the side chain atoms about the C_alpha position in the Ca-SCM model.  This new functionality
#   will allow users to have intrinsically unstructured Ca-SC segments attached to Go-based segments.
#
=pod
Example 1: The control file for creating a Calpha-SCM with non-go dihedrals
and charged side chains.
##############################################################################
# Control file for creating go model
#

##################################################################################
=cut

## BEGIN: read command line ##
if($#ARGV+1 < 1)
{
  print "USAGE: perl go_model.cntrl (Optional command line arguments)\n";
  die "ERROR: 1 command line arguement is required\n";
}
$cntrl = $ARGV[0];
### END: read command line ###

## BEGIN: defaults ##
$prefix = "none"; # output prefix
$charmm = "none"; # path to charmm executable
@pdb = ();
$offset = 0;
$nscal = 1.0;
$fnn = 1.0;
@pot = ();
@ndomain=(0);
$bondlength_go[0] = 1; # 1 = use the Calpha-Calpha bond lenth in PDB, 0 = set it to 3.81 angstroms
$dihedral_go[0] = 1; # 1 = use dihedrals that are GO, 0 = dihedrals are NON-GO!
$casm = 1; # 1 = create calpha-scm, 0 = create calpha model only!
$charges[0] = 0; # 1 = add charges, 0 = don't add chages
$improperdihed_go[0] = 1; # 1 = use improper dehdral angles that are GO, 0 = impropers are Non-Go!
$angle_dw[0] = 0; # 1 = use double well angle potential, 0 = don't use it.
$heav_cut = 4.5; # angstroms, Definition of cutoff for sidechain heavy atoms contact
$int_scfac = $nscal; # default
$int_file = "none";
$interact_scale = "undefined";
$mass_offset = 0; # Start labeling atom-id's as $mass_offset+1
@cap = ("A"); # C-alpha atom prefix in top, prm, psf files
@scp = ("B"); # side chain prefix
### END: defaults ###

## BEGIN: Parse the control file for information##
open(IN,"$cntrl") or die "ERROR: file $cntrl does not exist\n";
while(<IN>)
{
  if($_ !~ m/#/ and $_ =~ m/\w+/) # if the input line does not contain a pound sign, and has more than just white spaces
  {
    ($_) = cleanws($_);
    @dat = split(/=/,$_);
    ($dat[0]) = cleanws(lc($dat[0]));
    ($dat[1]) = cleanws($dat[1]);
    @dat2=split(/\s+/,$dat[1]);
    if($dat[0] =~ m/prefix/) {$prefix = $dat[1]}
    elsif($dat[0] =~ m/pdb/)
    {
      for($i=0;$i<=$#dat2;$i++) {push(@pdb,$dat2[$i])} # read in the list of PDB file names here.
    }
    elsif($dat[0] =~ m/charmm/) {$charmm = $dat[1]}
    elsif($dat[0] =~ m/nscal/) {$nscal = $dat[1]}
    elsif($dat[0] =~ m/fnn/) {$fnn = $dat[1]}
    elsif($dat[0] =~ m/pot/) 
    {
      for($i=0;$i<=$#dat2;$i++) {push(@pot,$dat2[$i])} # read in a list of potentials, it must match the number of PDBs
    }
    elsif($dat[0] =~ m/ndomain/) 
    {
      @ndomain = ();
      for($i=0;$i<=$#dat2;$i++) {$ndomain[$i]=$dat2[$i]}
    }
    elsif($dat[0] =~ m/bondlength_go/)
    {
      @bondlength_go = ();
      for($i=0;$i<=$#dat2;$i++) {push(@bondlength_go,$dat2[$i])} # read in a list of bond potentials, i.e., PDB segment is Go or Non-go
    }
    elsif($dat[0] =~ m/dihedral_go/)
    {
      @dihedral_go = ();
      for($i=0;$i<=$#dat2;$i++) {push(@dihedral_go,$dat2[$i])} # read in a list of dihedral potentials, i.e., PDB segment is Go or Non-go
    }
    elsif($dat[0] =~ m/charges/)  
    {
      @charges = ();
      for($i=0;$i<=$#dat2;$i++) {push(@charges,$dat2[$i])} # read in charge status for PDB   
    }
    elsif($dat[0] =~ m/angle_dw/)    
    {
      @angle_dw = ();
      for($i=0;$i<=$#dat2;$i++) {push(@angle_dw,$dat2[$i])} # read in a list of angle double-well potentials, i.e., PDB segment uses it or not   
    }
    elsif($dat[0] =~ m/improperdihed_go/)
    {
      @improperdihed_go = ();
      for($i=0;$i<=$#dat2;$i++) {push(@improperdihed_go,$dat2[$i])} # read in a list of improper dihedral potentials, i.e., PDB segment is Go or Non-go
    }
    elsif($dat[0] =~ m/casm/) {$casm = $dat[1]}
    elsif($dat[0] =~ m/int_scfac/) {$int_scfac = $dat[1]}
    elsif($dat[0] =~ m/int_file/) {$int_file = $dat[1]}
    elsif($dat[0] =~ m/interact_scale/) {$interact_scale = $dat[1]}
    elsif($dat[0] =~ m/mass_offset/) {$mass_offset = $dat[1]}
    elsif($dat[0] =~ m/sc_name/) 
    {
      @scp=();
      for($i=0;$i<=$#dat2;$i++) {push(@scp,$dat2[$i])}
    }
    elsif($dat[0] =~ m/ca_name/) 
    {
      @cap=();
      for($i=0;$i<=$#dat2;$i++) {push(@cap,$dat2[$i])}
    }
    else 
    {
      $error_unrecognized_command = "true";
      $command = $_;
    }
  }
}
close IN;

if($#ARGV+1 > 1)
{
  print "Additional command line arguemnts supplied,\n";
  print "commands in control file may be overridden.\n";
  for($i=1;$i<=$#ARGV;$i++)
  {
    ($ARGV[$i]) = cleanws($ARGV[$i]);
    @dat = split(/=/,$ARGV[$i]);
    if($dat[0] =~ m/prefix/) {$prefix = $dat[1]}
    elsif($dat[0] =~ m/pdb/) {push(@pdb,$dat[1])}
    elsif($dat[0] =~ m/charmm/) {$charmm = $dat[1]}
    elsif($dat[0] =~ m/nscal/) {$nscal = $dat[1]}
    elsif($dat[0] =~ m/fnn/) {$fnn = $dat[1]}
    elsif($dat[0] =~ m/pot/) {push(@pot,$dat[1])}
    elsif($dat[0] =~ m/ndomain/) {$ndomain = $dat[1]}
    elsif($dat[0] =~ m/casm/) {$casm = $dat[1]}
    elsif($dat[0] =~ m/bondlength_go/) {push(@bondlength_go,$dat[1])}
    elsif($dat[0] =~ m/dihedral_go/) {push(@dihedral_go,$dat[1])}
    elsif($dat[0] =~ m/charges/) {push(@charges,$dat[1])}
    elsif($dat[0] =~ m/angle_dw/) {push(@angle_dw,$dat[1])}
    elsif($dat[0] =~ m/improperdihed_go/) {push(@improperdihed_go,$dat[1])}
    elsif($dat[0] =~ m/int_scfac/) {$int_scfac = $dat[1]}
    elsif($dat[0] =~ m/int_file/) {$int_file = $dat[1]}
    elsif($dat[0] =~ m/interact_scale/) {$interact_scale = $dat[1]}
    elsif($dat[0] =~ m/mass_offset/) {$mass_offset = $dat[1]}
    elsif($dat[0] =~ m/sc_name/) {@scp = $dat[1]}
    elsif($dat[0] =~ m/ca_name/) {@cap = $dat[1]}
    else
    {
      $error_unrecognized_command = "true";
      $command = $_;
    }
  }
}

print "\n##############################\n";
print "Build CG Model: Version $version\n";
print "Edward P. O'Brien Jr.\n";
print "Dept. of Chemistry\n";
print "Penn State University\n";
print "##############################\n\n";

print "Default values:\n";
print "pdb = @pdb\n";
print "prefix = $prefix\n";
print "charmm = $charmm\n";
print "casm = $casm\n";
print "bondlength_go = @bondlength_go\n";
print "dihedral_go = @dihedral_go\n";
print "improperdihed_go = @improperdihed_go\n";
print "charges = @charges\n";
print "nscal = $nscal\n";
print "fnn = $fnn\n";
print "pot = @pot\n";
print "ndomain = @ndomain\n";
print "angle_dw = @angle_dw\n";
print "mass_offset = $mass_offset\n";
print "sc_name = @scp\n";
print "ca_name = @cap\n";
print "interact_scale = $interact_scale\n";

if($#pdb != $#pot) {die "ERROR: pdb != pot, $#pdb != $#pot\n"}
@pot_list=();
for($i=0;$i<=$#pot;$i++)
{
  if(uc($pot[$i]) =~ m/GENERIC/)
  {
    @dat=split(/\-/,$pot[$i]);
    if($#dat+1==1) {die "ERROR: Generic potential keyword must be invoked as 'generic-bt'\n"}
    else
    {
      if(($dat[1] ne "bt")and($dat[1] ne "mj")and($dat[1] ne "kgs"))
      {
        die "ERROR: You can only invoke Generic potential keyword as 'generic-bt' or 'generic-mj' or 'generic-kgs'\n";
      }
      else{push(@pot_list,uc($dat[1]))}      
    }
  }
  else {push(@pot_list,uc($pot[$i]))}
}
my %string = map { $_, 1 } @pot_list;
if (keys %string == 1) {$final_pot=$pot_list[0]} # all equal
else {die "ERROR: The same non-bonded potential must be used throughout the entire protein; check 'pot' keywords! @pot_list\n"}

if($#pdb != $#cap)
{
  if($#cap+1==1) 
  {
    $temp=$cap[0];
    $temp2=$scp[0];
    for($k=0;$k<=$#pdb;$k++) 
    {
      $cap[$k]=$temp;
      $scp[$k]=$temp2;
    }
  }
  else
  {
    die "ERROR: Your number of PDB files does not match your number of CA and SC names defined in your CNTRL file\n"; 
  }
}

if($#pdb+1 == 1 and $interact_scale ne "undefined") {die "ERROR: You cannot invoke the interact_scale keyword in CNTRL file with only one PDB file (@pdb) defined\n"}
### END: get info from control file ###

## BEGIND: Conditional Defaults ##
if($casm==1) 
{
  $ene_bsc = 0.03;
  $single_hbond_ene = 0.50; # energy of a hydrogen bond for everthing but helices
  $single_hbond_ene_helix = 0.50; # energy of a hydrogen bond in a helix
  } # energy of a backbone-sidechain native contact
  else
  {
    $ene_bsc = 0.37;  
  $single_hbond_ene = 0.75; # energy of a hydrogen bond for everthing but helices
  $single_hbond_ene_helix = 0.75; # energy of a hydrogen bond in a helix
}

my @dom_nscal = (); #add by Yang

if($int_file ne "none")
{
  open(IN,"$int_file") or die "ERROR: File $int_file does not exist\n";
  while(<IN>)
  {
    chomp($_);
    ($_) = cleanws($_);
    @dat = split(/\=/,$_);
    (@junk) = cleanws($dat[1]);
    @dat2 = split(/\s+/,$junk[0]);
    if($_ =~ m/scale factor/)
    {
      $int_scfac = $dat2[0];
      push@dom_nscal, $dat2[0]; #add by Yang
    }
    if($_ =~ m/domain/)
    {
      $ndomain++;
      @dat3 = split(/\-/,$dat2[0]);
      $dom[$ndomain][1] = $dat3[0];
      $dom[$ndomain][2] = $dat3[1];
      if($dom[$ndomain][1] > $dom[$ndomain][2]) 
      {
        die "ERROR: When defining the domains in the interface file $dom[$ndomain][1] is Greater than $dom[$ndomain][2]!\n";
      }
    }
  }
  close IN;
  print "$ndomain domain(s) defined in the Interface file $int_file\n";
  print "Domain information:\n";
  for($i=1;$i<=$ndomain;$i++)
  {
    print "Domain $i: $dom[$i][1] to $dom[$i][2]\n";
  }
}

## Loop-up table for uniquely indentifying residues #
@aa = ("GLY","ALA","VAL","LEU","ILE","MET","PHE","PRO","SER","THR","CYS","ASN","GLN","TYR","TRP","ASP","GLU","HIS","LYS","ARG");
if($#aa+1 != 20) {die "ERROR\n"}
for($i=0;$i<=$#aa;$i++)
{
  $num = $i+1;
  $res{$aa[$i]} = $num;
  $n2res{$num} = $aa[$i];
}

%Mass =
(
  "N" => 14.0067,
  "H" => 1.00794,
  "C" => 12.011,
  "O" => 15.9994,
  "S" => 32.06,
  );
# number of heavy atoms in sidechains
%refNscat =
(
  "ALA" => 1,
  "CYS" => 2,
  "ASP" => 4,
  "GLU" => 5,
  "PHE" => 7,
  "GLY" => 0,
  "HIS" => 6,
  "HSD" => 6,
  "HSE" => 6,
  "HSP" => 6,
  "ILE" => 4,
  "LYS" => 5,
  "LEU" => 4,
  "MET" => 4,
  "ASN" => 4,
  "PRO" => 3,
  "GLN" => 5,
  "ARG" => 7,
  "SER" => 2,
  "THR" => 3,
  "VAL" => 3,
  "TRP" => 10,
  "TYR" => 8
  );
# charges on side chains at pH 7
%refcharge =
(
  "ALA" => 0.0,
  "CYS" => 0.0,
  "ASP" => -1.0,
  "GLU" => -1.0,
  "PHE" => 0.0,
  "GLY" => 0.0,
  "HIS" => 0.0,
  "HSD" => 0.0,
  "HSE" => 0.0,
  "HSP" => 0.0,
  "ILE" => 0.0,
  "LYS" => 1.0,
  "LEU" => 0.0,
  "MET" => 0.0,
  "ASN" => 0.0,
  "PRO" => 0.0,
  "GLN" => 0.0,
  "ARG" => 1.0,
  "SER" => 0.0,
  "THR" => 0.0,
  "VAL" => 0.0,
  "TRP" => 0.0,
  "TYR" => 0.0
  );

# Generic C_alpha side-chain center of mass distance
%lbs_nongo =
(
  "ASP" => 2.46916481058687,
  "PRO" => 1.87381801537346,
  "LYS" => 3.49738414814426,
  "ILE" => 2.25260184847053,
  "TRP" => 3.58251993741888,
  "CYS" => 2.06666004558289,
  "HSD" => 3.15209719417679,
  "PHE" => 3.38385541816659,
  "HSP" => 3.15209719417679,
  "GLN" => 3.08654121335,
  "SER" => 1.89840600762153,
  "ASN" => 2.46916481058687,
  "VAL" => 1.93953811063784,
  "LEU" => 2.56580983973678,
  "TYR" => 3.38981664391425,
  "GLU" => 3.07971386504681,
  "ARG" => 3.39687572938579,
  "THR" => 1.931721703272,
  "ALA" => 1.51146031725997,
  "MET" => 2.95389402456081,
  "HIS" => 3.15209719417679,
  "HSE" => 3.15209719417679
  );

%improper_nongo =
(
  "ASP" => 14.655341300544,
  "PRO" => 26.763068425539,
  "LYS" => 12.765248692601,
  "ILE" => 13.5446902008313,
  "TRP" => 11.4483488626106,
  "CYS" => 20.0484470024042,
  "HSD" => 14.9962640689562,
  "PHE" => 10.9217771918902,
  "HSP" => 14.9962640689562,
  "GLN" => 17.3050853491068,
  "SER" => 20.1390130256255,
  "ASN" => 14.655341300544,
  "VAL" => 13.3216022614598,
  "LEU" => 11.8137180266206,
  "TYR" => 12.2715081962165,
  "GLU" => 15.4130821146834,
  "ARG" => 15.5451613009777,
  "THR" => 16.2956083930276,
  "ALA" => 16.8418866013662,
  "MET" => 12.7046284165739,
  "HIS" => 14.9962640689562,
  "HSE" => 14.9962640689562
  );

%ang_sb_nongo =
(
  "ASP" => 120.380153696218,
  "PRO" => 125.127927161651,
  "LYS" => 119.523270610009,
  "ILE" => 118.791108398805,
  "TRP" => 130.018548241749,
  "CYS" => 110.512719347428,
  "HSD" => 116.815900172681,
  "PHE" => 122.937540996701,
  "HSP" => 116.815900172681,
  "GLN" => 116.182123224059,
  "SER" => 107.971234136647,
  "ASN" => 120.380153696218,
  "VAL" => 112.877421898116,
  "LEU" => 123.32179171436,
  "TYR" => 116.783314494739,
  "GLU" => 116.659068554985,
  "ARG" => 119.709740783191,
  "THR" => 111.719883260793,
  "ALA" => 108.623605160075,
  "MET" => 116.636559053295,
  "HIS" => 116.815900172681,
  "HSE" => 116.815900172681
  );

%ang_bs_nongo =
(
  "ASP" => 116.629356207687,
  "PRO" => 79.4932105625367,
  "LYS" => 119.779735484239,
  "ILE" => 116.923861483529,
  "TRP" => 100.858690902849,
  "CYS" => 114.816253227757,
  "HSD" => 115.848569293979,
  "PHE" => 112.804608190743,
  "HSP" => 115.848569293979,
  "GLN" => 119.106753006548,
  "SER" => 116.361829754186,
  "ASN" => 116.629356207687,
  "VAL" => 121.299281732077,
  "LEU" => 117.587011217416,
  "TYR" => 116.72484692836,
  "GLU" => 119.507585037498,
  "ARG" => 117.532816176021,
  "THR" => 117.044133956143,
  "ALA" => 120.747734648009,
  "MET" => 123.234171432545,
  "HIS" => 115.848569293979,
  "HSE" => 115.848569293979
  );

# segment id relationships
@alphabet = ('A'..'Z');
$nseg = 0;
foreach $letter (@alphabet)
{
  $nseg++;
  $segid2num{"$letter"}=$nseg;   
}
# Miyazawa-Jernigan statistical potential values
# OR Bentancourt-Thirumalai potential
# get contact potential minimum energies
$dir = getcwd;
@dat = split(/\//,$dir);
$root = join('/',@dat[0..4]);
#if(uc($final_pot) =~ m/MJ/) {$miya = "$root/software/shared_files/mj_contact_potential.dat"}
#elsif(uc($final_pot) =~ m/KGS/) {$miya = "$root/software/shared_files/kgs_contact_potential.dat"}
#elsif(uc($final_pot) =~ m/BT/) {$miya = "$root/software/shared_files/bt_contact_potential.dat"}
#elsif(uc($final_pot) =~ m/CT/) {$miya = "$root/software/shared_files/mj_contact_potential.dat"} # Dirty, temporary solution
#elsif(uc($final_pot) =~ m/GENERIC/) {$miya = "$root/software/shared_files/mj_contact_potential.dat"} # Dirty, temporary solution
#else {die "ERROR: Unrecognized force-field $final_pot\n";}
if(uc($final_pot) =~ m/MJ/) {$miya = "$ENV{CG_MODEL_HOME}/shared_files/mj_contact_potential.dat"}
elsif(uc($final_pot) =~ m/KGS/) {$miya = "$ENV{CG_MODEL_HOME}/shared_files/kgs_contact_potential.dat"}
elsif(uc($final_pot) =~ m/BT/) {$miya = "$ENV{CG_MODEL_HOME}/shared_files/bt_contact_potential.dat"}
elsif(uc($final_pot) =~ m/CT/) {$miya = "$ENV{CG_MODEL_HOME}/shared_files/mj_contact_potential.dat"} # Dirty, temporary solution
elsif(uc($final_pot) =~ m/GENERIC/) {$miya = "$ENV{CG_MODEL_HOME}/shared_files/mj_contact_potential.dat"} # Dirty, temporary solution
else {die "ERROR: Unrecognized force-field $final_pot\n";}
open(IN,"$miya") or die "ERROR: file $miya does not exist\n";
$nrows = 0;
$avgmj = 0;
$nmj = 0;
while(<IN>)
{
  if($_ !~ m/#/)
  {
    @dat = split(/\s+/,$_);
    if($_ =~ m/AA/)
    {
      for($i=1;$i<=$#dat;$i++)
      {
        $vec[$i] = $res{uc($dat[$i])};
      }
      if($#dat != 20) {die "ERROR: missing residues in file $miya\n"}
    }
    else
    {
      $nrows++;
      $tc = 1;
      for($i=0;$i<=$#dat;$i++)
      {
        if(uc($final_pot) =~ m/MJ/)
        {
          $eps[$vec[$nrows]][$vec[$tc]] = $nscal*abs($dat[$i]-1.2);
          $eps[$vec[$tc]][$vec[$nrows]] = $nscal*abs($dat[$i]-1.2);
          $avg_mj += $nscal*abs($dat[$i]-1.2);
        }
        elsif(uc($final_pot) =~ m/BT/)
        {
          $eps[$vec[$nrows]][$vec[$tc]] = $nscal*abs($dat[$i]-0.6);
          $eps[$vec[$tc]][$vec[$nrows]] = $nscal*abs($dat[$i]-0.6);
          $avg_mj += $nscal*abs($dat[$i]-0.6);
        }
        elsif(uc($final_pot) =~ m/KGS/)
        {
          $eps[$vec[$nrows]][$vec[$tc]] = $nscal*abs($dat[$i]-1.8);
          $eps[$vec[$tc]][$vec[$nrows]] = $nscal*abs($dat[$i]-1.8);
          $avg_mj += $nscal*abs($dat[$i]-0.6);
        }
        elsif(uc($final_pot) =~ m/CT/)
        {
          $eps[$vec[$nrows]][$vec[$tc]] = $nscal;
          $eps[$vec[$tc]][$vec[$nrows]] = $nscal;
          $avg_mj += $nscal;
        }
        elsif(uc($final_pot) =~ m/GENERIC/) {}
        else {die "ERROR: unrecognized force field: $final_pot\n"}
        $tc++;
        $nmj++;
      }
      if($nrows > 20) {die "ERROR 2: missing residues in file $miya: $nrows \n"}
      if($#dat + 1 != $nrows) {die "ERROR 3: missing residues in file $miya, $#dat + 1 != $nrows \n"}
    }
  }
}
close IN;
$avg_mj = $avg_mj/$nmj;
print "The average $final_pot interaction energy is $avg_mj\n";

# Read in the generic backbone dihedral potential of CL Brooks if NON-GO dihedrals
# requested by user.
print "Reading in dihedrals for each PDB\n";
if($dihedral_go[$np] == 0) 
{
 #$file = "$root/software/shared_files/karanicolas_dihe_parm.dat";
 $file = "$ENV{CG_MODEL_HOME}/shared_files/karanicolas_dihe_parm.dat";
 open(IN,"$file") or die "ERROR: file $file does not exist\n";
 while(<IN>)
 {
   chomp($_);
   @dat = split(/\s+/,$_);
   $r1 = uc($dat[0]);
   $r2 = uc($dat[1]);
   if($r1 ne $r1_old or $r2 ne $r2_old) {$nphi=0}
   $nphi++;
   $dihedb_nongo[$res{$r1}][$res{$r2}][$nphi][1] = 0.756*$dat[2]; # We scale to match the values from the MMTSB website
   $dihedb_nongo[$res{$r1}][$res{$r2}][$nphi][2] = $dat[3];
   $dihedb_nongo[$res{$r1}][$res{$r2}][$nphi][3] = $dat[4];
   $r1_old = $r1;
   $r2_old = $r2;
   if($nphi > 4) {die "ERROR: nphi = $nphi upon reading in generic dihedral file\n $_\n"}
 }
}

# mass of amino acids
# UNSURE! about pro, arg, his and cys weights
%aaSCmass = (
  "ALA" => 71.000000,
  "CYS" => 114.000000,
  "ASP" => 114.000000,
  "GLU" => 128.000000,
  "PHE" => 147.000000,
  "GLY" => 57.000000,
  "HIS" => 114.000000,
  "HSD" => 114.000000,
  "HSE" => 114.000000,
  "HSP" => 114.000000,
  "ILE" => 113.000000,
  "LYS" => 128.000000,
  "LEU" => 113.000000,
  "MET" => 131.000000,
  "ASN" => 114.000000,
  "PRO" => 114.000000,
  "GLN" => 128.000000,
  "ARG" => 114.000000,
  "SER" => 87.000000,
  "THR" => 101.000000,
  "VAL" => 99.000000,
  "TRP" => 186.000000,
  "TYR" => 163.000000
  );

# vdw radius of sidechains
%rvdw =
(
  "ALA" => 2.51958406732374,
  "CYS" => 2.73823091624513,
  "ASP" => 2.79030096923572,
  "GLU" => 2.96332591119925,
  "PHE" => 3.18235414984794,
  "GLY" => 2.25450393833984,
  "HIS" => 3.04273820988499,
  "HSD" => 3.04273820988499,
  "HSE" => 3.04273820988499,
  "HSP" => 3.04273820988499,
  "ILE" => 3.09345983013354,
  "LYS" => 3.18235414984794,
  "LEU" => 3.09345983013354,
  "MET" => 3.09345983013354,
  "ASN" => 2.84049696898525,
  "PRO" => 2.78004241717965,
  "GLN" => 3.00796101305807,
  "ARG" => 3.28138980397453,
  "SER" => 2.59265585208464,
  "THR" => 2.81059478021734,
  "VAL" => 2.92662460060742,
  "TRP" => 3.38869998431408,
  "TYR" => 3.22881842919248
  );

# Loop over each PDB file read in
for($np=0;$np<=$#pdb;$np++)
{
  print "Reading in PDB file $pdb[$np]\n";
  $nca[$np]=0;
  $nsc=0;
  @nb=();
  # (1) format PDB using MMTSB
  ##create output file name
  ($fileBaseName, $dirName, $fileExtension) = fileparse($pdb[$np], ('\.pdb') );
  $fileBaseName = lc($fileBaseName);
  $out = "$fileBaseName" . "_heavy_atoms.pdb";
  system "convpdb.pl -nsel heavy $pdb[$np] > temp.pdb";
  system "convpdb.pl -renumber 1 -nohetero -out charmm27 temp.pdb > $out";

  # (2) extract coordinates on Calpha and
  # sequence information (used for assigning dihedrals)
  open(IN,"$out") or die "ERROR: file $out does not exist\n";
  $hat=0;
  $multi_seg = 0; # 0 -> only one segment in PDB, 1 -> more than one segment in PDB
  $res_num = 0;
  while(<IN>)
  {
    if($_ =~ m/ATOM/) # Only read in lines that start with 'ATOM'
    {
      $hat++;
      @dat = split(/\s+/,$_);
      if($hat == 1) # run this test only once on the first atom's coordinates
      {
        $xtest = $dat[6];
        $ytest = $dat[7];
        $ztest = $dat[8];
      }
      if(int($xtest)!=$xtest and int($ytest)!=$ytest and int($ztest)!=$ztest) # Good, they're all floating points
      {
        $atomnam[$hat] = $dat[2];
        @chars=split(//,$dat[2]);
        $atom[$hat] = $chars[0];
        $resnam[$np][$hat] = $dat[3];
        $segnam[$hat] = $dat[4];
        $segid[$hat] = $segid2num{"$segnam[$hat]"};
        $junk[$hat] = $dat[5];
        if($junk[$hat] != $junk[$hat-1]) {$res_num++}
        $resnum[$hat] = $res_num;
        $resnam_num[$np][$res_num] = $resnam[$np][$hat];
        $map_pdbnum_newnum[$segid[$hat]][$dat[5]] = $res_num;
        $x[$hat] = $dat[6];
        $y[$hat] = $dat[7];
        $z[$hat] = $dat[8];
        print "$atom[$hat], $resnam[$np][$hat], $segnam[$hat], $atomnam[$hat], $junk[$hat], $res_num, $x[$hat], $y[$hat], $z[$hat]\n";
      }
      else # Assume xyz's are shifted by one index
      {
        $atomnam[$hat] = $dat[2];
        @chars=split(//,$dat[2]);
        $atom[$hat] = $chars[0];
        $resnam[$np][$hat] = $dat[3];
        $segnam[$hat] = $dat[-1];
        $segid[$hat] = $segid2num{"$segnam[$hat]"};
        $junk[$hat] = $dat[4];
        if($junk[$hat] != $junk[$hat-1]) {$res_num++}
        print "$atom[$hat], $resnam[$np][$hat], $segnam[$hat], $atomnam[$hat], $junk[$hat], $res_num\n";
        $resnum[$hat] = $res_num;
        $resnam_num[$np][$res_num] = $resnam[$np][$hat];
        $map_pdbnum_newnum[$segid[$hat]][$dat[4]] = $res_num;
        $x[$hat] = $dat[5];
        $y[$hat] = $dat[6];
        $z[$hat] = $dat[7];
      }

      # Determine if more than one protein segment is present in PDB
      if($hat > 1 and ($segnam[$hat] ne $segnam[$hat-1])) {$multi_seg=1}
      # ID backbone, and SC atoms
      # and in particular Calpha
      if($atomnam[$hat] eq 'N' or $atomnam[$hat] eq 'CA' or $atomnam[$hat] eq 'C' or $atomnam[$hat] eq 'O' or $atomnam[$hat] eq 'OT1' or $atomnam[$hat] eq 'OT2') #or ($atomnam[$hat] eq 'CD' and $resnam[$np][$hat] eq "PRO"))
      {
        $idbb[$hat] = 1;
        if($attype[$hat-1] eq 's') 
        {
          $send[$resnum[$hat-1]] = $hat-1;
        }
        $attype[$hat] = "b";
        $nb[$resnum[$hat]]++;
        if($atomnam[$hat] eq 'CA') 
        {
          $nca[$np]++;
          $ca[$nca[$np]] = $resnum[$hat]; 
          $xb[$nca[$np]] = $x[$hat];
          $yb[$nca[$np]] = $y[$hat];
          $zb[$nca[$np]] = $z[$hat];
          $xb2[$np][$nca[$np]] = $x[$hat];
          $yb2[$np][$nca[$np]] = $y[$hat];
          $zb2[$np][$nca[$np]] = $z[$hat];
          $seg[$np][$nca[$np]] = $segid[$hat];
          if($resnam[$np][$hat] eq 'GLY') 
          {
            $ngly++;
            $nsc++;
            $gly[$nsc] = 1; # yes,  glycine
            $ns[$np][$resnum[$hat]] = 0;
          }
        }
      }
      elsif($resnam[$np][$hat] ne 'GLY') 
      {
        $attype[$hat] = "s";
        $ns[$np][$resnum[$hat]]++;
        if($ns[$np][$resnum[$hat]] == 1) 
        {
          $nsc++;
          $gly[$nsc] = 0; # no, not glycine
          $sstrt[$resnum[$hat]] = $hat;
        }
      }
    }
  }  
  close IN;

  for($i=1;$i<=$nca[$np];$i++) {print "$i: $gly[$i], $ns[$np][$i]\n"}
  if($attype[$hat] eq 's')
  {
    $send[$resnum[$hat]] = $hat;
  }
  if($multi_seg == 1) {$maxres[$np] = $nca[$np]}
  else {$maxres[$np] = $resnum[$hat]}
  $maxat = $hat;
  print "number of glycines = $ngly\n";
  print "maxres = $maxres[$np], maxatoms = $maxat\n";

  # Identify what domain a residue is in
  if($ndomain[$np] > 0)
  {
    for($i=1;$i<=$maxres[$np];$i++)
    {
      for($j=1;$j<=$ndomain[$np];$j++)
      {      
        if($i >= $dom[$j][1] and $i <= $dom[$j][2])
        {
          $id_dom[$np][$i] = $j;
          $check[$i]++;
          if($check[$i] > 1) {die "ERROR: Residue $i has been found in more than one domain definition! $i >= $dom[$ndomain[$np]][1] and $i <= $dom[$ndomain[$np]][2]\n"} 
        }
      }
      if($check[$i] == 0) # residue i not in any domain defined in the interface file
      {
        $id_dom[$np][$i] = $ndomain[$np]+1;
      }
    }
  }

  ## quality check, nca = nres?
  if($nca[$np] != $maxres[$np]) {die "ERROR: nca != maxres, $nca[$np] != $maxres[$np]\n"}
  ## quality check, any backbone atoms missing?
  for($i=1;$i<=$maxres[$np];$i++)
  {
    if($nb[$i] != 4 and $i !=$maxres[$np]) 
    {
      die "ERROR: In pdb $pdb[$np] (# $np) the number of backbone atoms in residue $i incorrect: $nb[$i] != 4\n";
    }
    elsif(($nb[$i] != 4 and $nb[$i] != 5) and $i == $maxres[$np])
    {
      die "ERROR: In pdb $pdb[$np] the number of backbone atoms in residue $i incorrect: $nb[$i] != 4\n";
    }
  }
  ## quality check, any heavy atoms missing?
  for($i=1;$i<=$maxres[$np];$i++)
  {
    if($refNscat{uc($resnam_num[$np][$i])} != $ns[$np][$i]) {die "ERROR 1: number of heavy atoms in sidechain $i ($resnam_num[$np][$i]) incorrect: $refNscat{uc($resnam_num[$np][$i])} != $ns[$np][$i]\n"}
  }

  # (3) get Sidechain center-of-mass
  if($casm == 1) # Only perform this QC and calculation if we are using the calpha_side chain model
  {
    for($i=1;$i<=$maxres[$np];$i++)
    {
      if($ns[$np][$i] > 0)
      {
        if(defined($sstrt[$i]) and defined($send[$i])) 
        {
          &Determine_CM($sstrt[$i],$send[$i]);
          $xs[$i] = $Rcm[0];
          $ys[$i] = $Rcm[1];
          $zs[$i] = $Rcm[2];
          $xs2[$np][$i] = $Rcm[0];
          $ys2[$np][$i] = $Rcm[1];
          $zs2[$np][$i] = $Rcm[2];
        }
        else {die "ERROR: start and finish of sidechain $i not properly defined: $sstrt[$i],$send[$i])\n"}
      }
    }
  }

  # (3) Compute bond lengths, and angles, diherals, and sidechain sizes
  ##bond lengths, and sidechain sizes
  for($i=1;$i<=$nca[$np];$i++)
  {
    #ca-ca bond length
    if($i != $nca[$np]) 
    { 
      if($seg[$np][$i] == $seg[$np][$i+1])
      {
        if($bondlength_go[$np] == 1) {$lbb[$np][$i] = distance($xb[$i],$yb[$i],$zb[$i],$xb[$i+1],$yb[$i+1],$zb[$i+1])}
        else {$lbb[$np][$i] = 3.81}
        if($lbb[$np][$i] < 2.8 or $lbb[$np][$i] > 4) {die "ERROR: Ca-Ca length abnormal between residues $i and $i+1: $lbb[$np][$i]: $xb[$i],$yb[$i],$zb[$i],$xb[$i+1],$yb[$i+1],$zb[$i+1]\n"}
      }
    }
    #ca-cb bond length
    if($gly[$i] == 0 and $casm==1)
    {
      if($bondlength_go[$np]==1) {$lbs[$np][$i] = distance($xb[$i],$yb[$i],$zb[$i],$xs[$i],$ys[$i],$zs[$i])}
      else { $lbs[$np][$i] = $lbs_nongo{uc($resnam_num[$np][$i])} }
      if($lbs[$np][$i] > 7) {die "ERROR: Ca-Cb length abnormal for pdb = $np, residue $i: $lbs[$np][$i] = $xb[$i],$yb[$i],$zb[$i],$xs[$i],$ys[$i],$zs[$i]\n"}
      elsif($lbs[$np][$i] == 0) {die "ERROR: Ca-Cb length abnormal for residue $i equals 0\n"}

      #???????????
      # effective rmin values of repulsive non-native sidechain contacts
      $t1 = abs($eps[$res{$resnam_num[$np][$i]}][$res{$resnam_num[$np][$i]}]);
      $t2 = ($t1*(2*$lbs[$np][$i]*2**(1/6))**12/(1e-12))**(1/12);
      $rmin[$i] = $t2;
    }
  }
  ## bond angles
  # for backbone
  for($i=1;$i<=$nca[$np]-2;$i++)
  {
    if($angle_dw[$np] == 0) {$ang_b[$np][$i] = Angle($xb[$i],$yb[$i],$zb[$i],$xb[$i+1],$yb[$i+1],$zb[$i+1],$xb[$i+2],$yb[$i+2],$zb[$i+2])}
  }

  # for side chains
  if($casm==1)
  {
    for($i=1;$i<=$nca[$np]-1;$i++)
    {
      if($ns[$np][$i] > 0)
      {
        if($angle_dw[$np] == 0)
        {
          $ang_sb[$np][$i] = Angle($xs[$i],$ys[$i],$zs[$i],$xb[$i],$yb[$i],$zb[$i],$xb[$i+1],$yb[$i+1],$zb[$i+1]);
        }
        else
        {
          $ang_sb[$np][$i] = $ang_sb_nongo{uc($resnam_num[$np][$i])}; # use the generic parameters
        }
      }

      if($ns[$np][$i+1] > 0)
      {
        if($angle_dw[$np] == 0)
        {              
          $ang_bs[$np][$i] = Angle($xb[$i],$yb[$i],$zb[$i],$xb[$i+1],$yb[$i+1],$zb[$i+1],$xs[$i+1],$ys[$i+1],$zs[$i+1]);
        }               
        else              
        {               
          $ang_bs[$np][$i] = $ang_bs_nongo{uc($resnam_num[$np][$i+1])}; # use the generic parameters
        }
      }
    }
  }
  if($dihedral_go[$np] == 1) # Use go dihedrals for the backbone
  {
    ## backbone dihedrals
    for($i=1;$i<=$nca[$np]-3;$i++) { $dihedb[$np][$i] = get_Torsion($i,$i+1,$i+2,$i+3,b) }
  }
  if($casm==1)
  {
    ## sidechain/backbone dihedrals - to maintain chirality
    for($i=2;$i<=$nca[$np]-1;$i++) 
    {
      if($ns[$np][$i] > 0 and $casm==1) 
      {
        # ordering is important!
        if($improperdihed_go[$np]==1) {$diheds[$np][$i] = get_Torsion($i,$i-1,$i+1,$i,"s")} # bi, bi-1, bi+1, si 
        else 
        {
          $diheds[$np][$i] = $improper_nongo{uc($resnam_num[$np][$i])}; #} # use transferable improper dihedral
          print "Here $np, $i, $diheds[$np][$i] = $improper_nongo{uc($resnam_num[$np][$i])}\n";
        }
      }
    }
  }

  if(uc($pot[$np]) !~ m/GENERIC/)
  {
    # (4) Compute native contacts between side-chains
    print "Determining native contacts\n";
    for($i=1;$i<=$nca[$np];$i++)
    {
      for($j=1;$j<=$nca[$np];$j++) {$distances[$np][$i][$j] = distance($xb[$i],$yb[$i],$zb[$i],$xb[$j],$yb[$j],$zb[$j])}
    }
    print "Finished calculating distance matrix\n";
    for($i=1;$i<=$nca[$np]-3;$i++) # sidechain 1
    {
      for($j=$i+3;$j<=$nca[$np];$j++) # sidechain 2
      {
        $dij = distance($xb[$i],$yb[$i],$zb[$i],$xb[$j],$yb[$j],$zb[$j]);
        if($gly[$i] == 0 and $gly[$j] == 0)
        {
          $found=0;
          for($k=$sstrt[$i];$k<=$send[$i];$k++)
          {
            for($l=$sstrt[$j];$l<=$send[$j];$l++)
            {
              $dij = distance($x[$k],$y[$k],$z[$k],$x[$l],$y[$l],$z[$l]);
              if($dij <= $heav_cut and $found == 0)
              {
                $found = 1;
                $nnatsc[$np]++;
                $natss[$i][$j] = 1;
                $nsc1[$np][$nnatsc[$np]] = $i;
                $nsc2[$np][$nnatsc[$np]] = $j;
                $ncsc[$i]++;
                $ncsc[$j]++;
                if($casm==1) {$dij = distance($xs[$i],$ys[$i],$zs[$i],$xs[$j],$ys[$j],$zs[$j])}
                else {$dij = distance($xb[$i],$yb[$i],$zb[$i],$xb[$j],$yb[$j],$zb[$j])}
                $natdist[$np][$nnatsc[$np]] = $dij;
                $native[$np][$i][$j]=1;
              }
            }
          }
          if($found == 0)  
          {
            $nonnatsc++;
            $natss[$i][$j] = 0;
          }
        } 
      }
    }

    # (4) Compute native contacts between backbone and side-chains
    print "Determining Backbone - Side chain contacts\n";
    for($i=1;$i<=$nca[$np];$i++) # sidechain 1
    {
      for($j=1;$j<=$hat;$j++) # loop through all heavy atoms
      {
        $resid = $resnum[$j]; # id residue of heavy atom
        if($i < $resid - 2 or $i > $resid + 2)
        {
          if($gly[$i] == 0 and $idbb[$j] == 1 and $found[$resid][$i] != 1) # is it a backbone atom?
          {
            for($k=$sstrt[$i];$k<=$send[$i];$k++) # loop through sidechain $i's heavy atoms
            {
              $dij = distance($x[$k],$y[$k],$z[$k],$x[$j],$y[$j],$z[$j]); # compute the distance
              if($dij <= $heav_cut and $found[$resid][$i] != 1)
              {
                $found[$resid][$i] = 1;
                if($casm==1) {$dij = distance($xs[$i],$ys[$i],$zs[$i],$xb[$resid],$yb[$resid],$zb[$resid])}
                else {$dij = distance($xb[$i],$yb[$i],$zb[$i],$xb[$resid],$yb[$resid],$zb[$resid])}
                $nbsc[$np]++;
                $natbsc[$np][$nbsc[$np]][1] = $resid; # backbone ID
                $natbsc[$np][$nbsc[$np]][2] = $i; # sidechain ID
                $natdist_bsc[$np][$nbsc[$np]] = $dij;
                $native[$np][$i][$resid]=1;
              }      
            }
          }
        }
      }
    }

    print "# nat sc-sc contacts $nnatsc[$np], # nat sc-bb contacts $nbsc[$np], and  # non-nat sc-sc $nonnatsc[$np]\n";
    &visualwritepdb;

    open(CM,">$fileBaseName" . "_NumContacts_perSC.dat");
    for($i=1;$i<=$nca[$np];$i++)
    {
      print CM "$i $ncsc[$i]\n";
    }
    close CM;

    # (5) Determine hydrogen bonds that are present using STRIDE,
    #   and assign to Calpha-Calpha pairs. Also secondary structural elements
    #   within the native structure.   
    print "Determining the presence of hydrogen bonds using STRIDE\n";
    $nhbond[$np]=0;
    if($maxres[$np] > 4 and uc($pot[$np]) !~ m/GENERIC/)
    {
      &Stride; # Analyze structure with stride
      if($#D < 1)
      {
        # This is probably because there were no H-bonds in the structure, Stride dies in this case
        # Remedy it by adding a ghost structure with a H-bond
        die "Error: File had no info: nconverted = $j, $#D\n";
      } 
      $nr = 0;
      foreach $line (@D)
      {
        # Determine location and identity of secondary structural elements
        if($line =~ m/ASG /)
        {
          $nr++;
          # Determine per residue information
          $helical[$np][$nr]=0;
          @dat = split(/\s+/,$line);
          $res[$nr][1] = $dat[1];
          $res[$nr][2] = $dat[3];
          $res[$nr][3] = $dat[4];
          $res[$nr][4] = $nam2id{$dat[5]};
          $o2n{$dat[3]} = $dat[4];
          if($dat[6] =~ /Helix/) {$helical[$np][$nr]=1}
        }
        ## ASSUME!!! single chain, globular protein
        if( $line =~ m/^ACC / or $line =~ m/^DNR / )
        {
          # Get H-bonding info
          $RESnum1 = substr($line,11,4);
          $RESnum2 = substr($line,31,4);
          $segment1 = substr($line,8,2);
          $segment2 = substr($line,28,2);
          # Remove white spaces at begining of string
          $RESnum1 =~ s/^\s+//;
          $RESnum2 =~ s/^\s+//;
          $segment1 =~ s/^\s+//;
          $segment2 =~ s/^\s+//;
          # Remove white spaces at end of string
          $RESnum1 =~ s/\s+$//;
          $RESnum2 =~ s/\s+$//;
          $segment1 =~ s/\s+$//;
          $segment2 =~ s/\s+$//;
          # Renumber the chain and residue ids, $nSnap was defined earlier
          $nRes1 = $map_pdbnum_newnum[$segid2num{$segment1}][$RESnum1];
          $nRes2 = $map_pdbnum_newnum[$segid2num{$segment2}][$RESnum2];
          $native[$np][$nRes1][$nRes2]=1;
          $native[$np][$nRes2][$nRes1]=1;
          if($nRes1 > $maxres[$np] or $nRes2 > $maxres[$np]) {die "ERROR: here\n"}
          if($segment1 eq $segment2)
          {
            if($nRes1 < $nRes2) 
            {
              if(defined($hbond[$nRes1][$nRes2])) 
              {
                if($helical[$np][$nRes1] == 1 and $helical[$np][$nRes2] == 1 ) {$hbond_ene[$np][$nhbond[$np]] = 2*$single_hbond_ene_helix} # helical hbond
                else {$hbond_ene[$np][$nhbond[$np]] = 2*$single_hbond_ene} # not helical hbond
              }
              else 
              {
                $nhbond[$np]++;
                $hbond[$nRes1][$nRes2] = 1;
                $hbond_at[$np][$nhbond[$np]][1] = $nRes1;
                $hbond_at[$np][$nhbond[$np]][2] = $nRes2;
                if($helical[$np][$nRes1] == 1 and $helical[$np][$nRes2] == 1 ) {$hbond_ene[$np][$nhbond[$np]] = $single_hbond_ene_helix} # helical hbond
                else {$hbond_ene[$np][$nhbond[$np]] = $single_hbond_ene} # not helical hbond
              }
            }
          }
          else  
          {
            $nhbond[$np]++;
            $hbond[$nRes1][$nRes2] = 1;
            $hbond_at[$np][$nhbond[$np]][1] = $nRes1;
            $hbond_at[$np][$nhbond[$np]][2] = $nRes2;
            $hbond_ene[$np][$nhbond[$np]] = $single_hbond_ene;
          }
        }
      }  
      for($i=1;$i<=$nhbond[$np];$i++) {print "$i $hbond_ene[$np][$i], $hbond_at[$np][$i][1] $hbond_at[$np][$i][2]\n"}
      print "# of unique Hbonds $nhbond[$np]\n";
    }
    for($i=1;$i<=$nhbond[$np];$i++)
    {
      $n1 = $hbond_at[$np][$i][1];
      $n2 = $hbond_at[$np][$i][2];
      $dhb[$np][$i] = distance($xb[$n1],$yb[$n1],$zb[$n1],$xb[$n2],$yb[$n2],$zb[$n2]);
      print "$hbond_at[$np][$i][1],$hbond_at[$np][$i][2] $dhb[$np][$i]\n";
    }
  }

  # () create the new residue names
  for($i=1;$i<=$nca[$np];$i++)
  {
    if(uc($pot[$np]) !~ m/GENERIC/)
    {
      $resnam[$np][$i] = "G" . "$i";
      $newnam[$np][$i][1] = "$cap[$np]" . "$i";
      $newnam[$np][$i][2] = "$scp[$np]" . "$i";
    }
    elsif($interact_scale ne "undefined")
    { 
      $j=$i;
      $resnam[$np][$i] = "G" . "$j";
      $newnam[$np][$i][1] = "$cap[$np]" . "$j";
      $newnam[$np][$i][2] = "$scp[$np]" . "$j";
    }
    else
    { 
      $j=$i+$offset;
      $resnam[$np][$i] = "U" . "$j";
      $newnam[$np][$i][1] = "UA" . "$j";
      $newnam[$np][$i][2] = "UB" . "$j";
    }
  }

  if(uc($pot[$np]) =~ m/GENERIC/) {$offset=$nca[$np]}
  } # End loop over PDBs

# (6) print out new Ca-Cb crd file
if(uc($pot[$np]) !~ m/GENERIC/) {&writecrd}
# (7) print out new Ca-Cb pdb file
# (8) print out new Ca-Cb top file
if($prefix eq "none") 
{
  if($casm==1) {$out = "$fileBaseName" . "_ca-cb.top"}
  else {$out = "$fileBaseName" . "_ca.top"}
  $top_file = $out;
}
else
{
  if($casm==1) {$out = "$prefix" . "_ca-cb.top"}
  else {$out = "$prefix" . "_ca.top"}
  $top_file = $out;
}
open(OUT,">$out");
print OUT "* This CHARMM .top file describes a Ca-Cb Go model of $fileBaseName\n";
print OUT "*\n";
print OUT "20 1\n";
# create NEW topology file
## mass section
print OUT "! backbone masses\n";
$cnt = $mass_offset;
for($np=0;$np<=$#pdb;$np++)
{
  for($i=1;$i<=$nca[$np];$i++)
  {
    $cnt++;
    if($casm==1)
    {
      printf OUT "MASS %-4s%-9s%-10.6f\n",$cnt,$newnam[$np][$i][1],$aaSCmass{GLY};
      if($ns[$np][$i] > 0)
      {
        $cnt++;
        $mas = $aaSCmass{$resnam_num[$np][$i]} - $aaSCmass{GLY};
        printf OUT "MASS %-4s%-9s%-10.6f\n",$cnt,$newnam[$np][$i][2],$mas;
      }
    }
    else
    {
      printf OUT "MASS %-4s%-9s%-10.6f\n",$cnt,$newnam[$np][$i][1],$aaSCmass{$resnam_num[$np][$i]}; 
    }
  }
}
print OUT "\n";
print OUT "DECL +A\n";
print OUT "DECL -A\n";
print OUT "DECL #A";

## residue topology section
for($np=0;$np<=$#pdb;$np++)
{
  for($i=1;$i<=$maxres[$np];$i++)
  {
    $newres = $resnam[$np][$i];
    if($ns[$np][$i] > 0 and $casm==1)
    {
      print OUT "\n";
      if($charges[$np] == 0) {print OUT "RESI $newres     0.0\n"}
      else {printf OUT "RESI $newres    %2.1f\n",$refcharge{uc($resnam_num[$np][$i])}}
      print OUT "GROUP\n";
      print OUT "Atom A $newnam[$np][$i][1]     0.0\n";
      if($charges[$np] == 0) {print OUT "Atom B $newnam[$np][$i][2]     0.0\n"}
      else {printf OUT "Atom B $newnam[$np][$i][2]     %2.1f\n",$refcharge{uc($resnam_num[$np][$i])}}
      print OUT "Bond A B  A +A\n";
      print OUT "Angle -A A B  B A +A  -A A +A\n";
      print OUT "DIHE -A A +A #A\n";
      print OUT "IMPH A -A +A B\n";
    }
    else
    {
      print OUT "\n";
      if($charges[$np] == 0) {print OUT "RESI $newres     0.0\n"}
      else {printf OUT "RESI $newres     %2.1f\n",$refcharge{uc($resnam_num[$np][$i])}}
      print OUT "GROUP\n";
      if($charges[$np] == 0) {print OUT "Atom A $newnam[$np][$i][1]     0.0\n"}
      else {printf OUT "Atom A $newnam[$np][$i][1]     %2.1f\n",$refcharge{uc($resnam_num[$np][$i])}}
      print OUT "Bond A +A\n";
      print OUT "Angle -A A +A\n";
      print OUT "DIHE -A A +A #A\n";
    }
  }
}
print OUT "\nEND\n";
close OUT;

# (9) print out new Ca-Cb param file
for($i=0;$i<=$#pot;$i++) 
{
# if(uc($pot[$i]) !~ m/GENERIC/) {$temp=$pot[$i]}
$temp=$pot[$i]
}
if($prefix eq "none")
{
  $prefix = "$fileBaseName";
  if($ndomain[$np] == 0)
  {
    if(uc($temp) =~ m/MJ/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_mj.prm";}
    elsif(uc($temp) =~ m/KGS/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_kgs.prm";}
    elsif(uc($temp) =~ m/BT/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_bt.prm";}
    elsif(uc($temp) =~ m/CT/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_ct.prm";}
    elsif(uc($temp) =~ m/GENERIC/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_generic.prm";}
    else {die "ERROR 1: Unrecognized force-field $temp\n";}
  }
  else
  {
    if(uc($temp) =~ m/MJ/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_mj.prm";}
    elsif(uc($temp) =~ m/KGS/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_kgs.prm"}
    elsif(uc($temp) =~ m/BT/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_bt.prm"}
    elsif(uc($temp) =~ m/CT/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_ct.prm"}
    elsif(uc($temp) =~ m/GENERIC/) {$out = "$fileBaseName" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_generic.prm"}
    else {die "ERROR 2: Unrecognized force-field $temp\n";}
  }
}
else
{
  if($ndomain[$np] == 0)
  {
    if(uc($temp) =~ m/MJ/) {$out = "$prefix" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_mj.prm";}
    elsif(uc($temp) =~ m/KGS/) {$out = "$prefix" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_kgs.prm";}
    elsif(uc($temp) =~ m/BT/) {$out = "$prefix" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_bt.prm";}
    elsif(uc($temp) =~ m/CT/) {$out = "$prefix" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_ct.prm";}
    elsif(uc($temp) =~ m/GENERIC/) {$out = "$prefix" ."_nscal" . "$nscal" . "_fnn" . "$fnn" . "_go" . "_generic.prm";}
    else {die "ERROR 3: Unrecognized force-field $temp\n";}
  }
  else
  {
    if(uc($temp) =~ m/MJ/) {$out = "$prefix" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_mj.prm";}
    elsif(uc($temp) =~ m/KGS/) {$out = "$prefix" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_kgs.prm"}
    elsif(uc($temp) =~ m/BT/) {$out = "$prefix" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_bt.prm"}
    elsif(uc($temp) =~ m/CT/) {$out = "$prefix" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_ct.prm"}
    elsif(uc($temp) =~ m/GENERIC/) {$out = "$prefix" ."_nscal" . "$nscal" . "_interface" . "$int_scfac" . "_fnn" . "$fnn" . "_go" . "_generic.prm"}
    else {die "ERROR 4: Unrecognized force-field $temp\n";}
  }
}
$prm_file = $out;
open(OUT,">$out");
# create NEW parameter file
print OUT "* This CHARMM .param file describes a Go model of $fileBaseName\n";
print OUT "*\n";
print OUT "\n";
print OUT "ATOM\n";
$cnt = $mass_offset;
for($np=0;$np<=$#pdb;$np++)
{
  for($i=1;$i<=$nca[$np];$i++)
  {
    $cnt++;
    if($casm==1)
    {
      printf OUT "MASS %-5s %-9s%-10.6f\n",$cnt,$newnam[$np][$i][1],$aaSCmass{GLY};
      if($ns[$np][$i] > 0)
      {
        $cnt++;
        $mas = $aaSCmass{$resnam_num[$np][$i]} - $aaSCmass{GLY};
        printf OUT "MASS %-5s %-9s%-10.6f\n",$cnt,$newnam[$np][$i][2],$mas;
      }
    }
    else
    {
      printf OUT "MASS %-5s %-9s%-10.6f\n",$cnt,$newnam[$np][$i][1],$aaSCmass{$resnam_num[$np][$i]};
    }
  }
}

print OUT "\n";
## bond section
print OUT "BOND\n";
$kb = 50;
for($np=0;$np<=$#pdb;$np++)
{
  for($i=1;$i<=$maxres[$np];$i++)
  {
    if($i<=$maxres[$np]-1) # print out backbone bonding
    { 
      if($seg[$np][$i] == $seg[$np][$i+1] and $bondlength_go[$np] == 1) {printf OUT "%-8s%-10s%-12.6f%-9.6f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$kb,$lbb[$np][$i]}
      elsif($seg[$np][$i] == $seg[$np][$i+1] and $bondlength_go[$np] == 0) {printf OUT "%-8s%-10s%-12.6f%-9.6f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$kb,3.81}
    }
    if($ns[$np][$i] > 0 and $casm==1) # print out backbone-sidechain bonding
    {
      if($bondlength_go[$np] == 1) {printf OUT "%-8s%-10s%-12.6f%-9.6f\n",$newnam[$np][$i][1],$newnam[$np][$i][2],$kb,$lbs[$np][$i]}
      elsif($bondlength_go[$np] == 0) {printf OUT "%-8s%-10s%-12.6f%-9.6f\n",$newnam[$np][$i][1],$newnam[$np][$i][2],$kb,$lbs_nongo{uc($resnam_num[$np][$i])} }
    }
  }
  # Add bond between last and first residues of PDB files
  if($#pdb>0 and $np!=$#pdb) 
  {
    printf OUT "%-8s%-10s%-12.6f%-9.6f\n",$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$kb,3.81;
  }
}
print OUT "\n";
## angle section
print OUT "ANGLE\n";
$ka = 30;
for($np=0;$np<=$#pdb;$np++)
{
  for($i=1;$i<=$maxres[$np]-1;$i++)
  {
    if($i<=$maxres[$np]-2) # Ca,i - Ca,i+1 - Ca,i+2
    {
      if($seg[$np][$i] == $seg[$np][$i+1] and $seg[$np][$i] == $seg[$np][$i+2])
      {
        if($angle_dw[$np] == 0) {printf OUT "%-8s%-8s%-10s%11.6f%11.6f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+2][1],$ka,$ang_b[$np][$i]}
        else {printf OUT "%-8s%-8s%-10s  106.4 91.7 26.3 130.0 0.1 4.3\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+2][1]}
      }
    }
    if($ns[$np][$i] > 0 and $casm==1) # CB,i - Ca,i - Ca,i+1
    {
      if($seg[$np][$i] == $seg[$np][$i+1])
      {
        printf OUT "%-8s%-8s%-10s%11.6f%11.6f\n",$newnam[$np][$i][2],$newnam[$np][$i][1],$newnam[$np][$i+1][1],$ka,$ang_sb[$np][$i];
      }
    }
    if($ns[$np][$i+1] > 0 and $casm==1) # Ca,i - Ca,i+1 - Cb,i+1
    {
      if($seg[$np][$i] == $seg[$np][$i+1])
      {
        printf OUT "%-8s%-8s%-10s%11.6f%11.6f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+1][2],$ka,$ang_bs[$np][$i];
      }
    }
  }
  # Add angle energy terms between last 2 and first 2 residues of PDB files
  if($#pdb>0 and $np!=$#pdb)
  {
    if($nca[$np]>=2)
    {
      # Ca,i-1 - Ca,i - Ca,i+1
      if($angle_dw[$np] == 0) {printf OUT "%-8s%-8s%-10s%11.6f%11.6f\n",$newnam[$np][$maxres[$np]-1][1],$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$ka,100}
      else {printf OUT "%-8s%-8s%-10s  106.4 91.7 26.3 130.0 0.1 4.3\n",$newnam[$np][$maxres[$np]-1][1],$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1]}

      if($casm == 1)
      {
        # Cb,i - Ca,i - Ca,1
        printf OUT "%-8s%-8s%-10s%11.6f%11.6f\n",$newnam[$np][$maxres[$np]][2],$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$ka,$ang_sb_nongo{uc($resnam_num[$np][$maxres[$np]])};
      }
    }

    if($nca[$np]>=1)
    {
      # Ca,i - Ca,i+1 - Ca,i+2
      if($angle_dw[$np] == 0) {printf OUT "%-8s%-8s%-10s%11.6f%11.6f\n",$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$newnam[$np+1][2][1],$ka,100}
      else {printf OUT "%-8s%-8s%-10s  106.4 91.7 26.3 130.0 0.1 4.3\n",$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$newnam[$np+1][2][1]}

      if($casm == 1)
      {
        # Ca,i - Ca,1 - Cb,1
        printf OUT "%-8s%-8s%-10s%11.6f%11.6f\n",$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$newnam[$np+1][1][2],$ka,$ang_bs_nongo{uc($resnam_num[$np+1][1])};
      }
    }
  }
}
print OUT "\n";
## dihedral section
print OUT "DIHEDRAL\n";
print OUT "! backbone dihedrals\n";
for($np=0;$np<=$#pdb;$np++)
{
  if($dihedral_go[$np] == 1) # Use Go backbone dihedral angles 
  {
    for($i=1;$i<=$maxres[$np]-3;$i++)
    {
      if($seg[$np][$i] == $seg[$np][$i+1] and $seg[$np][$i] == $seg[$np][$i+2] and $seg[$np][$i] == $seg[$np][$i+3])
      {
        $delta = 1*$dihedb[$np][$i]-180;
        if($casm==1)
        {
          if($helical[$np][$i+1] == 1 and $helical[$np][$i+2] == 1) {$kd = 0.30} # helical
          else {$kd = 0.55} # not helical
          printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+2][1],$newnam[$np][$i+3][1],$kd,1,$delta;
          $delta = 3*$dihedb[$np][$i]-180;
          if($helical[$np][$i+1] == 1 and $helical[$np][$i+2] == 1) {$kd = 0.15} # helical
          else {$kd = 0.275} # not helical
          printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+2][1],$newnam[$np][$i+3][1],$kd,3,$delta;
        }
        else 
        {
          if($helical[$np][$i+1] == 1 and $helical[$np][$i+2] == 1) {$kd = 0.75} # helical
          else {$kd = 0.75} # not helical
          printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+2][1],$newnam[$np][$i+3][1],$kd,1,$delta;
          $delta = 3*$dihedb[$np][$i]-180;
          if($helical[$np][$i+1] == 1 and $helical[$np][$i+2] == 1) {$kd = 0.275} # helical
          else {$kd = 0.275} # not helical
          printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+2][1],$newnam[$np][$i+3][1],$kd,3,$delta;
        }
      }
    }
  }
  else # Use Non-go dihedrals
  {
    for($i=1;$i<=$maxres[$np]-3;$i++)
    {
      if($seg[$np][$i] == $seg[$np][$i+1] and $seg[$np][$i] == $seg[$np][$i+2] and $seg[$np][$i] == $seg[$np][$i+3])
      {
        for($j=1;$j<=4;$j++) 
        {
          $kd=$dihedb_nongo[$res{$resnam_num[$np][$i+1]}][$res{$resnam_num[$np][$i+2]}][$j][1];
          $period=$dihedb_nongo[$res{$resnam_num[$np][$i+1]}][$res{$resnam_num[$np][$i+2]}][$j][2];
          $delta=$dihedb_nongo[$res{$resnam_num[$np][$i+1]}][$res{$resnam_num[$np][$i+2]}][$j][3];
          printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$i][1],$newnam[$np][$i+1][1],$newnam[$np][$i+2][1],$newnam[$np][$i+3][1],$kd,$period,$delta;
        }
      }
    }
  }

  # Add dihedral energy terms between last 3 and first 3 residues of PDB files
  if($#pdb>0 and $np!=$#pdb)
  {
    # Ca-2, Ca-1, Ca, Ca+1
    if($nca[$np]>=3)
    {
      for($j=1;$j<=4;$j++)
      {
        $kd=    $dihedb_nongo[$res{$resnam_num[$np][$maxres[$np]-1]}][$res{$resnam_num[$np][$maxres[$np]]}][$j][1];
        $period=$dihedb_nongo[$res{$resnam_num[$np][$maxres[$np]-1]}][$res{$resnam_num[$np][$maxres[$np]]}][$j][2];
        $delta= $dihedb_nongo[$res{$resnam_num[$np][$maxres[$np]-1]}][$res{$resnam_num[$np][$maxres[$np]]}][$j][3];
        printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$maxres[$np]-2][1],$newnam[$np][$maxres[$np]-1][1],$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$kd,$period,$delta;
      }
    }
    # Ca-1, Ca, Ca+1, Ca+2
    if($nca[$np]>=2)
    {
      for($j=1;$j<=4;$j++)
      {
        $kd=    $dihedb_nongo[$res{$resnam_num[$np][$maxres[$np]]}][$res{$resnam_num[$np+1][1]}][$j][1];
        $period=$dihedb_nongo[$res{$resnam_num[$np][$maxres[$np]]}][$res{$resnam_num[$np+1][1]}][$j][2];
        $delta= $dihedb_nongo[$res{$resnam_num[$np][$maxres[$np]]}][$res{$resnam_num[$np+1][1]}][$j][3];
        printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$maxres[$np]-1][1],$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$newnam[$np+1][2][1],$kd,$period,$delta;
      }
    }
    # Ca, Ca+1, Ca+2, Ca+3
    if($nca[$np]>=1)
    {
      for($j=1;$j<=4;$j++)
      {
        $kd=    $dihedb_nongo[$res{$resnam_num[$np+1][1]}][$res{$resnam_num[$np+1][2]}][$j][1];
        $period=$dihedb_nongo[$res{$resnam_num[$np+1][1]}][$res{$resnam_num[$np+1][2]}][$j][2];
        $delta= $dihedb_nongo[$res{$resnam_num[$np+1][1]}][$res{$resnam_num[$np+1][2]}][$j][3];
        printf OUT "%-5s %-5s %-5s %-7s%-10.6f%-3s%-10.5f\n",$newnam[$np][$maxres[$np]][1],$newnam[$np+1][1][1],$newnam[$np+1][2][1],$newnam[$np+1][3][1],$kd,$period,$delta;
      }
    }
  }
}

print OUT "\nIMPHI\n";
print OUT "! sidechain improper dihedrals to maintain chirality\n";
if($casm == 1)
{
  for($np=0;$np<=$#pdb;$np++)
  {
    for($i=2;$i<=$maxres[$np]-1;$i++)
    {
      if($seg[$np][$i] == $seg[$np][$i-1] and $seg[$np][$i] == $seg[$np][$i+1])
      {
        $delta = $diheds[$np][$i]+180;
        $kd = 20*abs($avg_mj);
        if($ns[$np][$i] > 0) {printf OUT "%-5s %-5s %-5s %-7s%-10.6f%3s%10.5f\n",$newnam[$np][$i][1],$newnam[$np][$i-1][1],$newnam[$np][$i+1][1],$newnam[$np][$i][2],$kd,1,$delta }
      }
    }

    # Add improper dihedral energy terms at the last 1 and first 1 residues of PDB files
    if($#pdb>0 and $np!=$#pdb)
    {
      if($nca[$np]>=2)
      {
        # Ca,max - Ca,max-1 - Ca,1 - Cb,max
        $delta = $diheds[$np][$i]+180;
        $kd = 20*abs($avg_mj);
        if($ns[$np][$maxres[$np]] > 0) {printf OUT "%-5s %-5s %-5s %-7s%-10.6f%3s%10.5f\n",$newnam[$np][$maxres[$np]][1],$newnam[$np][$maxres[$np]-1][1],$newnam[$np+1][1][1],$newnam[$np][$maxres[$np]][2],$kd,1,$delta }      
      }
      if($nca[$np]>=2)
      {
        # Ca,i - Ca,i-1 - Ca,i+1 - Cb,i
        $delta = $diheds[$np][$i]+180;
        $kd = 20*abs($avg_mj);
        if($ns[$np+1][1] > 0) {printf OUT "%-5s %-5s %-5s %-7s%-10.6f%3s%10.5f\n",$newnam[$np+1][1][1],$newnam[$np][$maxres[$np]][1],$newnam[$np+1][2][1],$newnam[$np+1][1][2],$kd,1,$delta }      
      }
    }               
  }
}
print OUT "\n";
## nonbonded section
# I'm not sure whether to use NBXMOD 3 or 4
print OUT "NONBONDED NBXMOD 3 ATOM CDIEL SWITCH VATOM VDISTANCE VSWITCH -\n";
print OUT "CUTNB 32 CTOFNB 20 CTONNB 18 EPS 78.5 WMIN 1.5 E14FAC 1.0\n";
print OUT "!atom           e_min   r_min/2\n";
# if using the C-alpha only model do some preprocessing to determine the collision
# diameter of non-native interactions according to the Karanacolis-Brooks
# algorightm
for($np=0;$np<=$#pdb;$np++)
{
  if($casm != 1)
  {
    if(uc($pot[$np]) =~ m/GENERIC/)
    {
      for($i=1;$i<=$nca[$np];$i++)
      {
        $sigmin[$np][$i]=2*$rvdw{$resnam_num[$np][$i]};
      }
    }
    else
    {
      # determine the collision diameter
      for($i=1;$i<=$nca[$np]-3;$i++)
      { 
        $sigmin[$np][$i]=1000000;
        for($j=$i+3;$j<=$nca[$np];$j++)
        {
          if($native[$np][$i][$j] != 1)
          {
            if($distances[$np][$i][$j] < $sigmin[$np][$i]) {$sigmin[$np][$i]=$distances[$np][$i][$j]}
          }
        }
        # Check that the collision diameter is not outrageously large 
        if($sigmin[$np][$i] > 3*$rvdw{$resnam_num[$np][$i]}) 
        {
          print "Adjusting collision diameter for residue $i ($resnam_num[$np][$i]) from $sigmin[$np][$i] to 3*$rvdw{$resnam_num[$np][$i]}\n";  
          $sigmin[$np][$i]=3*$rvdw{$resnam_num[$np][$i]};
        }
      }
      # get last three collision diameters
      for($i=$nca[$np];$i>$nca[$np]-3;$i--)
      {
        $sigmin[$np][$i]=1000000;
        for($j=$i-3;$j>=1;$j--)
        {
          if($native[$np][$i][$j] != 1)
          {
            if($distances[$np][$i][$j] < $sigmin[$np][$i]) {$sigmin[$np][$i]=$distances[$np][$i][$j]}
          }
        }
        # Check that the collision diameter is not outrageously large 
        if($sigmin[$np][$i] > 3*$rvdw{$resnam_num[$np][$i]}) 
        {
          print "Adjusting collision diameter for residue $i ($resnam_num[$np][$i]) from $sigmin[$np][$i] to 3*$rvdw{$resnam_num[$np][$i]}\n";
          $sigmin[$np][$i]=3*$rvdw{$resnam_num[$np][$i]};
        }
      }
    }
  } 

  for($i=1;$i<=$nca[$np];$i++)
  {
    if($casm==1) 
    {
      $eps2 = -1e-12; #!!!! SYSTem dependent !!!!!!!!
      $rmin2 = 20.0;
      printf OUT "%-9s%-5.1f$eps2    %-10.6f\n",$newnam[$np][$i][1],0.0,$rmin2;
      $t1 = 1;
      $t2 = ($t1*(2*$rvdw{$resnam_num[$np][$i]}*2**(1/6))**12/(1e-12))**(1/12);
      $temp = $fnn*$t2/2;
      if($ns[$np][$i] > 0) { printf OUT "%-9s%-5.1f$eps2    %-10.6f\n",$newnam[$np][$i][2],0.0,$temp}
    }
    else
    {
      $eps2 = -0.000132;
      $rmin2 = $sigmin[$np][$i]*2**(1/6)/2;
      $temp = $fnn*$rmin2;
      printf OUT "%-9s%-5.1f$eps2    %-10.6f\n",$newnam[$np][$i][1],0.0,$temp;
    }
  }
}
print OUT "\n";
## NBFIX section
print OUT "NBFIX\n";
### native side-chain pairs and backbone Hbonding
for($np=0;$np<=$#pdb;$np++)
{
  if($casm==1) 
  {
    print OUT "! b-b due to Hbonding\n";
    for($i=1;$i<=$nhbond[$np];$i++)
    {
      if($ndomain[$np] == 0) # If interface not defined
      {
        printf OUT "%-8s%-11s%-11.6f%-11.6f\n",$newnam[$np][$hbond_at[$np][$i][1]][1],$newnam[$np][$hbond_at[$np][$i][2]][1],-$hbond_ene[$np][$i],$dhb[$np][$i];
        $totene_bb += $hbond_ene[$np][$i];
      }
      else # If interface is defined
      {
        if($id_dom[$np][$hbond_at[$np][$i][1]] == $id_dom[$np][$hbond_at[$np][$i][2]]) # residues within same domain
        {
          printf OUT "%-8s%-11s%-11.6f%-11.6f\n",$newnam[$np][$hbond_at[$np][$i][1]][1],$newnam[$np][$hbond_at[$np][$i][2]][1],-$hbond_ene[$np][$i],$dhb[$np][$i];
          $totene_bb += $hbond_ene[$np][$i];
        }
        else # residues in differt domain
        {
          $rescale = $int_scfac;
          printf OUT "%-8s%-11s%-11.6f%-11.6f ! Interface between domains $id_dom[$np][$hbond_at[$np][$i][1]], $id_dom[$np][$hbond_at[$np][$i][2]]\n",$newnam[$np][$hbond_at[$np][$i][1]][1],$newnam[$np][$hbond_at[$np][$i][2]][1],-$rescale,$dhb[$np][$i];
          $totene_bb += $hbond_ene[$np][$i];   
        }
      }
    }

    print OUT "! native side-chain interactions\n";
    if(uc($pot[$np]) =~ m/GENERIC/) # C-alpha - side chain model Generic non-bond interactions
    {
      print OUT "!Generic interactions between unstructured portions of this protein, from PDB $pdb[$np]\n";
      # Print out NBFIX energy values
      for($i=1;$i<=$nca[$np]-3;$i++)
      {
        for($j=$i+3;$j<=$nca[$np];$j++)
        {
          $temp=$rvdw{$resnam_num[$np][$i]} + $rvdw{$resnam_num[$np][$j]};
          $ene=(0.3/10)*$eps[$res{$resnam_num[$np][$i]}][$res{$resnam_num[$np][$j]}];
          printf OUT "%-8s%-11s%-11.6f%-11.6f\n",$newnam[$np][$i][1],$newnam[$np][$j][1],-$ene,$temp;
        }
      }
    }
    else # Go non-bond interactions 
    {
      for($i=1;$i<=$nnatsc[$np];$i++)
      {
        $temp = $natdist[$np][$i];
        if($eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}] == 0)
        {
          print "$i, $nsc1[$np][$i]-$nsc2[$np][$i], $resnam_num[$np][$nsc1[$np][$i]]-$resnam_num[$np][$nsc2[$np][$i]],$res{$resnam_num[$np][$nsc1[$np][$i]]}-$res{$resnam_num[$np][$nsc2[$np][$i]]}\n";
          die "ERROR 1: Well depth equal to zero!!! $eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}]\n";
        }
        if($ndomain[$np] == 0) # If interface not defined
        {
          printf OUT "%-8s%-11s%-13.6f%-11.6f\n",$newnam[$np][$nsc1[$np][$i]][2],$newnam[$np][$nsc2[$np][$i]][2],-$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}],$temp;
          $totene_sc += $eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}];
        }
        else # If interface is defined
        {
          if($id_dom[$np][$nsc1[$np][$i]] == $id_dom[$np][$nsc2[$np][$i]]) # residues within same domain
          {
            printf OUT "%-8s%-11s%-13.6f%-11.6f\n",$newnam[$np][$nsc1[$np][$i]][2],$newnam[$np][$nsc2[$np][$i]][2],-$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}],$temp;
            $totene_sc += $eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}];  
          }
          else # residues in differt domain
          {
            $rescale = $eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}]/$nscal;
            $rescale2 = $int_scfac*$rescale;
            printf OUT "%-8s%-11s%-13.6f%-11.6f ! Interface between domains $id_dom[$np][$nsc1[$np][$i]], $id_dom[$np][$nsc2[$np][$i]]\n",$newnam[$np][$nsc1[$np][$i]][2],$newnam[$np][$nsc2[$np][$i]][2],-$rescale2,$temp;
            $totene_sc += $rescale2;
          }
        }
      }

      print OUT "! backbone-sidechain interactions\n";
      for($i=1;$i<=$nbsc[$np];$i++)
      {
        if($ndomain[$np] == 0) # If interface not defined
        {
          printf OUT "%-8s%-11s%-11.6f%-11.6f\n",$newnam[$np][$natbsc[$np][$i][1]][1],$newnam[$np][$natbsc[$np][$i][2]][2],-$ene_bsc,$natdist_bsc[$np][$i];
          $totene_bsc += $ene_bsc;
        }
        else # If interface is defined
        {
          if($id_dom[$np][$natbsc[$np][$i][1]] == $id_dom[$np][$natbsc[$np][$i][2]]) # residues within same domain
          { 
            printf OUT "%-8s%-11s%-11.6f%-11.6f\n",$newnam[$np][$natbsc[$np][$i][1]][1],$newnam[$np][$natbsc[$np][$i][2]][2],-$ene_bsc,$natdist_bsc[$np][$i];
            $totene_bsc += $ene_bsc;
          }
          else # residues in differt domain
          {
            $rescale = $int_scfac*$ene_bsc;
            printf OUT "%-8s%-11s%-11.6f%-11.6f ! Interface between domains $id_dom[$np][$natbsc[$np][$i][1]] and $id_dom[$np][$natbsc[$np][$i][2]]\n",$newnam[$np][$natbsc[$np][$i][1]][1],$newnam[$np][$natbsc[$np][$i][2]][2],-$rescale,$natdist_bsc[$np][$i];
            $totene_bsc += $rescale;
          }
        }
      }
    }

    print OUT "\n";
    print OUT "! $totene_bb, $totene_sc, $totene_bsc\n";
  }
  else
  {
    if(uc($pot[$np]) !~ m/GENERIC/) # C-alpha model
    {
      print OUT "! b-b due to Hbonding plus native side-chain interactions plus backbone-sidechain interactions, from PDB $pdb[$np]\n";
      # Add up non-bonded energies
      @energy=();
      # loop through hydrogen bonds
      for($i=1;$i<=$nhbond[$np];$i++)
      {
        $distance[$np][$hbond_at[$np][$i][1]][$hbond_at[$np][$i][2]] = $dhb[$np][$i];
        if($ndomain[$np] == 0) # If interface not defined
        {
          $energy[$np][$hbond_at[$np][$i][1]][$hbond_at[$np][$i][2]]+=$hbond_ene[$np][$i];
          $energy[$np][$hbond_at[$np][$i][2]][$hbond_at[$np][$i][1]]+=$hbond_ene[$np][$i];
        }
        else #add by Yang
        {
          $energy[$np][$hbond_at[$np][$i][1]][$hbond_at[$np][$i][2]]+=$hbond_ene[$np][$i];
          $energy[$np][$hbond_at[$np][$i][2]][$hbond_at[$np][$i][1]]+=$hbond_ene[$np][$i];
        }
      }
      # loop through sc-sc interactions
      for($i=1;$i<=$nnatsc[$np];$i++)
      {
        $distance[$np][$nsc1[$np][$i]][$nsc2[$np][$i]] = $natdist[$np][$i];
        if($eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}] == 0)
        {
          print "$i, $nsc1[$np][$i]-$nsc2[$np][$i], $resnam_num[$np][$nsc1[$np][$i]]-$resnam_num[$np][$nsc2[$np][$i]],$res{$resnam_num[$np][$nsc1[$np][$i]]}-$res{$resnam_num[$np][$nsc2[$np][$i]]}\n";
          die "ERROR 2: Well depth equal to zero!!! $eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}]\n";
        }
        if($ndomain[$np] == 0) # If interface not defined
        {
          $energy[$np][$nsc1[$np][$i]][$nsc2[$np][$i]]+=$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}];
          $energy[$np][$nsc2[$np][$i]][$nsc1[$np][$i]]+=$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}];
        }
        else #add by Yang
        {
          if($id_dom[$np][$nsc1[$np][$i]] == $id_dom[$np][$nsc2[$np][$i]]) # residues within same domain
          {
            $energy[$np][$nsc1[$np][$i]][$nsc2[$np][$i]]+=$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}]/$nscal*$dom_nscal[$id_dom[$np][$nsc2[$np][$i]]-1];
            $energy[$np][$nsc2[$np][$i]][$nsc1[$np][$i]]+=$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}]/$nscal*$dom_nscal[$id_dom[$np][$nsc2[$np][$i]]-1];
          }
          else # residues in different domain
          {
            my ($interface_id_1, $interface_id_2);
            if($id_dom[$np][$nsc1[$np][$i]] < $id_dom[$np][$nsc2[$np][$i]])
            {
              $interface_id_1 = $id_dom[$np][$nsc1[$np][$i]];
              $interface_id_2 = $id_dom[$np][$nsc2[$np][$i]];
            }
            else
            {
              $interface_id_1 = $id_dom[$np][$nsc2[$np][$i]];
              $interface_id_2 = $id_dom[$np][$nsc1[$np][$i]];
            }
            my $num = $interface_id_2 - $interface_id_1 + ($interface_id_1 - 1) * $ndomain[$np] - $interface_id_1 * ($interface_id_1 - 1) / 2;
            my $i_nscal = $dom_nscal[$ndomain[$np] + $num - 1];

            $energy[$np][$nsc1[$np][$i]][$nsc2[$np][$i]]+=$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}]/$nscal*$i_nscal;
            $energy[$np][$nsc2[$np][$i]][$nsc1[$np][$i]]+=$eps[$res{$resnam_num[$np][$nsc1[$np][$i]]}][$res{$resnam_num[$np][$nsc2[$np][$i]]}]/$nscal*$i_nscal;
          }
        }
      }
      # loop through bb-sc interactions 
      for($i=1;$i<=$nbsc[$np];$i++)
      {
        $distance[$np][$natbsc[$np][$i][1]][$natbsc[$np][$i][2]]=$natdist_bsc[$np][$i];
        $distance[$np][$natbsc[$np][$i][2]][$natbsc[$np][$i][1]]=$natdist_bsc[$np][$i];
        if($ndomain[$np] == 0) # If interface not defined
        {
          $energy[$np][$natbsc[$np][$i][1]][$natbsc[$np][$i][2]]+=$ene_bsc;
          $energy[$np][$natbsc[$np][$i][2]][$natbsc[$np][$i][1]]+=$ene_bsc;
        }
        else #add by Yang
        {
          $energy[$np][$natbsc[$np][$i][1]][$natbsc[$np][$i][2]]+=$ene_bsc;
          $energy[$np][$natbsc[$np][$i][2]][$natbsc[$np][$i][1]]+=$ene_bsc;
        }
      }
      # Print out NBFIX energy values
      for($i=1;$i<=$nca[$np]-2;$i++)
      {
        for($j=$i+2;$j<=$nca[$np];$j++)
        {
          my $comment = ""; #add by Yang
          if($id_dom[$np][$i] eq $id_dom[$np][$j])
          {
            $comment = "! in Domain ".$id_dom[$np][$i];
          }
          else
          {
            $comment = "! in Interface ".$id_dom[$np][$i]." | ".$id_dom[$np][$j];
          }
          if($energy[$np][$i][$j]>0) {printf OUT "%-8s%-11s%-11.6f%-11.6f%s\n",$newnam[$np][$i][1],$newnam[$np][$j][1],-$energy[$np][$i][$j],$distance[$np][$i][$j],$comment}
        }
      }
    }
    else
    {
      print OUT "!Generic interactions between unstructured portions of this protein, from PDB $pdb[$np]\n";
      # Print out NBFIX energy values
      for($i=1;$i<=$nca[$np]-3;$i++)
      {
        for($j=$i+3;$j<=$nca[$np];$j++)
        {
          $temp=$rvdw{$resnam_num[$np][$i]} + $rvdw{$resnam_num[$np][$j]};
          $ene=(0.3/10)*$eps[$res{$resnam_num[$np][$i]}][$res{$resnam_num[$np][$j]}];
          printf OUT "%-8s%-11s%-11.6f%-11.6f\n",$newnam[$np][$i][1],$newnam[$np][$j][1],-$ene,$temp;
        }
      }
    }
  }
}
print OUT "\n";
print OUT "END\n";
close OUT;

$out="$prefix_"."$interact_scale.dat";
open(OUT,">$out");
# If the user wants the generic interaction between two different pdb files
if($interact_scale ne "undefined")
{
  for($np=0;$np<=$#pdb-1;$np++)
  {
    for($np2=$np+1;$np2<=$#pdb;$np2++)
    {
      if($cascm==1) {}
      else 
      {
        printf OUT "! The interaction energies between PDB files $pdb[$np] and $pdb[$np2] with interact_scale=$interact_scale\n";
        for($i=1;$i<=$nca[$np];$i++)
        {
          for($j=1;$j<=$nca[$np2];$j++)
          {
            $temp=$rvdw{$resnam_num[$np][$i]} + $rvdw{$resnam_num[$np2][$j]};
            $ene=$interact_scale*$eps[$res{$resnam_num[$np][$i]}][$res{$resnam_num[$np2][$j]}];
            printf OUT "%-8s%-11s%-11.6f%-11.6f\n",$newnam[$np][$i][1],$newnam[$np2][$j][1],-$ene,$temp;
          }
        }
      } 
    }
  }
}
print OUT "\nEND\n";
close OUT;

# (10) Write out PSF script and file if CHARMM path is defined
&write_psf;
#if($charmm ne "none") {system("$charmm < create_psf.inp")}

# (11) Minimize the starting structure to have reasonable bond lengths, etc.
#&write_minimization_script;
#if($charmm ne "none") {system("$charmm < minimize.inp")}

##################### SUBROUTINES ######################
sub cleanws
{
  my @a = @_;
  my($i);

  for($i=0;$i<=$#a;$i++)
  {
    # Remove white spaces at begining of string
    $a[$i] =~ s/^\s+//;
    # Remove white spaces at end of string
    $a[$i] =~ s/\s+$//;
  }
  return @a;
}

sub Determine_CM
{
  my ($strt,$end) = @_;
  my $m,$nat;
  # Determine Center of Mass
  @Rcm = 0;
  $SumMass = 0;
  $nat = 0;
  for($m=$strt;$m<=$end;$m++)
  {
    if($atomnam[$m] ne 'N' and $atomnam[$m] ne 'CA' and $atomnam[$m] ne 'C' and $atomnam[$m] ne 'O' and $atomnam[$m] ne 'OT1' and $atomnam[$m] ne 'OT2')
    {  
      $nat++;
      $Rcm[0] = $Rcm[0] + $Mass{$atom[$m]}*$x[$m];
      $Rcm[1] = $Rcm[1] + $Mass{$atom[$m]}*$y[$m];
      $Rcm[2] = $Rcm[2] + $Mass{$atom[$m]}*$z[$m];
      $SumMass = $SumMass + $Mass{$atom[$m]};
    }
  }
  $Rcm[0] = $Rcm[0]/$SumMass;
  $Rcm[1] = $Rcm[1]/$SumMass;
  $Rcm[2] = $Rcm[2]/$SumMass;
  if($refNscat{uc($resnam[$np][$end])} != $nat) {die "ERROR 2: number of heavy atoms in sidechain $resnum[$end] ($resnam[$np][$end]) incorrect in CM subroutine: $refNscat{uc($resnam[$np][$end])} != $nat\n"}
}

sub Stride
{
  @D = ();
  if( -e $out )
  {
    if( -s $out  )
    {
      $cmd = "stride $out -h |"; # calc main-chain H-bonds
      open(D, "$cmd") or die "Can't run program: $!\n";
      @D = <D>;
      close(D);
    }
    else { die "Error: File size is zero for $out\n" }
  }
  else { die "Error2: File $out does not exist\n" }
}

sub distance
{
  my ($x1,$y1,$z1,$x2,$y2,$z2) = @_;
  return (($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2)**0.5;
}

sub Angle
{
  my ($xp1,$yp1,$zp1,$xc,$yc,$zc,$xp2,$yp2,$zp2) = @_;

  @Dr1 = ($xc - $xp2,$yc - $yp2,$zc - $zp2);
  $CP1dist = ($Dr1[0]**2+$Dr1[1]**2+$Dr1[2]**2)**0.5;

  @Dr2 = ($xc - $xp1,$yc - $yp1,$zc - $zp1);
  $CP2dist = ($Dr2[0]**2+$Dr2[1]**2+$Dr2[2]**2)**0.5;

  @Dr1 = ($Dr1[0]/$CP1dist,$Dr1[1]/$CP1dist,$Dr1[2]/$CP1dist);
  @Dr2 = ($Dr2[0]/$CP2dist,$Dr2[1]/$CP2dist,$Dr2[2]/$CP2dist);

  return acos( $Dr1[0]*$Dr2[0] + $Dr1[1]*$Dr2[1] + $Dr1[2]*$Dr2[2] )*180/pi;
}

sub visualwritepdb
{
  my $i,$j,$k;

  $seg_name = "A";

  $count = 0;
  open(OUT,">temp_visual.pdb");
  for($i=1;$i<=$nca;$i++)
  {
    $count++;
    $X = sprintf("%7.3f",$xb[$i]);
    $Y = sprintf("%7.3f",$yb[$i]);
    $Z = sprintf("%7.3f",$zb[$i]);
    $Xf = sprintf("%7.7s",$X);
    $Yf = sprintf("%7.7s",$Y);
    $Zf = sprintf("%7.7s",$Z);
    $n1 = "A" ;
    $atom_name = $n1;
    printf OUT "%3s %6s %4s %4s %5s %10s %7s %7s %5.2f %5.2f %8s\n",ATOM,$count,$n1,$atom_name,$count,$Xf,$Yf,$Zf,1,0,$seg_name;

    if($gly[$i] == 0)
    { 
      $X = sprintf("%7.3f",$xs[$i]);
      $Y = sprintf("%7.3f",$ys[$i]);
      $Z = sprintf("%7.3f",$zs[$i]);
      $Xf = sprintf("%7.7s",$X);
      $Yf = sprintf("%7.7s",$Y);
      $Zf = sprintf("%7.7s",$Z);
      $n1 = "B";
      $atom_name = $n1;
      printf OUT "%4s %6s %4s %4s %5s %10s %7s %7s %5.2f %5.2f %8s\n",ATOM,$count,$n1,$atom_name,$count,$Xf,$Yf,$Zf,1,0,$seg_name; 
    }
  }
  close OUT;
}

sub writecrd
{
  my $i,$j,$k,$out;

  $seg_name = "A";
  $count = 0;
  if($prefix eq "none")
  {
    $prefix = "$fileBaseName";
  }

  if($casm==1) {$out = "$prefix" . "_ca-cb.cor"}
  else {$out = "$prefix" . "_ca.cor"}
  open(OUT,">$out");
  print OUT "*\n";

  for($np=0;$np<=$#pdb;$np++) {$total_nres+=$maxres[$np]}

  if($casm==1) {$temp = 2*$total_nres - $ngly}
  else {$temp=$total_nres}

  print OUT "  $temp\n";
  if($casm==1) {$out = "$prefix" . "_ca-cb.seq"}
  else {$out = "$prefix" . "_ca.seq"}
  open(OUT2,">$out");

  for($np=0;$np<=$#pdb;$np++)
  {
    for($i=1;$i<=$nca[$np];$i++)
    {
      $count++;
      $X = sprintf("%7.3f",$xb2[$np][$i]);
      $Y = sprintf("%7.3f",$yb2[$np][$i]);
      $Z = sprintf("%7.3f",$zb2[$np][$i]);
      $Xf = sprintf("%7.7s",$X);
      $Yf = sprintf("%7.7s",$Y);
      $Zf = sprintf("%7.7s",$Z);
      printf OUT "%5s%5s %-4s %-3s%11.5f%10.5f%10.5f%2s%5s%13.5f\n",$count,$count,$resnam[$np][$i],A,$Xf,$Yf,$Zf,$seg_name,1,0;
      printf OUT2 " $resnam[$np][$i]";
      if($ns[$np][$i] > 0 and $casm==1)
      {
        $X = sprintf("%7.3f",$xs2[$np][$i]);
        $Y = sprintf("%7.3f",$ys2[$np][$i]);
        $Z = sprintf("%7.3f",$zs2[$np][$i]);
        $Xf = sprintf("%7.7s",$X);
        $Yf = sprintf("%7.7s",$Y);
        $Zf = sprintf("%7.7s",$Z);
        printf OUT "%5s%5s %-4s %-3s%11.5f%10.5f%10.5f%2s%5s%13.5f\n",$count,$count,$resnam[$np][$i],B,$Xf,$Yf,$Zf,$seg_name,1,0;
      }
    }
  }
  close OUT;
  close OUT2;
}

sub write_psf
{
  my $i;
  my $out = "create_psf.inp";
  if($casm == 1) {$out2="$prefix"."_ca-cb"}
  else {$out2="$prefix"."_ca"}
  open(OUT,">$out");
  $total_nres=0;
  for($np=0;$np<=$#pdb;$np++) {$total_nres+=$maxres[$np]}
  print OUT "* CREATE PSF
  *   by Edward P. O'Brien Jr
  *

  prnlev 5

  ! file handle for input/output
  set out $out2
  ! segment name
  set seg A

  ! read parameter and topology files
  open unit 10 read form name $top_file
  read rtf unit 10 card
  close unit 10

  open unit 10 read form name $prm_file
  read param unit 10 card
  close unit 10

  READ SEQUENCE CARDS
  * CG model
  *
  $total_nres\n";

  $ncnt=0;
  for($np=0;$np<=$#pdb;$np++)
  {
    for($i=1;$i<=$nca[$np];$i++)
    {
      $ncnt++;
      print OUT "$resnam[$np][$i] ";
      if($ncnt % 25 == 0 ) {print OUT "\n"}
    }
  }
  print OUT "\n";
  print OUT '
  ! This creates the psf parameters
  GENERATE @seg FIRS none LAST none SETUP
  ! Build up coordinate from IC parameter
  IC PARA

  ! This writes the psf parameters
  open write unit 10 card name @out.psf
  write psf card unit 10
  close unit 10

  STOP';
}


sub write_minimization_script
{
  my $i;
  my $out = "minimize.inp";
  if($casm == 1) 
  {
    $out2="$prefix"."_ca-cb_mini.cor";
    $strt="$prefix"."_ca-cb.cor";
    $psf="$prefix"."_ca-cb.psf";
  }
  else 
  {
    $out2="$prefix"."_ca_mini.cor";
    $strt="$prefix"."_ca.cor";
    $psf="$prefix"."_ca.psf";
  }
  open(OUT,">$out");
  $total_nres=0;
  for($np=0;$np<=$#pdb;$np++) {$total_nres+=$maxres[$np]}
  print OUT "* MINIMIZE STRUCTURE
  *   by Edward P. O'Brien Jr
  *

  prnlev 5

  read rtf card name $top_file
  read param card name $prm_file
  read psf card name $psf
  read coor card name $strt

  mini sd nstep 10000

  write coor card name $out2

  STOP";
}

sub get_Torsion
{
  my($r1,$r2,$r3,$r4,$d) = @_;
  my @vec0,@vec1,@vec2,@vec3;
  # define vectors
  if(lc($d) eq 'b')
  {
    @vec0 = ($xb2[$np][$r1],$yb2[$np][$r1],$zb2[$np][$r1]);
    @vec1 = ($xb2[$np][$r2],$yb2[$np][$r2],$zb2[$np][$r2]);
    @vec2 = ($xb2[$np][$r3],$yb2[$np][$r3],$zb2[$np][$r3]);
    @vec3 = ($xb2[$np][$r4],$yb2[$np][$r4],$zb2[$np][$r4]);
  }
  elsif(lc($d) eq 's')
  {
    @vec0 = ($xb2[$np][$r1],$yb2[$np][$r1],$zb2[$np][$r1]);
    @vec1 = ($xb2[$np][$r2],$yb2[$np][$r2],$zb2[$np][$r2]);
    @vec2 = ($xb2[$np][$r3],$yb2[$np][$r3],$zb2[$np][$r3]);
    @vec3 = ($xs2[$np][$r4],$ys2[$np][$r4],$zs2[$np][$r4]);
  }
  else {die "ERROR: call to get_Torsion subroutine screwed up!\n"}

  # get difference vectors
  my @dv0 = VecSub(@vec0,@vec1);
  my @dv1 = VecSub(@vec1,@vec2);
  my @dv2 = VecSub(@vec2,@vec3);

  # get cross product
  @crossproduct_result0 = CrossProduct(@dv0,@dv1);
  @crossproduct_result1 = CrossProduct(@dv1,@dv2);

  # get magnitude of vectors
  $magnitude_result0 = Magnitude(@crossproduct_result0);
  $magnitude_result1 = Magnitude(@crossproduct_result1);

  # normalize vectors
  my @n0 = ($crossproduct_result0[0]/$magnitude_result0,$crossproduct_result0[1]/$magnitude_result0,$crossproduct_result0[2]/$magnitude_result0);
  my @n1 = ($crossproduct_result1[0]/$magnitude_result1,$crossproduct_result1[1]/$magnitude_result1,$crossproduct_result1[2]/$magnitude_result1);

  # Method for assigning + and - angles
  @crossproduct_result3 = CrossProduct(@n0,@n1);

  $dotproduct_result0 = DotProduct(@crossproduct_result3,@dv1);

  if( $dotproduct_result0 < 0 ) { $prefactor = 1 }
  else { $prefactor = -1 }

  # do the dot product and calc angle
  $dotproduct_result1 = DotProduct(@n0,@n1);
  return $prefactor*acos($dotproduct_result1)*180/pi;
}

sub DotProduct
#returns the dot product of a vector
{
#        shift @_;
my ($x1,$y1,$z1,$x2,$y2,$z2) = @_[0..5];
return ($x1*$x2)+($y1*$y2)+($z1*$z2);
}

sub CrossProduct
#returns the cross product of two vectors - ORDER MATTERS!!!!!
{
  my ($x,$y,$z);
#        shift @_;
my ($a1,$a2,$a3,$b1,$b2,$b3) = @_[0..5];
$x = (($a2*$b3) - ($a3*$b2));
$y = (($a3*$b1) - ($a1*$b3));
$z = (($a1*$b2) - ($a2*$b1));
return ($x,$y,$z);
}

sub Magnitude
#returns the magnitude of a vector
{
#        shift @_;
my ($x,$y,$z) = @_[0..2];
return (($x*$x)+($y*$y)+($z*$z))**(0.5)
}

sub VecSub
#subtracts vector b from vector a
{
#        shift @_;
my ($a1,$a2,$a3,$b1,$b2,$b3) = @_[0..5];
return ($a1-$b1,$a2-$b2,$a3-$b3);
}



