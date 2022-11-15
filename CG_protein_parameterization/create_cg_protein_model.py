#!/usr/bin/env python3
try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
from sys import stdout, exit, stderr
import getopt, os, time, random, math, traceback, io, sys, string
import parmed as pmd
import numpy as np

sys.setrecursionlimit(int(1e6))

usage = '\nUsage: create_cg_protein_model.py\n' \
        '       --ctrlfile | -f <model.ctrl> Control file for creating cg protein model\n'\
        '       [--help | -h] Print this information\n\n'\

######################## Data #########################
## Loop-up table for uniquely indentifying residues #
aa = ["GLY","ALA","VAL","LEU","ILE","MET","PHE","PRO","SER","THR","CYS","ASN","GLN","TYR","TRP","ASP","GLU","HIS","LYS","ARG"]
if len(aa) != 20:
    print('ERROR')
    sys.exit()
res2n = {}
n2res = {}
for i, a in enumerate(aa):
    res2n[a] = i
    n2res[i] = a

Mass = {"N": 14.0067,
        "H": 1.00794,
        "C": 12.011,
        "O": 15.9994,
        "S": 32.06,}

# number of heavy atoms in sidechains
refNscat = {"ALA": 1,
            "CYS": 2,
            "ASP": 4,
            "GLU": 5,
            "PHE": 7,
            "GLY": 0,
            "HIS": 6,
            "HSD": 6,
            "HSE": 6,
            "HSP": 6,
            "ILE": 4,
            "LYS": 5,
            "LEU": 4,
            "MET": 4,
            "ASN": 4,
            "PRO": 3,
            "GLN": 5,
            "ARG": 7,
            "SER": 2,
            "THR": 3,
            "VAL": 3,
            "TRP": 10,
            "TYR": 8}

# charges on side chains at pH 7
refcharge = {"ALA": 0.0,
             "CYS": 0.0,
             "ASP": -1.0,
             "GLU": -1.0,
             "PHE": 0.0,
             "GLY": 0.0,
             "HIS": 0.0,
             "HSD": 0.0,
             "HSE": 0.0,
             "HSP": 0.0,
             "ILE": 0.0,
             "LYS": 1.0,
             "LEU": 0.0,
             "MET": 0.0,
             "ASN": 0.0,
             "PRO": 0.0,
             "GLN": 0.0,
             "ARG": 1.0,
             "SER": 0.0,
             "THR": 0.0,
             "VAL": 0.0,
             "TRP": 0.0,
             "TYR": 0.0}

# Generic C_alpha side-chain center of mass distance
lbs_nongo = {"ASP": 2.46916481058687,
             "PRO": 1.87381801537346,
             "LYS": 3.49738414814426,
             "ILE": 2.25260184847053,
             "TRP": 3.58251993741888,
             "CYS": 2.06666004558289,
             "HSD": 3.15209719417679,
             "PHE": 3.38385541816659,
             "HSP": 3.15209719417679,
             "GLN": 3.08654121335,
             "SER": 1.89840600762153,
             "ASN": 2.46916481058687,
             "VAL": 1.93953811063784,
             "LEU": 2.56580983973678,
             "TYR": 3.38981664391425,
             "GLU": 3.07971386504681,
             "ARG": 3.39687572938579,
             "THR": 1.931721703272,
             "ALA": 1.51146031725997,
             "MET": 2.95389402456081,
             "HIS": 3.15209719417679,
             "HSE": 3.15209719417679}

improper_nongo = {"ASP": 14.655341300544,
                  "PRO": 26.763068425539,
                  "LYS": 12.765248692601,
                  "ILE": 13.5446902008313,
                  "TRP": 11.4483488626106,
                  "CYS": 20.0484470024042,
                  "HSD": 14.9962640689562,
                  "PHE": 10.9217771918902,
                  "HSP": 14.9962640689562,
                  "GLN": 17.3050853491068,
                  "SER": 20.1390130256255,
                  "ASN": 14.655341300544,
                  "VAL": 13.3216022614598,
                  "LEU": 11.8137180266206,
                  "TYR": 12.2715081962165,
                  "GLU": 15.4130821146834,
                  "ARG": 15.5451613009777,
                  "THR": 16.2956083930276,
                  "ALA": 16.8418866013662,
                  "MET": 12.7046284165739,
                  "HIS": 14.9962640689562,
                  "HSE": 14.9962640689562}

ang_sb_nongo = {"ASP": 120.380153696218,
                "PRO": 125.127927161651,
                "LYS": 119.523270610009,
                "ILE": 118.791108398805,
                "TRP": 130.018548241749,
                "CYS": 110.512719347428,
                "HSD": 116.815900172681,
                "PHE": 122.937540996701,
                "HSP": 116.815900172681,
                "GLN": 116.182123224059,
                "SER": 107.971234136647,
                "ASN": 120.380153696218,
                "VAL": 112.877421898116,
                "LEU": 123.32179171436,
                "TYR": 116.783314494739,
                "GLU": 116.659068554985,
                "ARG": 119.709740783191,
                "THR": 111.719883260793,
                "ALA": 108.623605160075,
                "MET": 116.636559053295,
                "HIS": 116.815900172681,
                "HSE": 116.815900172681}

ang_bs_nongo = {"ASP": 116.629356207687,
                "PRO": 79.4932105625367,
                "LYS": 119.779735484239,
                "ILE": 116.923861483529,
                "TRP": 100.858690902849,
                "CYS": 114.816253227757,
                "HSD": 115.848569293979,
                "PHE": 112.804608190743,
                "HSP": 115.848569293979,
                "GLN": 119.106753006548,
                "SER": 116.361829754186,
                "ASN": 116.629356207687,
                "VAL": 121.299281732077,
                "LEU": 117.587011217416,
                "TYR": 116.72484692836,
                "GLU": 119.507585037498,
                "ARG": 117.532816176021,
                "THR": 117.044133956143,
                "ALA": 120.747734648009,
                "MET": 123.234171432545,
                "HIS": 115.848569293979,
                "HSE": 115.848569293979}

# segment id relationships
alphabet = list(map(chr, range(ord('A'), ord('Z')+1)))
segid2num = {}
for nseg, letter in enumerate(alphabet):
    segid2num[letter] = nseg
    
# mass of amino acids
# UNSURE! about pro, arg, his and cys weights
aaSCmass = {"ALA": 71.000000,
            "CYS": 114.000000,
            "ASP": 114.000000,
            "GLU": 128.000000,
            "PHE": 147.000000,
            "GLY": 57.000000,
            "HIS": 114.000000,
            "HSD": 114.000000,
            "HSE": 114.000000,
            "HSP": 114.000000,
            "ILE": 113.000000,
            "LYS": 128.000000,
            "LEU": 113.000000,
            "MET": 131.000000,
            "ASN": 114.000000,
            "PRO": 114.000000,
            "GLN": 128.000000,
            "ARG": 114.000000,
            "SER": 87.000000,
            "THR": 101.000000,
            "VAL": 99.000000,
            "TRP": 186.000000,
            "TYR": 163.000000}

# vdw radius of sidechains
rvdw = {"ALA": 2.51958406732374,
        "CYS": 2.73823091624513,
        "ASP": 2.79030096923572,
        "GLU": 2.96332591119925,
        "PHE": 3.18235414984794,
        "GLY": 2.25450393833984,
        "HIS": 3.04273820988499,
        "HSD": 3.04273820988499,
        "HSE": 3.04273820988499,
        "HSP": 3.04273820988499,
        "ILE": 3.09345983013354,
        "LYS": 3.18235414984794,
        "LEU": 3.09345983013354,
        "MET": 3.09345983013354,
        "ASN": 2.84049696898525,
        "PRO": 2.78004241717965,
        "GLN": 3.00796101305807,
        "ARG": 3.28138980397453,
        "SER": 2.59265585208464,
        "THR": 2.81059478021734,
        "VAL": 2.92662460060742,
        "TRP": 3.38869998431408,
        "TYR": 3.22881842919248}

######################### Functions ###########################
# generate charmm .psf
def create_psf(struct, ca_list, name):
    # creat backbone bonds
    for i in range(len(ca_list)-1):
        segid_list = [ca_list[i+j].residue.segid for j in range(2)]
        segid_list = list(set(segid_list))
        if len(segid_list) == 1:
            struct.bonds.append(pmd.topologyobjects.Bond(ca_list[i], ca_list[i+1]))
    # creat backbone-sidechain bonds if exist
    for ca_atom in ca_list:
        if len(ca_atom.residue.atoms) > 1:
            b_bead = ca_atom.residue.atoms[1]
            struct.bonds.append(pmd.topologyobjects.Bond(ca_atom, b_bead))
    # create Angles
    for atm in struct.atoms:
        bond_list = atm.bond_partners
        if len(bond_list) > 1:
            for i in range(len(bond_list)-1):
                for j in range(i+1, len(bond_list)):
                    struct.angles.append(pmd.topologyobjects.Angle(bond_list[i], atm, bond_list[j]))
    # create Dihedrals
    for i in range(len(ca_list)-3):
        segid_list = [ca_list[i+j].residue.segid for j in range(4)]
        segid_list = list(set(segid_list))
        if len(segid_list) == 1:
            struct.dihedrals.append(pmd.topologyobjects.Dihedral(ca_list[i], ca_list[i+1], ca_list[i+2], ca_list[i+3]))
    # create Impropers
    for i in range(1, len(ca_list)-1):
        segid_list = [ca_list[i+j-1].residue.segid for j in range(3)]
        segid_list = list(set(segid_list))
        if len(segid_list) == 1 and len(ca_list[i].residue.atoms) > 1:
            b_bead = ca_list[i].residue.atoms[1]
            struct.impropers.append(pmd.topologyobjects.Improper(ca_list[i], ca_list[i-1], ca_list[i+1], b_bead))
    struct.save(name+'.psf', overwrite=True, vmd=False)
# END generate charmm .psf

# generate charmm .top
def Create_rtf(struct, out_name):
    global pdbfile, casm
    fo = open(out_name+'.top', 'w')
    if casm == 1:
        fo.write('* This CHARMM .top file describes a Ca-Cb Go model of %s\n*\n20 1\n'%pdbfile)
    else:
        fo.write('* This CHARMM .top file describes a Ca Go model of %s\n*\n20 1\n'%pdbfile)
    # MASS section
    fo.write('! backbone masses\n')
    for idx, atm in enumerate(struct.atoms):
        fo.write('MASS %-4s %-8s %.6f\n'%(str(idx+1), atm.type, atm.mass))
    fo.write('\n')
    fo.write('DECL +%s\n'%struct[0].name)
    fo.write('DECL -%s\n'%struct[0].name)
    fo.write('DECL #%s\n'%struct[0].name)
    # residue section
    for res in struct.residues:
        res_charge = 0
        for atm in res.atoms:
            res_charge += atm.charge
        fo.write('RESI %-6s %.1f\n'%(res.name, res_charge))
        fo.write('GROUP\n')
        for atm in res.atoms:
            fo.write('ATOM %s %-6s %.1f\n'%(atm.name, atm.type, atm.charge))
        if casm == 1 and len(res.atoms) != 1:
            fo.write("Bond %s %s  %s +%s\n"%(res.atoms[0].name, res.atoms[1].name, 
                                             res.atoms[0].name, res.atoms[0].name))
            fo.write("Angle -%s %s %s  %s %s +%s  -%s %s +%s\n"%(res.atoms[0].name, res.atoms[0].name, res.atoms[1].name, 
                                                                 res.atoms[1].name, res.atoms[0].name, res.atoms[0].name, 
                                                                 res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
            fo.write("DIHE -%s %s +%s #%s\n"%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
            fo.write("IMPH %s -%s +%s %s\n\n"%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name, res.atoms[1].name))
        else:
            fo.write('Bond %s +%s\n'%(res.atoms[0].name, res.atoms[0].name))
            fo.write('Angle -%s %s +%s\n'%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
            fo.write('DIHE -%s %s +%s #%s\n\n'%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
    # end section
    fo.write('END\n')
    fo.close()
# END generate charmm .top

def calc_distance(atom_1, atom_2):
    dist = ((atom_1.xx - atom_2.xx)**2 + (atom_1.xy - atom_2.xy)**2 + (atom_1.xz - atom_2.xz)**2)**0.5
    return dist

def cg_energy_minimization(cor, prefix, prm_file):
    temp = 310
    np = '1'
    timestep = 0.015*picoseconds
    fbsolu = 0.05/picosecond
    temp = temp*kelvin

    psf_pmd = pmd.charmm.CharmmPsfFile(prefix+'.psf')
    psf = CharmmPsfFile(prefix+'.psf')
    top = psf.topology
    os.system('parse_cg_cacb_prm.py -p '+prm_file+' -t '+prefix+'.top')
    name = prm_file.split('.prm')[0]
    forcefield = ForceField(name+'.xml')
    
    template_map = {}
    for chain in top.chains():
        for res in chain.residues():
            template_map[res] = res.name
                
    
    system = forcefield.createSystem(top, nonbondedCutoff=2.0*nanometer, 
                                     constraints=None, removeCMMotion=False, ignoreExternalBonds=True,
                                     residueTemplates=template_map)
    custom_nb_force = system.getForce(4)
    custom_nb_force.setUseSwitchingFunction(True)
    custom_nb_force.setSwitchingDistance(1.8*nanometer)
    custom_nb_force.setNonbondedMethod(custom_nb_force.CutoffNonPeriodic)
    
    # add position restraints
    force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    system.addForce(force)
    # END add position restraints
    
    # add position restraints for CA
    force = system.getForces()[-1]
    k = 100*kilocalorie/mole/angstrom**2
    for atm in top.atoms():
        if atm.name == 'A':
            force.addParticle(atm.index, (k, cor[atm.index][0], cor[atm.index][1], cor[atm.index][2]))
    
    integrator = LangevinIntegrator(temp, fbsolu, timestep)
    integrator.setConstraintTolerance(0.00001)
    # prepare simulation
    platform = Platform.getPlatformByName('CPU')
    properties = {'Threads': np}
    simulation = Simulation(top, system, integrator, platform, properties)
    simulation.context.setPositions(cor)
    simulation.context.setVelocitiesToTemperature(temp)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
    getEnergyDecomposition(stdout, simulation.context, system)
    print('   Potential energy before minimization: %.4f kcal/mol'%energy)
    simulation.minimizeEnergy(tolerance=0.1*kilocalories_per_mole)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
    getEnergyDecomposition(stdout, simulation.context, system)
    print('   Potential energy after minimization: %.4f kcal/mol'%energy)
    current_cor = simulation.context.getState(getPositions=True).getPositions()
    return current_cor

# remove bond constraints of 0 mass atoms
def rm_cons_0_mass(system):
    tag = 0
    while tag == 0 and system.getNumConstraints() != 0:
        for i in range(system.getNumConstraints()):
            con_i = system.getConstraintParameters(i)[0]
            con_j = system.getConstraintParameters(i)[1]
            mass_i = system.getParticleMass(con_i).value_in_unit(dalton)
            mass_j = system.getParticleMass(con_j).value_in_unit(dalton)
            if mass_i == 0 and mass_j == 0:
                system.removeConstraint(i)
                #print('Constraint %d is removed, range is %d'%(i, system.getNumConstraints()))
                tag = 0
                break
            elif mass_i == 0 or mass_j == 0:
                system.removeConstraint(i)
                #print('Constraint %d is removed, range is %d'%(i, system.getNumConstraints()))
                system.getForce(0).addBond(con_i, con_j, 3.81*angstroms, 50*kilocalories/mole/angstroms**2)
                tag = 0
                break
            else:
                tag = 1
# END remove bond constraints of 0 mass atoms

# energy decomposition 
def forcegroupify(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        f = str(type(force))
        s = f.split('\'')
        f = s[1]
        s = f.split('.')
        f = s[-1]
        forcegroups[i] = f
    return forcegroups

def getEnergyDecomposition(handle, context, system):
    forcegroups = forcegroupify(system)
    energies = {}
    for i, f in forcegroups.items():
        try:
            states = context.getState(getEnergy=True, groups={i})
        except ValueError as e:
            print(str(e))
            energies[i] = Quantity(np.nan, kilocalories/mole)
        else:
            energies[i] = states.getPotentialEnergy()
    results = energies
    handle.write('    Potential Energy:\n')
    for idd in energies.keys():
        handle.write('      %s: %.4f kcal/mol\n'%(forcegroups[idd], energies[idd].value_in_unit(kilocalories/mole))) 
    return results

##################################### MAIN #######################################
ctrlfile = ''

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:", ["ctrlfile="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-f", "--ctrlfile"):
        ctrlfile = arg

## Check dependency installation ##
if os.popen('stride 2>&1').readlines()[0].strip().endswith('command not found'):
    print('Error: Essential software "stride" is not installed.\nPlease install stride before coarse-graining.')
    sys.exit()

## BEGIN: defaults ##
pdbfile = ''
nscal_0 = '1'
nscal = float(nscal_0)
fnn_0 = '1'
fnn = float(fnn_0)
potential_name = "MJ"
casm = 1 # 1 = create calpha-scm, 0 = create calpha model only!
heav_cut = 4.5 # angstroms, Definition of cutoff for sidechain heavy atoms contact
domain_file = "None"
ca_prefix = "A" # C-alpha atom prefix in top, prm, psf files
sc_prefix = "B" # side chain prefix
### END: defaults ###

## BEGIN: Parse the control file for information##
if not os.path.exists(ctrlfile):
    print('Error: cannot find control file ' + ctrlfile + '.')
    sys.exit()

file_object = open(ctrlfile,'r')
try:
    for line in file_object:
        line = line.strip()
        if not line:
            # This is a blank line
            continue
        if line.startswith('#'):
            # This is a comment line
            continue
        if line.startswith('pdbfile'): # read in the PDB file name here.
            words = line.split('=')
            pdbfile = words[1].strip()
            continue
        if line.startswith('nscal'):
            words = line.split('=')
            nscal_0 = words[1].strip()
            nscal = float(nscal_0)
            continue
        if line.startswith('fnn'):
            words = line.split('=')
            fnn_0 = words[1].strip()
            fnn = float(fnn_0)
            continue
        if line.startswith('potential_name'): # read in the potential name
            words = line.split('=')
            potential_name = words[1].strip()
            continue
        if line.startswith('casm'):
            words = line.split('=')
            casm = int(words[1].strip())
            continue
        if line.startswith('domain_file'):
            words = line.split('=')
            domain_file = words[1].strip()
            continue
        if line.startswith('sc_prefix'):
            words = line.split('=')
            sc_prefix = words[1].strip()
            continue
        if line.startswith('ca_prefix'):
            words = line.split('=')
            ca_prefix = words[1].strip()
            continue
finally:
     file_object.close()

print("\n##########################################")
print("# Build CG Protein Model: Python version #")
print("#   Yang Jiang & Edward P. O'Brien Jr.   #")
print("#            Dept. of Chemistry          #")
print("#          Penn State University         #")
print("##########################################\n")

print("Configuration:")
print("pdbfile = %s"%pdbfile)
print("casm = %d"%casm)
print("nscal = %.4f"%nscal)
print("fnn = %.4f"%fnn)
print("potential_name = %s"%potential_name)
print("domain_file = %s"%domain_file)
print("sc_prefix = %s"%sc_prefix)
print("ca_prefix = %s"%ca_prefix)
print("")

if domain_file != "None":
    nscal_0 = '1'
    nscal = 1
    print('domain_file is defined, nscal will be ignored.\n')

if casm != 0 and casm != 1:
    print('ERROR: casm can only be either 0 (ca model) or 1 (ca-sidechain model).')

if potential_name.upper().startswith('GENERIC'):
    words = potential_name.split('-')
    if len(words) == 1:
        print("ERROR: Generic potential keyword must be invoked as 'generic-bt'")
        sys.exit()
    else:
        if words[-1].upper() != 'BT' and words[-1].upper() != 'MJ' and words[-1].upper() != 'KGS':
            print("ERROR: You can only invoke Generic potential keyword as 'generic-bt' or 'generic-mj' or 'generic-kgs'")
            sys.exit()
        else:
            potential_name = potential_name.upper()
            print("ERROR: The generic potential is not supported in this version.\nCoarse-graining terminated.")
            sys.exit()
else:
    potential_name = potential_name.upper()
### END: get info from control file ###

## BEGIND: Conditional Defaults ##
if casm == 1:  
    ene_bsc = 0.37 # energy of a backbone-sidechain native contact (0.03 in old version)
    single_hbond_ene = 0.75 # energy of a hydrogen bond for everthing but helices (0.50 in old version)
    single_hbond_ene_helix = 0.75 # energy of a hydrogen bond in a helix (0.50 in old version)
    bondlength_go = 0 # non-Go bond length
    angle_dw = 0 # Go angle potential
    dihedral_go = 1 # Go dihedral potential
    improperdihed_go = 1 # Go improper dihedral potential
    
else: 
    ene_bsc = 0.37;  
    single_hbond_ene = 0.75; # energy of a hydrogen bond for everthing but helices
    single_hbond_ene_helix = 0.75; # energy of a hydrogen bond in a helix
    bondlength_go = 0 # non-Go bond length
    angle_dw = 1 # double-well angle potential
    dihedral_go = 0 # non-Go dihedral potential
    improperdihed_go = 0 # non-Go improper dihedral potential

# read domain nscal values if domain is defined
dom_nscal = []
ndomain = 0
dom = []
if domain_file != "None":
    if not os.path.exists(domain_file):
        print("ERROR: File %s does not exist"%domain_file)
        sys.exit()
    f = open(domain_file)
    lines = f.readlines()
    f.close()
    for line in lines:
        line = line.strip()
        if line.startswith('scale factor'):
            words = line.split('=')
            dom_nscal.append(float(words[-1]))
        if line.startswith('domain'):
            ndomain += 1
            words = line.split('=')[-1].split('-')
            words = [int(w) for w in words]
            dom.append(words)
            if words[0] > words[1]:
                print("ERROR: When defining the domains in the interface file, index %d is Greater than %d!"%(words[0], words[1]))
                sys.exit()
    print('%d domain(s) defined in the Domain file %s'%(ndomain, domain_file))
    if ndomain == 0:
        print("ERROR: No domain definitions were read. Check the domain definition file!")
        sys.exit()
    print("Domain information:")
    for i, d in enumerate(dom):
        print("Domain %d: %d to %d"%(i+1, d[0], d[1]))
    print("")
    if len(dom_nscal) != (1+ndomain)*ndomain/2:
        print("ERROR: Incorrect number of interfaces assigned. (%d, should be %d)"%(len(dom_nscal)-ndomain, (ndomain-1)*ndomain/2))
        sys.exit()
# END read domain nscal values if domain is defined

# initialize nonbonding potential
root_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
if potential_name.startswith('MJ'):
    miya = root_dir+"/shared_files/mj_contact_potential.dat"
elif potential_name.startswith('KGS'):
    miya = root_dir+"/shared_files/kgs_contact_potential.dat"
elif potential_name.startswith('BT'):
    miya = root_dir+"/shared_files/bt_contact_potential.dat"
else:
    print("ERROR: Unrecognized force-field %s"%potential_name)
    sys.exit()

eps = np.zeros((20,20))

f = open(miya)
lines = f.readlines()
f.close()
nrows = 0
avg_mj = 0
nmj = 0
for line in lines:
    line = line.strip()
    if line.startswith('#'):
        continue
    if line.startswith('AA'):
        words = line.split()
        vec = []
        for w in words[1:]:
            vec.append(res2n[w.upper()])
        if len(vec) != 20:
            print("ERROR: missing residues in file %s"%miya)
            sys.exit()
    else:
        words = line.split()
        for tc, w in enumerate(words):
            w = float(w)
            if potential_name.startswith('MJ'):
                eps[vec[nrows]][vec[tc]] = nscal * abs(w-1.2)
                eps[vec[tc]][vec[nrows]] = nscal * abs(w-1.2)
                avg_mj += nscal * abs(w-1.2)
            elif potential_name.startswith('BT'):
                eps[vec[nrows]][vec[tc]] = nscal * abs(w-0.6)
                eps[vec[tc]][vec[nrows]] = nscal * abs(w-0.6)
                avg_mj += nscal * abs(w-0.6)
            elif potential_name.startswith('KGS'):
                eps[vec[nrows]][vec[tc]] = nscal * abs(w-1.8)
                eps[vec[tc]][vec[nrows]] = nscal * abs(w-1.8)
                avg_mj += nscal * abs(w-1.8)
            nmj += 1
        nrows += 1
        if nrows > 20:
            print("ERROR 2: missing residues in file %s: %d"%(miya, nrows))
            sys.exit()
        if len(words) != nrows:
            print("ERROR 3: missing residues in file %s, %d != %d"%(miya, len(words), nrows))
            sys.exit()

avg_mj = avg_mj/nmj
print("The average %s interaction energy is %.4f\n"%(potential_name, avg_mj))
# END initialize nonbonding potential

# Read in the generic backbone dihedral potential of CL Brooks if NON-GO dihedrals
# requested by user.
if dihedral_go == 0:
    dihedb_nongo = [[[] for j in range(20)] for i in range(20)]
    f = open(root_dir+"/shared_files/karanicolas_dihe_parm.dat")
    lines = f.readlines()
    f.close()
    nphi = 0
    r1_old = None
    r2_old = None
    for line in lines:
        line = line.strip()
        dat = line.split()
        r1 = dat[0].upper()
        r2 = dat[1].upper()
        if r1 != r1_old or r2 != r2_old:
            nphi = 0
        dihedb_nongo[res2n[r1]][res2n[r2]].append([0.756*float(dat[2]), int(dat[3]), float(dat[4])])
        nphi += 1
        r1_old = r1
        r2_old = r2
        if nphi > 4:
            print("ERROR: nphi = %d upon reading in generic dihedral file"%nphi)
            print(line)
            sys.exit()
# END Read in the generic backbone dihedral potential

resname_prefix = 'G'
atomname_prefix = ''
        
# Read PDB file
cg_structure = pmd.Structure()
print("Reading in PDB file %s"%pdbfile)

struct = pmd.load_file(pdbfile)
sel_idx = np.zeros(len(struct.atoms))
for idx, res in enumerate(struct.residues):
    res.number = idx+1
    if res.name in aa:
        for atm in res.atoms:
            if atm.element != 1:
                sel_idx[atm.idx] = 1
heavy_protein = struct[sel_idx]

for idx, res in enumerate(heavy_protein.residues):
    num_backbone = 0
    num_sidechain = 0
    for atm in res.atoms:
        if atm.name in ['C', 'N', 'O', 'CA']:
            num_backbone += 1
        elif atm.name != 'OXT':
            num_sidechain += 1
    if num_backbone != 4:
        print("ERROR: In pdb the number of backbone atoms in residue %d is incorrect: %d != 4"%(idx+1, num_backbone))
        sys.exit()
    if num_sidechain != refNscat[res.name]:
        print("ERROR: In pdb the number of sidechain atoms in residue %d is incorrect: %d != %d"%(idx+1, num_sidechain, refNscat[res.name]))
        sys.exit()

idx_atm = 0
ca_list = []
chain_id_list = []
for res in heavy_protein.residues:
    if not res.chain in chain_id_list:
        chain_id_list.append(res.chain)
if len(chain_id_list) > len(alphabet):
    print('ERROR: The number of chains in pdb file (%d) exceeds the maximum (%d)'%(len(chain_id_list), len(alphabet)))
    sys.exit()
resid = 0
chainid = chain_id_list[0]
for idx, res in enumerate(heavy_protein.residues):
    if res.segid == '':
        segid = alphabet[chain_id_list.index(res.chain)]
    else:
        segid = res.segid
    
    if res.chain != chainid:
        chainid = res.chain
        resid = 1
    else:
        resid += 1
    
    SC_Mass = aaSCmass[res.name] - aaSCmass['GLY']
    CA_Mass = aaSCmass['GLY']
    SC_COM = np.zeros(3)
    CA_COM = np.zeros(3)
    sum_SC_Mass = 0
    
    for atm in res.atoms:
        if atm.name not in ['C', 'N', 'O', 'CA', 'OXT']:
            sum_SC_Mass += atm.mass
            SC_COM += atm.mass * np.array([atm.xx, atm.xy, atm.xz])
        elif atm.name == 'CA':
            CA_COM[0] = atm.xx
            CA_COM[1] = atm.xy
            CA_COM[2] = atm.xz
    if sum_SC_Mass == 0:
        is_gly = True
    else:
        is_gly = False
        SC_COM /= sum_SC_Mass
    
    if casm == 0:
        cg_atm = pmd.topologyobjects.Atom(name=atomname_prefix+ca_prefix, 
                                          type=ca_prefix+str(idx+1), charge=refcharge[res.name], 
                                          mass=aaSCmass[res.name], number=idx_atm+1)
        cg_atm.xx = CA_COM[0]
        cg_atm.xy = CA_COM[1]
        cg_atm.xz = CA_COM[2]
        cg_structure.add_atom(cg_atm, resname_prefix+str(idx+1), resid, segid=segid, chain=res.chain)
        idx_atm += 1
        ca_list.append(cg_atm)
    else:
        ca_atm = pmd.topologyobjects.Atom(name=atomname_prefix+ca_prefix, 
                                          type=ca_prefix+str(idx+1), charge=0.0, 
                                          mass=CA_Mass, number=idx_atm+1)
        ca_atm.xx = CA_COM[0]
        ca_atm.xy = CA_COM[1]
        ca_atm.xz = CA_COM[2]
        cg_structure.add_atom(ca_atm, resname_prefix+str(idx+1), resid, segid=segid, chain=res.chain)
        idx_atm += 1
        ca_list.append(ca_atm)
        
        if not is_gly:
            sc_atm = pmd.topologyobjects.Atom(name=atomname_prefix+sc_prefix, 
                                              type=sc_prefix+str(idx+1), charge=refcharge[res.name], 
                                              mass=SC_Mass, number=idx_atm+1)
            sc_atm.xx = SC_COM[0]
            sc_atm.xy = SC_COM[1]
            sc_atm.xz = SC_COM[2]
            cg_structure.add_atom(sc_atm, resname_prefix+str(idx+1), resid, segid=segid, chain=res.chain)
            idx_atm += 1

# Assign domain id to atom
if ndomain != 0:
    print('Assign domain id to each atom')
    id_domain = []
    for atm in cg_structure.atoms:
        res_id = atm.residue.idx+1
        found = False
        for i, di in enumerate(dom):
            if res_id >= di[0] and res_id <= di[1]:
                id_domain.append(i)
                found = True
                break
        if not found:
            print('ERROR: %s is not located in any domain.'%atm)
            sys.exit()
    print('')

# Write psf, cor and top
output_prefix = pdbfile.strip().split('/')[-1].split('.pdb')[0]
if casm == 1:
    output_prefix += '_ca-cb'
else:
    output_prefix += '_ca'
print('Create psf')
create_psf(cg_structure, ca_list, output_prefix)
print('Create cor')
cg_structure.save(output_prefix+'.cor', overwrite=True, format='charmmcrd')
print('Create top\n')
Create_rtf(cg_structure, output_prefix)
    
# Prepare FF parameters
print("Determining native contacts")
dist_map = np.zeros((len(cg_structure.atoms), len(cg_structure.atoms)))
for idx_1, atm_1 in enumerate(cg_structure.atoms):
    for idx_2, atm_2 in enumerate(cg_structure.atoms):
        dist_map[idx_1, idx_2] = calc_distance(atm_1, atm_2)
print("Finished calculating distance matrix")

## Compute native contacts between side-chains
print("Determining side-chains - side-chains contacts")
native_ss_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
for i in range(len(cg_structure.residues)-3):
    res_1 = heavy_protein.residues[i]
    for j in range(i+3, len(cg_structure.residues)): # separate by 2 residues
        res_2 = heavy_protein.residues[j]
        found = False
        for atm_1 in res_1.atoms:
            for atm_2 in res_2.atoms:
                if not atm_1.name in ['C', 'N', 'O', 'CA', 'OXT'] and not atm_2.name in ['C', 'N', 'O', 'CA', 'OXT']:
                    dij = calc_distance(atm_1, atm_2)
                    if dij <= heav_cut:
                        native_ss_map[i,j] = 1
                        native_ss_map[j,i] = 1
                        found = True
                        break
            if found:
                break
## Compute native contacts between backbone and side-chains
print("Determining backbone - side-chains contacts")
native_bsc_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
for i in range(len(cg_structure.residues)):
    res_1 = heavy_protein.residues[i]
    for j in range(len(cg_structure.residues)):
        res_2 = heavy_protein.residues[j]
        if i < j-2 or i > j+2: # separate by 2 residues
            found = False
            for atm_1 in res_1.atoms:
                for atm_2 in res_2.atoms:
                    if atm_1.name in ['C', 'N', 'O', 'CA', 'OXT'] and atm_2.name not in ['C', 'N', 'O', 'CA', 'OXT']:
                        dij = calc_distance(atm_1, atm_2)
                        if dij <= heav_cut:
                            native_bsc_map[i,j] = 1
                            found = True
                            break
                if found:
                    break
print('# nat sc-sc contacts %d, # nat bb-sc contacts %d, and  # non-nat sc-sc %d'%(np.sum(native_ss_map)/2, 
      np.sum(native_bsc_map), (len(cg_structure.residues)-3)*(len(cg_structure.residues)-2)/2 - np.sum(native_ss_map)/2))
      
## Determine hydrogen bonds that are present using STRIDE,
## and assign to Calpha-Calpha pairs. Also secondary structural elements
## within the native structure.
print("Determining the presence of hydrogen bonds using STRIDE")
native_hb_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
helical_list = np.zeros(len(cg_structure.residues))
hb_ene_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
screen_out = os.popen('stride -h %s'%pdbfile).readlines()
for line in screen_out:
    line = line.strip()
    resid = 0
    if line.startswith('ASD '):
        if 'Helix' in line.split()[6]:
            helical_list[resid] = 1
        resid += 1
    if line.startswith('ACC ') or line.startswith('DNR '):
        # Get H-bonding info
        resid_1 = int(line[16:20])+1
        resid_2 = int(line[36:40])+1
        chainid_1 = line[8:10].strip()
        if chainid_1 == '-':
            chainid_1 = ''
        chainid_2 = line[28:30].strip()
        if chainid_2 == '-':
            chainid_2 = ''
        found = [0, 0]
        for idx, res in enumerate(cg_structure.residues):
            if res.number == resid_1 and res.chain == chainid_1:
                idx_1 = idx
                found[0] = 1
            elif res.number == resid_2 and res.chain == chainid_2:
                idx_2 = idx
                found[1] = 1
            if sum(found) == 2:
                break
        if sum(found) != 2:
            print("ERROR: Cannot find residue in parmed structure according to the Hbond info.\n  %s"%line)
            sys.exit()
        if chainid_1 == chainid_2:
            if idx_1 < idx_2:
                if native_hb_map[idx_1, idx_2] == 1:
                    if helical_list[idx_1] == 1 and helical_list[idx_2] == 1:
                        hb_ene_map[idx_1, idx_2] = 2*single_hbond_ene_helix
                        hb_ene_map[idx_2, idx_1] = 2*single_hbond_ene_helix
                    else:
                        hb_ene_map[idx_1, idx_2] = 2*single_hbond_ene
                        hb_ene_map[idx_2, idx_1] = 2*single_hbond_ene
                else:
                    native_hb_map[idx_1, idx_2] = 1
                    native_hb_map[idx_2, idx_1] = 1
                    if helical_list[idx_1] == 1 and helical_list[idx_2] == 1:
                        hb_ene_map[idx_1, idx_2] = single_hbond_ene_helix
                        hb_ene_map[idx_2, idx_1] = single_hbond_ene_helix
                    else:
                        hb_ene_map[idx_1, idx_2] = single_hbond_ene
                        hb_ene_map[idx_2, idx_1] = single_hbond_ene
        else:
            native_hb_map[idx_1, idx_2] = 1
            native_hb_map[idx_2, idx_1] = 1
            hb_ene_map[idx_1, idx_2] = single_hbond_ene
            hb_ene_map[idx_2, idx_1] = single_hbond_ene
num_hb = 0
for i in range(len(cg_structure.residues)-1):
    for j in range(i+1, len(cg_structure.residues)):
        if native_hb_map[i,j] == 1:
            num_hb += 1
            #print('%d %.4f, %d %d'%(num_hb, hb_ene_map[i,j], i+1, j+1))
print('# of unique Hbonds %d'%num_hb)
native_contact_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
for i in range(len(cg_structure.residues)):
    for j in range(len(cg_structure.residues)):
        if native_ss_map[i,j] == 1 or native_bsc_map[i,j] == 1 or native_hb_map[i,j] == 1:
            native_contact_map[i,j] == 1

## Write prm file ##
print('\nCreate prm\n')
prmfile = pdbfile.strip().split('/')[-1].split('.pdb')[0] + '_nscal' + nscal_0 + '_fnn' + fnn_0 + '_go_' + potential_name.lower() + '.prm'
f = open(prmfile, 'w')
f.write('* This CHARMM .param file describes a Go model of %s\n'%(pdbfile.split('/')[-1]))
f.write('*\n\n')
# Atomic mass
f.write('ATOM\n')
for idx, atm in enumerate(cg_structure.atoms):
    f.write('MASS %-5s %-8s %-10.6f\n'%(str(idx+1), atm.type, atm.mass))
f.write('\n')
# Bond section (non-go bondlength for both models)
f.write('BOND\n')
kb = 50.0
for idx, bond in enumerate(cg_structure.bonds):
    if bondlength_go == 0:
        if bond.atom2.name == (atomname_prefix+sc_prefix):
            res_idx = bond.atom2.residue.idx
            bond_length = lbs_nongo[heavy_protein.residues[res_idx].name]
            f.write('%-8s%-10s%-12.6f%-9.6f\n'%(bond.atom1.type, bond.atom2.type, kb, bond_length))
        else:
            f.write('%-8s%-10s%-12.6f%-9.6f\n'%(bond.atom1.type, bond.atom2.type, kb, 3.81))
    else:
        f.write('%-8s%-10s%-12.6f%-9.6f\n'%(bond.atom1.type, bond.atom2.type, kb, bond.measure()))
f.write('\n')
# Angle section
f.write('ANGLE\n')
ka = 30.0
for idx, angle in enumerate(cg_structure.angles):
    if angle_dw == 0:
        f.write('%-8s%-8s%-10s%11.6f%11.6f\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type, 
                                               ka, angle.measure()))
    else:
        if angle.atom1 == (atomname_prefix+sc_prefix):
            res_idx = angle.atom1.residue.idx
            angle_value = ang_sb_nongo[heavy_protein.residues[res_idx].name]
            f.write('%-8s%-8s%-10s%11.6f%11.6f\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type, 
                                                   ka, angle_value))
        elif angle.atom3 == (atomname_prefix+sc_prefix):
            res_idx = angle.atom3.residue.idx
            angle_value = ang_bs_nongo[heavy_protein.residues[res_idx].name]
            f.write('%-8s%-8s%-10s%11.6f%11.6f\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type, 
                                                   ka, angle_value))
        else:
            f.write('%-8s%-8s%-10s  106.4 91.7 26.3 130.0 0.1 4.3\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type))
f.write('\n')
# Dihedral section
f.write('DIHEDRAL\n')
f.write('! backbone dihedrals\n')
for idx, dihedral in enumerate(cg_structure.dihedrals):
    if dihedral_go == 1: # Use Go backbone dihedral angles
        delta = 1*dihedral.measure()-180
        if casm == 1:
            if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                kd = 0.30
            else: # not helical
                kd = 0.55
            f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                               dihedral.atom4.type, kd, 1, delta))
            delta = 3*dihedral.measure()-180
            if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                kd = 0.15
            else: # not helical
                kd = 0.275
            f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                               dihedral.atom4.type, kd, 3, delta))
        else:
            if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                kd = 0.75
            else: # not helical
                kd = 0.75
            f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                               dihedral.atom4.type, kd, 1, delta))
            delta = 3*dihedral.measure()-180
            if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                kd = 0.275
            else: # not helical
                kd = 0.275
            f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                               dihedral.atom4.type, kd, 3, delta))
    else: # Use Non-go dihedrals
        for i in range(4):
            res_idx_1 = dihedral.atom2.residue.idx
            res_idx_2 = dihedral.atom3.residue.idx
            [kd, period, delta] = dihedb_nongo[res2n[heavy_protein.residues[res_idx_1].name]][res2n[heavy_protein.residues[res_idx_2].name]][i]
            f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                               dihedral.atom4.type, kd, period, delta))
f.write('\n')
# Improper dihedral section
f.write('IMPHI\n')
f.write('! sidechain improper dihedrals to maintain chirality\n')
if casm == 1:
    for idx, improper in enumerate(cg_structure.impropers):
        if improperdihed_go == 1:
            angle = improper.measure()
        else:
            res_idx = improper.atom1.residue.idx
            angle = improper_nongo[heavy_protein.residues[res_idx].name] # use transferable improper dihedral
        delta = angle + 180
        kd = 20*abs(avg_mj)
        f.write('%-5s %-5s %-5s %-7s%.6f %-3d%-10.5f\n'%(improper.atom1.type, improper.atom2.type, improper.atom3.type, 
                                                           improper.atom4.type, kd, 1, delta))
f.write('\n')

## nonbonded section
f.write('NONBONDED NBXMOD 3 ATOM CDIEL SWITCH VATOM VDISTANCE VSWITCH -\n')
f.write('CUTNB 32 CTOFNB 20 CTONNB 18 EPS 78.5 WMIN 1.5 E14FAC 1.0\n')
f.write('!atom           e_min   r_min/2\n')
# if using the C-alpha only model do some preprocessing to determine the collision
# diameter of non-native interactions according to the Karanacolis-Brooks
# algorithm
if casm != 1:
    sigmin = 1000000*np.ones(len(cg_structure.residues))
    if potential_name.startswith('GENERIC'):
        for idx, res in enumerate(cg_structure.residues):
            sigmin[idx] = 2*rvdw[heavy_protein.residues[idx].name]
    else:
        # determine the collision diameter
        for i in range(len(cg_structure.residues)):
            for j in range(len(cg_structure.residues)):
                if native_contact_map[i,j] != 1 and (j < i-2 or j > i+2):
                    if dist_map[i,j] < sigmin[i]:
                        sigmin[i] = dist_map[i,j]
    for idx, atm in enumerate(cg_structure.atoms):
        eps2 = -0.000132
        rmin2 = sigmin[idx]*2**(1/6)/2
        temp = fnn*rmin2
        f.write("%-9s%-5.1f%-9.6f    %-10.6f\n"%(atm.type, 0.0, eps2, temp))
else:
    eps2 = '-1e-12' #!!!! SYSTem dependent !!!!!!!!
    rmin2 = 20.0
    for idx, atm in enumerate(cg_structure.atoms):
        if atm.name == (atomname_prefix+ca_prefix):
            f.write("%-9s%-5.1f%-s    %-10.6f\n"%(atm.type, 0.0, eps2, rmin2))
        else:
            t1 = 1
            t2 = (t1*(2*rvdw[heavy_protein.residues[atm.residue.idx].name]*2**(1/6))**12/(1e-12))**(1/12)
            temp = fnn*t2/2
            f.write("%-9s%-5.1f%-s    %-10.6f\n"%(atm.type, 0.0, eps2, temp))
f.write('\n')
## NBFIX section
f.write('NBFIX\n')
### native side-chain pairs and backbone Hbonding
if casm == 1:
    f.write('! b-b due to Hbonding\n')
    totene_bb = 0
    for i in range(len(cg_structure.residues)-1):
        for j in range(i+1, len(cg_structure.residues)):
            if native_hb_map[i,j] == 1:
                atm_i = cg_structure.residues[i].atoms[0]
                atm_j = cg_structure.residues[j].atoms[0]
                comment = ''
                if ndomain == 0: # No domain defined
                    ene = hb_ene_map[i,j]
                else: # Domain defined
                    if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                        di = id_domain[atm_i.idx]
                        comment = '! in Domain %d'%(di+1)
                        ene = hb_ene_map[i,j]
                    else: # in the interface
                        di = id_domain[atm_i.idx]
                        dj = id_domain[atm_j.idx]
                        comment = '! in Interface %d | %d'%(di+1, dj+1)
                        ii =  int((2*ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                        #ene = dom_nscal[ii] # ??? Use nscal at interface
                        ene = hb_ene_map[i,j] # ??? Use the same energy
                f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
                totene_bb += ene
    
    totene_sc = 0
    totene_bsc = 0
    if potential_name.startswith('GENERIC'): # C-alpha - side chain model Generic non-bond interactions
        f.write('!Generic interactions between unstructured portions of this protein\n')
        # Print out NBFIX energy values
        for i in range(len(cg_structure.residues)-3):
            resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
            for j in range(i+3, len(cg_structure.residues)):
                resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                atm_i = cg_structure.residues[i].atoms[1] # ??? should be side-chain
                atm_j = cg_structure.residues[j].atoms[1] # ??? should be side-chain
                temp = rvdw[resname_1] + rvdw[resname_2]
                ene=(0.3/10)*eps[res2n[resname_1]][res2n[resname_2]]
                f.write('%-8s%-11s%-13.6f%-11.6f\n'%(atm_i.type, atm_j.type, -ene, temp))
    else: # Go non-bond interactions 
        f.write('! native side-chain interactions\n')
        for i in range(len(cg_structure.residues)-1):
            resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
            for j in range(i+1, len(cg_structure.residues)):
                resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                if native_ss_map[i,j] == 1:
                    atm_i = cg_structure.residues[i].atoms[1]
                    atm_j = cg_structure.residues[j].atoms[1]
                    if eps[res2n[resname_1]][res2n[resname_2]] == 0:
                        print('ERROR 1: Well depth equal to zero!!! %s - %s'%(resname_1, resname_2))
                        sys.exit()
                    comment = ''
                    if ndomain == 0: # No domain defined
                        ene = eps[res2n[resname_1]][res2n[resname_2]]
                    else: # If domain is defined
                        if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                            di = id_domain[atm_i.idx]
                            comment = '! in Domain %d'%(di+1)
                            ene = eps[res2n[resname_1]][res2n[resname_2]] * dom_nscal[di] 
                        else: # in the interface
                            di = id_domain[atm_i.idx]
                            dj = id_domain[atm_j.idx]
                            comment = '! in Interface %d | %d'%(di+1, dj+1)
                            ii =  int((2*ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                            ene = eps[res2n[resname_1]][res2n[resname_2]] * dom_nscal[ii] 
                    f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
                    totene_sc += ene
                    
        f.write('! backbone-sidechain interactions\n')
        for i in range(len(cg_structure.residues)): 
            resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
            for j in range(len(cg_structure.residues)):
                resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                if native_bsc_map[i,j] == 1:
                    atm_i = cg_structure.residues[i].atoms[0] # backbone
                    atm_j = cg_structure.residues[j].atoms[1] # sidechain
                    comment = ''
                    if ndomain == 0: # No domain defined
                        ene = ene_bsc
                    else: # If domain is defined
                        if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                            di = id_domain[atm_i.idx]
                            comment = '! in Domain %d'%(di+1)
                            ene = ene_bsc
                        else: # in the interface
                            di = id_domain[atm_i.idx]
                            dj = id_domain[atm_j.idx]
                            comment = '! in Interface %d | %d'%(di+1, dj+1)
                            ii =  int((2*ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                            ene = ene_bsc * dom_nscal[ii] # Rescaled energy
                    f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
                    totene_bsc += ene
                    
    f.write('\n')
    f.write('! %.4f, %.4f, %.4f\n'%(totene_bb, totene_sc, totene_bsc))
else:
    if not potential_name.startswith('GENERIC'): # C-alpha model
        f.write('! b-b due to Hbonding plus native side-chain interactions plus backbone-sidechain interactions\n')
        # Add up non-bonded energies
        for i in range(len(cg_structure.residues)-1): 
            resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
            for j in range(i+1, len(cg_structure.residues)):
                resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                atm_i = cg_structure.residues[i].atoms[0]
                atm_j = cg_structure.residues[j].atoms[0]
                ene = 0
                # hydrogen bonds
                if native_hb_map[i,j] == 1:
                    if ndomain == 0: # No domain defined
                        ene += hb_ene_map[i,j]
                    else: # Domain defined
                        if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                            di = id_domain[atm_i.idx]
                            ene += hb_ene_map[i,j]
                        else: # in the interface
                            di = id_domain[atm_i.idx]
                            dj = id_domain[atm_j.idx]
                            ii =  int((2*ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                            ene += hb_ene_map[i,j] # Use the same energy
                # sc-sc interactions
                if native_ss_map[i,j] == 1:
                    if eps[res2n[resname_1]][res2n[resname_2]] == 0:
                        print('ERROR 1: Well depth equal to zero!!! %s - %s'%(resname_1, resname_2))
                        sys.exit()
                    if ndomain == 0: # No domain defined
                        ene += eps[res2n[resname_1]][res2n[resname_2]]
                    else: # Domain defined
                        if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                            di = id_domain[atm_i.idx]
                            ene += eps[res2n[resname_1]][res2n[resname_2]] * dom_nscal[di]
                        else: # in the interface
                            di = id_domain[atm_i.idx]
                            dj = id_domain[atm_j.idx]
                            ii =  int((2*ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                            ene += eps[res2n[resname_1]][res2n[resname_2]] * dom_nscal[ii] 
                # b-sc interactions
                if native_bsc_map[i,j] == 1:
                    if ndomain == 0: # No domain defined
                        ene += ene_bsc
                    else: # Domain defined
                        if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                            di = id_domain[atm_i.idx]
                            ene += ene_bsc
                        else: # in the interface
                            di = id_domain[atm_i.idx]
                            dj = id_domain[atm_j.idx]
                            ii =  int((2*ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                            ene += ene_bsc # Use the same energy
                if native_bsc_map[j,i] == 1:
                    if ndomain == 0: # No domain defined
                        ene += ene_bsc
                    else: # Domain defined
                        if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                            di = id_domain[atm_i.idx]
                            ene += ene_bsc
                        else: # in the interface
                            di = id_domain[atm_i.idx]
                            dj = id_domain[atm_j.idx]
                            ii =  int((2*ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                            ene += ene_bsc # Use the same energy
                
                # Write NBFIX
                if ene != 0:
                    comment = ''
                    if ndomain != 0:
                        if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                            di = id_domain[atm_i.idx]
                            comment = '! in Domain %d'%(di+1)
                        else: # in the interface
                            di = id_domain[atm_i.idx]
                            dj = id_domain[atm_j.idx]
                            comment = '! in Interface %d | %d'%(di+1, dj+1)
                    f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
    else:
        f.write('!Generic interactions between unstructured portions of this protein\n')
        # Print out NBFIX energy values
        for i in range(len(cg_structure.residues)-3):
            resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
            for j in range(i+3, len(cg_structure.residues)):
                resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                atm_i = cg_structure.residues[i].atoms[0]
                atm_j = cg_structure.residues[j].atoms[0]
                temp = rvdw[resname_1] + rvdw[resname_2]
                ene=(0.3/10)*eps[res2n[resname_1]][res2n[resname_2]]
                f.write('%-8s%-11s%-13.6f%-11.6f\n'%(atm_i.type, atm_j.type, -ene, temp))
f.write('\nEND\n')
f.close()

print('All done.')
