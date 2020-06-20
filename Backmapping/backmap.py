#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr
import getopt, os, time, random, math, traceback, io, sys
import parmed as pmd
import numpy as np
import mdtraj as mdt

usage = '\nUsage: python backmap.py\n' \
        '       --aa_pdb | -i <xxx.pdb> initial pdb file used to create the target C-alpha CG protein\n'\
        '       --cg_pdb | -c <xxx.pdb> pdb file for the target C-alpha CG protein\n'\
        '       [-- nproc | -n <NUM>] number of CPU processors used to run energy minimizations\n'\
        '       [--help | -h] Print this information\n\n'

def clean_pdb(pdb, out_dir):
    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP'];
    name = pdb.split('/')[-1].split('.pdb')[0]
    struct = pmd.load_file(pdb)
    sel_idx = np.zeros(len(struct.atoms))
    for idx, res in enumerate(struct.residues):
        res.number = idx+1
        if res.name in AA_name_list:
            for atm in res.atoms:
                sel_idx[atm.idx] = 1
    struct[sel_idx].save(out_dir+'/'+name+'_clean.pdb', overwrite=True)
    return name+'_clean.pdb'
    
def create_psf(name):
    segid = 'A'
    parm = pmd.charmm.CharmmParameterSet(name+'.top')
    f = open(name+'.seq','r')
    seq = f.readlines()[0].strip().split()
    f.close()
    struct = pmd.Structure()
    for resname in seq:
        struct += parm.residues[resname].to_structure()
    ca_list = []
    for atm in struct.atoms:
        atm.mass = parm.atom_types[atm.type].mass
        if atm.name == 'A':
            ca_list.append(atm)
    # creat backbond bonds
    for i in range(len(ca_list)-1):
        struct.bonds.append(pmd.topologyobjects.Bond(ca_list[i], ca_list[i+1]))
    # create Angles
    for atm in struct.atoms:
        bond_list = atm.bond_partners
        if len(bond_list) > 1:
            for i in range(len(bond_list)-1):
                for j in range(i+1, len(bond_list)):
                    struct.angles.append(pmd.topologyobjects.Angle(bond_list[i], atm, bond_list[j]))
    # create Dihedrals
    for i in range(len(ca_list)-3):
        struct.dihedrals.append(pmd.topologyobjects.Dihedral(ca_list[i], ca_list[i+1], ca_list[i+2], ca_list[i+3]))
    # create Impropers
    for i in range(1, len(ca_list)-1):
        if len(ca_list[i].residue.atoms) > 1:
            b_bead = ca_list[i].residue.atoms[1]
            struct.impropers.append(pmd.topologyobjects.Improper(ca_list[i], ca_list[i-1], ca_list[i+1], b_bead))
    for res in struct.residues:
        res.segid = segid
    struct.save(name+'.psf', overwrite=True)

def create_cg_model(pdb):
    os.system("mkdir create_model")
    os.chdir("create_model")

    fo = open("go_model.cntrl", "w")
    fo.write("pdbfile = ../%s\n"%pdb)
    fo.write("nscal = %.1f\n"%10)
    fo.write("potential_name = mj\n")
    fo.write("casm = 1\n")
    fo.close()
  
    os.system("create_cg_protein_model.py -f go_model.cntrl > go_model.log 2>&1");
    
    name = pdb.split('.pdb')[0].lower()
    prefix = name+'_ca-cb'
    prm_name = name + '_nscal10.0_fnn1_go_mj.prm'
    
    if os.path.exists(prefix+'.psf'):
        os.system('cp *.psf ../')
        os.system('cp *.cor ../')
        os.system('cp *.top ../')
        os.system('cp *.prm ../')
        os.chdir('../')
    else:
        print("Error: failed to create CG model from %s\n\n"%pdb)
        sys.exit()
    return (prefix, prm_name)
    
def add_sc_beads(cg_pdb, cacb_struct):
    cor = pmd.load_file(cg_pdb)
    cor = cor.coordinates
    new_cacb_struct = cacb_struct.copy(pmd.Structure)
    idx = 0
    for res in new_cacb_struct.residues:
        res.atoms[0].xx = cor[idx,0]
        res.atoms[0].xy = cor[idx,1]
        res.atoms[0].xz = cor[idx,2]
        if len(res.atoms) > 1:
            cor1 = cacb_struct.coordinates[res.atoms[0].idx,:]
            cor2 = cacb_struct.coordinates[res.atoms[1].idx,:]
            bond_length = np.sum((cor1-cor2)**2)**0.5
            res.atoms[1].xx = cor[idx,0] + bond_length
            res.atoms[1].xy = cor[idx,1]
            res.atoms[1].xz = cor[idx,2]
        idx += 1
    return new_cacb_struct

def cacb_energy_minimization(cor, prefix, prm_file):
    global nproc
    temp = 310
    timestep = 0.015*picoseconds
    fbsolu = 0.05/picosecond
    temp = temp*kelvin

    psf_pmd = pmd.charmm.CharmmPsfFile(prefix+'.psf')
    psf = CharmmPsfFile(prefix+'.psf')
    top = psf.topology
    os.system('parse_cg_cacb_prm.py -p '+prm_file+' -t '+prefix+'.top')
    name = prm_file.split('.prm')[0]
    forcefield = ForceField(name+'.xml')
    
    # re-name residues that are changed by openmm
    for resid, res in enumerate(top.residues()):
        if res.name != psf_pmd.residues[resid].name:
            res.name = psf_pmd.residues[resid].name
    
    template_map = {}
    for chain in top.chains():
        for res in chain.residues():
            template_map[res] = res.name
                
    
    system = forcefield.createSystem(top, nonbondedCutoff=2.0*nanometer, switchDistance=1.8*nanometer, 
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
    properties = {'Threads': nproc}
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

def Call_PD2(mini_pdb):
    min_opt = 500
    print("-> Calling PD2 to reconstruct backbone")
    PD2_dir = os.popen('which pd2_ca2main').readlines()[0].strip()
    if PD2_dir.startswith('which: no '):
        print("Error: cannot find pd2_ca2main in your PATH\n")
        sys.exit()
    else:
        words = PD2_dir.split('/bin')
        PD2_database_dir = words[0] + "/database/"
        
    output_from_PD2 = mini_pdb.split('.pdb')[0] + "_PD2_min.pdb"

    os.system("pd2_ca2main --database "+PD2_database_dir+" --ca2main:new_fixed_ca --ca2main:bb_min_steps "
              +str(min_opt)+" -i "+mini_pdb+" -o "+output_from_PD2+" > PD2_log.log 2>&1")

    f = open("PD2_log.log", 'r')
    screen = ""
    for line in f.readlines():
        screen += line
    if screen.find('Validating output PDB:') != -1:
        if screen.find('output PDB seems OK') != -1:
            print("   Backbone reconstruction completed.")
        else:
            print("   Warning: Backbone reconstruction completed with a bad structure.")
    else:
        print("Error: backbone reconstruction failed. Check PD2_log.log to see details.\n")
        sys.exit()

    return output_from_PD2

def Call_Pulchra(PD2_rebult_pdb):
    print("-> Calling pulchra to reconstruct all-atom PDB")

    os.system("pulchra -v -g -b -q "+PD2_rebult_pdb+" > pulchra.log")

    pdb_code = PD2_rebult_pdb.split('.pdb')[0]
    old_name = pdb_code + ".rebuilt.pdb"
    new_name = pdb_code + "_pulchra.pdb"
    os.system("mv "+old_name+" "+new_name)

    print("   Reconstructed all-atom PDB "+new_name)

    return new_name

def GBSA_minimization(input_pdb, maxcyc):
    ncyc = int(maxcyc/2)
    pdb_code = input_pdb.split('.pdb')[0]

    print("-> Running Amber energy minimization for "+str(maxcyc)+" steps")

    tem_dir = "Amber_min_"+str(maxcyc)
    os.system('mkdir '+tem_dir)

    f = open(tem_dir+"/leap.in", 'w')
    f.write("""source leaprc.protein.ff14SB
a = loadpdb """+input_pdb+"""
saveamberparm a """+tem_dir+"""/protein.top """+tem_dir+"""/protein.crd
quit\n""")
    f.close()

    os.system("tleap -s -f "+tem_dir+"/leap.in > /dev/null")

    f = open(tem_dir+"/min1.in", 'w')
    f.write("""min1
 &cntrl
  imin=1,
  maxcyc="""+str(maxcyc)+""",
  ncyc   = """+str(ncyc)+""",
  ntb    = 0,
  igb    = 1, intdiel = 1, extdiel = 78.5,
  gbsa = 1,
  cut    = 999,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=100.0,
 /
""")
    f.close()

    os.system("mpirun -np 4 sander.MPI  -O -i "+tem_dir+"/min1.in -p "+tem_dir+"/protein.top -c "+tem_dir+
              "/protein.crd -r "+tem_dir+"/min1.rst -o "+tem_dir+"/min1.out -ref "+tem_dir+"/protein.crd -inf "+tem_dir+"/mdinfo")
  
    energy = "NaN"
    f = open(tem_dir+"/min1.out")
    lines = f.readlines()
    idx = 0
    while idx < len(lines):
        line = lines[idx]
        if line.find("FINAL RESULTS") != -1:
            idx += 1
            while idx < len(lines):
                line = lines[idx]
                if line.find("   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER") != -1:
                    line = lines[idx+1]
                    energy = line.strip().split()[1]
                    break
                idx += 1
        idx += 1
    f.close()

    if energy == "NaN":
        print("   Warning: Failed to run Amber Energy Minimization\n   Really a bad structure")
    else:
        print("   Done")
    
    struct = pmd.load_file(tem_dir+'/protein.top')
    cor = pmd.load_file(tem_dir+'/min1.rst')
    struct.coordinates = cor.coordinates
    struct.save(pdb_code+'_min'+str(maxcyc)+'.pdb', overwrite=True)

    return (energy, pdb_code+'_min'+str(maxcyc)+'.pdb')

def Amber_water_minimization(input_pdb, maxcyc):
    ncyc = int(maxcyc/2)
    pdb_code = input_pdb.split('.pdb')[0]

    print("-> Running Amber energy minimization for "+str(maxcyc)+" steps")

    tem_dir = "Amber_min_"+str(maxcyc)
    os.system('mkdir '+tem_dir)
    
    f = open(tem_dir+"/leap.in", 'w')
    f.write("""source leaprc.protein.ff14SB
source leaprc.water.tip3p
a = loadpdb """+input_pdb+"""
solvatebox a TIP3PBOX 10.0
addions a Na+ 0.0
saveamberparm a """+tem_dir+"""/protein.top """+tem_dir+"""/protein.crd
quit\n""")
    f.close()

    os.system("tleap -s -f "+tem_dir+"/leap.in > /dev/null")

    f = open(tem_dir+"/min1.in", 'w')
    f.write("""min1
 &cntrl
  imin=1,
  maxcyc="""+str(maxcyc)+""",
  ncyc   = """+str(ncyc)+""",
  cut=8.0, ntb=1, ntc=1, ntf=1,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=100.0,
 /
""")
    f.close()

    os.system("mpirun -np 4 pmemd.MPI  -O -i "+tem_dir+"/min1.in -p "+tem_dir+"/protein.top -c "+tem_dir+
              "/protein.crd -r "+tem_dir+"/min1.rst -o "+tem_dir+"/min1.out -ref "+tem_dir+"/protein.crd -inf "+tem_dir+"/mdinfo")
  
    energy = "NaN"
    f = open(tem_dir+"/min1.out")
    lines = f.readlines()
    idx = 0
    while idx < len(lines):
        line = lines[idx]
        if line.find("FINAL RESULTS") != -1:
            idx += 1
            while idx < len(lines):
                line = lines[idx]
                if line.find("   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER") != -1:
                    line = lines[idx+1]
                    energy = line.strip().split()[1]
                    break
                idx += 1
        idx += 1
    f.close()

    if energy == "NaN":
        print("   Warning: Failed to run Amber Energy Minimization\n   Really a bad structure")
    else:
        print("   Done")
    
    struct = pmd.load_file(tem_dir+'/protein.top')
    cor = pmd.load_file(tem_dir+'/min1.rst')
    struct.coordinates = cor.coordinates
    struct['(!@/H)&(!:WAT,Na+)'].save(pdb_code+'_min'+str(maxcyc)+'.pdb', overwrite=True)

    return (energy, pdb_code+'_min'+str(maxcyc)+'.pdb')

def Amber_vacuum_minimization(input_pdb, maxcyc):
    ncyc = int(maxcyc/2)
    pdb_code = input_pdb.split('.pdb')[0]

    print("-> Running Amber energy minimization for "+str(maxcyc)+" steps")

    tem_dir = "Amber_min_"+str(maxcyc)
    os.system('mkdir '+tem_dir)
    
    f = open(tem_dir+"/leap.in", 'w')
    f.write("""source leaprc.protein.ff14SB
a = loadpdb """+input_pdb+"""
saveamberparm a """+tem_dir+"""/protein.top """+tem_dir+"""/protein.crd
quit\n""")
    f.close()

    os.system("tleap -s -f "+tem_dir+"/leap.in > /dev/null")

    f = open(tem_dir+"/min1.in", 'w')
    f.write("""min1
 &cntrl
  imin=1,
  maxcyc="""+str(maxcyc)+""",
  ncyc   = """+str(ncyc)+""",
  cut=999, ntb=0,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=100.0,
 /
""")
    f.close()

    os.system("sander -O -i "+tem_dir+"/min1.in -p "+tem_dir+"/protein.top -c "+tem_dir+
              "/protein.crd -r "+tem_dir+"/min1.rst -o "+tem_dir+"/min1.out -ref "+tem_dir+"/protein.crd -inf "+tem_dir+"/mdinfo")
  
    energy = "NaN"
    f = open(tem_dir+"/min1.out")
    lines = f.readlines()
    idx = 0
    while idx < len(lines):
        line = lines[idx]
        if line.find("FINAL RESULTS") != -1:
            idx += 1
            while idx < len(lines):
                line = lines[idx]
                if line.find("   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER") != -1:
                    line = lines[idx+1]
                    energy = line.strip().split()[1]
                    break
                idx += 1
        idx += 1
    f.close()

    if energy == "NaN":
        print("   Warning: Failed to run Amber Energy Minimization\n   Really a bad structure")
    else:
        print("   Done")
    
    struct = pmd.load_file(tem_dir+'/protein.top')
    cor = pmd.load_file(tem_dir+'/min1.rst')
    struct.coordinates = cor.coordinates
    struct['!@/H'].save(pdb_code+'_min'+str(maxcyc)+'.pdb', overwrite=True)

    return (energy, pdb_code+'_min'+str(maxcyc)+'.pdb')

def OpenMM_vacuum_minimization(input_pdb, maxcyc):
    global nproc
    pdb_code = input_pdb.split('.pdb')[0]

    print("-> Running all-atom energy minimization for %d steps in vacuum via OpenMM"%maxcyc)

    platform = Platform.getPlatformByName('CPU')
    properties = {'Threads': nproc}

    forcefield = ForceField('amber14-all.xml')
    pdb = pdbfile.PDBFile(input_pdb)

    # Check if the end residue has missing OXT atom and add if needed
    for chain in pdb.topology.chains():
        end_res = list(chain.residues())[-1]
        found = False
        for atom in end_res.atoms():
            if atom.name == 'OXT':
                found = True
            elif atom.name == 'C':
                C_atom = atom
            elif atom.name == 'CA':
                CA_atom = atom
            elif atom.name == 'O':
                O_atom = atom
        C_position = np.array(pdb.positions[C_atom.index].value_in_unit(nanometer))
        CA_position = np.array(pdb.positions[CA_atom.index].value_in_unit(nanometer))
        O_position = np.array(pdb.positions[O_atom.index].value_in_unit(nanometer))
        if not found:
            new_atom = pdb.topology.addAtom('OXT', element.oxygen, end_res)
            pdb.topology.addBond(C_atom, new_atom)
            new_position = np.dot(rotation_matrix(C_position-CA_position, np.pi), O_position-C_position) + C_position
            new_position = Quantity(value=Vec3(x=new_position[0], y=new_position[1], z=new_position[2]), unit=nanometer)
            pdb.positions.insert(O_atom.index+1, new_position)

    model = modeller.Modeller(pdb.topology, pdb.positions)
    model.addHydrogens(forcefield=forcefield, pH=7.0, variants=None, platform=platform)

    top = model.topology
    structure = pmd.openmm.load_topology(top)
    cor = model.positions
    #structure.positions = cor
    #structure.save('111.pdb', overwrite=True)
    
    system = forcefield.createSystem(top, nonbondedMethod=NoCutoff, constraints=None)
    
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
    k = 500*kilocalorie/mole/angstrom**2
    for atm in top.atoms():
        if atm.name == 'CA':
            force.addParticle(atm.index, (k, cor[atm.index][0], cor[atm.index][1], cor[atm.index][2]))
    
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator.setConstraintTolerance(0.00001)
    simulation = Simulation(top, system, integrator, platform, properties)
    simulation.context.setPositions(cor)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
    getEnergyDecomposition(stdout, simulation.context, system)
    print('   Potential energy before minimization: %.4f kcal/mol'%energy)
    simulation.minimizeEnergy(maxIterations=maxcyc)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
    getEnergyDecomposition(stdout, simulation.context, system)
    print('   Potential energy after minimization: %.4f kcal/mol'%energy)
    current_cor = simulation.context.getState(getPositions=True).getPositions()
    
    structure.positions = current_cor
    structure['!@/H'].save(pdb_code+'_OpenMM_min.pdb', overwrite=True)
    return pdb_code+'_OpenMM_min.pdb'
    
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

##################################### MAIN #######################################
if len(sys.argv) == 1:
    print(usage)
    sys.exit()

nproc = None
try:
    opts, args = getopt.getopt(sys.argv[1:],"h:i:c:n:", ["help", "aa_pdb=", "cg_pdb=", "nproc="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt in ("-i", "--aa_pdb"):
        aa_pdb = arg
    elif opt in ("-c", "--cg_pdb"):
        cg_pdb = arg
    elif opt in ("-n", "--nproc"):
        nproc = arg
        
if nproc == None:
    nproc = '1'

print("-> Cleaning PDB file %s"%aa_pdb)
name = cg_pdb.split('/')[-1].split('.pdb')[0]
work_dir = 'rebuild_'+name
os.system('mkdir '+work_dir)
aa_clean_pdb = clean_pdb(aa_pdb, work_dir)
os.chdir(work_dir)
print('   Done')

# buld ca-cb model
print("-> Building ca-cb model for %s"%aa_clean_pdb)
(prefix, prm_file) = create_cg_model(aa_clean_pdb)
print('   Done')

cacb_struct = pmd.load_file(prefix+'.psf')
cacb_cor = pmd.load_file(prefix+'.cor')
cacb_struct.coordinates = cacb_cor.coordinates

# add SC beads to cg pdb
print("-> Adding side chain beads")
target_name = name
cg_sc_struct = add_sc_beads('../'+cg_pdb, cacb_struct)
print('   Done')

# run energy minimization for cacb model
print("-> Running energy minimization for ca-cb model")
cg_sc_min_cor = cacb_energy_minimization(cg_sc_struct.positions, prefix, prm_file)
aa_pdb_struct = pmd.load_file(aa_clean_pdb)
for index in range(len(aa_pdb_struct.residues)):
    res_name = aa_pdb_struct.residues[index].name
    cg_sc_struct.residues[index].name = res_name
for atm in cg_sc_struct.atoms:
    if atm.name == 'A':
        atm.name = ' CA'
    elif atm.name == 'B':
        atm.name = ' SC'
cg_sc_struct.positions = cg_sc_min_cor
cg_sc_struct.save(target_name+'_mini.pdb', overwrite=True)
print('   Done')

output_from_PD2 = Call_PD2(target_name+'_mini.pdb')

output_from_Pultra = Call_Pulchra(output_from_PD2)

try:
    rec_pdb = OpenMM_vacuum_minimization(output_from_Pultra, 50)
    os.system('cp '+rec_pdb+' ../'+target_name+'_rebuilt.pdb')
except Exception as e:
    traceback.print_exc()
    print('Failed to run OpenMM minimization. Use Pulchra result instead.')
    os.system('cp '+output_from_Pultra+' ../'+target_name+'_rebuilt.pdb')

os.chdir('../')
