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
import os, time, traceback
import parmed as pmd
import numpy as np
import mdtraj as mdt

###### convert time seconds to hours ######
def convert_time(seconds):
    return seconds/3600
###### END convert time seconds to hours ######

###### calculate native contact fraction ######
def calc_Q(current_cor):
    global native_contact_map, native_distance_map, native_contact_num, sec_strc_def, sdist
    current_contact_num = 0
    for i in range(len(current_cor)-4):
        tag_i = 0
        for rs in sec_strc_def:
            if i >= rs[0]-1 and i <= rs[1]-1:
                tag_i = 1
                break
        if tag_i == 0:
            continue
        for j in range(i+4, len(current_cor)):
            tag_j = 0
            for rs in sec_strc_def:
                if j >= rs[0]-1 and j <= rs[1]-1:
                    tag_j = 1
                    break
            if tag_j == 0:
                continue
            if native_contact_map[i][j] == 1:
                dist = pow(pow(current_cor[i][0] - current_cor[j][0], 2) + pow(current_cor[i][1] - current_cor[j][1], 2)
                    + pow(current_cor[i][2] - current_cor[j][2], 2), 0.5)
                if dist <= sdist * native_distance_map[i][j]:
                    current_contact_num += 1
    return current_contact_num / native_contact_num

###### END calculate native contact fraction ######

###### Q mod filter ######
def calc_Q_mod(Q_ts):
    edges = np.arange(0, 1.02, 0.02)
    (N, be) = np.histogram(Q_ts, bins=edges)
    idx = np.argwhere(N == np.max(N))[0]
    Q_mod = (edges[idx]+edges[idx+1])/2
    return Q_mod
###### END Q mod filter ######

############## MAIN #################
psffile = sys.argv[1]
corfile = sys.argv[2]
prmfile = sys.argv[3]
temp = float(sys.argv[4])
ppn = sys.argv[5]
outname = sys.argv[6]
rand = int(sys.argv[7])
vecfile = sys.argv[8]
sim_step = int(sys.argv[9])
secondary_structure_def = sys.argv[10]
Q_threshold = float(sys.argv[11])
native_cor = sys.argv[12]
if len(sys.argv) == 14:
    use_gpu = int(sys.argv[13])
else:
    use_gpu = -1
cpfile = outname+'.ncrst'

timestep = 0.015*picoseconds
fbsolu = 0.05/picosecond
temp = temp*kelvin
nsteps_save = 5000
half_window = 100

dist_cutoff = 8 # distance cutoff for finding native contact
sdist = 1.2 # multiple factor of native distance to determine native contact in trajectory
fold_nframe = 100 # number of frames to determine folding status

current_step = 0
folding_tag = 0
if_restart = 0
Q_list = []
if os.path.exists(outname+'.out'):
    ff = open(outname+'.out')
    lines = ff.readlines()
    last_line = lines[-1].strip()
    ff.close()
    if last_line.startswith('Done'):
        print('All Done.')
        exit()
    elif not last_line.startswith('Time') and os.path.getsize(outname+'.out') != 0:
        if_restart = 1
        current_step = int(last_line.split()[1])
        Q_list = [float(l.strip().split()[2]) for l in lines[-min(2*half_window+1, len(lines)-1):]]
        q_list = [l.strip().split()[3] for l in lines[-min(fold_nframe, len(lines)-1):]]
        for i in range(1,len(q_list)+1):
            if q_list[-i] != 'NA':
                if float(q_list[-i]) >= Q_threshold:
                    folding_tag += 1
                else:
                    break
            else:
                break
else:
    f = open(outname+'.out', 'w')
    f.write('%10s %20s %10s %10s\n'%('Time(ns)', 'Steps', 'Q_tot', 'Q_mod'))
    f.close()

### contact map and distance map for start structure ###
native_cor = CharmmCrdFile(native_cor)
native_cor = native_cor.positions.value_in_unit(angstrom)
native_contact_map = [[0 for j in range(len(native_cor))] for i in range(len(native_cor))]
native_distance_map = [[0 for j in range(len(native_cor))] for i in range(len(native_cor))]
sec_strc_def = []
sec_def_object = open(secondary_structure_def,'r')
for line in sec_def_object.readlines():
    line = line.strip()
    if line != '':
        words = line.split()
        sec_strc_def.append([int(words[1]), int(words[2])])
for i in range(len(native_cor)-4):
    for j in range(i+4, len(native_cor)):
        dist = pow(pow(native_cor[i][0] - native_cor[j][0], 2) + pow(native_cor[i][1] - native_cor[j][1], 2)
            + pow(native_cor[i][2] - native_cor[j][2], 2), 0.5)
        if dist <= dist_cutoff:
            native_contact_map[i][j] = 1
            native_distance_map[i][j] = dist
            native_contact_map[j][i] = 1
            native_distance_map[j][i] = dist

sec_def_object.close()
native_contact_num = 0
for i in range(len(native_cor)-4):
    tag_i = 0
    for rs in sec_strc_def:
        if i >= rs[0]-1 and i <= rs[1]-1:
            tag_i = 1
            break
    if tag_i == 0:
        continue
    for j in range(i+4, len(native_cor)):
        tag_j = 0
        for rs in sec_strc_def:
            if j >= rs[0]-1 and j <= rs[1]-1:
                tag_j = 1
                break
        if tag_j == 0:
            continue
        if native_contact_map[i][j] == 1:
            native_contact_num += 1
### END contact map and distance map for start structure ###

psf = CharmmPsfFile(psffile)
psf_pmd = pmd.load_file(psffile)
cor = CharmmCrdFile(corfile)
forcefield = ForceField(prmfile)
top = psf.topology
# re-name residues that are changed by openmm
for resid, res in enumerate(top.residues()):
    if res.name != psf_pmd.residues[resid].name:
        res.name = psf_pmd.residues[resid].name
templete_map = {}
for chain in top.chains():
    for res in chain.residues():
        templete_map[res] = res.name
system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=2.0*nanometer, 
        constraints=AllBonds, removeCMMotion=False, ignoreExternalBonds=True, 
        residueTemplates=templete_map)
custom_nb_force = system.getForce(4)
custom_nb_force.setUseSwitchingFunction(True)
custom_nb_force.setSwitchingDistance(1.8*nanometer)
integrator = LangevinIntegrator(temp, fbsolu, timestep)
integrator.setConstraintTolerance(0.00001)
integrator.setRandomNumberSeed(rand)

# prepare simulation
if use_gpu == -1:
    properties = {'Threads': ppn}
    platform = Platform.getPlatformByName('CPU')
else:
    #dev_index = dev_index_list[int(multiprocessing.current_process().name.split('-')[-1])-1]
    dev_index = use_gpu
    properties = {'CudaPrecision': 'mixed'}
    properties["DeviceIndex"] = "%d"%(dev_index);
    platform = Platform.getPlatformByName('CUDA')

# Attempt of creating the simulation object (sometimes fail due to CUDA environment)
i_attempt = 0
while True:
    try:
        simulation = Simulation(top, system, integrator, platform, properties)
    except Exception as e:
        print('Error occurred at attempt %d...'%(i_attempt+1))
        traceback.print_exc()
        i_attempt += 1
        continue
    else:
        break

if if_restart != 0:
    try:
        rst = pmd.load_file(cpfile)
        simulation.context.setPositions(rst.coordinates[0]*angstrom)
        simulation.context.setVelocities(rst.velocities[0]*angstrom/picosecond)
    except Exception as e:
        print(e)
        print('Warning: Fail to load checkpoint, use the last frame and random velocity instead.')
        dcd_traj = mdt.load(outname+'.dcd', top=psffile)
        dcd_traj[-1].save(outname+'.pdb')
        current_cor = PDBFile(outname+'.pdb')
        os.system('rm -f '+outname+'.pdb')
        simulation.context.setPositions(current_cor.getPositions())
        simulation.context.setVelocitiesToTemperature(temp)
else:
    vec = []
    fo = open(vecfile, 'r')
    for line in fo:
        vec_list = line.strip().split()
        vec_quantity = []
        for v in vec_list:
            vec_quantity.append(float(v)*angstroms/picoseconds)
        vec.append(vec_quantity)
    simulation.context.setPositions(cor.positions)
    simulation.context.setVelocities(vec)

# append reporters
simulation.reporters = []
simulation.reporters.append(pmd.openmm.reporters.RestartReporter(cpfile, nsteps_save, netcdf=True))
if if_restart != 0:
    simulation.reporters.append(DCDReporter(outname+'.dcd', nsteps_save, append=True))
else:
    simulation.reporters.append(DCDReporter(outname+'.dcd', nsteps_save, append=False))

# run production simulation
start_time = time.time()
nframe = 0
while True:
    simulation.step(nsteps_save)
    nframe += 1
    if if_restart == 0:
        time_id = nsteps_save * nframe * timestep.value_in_unit(nanosecond)
        step_id = nsteps_save * nframe
    else:
        time_id = (current_step + nsteps_save * nframe) * timestep.value_in_unit(nanosecond)
        step_id = current_step + nsteps_save * nframe
    current_cor = simulation.context.getState(getPositions=True).getPositions().value_in_unit(angstrom)
    Q = calc_Q(current_cor)
    if len(Q_list) < 2*half_window+1:
        Q_list.append(Q)
    else:
        Q_list.pop(0)
        Q_list.append(Q)
        
    if len(Q_list) == 2*half_window+1:
        Q_mod = calc_Q_mod(np.array(Q_list))
        if Q_mod >= Q_threshold:
            folding_tag += 1
        else:
            folding_tag = 0
        Q_mod = '%10.4f'%Q_mod
    else:
        Q_mod = '%10s'%'NA'
    f = open(outname+'.out', 'a')
    f.write('%10.3f %20d %10.4f %s\n'%(time_id, step_id, Q, Q_mod))
    f.close()
    end_time = time.time()
    time_cost = convert_time(end_time - start_time)
    if folding_tag == fold_nframe:
        f = open(outname+'.out', 'a')
        f.write('Done.\n')
        f.close()
        break
    elif step_id >= sim_step:
        break
