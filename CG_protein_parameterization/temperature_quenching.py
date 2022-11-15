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
import getopt, os, time, multiprocessing, random, math
import parmed as pmd

usage = '\nUsage: python temperature_quenching.py\n' \
        '       --ctrlfile | -f <TQ.ctrl> Control file for temperature quenching\n'\
        '       [--help | -h] Print this information\n\n'\
        ' Example of cntrol file:\n'\
        '  tpn = 20\n'\
        '  ppn = 1\n'\
        '  num_traj = 500\n'\
        '  psf = setup/1shf_clean_ca.psf\n'\
        '  top = setup/1shf_clean_ca.top\n'\
        '  param = setup/1shf_clean_nscal1.4_fnn1_go_bt.prmn\n'\
        '  starting_strucs = setup/1shf_clean_ca.cor\n'\
        '  nsteps_equil = 40000000\n'\
        '  temp_equil = 1000\n'\
        '  temp_prod = 310\n'\
        '  Q_threshold = 0.75\n'\
        '  secondary_structure_def = setup/secondary_struc_defs.txt\n'\
        '  log_file = info.log\n'\
        '  restart = 0\n'

###### run Langevin Dynamics ######
def run_TQ_LD(index, rand):
    global nsteps_equil, temp_equil, temp_prod, tag_restart_equil, tag_restart_prod
    global psf, psf_pmd, forcefield, templete_map, start_cor, properties, current_step
    global nsteps_save, timestep, Q_threshold, fold_nframe, folding_array

    fbsolu = 0.05/picosecond
    nonbond_cutoff = 2.0*nanometer
    switch_cutoff = 1.8*nanometer
    constraint_tolerance = 0.00001
    top = psf.topology
    # re-name residues that are changed by openmm
    for resid, res in enumerate(top.residues()):
        if res.name != psf_pmd.residues[resid].name:
            res.name = psf_pmd.residues[resid].name
    platform = Platform.getPlatformByName('CPU')
    system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=nonbond_cutoff, 
        constraints=AllBonds, removeCMMotion=False, ignoreExternalBonds=True,
        residueTemplates=templete_map)
    
    custom_nb_force = system.getForce(4)
    custom_nb_force.setUseSwitchingFunction(True)
    custom_nb_force.setSwitchingDistance(switch_cutoff)
    
    if tag_restart_prod[index-1] == 1:
        integrator = LangevinIntegrator(temp_prod, fbsolu, timestep)
        integrator.setConstraintTolerance(constraint_tolerance)
        simulation = Simulation(top, system, integrator, platform, properties)
        with open('traj/'+str(index)+'.chk', 'rb') as f: simulation.context.loadCheckpoint(f.read())
        simulation.reporters = []
        simulation.reporters.append(CheckpointReporter('traj/'+str(index)+'.chk', nsteps_save))
        if current_step[index-1] == 0:
            simulation.reporters.append(DCDReporter('traj/'+str(index)+'_prod.dcd', nsteps_save, append=False))
        else:
            simulation.reporters.append(DCDReporter('traj/'+str(index)+'_prod.dcd', nsteps_save, append=True))
    else:
        integrator = LangevinIntegrator(temp_equil, fbsolu, timestep)
        integrator.setConstraintTolerance(constraint_tolerance)
        if tag_restart_equil[index-1] == 1:
            simulation = Simulation(top, system, integrator, platform, properties)
            with open('traj/'+str(index)+'.chk', 'rb') as f: simulation.context.loadCheckpoint(f.read())
        else:
            integrator.setRandomNumberSeed(rand)
            simulation = Simulation(top, system, integrator, platform, properties)
            simulation.context.setPositions(start_cor.positions)
            simulation.context.setVelocitiesToTemperature(temp_prod) 
        simulation.reporters = []
        simulation.reporters.append(StateDataReporter('output/'+str(index)+'_equil.out', nsteps_save, step=True,
            potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
            speed=True, totalSteps=nsteps_equil-current_step[index-1], separator='\t'))
        simulation.reporters.append(CheckpointReporter('traj/'+str(index)+'.chk', nsteps_save))
        simulation.reporters.append(PDBReporter('traj/'+str(index)+'_equil.pdb', nsteps_equil-current_step[index-1]))
        simulation.step(nsteps_equil - current_step[index-1])
        pdb = pmd.load_file('traj/'+str(index)+'_equil.pdb')
        pdb.save('traj/'+str(index)+'_equil.cor', format='charmmcrd', overwrite=True)
        os.remove('traj/'+str(index)+'_equil.pdb')

        integrator = LangevinIntegrator(temp_prod, fbsolu, timestep)
        integrator.setConstraintTolerance(constraint_tolerance)
        simulation = Simulation(top, system, integrator, platform, properties)
        with open('traj/'+str(index)+'.chk', 'rb') as f: simulation.context.loadCheckpoint(f.read())
        simulation.reporters = []
        simulation.reporters.append(CheckpointReporter('traj/'+str(index)+'.chk', nsteps_save))
        simulation.reporters.append(DCDReporter('traj/'+str(index)+'_prod.dcd', nsteps_save, append=False))
    
    folding_tag = folding_array[index-1]
    nframe = 0
    while folding_tag == 0:
        simulation.step(nsteps_save)
        nframe += 1
        if tag_restart_prod[index-1] == 0:
            time_id = nsteps_save * nframe * timestep.value_in_unit(nanosecond)
            step_id = nsteps_save * nframe
        else:
            time_id = (current_step[index-1] + nsteps_save * nframe) * timestep.value_in_unit(nanosecond)
            step_id = current_step[index-1] + nsteps_save * nframe
        current_cor = simulation.context.getState(getPositions=True).getPositions().value_in_unit(angstrom)
        Q = calc_Q(current_cor)
        f = open('output/'+str(index)+'_prod.out', 'a')
        f.write('%10.3f %20d %10.3f\n'%(time_id, step_id, Q))
        f.close()
        results = updat_Q(index)
        folding_tag = results[0]

###### END run Langevin Dynamics ######

###### convert time seconds to hours ######
def convert_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return ("%d:%02d:%02d"%(h, m, s))
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

###### update Q_array and folding_array ######
def updat_Q(index):
    global Q_threshold, fold_nframe
    q_array = []
    last_line_list = os.popen('tail -n '+str(fold_nframe)+' output/'+str(index)+'_prod.out').readlines()
    if len(last_line_list) == fold_nframe:
        tag = 1
    else:
        tag = 0
    for line in last_line_list:
        line = line.strip()
        words = line.split()
        q_array.append(float(words[2]))
        if float(words[2]) < Q_threshold:
            tag = 0
    return [tag, q_array]
###### END update Q_array and folding_array ######

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

tpn = 20 # total number of processors
ppn = 1 # number of processors for each trajectory
num_traj = 500 # total number of trajectories
nsteps_equil = 4000000 # number of steps in equilibrium simulation for each exchange, 60 ns
temp_equil = 1000.0 * kelvin # temperature of equilibrium simulation
temp_prod = 310.0 * kelvin # temperature of equilibrium simulation
restart = 0 # flag of whether run restart
log_file = 'info.log' # log file name
psf = '' # Charmm psf file for CG model
top = '' # Charmm top file for CG model
param = '' # Charmm prm file for CG model
starting_strucs = '' # starting structures (Charmm cor file)
nsteps_save = 1000 # save chk, out and dcd each nsteps_save steps
Q_threshold = 0 # threshold for determine folding status
secondary_structure_def = '' # secondary structure defination file
timestep = 0.015*picoseconds # time step for simulation
fold_time = 150*picoseconds # time to determine folding status
dist_cutoff = 8 # distance cutoff for finding native contact
sdist = 1.2 # multiple factor of native distance to determine native contact in trajectory
fold_nframe = int(fold_time/nsteps_save/timestep) # number of frames to determine folding status
sleep_time = 5 # how ofen (seconds) the main process check and write the log file

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
        if line.startswith('tpn'):
            words = line.split('=')
            tpn = int(words[1].strip())
            continue
        if line.startswith('ppn'):
            words = line.split('=')
            ppn = int(words[1].strip())
            continue
        if line.startswith('num_traj'):
            words = line.split('=')
            num_traj = int(words[1].strip())
            continue
        if line.startswith('nsteps_equil'):
            words = line.split('=')
            nsteps_equil = int(words[1].strip())
            continue
        if line.startswith('temp_equil'):
            words = line.split('=')
            temp_equil = float(words[1].strip()) * kelvin
            continue
        if line.startswith('temp_prod'):
            words = line.split('=')
            temp_prod = float(words[1].strip()) * kelvin
            continue
        if line.startswith('restart'):
            words = line.split('=')
            restart = int(words[1].strip())
            continue
        if line.startswith('log_file'):
            words = line.split('=')
            log_file = words[1].strip()
            continue
        if line.startswith('psf'):
            words = line.split('=')
            psf = words[1].strip()
            continue
        if line.startswith('top'):
            words = line.split('=')
            top = words[1].strip()
            continue
        if line.startswith('param'):
            words = line.split('=')
            param = words[1].strip()
            continue
        if line.startswith('starting_strucs'):
            words = line.split('=')
            starting_strucs = words[1].strip()
            continue
        if line.startswith('Q_threshold'):
            words = line.split('=')
            Q_threshold = float(words[1].strip())
            continue
        if line.startswith('secondary_structure_def'):
            words = line.split('=')
            secondary_structure_def = words[1].strip()
            continue
finally:
     file_object.close()

###### Check Control Parameters ######
if tpn == 0:
    print('Error: no processor number specified.')
    sys.exit()
if ppn == 0:
    print('Error: no processor number for each trajectory specified.')
    sys.exit()
if tpn%ppn != 0:
    print('Error: tpn is not the multiple of ppn.')
    sys.exit()
if num_traj == 0:
    print('Error: no trajectory number specified.')
    sys.exit()
if psf == '':
    print('Error: no Charmm psf file specified.')
    sys.exit()
if top == '':
    print('Error: no Charmm top file specified.')
    sys.exit()
if param == '':
    print('Error: no Charmm prm file specified.')
    sys.exit()
if starting_strucs == '':
    print('Error: no starting structures specified.')
    sys.exit()
if restart != 0 and restart != 1:
    print('Error: restart can only be 0 (no restart) or 1 (restart).')
    sys.exit()
if Q_threshold == 0:
    print('Error: no Q_threshold specified.')
    sys.exit()
if secondary_structure_def == '':
    print('Error: no secondary_structure_def specified.')
    sys.exit()
current_step = [0 for i in range(num_traj)]
tag_restart_equil = [0 for i in range(num_traj)]
tag_restart_prod = [0 for i in range(num_traj)]
Q_array = [[0 for j in range(fold_nframe)] for i in range(num_traj)]
folding_array = [0 for i in range(num_traj)]
if restart == 1:
    for i in range(1, num_traj + 1):
        if os.path.exists('output/'+str(i)+'_prod.out'):
            last_line = os.popen('tail -n 1 output/'+str(i)+'_prod.out').read().strip()
            words = last_line.split()
            current_step[i-1] = int(words[1])
            tag_restart_prod[i-1] = 1
            results = updat_Q(i)
            folding_array[i-1] = results[0]
            Q_array[i-1] = results[1]
        elif os.path.exists('output/'+str(i)+'_equil.out'):
            last_line = os.popen('tail -n 1 output/'+str(i)+'_equil.out').read().strip()
            words = last_line.split()
            current_step[i-1] = int(words[1])
            tag_restart_equil[i-1] = 1
        else:
            tag_restart_equil[i-1] = 0
###### END Check Control Parameters ######

###### Setup writing log files ######
os.system('parse_cg_prm.py -t '+top+' -p '+param)
xml_param = param.split('.prm')
xml_param = xml_param[0]+'.xml'

log_file_object = open(log_file,'w')
log_head = ''
log_head += 'Temperature Quenching for CG Model using OpenMM\nAuthor: Yang Jiang; Ed O\'Brien.\n'
log_head += 'Start at '+time.asctime(time.localtime(time.time()))+'\n'
log_head += 'Total number of processors: '+str(tpn)+'\n'
log_head += 'Number of processors for each trajectory: '+str(ppn)+'\n'
log_head += 'Number of trajectories: '+str(num_traj)+'\n'
log_head += 'Number of steps in equilibrium for each trajectory: '+str(nsteps_equil)+'\n'
log_head += 'Temperature of equilibrium simulation: '+str(temp_equil)+'\n'
log_head += 'Temperature of production simulation: '+str(temp_prod)+'\n'
log_head += 'Log file name: '+str(log_file)+'\n'
log_head += 'Charmm psf file: '+str(psf)+'\n'
log_head += 'Charmm top file: '+str(top)+'\n'
log_head += 'Charmm prm file: '+str(param)+'\n'
log_head += 'Secondary structure defination file: '+str(secondary_structure_def)+'\n'
log_head += 'Q threshold: '+str(Q_threshold)+'\n'
if restart == 0:
    log_head += 'No restart requested\n'
else:
    log_head += 'Restart requested\n'
log_head += 'Starting structure: '+starting_strucs+'\n'
log_head += 'File save steps: '+str(nsteps_save)+'\n'
log_head += 'Time step: '+str(timestep)+'\n'

if not os.path.exists('output'):
    os.mkdir('output')
elif restart == 0:
    os.system('rm -rf output')
    os.mkdir('output')
if not os.path.exists('traj'):
    os.mkdir('traj')
elif restart == 0:
    os.system('rm -rf traj')
    os.mkdir('traj')

###### END Setup writing log files ######

###### Temperature Quenching ######
nprocess = int(tpn/ppn)
psf_pmd = pmd.charmm.psf.CharmmPsfFile(psf)
psf = CharmmPsfFile(psf)
forcefield = ForceField(xml_param)
top = psf.topology
start_cor = CharmmCrdFile(starting_strucs)
# re-name residues that are changed by openmm
for resid, res in enumerate(top.residues()):
    if res.name != psf_pmd.residues[resid].name:
        res.name = psf_pmd.residues[resid].name
templete_map = {}
for chain in top.chains():
    for res in chain.residues():
        templete_map[res] = res.name
properties = {'Threads': str(ppn)}

### contact map and distance map for start structure ###
native_cor = start_cor.positions.value_in_unit(angstrom)
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
        else:
            native_contact_map[i][j] = 0
            native_distance_map[i][j] = 0
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

log_head += 'Native contacts found in %s: %d\n'%(starting_strucs, native_contact_num)
log_head += '\n%10s %20s %20s %10s %10s %10s %10s\n'%('SIM_ID', 'START_STEP', 'SIM_STATUS', 'Q_LAST', 'IF_FOLDED' ,'TIME_USED', 'SPEED')
log_output = log_head
start_str = []
for i in range(1, num_traj + 1):
    if tag_restart_prod[i-1] == 0:
        start_str.append('equil@'+str(current_step[i-1]+1))
    else:
        start_str.append('prod@'+str(current_step[i-1]+1))
    log_output += '%10s %20s %20s %10s %10s %10s %10s\n'%(str(i), start_str[i-1], 'wait', '--', '--', '--', '--')
log_file_object.write(log_output)
log_file_object.close()

pool = multiprocessing.Pool(nprocess)
for i in range(1, num_traj + 1):
    rand = random.randint(10,1000000000)
    pool.apply_async(run_TQ_LD, (i, rand))

start_time = [time.time() for i in range(num_traj)]
end_time = [time.time() for i in range(num_traj)]
while True:
    time.sleep(sleep_time)
    log_file_object = open(log_file,'w')
    log_output = log_head
    for i in range(1, num_traj + 1):
        if os.path.exists('output/'+str(i)+'_equil.out'):
            if not os.path.exists('output/'+str(i)+'_prod.out') and tag_restart_prod[i-1] == 0:
                last_line = os.popen('tail -n 1 output/'+str(i)+'_equil.out').read().strip()
                words = last_line.split()
                if len(words) != 6:
                    log_output += '%10s %20s %20s %10s %10s %10s %10s\n'%(str(i), start_str[i-1], 'wait', '--', '--', '--', '--')
                    start_time[i-1] = time.time()
                else:
                    speed = str(float(words[4]))
                    end_time[i-1] = time.time()
                    log_output += '%10s %20s %20s %10s %10s %10s %10s\n'%(str(i), start_str[i-1], 'equil@'+str(current_step[i-1]+int(words[1])), 
                        '--', '--', convert_time(end_time[i-1] - start_time[i-1]), speed)
            elif os.path.exists('output/'+str(i)+'_prod.out'):
                last_line = os.popen('tail -n 1 output/'+str(i)+'_prod.out').read().strip()
                words = last_line.split()
                if words[1] == current_step[i-1]:
                    log_output += '%10s %20s %20s %10s %10s %10s %10s\n'%(str(i), start_str[i-1], 'wait', '--', '--', '--', '--')
                    start_time[i-1] = time.time()
                else:
                    results = updat_Q(i)
                    folding_array[i-1] = results[0]
                    Q_array[i-1] = results[1]
                    if results[0] == 1:
                        speed = (float(words[1]) - current_step[i-1]) / (end_time[i-1] - start_time[i-1])
                        log_output += '%10s %20s %20s %10.3f %10s %10s %10.1f\n'%(str(i), start_str[i-1], 'Done@'+words[1], 
                            Q_array[i-1][-1], str(folding_array[i-1]), convert_time(end_time[i-1] - start_time[i-1]), speed)
                    else:
                        end_time[i-1] = time.time()
                        speed = (float(words[1]) - current_step[i-1]) / (end_time[i-1] - start_time[i-1])
                        log_output += '%10s %20s %20s %10.3f %10s %10s %10.1f\n'%(str(i), start_str[i-1], 'prod@'+words[1], 
                            Q_array[i-1][-1], str(folding_array[i-1]), convert_time(end_time[i-1] - start_time[i-1]), speed)
            else:
                log_output += '%10s %20s %20s %10s %10s %10s %10s\n'%(str(i), start_str[i-1], 'wait', '--', '--', '--', '--')
                start_time[i-1] = time.time()
        else:
            log_output += '%10s %20s %20s %10s %10s %10s %10s\n'%(str(i), start_str[i-1], 'wait', '--', '--', '--', '--')
            start_time[i-1] = time.time()
    log_file_object.write(log_output)
    log_file_object.close()
    
    if len(pool._cache) == 0:
        break

pool.close()
pool.join()
