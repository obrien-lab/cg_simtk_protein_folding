#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr
import getopt, os, time, multiprocessing, random, math
import parmed as pmd
import mdtraj

usage = '\nUsage: python parallel_temperature_REX.py\n' \
        '       --ctrlfile | -f <REX.ctrl> Control file for temperature replica exchange\n'\
        '       [--help | -h] Print this information\n\n'\
        ' Example of cntrol file:\n'\
        '  nwindows = 12\n'\
        '  temps = 280 291 304 319 335 339 344 349 354 364 375 400\n'\
        '  tpn = 12\n'\
        '  psf = setup/1shf_clean_ca.psf\n'\
        '  top = setup/1shf_clean_ca.top\n'\
        '  param = setup/1shf_clean_nscal1.4_fnn1_go_bt.prmn\n'\
        '  nexch_equil = 2\n'\
        '  nsteps_equil = 10000\n'\
        '  nexch_prod = 100000\n'\
        '  nsteps_prod = 3000\n'\
        '  log_file = info.log\n'\
        '  ene_file_prefix = ene\n'\
        '  accp_file_prefix = stats\n'\
        '  starting_strucs_t1 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t2 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t3 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t4 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t5 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t6 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t7 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t8 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t9 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t10 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t11 = setup/1shf_clean_ca.cor\n'\
        '  starting_strucs_t12 = setup/1shf_clean_ca.cor\n'

###### run Langevin Dynamics ######
def run_REX_LD(psf_file, psf, forcefield, templete_map, cor, temp, strtemp, outname, properties, simulation_steps, trajname, rand, window, return_dict):
    start_time = time.time()
    timestep = 0.015*picoseconds
    fbsolu = 0.05/picosecond
    nonbond_cutoff = 2.0*nanometer
    switch_cutoff = 1.8*nanometer
    constraint_tolerance = 0.00001
    psf_pmd = pmd.charmm.psf.CharmmPsfFile(psf_file)
    top = psf.topology
    # re-name residues that are changed by openmm
    for resid, res in enumerate(top.residues()):
        if res.name != psf_pmd.residues[resid].name:
            res.name = psf_pmd.residues[resid].name
    platform = Platform.getPlatformByName('CPU')
    system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=nonbond_cutoff, switchDistance=switch_cutoff, 
        constraints=AllBonds, removeCMMotion=True, ignoreExternalBonds=True,
        residueTemplates=templete_map)
    # must set to use switching function explicitly for CG Custom Nonbond Force #
    custom_nb_force = system.getForce(4)
    custom_nb_force.setUseSwitchingFunction(True)
    custom_nb_force.setSwitchingDistance(switch_cutoff)
    # End set to use switching function explicitly for CG Custom Nonbond Force #
    integrator = LangevinIntegrator(temp, fbsolu, timestep)
    integrator.setConstraintTolerance(constraint_tolerance)
    integrator.setRandomNumberSeed(rand)
    simulation = Simulation(top, system, integrator, platform, properties)
    simulation.context.setPositions(cor.positions)
    simulation.context.setVelocitiesToTemperature(strtemp)
    simulation.reporters = []
    #simulation.reporters.append(PDBReporter(outname+'.pdb', simulation_steps))
    if trajname != '':
        words = outname.split('/')
        words = words[1].split('_')
        if words[1] == '1':
            simulation.reporters.append(DCDReporter(trajname, simulation_steps, append=False))
        else:
            simulation.reporters.append(DCDReporter(trajname, simulation_steps, append=True))
    simulation.step(simulation_steps)
    #pdb = pmd.load_file(outname+'.pdb')
    #pdb.save(outname, format='charmmcrd', overwrite=True)
    #os.remove(outname+'.pdb')
    current_cor = simulation.context.getState(getPositions=True).getPositions()
    psf_pmd.positions = current_cor
    psf_pmd.save(outname, format='charmmcrd', overwrite=True)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    energy = energy.value_in_unit(kilocalorie/mole)
    end_time = time.time()
    return_dict[window] = [energy, end_time-start_time]
###### END run Langevin Dynamics ######

###### Temperature Swap ######
def Temperature_Swap(windows, temps, energy):
    kb = 1.9872/1000 # kcal/mol
    exch_map = []
    for i in range(windows[-1]+1):
        exch_map.append(i)
    for i in range(1, len(windows), 2):
        temp_1 = temps[windows[i-1]].value_in_unit(kelvin)
        temp_2 = temps[windows[i]].value_in_unit(kelvin)
        energy_1 = energy[windows[i-1]]
        energy_2 = energy[windows[i]]
        de = energy_2 - energy_1
        dt = 1 / temp_1 - 1 / temp_2
        exp = math.exp(-dt*de/kb)
        rand = random.random()
        if (exp < 1 and rand < exp) or (exp >= 1): # Accept
            exch_map[windows[i-1]] = windows[i]
            exch_map[windows[i]] = windows[i-1]
        else:
            exch_map[windows[i-1]] = windows[i-1]
            exch_map[windows[i]] = windows[i]
    return exch_map
###### END Temperature Swap ######

###### convert time seconds to hours ######
def convert_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return ("%d:%02d:%02d"%(h, m, s))
###### END convert time seconds to hours ######

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

nwin = 0 # window number
temps = [] # window temperature
tpn = 0 # total number of processors
ppn = 0 # number of processors for each window
nexch_equil = 2 # number of exchanges in equilibrium
nsteps_equil = 10000 # number of steps in equilibrium simulation for each exchange
nexch_prod = 100000 # number of exchanges in production
nsteps_prod = 3000 # number of steps in production simulation for each exchange
restart = 0 # flag of whether run restart
nsteps_start = 1 # start step for restart
log_file = 'info.log' # log file name
ene_file_prefix = 'ene' # energy file prefix
accp_file_prefix = 'stats' # acceptance file prefix
psf = '' # Charmm psf file for CG model
top = '' # Charmm top file for CG model
param = '' # Charmm prm file for CG model
starting_strucs = [] # starting structures (Charmm cor file)

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
        if line.startswith('nwindows'):
            words = line.split()
            nwin = int(words[2])
            continue
        if line.startswith('temps'):
            words = line.split()
            for i in range(2, len(words), 1):
                temps.append(float(words[i]))
            continue
        if line.startswith('tpn'):
            words = line.split()
            tpn = int(words[2])
            continue
        if line.startswith('nexch_equil'):
            words = line.split()
            nexch_equil = int(words[2])
            continue
        if line.startswith('nsteps_equil'):
            words = line.split()
            nsteps_equil = int(words[2])
            continue
        if line.startswith('nexch_prod'):
            words = line.split()
            nexch_prod = int(words[2])
            continue
        if line.startswith('nsteps_prod'):
            words = line.split()
            nsteps_prod = int(words[2])
            continue
        if line.startswith('restart'):
            words = line.split()
            restart = int(words[2])
            continue
        if line.startswith('nsteps_start'):
            words = line.split()
            nsteps_start = int(words[2])
            continue
        if line.startswith('log_file'):
            words = line.split()
            log_file = words[2]
            continue
        if line.startswith('ene_file_prefix'):
            words = line.split()
            ene_file_prefix = words[2]
            continue
        if line.startswith('accp_file_prefix'):
            words = line.split()
            accp_file_prefix = words[2]
            continue
        if line.startswith('psf'):
            words = line.split()
            psf = words[2]
            continue
        if line.startswith('top'):
            words = line.split()
            top = words[2]
            continue
        if line.startswith('param'):
            words = line.split()
            param = words[2]
            continue
        if line.startswith('starting_strucs'):
            words = line.split()
            starting_strucs.append(words[2])
            continue
finally:
     file_object.close()

###### Check Control Parameters ######
if nwin != len(temps):
    print('Error: window number and temperature number mismatch.')
    sys.exit()
if nwin == 0:
    print('Error: no window number specified.')
    sys.exit()
if len(temps) == 0:
    print('Error: no temperature specified.')
    sys.exit()
if tpn == 0:
    print('Error: no processor number specified.')
    sys.exit()
if tpn % nwin == 0:
    ppn = int(tpn / nwin)
else:
    print('Error: total processor number are not multiple of window number.')
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
if len(starting_strucs) == 0 and restart == 0:
    print('Error: no starting structures specified.')
    sys.exit()
if len(starting_strucs) != nwin and restart == 0:
    print('Error: window number and structure number mismatch.')
    sys.exit()
if restart == 1:
    nexch_equil = 0
    starting_strucs = []
    n_frame_traj = []
    n_frame_ene = []
    for i in range(nwin):
        traj = mdtraj.load_dcd('aa'+str(i+1)+'/mc1.dcd', top=psf)
        n_frame_traj.append(traj.n_frames)
        n_frame_ene.append(int(os.popen('cat aa'+str(i+1)+'/'+ene_file_prefix+'_1.log | wc -l').read().strip()))
    min_n_frame_traj = min(n_frame_traj)
    min_n_frame_ene = min(n_frame_ene)
    min_n_frame = min(min_n_frame_traj, min_n_frame_ene)
    if min_n_frame < nsteps_start-1:
        print('Warning: min frame number is less than nsteps_start-1; '+str(min_n_frame)+' < '+str(nsteps_start-1))
        nsteps_start = min_n_frame+1
        print('Adjust nsteps_start = '+str(nsteps_start))
    for i in range(nwin):
        if n_frame_traj[i] >= nsteps_start-1:
            traj = mdtraj.load_dcd('aa'+str(i+1)+'/mc1.dcd', top=psf)
            new_traj = traj[0:nsteps_start-1]
            new_traj.save('aa'+str(i+1)+'/mc1.dcd', force_overwrite=True)
            new_traj[-1].save('aa'+str(i+1)+'/1_'+str(new_traj.n_frames)+'_prod.pdb', force_overwrite=True)
            pdb = pmd.load_file('aa'+str(i+1)+'/1_'+str(new_traj.n_frames)+'_prod.pdb')
            pdb.save('aa'+str(i+1)+'/1_'+str(new_traj.n_frames)+'_prod.cor', format='charmmcrd', overwrite=True)
            os.remove('aa'+str(i+1)+'/1_'+str(new_traj.n_frames)+'_prod.pdb')
        if n_frame_ene[i] > nsteps_start-1:
            f = open('aa'+str(i+1)+'/'+ene_file_prefix+'_1.log', 'r')
            f_list = f.readlines()
            f.close()
            fo = open('aa'+str(i+1)+'/'+ene_file_prefix+'_1.log', 'w')
            for line in range(nsteps_start-1):
                fo.write(f_list[line])
            fo.close()
            #os.system('sed \''+str(nsteps_start)+',$d\' -i aa'+str(i+1)+'/'+ene_file_prefix+'_1.log')
        starting_strucs.append('aa'+str(i+1)+'/1_'+str(nsteps_start-1)+'_prod.cor')
    log_file = log_file.split('.')[0] + '_r_' + str(nsteps_start) + '.log'
    accp_file_prefix = accp_file_prefix + '_r_' + str(nsteps_start)
elif restart == 0:
    nsteps_start = 1
    os.system('rm -rf aa*')
    os.system('rm -rf logs')
else:
    print('Error: restart can only be 0 (no restart) or 1 (restart).')
    sys.exit()
###### END Check Control Parameters ######

###### Setup writing log files ######

os.system('parse_cg_prm.py -t '+top+' -p '+param)
xml_param = param.split('.prm')
xml_param = xml_param[0]+'.xml'

ene_file = ene_file_prefix+'_1.log'
accp_file = accp_file_prefix+'-tswap-1-1.log'

log_file_object = open(log_file,'w')
log_file_object.write('Parallel Temperature Replica Exchange for CG Model using OpenMM\nAuthor: Yang Jiang; Ed O\'Brien.\n')
log_file_object.write('Start at '+time.asctime(time.localtime(time.time()))+'\n')
log_file_object.write('Number of windows: '+str(nwin)+'\n')
log_file_object.write('Temperatures: ')
for i in range(len(temps)):
    log_file_object.write(str(temps[i])+' ')
log_file_object.write('\n')
log_file_object.write('Total number of processors: '+str(tpn)+'\n')
log_file_object.write('Number of processors for each window: '+str(ppn)+'\n')
log_file_object.write('Number of exchanges in equilibrium: '+str(nexch_equil)+'\n')
log_file_object.write('Number of steps in equilibrium for each exchange: '+str(nsteps_equil)+'\n')
log_file_object.write('Number of exchanges in production: '+str(nexch_prod)+'\n')
log_file_object.write('Number of steps in production for each exchange: '+str(nsteps_prod)+'\n')
log_file_object.write('Log file name: '+str(log_file)+'\n')
log_file_object.write('Eenergy file name: '+ene_file+'\n')
log_file_object.write('Acceptance file name: '+accp_file+'\n')
log_file_object.write('Charmm psf file: '+str(psf)+'\n')
log_file_object.write('Charmm top file: '+str(top)+'\n')
log_file_object.write('Charmm prm file: '+str(param)+'\n')
if restart == 0:
    log_file_object.write('No restart requested\n')
else:
    log_file_object.write('Restart requested from '+str(nsteps_start)+'\n')
for i in range(len(starting_strucs)):
    log_file_object.write('Starting structure for T'+str(i+1)+': '+starting_strucs[i]+'\n')
log_file_object.close()

if not os.path.exists('logs'):
    os.mkdir('logs')
accp_file_object = open('logs/'+accp_file,'w')
accp_file_object.close()

for i in range(nwin):
    if not os.path.exists('aa'+str(i+1)):
        os.mkdir('aa'+str(i+1))
###### END Setup writing log files ######

###### Temperature REX ######
energy = []
process_time = []
accp = []
cor_list = []
exch_map = []
window_track = []
nexch = []
#energy_file_object = []
for i in range(nwin):
    temps[i] = temps[i]*kelvin
    cor_list.append(starting_strucs[i])
    exch_map.append(i)
    window_track.append(i)
    accp.append(0)
    nexch.append(0)
psf_file = psf
psf = CharmmPsfFile(psf)
psf_pmd = pmd.charmm.psf.CharmmPsfFile(psf_file)
forcefield = ForceField(xml_param)
top = psf.topology
# re-name residues that are changed by openmm
for resid, res in enumerate(top.residues()):
    if res.name != psf_pmd.residues[resid].name:
        res.name = psf_pmd.residues[resid].name
templete_map = {}
for chain in top.chains():
    for res in chain.residues():
        templete_map[res] = res.name
properties = {'Threads': str(ppn)}

###### equil phase ######
for i in range(nexch_equil):
    energy = []
    process_time = []
    process_pool = []
    start_time = time.time()
    return_dict = multiprocessing.Manager().dict()
    log_file_object = open(log_file,'a')
    log_file_object.write('EQUIL '+str(i+1)+': ')
    for window in range(nwin):
#        if window == nwin-1:
#            log_file_object.write(str(window_track[window]+1)+'| ')
#        else:
#            log_file_object.write(str(window_track[window]+1)+'||| ')
        
        strtemp = temps[exch_map[window]]
        cor = CharmmCrdFile(cor_list[window])
        outname = 'aa'+str(window+1)+'/1_'+str(i+1)+'_equil.cor'
        rand = random.randint(10,1000000000)
        p = multiprocessing.Process(target=run_REX_LD, args=(psf_file, psf, forcefield, templete_map, cor, temps[window], 
            strtemp, outname, properties, nsteps_equil, '', rand, window, return_dict))
        p.daemon = True
        process_pool.append(p)
    for window in range(nwin):
        process_pool[window].start()
    for window in range(nwin):
        process_pool[window].join()
    for window in range(nwin):
        results = return_dict[window]
        energy.append(results[0])
        process_time.append(results[1])
    end_time = time.time()
    cost_time = end_time-start_time
    exch_map = Temperature_Swap(range(i%2, nwin), temps, energy)
    log_file_object.write('%s %s %.2f %.2f %s\n'%(time.strftime('%H:%M:%S', time.localtime(start_time)), 
        time.strftime('%H:%M:%S', time.localtime(end_time)), min(process_time), max(process_time), 
        convert_time((nexch_equil-i-1)*cost_time)))
    log_file_object.close()

    # update cor_list and window_track
    new_window_track = [0 for window in range(nwin)]
    for window in range(nwin):
        cor_list[window] = 'aa'+str(exch_map[window]+1)+'/1_'+str(i+1)+'_equil.cor'
        new_window_track[window] = window_track[exch_map[window]]
    window_track = new_window_track

###### prod phase ######
for i in range(nwin):
    window_track[i] = i
for i in range(nsteps_start-1, nexch_prod):
    energy = []
    process_time = []
    process_pool = []
    start_time = time.time()
    return_dict = multiprocessing.Manager().dict()
    log_file_object = open(log_file,'a')
    log_file_object.write('PROD '+str(i+1)+': ')
    for window in range(nwin):
        if window == nwin-1:
            log_file_object.write(str(window_track[window]+1)+'| ')
        else:
            log_file_object.write(str(window_track[window]+1)+'||| ')
        
        strtemp = temps[exch_map[window]]
        cor = CharmmCrdFile(cor_list[window])
        outname = 'aa'+str(window+1)+'/1_'+str(i+1)+'_prod.cor'
        rand = random.randint(10,1000000000)
        p = multiprocessing.Process(target=run_REX_LD, args=(psf_file, psf, forcefield, templete_map, cor, temps[window], 
            strtemp, outname, properties, nsteps_prod, 'aa'+str(window+1)+'/mc1.dcd', rand, window, return_dict))
        p.daemon = True
        process_pool.append(p)
    for window in range(nwin):
        process_pool[window].start()
    for window in range(nwin):
        process_pool[window].join()
    for window in range(nwin):
        results = return_dict[window]
        energy.append(results[0])
        f = open('aa'+str(window+1)+'/'+ene_file, 'a')
        f.write('%.6f\n'%(results[0]))
        f.close()
        process_time.append(results[1])
    end_time = time.time()
    cost_time = end_time-start_time
    exch_map = Temperature_Swap(range(i%2, nwin), temps, energy)
    log_file_object.write('%s %s %.2f %.2f %s\n'%(time.strftime('%H:%M:%S', time.localtime(start_time)), 
        time.strftime('%H:%M:%S', time.localtime(end_time)), min(process_time), max(process_time), 
        convert_time((nexch_prod-i-1)*cost_time)))
    log_file_object.close()

    # update cor_list and window_track
    new_window_track = [0 for window in range(nwin)]
    for window in range(nwin):
        cor_list[window] = 'aa'+str(exch_map[window]+1)+'/1_'+str(i+1)+'_prod.cor'
        new_window_track[window] = window_track[exch_map[window]]
        if i > 0:
            os.remove('aa'+str(window+1)+'/1_'+str(i)+'_prod.cor')
    window_track = new_window_track

    # print acceptence ratio
    accp_file_object = open('logs/'+accp_file,'a')
    accp_file_object.write('T-EXCHNG %d: '%(i+1))
    if i%2 == 0:
        if nwin%2 == 0:
            for window in range(nwin):
                nexch[window] += 1
        else:
            for window in range(nwin-1):
                nexch[window] += 1
    else:
        if (nwin-1)%2 == 0:
            for window in range(1, nwin):
                nexch[window] += 1
        else:
            for window in range(1, nwin-1):
                nexch[window] += 1
    tag = -1
    for window in range(nwin):     
        if window != exch_map[window]:
            accp[window] += 1
            if window != tag:
                accp_file_object.write('%.1f <=> %.1f '%(temps[window].value_in_unit(kelvin), 
                    temps[exch_map[window]].value_in_unit(kelvin)))
            tag = exch_map[window]
    for window in range(nwin):
    	if nexch[window] == 0:
    		accp_file_object.write('%.2f '%(0.0))
    	else:
    		accp_file_object.write('%.2f '%(accp[window]/nexch[window]))
    accp_file_object.write('\n')
    accp_file_object.close()
