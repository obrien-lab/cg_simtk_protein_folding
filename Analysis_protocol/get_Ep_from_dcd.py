#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr
import getopt, os, time, multiprocessing, random, math
import parmed as pmd
import mdtraj

usage = 'Usage: python get_Ep_from_dcd.py\n' \
        '              --psffile | -p <CG.psf> Charmm psf file for CG model\n'\
        '              --dcdfile | -c <CG.dcd> Charmm dcd file for CG model\n'\
        '              --xmlfile | -x <CG.xml> OpenMM force filed xml file for CG model\n'\
        '              --outfile | -o <Ep.dat> Output file name\n'\
        '              --mask | -m <MASK> Amber type atom mask for selection\n'

def my_progress_bar(total_step, current_step, used_time):
	barLength = 50
	if current_step == 0:
		elapsed_time = 0
	else:
		elapsed_time = used_time / current_step * (total_step - current_step)
	elapsed_time = convert_time(elapsed_time)
	progress = current_step / total_step
	block = int(round(barLength*progress))
	text = '\r[%s] %.1f%% %s'%('#'*block + '-'*(barLength - block), progress*100, elapsed_time)
	stdout.write(text)
	stdout.flush()

def convert_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return ("%d:%02d:%02d"%(h, m, s))

psffile = None
dcdfile = None
xmlfile = None
outfile = 'Ep.dat'
mask = None

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hp:c:x:o:m:",["psffile=", "dcdfile=", "xmlfile=", "outfile=", "mask="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-p", "--psffile"):
        psffile = arg
    elif opt in ("-c", "--dcdfile"):
        dcdfile = arg
    elif opt in ("-x", "--xmlfile"):
        xmlfile = arg
    elif opt in ("-o", "--outfile"):
        outfile = arg
    elif opt in ("-m", "--mask"):
        mask = arg

timestep = 0.015*picoseconds
fbsolu = 0.05/picosecond
temp = 310*kelvin
nonbond_cutoff = 2.0*nanometer
switch_cutoff = 1.8*nanometer

traj = mdtraj.load_dcd(dcdfile, top=psffile)
print('Total number of frames: %d'%traj.n_frames)

psf = pmd.charmm.psf.CharmmPsfFile(psffile)
forcefield = ForceField(xmlfile)

pos = traj.openmm_positions(0)
psf.positions = pos
new_psf = psf[mask]

new_psf.save('tmp.psf', overwrite=True)
new_psf = CharmmPsfFile('tmp.psf')
psf_pmd = pmd.load_file('tmp.psf')
os.system('rm -f tmp.psf')
top = new_psf.topology

# re-name residues that are changed by openmm
for resid, res in enumerate(top.residues()):
    if res.name != psf_pmd.residues[resid].name:
        res.name = psf_pmd.residues[resid].name

template_map = {}
for chain in top.chains():
    for res in chain.residues():
        template_map[res] = res.name
system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
	nonbondedCutoff=nonbond_cutoff, switchDistance=switch_cutoff, 
	constraints=AllBonds, removeCMMotion=True, ignoreExternalBonds=True,
	residueTemplates=template_map)

# must set to use switching function explicitly for CG Custom Nonbond Force #
custom_nb_force = system.getForce(4)
custom_nb_force.setUseSwitchingFunction(True)
custom_nb_force.setSwitchingDistance(switch_cutoff)
# End set to use switching function explicitly for CG Custom Nonbond Force #
integrator = LangevinIntegrator(temp, fbsolu, timestep)
integrator.setConstraintTolerance(0.00001)

# prepare simulation
platform = Platform.getPlatformByName('CPU')
properties = {'Threads': '1'}
#properties = {'Precision': 'double'}
simulation = Simulation(top, system, integrator, platform, properties)

fo = open(outfile, 'w')
fo.close()

my_progress_bar(traj.n_frames, 0, 0)

start_time = time.time()

for frame in range(traj.n_frames):
	pos = traj.openmm_positions(frame)
	psf.positions = pos
	new_psf_1 = psf[mask]
	pos = new_psf_1.positions

	simulation.context.setPositions(pos)

	energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
	fo = open(outfile, 'a')
	fo.write('%.6f\n'%energy)
	fo.close()
	end_time = time.time()
	my_progress_bar(traj.n_frames, frame+1, end_time-start_time)

print('\n')
