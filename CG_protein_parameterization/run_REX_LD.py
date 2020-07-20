#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr
import getopt, os, time
import parmed as pmd

usage = 'Usage: python run_REX_LD.py\n' \
        '       --psffile | -p <CG.psf> Charmm psf file\n'\
        '       --corfile | -c <CG.cor> Charmm cor file\n'\
        '       --prmfile | -x <CG.xml> OpenMM xml forcefield file\n'\
        '       --temperature | -t <xxx Kelvin> Simulation temperature in Kelvin\n'\
        '       --strtemp | -b <xxx Kelvin> Start temperature to generate velocity\n'\
        '       --totalstep | -s <xxx Steps> Total simulation steps\n'\
        '       --processnum | -n <xxx Processors> Processor number for this run in parallel\n'\
        '       --outname | -o <xxx.dcd> Trajectory file name\n'\
        '       [--help | -h] Print this information\n'

psffile = None
corfile = None
prmfile = None
temp = None
strtemp = None
simulation_steps = None
np = None
outname = None

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hp:c:x:t:s:n:o:b:",
    	["psffile=", "corfile=", "prmfile=", "temperature=", "totalstep=", "processnum=", "outname=", "strtemp="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-p", "--psffile"):
        psffile = arg
    elif opt in ("-c", "--corfile"):
        corfile = arg
    elif opt in ("-x", "--prmfile"):
        prmfile = arg
    elif opt in ("-t", "--temperature"):
        temp = float(arg)
    elif opt in ("-b", "--strtemp"):
        strtemp = float(arg)
    elif opt in ("-s", "--totalstep"):
        simulation_steps = int(arg)
    elif opt in ("-n", "--processnum"):
        np = arg
    elif opt in ("-o", "--outname"):
        outname = arg


timestep = 0.015*picoseconds
fbsolu = 0.05/picosecond
temp = temp*kelvin

psf = CharmmPsfFile(psffile)
cor = CharmmCrdFile(corfile)
forcefield = ForceField(prmfile)
top = psf.topology
templete_map = {}
for chain in top.chains():
    for res in chain.residues():
        templete_map[res] = res.name
system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=2.0*nanometer, switchDistance=1.8*nanometer, 
        constraints=AllBonds, removeCMMotion=True, ignoreExternalBonds=True, 
        residueTemplates=templete_map)
custom_nb_force = system.getForce(4)
custom_nb_force.setUseSwitchingFunction(True)
custom_nb_force.setSwitchingDistance(1.8*nanometer)

integrator = LangevinIntegrator(temp, fbsolu, timestep)
integrator.setConstraintTolerance(0.00001)

# prepare simulation
platform = Platform.getPlatformByName('CPU')
properties = {'Threads': np}
simulation = Simulation(top, system, integrator, platform, properties)
simulation.context.setPositions(cor.positions)
simulation.context.setVelocitiesToTemperature(strtemp)

# append reporters
simulation.reporters = []
simulation.reporters.append(DCDReporter(outname+'.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=simulation_steps, separator='\t'))

# run production simulation
print('Running Production...')
start_time = time.time()
simulation.step(simulation_steps)
current_cor = simulation.context.getState(getPositions=True).getPositions()
psf_pmd = pmd.load_file(psffile)
psf_pmd.positions = current_cor
psf_pmd.save(outname+'.cor', format='charmmcrd', overwrite=True)
end_time = time.time()
speed = simulation_steps/(end_time-start_time)
print('Speed: %d steps/second'%(int(speed)))
print('Done!')

