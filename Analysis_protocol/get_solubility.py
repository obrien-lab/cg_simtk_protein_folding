#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np
import mdtraj as mdt
import parmed as pmd

################################# Arguments ###################################
# Default parameters
traj_file = ''
psf_file = ''
native_AA_pdb = ''
num_samples = 5
n_boot = int(1e4)
fasta_seq = ''
disordered_pdb = ''
# This is the trimer interface region should be removed from propensity estimation
offset_sel = ''
agg_prop_file = 'agg_pro.dat'
deg_sel = ':ILE,VAL,LEU,PHE,CYS,MET,ALA,GLY,TRP'

# read control file
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
        if line.startswith('traj_file'):
            words = line.split('=')
            traj_file = words[1].strip()
            continue
        if line.startswith('psf_file'):
            words = line.split('=')
            psf_file = words[1].strip()
            continue
        if line.startswith('native_AA_pdb'):
            words = line.split('=')
            native_AA_pdb = words[1].strip()
            continue
        if line.startswith('num_samples'):
            words = line.split('=')
            num_samples = int(words[1].strip())
            continue
        if line.startswith('n_boot'):
            words = line.split('=')
            n_boot = int(words[1].strip())
            continue
        if line.startswith('fasta_seq'):
            words = line.split('=')
            fasta_seq = words[1].strip()
            continue
        if line.startswith('disordered_pdb'):
            words = line.split('=')
            disordered_pdb = words[1].strip()
            continue
        if line.startswith('offset_sel'):
            words = line.split('=')
            offset_sel = words[1].strip()
            continue
        if line.startswith('agg_prop_file'):
            words = line.split('=')
            agg_prop_file = words[1].strip()
            continue
        if line.startswith('deg_sel'):
            words = line.split('=')
            deg_sel = words[1].strip()
            continue
        
finally:
     file_object.close()

################################# Functions ###################################
def bootstrap(boot_fun, data, n_time):
    idx_list = np.arange(len(data))
    if len(data.shape) == 1:
        boot_stat = np.zeros(n_time)
    elif len(data.shape) == 2:
        boot_stat = np.zeros((n_time, data.shape[1]))
    else:
        print('bootstrap: Can only handle 1 or 2 dimentional data')
        sys.exit()
        
    for i in range(n_time):
        sample_idx_list = np.random.choice(idx_list, len(idx_list))
        if len(data.shape) == 1:
            new_data = data[sample_idx_list]
            boot_stat[i] = boot_fun(new_data)
        elif len(data.shape) == 2:
            new_data = data[sample_idx_list, :]
            for j in range(data.shape[1]):
                boot_stat[i,j] = boot_fun(new_data[:,j])
            
    return boot_stat

################################## MAIN #######################################
struct = pmd.load_file(psf_file)
ext_struct = pmd.load_file(disordered_pdb)
for idx, res in enumerate(ext_struct.residues):
    struct.residues[idx].name = res.name
# Estimate Aggregation site
agg_sel = ':'
f = open('agg_pro.dat')
lines = f.readlines()
f.close()
for li, line in enumerate(lines):
    if line.startswith('HITS'):
        break
for line in lines[li+1:]:
    line = line.strip()
    if line.startswith('>--->'):
        agg_sel += line.split(':')[1]
        break
agg_sel += offset_sel
print('Agg_sel: %s'%agg_sel)
agg_idx = list(pmd.amber.AmberMask(struct, agg_sel).Selected())
# Estimate Degradation site
deg_sel += offset_sel
print('Deg_sel: %s'%deg_sel)
deg_idx = list(pmd.amber.AmberMask(struct, deg_sel).Selected())
# Estimate Chaperone binding site
cb_sel = ':'
lines = os.popen('ChaperISM.py '+fasta_seq+' -qt').readlines()
for li, line in enumerate(lines):
    if line.startswith(' POSITION '):
        break
cb_site_list = []
for line in lines[li+2:len(lines)-1]:
    words = line.strip().split()
    if words[-1] == '*':
        #cb_site_list += list(np.arange(int(words[0])+1, int(words[0])+8))
        cb_site_list.append(int(words[0])+4)
cb_site_list = np.unique(cb_site_list)
cb_site_list = [str(c) for c in cb_site_list]
cb_sel += ','.join(cb_site_list)
cb_sel += offset_sel
print('CB_sel: %s'%cb_sel)
cb_idx = list(pmd.amber.AmberMask(struct, cb_sel).Selected())

# Calculate SASA
traj = mdt.load(traj_file, top=psf_file)
SASA_list = []
for frame_id in range(traj.n_frames):
#for frame_id in range(5):
    SASA_0 = []
    state_id = int(frame_id / num_samples)
    rep_id = frame_id - state_id * num_samples
    os.system('mkdir tmp')
    os.chdir('tmp')
    traj[frame_id].save('tmp.pdb', force_overwrite=True)
    os.system('backmap.py -i '+native_AA_pdb+' -c tmp.pdb > /dev/null 2>&1')
    if not os.path.exists('tmp_rebuilt.pdb'):
        SASA_0 = [np.nan, np.nan, np.nan]
        print('Failed to backmap')
    else:
        struct = mdt.load('tmp_rebuilt.pdb')
        sasa = mdt.shrake_rupley(struct, mode='residue')[0]
        agg_sasa = np.sum(sasa[agg_idx])
        deg_sasa = np.sum(sasa[deg_idx])
        cb_sasa = np.sum(sasa[cb_idx])
        SASA_0 = [agg_sasa, deg_sasa, cb_sasa]
    SASA_list.append(SASA_0)
    os.chdir('../')
    os.system('rm -rf tmp/')
    print('Frame %d done'%frame_id)
struct = mdt.load(disordered_pdb)
sasa = mdt.shrake_rupley(struct, mode='residue')[0]
agg_sasa = np.sum(sasa[agg_idx])
deg_sasa = np.sum(sasa[deg_idx])
cb_sasa = np.sum(sasa[cb_idx])
SASA_list = [[agg_sasa, deg_sasa, cb_sasa]] + SASA_list
SASA_list = np.array(SASA_list)

n_states = int((len(SASA_list)-1) / num_samples)
SASA_avg = np.zeros((n_states+1, SASA_list.shape[1]))
SASA_std = np.zeros((n_states+1, SASA_list.shape[1]))
SASA_avg[0,:] = SASA_list[0,:]
SASA_std[0,:] = np.zeros(SASA_list.shape[1])
for i in range(n_states):
    SASA_avg[i+1,:] = np.mean(SASA_list[i*num_samples+1:(i+1)*num_samples+1,:], axis=0)
    SASA_std[i+1,:] = np.std(SASA_list[i*num_samples+1:(i+1)*num_samples+1,:], axis=0)

fo = open('solubility.dat','w')
fo.write('%6s %10s %10s %10s %10s %10s %10s %10s %15s %10s %15s %10s %15s %10s %10s %15s\n'%(
                'State', 'Agg_SA_avg', 'Agg_SA_std', 'Deg_SA_avg', 'Deg_SA_std', 'CB_SA_avg', 
                'CB_SA_std', 'Agg_Prop', 'Agg_95CI', 'Deg_Prop', 'Deg_95CI', 'CB_Prop', 'CB_95CI', 
                'Tot_Prop', 'Sol_Perc', 'Sol_95CI'))
agg_prop = SASA_avg[0,0]/SASA_avg[-1,0]-1
agg_lb = np.nan
agg_ub = np.nan
deg_prop = SASA_avg[0,1]/SASA_avg[-1,1]-1
deg_lb = np.nan
deg_ub = np.nan
cb_prop = SASA_avg[0,2]/SASA_avg[-1,2]-1
cb_lb = np.nan
cb_ub = np.nan
tot_prop_0 = agg_prop + deg_prop - cb_prop
sol_lb = np.nan
sol_ub = np.nan
fo.write('%6s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.4f %15s %10.4f %15s %10.4f %15s %10.4f %10.4f %15s\n'%(
                'ext', SASA_avg[0,0], SASA_std[0,0], SASA_avg[0,1], SASA_std[0,1], SASA_avg[0,2], SASA_std[0,2], 
                agg_prop, '[%.4f,%.4f]'%(agg_lb, agg_ub), deg_prop, '[%.4f,%.4f]'%(deg_lb, deg_ub), 
                cb_prop, '[%.4f,%.4f]'%(cb_lb, cb_ub), tot_prop_0, 0, '[%.4f,%.4f]'%(sol_lb, sol_ub)))

SASA_boot = []
for i in range(1,n_states+1):
    data = SASA_list[(i-1)*num_samples+1:i*num_samples+1,:]
    SASA_boot.append(bootstrap(np.mean, data, n_boot))

for i in range(1,n_states+1):
    prop_boot = SASA_boot[i-1]/SASA_boot[-1]-1
    prop_boot[prop_boot<0] = 0
    
    agg_prop = SASA_avg[i,0]/SASA_avg[-1,0]-1
    if agg_prop < 0:
        print('Correct state %d aggregation propensity %.4f to be zeros'%(i, agg_prop))
        agg_prop = 0
    agg_lb = np.percentile(prop_boot[:,0], 2.5)
    agg_ub = np.percentile(prop_boot[:,0], 97.5)
        
    deg_prop = SASA_avg[i,1]/SASA_avg[-1,1]-1
    if deg_prop < 0:
        print('Correct state %d degradation propensity %.4f to be zeros'%(i, deg_prop))
        deg_prop = 0
    deg_lb = np.percentile(prop_boot[:,1], 2.5)
    deg_ub = np.percentile(prop_boot[:,1], 97.5)
    
    cb_prop = SASA_avg[i,2]/SASA_avg[-1,2]-1
    if cb_prop < 0:
        print('Correct state %d Hsp70 binding propensity %.4f to be zeros'%(i, cb_prop))
        cb_prop = 0
    cb_lb = np.percentile(prop_boot[:,2], 2.5)
    cb_ub = np.percentile(prop_boot[:,2], 97.5)
    
    tot_prop = agg_prop + deg_prop - cb_prop
    if tot_prop < 0:
        print('Correct state %d total propensity %.4f to be zeros'%(i, tot_prop))
        tot_prop = 0
        
    sol_perc = 1 - tot_prop / tot_prop_0
    sol_boot = 1 - (prop_boot[:,0] + prop_boot[:,1] - prop_boot[:,2]) / tot_prop_0
    sol_lb = np.percentile(sol_boot, 2.5)
    sol_ub = np.percentile(sol_boot, 97.5)
    
    fo.write('%6d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.4f %15s %10.4f %15s %10.4f %15s %10.4f %10.4f %15s\n'%(
                    i, SASA_avg[i,0], SASA_std[i,0], SASA_avg[i,1], SASA_std[i,1], SASA_avg[i,2], SASA_std[i,2], 
                    agg_prop, '[%.4f,%.4f]'%(agg_lb, agg_ub), deg_prop, '[%.4f,%.4f]'%(deg_lb, deg_ub), 
                    cb_prop, '[%.4f,%.4f]'%(cb_lb, cb_ub), tot_prop, sol_perc, '[%.4f,%.4f]'%(sol_lb, sol_ub)))
fo.close()
    

