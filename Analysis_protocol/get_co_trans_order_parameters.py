#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
import parmed as pmd
import mdtraj as mdt

################################# Arguments ###################################
# Default parameters
n_traj = 100
mutant_type_list = ['fast', 'slow']
prefix_dir = ''

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
        if line.startswith('n_traj'):
            words = line.split('=')
            n_traj = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('prefix_dir'):
            words = line.split('=')
            prefix_dir = words[1].strip()
            continue
finally:
     file_object.close()

################################# Functions ###################################    
def get_Qbb_T():
    global n_traj, co_traj_dir, mutant_type
    print('QBB:')
    qbb_list = [];
    T_list = [];
    for i in range(n_traj):
        fo = open(co_traj_dir+'qbb_raw_'+str(i+1)+'_t.dat', 'r')
        C = fo.readlines()
        fo.close()
        C = [C[k].strip().split() for k in range(1,len(C))]
        C = np.array(C, dtype=np.float32)
        idx = np.where(C[:,2] < 4)
        Q_ts = C[idx[0],-1].reshape((len(C[idx[0],-1]),1))
        T = C[idx[0],:3]
        qbb_list.append(Q_ts)
        T_list.append(T)
        print('  %s Traj %d Done'%(mutant_type, i+1)) 
    np.save('%s_QBB.npy'%mutant_type, qbb_list)
    np.save('%s_T.npy'%mutant_type, T_list)
    
def get_entanglement():
    global n_traj, co_traj_dir, mutant_type
    print('G:')
    G_list = [];
    for i in range(n_traj):
        fo = open(co_dir+'analysis/G_full_vs_T/G_'+str(i+1)+'_t.dat', 'r')
        C = fo.readlines()
        fo.close()
        C = [C[k].strip().split() for k in range(1,len(C))]
        C = np.array(C)
        idx = np.where(C[:,2].astype(int) < 4)
        G_ts = C[idx[0],4:].astype(np.float32)
        G_list.append(G_ts)
        print('  %s Traj %d Done'%(mutant_type, i+1))  
    np.save('%s_Entanglement.npy'%mutant_type, G_list)
    
################################## MAIN #######################################
for mutant_type in mutant_type_list:
    co_dir = prefix_dir+'/continuous_synthesis/'+mutant_type+('/1-%d/'%(n_traj))
    co_traj_dir = co_dir+'analysis/qbb_full_vs_T/'

    # Get Qbb trajectory
    get_Qbb_T()

    # Get Entanglement trajectory
    get_entanglement()
