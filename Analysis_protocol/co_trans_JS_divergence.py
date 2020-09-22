#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np
from scipy.stats import entropy

################################# Arguments ###################################
# Default parameters
n_traj = 100
mutant_type_list = ['fast', 'slow']
co_msm_data_file = './msm_data.npz'
co_TS_data_dir = '../'
num_sample = 1
state_type = 'micro' # can be 'micro' or 'meta'

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
        if line.startswith('co_msm_data_file'):
            words = line.split('=')
            co_msm_data_file = words[1].strip()
            continue
        if line.startswith('num_sample'):
            words = line.split('=')
            num_sample = int(words[1].strip())
            continue
        if line.startswith('state_type'):
            words = line.split('=')
            state_type = words[1].strip()
            continue
finally:
     file_object.close()

################################# Functions ###################################


################################## MAIN #######################################
# Read co-trans dtrajs
if state_type == 'micro':
    dtrajs = np.load(co_msm_data_file, allow_pickle=True)['dtrajs']
elif state_type == 'meta':
    dtrajs = np.load(co_msm_data_file, allow_pickle=True)['meta_dtrajs']
else:
    print('Error: state_type can only be "micro" or "meta"')
    sys.exit()

meta_set = np.load(co_msm_data_file, allow_pickle=True)['meta_set']

n_states = 0
for i in range(len(dtrajs)):
    if np.max(dtrajs[i]) > n_states:
        n_states = np.max(dtrajs[i])
n_states += 1
print('Number of %sstates: %d'%(state_type, n_states))

mtype2trajid = [np.arange(n_traj*i_ax, n_traj*(i_ax+1)).astype(int) for i_ax, mutant_type in enumerate(mutant_type_list)]

# analysis MSM for each mutant
P_list = []
for i_ax, mutant_type in enumerate(mutant_type_list):
    ts = np.load('%s/%s_T.npy'%(co_TS_data_dir, mutant_type), allow_pickle=True)
    max_len = int(ts[0][-1,1])
    dtrajs_0 = dtrajs[mtype2trajid[i_ax]]
    P_list_0 = np.zeros((max_len, n_states))
    for ncl in range(max_len):
        for idx_0, md in enumerate(dtrajs_0):
            idx_list = np.where(ts[idx_0][:,1] == ncl+1)[0]
            idx_list = idx_list[max([-num_sample, -len(idx_list)]):]
            for mmd in md[idx_list]:
                P_list_0[ncl, mmd] += 1
        P_list_0[ncl, :] = P_list_0[ncl, :] / (num_sample*len(dtrajs_0))
    P_list.append(P_list_0)

# Jensen-Shannon divergence
JS_list = []
if state_type == 'micro':
    for ms in meta_set:
        P_list_0 = []
        for i_ax, mutant_type in enumerate(mutant_type_list):
            P = np.copy(P_list[i_ax][:,ms])
            for i in range(len(P)):
                if np.sum(P[i,:]) != 0:
                    P[i,:] = P[i,:] / np.sum(P[i,:])
                else:
                    P[i,0] = 1
            P_list_0.append(P)
        M_0 = 0.5 * (P_list_0[0] + P_list_0[1])
        JS_list.append(0.5 * (entropy(P_list_0[0], M_0, axis=1) + entropy(P_list_0[1], M_0, axis=1)))
M = 0.5 * (P_list[0] + P_list[1])
JS_list.append(0.5 * (entropy(P_list[0], M, axis=1) + entropy(P_list[1], M, axis=1)))
JS_list = np.array(JS_list).T
    
fo = open('JS_div_%s.dat'%state_type, 'w')
fo.write('%10s '%('NC_len'))
for j in range(JS_list.shape[1]-1):
    fo.write('%10s '%('C%d'%(j+1)))
fo.write('%10s\n'%('JSD'))
for i in range(max_len):
    fo.write('%10d '%(i+1))
    for j in range(JS_list.shape[1]):
        fo.write('%10.4f '%(JS_list[i,j]))
    fo.write('\n')
fo.close()
