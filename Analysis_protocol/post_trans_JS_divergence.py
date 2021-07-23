#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np
from scipy.stats import entropy

################################# Arguments ###################################
# Default parameters
dt = 0.015/1000
nsave = 5000
alpha = 4331293.0
n_window = 200
n_traj = 100
num_rep = 10
mutant_type_list = ['fast', 'slow']
msm_data_file = './msm_data.npz'
state_type = 'micro' # can be 'micro' or 'meta'

# read control file
ctrlfile = ''

usage = '''post_trans_JS_divergence.py
  --ctrlfile | -f <jsd.ctrl> Control file for JSD estimation.
  [--help | -h] Print help information
'''

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
        if line.startswith('dt'):
            words = line.split('=')
            dt = float(words[1].strip())
            continue
        if line.startswith('nsave'):
            words = line.split('=')
            nsave = int(words[1].strip())
            continue
        if line.startswith('alpha'):
            words = line.split('=')
            alpha = float(words[1].strip())
            continue
        if line.startswith('n_window'):
            words = line.split('=')
            n_window = int(words[1].strip())
            continue
        if line.startswith('n_traj'):
            words = line.split('=')
            n_traj = int(words[1].strip())
            continue
        if line.startswith('num_rep'):
            words = line.split('=')
            num_rep = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('msm_data_file'):
            words = line.split('=')
            msm_data_file = words[1].strip()
            continue
        if line.startswith('state_type'):
            words = line.split('=')
            state_type = words[1].strip()
            continue
finally:
     file_object.close()
     
dt = dt*nsave*alpha/1e9 # in seconds

################################# Functions ###################################


################################## MAIN #######################################
# Read post-trans dtrajs
if state_type == 'micro':
    dtrajs = np.load(msm_data_file, allow_pickle=True)['dtrajs']
elif state_type == 'meta':
    dtrajs = np.load(msm_data_file, allow_pickle=True)['meta_dtrajs']
else:
    print('Error: state_type can only be "micro" or "meta"')
    sys.exit()

meta_set = np.load(msm_data_file, allow_pickle=True)['meta_set']

n_states = 0
for i in range(len(dtrajs)):
    if np.max(dtrajs[i]) > n_states:
        n_states = np.max(dtrajs[i])
n_states += 1
print('Number of %sstates: %d'%(state_type, n_states))

max_T_len = 0
for i in range(len(dtrajs)):
    if dtrajs[i].shape[0] > max_T_len:
        max_T_len = dtrajs[i].shape[0]

for dtraj in dtrajs:
    (N, be) = np.histogram(dtraj[-n_window:], bins=np.arange(-0.5, n_states, 1))
    dtraj_last = np.argwhere(N == np.max(N))[0][0]
    

mtype2trajid = [np.arange(n_traj*num_rep*i_ax, n_traj*num_rep*(i_ax+1)).astype(int) for i_ax, mutant_type in enumerate(mutant_type_list)]

# analysis MSM for each mutant
P_list = []
for i_ax, mutant_type in enumerate(mutant_type_list):
    dtrajs_0 = dtrajs[mtype2trajid[i_ax]]
    P_list_0 = np.zeros((max_T_len, n_states))
    for i in range(len(dtrajs_0)):
        (N, be) = np.histogram(dtrajs_0[i][-n_window:], bins=np.arange(-0.5, n_states, 1))
        dtraj_last = np.argwhere(N == np.max(N))[0][0]
        for j in range(max_T_len):
            if j >= len(dtrajs_0[i]): 
                state_0 = dtraj_last
            else:
                state_0 = dtrajs_0[i][j]
            P_list_0[j,state_0] += 1
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
fo.write('%10s '%('Time(s)'))
for j in range(JS_list.shape[1]-1):
    fo.write('%10s '%('P%d'%(j+1)))
fo.write('%10s\n'%('JSD'))
for i in range(max_T_len):
    fo.write('%10.4f '%((i+1)*dt))
    for j in range(JS_list.shape[1]):
        fo.write('%10.4f '%(JS_list[i,j]))
    fo.write('\n')
fo.close()
