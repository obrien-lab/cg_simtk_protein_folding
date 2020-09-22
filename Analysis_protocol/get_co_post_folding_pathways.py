#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np

################################# Arguments ###################################
# Default parameters
n_traj = 100
num_rep = 10
mutant_type_list = ['fast', 'slow']
n_window = 200
co_msm_data_file = ''
post_msm_data_file = ''

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
        if line.startswith('num_rep'):
            words = line.split('=')
            num_rep = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('n_window'):
            words = line.split('=')
            n_window = int(words[1].strip())
            continue
        if line.startswith('co_msm_data_file'):
            words = line.split('=')
            co_msm_data_file = words[1].strip()
            continue
        if line.startswith('post_msm_data_file'):
            words = line.split('=')
            post_msm_data_file = words[1].strip()
            continue
finally:
     file_object.close()

################################# Functions ###################################
def get_pathways(meta_dtrajs, start_states, end_states):
    pathways = {}
    states_on_pathway = []
    start_states = [str(s+1) for s in start_states]
    end_states = [str(s+1) for s in end_states]

    for traj_idx, md in enumerate(meta_dtrajs):
        path = []
        for idx, mdi in enumerate(md):
            if str(mdi+1) in start_states:
                path.append(str(mdi+1))
                break
        if idx == len(md)-1:
            continue
            
        for mdi in md[idx+1:]:
            tag_find = False
            for pi in range(len(path)):
                if path[pi] == str(mdi+1):
                    path = path[0:pi+1]
                    tag_find = True
                    break
            if not tag_find:
                path.append(str(mdi+1))
        
        if path[-1] not in end_states:
            continue      
        
        for p in path:
            if not int(p) in states_on_pathway:
                states_on_pathway.append(int(p))
        
        path = ' -> '.join(path)
        if not path in pathways.keys():
            pathways[path] = 1
        else:
            pathways[path] += 1
    
    tot_num = 0
    for path in pathways.keys():
        tot_num += pathways[path]
    for path in pathways.keys():
        pathways[path] /= tot_num

    sort_pathways = sorted(pathways.items(), key=lambda x: x[1], reverse=True)
    
    states_on_pathway.sort()
    return [sort_pathways, states_on_pathway]

################################## MAIN #######################################
# Read co-trans metastable dtrajs
co_meta_dtrajs = np.load(co_msm_data_file, allow_pickle=True)['meta_dtrajs']
co_n_states = 0
for i in range(len(co_meta_dtrajs)):
    if np.max(co_meta_dtrajs[i]) > co_n_states:
        co_n_states = np.max(co_meta_dtrajs[i])
co_n_states += 1
# Read post-trans metastable dtrajs
post_meta_dtrajs = np.load(post_msm_data_file, allow_pickle=True)['meta_dtrajs']
max_T_len = 0
post_n_states = 0
for i in range(len(post_meta_dtrajs)):
    if post_meta_dtrajs[i].shape[0] > max_T_len:
        max_T_len = post_meta_dtrajs[i].shape[0]
    if np.max(post_meta_dtrajs[i]) > post_n_states:
        post_n_states = np.max(post_meta_dtrajs[i])
post_n_states += 1
post_meta_dtrajs_extended = []
for i in range(len(post_meta_dtrajs)):
    (N, be) = np.histogram(post_meta_dtrajs[i][-n_window:], bins=np.arange(-0.5, post_n_states, 1))
    meta_dtraj_last = np.argwhere(N == np.max(N))[0][0]
    mde = []
    for j in range(max_T_len):
        if j >= len(post_meta_dtrajs[i]): 
            state_0 = meta_dtraj_last
        else:
            state_0 = post_meta_dtrajs[i][j]
        mde.append(state_0)
    post_meta_dtrajs_extended.append(mde)

print('Total number of states: %d'%(co_n_states + post_n_states))

# combine dtrajs
meta_dtrajs = []
for i in range(len(co_meta_dtrajs)):
    for j in range(num_rep):
        md = np.hstack((co_meta_dtrajs[i], post_meta_dtrajs_extended[i*num_rep+j] + co_n_states))
        meta_dtrajs.append(md)
meta_dtrajs = np.array(meta_dtrajs)

mtype2trajid = [np.arange(n_traj*num_rep*i_ax, n_traj*num_rep*(i_ax+1)).astype(int) for i_ax, mutant_type in enumerate(mutant_type_list)]

# analysis MSM for each mutant
fo = open('pathways.dat', 'w')
for i_ax, mutant_type in enumerate(mutant_type_list):    
    # flux analysis
    A = [0]
    B = [md[-1] for md in meta_dtrajs[mtype2trajid[i_ax]]]
    B = np.unique(B)
    
    [pathways, state_on_pathway] = get_pathways(meta_dtrajs[mtype2trajid[i_ax]], A, B)
    
    fo.write('%s pathways:\n'%mutant_type)
    fo.write('%-12s %s\n'%('percentage','path'))
    fo.write('-------------------------------------\n')
    for path in pathways:
        fo.write('%-12s %-s\n'%('%.2f %%'%(path[1] * 100), path[0]))
    fo.write('\n')
