#!/usr/bin/env python3
import sys, getopt, math, os, time, multiprocessing, traceback
import numpy as np

################################# Arguments ###################################
# Default values
end_t = 60 # in seconds
dt = 0.015/1000
nsave = 5000
alpha = 4331293.0
n_traj = 100
mutant_type_list = ['fast', 'slow']
n_rep = 10
msm_data_file = './msm_data.npz'
n_boot = 10000
num_proc = 20
num_points_plot = 1000
if_boot = True

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
        if line.startswith('end_t'):
            words = line.split('=')
            end_t = float(words[1].strip())
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
        if line.startswith('n_traj'):
            words = line.split('=')
            n_traj = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('n_rep'):
            words = line.split('=')
            n_rep = int(words[1].strip())
            continue
        if line.startswith('msm_data_file'):
            words = line.split('=')
            msm_data_file = words[1].strip()
            continue
        if line.startswith('n_boot'):
            words = line.split('=')
            n_boot = int(words[1].strip())
            continue
        if line.startswith('num_proc'):
            words = line.split('=')
            num_proc = int(words[1].strip())
            continue
        if line.startswith('num_points_plot'):
            words = line.split('=')
            num_points_plot = int(words[1].strip())
            continue
        if line.startswith('if_boot'):
            words = line.split('=')
            if_boot = int(words[1].strip())
            if if_boot == 1:
                if_boot = True
            elif if_boot == 0:
                if_boot = False
            else:
                print('Error: if_boot can only be either 0 or 1.')
                sys.exit()
            continue
finally:
     file_object.close()

dt = dt*nsave*alpha/1e9 # in seconds

################################# Functions ###################################
def boot_fun(data, n_states, job_id):
    PPT = np.zeros((data.shape[1], n_states))
    for md in data:
        for i in range(data.shape[1]):
            PPT[i,md[i]]+=1
    PPT /= len(data)
    # print('Bootstrapping done for %d'%(job_id+1))
    return PPT

def bootstrap(boot_fun, data, n_boot):
    global num_proc, n_states
    
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    
    idx_list = np.arange(len(data))
    
    start_time = time.time()
    print('start bootsrapping')
    for i in range(n_boot):
        sample_idx_list = np.random.choice(idx_list, len(idx_list))
        new_data = data[sample_idx_list,:]
        pool_list.append(pool.apply_async(boot_fun, (new_data, n_states, i, )))
    pool.close()
    pool.join()
    boot_stat = [p.get() for p in pool_list]
    used_time = time.time() - start_time
    print('%.2fs'%used_time)
    return boot_stat

################################## MAIN #######################################
npzfile = np.load(msm_data_file, allow_pickle=True)

max_T_len = int(np.ceil(end_t/dt))
interval = int(max_T_len / num_points_plot)
sample_idx = [max_T_len-1-i*interval for i in range(int(max_T_len/interval), -1, -1)]
if sample_idx[0] != 0:
    sample_idx = [0] + sample_idx

md_list = npzfile['meta_dtrajs']
n_states = 0
tag_error = False
for i, md in enumerate(md_list):
    if n_states < np.max(md):
        n_states = np.max(md)
    if len(md) < max_T_len:
        print("Error: Traj #%d stopped early"%(i+1))
        tag_error = True
n_states += 1   

if tag_error:
    sys.exit()

si = 0
MSTS_list = []
boot_stat_list = []
for i_ax, mutant_type in enumerate(mutant_type_list):
    ei = si+n_traj*n_rep
    meta_dtrajs = md_list[np.arange(si, ei)]
    meta_dtrajs = meta_dtrajs[:,sample_idx]
    si = ei
    
    # MSTS
    PPT = np.zeros((len(sample_idx), n_states))
    t_span = (np.array(sample_idx)+1)*dt
    for md in meta_dtrajs:
        for i in range(len(sample_idx)):
            PPT[i,md[i]]+=1
    PPT /= len(meta_dtrajs)
    fff = open('%s_MSTS.dat'%mutant_type, 'w')
    for i in range(PPT.shape[0]):
        fff.write('%8.4f '%t_span[i])
        for j in range(n_states):
            fff.write('%6.4f '%PPT[i,j])
        fff.write('\n')
    fff.close()
    MSTS_list.append(PPT)
    
    # bootstrap
    if if_boot:
        boot_stat = bootstrap(boot_fun, meta_dtrajs, n_boot)
    else:
        all_data = np.load('MSTS_data.npz', allow_pickle=True)
        boot_stat = all_data['boot_stat_list'][i_ax]

np.savez('MSTS_data.npz', 
         t_span = t_span,
         MSTS_list = MSTS_list,
         boot_stat_list = boot_stat_list)
