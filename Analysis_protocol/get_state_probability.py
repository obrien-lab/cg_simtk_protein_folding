#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
from scipy.integrate import solve_ivp
import msmtools

################################# Arguments ###################################
# Default values
end_t = 60 # in seconds
dt = 0.015/1000
nsave = 5000
alpha = 4331293.0
n_window = 200
n_traj = 100
mutant_type_list = ['fast', 'slow']
lag_t = 1
start_idx = 1
end_idx = 10
msm_data_file = './msm_data.npz'
prefix_dir = ''
n_boot = 10000
num_proc = 20
if_extend_dtraj = True
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
        if line.startswith('n_window'):
            words = line.split('=')
            n_window = int(words[1].strip())
            continue
        if line.startswith('n_traj'):
            words = line.split('=')
            n_traj = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('lag_t'):
            words = line.split('=')
            lag_t = int(words[1].strip())
            continue
        if line.startswith('start_idx'):
            words = line.split('=')
            start_idx = int(words[1].strip())
            continue
        if line.startswith('end_idx'):
            words = line.split('=')
            end_idx = int(words[1].strip())
            continue
        if line.startswith('msm_data_file'):
            words = line.split('=')
            msm_data_file = words[1].strip()
            continue
        if line.startswith('prefix_dir'):
            words = line.split('=')
            prefix_dir = words[1].strip()
            continue
        if line.startswith('n_boot'):
            words = line.split('=')
            n_boot = int(words[1].strip())
            continue
        if line.startswith('num_proc'):
            words = line.split('=')
            num_proc = int(words[1].strip())
            continue
        if line.startswith('if_extend_dtraj'):
            words = line.split('=')
            if_extend_dtraj = int(words[1].strip())
            if if_extend_dtraj == 1:
                if_extend_dtraj = True
            elif if_extend_dtraj == 0:
                if_extend_dtraj = False
            else:
                print('Error: if_extend_dtraj can only be either 0 or 1.')
                sys.exit()
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

# Building Master Equations
def dP(t, P, M):
    return np.dot(M,P)

def Master_equation(t_span, K, P0, t_eval):
    sol = solve_ivp(dP, t_span, P0, t_eval=t_eval, args=(K,))
    return sol
    
def estimate_rate_matrix_2(C_matrix, dt):
    n_state = len(C_matrix)
    ksum_matrix = np.zeros((n_state, n_state))
    P_matrix = np.zeros((n_state, n_state))
    for i in range(n_state):
        if C_matrix[i,i] != 0:
            ksum_matrix[i,i] = -math.log(C_matrix[i,i]/np.sum(C_matrix[i,:]))/dt
        else:
            ksum_matrix[i,i] = 2/dt # assume the mean dwell time of this state is dt/2
        for j in range(n_state):
            if i != j:
                P_matrix[j,i] = C_matrix[i,j]
        if np.sum(P_matrix[:,i]) != 0:
            P_matrix[:,i] /= np.sum(P_matrix[:,i])
    
    rate_matrix = np.dot(P_matrix, ksum_matrix)
    for i in range(n_state):
        rate_matrix[i,i] = -np.sum(rate_matrix[:,i])
        
    return rate_matrix

def exp_fun(x, k):
    return np.exp(-k*x)

def boot_fun(data, d0, dt, end_t, i):
    rate_matrix = estimate_rate_matrix_2(data, dt)
    n_states = rate_matrix.shape[0]
    t_span = [0, end_t]
    t_eval = np.linspace(0, end_t, 1000)
    P0 = np.zeros(n_states)
    for md in d0:
        P0[d0] += 1
    P0 = P0 / np.sum(P0)
    sol = Master_equation(t_span, rate_matrix, P0, t_eval)
    #print('Bootstrapping done for %d'%(i+1))
    return sol

def bootstrap(boot_fun, data, d0, n_boot, dt, end_t):
    global num_proc
    
    n_states = len(data[0])
    
    pool = multiprocessing.Pool(num_proc)
    pool_list = []
    
    idx_list = np.arange(len(data))
    
    start_time = time.time()
    print('start bootsrapping')
    for i in range(n_boot):
        sample_idx_list = np.random.choice(idx_list, len(idx_list))
        new_data = np.zeros((n_states,n_states))
        for idx in sample_idx_list:
            new_data += data[idx]
        pool_list.append(pool.apply_async(boot_fun, (new_data, d0[sample_idx_list], dt, end_t, i, )))
    pool.close()
    pool.join()
    boot_stat = [p.get() for p in pool_list]
    used_time = time.time() - start_time
    print('%.2fs'%used_time)
    return boot_stat

################################## MAIN #######################################
npzfile = np.load(msm_data_file, allow_pickle=True)

C_matrix_list = []
METS_list = []
METS_boot_list = []
si = 0
for i_ax, mutant_type in enumerate(mutant_type_list):
    po_dir = prefix_dir+'/post_translation/'+mutant_type+'/1-100/'
    ei = si+n_traj*(end_idx-start_idx+1)
    meta_dtrajs = npzfile['meta_dtrajs'][np.arange(si, ei)]
    si = ei
    n_states = 0
    for md in meta_dtrajs:
        if n_states < np.max(md):
            n_states = np.max(md)
    n_states += 1
    
    max_T_len = 0
    for i in range(len(meta_dtrajs)):
        if meta_dtrajs[i].shape[0] > max_T_len:
            max_T_len = meta_dtrajs[i].shape[0]
    
    if if_extend_dtraj:
        C_list = []
        for i in range(len(meta_dtrajs)):
            (N, be) = np.histogram(meta_dtrajs[i][-n_window:], bins=np.arange(-0.5, n_states, 1))
            meta_dtraj_last = np.argwhere(N == np.max(N))[0][0]
            mde = []
            for j in range(max_T_len):
                if j >= len(meta_dtrajs[i]): 
                    state_0 = meta_dtraj_last
                else:
                    state_0 = meta_dtrajs[i][j]
                mde.append(state_0)
            
            C = msmtools.estimation.count_matrix(mde, lag_t)
            C = C.toarray()
            if len(C) < n_states:
                C1 = np.zeros((n_states,n_states))
                C1[:len(C),:len(C)] = C
                C_list.append(C1)
            else:
                C_list.append(C)
            #print('Done for traj %d'%(i+1))
    else:
        all_data = np.load('METS_data_lag%d.npz'%lag_t, allow_pickle=True)
        C_list = all_data['C_matrix_list'][i_ax]
    C_matrix_list.append(C_list)
    
    if not if_boot:
        all_data = np.load('METS_data_lag%d.npz'%lag_t, allow_pickle=True)
        sol = all_data['METS_list'][i_ax]
        t_span = sol.t
        PP = sol.y
        boot_stat = list(all_data['METS_boot_list'][i_ax])
        METS_list.append(sol)
        METS_boot_list.append(boot_stat)
    else:
        C_matrix = np.zeros((n_states, n_states))
        for C in C_list:
            C_matrix += C
        rate_matrix = estimate_rate_matrix_2(C_matrix, dt*lag_t)
        
        t_span = [0, end_t]
        t_eval = np.linspace(0, end_t, 1000)
        P0 = np.zeros(n_states)
        for md in meta_dtrajs:
            P0[md[0]] += 1
        P0 = P0 / np.sum(P0)
        sol = Master_equation(t_span, rate_matrix, P0, t_eval)
        t_span = sol.t
        PP = sol.y
        
        METS_list.append(sol)
        
        # bootstrap
        d0 = np.array([mde[0] for mde in meta_dtrajs])
        boot_stat = bootstrap(boot_fun, C_list, d0, n_boot, dt*lag_t, end_t)
        METS_boot_list.append(boot_stat)
        
    fff = open('%s_METS_lag%d.dat'%(mutant_type, lag_t), 'w')
    for i in range(t_span.shape[0]):
        fff.write('%8.4f '%(t_span[i]))
        for j in range(n_states):
            fff.write('%6.4f '%abs(PP[j,i]))
        fff.write('\n')
    fff.close()

np.savez('METS_data_lag%d.npz'%lag_t, 
         C_matrix_list = C_matrix_list,
         METS_list = METS_list,
         METS_boot_list = METS_boot_list)
