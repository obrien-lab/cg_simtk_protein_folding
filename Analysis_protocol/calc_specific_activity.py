#!/usr/bin/env python3
import sys, math, getopt, os, time, traceback, multiprocessing
import numpy as np
from scipy.integrate import solve_ivp

########## Functions ############
def argmedian(data):
    data = np.array(data)
    if len(data)%2 == 0:
        idx = [int(len(data)/2)-1, int(len(data)/2)]
    else:
        idx = [int(len(data)/2)]
    return np.argsort(data)[idx]

def boot_fun(data):
    global RT
    return(np.median(np.exp(-data/RT)))

def bootstrap(boot_fun, data, n_time):
    idx_list = np.arange(len(data))
    if len(data.shape) == 1:
        boot_stat = np.zeros(n_time)
    elif len(data.shape) == 2:
        boot_stat = np.zeros((data.shape[0] ,n_time))
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
        
def get_dG(pmf_file):
    f = open(pmf_file)
    lines = f.readlines()
    f.close()
    G_list = [float(l.strip().split()[1]) for l in lines if not 'inf' in l]
    G_list = np.array(G_list)
    G_max = np.max(G_list)
    idx = np.argwhere(G_list==G_max)
    idx = idx[0][0]
    G_min = np.min(G_list[:idx+1])
    G_max_err = float(lines[idx].strip().split()[-1])
    idx = np.argwhere(G_list==G_min)
    idx = idx[0][0]
    G_min_err = float(lines[idx].strip().split()[-1])
    G_err = (G_max_err**2 + G_min_err**2)**0.5
    return G_max-G_min, G_err


def permutation_test(perm_stat_fun, data_1, data_2, num_perm, perm_fun):
    global nproc
    combined_data = np.array(list(data_1) + list(data_2))
    t0 = perm_stat_fun(data_1) - perm_stat_fun(data_2)
    
    pool = multiprocessing.Pool(nproc)
    pool_list = []
    start_time = time.time()
    print('start permutation test')
    for i in range(num_perm):
        perm_idx_list = np.random.permutation(np.arange(len(combined_data)))
        pool_list.append(pool.apply_async(perm_fun, (perm_stat_fun, perm_idx_list, combined_data, len(data_1))))
    pool.close()
    pool.join()
    t_dist = [p.get() for p in pool_list]
    p = 0
    for t in t_dist:
        if t >= t0:
            p += 1
    p = (p+1)/(num_perm+1)
    used_time = time.time() - start_time
    print('%.2fs'%used_time)
    return p
    
def perm_fun(perm_stat_fun, perm_idx_list, combined_data, length_1):
    d_1 = perm_stat_fun(combined_data[perm_idx_list[:length_1]])
    d_2 = perm_stat_fun(combined_data[perm_idx_list[length_1:]])
    #print('OK')
    return (d_1-d_2)

def perm_stat_fun(data):
    global dt, end_t, k_median, sol_weight
    C_matrix = data[0][1]
    d0 = [data[0][0]]
    for i in range(1, len(data)):
        C_matrix += data[i][1]
        d0.append(data[i][0])
    rate_matrix = estimate_rate_matrix_2(C_matrix, dt)
    n_states = rate_matrix.shape[0]
    t_span = [0, end_t]
    t_eval = np.linspace(0, end_t, 1000)
    P0 = np.zeros(n_states)
    for md in d0:
        P0[d0] += 1
    P0 = P0 / np.sum(P0)
    sol = Master_equation(t_span, rate_matrix, P0, t_eval)
    w = sol.y[:,-1] * sol_weight
    w = w / np.sum(w)
    P = np.sum(w*k_median)
    return P
    
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

#################################
end_t = 60 # in seconds
dt = 0.015/1000
nsave = 5000
alpha = 4331293.0
nproc = 20
mutant_type_list = ['fast', 'slow']
num_rep = 5
n_perm = int(1e6)
perm_type = 0 #0: simple permutation; 1: permutate trajs
T = 310
n_boot = 10000
QMMM_sample_dir = ''
msm_data_file = ''
mets_data_file = ''
sol_weight = []

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
        if line.startswith('nproc'):
            words = line.split('=')
            nproc = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('num_rep'):
            words = line.split('=')
            num_rep = int(words[1].strip())
            continue
        if line.startswith('n_perm'):
            words = line.split('=')
            n_perm = int(words[1].strip())
            continue
        if line.startswith('perm_type'):
            words = line.split('=')
            perm_type = int(words[1].strip())
            continue
        if line.startswith('T'):
            words = line.split('=')
            T = float(words[1].strip())
            continue
        if line.startswith('n_boot'):
            words = line.split('=')
            n_boot = int(words[1].strip())
            continue
        if line.startswith('QMMM_sample_dir'):
            words = line.split('=')
            QMMM_sample_dir = words[1].strip()
            continue
        if line.startswith('msm_data_file'):
            words = line.split('=')
            msm_data_file = words[1].strip()
            continue
        if line.startswith('mets_data_file'):
            words = line.split('=')
            mets_data_file = words[1].strip()
            continue
        if line.startswith('sol_weight'):
            words = line.split('=')
            words = words[1].strip().split()
            sol_weight = [float(w) for w in words]
            continue
finally:
     file_object.close()    

dt = dt*nsave*alpha/1e9 # in seconds
R = 8.3145*1e-3/4.184 # kcal/mol/K
RT = R*T
plot_num = 1000

mets_data = np.load(mets_data_file, allow_pickle=True)
n_states = mets_data['METS_list'][0].y.shape[0]

if len(sol_weight) == 0:
    sol_weight = np.ones(n_states)
elif len(sol_weight) != n_states:
    print('Error: the length of sol_weight (%d) does not equal to the number of metastable states (%d).'%(len(sol_weight), n_states))
    sys.exit()
    
if perm_type == 0:
    print('Use a simple permutation test for the last SA values.\nThis will only involve the final distribution of SA in the hypothesis test and error estimation.')
else:
    print('Use a permutation test of shuffling the transition count matrices.\nThis will involve the transition kinetics in the hypothesis test and error estimation.')

#######################################
work_dir = os.popen('pwd').readline().strip()
data = {}
data_err = {}
os.chdir(QMMM_sample_dir)
state_list = os.popen('ls -d */').readlines()
state_list = [int(s.strip()[:-1]) for s in state_list]
state_list.sort()

for state_id in state_list:
    rep_list = os.popen('ls %d/'%state_id).readlines()
    data[state_id] = np.Inf*np.ones(num_rep)
    data_err[state_id] = np.Inf*np.ones(num_rep)
    for rep_dir in rep_list:
        rep_dir = rep_dir.strip()
        if not os.path.exists('%d/%s/dE_dG.dat'%(state_id, rep_dir)):
            continue
        pmf_file = '%d/%s/US/pmf.dat'%(state_id, rep_dir)
        dG, dG_err = get_dG(pmf_file)
        data[state_id][int(rep_dir)-1] = dG
        data_err[state_id][int(rep_dir)-1] = dG_err
        
os.chdir(work_dir)

for key in data.keys():
    data[key] = data[key][np.where(np.invert(np.isinf(data[key])))]
    data_err[key] = data_err[key][np.where(np.invert(np.isinf(data_err[key])))]

k_median = np.zeros(n_states)
f_median = open('dG_median.dat', 'w')
for state_id in state_list:
    data_0 = data[state_id]
    data_0_err = data_err[state_id]
    if len(data_0) != 0:
        dG_median = np.median(data_0)
        dG_median_err = data_0_err[argmedian(data_0)]
        if len(dG_median_err) == 1:
            dG_median_err = dG_median_err[0]
        else:
            dG_median_err = np.sum(dG_median_err**2)**0.5/2
        f_median.write('%.4f +/- %.4f '%(dG_median, dG_median_err))
        k_median[state_id-1] = np.exp(-dG_median/RT)
f_median.write('\n')
f_median.close()

perm_data_list = []
ts_list = []
sa_list = []
sa_boot_list = []
for idx, mutant_type in enumerate(mutant_type_list):
    w = mets_data['METS_list'][idx].y.T
    ts = mets_data['METS_list'][idx].t
    w = w * sol_weight
    w = w / (np.sum(w, axis=1).reshape(len(w),1))
    ts_list.append(ts)
    
    sa_list.append(np.dot(w, k_median))
    perm_data = []
    n_traj = len(mets_data['C_matrix_list'][idx])
    n_w = w[-1]*n_traj
    n_w = [round(nw) for nw in n_w]
    n_w[0] = n_traj - sum(n_w[1:])
    for iidx, nw in enumerate(n_w):
        for i in range(int(nw)):
            perm_data.append(k_median[iidx])
    perm_data_list.append(perm_data)
    
    METS_boot_data = mets_data['METS_boot_list'][idx]
    sa_boot = np.zeros((w.shape[0], len(METS_boot_data)))
    for i in range(len(METS_boot_data)):
        w = METS_boot_data[i].y.T * sol_weight
        w = w / (np.sum(w, axis=1).reshape(len(w),1))
        sa_boot[:,i] = np.dot(w, k_median)
    sa_boot_list.append(sa_boot)    

# normalize using max of the end values of fast and slow
sa_max = max([sa_list[0][-1], sa_list[1][-1]])
if sa_list[0][-1] > sa_list[1][-1]:
    idx_0 = 0
    idx_1 = 1
else:
    idx_0 = 1
    idx_1 = 0
norm_sa_list = [sa/sa_max for sa in sa_list]

# normalized by the last boot value of the mutant with larger mean SA
n_boot = sa_boot_list[0].shape[1]
norm_sa_boot_list = [sa/sa_boot_list[idx_0][-1,:] for sa in sa_boot_list]

# permutation test for significant difference between un-normalized fast SA and slow SA at the end of time
if perm_type == 0:
    p = permutation_test(np.mean, perm_data_list[idx_0], perm_data_list[idx_1], n_perm, perm_fun)
else:
    C_list = mets_data['C_matrix_list']
    n_traj = len(C_list[0])
    msm_data = np.load(msm_data_file, allow_pickle=True)
    d0_list = []
    d0_list.append([md[0] for md in msm_data['meta_dtrajs'][:n_traj]])
    d0_list.append([md[0] for md in msm_data['meta_dtrajs'][n_traj:]])
    d_1 = [[d0_list[idx_0][i], C_list[idx_0][i]] for i in range(len(d0_list[idx_0]))]
    d_2 = [[d0_list[idx_1][i], C_list[idx_1][i]] for i in range(len(d0_list[idx_1]))]
    p = permutation_test(perm_stat_fun, d_1, d_2, n_perm, perm_fun)
    
print('Done permutation test')

x = []
y = [np.zeros(plot_num) for mutant_type in mutant_type_list]
y_lb = [np.zeros(plot_num) for mutant_type in mutant_type_list]
y_ub = [np.zeros(plot_num) for mutant_type in mutant_type_list]
for idx, mutant_type in enumerate(mutant_type_list):
    ts_idx = np.linspace(0, len(ts_list[idx])-1, plot_num)
    ts_idx = ts_idx.astype(int)
    
    fo = open('%s_SA.dat'%mutant_type, 'w')
    iidx = 0
    for i, t in enumerate(ts_list[idx]):
        sa_mean = sa_list[idx][i]
        sa_boot = sa_boot_list[idx][i,:]
        sa_lb = np.percentile(sa_boot, 2.5)
        sa_ub = np.percentile(sa_boot, 97.5)
        norm_sa_mean = norm_sa_list[idx][i]
        norm_sa_boot = norm_sa_boot_list[idx][i,:]
        norm_sa_lb = np.percentile(norm_sa_boot, 2.5)
        norm_sa_ub = np.percentile(norm_sa_boot, 97.5)
        fo.write('%8.4f %.4e [%.4e, %.4e] %.4f [%.4f, %.4f]\n'%(t, sa_mean, 
                 sa_lb, sa_ub, norm_sa_mean, norm_sa_lb, norm_sa_ub))
        if i in ts_idx:
            y[idx][iidx] = norm_sa_mean
            y_lb[idx][iidx] = norm_sa_lb
            y_ub[idx][iidx] = norm_sa_ub
            iidx += 1
    fo.close()

    x.append(ts_list[idx][ts_idx])
    print('Done writing %s SA.'%mutant_type)

# estimate error bar
y_end = np.array([sa_list[idx][-1] for idx in range(len(mutant_type_list))])
if y_end[0] > y_end[1]:
    idx_0 = 0
    idx_1 = 1
else:
    idx_0 = 1
    idx_1 = 0
y_end = y_end / sa_list[idx_0][-1]

if perm_type == 1:
    y_boot_end = np.array([sa_boot_list[idx][-1,:] for idx in range(len(mutant_type_list))])
    y_boot_end_1 = y_boot_end / y_boot_end[idx_0,:] # normalized by the boot value of the mutant with larger mean SA
    y_boot_end_2 = y_boot_end / y_boot_end[idx_0,:]
else:
    y_boot_end = []
    for idx in range(len(mutant_type_list)):
        y_boot_end.append(bootstrap(np.mean, np.array(perm_data_list[idx]), n_perm))
    y_boot_end = np.array(y_boot_end)
    y_boot_end_2 = y_boot_end / y_boot_end[idx_0,:]

y_end_lb = np.percentile(y_boot_end_2, 2.5, axis=1)
y_end_ub = np.percentile(y_boot_end_2, 97.5, axis=1)

for idx, mutant_type in enumerate(mutant_type_list):
    print("%s SA = %.4e; 95%%CI = [%.4e, %.4e]"%(mutant_type, sa_list[idx][-1], np.percentile(y_boot_end, 2.5, axis=1)[idx], 
                                                 np.percentile(y_boot_end, 97.5, axis=1)[idx]))
print('%s normalized SA = %.4f; 95%%CI = [%.4f, %.4f]; p = %.4e'%(mutant_type_list[idx_1], y_end[idx_1],  y_end_lb[idx_1], y_end_ub[idx_1], p))
