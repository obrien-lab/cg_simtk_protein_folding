#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import pyemma as pem
import parmed as pmd
import mdtraj as mdt
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
n_cluster = 400
stride=10
n_large_states = 10
n_small_states = 2
lag_t = 1
start_idx = 1
end_idx = 10
sample_size = 5
native_AA_pdb = ''
prefix_dir = ''
visualiz_threshold = 0.02
if_cluster = True
if_visualize = True
if_sample = True

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
        if line.startswith('n_cluster'):
            words = line.split('=')
            n_cluster = int(words[1].strip())
            continue
        if line.startswith('stride'):
            words = line.split('=')
            stride = int(words[1].strip())
            continue
        if line.startswith('n_large_states'):
            words = line.split('=')
            n_large_states = int(words[1].strip())
            continue
        if line.startswith('n_small_states'):
            words = line.split('=')
            n_small_states = int(words[1].strip())
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
        if line.startswith('sample_size'):
            words = line.split('=')
            sample_size = int(words[1].strip())
            continue
        if line.startswith('native_AA_pdb'):
            words = line.split('=')
            native_AA_pdb = words[1].strip()
            continue
        if line.startswith('prefix_dir'):
            words = line.split('=')
            prefix_dir = words[1].strip()
            continue
        if line.startswith('visualiz_threshold'):
            words = line.split('=')
            visualiz_threshold = float(words[1].strip())
            continue
        if line.startswith('if_cluster'):
            words = line.split('=')
            if_cluster = int(words[1].strip())
            if if_cluster == 1:
                if_cluster = True
            elif if_cluster == 0:
                if_cluster = False
            else:
                print('Error: if_cluster can only be either 0 or 1.')
                sys.exit()
            continue
        if line.startswith('if_visualize'):
            words = line.split('=')
            if_visualize = int(words[1].strip())
            if if_visualize == 1:
                if_visualize = True
            elif if_visualize == 0:
                if_visualize = False
            else:
                print('Error: if_visualize can only be either 0 or 1.')
                sys.exit()
            continue
        if line.startswith('if_sample'):
            words = line.split('=')
            if_sample = int(words[1].strip())
            if if_sample == 1:
                if_sample = True
            elif if_sample == 0:
                if_sample = False
            else:
                print('Error: if_sample can only be either 0 or 1.')
                sys.exit()
            continue
finally:
     file_object.close()

dt = dt*nsave*alpha/1e9 # in seconds
################################# Functions ###################################
def standardize(data):
    data_con = data[0]
    for i in range(1, len(data)):
        data_con = np.vstack((data_con, data[i]))
    data_mean = np.mean(data_con, axis=0)
    data_std = np.std(data_con, axis=0)
    result = [(d - data_mean) / data_std for d in data]
    return [result, data_mean, data_std]

def unstandardize(data, data_mean, data_std):
    result = data * data_std + data_mean
    return result

# Building Master Equations
def dP(t, P, M):
    return np.dot(M,P)

def Master_equation(t_span, K, P0, t_eval):
    sol = solve_ivp(dP, t_span, P0, t_eval=t_eval, args=(K,))
    return sol
    
def estimate_rate_matrix(mutant_type, dtrajs, n_state, dt):    
    ksum_matrix = np.zeros((n_state, n_state))
    P_matrix = np.zeros((n_state, n_state))
    C_matrix = msmtools.estimation.count_matrix(dtrajs, 1)
    C_matrix = C_matrix.toarray()
    if len(C_matrix) != n_state:
        for i in range(len(C_matrix), n_state):
            C_matrix.append(np.zeros((1,C_matrix.shape[1])), axis=0)
            C_matrix.append(np.zeros((C_matrix.shape[0],1)), axis=1)
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
    # enforce native state to be a sink
    # for i in range(n_state):
    #     rate_matrix[i,-1] = 0
    for i in range(n_state):
        rate_matrix[i,i] = -np.sum(rate_matrix[:,i])
        
    return rate_matrix

def exp_fun(x, k):
    return np.exp(-k*x)

def calc_G_list(coor, sel, cutoff, terminal_cutoff):
    n_atom = coor.shape[0]
    # Generate contact matrix
    R = np.zeros((n_atom-1, 3))
    dR = np.zeros((n_atom-1, 3))
    for i in range(n_atom-1):
        R[i] = (coor[i, :] + coor[i+1, :])/2
        dR[i] = coor[i+1, :] - coor[i, :]
    M = np.zeros((n_atom-1, n_atom-1))
    for i in range(n_atom-2):
        for j in range(i+1, n_atom-1):
            v1 = (R[i] - R[j]) / np.sum((R[i] - R[j])**2)**(3/2)
            v2 = np.cross(dR[i], dR[j])
            M[i,j] = np.dot(v1, v2)
            M[j,i] = M[i,j]

    # Calculate G
    G_list = []
    for ii, r1_range in enumerate(sel):
        r1_i = r1_range[0]
        r1_j = r1_range[1]
        coor_1 = coor[r1_i, :]
        coor_2 = coor[r1_j, :]
        dist = np.sum((coor_1-coor_2)**2)**0.5
        if dist < cutoff:
            G0 = [0,0]
            for idx, r2_range in enumerate([[terminal_cutoff,r1_i-4], [r1_j+4, n_atom-terminal_cutoff-1]]):
                r2_i = r2_range[0]
                r2_j = r2_range[1]
                for r1 in range(r1_i, r1_j):
                    for r2 in range(r2_i, r2_j):
                        G0[idx] += M[r1, r2]
            G_list.append([G0[0]/4/3.14, G0[1]/4/3.14])
        else:
            G_list.append([np.nan, np.nan])
    return (M, G_list)

def gen_state_visualizion(state_id, psf, native_cor, state_cor, native_AA_pdb, if_entangled):
    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP']
    
    print('Generate visualization of state %d'%(state_id))
    os.system('mkdir state_struct')
    os.chdir('state_struct')
    
    cutoff = 0.8 # in nm
    min_interval = 4
    terminal_cutoff = 5
    
    ## Get native contacts ##
    struct = pmd.load_file(psf)
    coor = pmd.load_file(native_cor)
    struct.coordinates = coor.coordinates
    state_cor[0].save('tmp.pdb', force_overwrite=True)

    native_contact = []
    for i in range(0, len(struct.atoms)-min_interval):
        coor_1 = np.array([struct[i].xx, struct[i].xy, struct[i].xz])
        for j in range(i+min_interval, len(struct.atoms)):
            coor_2 = np.array([struct[j].xx, struct[j].xy, struct[j].xz])
            dist = np.sum((coor_1-coor_2)**2)**0.5
            if (dist < 10*cutoff):
                native_contact.append([i, j])
    
    _, G_native_list = calc_G_list(struct.coordinates/10, native_contact, cutoff, terminal_cutoff)
    
    M, G_list = calc_G_list(state_cor.xyz[0,:,:], native_contact, cutoff, terminal_cutoff)
    
    g_max = 0
    i_max = 0
    for i in range(len(native_contact)):
        if abs(G_list[i][0] + G_list[i][1]) - abs(G_native_list[i][0] + G_native_list[i][1]) > g_max and not np.isnan(G_list[i][1]):
            i_max = i
            g_max = abs(G_list[i][0] + G_list[i][1]) - abs(G_native_list[i][0] + G_native_list[i][1]) 
            
    idx_max = [native_contact[i_max][0]+1, native_contact[i_max][1]+1]
            
    min_dist = idx_max[1] - idx_max[0]
    idx_max_0 = [idx_max[0], idx_max[1]]
    for i, nc in enumerate(native_contact):
        if nc[0] > idx_max_0[0]-1 and nc[1] < idx_max_0[1]-1 and nc[1] - nc[0] < min_dist and nc[1] - nc[0] > 15 \
           and not np.isnan(G_list[i][1]) and abs(round(G_list[i][1])+round(G_list[i][0]))>abs(round(G_native_list[i][1])+round(G_native_list[i][0])) \
           and abs(round(G_list[i][1])+round(G_list[i][0])) == abs(round(G_list[i_max][1])+round(G_list[i_max][0])):
            idx_max = [nc[0]+1, nc[1]+1]
            min_dist = nc[1] - nc[0]
            i_max = i
    
    if abs(G_list[i_max][0]) > abs(G_list[i_max][1]):
        idx_thread = [terminal_cutoff+1, native_contact[i_max][0]-4+1]
        G_list_max = G_list[i_max][0]
    else:
        idx_thread = [native_contact[i_max][1]+4+1, len(struct.atoms)-terminal_cutoff]
        G_list_max = G_list[i_max][1]
    g_max = 0
    len_thread = 20
    idx_thread_max = [idx_thread[0], idx_thread[1]]    
    for i in range(idx_thread[0]-1, idx_thread[1]-len_thread+1):
        j = i+len_thread
        g = 0
        for r1 in range(native_contact[i_max][0], native_contact[i_max][1]):
            for r2 in range(i, j):
                g += M[r1, r2]
        if abs(g) > g_max:
            g_max = abs(g)
            idx_thread_max = [i+1, j]
    
    # backmap
    os.system('backmap.py -i '+native_AA_pdb+' -c tmp.pdb')
    os.system('mv tmp_rebuilt.pdb state_%d.pdb'%state_id)
    os.system('rm -f tmp.pdb')
    os.system('rm -rf ./rebuild_tmp/')
    
    pdb_struct = pmd.load_file(native_AA_pdb)
    idx_offset = 0
    for res in pdb_struct.residues:
        if res.name in AA_name_list:
            idx_offset = res.number - 1
            break
    
    f = open('vmd_%d.tcl'%state_id, 'w')
    if if_entangled:
        f.write('''display rendermode GLSL
axes location off

color Display {Background} white

mol new '''+native_AA_pdb+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 6
mol selection {resid '''+('%d'%(1+idx_offset))+''' to '''+('%d'%(len(struct.atoms)+idx_offset))+'''}
mol material AOChalky
mol addrep top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 4
mol selection {resid '''+('%d'%(idx_max[0]+idx_offset))+''' to '''+('%d'%(idx_max[1]+idx_offset))+'''}
mol material Opaque
mol addrep top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 12
mol selection {resid '''+('%d'%(idx_thread_max[0]+idx_offset))+''' to '''+('%d'%(idx_thread_max[1]+idx_offset))+'''}
mol material Opaque
mol addrep top
mol representation VDW 1.000000 12.000000
mol color Name
mol selection {not resid '''+('%d'%(1+idx_offset))+''' to '''+('%d'%(len(struct.atoms)+idx_offset))+''' and not water}
mol material Opaque
mol addrep top

mol new ./'''+('state_%d.pdb'%state_id)+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 6
mol selection {all}
mol material AOChalky
mol addrep top
mol representation NewCartoon 0.350000 10.000000 4.100000 0
mol color ColorID 1
mol selection {resid '''+('%d'%(idx_max[0]))+''' to '''+('%d'%(idx_max[1]))+'''}
mol material Opaque
mol addrep top
mol representation NewCartoon 0.350000 10.000000 4.100000 0
mol color ColorID 0
mol selection {resid '''+('%d'%(idx_thread_max[0]))+''' to '''+('%d'%(idx_thread_max[1]))+'''}
mol material Opaque
mol addrep top
mol representation VDW 1.000000 12.000000
mol color ColorID 3
mol selection {resid '''+('%d'%(idx_max[0]))+''' '''+('%d'%(idx_max[1]))+''' and name CA}
mol material Opaque
mol addrep top

set sel [atomselect top "resid '''+('%d'%(idx_max[0]))+''' '''+('%d'%(idx_max[1]))+''' and name CA"]
set idx [$sel get index]
topo addbond [lindex $idx 0] [lindex $idx 1]
mol representation Bonds 0.300000 12.000000
mol color ColorID 3
mol selection {resid '''+('%d'%(idx_max[0]))+''' '''+('%d'%(idx_max[1]))+''' and name CA}
mol material Opaque
mol addrep top

set sel1 [atomselect 0 "resid '''+('%d'%(1+idx_offset))+''' to '''+('%d'%(len(struct.atoms)+idx_offset))+''' and not (resid '''+('%d'%(idx_max[0]+idx_offset))+''' to '''+('%d'%(idx_max[1]+idx_offset))+''' '''+('%d'%(idx_thread_max[0]+idx_offset))+''' to '''+('%d'%(idx_thread_max[1]+idx_offset))+''') and name CA"]
set sel2 [atomselect 1 "resid 1 to '''+('%d'%len(struct.atoms))+''' and not (resid '''+('%d'%(idx_max[0]))+''' to '''+('%d'%(idx_max[1]))+''' '''+('%d'%(idx_thread_max[0]))+''' to '''+('%d'%(idx_thread_max[1]))+''') and name CA"]
set trans_mat [measure fit $sel1 $sel2]
set move_sel [atomselect 0 "all"]
$move_sel move $trans_mat
''')
    else:
        f.write('''display rendermode GLSL
axes location off

color Display {Background} white

mol new ./'''+('state_%d.pdb'%state_id)+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 10
mol selection {all}
mol material AOChalky
mol addrep top
''')
    f.close()
    os.chdir('../')

def get_co_po_dir(prefix_dir, mutant_type):
    global n_traj
    co_dir = prefix_dir+'/continuous_synthesis/'+mutant_type+('/1-%d/'%n_traj)
    po_dir = prefix_dir+'/post_translation/'+mutant_type+('/1-%d/'%n_traj)
    psf_file = os.popen('ls %s/setup/*_ca.psf'%po_dir).readlines()[0].strip()
    cor_file = os.popen('ls %s/setup/*_ca.cor'%po_dir).readlines()[0].strip()
    return (co_dir, po_dir, psf_file, cor_file)

################################## MAIN #######################################
if not if_cluster:
    npzfile = np.load('msm_data.npz', allow_pickle=True)

prefix_file = '_act'

G_list_0_list = []
cor_list = []
trajid_list = []
trajid2mtype = []
mtype2trajid = []
rate_matrix_list = []
xlim_list = []
ylim_list = []
fig_list = []
fig_tpt_list = []

# combine trajs and do clustering and PCCA++
for i_ax, mutant_type in enumerate(mutant_type_list):
    co_dir, po_dir, psf_file, cor_file = get_co_po_dir(prefix_dir, mutant_type)
    
    max_length = len(pmd.load_file(psf_file).atoms)
    
    # Get number of native contact
    cutoff = 0.8 # in nm
    min_interval = 4
    struct = pmd.load_file(psf_file)
    coor = pmd.load_file(cor_file)
    struct.coordinates = coor.coordinates
    native_contact = []
    for i in range(0, len(struct.atoms)-min_interval):
        coor_1 = np.array([struct[i].xx, struct[i].xy, struct[i].xz])
        for j in range(i+min_interval, len(struct.atoms)):
            coor_2 = np.array([struct[j].xx, struct[j].xy, struct[j].xz])
            dist = np.sum((coor_1-coor_2)**2)**0.5
            if (dist < 10*cutoff):
                native_contact.append([i, j])
    num_nc = len(native_contact)
    print('# of native contacts: %d'%num_nc)

    # Get Qbb trajectory
    qbb_list = list(np.load(mutant_type+'_QBB'+prefix_file+'.npy', allow_pickle=True))
    qbb_list = [q.reshape(q.shape[0],1) for q in qbb_list]
    
    # Get Entanglement trajectory
    G_list_0 = list(np.load(mutant_type+'_Entanglement.npy', allow_pickle=True))
    G_list_0_list += G_list_0
    G_list = [np.sum(g[:,1:6], axis=1)/num_nc for g in G_list_0]
    #G_list = [np.sum(g[:,1:3], axis=1)/num_nc for g in G_list_0]
    #G_list = [g[:,-1] for g in G_list_0]
    G_list = [g.reshape(g.shape[0],1) for g in G_list]
    
    for i in range(n_traj):
        for j in range(start_idx-1,end_idx):
            out_file = po_dir+str(i+1)+'/traj_'+str(i+1)+'_'+str(j+1)+'.out'
            ff = open(out_file)
            lines = ff.readlines()
            ff.close()
            if lines[-1].startswith('Done'):
                end_frame = int(int(lines[-2].strip().split()[1])/nsave)
            else:
                end_frame = int(int(lines[-1].strip().split()[1])/nsave)
            if len(qbb_list[(end_idx-start_idx+1)*i+j]) > end_frame:
                qbb_list[(end_idx-start_idx+1)*i+j] = qbb_list[(end_idx-start_idx+1)*i+j][:end_frame]
                qbb_list[(end_idx-start_idx+1)*i+j].reshape(qbb_list[(end_idx-start_idx+1)*i+j].shape[0],1)
            if len(G_list[(end_idx-start_idx+1)*i+j]) > end_frame:
                G_list[(end_idx-start_idx+1)*i+j] = G_list[(end_idx-start_idx+1)*i+j][:end_frame]
                G_list[(end_idx-start_idx+1)*i+j].reshape(G_list[(end_idx-start_idx+1)*i+j].shape[0],1)
    
    for i in range(len(qbb_list)):
        if len(qbb_list[i]) != len(G_list[i]):
            print('Traj %d_%d has mismatched data length (%d qbb_list vs. %d G_list)'%
                  (int(i/(end_idx-start_idx+1))+1, 
                   i-int(i/(end_idx-start_idx+1))*(end_idx-start_idx+1)+1, 
                   len(qbb_list[i]), len(G_list[i])))
            min_idx = min([len(qbb_list[i]), len(G_list[i])])
            qbb_list[i] = qbb_list[i][:min_idx]
            G_list[i] = G_list[i][:min_idx]
            
    trajid_list += [i for i in range(len(qbb_list))]
    mtype2trajid.append([i+len(cor_list) for i in range(len(qbb_list))])
    trajid2mtype += [i_ax for i in range(len(qbb_list))]
    cor_list += [np.hstack((qbb_list[i], G_list[i])) for i in range(len(qbb_list))]
cor_list = np.array(cor_list)
    
#Clustering
if if_cluster:
    std_cor_list, cor_mean, cor_std = standardize(cor_list)
    cluster = pem.coordinates.cluster_kmeans(std_cor_list, k=n_cluster, max_iter=5000, stride=stride)
    dtrajs = cluster.dtrajs
    center = unstandardize(cluster.clustercenters, cor_mean, cor_std)
else:
    dtrajs = list(npzfile['dtrajs'])
    center = npzfile['center']

# Get connective groups and build MSMs
c_matrix = msmtools.estimation.count_matrix(dtrajs, lag_t).toarray()
sub_groups = msmtools.estimation.connected_sets(c_matrix)
active_groups = []
for sg in sub_groups:
    for ssg in sg:
        tag_found = False
        for dtraj in dtrajs:
            if ssg in dtraj:
                tag_found = True
                break
        if not tag_found:
            break
    if tag_found:
        active_groups.append(sg)
print('Total number of active groups: %d'%(len(active_groups)))

msm_list = []        
for ag in active_groups:
    cm = msmtools.estimation.largest_connected_submatrix(c_matrix, lcc=ag)
    if len(cm) == 1:
        msm = None
    else:
        T = msmtools.estimation.transition_matrix(cm, reversible=True)
        msm = pem.msm.markov_model(T, dt_model=str(dt)+' s')
    msm_list.append(msm)

meta_dist = []
meta_set = []
eigenvalues_list = []
for idx_msm, msm in enumerate(msm_list):
    if idx_msm == 0:
        n_states = n_large_states
    else:
        n_states = n_small_states
    if msm == None:
        eigenvalues_list.append(None)
        dist = np.zeros(n_cluster)
        iidx = active_groups[idx_msm][0]
        dist[iidx] = 1.0
        meta_dist.append(dist)
        meta_set.append(active_groups[idx_msm])
    else:
        eigenvalues_list.append(msm.eigenvalues())
        # coarse-graining 
        while n_states > 1:
            tag_empty = False
            pcca = msm.pcca(n_states)
            for ms in msm.metastable_sets:
                if ms.size == 0:
                    tag_empty = True
                    break
            if not tag_empty:
                break
            else:
                n_states -= 1
                print('Reduced number of states to %d for active group %d'%(n_states, idx_msm+1))
        if n_states == 1:
            # use observation prob distribution for non-active set
            dist = np.zeros(n_cluster)
            for nas in active_groups[idx_msm]:
                for dtraj in dtrajs:
                    dist[nas] += np.count_nonzero(dtraj == nas)
            dist /= np.sum(dist)
            meta_dist.append(dist)
            meta_set.append(active_groups[idx_msm])
        else:
            for i, md in enumerate(msm.metastable_distributions):
                dist = np.zeros(n_cluster)
                s = np.sum(md[msm.metastable_sets[i]])
                set_0 = []
                for idx in msm.metastable_sets[i]:
                    iidx = active_groups[idx_msm][idx]
                    dist[iidx] = md[idx]
                    set_0.append(iidx)
                dist = dist / s
                meta_dist.append(dist)
                meta_set.append(set_0)
meta_dist = np.array(meta_dist)
meta_set = np.array(meta_set)

coarse_state_centers = center[meta_dist.argmax(1)]
cg_center_order_idx = np.argsort(coarse_state_centers[:,0])
micro_to_meta = np.zeros(n_cluster)
meta_set = meta_set[cg_center_order_idx]
meta_dist = meta_dist[cg_center_order_idx, :]
for idx, ms in enumerate(meta_set):
    for mms in ms:
        micro_to_meta[mms] = idx
meta_dtrajs = []
for traj in dtrajs:
    meta_traj = np.zeros(len(traj), dtype=int)
    for i, idx in enumerate(traj):
        meta_traj[i] = micro_to_meta[idx]
    meta_dtrajs.append(meta_traj)
meta_dtrajs = np.array(meta_dtrajs)

n_states = len(meta_set)
print('Total %d metastable states were grouped'%n_states)

if if_sample:
    cluster_indexes = pem.util.discrete_trajectories.index_states(dtrajs)
    if len(cluster_indexes) < n_cluster:
        cluster_indexes = list(cluster_indexes)
        for i in range(len(cluster_indexes), n_cluster):
            cluster_indexes.append(np.array([[]]))
        cluster_indexes = np.array(cluster_indexes)
    samples = pem.util.discrete_trajectories.sample_indexes_by_distribution(cluster_indexes, 
                                                                            meta_dist, 
                                                                            sample_size)
else:
    samples = npzfile['meta_samples']
meta_samples = samples

sampled_traj = None
visualiz_G = []
for i, meta_state in enumerate(samples):
    visualiz_G.append(cor_list[meta_state[0][0]][meta_state[0][1],1])
    for idx in meta_state:
        traj_idx = idx[0]
        frame_idx = idx[1]
        co_dir, po_dir, psf_file, cor_file = get_co_po_dir(prefix_dir, mutant_type_list[trajid2mtype[traj_idx]])
        traj_idx = trajid_list[traj_idx]
        traj_idx_1 = int(traj_idx / (end_idx-start_idx+1))
        traj_idx_2 = int(traj_idx - traj_idx_1 * (end_idx-start_idx+1))
        traj = mdt.load(co_dir+'traj/'+str(traj_idx_1+1)+'/rnc_l'+str(max_length)+'_ejection.dcd', top=co_dir+'traj/'+str(traj_idx_1+1)+'/rnc_l'+str(max_length)+'.psf')
        traj += mdt.load(co_dir+'traj/'+str(traj_idx_1+1)+'/rnc_l'+str(max_length)+'_dissociation.dcd', top=co_dir+'traj/'+str(traj_idx_1+1)+'/rnc_l'+str(max_length)+'.psf')
        sel = traj.topology.select('resid 0 to '+str(max_length-1))
        traj = traj.atom_slice(sel)
        traj += mdt.load(po_dir+str(traj_idx_1+1)+'/traj_'+str(traj_idx_1+1)+'_'+str(traj_idx_2+1)+'.dcd', top=psf_file)
        if sampled_traj is None:
            sampled_traj = traj[frame_idx]
        else:
            sampled_traj += traj[frame_idx]
    print('Get samples of metastable state %d'%(i+1))
sampled_traj = sampled_traj.center_coordinates()
sampled_traj = sampled_traj.superpose(sampled_traj)
sampled_traj.save('sampled_traj.dcd', force_overwrite=True)

if_entangled_list = [False for i in range(n_states)]
# analysis MSM for each mutant
for i_ax, mutant_type in enumerate(mutant_type_list):
    co_dir, po_dir, psf_file, cor_file = get_co_po_dir(prefix_dir, mutant_type)
    
    state_indices = []
    for state_id in range(0, n_states):
        state_indices.append([])
        for i_1, md in enumerate(meta_dtrajs[mtype2trajid[i_ax]]):
            for i_2, mdd in enumerate(md):
                if mdd == state_id:
                    state_indices[-1].append([mtype2trajid[i_ax][i_1], i_2])
    for state_id in range(1, n_states+1):
        if len(state_indices[state_id-1]) == 0:
            continue
        g_list = np.zeros(5)
        G_avg = 0
        for si in state_indices[state_id-1]:
            g_list += G_list_0_list[si[0]][si[1],1:6]
            G_avg += cor_list[si[0]][si[1],1]
        g_list /= len(state_indices[state_id-1])
        g_list /= g_list.sum()
        g_list *= 100
        G_avg /= len(state_indices[state_id-1])
        if g_list[0]+g_list[1] > 50 and G_avg > visualiz_threshold:
            if_entangled_list[state_id-1] = True
    
    max_T_len = 0
    for i in mtype2trajid[i_ax]:
        if meta_dtrajs[i].shape[0] > max_T_len:
            max_T_len = meta_dtrajs[i].shape[0]
    
    PPT = np.zeros((max_T_len, n_states))
    meta_dtrajs_extended = []
    for i in mtype2trajid[i_ax]:
        (N, be) = np.histogram(meta_dtrajs[i][-n_window:], bins=np.arange(-0.5, n_states, 1))
        meta_dtraj_last = np.argwhere(N == np.max(N))[0][0]
        for j in range(max_T_len):
            if j >= len(meta_dtrajs[i]): 
                state_0 = meta_dtraj_last
            else:
                state_0 = meta_dtrajs[i][j]
            PPT[j,state_0] += 1
        mde = []
        for j in range(max_T_len):
            if j >= len(meta_dtrajs[i]): 
                state_0 = meta_dtraj_last
            else:
                state_0 = meta_dtrajs[i][j]
            mde.append(state_0)
        meta_dtrajs_extended.append(mde)
        
    rate_matrix = estimate_rate_matrix(mutant_type, meta_dtrajs_extended, n_states, dt)
    rate_matrix_list.append(rate_matrix)
    
    fff = open('%s_MSTS.dat'%mutant_type, 'w')
    for i in range(PPT.shape[0]):
        PPT[i,:] = PPT[i,:] / np.sum(PPT[i,:])
        fff.write('%8.4f '%((i+1)*dt))
        for j in range(n_states):
            fff.write('%6.4f '%PPT[i,j])
        fff.write('\n')
    fff.close()
    
    # build master equation
    t_span = [0, end_t]
    t_eval = np.linspace(0, end_t, 1000)
    P0 = PPT[0,:]
    sol = Master_equation(t_span, rate_matrix, P0, t_eval)
    t_span = sol.t
    PP = sol.y
    
    fff = open('%s_METS.dat'%mutant_type, 'w')
    for i in range(t_span.shape[0]):
        fff.write('%8.4f '%(t_span[i]))
        for j in range(n_states):
            fff.write('%6.4f '%abs(PP[j,i]))
        fff.write('\n')
    fff.close()
    
np.savez('msm_data.npz', 
         dtrajs = dtrajs,
         center = center,
         eigenvalues_list = eigenvalues_list,
         meta_dtrajs = meta_dtrajs,
         meta_set = meta_set,
         meta_samples = meta_samples,
         rate_matrix_list = rate_matrix_list)
    
# Visualize
if if_visualize:
    os.system('rm -rf state_struct')
    for state_id in range(1, n_states+1):
        frame_id = (state_id-1)*sample_size
        gen_state_visualizion(state_id, psf_file, cor_file, sampled_traj[frame_id], native_AA_pdb, if_entangled_list[state_id-1])
