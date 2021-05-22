#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np
import pyemma as pem
import parmed as pmd
import mdtraj as mdt
import msmtools

################################# Arguments ###################################
# Default values
dt = 0.015/1000
nsave = 5000
alpha = 4331293.0
dt = dt*nsave*alpha/1e9 # in seconds
NCL_div = [1]
n_traj = 100
mutant_type_list = ['fast', 'slow']
n_cluster_per_div = 200
stride=10
n_large_states = 3
n_small_states = 2
lag_t = 1
sample_size = 5
ribo_psf_file = ''
ribo_cor_file = ''
native_AA_pdb = ''
prefix_dir = ''
visualiz_threshold = 0.01
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
        if line.startswith('NCL_div'):
            words = line.split('=')
            words = words[1].strip().split()
            NCL_div = [int(w) for w in words]
            continue
        if line.startswith('n_traj'):
            words = line.split('=')
            n_traj = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('n_cluster_per_div'):
            words = line.split('=')
            n_cluster_per_div = int(words[1].strip())
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
        if line.startswith('sample_size'):
            words = line.split('=')
            sample_size = int(words[1].strip())
            continue
        if line.startswith('ribo_psf'):
            words = line.split('=')
            ribo_psf_file = words[1].strip()
            continue
        if line.startswith('ribo_cor'):
            words = line.split('=')
            ribo_cor_file = words[1].strip()
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

################################# Functions ###################################
def standardize(data):
    data_con = data[0]
    for i in range(1, len(data)):
        data_con = np.vstack((data_con, data[i]))
    data_mean = np.mean(data_con, axis=0)
    data_std = np.std(data_con, axis=0)
    data_std[data_std==0] = 1
    result = [(d - data_mean) / data_std for d in data]
    return [result, data_mean, data_std]

def unstandardize(data, data_mean, data_std):
    result = data * data_std + data_mean
    return result

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
                #G0[idx] = abs(G0[idx])
            #G0 = sum(G0)/4/3.14
            G_list.append([G0[0]/4/3.14, G0[1]/4/3.14])
        else:
            G_list.append([np.nan, np.nan])
    return (M, G_list)
    
def gen_state_visualizion(state_id, native_psf, native_cor, state_psf, state_cor, native_struct, nc_length, if_entangled):
    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP']
    
    print('Generate visualization of state %d'%(state_id))
    #os.system('mkdir state_struct')
    os.chdir('state_struct')
    
    cutoff = 0.8 # in nm
    min_interval = 4
    terminal_cutoff = 5
    
    ## Get native contacts ##
    struct = pmd.load_file(native_psf)
    coor = pmd.load_file(native_cor)
    struct.coordinates = coor.coordinates
    struct = struct[':1-%d'%nc_length]
    
    sel = state_cor.topology.select("residue 1 to %d"%nc_length)
    state_cor[0].atom_slice(sel).save('tmp.pdb', force_overwrite=True)
    
    rnc_struct = pmd.load_file(state_psf)
    rnc_struct.coordinates = state_cor.xyz[0, :, :]*10
    L24_struct = rnc_struct['@B']
    
    ribo_struct_0 = ribo_struct['!@B']
    ribo_struct_0 = L24_struct + ribo_struct_0

    native_contact = []
    for i in range(0, len(struct.atoms)-min_interval):
        coor_1 = np.array([struct[i].xx, struct[i].xy, struct[i].xz])
        for j in range(i+min_interval, len(struct.atoms)):
            coor_2 = np.array([struct[j].xx, struct[j].xy, struct[j].xz])
            dist = np.sum((coor_1-coor_2)**2)**0.5
            if (dist < 10*cutoff):
                native_contact.append([i, j])
    
    _, G_native_list = calc_G_list(struct.coordinates/10, native_contact, cutoff, terminal_cutoff)
    
    M, G_list = calc_G_list(state_cor.xyz[0,np.arange(nc_length).astype(int),:], native_contact, cutoff, terminal_cutoff)
    
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
    native_struct[':1-%d'%nc_length].save('target.pdb')
    os.system('backmap.py -i target.pdb -c tmp.pdb')
    os.system('mv tmp_rebuilt.pdb state_%d.pdb'%state_id)
    os.system('rm -f tmp.pdb')
    os.system('rm -f target.pdb')
    os.system('rm -rf ./rebuild_tmp/')
    
    # combine backmapped NC and ribosome
    nc_struct = pmd.load_file('state_%d.pdb'%state_id)
    rnc_struct = nc_struct + ribo_struct_0
    for i in range(nc_length):
        rnc_struct.residues[i].segid = 'A'
    rnc_struct.save('state_%d.pdb'%state_id, charmm=True, overwrite=True)
    
    # Set visualization via VMD
    f = open('vmd_%d.tcl'%state_id, 'w')
    f.write('''package require topotools
display rendermode GLSL
axes location off
display projection Orthographic

color Display {Background} white
''')

    if if_entangled:
        f.write('''
mol new ./'''+('state_%d.pdb'%state_id)+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.600000 10.000000 3.000000 0
mol color ColorID 6
mol selection {segname A}
mol material AOEdgy
mol addrep top
mol representation NewCartoon 0.610000 10.000000 3.000000 0
mol color ColorID 1
mol selection {resid '''+('%d'%(idx_max[0]))+''' to '''+('%d'%(idx_max[1]))+'''}
mol material AOEdgy
mol addrep top
mol representation NewCartoon 0.610000 10.000000 3.000000 0
mol color ColorID 0
mol selection {resid '''+('%d'%(idx_thread_max[0]))+''' to '''+('%d'%(idx_thread_max[1]))+'''}
mol material AOEdgy
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
''')
    else:
        f.write('''
mol new ./'''+('state_%d.pdb'%state_id)+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.600000 10.000000 3.000000 0
mol color ColorID 10
mol selection {segname A}
mol material AOEdgy
mol addrep top
''')

    f.write('''
mol representation QuickSurf 2.000000 1.000000 3.000000 1.000000
mol color ColorID 8
mol selection {not segname A AtR PtR }
mol material Transparent
mol addrep top
mol representation QuickSurf 2.000000 1.000000 3.000000 1.000000
mol color ColorID 3
mol selection {segname AtR}
mol material Diffuse
mol addrep top
mol representation QuickSurf 2.000000 1.000000 3.000000 1.000000
mol color ColorID 3
mol selection {segname PtR}
mol material Diffuse
mol addrep top

set viewplist {}
set fixedlist {}
set viewpoints([molinfo top]) {{{1 0 0 -25.9619} {0 1 0 -9.60361} {0 0 1 -1.88204} {0 0 0 1}} {{0.579449 0.00158709 -0.815005 0} {0.130733 0.986868 0.0948694 0} {0.804455 -0.16152 0.571631 0} {0 0 0 1}} {{0.0169804 0 0 0} {0 0.0169804 0 0} {0 0 0.0169804 0} {0 0 0 1}} {{1 0 0 -0.98} {0 1 0 0.14} {0 0 1 -1.29} {0 0 0 1}}}
lappend viewplist [molinfo top]
set topmol [molinfo top]
# done with molecule 0
foreach v $viewplist {
  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)
}
foreach v $fixedlist {
  molinfo $v set fixed 1
}
unset viewplist
unset fixedlist
mol top $topmol
unset topmol
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

ribo_struct = pmd.load_file(ribo_psf_file)
ribo_cor = pmd.load_file(ribo_cor_file)
ribo_struct.coordinates = ribo_cor.coordinates

native_struct = pmd.load_file(native_AA_pdb)

G_list_0_list = []
NCL_list = []
T_list = []
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
    qbb_list = list(np.load(mutant_type+'_QBB.npy', allow_pickle=True))
    qbb_list = [q.reshape(q.shape[0],1) for q in qbb_list]
    
    # Get Entanglement trajectory
    G_list_0 = list(np.load(mutant_type+'_Entanglement.npy', allow_pickle=True))
    G_list_0_list += G_list_0
    G_list = [np.sum(g[:,:5], axis=1)/num_nc for g in G_list_0]
    G_list = [g.reshape(g.shape[0],1) for g in G_list]
    
    # NC length trajectories
    T_list_0 = list(np.load(mutant_type+'_T.npy', allow_pickle=True))
    T_list += T_list_0
    NCL_list += [t[:,1].reshape(len(t), 1) for t in T_list_0]
            
    trajid_list += [i for i in range(len(qbb_list))]
    mtype2trajid.append([i+len(cor_list) for i in range(len(qbb_list))])
    trajid2mtype += [i_ax for i in range(len(qbb_list))]
    cor_list += [np.hstack((qbb_list[i], G_list[i])) for i in range(len(qbb_list))]
cor_list = np.array(cor_list)

#Clustering
if if_cluster:
    dtrajs = [[] for i in range(len(cor_list))]
    center = np.zeros((1, cor_list[0].shape[1]))
    for i_div in range(len(NCL_div)):
        start_ncl = NCL_div[i_div]
        if i_div == len(NCL_div)-1:
            end_ncl = max_length
        else:
            end_ncl = NCL_div[i_div+1]-1
        cor_list_0 = []
        for i, cor in enumerate(cor_list):
            idx_1 = np.where(NCL_list[i][:,0] >= start_ncl)[0]
            idx_2 = np.where(NCL_list[i][idx_1,0] <= end_ncl)[0]
            idx_list = idx_1[idx_2]
            cor_list_0.append(cor[idx_list, :])
        
        std_cor_list, cor_mean, cor_std = standardize(cor_list_0)
        cluster = pem.coordinates.cluster_kmeans(std_cor_list, k=n_cluster_per_div, max_iter=5000, stride=stride)
        dtrajs_0 = cluster.dtrajs
        center_0 = unstandardize(cluster.clustercenters, cor_mean, cor_std)
        
        offset = n_cluster_per_div*i_div
        for i in range(len(cor_list)):
            dtrajs[i] += list(dtrajs_0[i] + offset)
        center = np.vstack((center, center_0))
    dtrajs = [np.array(d) for d in dtrajs]
    center = np.delete(center, 0, 0)
else:
    dtrajs = list(npzfile['dtrajs'])
    center = npzfile['center']

n_cluster = n_cluster_per_div*len(NCL_div)

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
    n_states = n_large_states
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
for i in range(2):
    coarse_state_centers[:,i] = np.round(coarse_state_centers[:,i], decimals=4)
cg_center_order_idx = np.argsort(coarse_state_centers[:,0])
print('metastable states centers:')
print(coarse_state_centers[cg_center_order_idx])
micro_to_meta = np.zeros(n_cluster)
meta_set = meta_set[cg_center_order_idx]
meta_dist = meta_dist[cg_center_order_idx, :]
div_to_meta = [[] for i in range(len(NCL_div))]
for idx, ms in enumerate(meta_set):
    for mms in ms:
        micro_to_meta[mms] = idx
    div_to_meta[int(ms[0] / n_cluster_per_div)].append(idx)
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
    
np.savez('msm_data.npz', 
         dtrajs = dtrajs,
         center = center,
         eigenvalues_list = eigenvalues_list,
         meta_dtrajs = meta_dtrajs,
         meta_set = meta_set,
         meta_samples = meta_samples)
    
# Visualize
if if_visualize:
    os.system('rm -rf state_struct')
    os.system('mkdir state_struct')
    for state_id in range(1, n_states+1):
        traj_idx_0 = meta_samples[state_id-1][0, 0]
        frame_idx = meta_samples[state_id-1][0, 1]
        G = cor_list[traj_idx_0][frame_idx, 1]
        if G >= visualiz_threshold:
            if_entangled = True
        else:
            if_entangled = False
        ncl = int(T_list[traj_idx_0][frame_idx, 1])
        co_dir, po_dir, native_psf_file, cor_file = get_co_po_dir(prefix_dir, mutant_type_list[trajid2mtype[traj_idx_0]])
        stage = T_list[traj_idx_0][frame_idx, 2] + 1
        if stage == 4:
            stage = 1
            ncl += 1
        traj_idx = trajid_list[traj_idx_0]
        psf_file = co_dir + '/traj/%d/rnc_l%d.psf'%(traj_idx+1, ncl)
        if frame_idx == len(NCL_list[traj_idx_0])-1:
            traj_file = co_dir + '/traj/%d/rnc_l%d_stage_%d_final.cor'%(traj_idx+1, ncl, stage)
            frame_idx = 0
        elif T_list[traj_idx_0][frame_idx, 2] != T_list[traj_idx_0][frame_idx+1, 2]:
            traj_file = co_dir + '/traj/%d/rnc_l%d_stage_%d_final.cor'%(traj_idx+1, ncl, stage)
            frame_idx = 0
        else:
            traj_file = co_dir + '/traj/%d/rnc_l%d_stage_%d.dcd'%(traj_idx+1, ncl, stage)
            frame_idx_0 = 0
            for i in range(frame_idx, -1, -1):
                if T_list[traj_idx_0][frame_idx, 2] != T_list[traj_idx_0][i, 2]:
                    frame_idx_0 = frame_idx - i - 1
                    break
                elif i == 0:
                    frame_idx_0 = frame_idx
                    break
            frame_idx = frame_idx_0
        if traj_file.endswith('.cor'):
            struct = pmd.load_file(psf_file)
            struct.coordinates = pmd.load_file(traj_file).coordinates
            struct.save('state_struct/tmp.pdb', overwrite=True)
            traj_cor = mdt.load('state_struct/tmp.pdb', top=psf_file)[frame_idx]
        else:
            traj_cor = mdt.load(traj_file, top=psf_file)[frame_idx]
                
        gen_state_visualizion(state_id, native_psf_file, cor_file, psf_file, traj_cor, native_struct, ncl, if_entangled)
