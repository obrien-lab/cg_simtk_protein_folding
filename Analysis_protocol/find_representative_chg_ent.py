#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback, glob, copy
import pickle
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage, cophenet
from scipy.spatial.distance import squareform
import parmed as pmd
import mdtraj as mdt
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
# matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.labelsize'] = 'small'
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 4
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['xtick.labelsize'] = 'x-small'
matplotlib.rcParams['ytick.labelsize'] = 'x-small'
matplotlib.rcParams['legend.fontsize'] = 'x-small'
matplotlib.rcParams['figure.dpi'] = 600

usage = '''
  Usage: python find_representative_chg_ent.py 
                --data | -d <.PKL FILES> for entanglement analysis
                --psf | -p <INPUT.PSF> for the protein structure 
                --traj_dir_prefix | -t <TRAJECTORY> directory prefix for 
                                       entanglement analysis.
                [--classify_key | -k] <KEY> for classifying changes of entanglements.
                                      Default "linking_number".
                [--cluster_data | -c] <CLUSTER.NPZ> for pre-calculated clusters.
                                      Default is None, do clustering from -d data.
                [--viz | -v] <0 or 1> 0: Do not generate VMD visualization 
                             scripts; 1: Generate VMD visualization scripts; 
                             Default 1.
                [--top_struct | -s] <CUTOFF> for unique entangled structure visualization.
                                    An integer >= 1 means the number of top structures;
                                    A number < 1 means the probability cutoff 
                                    (skip lower probability);
                [--backmap | -b] <0 or 1> 0: Do not backmap; 1: Backmap;
                                 Default 1. Invalid when --viz=0.
                [--ref_AA_pdb | -r] <REFERENCE.PDB> for backmapping; 
                                    Invalid when --viz=0. No default value.
'''

################### Functions #########################
def save_pickle(filename, mode, data, protocol=4):
    with open(filename, mode) as fh:
        pickle.dump(data, fh, protocol=protocol)

def load_pickle(filename):
    data = []
    with open(filename, 'rb') as fr:
        try:
            while True:
                data.append(pickle.load(fr))
        except EOFError:
            pass
    data_dict = {}
    for da in data:
        data_dict.update(da)
    return data_dict

def pdist_loop_overlap(data_array_1, data_array_2):
    M_1 = np.repeat(data_array_1.reshape((data_array_1.shape[0], 1, data_array_1.shape[1])), data_array_2.shape[0], axis=1)
    M_2 = np.repeat(data_array_2.reshape((1, data_array_2.shape[0], data_array_2.shape[1])), data_array_1.shape[0], axis=0)
    M = np.concatenate((M_1, M_2), axis=-1)
    del M_1, M_2
    dist_M = (np.max(M, axis=-1) - np.min(M, axis=-1) + 1) / (M[:,:,1]-M[:,:,0]+M[:,:,3]-M[:,:,2])
    # Make distance between inclusive loops to be minimum
    dist_M[(M[:,:,2]-M[:,:,0])*(M[:,:,1]-M[:,:,3]) >= 0] = 0.5 
    return dist_M

def pdist_thread_deviation(data_array_1, data_array_2):
    M_1 = np.repeat(data_array_1.reshape((data_array_1.shape[0], 1, data_array_1.shape[1])), data_array_2.shape[0], axis=1)
    M_2 = np.repeat(data_array_2.reshape((1, data_array_2.shape[0], data_array_2.shape[1])), data_array_1.shape[0], axis=0)
    M = np.concatenate((M_1, M_2), axis=-1)
    del M_1, M_2
    dist_M = np.abs((M[:,:,2:4]-M[:,:,0:2]))
    dist_M[M[:,:,2:4]*M[:,:,0:2] < 0] = 10 # Make distance between no crossing and crossings to be small
    del M
    dist = np.max(dist_M, axis=-1)
    return dist
    
def pdist_cross_contamination(data_array_1, data_array_2):
    # data looks like [nc_1, nc_2, cross_N1, cross_N2, ..., cross_C1, cross_C2, ...]
    M_1 = np.repeat(data_array_1.reshape((data_array_1.shape[0], 1, data_array_1.shape[1])), data_array_2.shape[0], axis=1)
    M_2 = np.repeat(data_array_2.reshape((1, data_array_2.shape[0], data_array_2.shape[1])), data_array_1.shape[0], axis=0)
    M = np.concatenate((M_1, M_2), axis=-1)
    del M_1, M_2

    # Distance for cross_2 contaminate loop_1
    idx_array_1 = np.zeros((data_array_2.shape[1]-2,2), dtype=int)
    idx_array_1[:,0] = 1
    idx_array_1[:,1] = np.arange(data_array_1.shape[1]+2, data_array_1.shape[1]+data_array_2.shape[1], dtype=int)
    idx_array_2 = np.zeros((data_array_2.shape[1]-2,2), dtype=int)
    idx_array_2[:,0] = np.arange(data_array_1.shape[1]+2, data_array_1.shape[1]+data_array_2.shape[1], dtype=int)
    idx_array_2[:,1] = 0
    L = (M[:,:,1]-M[:,:,0]).reshape((M.shape[0], M.shape[1], 1))

    dist_M_1 = np.min(M[:,:,idx_array_1]-M[:,:,idx_array_2], axis=-1) / L
    dist_M_1[dist_M_1 <= 0] = 0
    dist_M_1[dist_M_1 >= 1] = 0
    
    # Distance for cross_1 contaminate loop_2
    idx_array_1 = np.zeros((data_array_1.shape[1]-2,2), dtype=int)
    idx_array_1[:,0] = data_array_1.shape[1]+1
    idx_array_1[:,1] = np.arange(2, data_array_1.shape[1], dtype=int)
    idx_array_2 = np.zeros((data_array_1.shape[1]-2,2), dtype=int)
    idx_array_2[:,0] = np.arange(2, data_array_1.shape[1], dtype=int)
    idx_array_2[:,1] = data_array_1.shape[1]
    L = (M[:,:,data_array_1.shape[1]+1]-M[:,:,data_array_1.shape[1]]).reshape((M.shape[0], M.shape[1], 1))

    dist_M_2 = np.min(M[:,:,idx_array_1]-M[:,:,idx_array_2], axis=-1) / L
    dist_M_2[dist_M_2 <= 0] = 0
    dist_M_2[dist_M_2 >= 1] = 0
    
    del M
    dist_M = np.max(np.concatenate((dist_M_1, dist_M_2), axis=-1), axis=-1)
    return dist_M

def agglomerative_clustering(dist, cluster_method, cluster_dist_cutoff, num_perm):
    min_SSDIFN = np.inf
    best_Z = None
    best_perm_idx_list = None
    pdist = squareform(dist, checks=False)
    if np.sum(pdist**2) == 0:
        best_Z = linkage(pdist, method=cluster_method)
        best_perm_idx_list = np.arange(dist.shape[0])
    else:
        # permuCLUSTER
        for idx_perm in range(np.max([1, num_perm])):
            perm_idx_list = np.random.permutation(np.arange(dist.shape[0]))
            pdist = squareform(dist[perm_idx_list,:][:,perm_idx_list], checks=False)
            Z = linkage(pdist, method=cluster_method)
            cdist = cophenet(Z)
            SSDIFN = np.sum((pdist - cdist)**2)/np.sum(pdist**2)
            if SSDIFN < min_SSDIFN:
                min_SSDIFN = SSDIFN
                best_Z = Z
                best_perm_idx_list = perm_idx_list
    cluster_id_list = fcluster(best_Z, cluster_dist_cutoff, criterion='distance')
    backmap_list = np.zeros(len(best_perm_idx_list), dtype=int)
    for i, j in enumerate(best_perm_idx_list):
        backmap_list[j] = i
    cluster_id_list = cluster_id_list[backmap_list]
    return cluster_id_list

def do_clustering(map_list, chg_ent_fingerprint_list, key, pdist_fun, cluster_method, cluster_dist_cutoff, num_perm=100):
    data = []
    # Get max number of crossings
    max_n_cross = 1
    if 'cross_contamination' in key:
        for map_idx in map_list:
            fingerprint = chg_ent_fingerprint_list[map_idx[0]][map_idx[1]][tuple(map_idx[2:4])]
            cr = fingerprint['crossing_resid']
            for ci, c in enumerate(cr):
                if len(c) > max_n_cross:
                    max_n_cross = len(c)
    # Prepare clustering data for distance calculation
    for map_idx in map_list:
        fingerprint = chg_ent_fingerprint_list[map_idx[0]][map_idx[1]][tuple(map_idx[2:4])]
        if 'crossing_resid' in key:
            ter_idx = int(key.split('_')[-1])
            mc = []
            cr = fingerprint['crossing_resid']
            ref_cr = fingerprint['ref_crossing_resid']
            for c in [ref_cr[ter_idx], cr[ter_idx]]:
                if len(c) == 0:
                    mc.append(-1)
                else:
                    mc.append(np.median(c))
            data.append(mc)
        elif 'native_contact' in key:
            nc = fingerprint['native_contact']
            data.append(nc)
        elif 'cross_contamination' in key:
            nc = fingerprint['native_contact']
            mc = []
            cr = fingerprint['crossing_resid']
            for ci, c in enumerate(cr):
                for cii in range(max_n_cross):
                    if cii >= len(c):
                        mc.append(-1)
                    else:
                        mc.append(c[cii])
            data.append(nc + mc)
        else:
            print('Error: Unknown key specified for do_clustering(), %s'%(key))
            sys.exit()
    data = np.array(data)
    # Reduce data size (saving memory usage) by combining duplicated data points
    reduced_data = np.unique(data, axis=0)
    reduced_data_map = np.array([np.all(data == d, axis=1).nonzero()[0].tolist() for d in reduced_data], dtype=object)
    
    # If all data are the same, group them into a single cluster and return
    if len(reduced_data) == 1:
        cluster_data = [map_list]
        return cluster_data
    
    # Chunk data to reduce memory usage if the expanded matrix occupy >= memory_cutoff
    if reduced_data.nbytes ** 2 >= memory_cutoff:
        n_chunk = int(np.ceil(reduced_data.nbytes / np.sqrt(memory_cutoff/2)))
        len_chunk = int(np.ceil(len(reduced_data) / n_chunk))
        dist = np.zeros((len(reduced_data), len(reduced_data)))
        for i in range(n_chunk):
            i_1 = i*len_chunk
            i_2 = np.min([(i+1)*len_chunk,len(reduced_data)])
            for j in range(n_chunk):
                j_1 = j*len_chunk
                j_2 = np.min([(j+1)*len_chunk,len(reduced_data)])
                dist[i_1:i_2, j_1:j_2] = pdist_fun(reduced_data[i_1:i_2], reduced_data[j_1:j_2])
    else:
        dist = pdist_fun(reduced_data, reduced_data)
    
    if 'cross_contamination' in key:
        # Do divisive clustering
        cluster_idx_mapping = [list(np.arange(dist.shape[0]))]
        while True:
            cluster_1 = cluster_idx_mapping[-1]
            cluster_2 = []
            cluster_0 = copy.deepcopy(cluster_1)
            for i in range(len(cluster_0)-1):
                rm_idx_list = np.where(dist[cluster_0[i],cluster_0[i+1:]] >= cluster_dist_cutoff)[0]
                if len(rm_idx_list) > 0:
                    cluster_1.remove(cluster_0[i])
                    cluster_2.append(cluster_0[i])
            if len(cluster_2) > 0:
                cluster_idx_mapping.append(cluster_2)
            else:
                break
        # Do agglomerative clustering
        if len(cluster_idx_mapping) > 1:
            dist_0 = np.zeros((len(cluster_idx_mapping),len(cluster_idx_mapping)))
            for i in range(len(dist_0)-1):
                for j in range(i+1, len(dist_0)):
                    dist_0[i,j] = np.max(dist[cluster_idx_mapping[i],:][:,cluster_idx_mapping[j]])
            cluster_id_list_0 = agglomerative_clustering(dist_0, cluster_method, cluster_dist_cutoff, num_perm)
        else:
            cluster_id_list_0 = np.array([1], dtype=int)
        cluster_id_list = np.zeros(dist.shape[0], dtype=int)
        for cluster_idx, mapping in enumerate(cluster_idx_mapping):
            cluster_id_list[mapping] = cluster_id_list_0[cluster_idx]
    else:
        # Do agglomerative clustering
        cluster_id_list = agglomerative_clustering(dist, cluster_method, cluster_dist_cutoff, num_perm)
    n_cluster = np.max(cluster_id_list)
    
    # Back-Mapping indices 
    cluster_data = []
    for cluster_id in range(n_cluster):
        idx = np.where(cluster_id_list == cluster_id+1)[0]
        idx_list = reduced_data_map[idx].tolist()
        idx_list_1 = []
        for i in idx_list:
            idx_list_1 += i
        idx_list_1 = sorted(idx_list_1)
        cluster_data.append(np.array(map_list)[idx_list_1].tolist())
    return cluster_data

def cluster_chg_ent(chg_ent_keyword_dict, chg_ent_fingerprint_list, cluster_method=['average', 'average', 'complete'], cluster_dist_cutoff=[20, 1.0, 0.1]):
    cluster_data_keys = sorted(list(chg_ent_keyword_dict.keys()))
    cluster_data = {key: [] for key in cluster_data_keys}
    cluster_tree = {key: [] for key in cluster_data_keys}
    for key in cluster_data_keys:
        map_list = chg_ent_keyword_dict[key]

        # if len(map_list) == 0:
            # continue
        # elif len(map_list) == 1:
            # cluster_data[key].append([map_list[0]])
            # print('Found 1 cluster(s) for %s'%(key))
            # continue
        
        # First clustering based on N-ter crossing residues 
        N_cr_cluster_data = do_clustering(map_list, chg_ent_fingerprint_list, 'crossing_resid_0', 
                                          pdist_thread_deviation, cluster_method[0], cluster_dist_cutoff[0])
        
        backtrace_idx_list = []
        idx_1 = 0
        idx_2 = 0
        idx_3 = 0
        # Second clustering based on C-ter crossing residues 
        for map_list_N_cr in N_cr_cluster_data:
            # if len(map_list_N_cr) == 0:
                # continue
            # elif len(map_list_N_cr) == 1:
                # cluster_data[key].append([map_list_N_cr[0]])
                # continue
            C_cr_cluster_data = do_clustering(map_list_N_cr, chg_ent_fingerprint_list, 'crossing_resid_1', 
                                              pdist_thread_deviation, cluster_method[0], cluster_dist_cutoff[0])
        
            # Third clustering based on loop
            for map_list_C_cr in C_cr_cluster_data:
                # if len(map_list_C_cr) == 0:
                    # continue
                # elif len(map_list_C_cr) == 1:
                    # cluster_data[key].append([map_list_C_cr[0]])
                    # continue
                nc_cluster_data = do_clustering(map_list_C_cr, chg_ent_fingerprint_list, 'native_contact', 
                                                pdist_loop_overlap, cluster_method[1], cluster_dist_cutoff[1])
                
                # Fourth clustering based on cross contamination
                for map_list_nc in nc_cluster_data:
                    # if len(map_list_nc) == 0:
                        # continue
                    # elif len(map_list_nc) == 1:
                        # cluster_data[key].append([map_list_nc[0]])
                        # continue
                    final_cluster_data = do_clustering(map_list_nc, chg_ent_fingerprint_list, 'cross_contamination', 
                                                       pdist_cross_contamination, cluster_method[2], cluster_dist_cutoff[2])
                    for final_cluster in final_cluster_data:
                        cluster_data[key].append(final_cluster)
                        backtrace_idx_list.append([idx_1, idx_2, idx_3])
                    idx_3 += 1
                idx_2 += 1
            idx_1 += 1
        
        backtrace_idx_list = np.array(backtrace_idx_list, dtype=int)
        cluster_tree[key] = [[np.where(backtrace_idx_list[:,i]==j)[0].tolist() for j in range(backtrace_idx_list[:,i].max()+1)] for i in range(backtrace_idx_list.shape[1])]

        n_cluster = len(cluster_data[key])
        print('Found %d cluster(s) for %s'%(n_cluster, key))
    
    return (cluster_data, cluster_tree)

def find_representative_entanglement(cluster_data, ent_cluster_idx_map):
    # most probable loop midpoint
    rep_ent_list = []
    for [key, idx] in ent_cluster_idx_map:
        cluster = np.array(cluster_data[key][idx])
        loop_midpoint_list = np.mean(cluster[:,2:4], axis=1)
        loop_max = np.max(cluster[:,3])
        loop_min = np.min(cluster[:,2])
        # mode
        bins = np.arange(loop_min, loop_max, 5)
        if len(bins) == 1:
            bins = np.arange(loop_min, loop_max, 1)
        hist, edges = np.histogram(loop_midpoint_list, bins=bins)
        idx = np.argmax(hist)
        min_loop_len = 1e6
        idx0 = np.where(loop_midpoint_list >= edges[idx])[0]
        idx1 = np.where(loop_midpoint_list[idx0] < edges[idx+1])[0]
        for iidx in idx0[idx1]:
            if cluster[iidx,3]-cluster[iidx,2] < min_loop_len:
                rep_ent = cluster[iidx]
                min_loop_len = rep_ent[3] - rep_ent[2]
        rep_ent_list.append(rep_ent)
    return rep_ent_list

def gen_state_visualizion(state_id, psf, state_cor, native_AA_pdb, rep_ent_dict, if_backmap=True, pulchra_only=False):
    def idx2sel(idx_list):
        if len(idx_list) == 0:
            return ''
        else:
            sel = 'index'
            idx_0 = idx_list[0]
            idx_1 = idx_list[0]
            sel_0 = ' %d'%idx_0
            for i in range(1, len(idx_list)):
                if idx_list[i] == idx_list[i-1] + 1:
                    idx_1 = idx_list[i]
                else:
                    if idx_1 > idx_0:
                        sel_0 += ' to %d'%idx_1
                    sel += sel_0
                    idx_0 = idx_list[i]
                    idx_1 = idx_list[i]
                    sel_0 = ' %d'%idx_0
            if idx_1 > idx_0:
                sel_0 += ' to %d'%idx_1
            sel += sel_0
            return sel

    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP']
    protein_colorid_list = [6, 6]
    loop_colorid_list = [1, 1]
    thread_colorid_list = [0, 0]
    nc_colorid_list = [3, 3]
    crossing_colorid_list = [8, 8]
    thread_cutoff=3
    terminal_cutoff=3
    
    print('Generate visualization of state %d'%(state_id))

    struct = pmd.load_file(psf)
    struct.coordinates = state_cor

    # backmap
    if if_backmap:
        if pulchra_only:
            pulchra_only = '1'
        else:
            pulchra_only = '0'
        struct.save('tmp.pdb', overwrite=True)
        os.system('backmap.py -i '+native_AA_pdb+' -c tmp.pdb -p '+pulchra_only)
        os.system('mv tmp_rebuilt.pdb state_%d.pdb'%state_id)
        os.system('rm -f tmp.pdb')
        os.system('rm -rf ./rebuild_tmp/')
    else:
        struct.save('state_%d.pdb'%state_id, overwrite=True)
    
    ref_struct = pmd.load_file(native_AA_pdb)
    current_struct = pmd.load_file('state_%d.pdb'%state_id)

    if len(list(rep_ent_dict.keys())) == 0:
        # no change of entaglement
        f = open('vmd_s%d_none.tcl'%(state_id), 'w')
        f.write('# Entanglement type: no change\n')
        f.write('''display rendermode GLSL
axes location off

color Display {Background} white

mol new ./'''+('state_%d.pdb'%state_id)+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID '''+str(protein_colorid_list[1])+'''
mol selection {all}
mol material AOChalky
mol addrep top
''')
        f.close()

    # Create vmd script for each type of change
    for ent_code, rep_ent_list in rep_ent_dict.items():
        pmd_struct_list = [ref_struct, current_struct]
        struct_dir_list = [native_AA_pdb, './state_%d.pdb'%state_id]
        key_prefix_list = ['ref_', '']
        repres_list = ['', '']
        align_sel_list = ['', '']
        
        vmd_script = '''# Entanglement type: '''+str(rep_ent_list[0]['type'])+'''
package require topotools
display rendermode GLSL
axes location off

color Display {Background} white

'''
        for struct_idx, pmd_struct in enumerate(pmd_struct_list):
            struct_dir = struct_dir_list[struct_idx]
            protein_colorid = protein_colorid_list[struct_idx]
            loop_colorid = loop_colorid_list[struct_idx]
            thread_colorid = thread_colorid_list[struct_idx]
            nc_colorid = nc_colorid_list[struct_idx]
            crossing_colorid = crossing_colorid_list[struct_idx]
            key_prefix = key_prefix_list[struct_idx]

            # Clean ligands
            clean_sel_idx = np.zeros(len(pmd_struct.atoms))
            for res in pmd_struct.residues:
                if res.name in AA_name_list:
                    for atm in res.atoms:
                        clean_sel_idx[atm.idx] = 1
            pmd_clean_struct = pmd_struct[clean_sel_idx]
            clean_idx_to_idx = np.where(clean_sel_idx == 1)[0]

            # vmd selection string for protein
            idx_list = []
            for res in pmd_struct.residues:
                if res.name in AA_name_list:
                    idx_list += [atm.idx for atm in res.atoms]
            vmd_sel = idx2sel(idx_list)

            repres = '''mol new '''+struct_dir+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID '''+str(protein_colorid)+'''
mol selection {'''+vmd_sel+'''}
mol material AOChalky
mol addrep top
'''
            align_sel = vmd_sel
            for chg_ent_fingerprint in rep_ent_list:
                nc = chg_ent_fingerprint[key_prefix+'native_contact']

                idx_list = []
                for res in pmd_clean_struct.residues:
                    if res.idx in nc:
                        idx_list += [atm.idx for atm in res.atoms if atm.name == 'CA']
                nc_sel = idx2sel(clean_idx_to_idx[idx_list])

                idx_list = []
                for res in pmd_clean_struct.residues:
                    if res.idx >= nc[0] and res.idx <= nc[1]:
                        idx_list += [atm.idx for atm in res.atoms]
                loop_sel = idx2sel(clean_idx_to_idx[idx_list])

                align_sel += ' and not (%s)'%loop_sel
                ref_coss_resid = chg_ent_fingerprint['ref_crossing_resid']
                cross_resid = chg_ent_fingerprint['crossing_resid']
                thread = []
                thread_sel_list = []
                for ter_idx in range(len(ref_coss_resid)):
                    thread_0 = []
                    resid_list = ref_coss_resid[ter_idx] + cross_resid[ter_idx]
                    if len(resid_list) > 0:
                        thread_0 = [np.min(resid_list)-5, np.max(resid_list)+5]
                        if ter_idx == 0:
                            thread_0[0] = np.max([thread_0[0], terminal_cutoff])
                            thread_0[1] = np.min([thread_0[1], nc[0]-thread_cutoff])
                        else:
                            thread_0[0] = np.max([thread_0[0], nc[1]+thread_cutoff])
                            thread_0[1] = np.min([thread_0[1], len(struct.atoms)-1-terminal_cutoff])
                        idx_list = []
                        for res in pmd_clean_struct.residues:
                            if res.idx >= thread_0[0] and res.idx <= thread_0[1]:
                                idx_list += [atm.idx for atm in res.atoms]
                        thread_0_sel = idx2sel(clean_idx_to_idx[idx_list])
                        thread_sel_list.append(thread_0_sel)
                        align_sel += ' and not (%s)'%thread_0_sel
                    else:
                        thread_sel_list.append('')
                    thread.append(thread_0)

                ln = chg_ent_fingerprint[key_prefix+'linking_number']
                cross = []
                for i in range(len(chg_ent_fingerprint[key_prefix+'crossing_resid'])):
                    cross.append([])
                    for j in range(len(chg_ent_fingerprint[key_prefix+'crossing_resid'][i])):
                        cross[-1].append(chg_ent_fingerprint[key_prefix+'crossing_pattern'][i][j]+str(chg_ent_fingerprint[key_prefix+'crossing_resid'][i][j]))
                repres += '# idx: native contact %s, linking number %s, crossings %s.\n'%(str(nc), str(ln), str(cross))
                repres +='''mol representation NewCartoon 0.350000 10.000000 4.100000 0
mol color ColorID '''+str(loop_colorid)+'''
mol selection {'''+loop_sel+'''}
mol material Opaque
mol addrep top
mol representation VDW 1.000000 12.000000
mol color ColorID '''+str(nc_colorid)+'''
mol selection {'''+nc_sel+'''}
mol material Opaque
mol addrep top
set sel [atomselect top "'''+nc_sel+'''"]
set idx [$sel get index]
topo addbond [lindex $idx 0] [lindex $idx 1]
mol representation Bonds 0.300000 12.000000
mol color ColorID '''+str(nc_colorid)+'''
mol selection {'''+nc_sel+'''}
mol material Opaque
mol addrep top
'''
                for ter_idx, thread_resid in enumerate(thread):
                    if len(thread_resid) == 0:
                        continue
                    repres += '''mol representation NewCartoon 0.350000 10.000000 4.100000 0
mol color ColorID '''+str(thread_colorid)+'''
mol selection {'''+thread_sel_list[ter_idx]+'''}
mol material Opaque
mol addrep top
'''
                    if len(chg_ent_fingerprint[key_prefix+'crossing_resid'][ter_idx]) > 0:
                        idx_list = []
                        for res in pmd_clean_struct.residues:
                            if res.idx in chg_ent_fingerprint[key_prefix+'crossing_resid'][ter_idx]:
                                idx_list += [atm.idx for atm in res.atoms if atm.name == 'CA']
                        crossing_sel = idx2sel(clean_idx_to_idx[idx_list])
                        repres += '''mol representation VDW 1.000000 12.000000
mol color ColorID '''+str(crossing_colorid)+'''
mol selection {'''+crossing_sel+'''}
mol material Opaque
mol addrep top
'''
                
            if struct_idx == 0:
                repres += '''mol representation VDW 1.000000 12.000000
mol color Name
mol selection {not ('''+vmd_sel+''') and not water}
mol material Opaque
mol addrep top
'''
            repres_list[struct_idx] = repres
            align_sel_list[struct_idx] = align_sel

        vmd_script += '\n'.join(repres_list)
        vmd_script += '''
set sel1 [atomselect 0 "'''+align_sel_list[0]+''' and name CA"]
set sel2 [atomselect 1 "'''+align_sel_list[1]+''' and name CA"]
set trans_mat [measure fit $sel1 $sel2]
set move_sel [atomselect 0 "all"]
$move_sel move $trans_mat
'''
        f = open('vmd_s%d_n%s_c%s.tcl'%(state_id, ent_code[0], ent_code[1]), 'w')
        f.write(vmd_script)
        f.close()
    
###################### MAIN ###########################
if len(sys.argv) == 1:
    print(usage)
    sys.exit()

pkl_file_path = None
psf_file = None
traj_dir_prefix = None
if_viz = 1
top_struct = 0.01
if_backmap = 1
pulchra_only = False
native_AA_pdb = None
npz_data_file = None
classify_key = ['linking_number']
cluster_method = ['average', 'average', 'complete']
# cluster_dist_cutoff = [20, 1.0, 0.6] # Allow contamination
cluster_dist_cutoff = [20, 1.0, 0.1] # No contamination
memory_cutoff = 6.4e10 # 64 Gb
max_plot_samples = 1000

try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "h:d:p:t:k:c:v:s:b:r:", 
                               ["help", "data=", "psf=", "traj_dir_prefix=", "classify_key=", "cluster_data=", "viz=", "top_struct=", "backmap=", "ref_AA_pdb="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt in ("-d", "--data"):
        pkl_file_path = os.path.abspath(arg)
    elif opt in ("-p", "--psf"):
        psf_file = os.path.abspath(arg)
    elif opt in ("-t", "--traj_dir_prefix"):
        traj_dir_prefix = os.path.abspath(arg)
    elif opt in ("-k", "--classify_key"):
        classify_key = arg.strip().split()
    elif opt in ("-c", "--cluster_data"):
        npz_data_file = os.path.abspath(arg)
    elif opt in ("-v", "--viz"):
        if_viz = int(arg)
    elif opt in ("-s", "--top_struct"):
        top_struct = float(arg)
    elif opt in ("-b", "--backmap"):
        if_backmap = int(arg)
    elif opt in ("-r", "--ref_AA_pdb"):
        native_AA_pdb = os.path.abspath(arg)

tag_err = False
if pkl_file_path is None and npz_data_file is None:
    print('Error: No pkl data file specified.')
    tag_err = True
if psf_file is None:
    print('Error: No psf file specified.')
    tag_err = True
if traj_dir_prefix is None:
    print('Error: No trajectory directory prefix specified.')
    tag_err = True
if if_viz != 0 and if_viz != 1:
    print('Error: -v must be either 0 or 1.')
    tag_err = True
if if_backmap != 0 and if_backmap != 1:
    print('Error: -b must be either 0 or 1.')
    tag_err = True
if native_AA_pdb is None and if_viz == 1:
    print('Error: -r must be specified when --viz=1.')
    tag_err = True
    
if tag_err:
    print(usage)
    sys.exit()

if if_viz == 1:
    if_viz = True
else:
    if_viz = False
if if_backmap == 1:
    if_backmap = True
else:
    if_backmap = False

if not pkl_file_path is None:
    ent_data_file_list = glob.glob(pkl_file_path)
    if len(ent_data_file_list) == 0:
        print('Error: No files found in the path "%s"'%pkl_file_path)
        sys.exit()

struct = pmd.load_file(psf_file)
n_residue = len(struct.residues)

if npz_data_file is None:
    # Classify changes of entanglement based on the keyword 
    # "[change_code, classify_key_1_N, classify_key_1_C, classify_key_2_N, classify_key_2_C, ...]"
    print('Reading pkl data and classify changes of entanglement...')
    chg_ent_fingerprint_list = []
    Q_list = []
    chg_ent_keyword_dict = {}
    chg_ent_keyword_list = []
    idx2frame = []
    idx2trajfile = []
    dtrajs = []
    # combined_traj = None
    for traj_idx, ent_data_file in enumerate(ent_data_file_list):
        proc_str = 'Processing %s (%6d / %6d)...'%(ent_data_file, traj_idx+1, len(ent_data_file_list))
        print(proc_str, end='\r')
        chg_ent_fingerprint_list.append({})
        Q_list.append({})
        ent_data = load_pickle(ent_data_file)
        frame_list = list(ent_data.keys())
        # sort frames, for old python version that doesn't support ordered dictionary
        frame_list = sorted([frame for frame in frame_list if frame != 'ref'])
        dtrajs.append([[] for frame in frame_list])
        idx2frame.append(frame_list)
        traj_file = traj_dir_prefix + '/' + ent_data_file.split('/')[-1].split('.')[0][:-4] + '.dcd'
        idx2trajfile.append(traj_file)
        # if combined_traj is None:
            # combined_traj = mdt.load(traj_file, top=psf_file)[frame_list]
        # else:
            # combined_traj += mdt.load(traj_file, top=psf_file)[frame_list]
        for frame in frame_list:
            chg_ent_fingerprint_list[-1][frame] = {}
            Q_list[-1][frame] = np.sum(list(ent_data[frame]['G_dict'].values())) / len(list(ent_data['ref']['ent_fingerprint'].keys())) / 2
            for nc, fingerprint in ent_data[frame]['chg_ent_fingerprint'].items():
                if fingerprint['type'] == ['no change', 'no change']: # No change of entanglement for both termini
                    continue
                chg_ent_fingerprint_list[-1][frame][nc] = fingerprint
                chg_ent_keyword = fingerprint['code'].copy()
                for ck in classify_key:
                    if type(fingerprint[ck]) == list:
                        chg_ent_keyword += fingerprint[ck]
                    else:
                        chg_ent_keyword += [fingerprint[ck]]
                chg_ent_keyword = str(chg_ent_keyword)
                if chg_ent_keyword not in chg_ent_keyword_list:
                    chg_ent_keyword_dict[chg_ent_keyword] = []
                    chg_ent_keyword_list.append(chg_ent_keyword)
                chg_ent_keyword_dict[chg_ent_keyword].append([traj_idx, frame] + list(nc))
        print(' '*len(proc_str), end='\r')
    print('%d data files have been read.'%len(ent_data_file_list))
    
    # # Calculate pairewise RMSD
    # combined_traj_idx_list = [[traj_idx, frame] for traj_idx, frame_list in enumerate(idx2frame) for frame in frame_list]
    # RMSD_array = np.zeros((len(combined_traj_idx_list), len(combined_traj_idx_list)))
    # for i in range(len(combined_traj_idx_list)-1):
        # rmsd_list = mdt.rmsd(combined_traj[i:], combined_traj, frame=i, precentered=True)
        # RMSD_array[i, i:] = rmsd_list
        # RMSD_array[i:, i] = rmsd_list

    # cluster changes of entanglements found in the trajectories
    print('Clustering changes of entanglement for %d keywords...'%(len(chg_ent_keyword_list)))
    ent_cluster_data, ent_cluster_tree = cluster_chg_ent(chg_ent_keyword_dict, chg_ent_fingerprint_list, cluster_method=cluster_method, 
                                                         cluster_dist_cutoff=cluster_dist_cutoff)
    chg_ent_keyword_list = sorted(chg_ent_keyword_list)
    # Save calculted data in case job is unexpectedly terminated
    np.savez('cluster_data_%s.npz'%('_'.join(classify_key)),
             chg_ent_fingerprint_list=chg_ent_fingerprint_list,
             Q_list=Q_list,
             chg_ent_keyword_dict=chg_ent_keyword_dict,
             chg_ent_keyword_list=chg_ent_keyword_list,
             idx2trajfile=idx2trajfile,
             idx2frame=idx2frame,
             # RMSD_array=RMSD_array,
             ent_cluster_data=ent_cluster_data,
             ent_cluster_tree=ent_cluster_tree)
else:
    print('Reading clustering data from %s...'%npz_data_file)
    npz_data = np.load(npz_data_file, allow_pickle=True)
    chg_ent_fingerprint_list = npz_data['chg_ent_fingerprint_list'].tolist()
    Q_list = npz_data['Q_list'].tolist()        
    chg_ent_keyword_dict = npz_data['chg_ent_keyword_dict'].item()
    chg_ent_keyword_list = npz_data['chg_ent_keyword_list'].tolist()
    idx2frame = npz_data['idx2frame'].tolist()
    idx2trajfile = npz_data['idx2trajfile'].tolist()
    dtrajs = [[[] for frame in chg_ent_fingerprint.keys()] for chg_ent_fingerprint in chg_ent_fingerprint_list]
    ent_cluster_data = npz_data['ent_cluster_data'].item()
    ent_cluster_tree = npz_data['ent_cluster_tree'].item()
    # # cluster changes of entanglements found in the trajectories
    # print('Clustering changes of entanglement for %d keywords...'%(len(chg_ent_keyword_list)))
    # ent_cluster_data, ent_cluster_tree = cluster_chg_ent(chg_ent_keyword_dict, chg_ent_fingerprint_list, cluster_method=cluster_method, 
                                                         # cluster_dist_cutoff=cluster_dist_cutoff)
    # combined_traj_idx_list = [[traj_idx, frame] for traj_idx, frame_list in enumerate(idx2frame) for frame in frame_list]
    # RMSD_array = npz_data['RMSD_array']

ent_cluster_idx_map = []
for ent_keyword in chg_ent_keyword_list:
    for i  in range(len(ent_cluster_data[ent_keyword])):
        ent_cluster_idx_map.append([ent_keyword, i])

# Print cluster tree
cluster_headers = ['After clustering on N crossing', 'After clustering on C crossing', 'After clustering on loop']
with open('cluster_tree_%s.dat'%('_'.join(classify_key)), 'w') as f:
    for ent_keyword in chg_ent_keyword_list:
        f.write(ent_keyword+'\n')
        clusters = ent_cluster_tree[ent_keyword]
        for i in range(len(clusters)):
            f.write(' '*4 + cluster_headers[i] + ':\n')
            for cluster in clusters[i]:
                f.write(' '*8 + '[')
                for ci, c in enumerate(cluster):
                    cluster_id = ent_cluster_idx_map.index([ent_keyword, c])+1
                    if ci == 0:
                        f.write('%d'%cluster_id)
                    else:
                        f.write(', %d'%cluster_id)
                f.write(']\n')
        f.write('\n')

# Find representative changes of entanglement in each cluster
print('Finding representative changes of entanglement...')
rep_chg_ent_list = find_representative_entanglement(ent_cluster_data, ent_cluster_idx_map)
# Create dataframe and save data
data = []
column_list = ['Keywords', 'Trajectory', 'Frame', 'Native Contact (Residue Index)',
               'Ref N-ter Crossing', 'Ref C-ter Crossing', 'N-ter Crossing', 'C-ter Crossing',
               'Ref N-ter GLN', 'Ref C-ter GLN', 'N-ter GLN', 'C-ter GLN',
               'Ref N-ter Linking Number', 'Ref C-ter Linking Number', 'N-ter Linking Number', 'C-ter Linking Number']
index_list = []
for state_id, rep_chg_ent in enumerate(rep_chg_ent_list):
    index_list.append(state_id+1)
    [traj_idx, frame_idx] = rep_chg_ent[:2]
    nc = tuple(rep_chg_ent[2:])
    keyword = ent_cluster_idx_map[state_id][0]
    chg_ent_fingerprint = chg_ent_fingerprint_list[traj_idx][frame_idx][nc]
    cross = []
    for i in range(len(chg_ent_fingerprint['crossing_resid'])):
        cross.append([])
        for j in range(len(chg_ent_fingerprint['crossing_resid'][i])):
            cross[-1].append(chg_ent_fingerprint['crossing_pattern'][i][j]+str(chg_ent_fingerprint['crossing_resid'][i][j]))
    ref_cross = []
    for i in range(len(chg_ent_fingerprint['ref_crossing_resid'])):
        ref_cross.append([])
        for j in range(len(chg_ent_fingerprint['ref_crossing_resid'][i])):
            ref_cross[-1].append(chg_ent_fingerprint['ref_crossing_pattern'][i][j]+str(chg_ent_fingerprint['ref_crossing_resid'][i][j]))
    GLN = chg_ent_fingerprint['GLN']
    ref_GLN = chg_ent_fingerprint['ref_GLN']
    LN = chg_ent_fingerprint['linking_number']
    ref_LN = chg_ent_fingerprint['ref_linking_number']

    data_0 = [keyword, idx2trajfile[traj_idx], frame_idx, nc,
              ref_cross[0], ref_cross[1], cross[0], cross[1],
              ref_GLN[0], ref_GLN[1], GLN[0], GLN[1],
              ref_LN[0], ref_LN[1], LN[0], LN[1]]
    data.append(data_0)
df = pd.DataFrame(data, columns=column_list, index=index_list)
df.to_csv('rep_chg_ent_%s.csv'%('_'.join(classify_key)), index_label='State ID')

# plot entanglement distribution
n_cluster = len(ent_cluster_idx_map)
fig = plt.figure(figsize=(np.max([6, 0.3*n_cluster]),5))
ax = fig.add_subplot(1,1,1)
window_width = 0.8
for state_id, [key, cluster_idx] in enumerate(ent_cluster_idx_map):
    cluster = ent_cluster_data[key][cluster_idx]
    nc_list = [c[2:4] for c in cluster]
    sort_index = [i for i, x in sorted(enumerate(nc_list), key=lambda x: (x[1][1]-x[1][0], x[1][0]))]
    sort_index = np.array(sort_index)
    if len(sort_index) <= max_plot_samples:
        plot_idx = np.arange(0, len(sort_index), 1, dtype=int)
    else:
        plot_idx = np.linspace(0, len(sort_index)-1, max_plot_samples, dtype=int)
    for idx, ci in enumerate(sort_index[plot_idx]):
        c = cluster[ci]
        nc = c[2:4]
        traj_idx = c[0]
        frame_idx = c[1]
        fingerprint = chg_ent_fingerprint_list[traj_idx][frame_idx][tuple(nc)]
        crossings = fingerprint['crossing_resid']
        ref_crossings = fingerprint['ref_crossing_resid']
        # plot loop
        x = state_id+1-window_width/2 + (idx+1)*window_width/(len(plot_idx)+1)
        ax.plot([x,x], nc, '-', color='tomato', linewidth=0.5, alpha=0.4)
        # plot crossings
        x = state_id+1-window_width/2 + (idx+1)*window_width/(len(plot_idx)+1)
        for ccr in ref_crossings:
            for cc in ccr:
                ax.plot([x, x], [cc-0.5, cc+0.5], '-', color='green', linewidth=0.5, alpha=0.4)
        for ccr in crossings:
            for cc in ccr:
                ax.plot([x, x], [cc-0.5, cc+0.5], '-', color='blue', linewidth=0.5, alpha=0.4)
ax.set_xticks(np.arange(1,n_cluster+1,1), np.arange(1,n_cluster+1,1))
ax.set_xlim([0, n_cluster+1])
ax.set_xlabel('Cluster')
ax.set_ylabel('Residue index')
fig.savefig('chg_ent_%s_distribution.pdf'%('_'.join(classify_key)), bbox_inches='tight')
del fig

print('Clustering structures with unique combinations of changes of entanglements...')
# Assign entanglement clusters (list of ent_cluster_idx) in discrete trajectories
for key, ent_clusters in ent_cluster_data.items():
    for i, ent_cluster in enumerate(ent_clusters):
        cluster_id = ent_cluster_idx_map.index([key, i])
        for chg_ent_keyword in ent_cluster:
            traj_idx, frame = chg_ent_keyword[:2]
            dtrajs[traj_idx][idx2frame[traj_idx].index(frame)].append(cluster_id)

# Strip same cluster ids in each frame
chg_ent_structure_keyword_list = []
for dtraj in dtrajs:
    for i, cluster_id_list in enumerate(dtraj):
        dtraj[i] = sorted(list(set(cluster_id_list)))
        if str(dtraj[i]) not in chg_ent_structure_keyword_list:
            chg_ent_structure_keyword_list.append(str(dtraj[i]))
chg_ent_structure_keyword_list = sorted(chg_ent_structure_keyword_list)

# cluster trajectory frames with different combinations of changes in entanglement
chg_ent_structure_cluster_data = {chg_ent_structure_keyword: [] for chg_ent_structure_keyword in chg_ent_structure_keyword_list}
for traj_idx, dtraj in enumerate(dtrajs):
    for frame_idx, cluster_id_list in enumerate(dtraj):
        chg_ent_structure_cluster_data[str(cluster_id_list)].append([traj_idx, idx2frame[traj_idx][frame_idx]])
Num_struct_list = [len(chg_ent_structure_cluster_data[keyword]) for keyword in chg_ent_structure_keyword_list]
sort_idx = np.argsort(-np.array(Num_struct_list, dtype=int))
sorted_chg_ent_structure_keyword_list = [chg_ent_structure_keyword_list[idx] for idx in sort_idx]
sorted_Num_struct_list = [Num_struct_list[idx] for idx in sort_idx]

print('Find representative combinations of changes of entanglements in structures...')
# Find representative changes of entanglement (minimal loop) in each frame
rep_chg_ent_dtrajs = []
for traj_idx, dtraj in enumerate(dtrajs):
    rep_chg_ent_dtrajs.append([])
    for frame_idx, cluster_id_list in enumerate(dtraj):
        frame_idx_0 = idx2frame[traj_idx][frame_idx]
        rep_chg_ent_dtrajs[-1].append({})
        for cluster_id in cluster_id_list:
            [ent_keyword, idx] = ent_cluster_idx_map[cluster_id]
            nc_list = []
            for element in ent_cluster_data[ent_keyword][idx]:
                if traj_idx == element[0] and frame_idx_0 == element[1]:
                    nc_list.append(tuple(element[2:]))
            rep_nc = nc_list[0]
            for nc in nc_list:
                if nc[1]-nc[0] < rep_nc[1]-rep_nc[0]:
                    rep_nc = nc
            rep_chg_ent_dtrajs[-1][-1][cluster_id] = chg_ent_fingerprint_list[traj_idx][frame_idx_0][rep_nc]
    
# Find representative structures (max Q) for each combination
rep_struct_data = {}
for keyword in sorted_chg_ent_structure_keyword_list:
    Q_0 = 0
    rep_struct_data[keyword] = chg_ent_structure_cluster_data[keyword][0]
    for [traj_idx, frame_idx] in chg_ent_structure_cluster_data[keyword]:
        Q = Q_list[traj_idx][frame_idx]
        if Q > Q_0:
            Q_0 = Q
            rep_struct_data[keyword] = [traj_idx, frame_idx]

# Create dataframe and save data
data = []
# column_list = ['Rep_chg_ents', 'Num of structures', 'Probability',
#                'Rep trajectory', 'Rep frame', 'Max Q', 'Median Q', 'Max RMSD (nm)', 'Median RMSD (nm)']
column_list = ['Rep_chg_ents', 'Num of structures', 'Probability',
               'Rep trajectory', 'Rep frame', 'Max Q', 'Median Q']
index_list = []
tot_num_frames = 0
for traj_idx, dtraj in enumerate(dtrajs):
    tot_num_frames += len(dtraj)
for state_id, keyword in enumerate(sorted_chg_ent_structure_keyword_list):
    index_list.append(state_id+1)
    Rep_chg_ents = str(list(np.array(eval(keyword))+1))
    Num = len(chg_ent_structure_cluster_data[keyword])
    Prob = Num / tot_num_frames
    Q_0_list = [Q_list[cd[0]][cd[1]] for cd in chg_ent_structure_cluster_data[keyword]]
    max_Q = np.max(Q_0_list)
    median_Q = np.median(Q_0_list)
    # idx_list = [combined_traj_idx_list.index(cd) for cd in chg_ent_structure_cluster_data[keyword]]
    # RMSD_list = squareform(RMSD_array[np.ix_(idx_list,idx_list)], checks=False)
    # max_RMSD = 0
    # median_RMSD = 0
    # if len(RMSD_list) > 0:
        # max_RMSD = np.max(RMSD_list)
        # median_RMSD = np.median(RMSD_list)
    # data_0 = [Rep_chg_ents, Num, Prob, idx2trajfile[rep_struct_data[keyword][0]], rep_struct_data[keyword][1], max_Q, median_Q, max_RMSD, median_RMSD]
    data_0 = [Rep_chg_ents, Num, Prob, idx2trajfile[rep_struct_data[keyword][0]], rep_struct_data[keyword][1], max_Q, median_Q]
    data.append(data_0)
df = pd.DataFrame(data, columns=column_list, index=index_list)
df.to_csv('chg_ent_struct_%s.csv'%('_'.join(classify_key)), index_label='State ID')

# Save data
np.savez('cluster_data_%s.npz'%('_'.join(classify_key)),
         chg_ent_fingerprint_list=chg_ent_fingerprint_list,
         Q_list=Q_list,
         chg_ent_keyword_dict=chg_ent_keyword_dict,
         chg_ent_keyword_list=chg_ent_keyword_list,
         idx2trajfile=idx2trajfile,
         idx2frame=idx2frame,
         # RMSD_array=RMSD_array,
         ent_cluster_data=ent_cluster_data,
         ent_cluster_tree=ent_cluster_tree,
         rep_chg_ent_list=rep_chg_ent_list,
         dtrajs=dtrajs,
         rep_chg_ent_dtrajs=rep_chg_ent_dtrajs,
         sorted_chg_ent_structure_keyword_list=sorted_chg_ent_structure_keyword_list,
         chg_ent_structure_cluster_data=chg_ent_structure_cluster_data,
         rep_struct_data=rep_struct_data)

if if_viz:
    # Generate visualiztion for representative changes of entanglement in each cluster
    print('Generate visualization for representative changes of entanglement...')
    os.system('mkdir viz_rep_chg_ent_%s'%('_'.join(classify_key)))
    for state_id, rep_chg_ent in enumerate(rep_chg_ent_list):
        state_cor = mdt.load(idx2trajfile[rep_chg_ent[0]], top=psf_file)[rep_chg_ent[1]].center_coordinates().xyz*10
        nc = (rep_chg_ent[2], rep_chg_ent[3])
        chg_ent_fingerprint = chg_ent_fingerprint_list[rep_chg_ent[0]][rep_chg_ent[1]][nc]
        rep_ent_dict = {tuple(chg_ent_fingerprint['code']): [chg_ent_fingerprint]}
        os.chdir('viz_rep_chg_ent_%s'%('_'.join(classify_key)))
        gen_state_visualizion(state_id+1, psf_file, state_cor, native_AA_pdb, rep_ent_dict, if_backmap=if_backmap, pulchra_only=pulchra_only)
        os.chdir('../')
    
    # Generate visualiztion for representative changes of entanglements in each structural cluster
    print('Generate visualization for unique entangled structures...')
    if top_struct >= 1:
        viz_dir = 'viz_chg_ent_struct_%s_%d'%('_'.join(classify_key), top_struct)
    else:
        viz_dir = 'viz_chg_ent_struct_%s_%.4f'%('_'.join(classify_key), top_struct)
    os.system('mkdir %s'%viz_dir)
    for state_id, keyword in enumerate(sorted_chg_ent_structure_keyword_list):
        if top_struct >= 1 and state_id >= top_struct:
            break
        elif top_struct < 1 and sorted_Num_struct_list[state_id]/tot_num_frames < top_struct:
            break
        [traj_idx, frame_idx] = rep_struct_data[keyword]
        state_cor = mdt.load(idx2trajfile[traj_idx], top=psf_file)[frame_idx].center_coordinates().xyz*10
        frame_idx_0 = idx2frame[traj_idx].index(frame_idx)
        rep_chg_ent_dict = rep_chg_ent_dtrajs[traj_idx][frame_idx_0]
        rep_ent_dict = {tuple(v['code']): [] for k, v in rep_chg_ent_dict.items()}
        for k, v in rep_chg_ent_dict.items():
            rep_ent_dict[tuple(v['code'])].append(v)
        os.chdir(viz_dir)
        gen_state_visualizion(state_id+1, psf_file, state_cor, native_AA_pdb, rep_ent_dict, if_backmap=if_backmap, pulchra_only=pulchra_only)
        os.chdir('../')
print('All done.')
