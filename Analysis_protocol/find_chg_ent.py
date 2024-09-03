#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback, copy
import pickle
import numpy as np
from scipy.spatial.distance import pdist, squareform
import parmed as pmd
import mdtraj as mdt
import topoly

start_time_0 = time.time()

usage = '''
  Usage: python find_chg_ent.py 
              --psf | -p <INPUT.PSF> for the protein structure
              --ref | -c <REFERENCE> coordinates file recognized by parmed
              --traj | -t <TRAJECTORY> for analysis, recognized by mdtraj
              [--begin | -b] <START STEP> to split trajectory. 
                                  Default is 0. Can have a negative value.
              [--end | -e] <END STEP> to split trajectory. 
                                Default is the end frame.
              [--stride | -s] <STRIDE> for frame selection. Default is 1.
              [--mask | -k] <AMBER MASK> for a subset of atoms to do the calculation. 
                                   Default is all atoms. Index starts from 1.
                                   Example MASK: 
                                   ':1-45' for selection from resid 1 to 45;
                                   ':1,3,5-9' for selection of resid 1, 3 
                                   and from 5 to 9.
                                   ':1-30@CA' for selection of CA atoms from
                                   resid 1 to 30.
              [--outdir | -o] <DIRECTORY> for the outputs. Default is the directory
                                          of your trajectories
              [--restart | -r] <0 or 1> 0: Do not restart calculation; 1: Restart 
                                        calculation. Default is 0.
              [--help | -h]\n
'''

################### Parameters ########################
def gen_chg_ent_code_dict():
    dict_1 = {'L': 'linking number',
              'C': 'linking chirality'}
    dict_2 = {'#': 'no change of',
              '+': 'gain of',
              '-': 'loss of',
              '~': 'switched',}
    chg_ent_code_dict = {}
    for chg_L in ['#', '+', '-']:
        for chg_C in ['#', '~']:
            code = 'L'+chg_L+'C'+chg_C
            if chg_L == '#' and chg_C == '#':
                string = 'no change'
            else:
                string = dict_2[chg_L] + ' ' + dict_1['L'] + ' & ' + dict_2[chg_C] + ' ' + dict_1['C']
            chg_ent_code_dict[code] = string
    return chg_ent_code_dict

chg_ent_code_dict = gen_chg_ent_code_dict()

################### Functions #########################
def calc_GLN_list(coor, sel, cutoff=8, thread_cutoff=4, terminal_cutoff=5):
    n_atom = coor.shape[0]
    # Generate contact matrix
    R = (coor[1:,:] + coor[:-1,:])/2
    dR = coor[1:,:] - coor[:-1,:]
    M_1 = np.repeat(R.reshape((R.shape[0], 1, R.shape[1])), R.shape[0], axis=1)
    M_2 = np.repeat(R.reshape((1, R.shape[0], R.shape[1])), R.shape[0], axis=0)
    v1 = (M_1 - M_2) / np.linalg.norm(M_1 - M_2, axis=-1, keepdims=True)**3
    v1[np.isnan(v1)] = 0
    M_1 = np.repeat(dR.reshape((dR.shape[0], 1, dR.shape[1])), dR.shape[0], axis=1)
    M_2 = np.repeat(dR.reshape((1, dR.shape[0], dR.shape[1])), dR.shape[0], axis=0)
    v2 = np.cross(M_1, M_2)
    M = np.sum(v1*v2, axis=-1)
    del M_1, M_2, v1, v2

    # Calculate GLN; [[N-ter_GLN, C_ter_GLN]]
    GLN_list = []
    for ii, r1_range in enumerate(sel):
        r1_i = r1_range[0]
        r1_j = r1_range[1]
        coor_1 = coor[r1_i, :]
        coor_2 = coor[r1_j, :]
        dist = np.sum((coor_1-coor_2)**2)**0.5
        G0 = [np.nan, np.nan]
        if dist <= cutoff:
            range_list = [[terminal_cutoff,r1_i-thread_cutoff], 
                          [r1_j+thread_cutoff, n_atom-terminal_cutoff-1]]
            for idx, r2_range in enumerate(range_list):
                r2_i = r2_range[0]
                r2_j = r2_range[1]
                if r2_j < 0:
                    r2_j = 0
                G0[idx] = np.sum(M[r1_i:r1_j, r2_i:r2_j])
        GLN_list.append([G0[0]/4/np.pi, G0[1]/4/np.pi])
    return (M, GLN_list)

def get_native_contacts(coor, ref_native_contact_list=None, min_interval=4, cutoff=8.0):
    n_atom = coor.shape[0]
    native_contact_list = []
    if ref_native_contact_list is None:
        for i in range(0, n_atom-min_interval):
            coor_1 = np.array(coor[i])
            for j in range(i+min_interval, n_atom):
                coor_2 = np.array(coor[j])
                dist = np.sum((coor_1-coor_2)**2)**0.5
                if (dist <= cutoff):
                    native_contact_list.append((i, j))
    else:
        for ref_nc in ref_native_contact_list:
            coor_1 = np.array(coor[ref_nc[0]])
            coor_2 = np.array(coor[ref_nc[1]])
            dist = np.sum((coor_1-coor_2)**2)**0.5
            if (dist <= cutoff):
                native_contact_list.append((i, j))
    return native_contact_list

def identify_entanglement_by_topoly(coor, native_contact_list, GLN_list, ref_ent_dict={}, cross_cutoff=10, thread_cutoff=4, terminal_cutoff=5, GLN_cutoff=0.5):
    adjacent_residue_cutoff = 2
    distance_map = squareform(pdist(coor, 'euclidean'))
    coor_list = coor.tolist()
    ent_native_contact_list = []
    for idx, nc in enumerate(native_contact_list):
        if np.any(np.abs(GLN_list[idx]) >= GLN_cutoff):
            ent_native_contact_list.append(nc)
        elif tuple(nc) in ref_ent_dict.keys() and np.sum(np.abs(ref_ent_dict[tuple(nc)]['linking_number'])) > 0:
            ent_native_contact_list.append(nc)
    results = topoly.lasso_type(coor_list, loop_indices=ent_native_contact_list, min_dist=[cross_cutoff, thread_cutoff, terminal_cutoff], more_info=True)
    nc_ent_fingerprint_dict = {}
    for idx, nc in enumerate(native_contact_list):
        fingerprint = {'GLN': GLN_list[idx],
                     'crossing_resid': [[], []],
                     'crossing_pattern': ['', ''],
                     'linking_number': [0, 0],
                     'native_contact': list(nc),
                     'surrounding_resid': [[],[]],}
        if nc in results.keys():
            crossingsN = results[nc]['crossingsN']
            crossingsC = results[nc]['crossingsC']
            for idx_ter, crossings in enumerate([crossingsN, crossingsC]):
                for c in crossings:
                    # Remove trivial crossings
                    if len(fingerprint['crossing_resid'][idx_ter]) > 0 and abs(int(c[1:])-fingerprint['crossing_resid'][idx_ter][-1]) <= adjacent_residue_cutoff:
                        # If same residue (or adjacent residue) pierces multiple surfaces
                        # !!! Be careful about "adjacent_residue_cutoff = 2", 2 used here may lose some true multiple-surface crossings in the result. !!!
                        # !!! One need to reduce it to 1 or even 0 when this situation happens. !!!
                        if fingerprint['crossing_pattern'][idx_ter][-1] == c[0]:
                            # Same piercing direction, skip this crossing
                            continue
                        else:
                            # Opposite piercing direction, crossings cancelled
                            del fingerprint['crossing_resid'][idx_ter][-1]
                            fingerprint['crossing_pattern'][idx_ter] = fingerprint['crossing_pattern'][idx_ter][:-1]
                    elif len(fingerprint['crossing_resid'][idx_ter]) > 0 and abs(int(c[1:])-fingerprint['crossing_resid'][idx_ter][-1]) < cross_cutoff:
                        if fingerprint['crossing_pattern'][idx_ter][-1] == c[0]:
                            # Same piercing direction, add this crossing
                            fingerprint['crossing_resid'][idx_ter].append(int(c[1:]))
                            fingerprint['crossing_pattern'][idx_ter] += c[0]
                        else:
                            # Opposite piercing direction, crossings cancelled
                            del fingerprint['crossing_resid'][idx_ter][-1]
                            fingerprint['crossing_pattern'][idx_ter] = fingerprint['crossing_pattern'][idx_ter][:-1]
                    else:
                        # add this crossing
                        fingerprint['crossing_resid'][idx_ter].append(int(c[1:]))
                        fingerprint['crossing_pattern'][idx_ter] += c[0]
                for cp in fingerprint['crossing_pattern'][idx_ter]:
                    if cp == '+':
                        fingerprint['linking_number'][idx_ter] += 1
                    elif cp == '-':
                        fingerprint['linking_number'][idx_ter] -= 1
                # Find surrounding resids
                dm = np.where(distance_map < 8.0, 1, 0)[fingerprint['crossing_resid'][idx_ter], :]
                surr_resid = list(np.where(np.sum(dm, axis=0) > 0)[0])
                fingerprint['surrounding_resid'][idx_ter] = surr_resid

        nc_ent_fingerprint_dict[nc] = fingerprint
    return results, nc_ent_fingerprint_dict

def identify_change_entanglement(ref_native_contact_list, nc_ent_fingerprint_dict, ref_nc_ent_fingerprint_dict, dGLN_cutoff=0.2):
    chg_ent_fingerprint = {}
    for idx, nc in enumerate(ref_native_contact_list):
        fingerprint = nc_ent_fingerprint_dict[nc]
        ref_fingerprint = ref_nc_ent_fingerprint_dict[nc]
        chg_ent = ['L#C#', 'L#C#']
        tag_modify = [False, False]
        if not np.any(np.isnan(fingerprint['GLN'])):
            for idx_ter in range(2):
                crossing_pattern = fingerprint['crossing_pattern'][idx_ter]
                ref_crossing_pattern = ref_fingerprint['crossing_pattern'][idx_ter]
                if '*' in crossing_pattern or '*' in ref_crossing_pattern:
                    # Skip the comparison
                    continue
                linking_number = fingerprint['linking_number'][idx_ter]
                GLN = fingerprint['GLN'][idx_ter]
                ref_linking_number = ref_fingerprint['linking_number'][idx_ter]
                ref_GLN = ref_fingerprint['GLN'][idx_ter]
                chg_pattern = ['#', '#'] # patterns of [linking number, chirality of linkage]
                if abs(linking_number) > abs(ref_linking_number):
                    if abs(GLN) - abs(ref_GLN) > dGLN_cutoff:
                        chg_pattern[0] = '+'
                    else:
                        tag_modify[idx_ter] = True
                elif abs(linking_number) < abs(ref_linking_number):
                    if abs(GLN) - abs(ref_GLN) < -dGLN_cutoff:
                        chg_pattern[0] = '-'
                    else:
                        tag_modify[idx_ter] = True
                if linking_number * ref_linking_number < 0 and GLN * ref_GLN < 0:
                    chg_pattern[1] = '~'
                chg_code = 'L%sC%s'%tuple(chg_pattern)
                chg_ent[idx_ter] = chg_code

        chg_ent_fingerprint[nc] = {'type': [chg_ent_code_dict[code] for code in chg_ent],
                                 'code': chg_ent,}
        for key in fingerprint.keys():
            chg_ent_fingerprint[nc][key] = copy.copy(fingerprint[key])
        for key in ref_fingerprint.keys():
            chg_ent_fingerprint[nc]['ref_'+key] = copy.copy(ref_fingerprint[key])
        # Modify fingerprints if needed
        for idx_ter, tag in enumerate(tag_modify):
            if tag:
                chg_ent_fingerprint[nc]['crossing_resid'][idx_ter] = chg_ent_fingerprint[nc]['ref_crossing_resid'][idx_ter]
                chg_ent_fingerprint[nc]['crossing_pattern'][idx_ter] = chg_ent_fingerprint[nc]['ref_crossing_pattern'][idx_ter]
                chg_ent_fingerprint[nc]['linking_number'][idx_ter] = chg_ent_fingerprint[nc]['ref_linking_number'][idx_ter]
                
    return chg_ent_fingerprint

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

def check_restart():
    global restart, begin, end
    num_header_lines = len(list(chg_ent_code_dict.keys()))+3
    if restart == 1:
        last_G_frame = np.nan
        last_ent_frame = np.nan
        if not os.path.exists(out_G_file) or os.path.getsize(out_G_file) == 0:
            print('Warning: No output file found. Automatically set restart = 0.')
            restart = 0
        elif not os.path.exists(out_data_file) or os.path.getsize(out_data_file) == 0:
            print('Warning: No ent data file found. Automatically set restart = 0.')
            restart = 0
        else:
            last_line_idx = np.nan
            last_G_frame = np.nan
            with open(out_G_file, 'r') as f:
                lines = f.readlines()
            if len(lines) <= num_header_lines:
                restart = 0
            elif len(lines) == num_header_lines+1 and not lines[-1].endswith('\n'):
                restart = 0
            else:
                if lines[-1].endswith('\n'):
                    last_line = lines[-1]
                    last_line_idx = len(lines)
                else:
                    last_line = lines[-2]
                    last_line_idx = len(lines)-1
                last_G_frame = int(last_line.strip().split()[0])
            
            try:
                ent_data = load_pickle(out_data_file)
            except Exception as e:
                traceback.print_exc()
                restart = 0
            else:
                ent_data_frame_list = list(ent_data.keys())
                # sort frames, for old python version that doesn't support ordered dictionary
                ent_data_frame_list = ['ref'] + sorted([e for e in ent_data_frame_list if e != 'ref'])
                if len(ent_data_frame_list) == 1:
                    restart = 0
                elif len(ent_data_frame_list) == 2 and len(list(ent_data[ent_data_frame_list[-1]].keys())) < len(list(ent_data[ent_data_frame_list[0]].keys())):
                    restart = 0
                else:
                    if len(list(ent_data[ent_data_frame_list[-1]].keys())) < len(list(ent_data[ent_data_frame_list[0]].keys())):
                        last_ent_frame = ent_data_frame_list[-2]
                    else:
                        last_ent_frame = ent_data_frame_list[-1]

            if restart == 1:
                last_frame = np.min([last_ent_frame, last_G_frame])
                begin = last_frame + stride
                print("Restart at frame %d"%begin)
                if end < begin:
                    print("Warning: end frame id %d < start frame id %d. No action took."%(end, begin))
                    sys.exit()
                with open(out_G_file, 'w') as f:
                    for line in lines[:num_header_lines]:
                        f.write(line)
                    for line in lines[num_header_lines:]:
                        frame = int(line.strip().split()[0])
                        if frame <= last_frame:
                            f.write(line)
                for frame in ent_data_frame_list:
                    if frame == 'ref':
                        save_pickle(out_data_file, 'wb', {frame: ent_data[frame]})
                    elif frame <= last_frame:
                        save_pickle(out_data_file, 'ab', {frame: ent_data[frame]})
    if restart == 0:
        if begin < 0:
            begin = n_frames + begin
        if begin < 0 or begin >= n_frames:
            print('Error: invalid start frame index')
            sys.exit()
        else:
            print("Start at frame %d"%begin)

def convert_time(seconds):
    if np.isnan(seconds):
        return "%6s:%2s:%2s"%('--', '--', '--')
    else:
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        return "%6d:%02d:%02d"%(h, m, s)

############################ MAIN #########################################
if len(sys.argv) == 1:
    print(usage)
    sys.exit()

psf_file = None
ref_cor_file = None
traj_file = None
begin = 0
end = None
stride = 1
mask = ':='
outdir = './'
restart = 0
try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "h:p:c:t:b:e:s:k:o:r:", 
                               ["help", "psf=", "ref=", "traj=", "begin=", "end=", "stride=", "mask=", "outdir=", "restart="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt in ("-p", "--psf"):
        psf_file = arg
    elif opt in ("-c", "--ref"):
        ref_cor_file = arg
    elif opt in ("-t", "--traj"):
        traj_file = arg
    elif opt in ("-b", "--begin"):
        begin = int(arg)
    elif opt in ("-e", "--end"):
        end = int(arg)
    elif opt in ("-s", "--stride"):
        stride = int(arg)
    elif opt in ("-k", "--mask"):
        mask = arg
    elif opt in ("-o", "--outdir"):
        outdir = arg
    elif opt in ("-r", "--restart"):
        restart = int(arg)

tag_err = False
if psf_file is None:
    print('Error: No psf file specified.')
    tag_err = True
if ref_cor_file is None:
    print('Error: No reference coordinates file specified.')
    tag_err = True
if traj_file is None:
    print('Error: No trajectory file specified.')
    tag_err = True
if restart != 0 and restart != 1:
    print('Error: -r must be either 0 or 1.')
    tag_err = True
    
if tag_err:
    print(usage)
    sys.exit()
    
# Load reference structure
struct = pmd.load_file(psf_file)
ref_coor = pmd.load_file(ref_cor_file).coordinates
struct.coordinates = ref_coor
# Load trajectory
traj = mdt.load(traj_file, top=psf_file).center_coordinates()
# Apply selection
sel_array = pmd.amber.AmberMask(struct, mask).Selection()
sel_array = [i for i, sel in enumerate(sel_array) if sel == 1]
if len(sel_array) == 0:
    print('Error: Found 0 atoms base on the mask "%s"'%mask)
    sys.exit()
else:
    print('%d atoms selected'%len(sel_array))
struct = struct[mask]
ref_coor = struct.coordinates
traj = traj.atom_slice(sel_array)
n_frames = len(traj.xyz)

out_G_file = outdir + '/' + traj_file.split('/')[-1].split('.')[0] + '_G.dat'
out_data_file = outdir + '/' + traj_file.split('/')[-1].split('.')[0] + '_ent.pkl'
if end is None or end > n_frames-1:
    end = n_frames-1
    print('Automatically change the end frame id to maximum frame %d'%(end))
# check restart
check_restart()

# Analyze reference structure
print('Analyzing reference structure... ', end='\r')
ref_native_contact_list = get_native_contacts(ref_coor)
ref_nnc = len(ref_native_contact_list)
ref_M, ref_GLN_list = calc_GLN_list(ref_coor, ref_native_contact_list)
if restart == 0:
    ref_topoly_results, ref_nc_ent_fingerprint_dict = identify_entanglement_by_topoly(ref_coor, ref_native_contact_list, ref_GLN_list, GLN_cutoff=0.0)
else:
    ref_nc_ent_fingerprint_dict = load_pickle(out_data_file)['ref']['ent_fingerprint']
print('Analyzing reference structure... Done.')
if restart == 0:
    with open(out_G_file, 'w') as f:
        for code_idx, ent_code in enumerate(sorted(list(chg_ent_code_dict.keys()), reverse=True)):
            ent_type = chg_ent_code_dict[ent_code]
            if ent_code == 'L#C#':
                f.write('# G%d: %s, loop formed & %s\n'%(code_idx+1, ent_code, ent_type))
            else:
                f.write('# G%d: %s, %s\n'%(code_idx+1, ent_code, ent_type))
        f.write('# G: Number of change of entanglement (G1+...+G%d) / (2 x Number of native contacts in reference structure)\n'%code_idx)
        f.write('# Number of native contact in the reference structure: %d\n'%(ref_nnc))
        f.write('%-12s'%'Frame')
        for code_idx, ent_code in enumerate(sorted(list(chg_ent_code_dict.keys()), reverse=True)):
            f.write(' %-5s'%('G%d'%code_idx))
        f.write(' %-6s\n'%'G')
    pickle_data = {'ref': {'ent_fingerprint': ref_nc_ent_fingerprint_dict,
                           'chg_ent_fingerprint': None,
                           'G_dict': None,
                           'G': None,}}
    save_pickle(out_data_file, 'wb', pickle_data)

# Loop through frames
start_time = time.time()
print('In total %d frames will be analyzed.'%(int((end-begin)/stride)+1))
for i_frame, coor in enumerate(traj.xyz[begin:end+1:stride,:,:]):
    frame_idx = begin + i_frame * stride
    # Progress
    progress = i_frame/(int((end-begin)/stride)+1)
    speed = progress/(time.time() - start_time)
    if speed == 0:
        remaining_time = np.nan
    else:
        remaining_time = (1-progress)/speed
    print('Analyzing frame %12d (%5.1f%%), %s remaining'%(frame_idx, progress*100, convert_time(remaining_time)), end='\r')
    # Analysis
    coor = coor * 10 # convert nanometer to Angstrom
    M, GLN_list = calc_GLN_list(coor, ref_native_contact_list)
    topoly_results, nc_ent_fingerprint_dict = identify_entanglement_by_topoly(coor, ref_native_contact_list, GLN_list, 
                                                                              ref_ent_dict=ref_nc_ent_fingerprint_dict, GLN_cutoff=0.5)
    chg_ent_fingerprint = identify_change_entanglement(ref_native_contact_list, nc_ent_fingerprint_dict, ref_nc_ent_fingerprint_dict)
    # Calculate G
    G_dict = {ent_code: 0 for ent_code in sorted(list(chg_ent_code_dict.keys()), reverse=True)}
    for nc, fingerprint in chg_ent_fingerprint.items():
        if not np.any(np.isnan(fingerprint['GLN'])):
            for chg_ent in fingerprint['code']:
                G_dict[chg_ent] += 1
    G = 0
    for ent_code in sorted(list(chg_ent_code_dict.keys()), reverse=True):
        if ent_code != 'L#C#':
            G += G_dict[ent_code]
    G = G / 2 / ref_nnc
    # Write data
    pickle_data = {frame_idx: {'ent_fingerprint': nc_ent_fingerprint_dict,
                               'chg_ent_fingerprint': chg_ent_fingerprint,
                               'G_dict': G_dict,
                               'G': G,}}
    with open(out_G_file, 'a') as f:
        f.write('%-12d'%(frame_idx))
        for code_idx, ent_code in enumerate(sorted(list(chg_ent_code_dict.keys()), reverse=True)):
            f.write(' %-5d'%G_dict[ent_code])
        f.write(' %-6.4f\n'%G)
    save_pickle(out_data_file, 'ab', pickle_data)
print('Analyzing frame %12d (%5.1f%%), %s remaining'%(frame_idx, 100.0, convert_time(0)))
print('All done. Elapsed time: %s'%(convert_time(time.time() - start_time_0)))
