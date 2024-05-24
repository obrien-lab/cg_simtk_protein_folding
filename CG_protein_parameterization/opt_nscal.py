#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Yang Jiang @ PSU
"""

import os, time, traceback, io, sys, getopt, multiprocessing, random
import parmed as pmd
import numpy as np

usage = '''
  Usage: python opt_nscal.py 
                --input | -i <input.pdb> for CG model creation
                --domain | -d <domain.dat> for domain defination
                [--temp | -t] <Temperature> in Kelvin
                [--ppn | -n] <number of CPUs> for each trajectory. Default 1.
                        if 0, use a single GPU for each trajectory.
                [--tpn | -p] <total number of CPUs>. Default 10.
                [--ntraj | -j] <number of trajectories>. Default 10.
                [--restart | -r] <0 or 1> restart optimization. Default 0, not restart. 
                [--nscal | -s] <nscal_level.dat> for nscal levels. Default values were 
                               obtained from a training set of 18 small single-domain 
                               proteins.
                [--casm | -c] <0 or 1> CG model type. Default 0, C-alpha model.
                              1, C-alpha side chain model.
                [--help | -h]
  
  Example domain.dat:
  1:96 302:350 a #Domain 1 is from resid 1 to 96 and 302 to 350, and in
                 #alpha-helix class
  97:155 b #Domain 2 is from resid 97 to 155 and in beta-sheet class
  156:301 c #Domain 3 is from resid 156 to 301 and in alpha-beta class
  
  Example nscal_level.dat:
  a 1.1954 1.4704 1.7453 2.0322 2.5044 1.7453
  b 1.4732 1.8120 2.1508 2.5044 2.5044 2.1508
  c 1.1556 1.4213 1.6871 1.9644 2.5044 1.6871
  i 1.2747 1.5679 1.8611 2.1670 2.5044 1.8611
'''

################# Arguments #################
if len(sys.argv) == 1:
    print(usage)
    sys.exit()

input_pdb = ""
dom_def_file = ""
temperature = 310 # Kelvin
ppn = 1
tpn = 10
ntraj = 10
restart = 0
nscal_level_file = ""
casm = 0 # use C-alpha model

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:d:t:n:p:j:r:s:c:", ["help", "input=", "domain=", "temp=", "ppn=", "tpn=", "ntraj=", "restart=", "nscal=", "casm="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt in ('-h', "--help"):
        print(usage)
        sys.exit()
    elif opt in ("-i", "--input"):
        input_pdb = arg
    elif opt in ("-d", "--domain"):
        dom_def_file = arg
    elif opt in ("-t", "--temp"):
        temperature = float(arg)
    elif opt in ("-n", "--ppn"):
        ppn = int(arg)
    elif opt in ("-p", "--tpn"):
        tpn = int(arg)
    elif opt in ("-j", "--ntraj"):
        ntraj = int(arg)
    elif opt in ("-r", "--restart"):
        restart = int(arg)
    elif opt in ("-s", "--nscal"):
        nscal_level_file = arg
    elif opt in ("-c", "--casm"):
        casm = int(arg)

has_error = False
if input_pdb == "":
    print('Error: You must specify the path to the all-atom pdb file for the CG model creation.')
    has_error = True
if dom_def_file == "":
    print('Error: You must specify the path to the domain definition file.')
    has_error = True
if ppn < 0:
    print('Error: ppn < 0.')
if tpn <=0:
    print('Error: tpn <= 0.')
if ntraj <=0:
    print('Error: ntraj <= 0.')
if ppn > tpn:
    print('Error: ppn > tpn.')
if restart != 0 and restart != 1:
    print('Error: restart must be either 0 or 1.')
if casm != 0 and casm != 1:
    print('Error: casm must be either 0 or 1.')
    
if has_error:
    print(usage)
    sys.exit()

if ppn == 0:
    use_gpu = 1
    ppn = 1
else:
    use_gpu = 0

nproc = int(tpn/ppn) # Total number of processors
nscal_set = {"a": [1.1954, 1.4704, 1.7453, 2.0322, 2.5044, 1.7453],
             "b": [1.4732, 1.8120, 2.1508, 2.5044, 2.5044, 2.1508],
             "c": [1.1556, 1.4213, 1.6871, 1.9644, 2.5044, 1.6871],
             "i": [1.2747, 1.5679, 1.8611, 2.1670, 2.5044, 1.8611]}
if nscal_level_file != "":
    nscal_set = parse_nscal_levels(nscal_level_file)
sim_step = 66666667 #for 1000 ns
Q_threshold = 0.6688
frame_threshold = 0.98
sleep_time = 10 # s

if restart == 1 and not os.path.exists('opt_nscal.log'):
    restart = 0

################### Functions #########################
def parse_nscal_levels(nscal_level_file):
    nscal_set = {"a": [],
                 "b": [],
                 "c": [],
                 "i": []}
    f = open(nscal_level_file)
    lines = f.readlines()
    f.close()
    for line in lines:
        words = line.strip().split()
        if words[0] in list(nscal_set.keys()):
            nscal_set[words[0]] = [float(n) for n in words[1:]]
    
    for (k, v) in nscal_set.items():
        if len(v) == 0:
            print('Error: No nscal levels found for class "%s" in the file %s.'%(k, nscal_level_file))
            sys.exit()
        elif len(v) != len(nscal_set["a"]):
            print("Error: Inconsistent number of nscal levels found in the file %s."%(nscal_level_file))
            sys.exit()
    
    return nscal_set

def parse_domain(dom_def_file):
    domain = []
    domian_class = {"a": "Alpha-helix",
                    "b": "Beta-sheet",
                    "c": "Alpha-Beta"}
    f = open(dom_def_file)
    lines = f.readlines()
    f.close()
    
    n = 0 # domain index
    start_min = np.inf
    for line in lines:
        line = line.strip()
        if line == "" or line.startswith("#"):
            continue
        data = line.split("#")[0].strip().split()
        res_range = data[:-1]
        sec_class = data[-1]
        out = "   Domain %d: "%(n+1)
        if sec_class not in list(domian_class.keys()):
            print("Error: Forgot to specify class as a, b or c in domain definition?")
            sys.exit()
        dom = {"range": [],
               "class": sec_class}
        for d in res_range:
            [start, end] = d.split(':')
            out += "%s ~ %s; "%(start, end)
            if int(start) < start_min:
                start_min = int(start)
            dom["range"].append([int(start), int(end)])
        print(out + "Class: %s"%domian_class[sec_class])
        domain.append(dom)
        n += 1
    
    offset = start_min - 1
    for d in domain:
        for r in d["range"]:
            r[0] -= offset
            r[1] -= offset
    
    ndomain = len(domain)
    for i in range(ndomain-1):
        for j in range(i+1, ndomain):
            interface = {"range": [i, j],
                         "class": "i"}
            domain.append(interface)
    return domain

def clean_pdb(input_pdb):
    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP'];
    print("-> Cleaning PDB file %s"%input_pdb)
    name = input_pdb.split('/')[-1].split('.pdb')[0]
    struct = pmd.load_file(input_pdb)
    sel_idx = np.zeros(len(struct.atoms))
    for idx, res in enumerate(struct.residues):
        res.number = idx+1
        if res.name in AA_name_list:
            for atm in res.atoms:
                sel_idx[atm.idx] = 1
    struct[sel_idx].save(name+'_clean.pdb', overwrite=True, altlocs="occupancy")
    print("   PDB file cleaned")
    return name+'_clean.pdb'
    
def creat_CG_model(pdb_file, casm, domain):
    if casm == 0:
        print("-> Creating Ca model")
        potential_name = "bt"
    else:
        print("-> Creating Ca-SCM model")
        potential_name = "mj"
    os.system('mkdir create_model')
    os.chdir('create_model')
    
    ndomain = 0
    for d in domain:
        if d['class'] != 'i':
            ndomain += 1
    dom_idx = 0
    dom_idx_map = []
    nscal_array = np.zeros((ndomain,ndomain))
    f = open('domain_def.dat', 'w')
    for d in domain:
        if d['class'] == 'i':
            nscal_array[d['range'][0], d['range'][1]] = d['nscal']
            nscal_array[d['range'][1], d['range'][0]] = d['nscal']
        else:
            nscal_array[dom_idx, dom_idx] = d['nscal']
            for a in d['range']:
                f.write('domain = %d-%d\n'%(a[0], a[1]))
                f.write('scale factor = %f\n'%d['nscal'])
                dom_idx_map.append(dom_idx)
            dom_idx += 1
    ndomain = len(dom_idx_map)
    for i in range(ndomain-1):
        for j in range(i+1, ndomain):
            f.write('scale factor = %f\n'%nscal_array[dom_idx_map[i],dom_idx_map[j]])
    f.close()
    
    f = open('go_model.cntrl', 'w')
    f.write('''pdbfile = ../'''+pdb_file+'''
casm = '''+str(casm)+'''
potential_name = '''+potential_name+'''
domain_file = domain_def.dat
''')
    f.close()
    os.system("create_cg_protein_model.py -f go_model.cntrl > go_model.log 2>&1")
    
    name = pdb_file.split('/')[-1].split('.pdb')[0]
    if casm == 0:
        prefix = name+'_ca'
    else:
        prefix = name+'_ca-cb'
    prm_name = name + '_nscal1_fnn1_go_' + potential_name + '.prm'
    
    if os.path.exists(prefix+'.psf'):
        os.system('cp *.psf ../')
        os.system('cp *.cor ../')
        os.system('cp *.top ../')
        os.system('cp *.prm ../')
        os.chdir('../')
    else:
        print("Error: failed to create CG model from %s\n"%pdb_file)
        sys.exit()
    return (prefix, prm_name)
    
def get_secondary_structure(pdb):
    print("-> Getting secondary structure information")
    
    screen_out = os.popen('stride '+pdb).readlines()
    if not screen_out[0].startswith('REM  -------------------- Secondary structure summary -------------------  ~~~~'):
        print(''.join(screen_out))
        sys.exit()

    pdb_struct = pmd.load_file(pdb)
    sec_ele_list = []
    for line in screen_out:
        line = line.strip()
        if line.startswith('LOC '):
            sec_name = line[5:17].strip()
            if 'Helix' in sec_name or 'Strand' in sec_name:
                chainid = line[28]
                start_resnum = int(line[21:27].strip())
                end_resnum = int(line[39:45].strip())
                length = end_resnum - start_resnum + 1
                if length >= 4:
                    start_resid = np.nan
                    end_resid = np.nan
                    for res in pdb_struct.residues:
                        if res.chain == chainid and res.number == start_resnum:
                            start_resid = res.idx
                        elif res.chain == chainid and res.number == end_resnum:
                            end_resid = res.idx
                        if not np.isnan(start_resid) and not np.isnan(end_resid):
                            break
                    if np.isnan(start_resid) or np.isnan(end_resid):
                        print('Error: Cannot find residue %d or %d in chain %s in %s.'%(start_resnum, end_resnum, chainid, pdb))
                        sys.exit()
                    sec_ele_list.append([start_resid+1, end_resid+1])
    sec_ele_list.sort(key=lambda x: x[0])
    
    f = open('setup/secondary_struc_defs.txt', 'w')
    for idx, sec_ele in enumerate(sec_ele_list):
        f.write('%d %d %d\n'%(idx+1, sec_ele[0], sec_ele[1]))
    f.close()
    print("   Done.")

def convert_time(t):
    if np.isnan(t):
        return str(t)
    h = int(t / 3600)
    m = int((t - h * 3600) / 60)
    s = int(t - h * 3600 - m * 60)
    return "%d:%d:%d"%(h,m,s)    

def run_simulation(idx, prefix, prm_name, rand):
    if use_gpu == 1:
        gpu = int(multiprocessing.current_process().name.split('-')[-1])-1-worker_idx
    else:
        gpu = -1
    os.system('post_trans_single_run_v2.py setup/%s.psf setup/%s.cor setup/%s %f %d %d %d %d ../setup/secondary_struc_defs.txt 1.1 setup/%s.cor %d'%(prefix, prefix, prm_name, temperature, ppn, idx+1, rand, sim_step, prefix, gpu))
    os.system('calc_native_contact_fraction.pl -i setup/%s.cor -d ../%s -s ../setup/secondary_struc_defs.txt -t %d.dcd -r 1'%(prefix, dom_def_file, idx+1))

########################### MAIN #########################################
os.system('mkdir setup')
domain = parse_domain(dom_def_file)
nscal_index_list = [0 for d in domain]
start_nscal_level = 0
n_calc_Q = -1
if restart == 0:
    fo = open('opt_nscal.log', 'w')
    fo.write('## Start at %s\n'%time.asctime(time.localtime(time.time())))
    fo.close()
else:
    f = open('opt_nscal.log')
    lines = f.readlines()
    f.close()
    tag_nscal = False
    dom_idx = 0
    for line in lines:
        if line.startswith('## Round '):
            start_nscal_level = int(line.strip().split(':')[0].strip().split()[-1])-1
            dom_idx = 0
        elif line.startswith('## Final nscal:'):
            print('Finished. No need to restart.')
            sys.exit()
        elif line.startswith('-> Set nscal as:'):
            tag_nscal = True
        elif line.startswith('-> Running simulations...'):
            tag_nscal = False
        elif line.startswith('-> Probability of domain stability:'):
            n_calc_Q += 1
        elif tag_nscal:
            nscal = line.strip().split('=')[-1].strip()
            old_nscal = "%.4f"%(nscal_set[domain[dom_idx]['class']][nscal_index_list[dom_idx]])
            if old_nscal != nscal:
                nscal_index_list[dom_idx] += 1
            new_nscal = "%.4f"%(nscal_set[domain[dom_idx]['class']][nscal_index_list[dom_idx]])
            if new_nscal != nscal:
                print('Error: nscal mismatch at round #%d, domain #%d level #%d'%(start_nscal_level+1, dom_idx+1, nscal_index_list[dom_idx]+1))
                sys.exit()
            dom_idx += 1
    f = open('opt_nscal.log', 'w')
    for idx, line in enumerate(lines):
        if line.startswith('## Round %d:'%(start_nscal_level+1)):
            break
        f.write(line)
    if n_calc_Q == start_nscal_level:
        for line in lines[idx+1:]:
            if line.startswith('-> Probability of domain stability:'):
                break
            f.write(line)
    f.close()

clean_pdb = clean_pdb(input_pdb)
get_secondary_structure(clean_pdb)

worker_idx = 0
for iteration in range(start_nscal_level, len(nscal_set['a'])):
    os.system('mkdir round_%d'%(iteration+1))
    os.chdir('round_%d'%(iteration+1))
    
    if n_calc_Q < start_nscal_level or iteration > start_nscal_level:
        os.system('mkdir setup')
        fo = open('../opt_nscal.log', 'a')
        fo.write("## Round %d:\n"%(iteration+1))
        fo.write("-> Set nscal as:\n")
        for dom_idx, dom in enumerate(domain):
            dom['nscal'] = nscal_set[dom['class']][nscal_index_list[dom_idx]]
            if dom['class'] == 'i':
                fo.write('   Interface %d|%d: nscal = %.4f\n'%(dom['range'][0]+1, dom['range'][1]+1, dom['nscal']))
            else:
                fo.write('   Domain %d: nscal = %.4f\n'%(dom_idx+1, dom['nscal']))
        fo.close()
        
        os.chdir('setup')
        (prefix, prm_name) = creat_CG_model("../../%s"%clean_pdb, casm, domain)
        os.system('parse_cg_prm.py -p %s -t %s'%(prm_name, prefix+'.top'))
        prm_name = prm_name.split('.prm')[0]+'.xml'
        os.chdir('../')
        
        fo = open('../opt_nscal.log', 'a')
        fo.write("-> Running simulations...\n")
        fo.close()
        
        start_time = [time.time() for i in range(ntraj)]
        start_step = [0 for i in range(ntraj)]
        for i in range(ntraj):
            if os.path.exists('%d.out'%(i+1)):
                info = os.popen('tail -n 1 %d.out'%(i+1)).readlines()[0]
                if not info.strip().startswith('Time') and not info.strip() == '':
                    start_step[i] = int(info.strip().split()[1])
        
        pool = multiprocessing.Pool(nproc)
        for i in range(ntraj):
            rand = int(random.random()*1e7)
            pool.apply_async(run_simulation, (i, prefix, prm_name, rand,))

        while True:
            time.sleep(sleep_time)
            
            f = open('simulation.log', 'w')
            f.write("%-20s %-20s %-20s %-20s %-20s\n"%("#Num", "Step", "Speed(Step/s)", "Used_Time", "Rest_Time"))
            for i in range(ntraj):
                if os.path.exists('qbb_%d.dat'%(i+1)):
                    current_step = 'Done'
                    speed = '--'
                    used_time = np.nan
                elif os.path.exists('%d.out'%(i+1)):
                    info = os.popen('tail -n 1 %d.out'%(i+1)).readlines()[0]
                    if not info.strip().startswith('Time') and not info.strip() == '':
                        current_step = int(info.strip().split()[1])
                    else:
                        current_step = 0
                    used_time = time.time() - start_time[i]
                    speed = (current_step - start_step[i])/used_time
                else:
                    current_step = 0
                    used_time = time.time() - start_time[i]
                    speed = (current_step - start_step[i])/used_time
                if speed == 0 or speed == '--':
                    rest_time = np.nan
                else:
                    rest_time = (sim_step - current_step)/speed
                    speed = int(speed)
                f.write("%-20d %-20s %-20s %-20s %-20s\n"%(i+1, str(current_step), str(speed), convert_time(used_time), convert_time(rest_time)))
            f.close()
            if len(pool._cache) == 0:
                break
        pool.close()
        pool.join()
        worker_idx += nproc
    
    # Get Q
    fo = open('../opt_nscal.log', 'a')
    fo.write("-> Probability of domain stability:\n")
    fo.close()
    if_stable = [[] for d in domain]
    stable_frac = [[] for d in domain]
    for i in range(ntraj):
        f = open('qbb_%d.dat'%(i+1))
        lines = f.readlines()
        f.close()
        num_fold_list = np.zeros(len(domain))
        for line in lines[1:]:
            qbb_list = line.strip().split()[:-1]
            for j, qbb in enumerate(qbb_list):
                if float(qbb) > Q_threshold or float(qbb) == -1:
                    num_fold_list[j] += 1
        frac = num_fold_list/(len(lines)-1)
        for j, f in enumerate(frac):
            stable_frac[j].append(f)
            if f >= frame_threshold:
                if_stable[j].append(True)
            else:
                if_stable[j].append(False)
    fo = open('../opt_nscal.log', 'a')
    tag_break = True
    for dom_idx, dom in enumerate(domain):
        if dom['class'] == 'i':
            fo.write('   Interface %d|%d: '%(dom['range'][0]+1, dom['range'][1]+1))
        else:
            fo.write('   Domain %d: '%(dom_idx+1))
        for frac in stable_frac[dom_idx]:
            fo.write('%.3f '%frac)
        if np.all(if_stable[dom_idx]):
            fo.write('stable.\n')
        else:
            fo.write('instable.\n')
            # update nscal index
            nscal_index_list[dom_idx] += 1
            tag_break = False
    fo.close()
    os.chdir('../')
    if tag_break:
        break

fo = open('opt_nscal.log', 'a')
fo.write('## Final nscal:\n')
for dom_idx, dom in enumerate(domain):
    if dom['class'] == 'i':
        fo.write('   Interface %d|%d: nscal = %.4f\n'%(dom['range'][0]+1, dom['range'][1]+1, dom['nscal']))
    else:
        fo.write('   Domain %d: nscal = %.4f\n'%(dom_idx+1, dom['nscal']))
fo.close()
