#!/usr/bin/env python3
import sys, getopt, math, os, multiprocessing, time, traceback
import numpy as np
import parmed as pmd
import mdtraj as mdt

################################# Arguments ###################################
# Default parameters
n_traj = 100
rep_per_traj = 10
mutant_type_list = ['fast', 'slow']
act_mask = ''
prefix_dir = ''

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
        if line.startswith('rep_per_traj'):
            words = line.split('=')
            rep_per_traj = int(words[1].strip())
            continue
        if line.startswith('mutant_type_list'):
            words = line.split('=')
            mutant_type_list = words[1].strip().split()
            continue
        if line.startswith('act_mask'):
            words = line.split('=')
            act_mask = words[1].strip()
            continue
        if line.startswith('prefix_dir'):
            words = line.split('=')
            prefix_dir = words[1].strip()
            continue
finally:
     file_object.close()

################################# Functions ###################################    
def get_Qbb_act():
    global n_traj, rep_per_traj, co_dir, max_length, cor_file, po_dir, act_mask, sec_def, dom_def, mutant_type
    print(mutant_type+' QBB_act:')
    qbb_list = [];
    for i in range(n_traj):
        traj = mdt.load(co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'_ejection.dcd', top=co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'.psf')
        traj += mdt.load(co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'_dissociation.dcd', top=co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'.psf')
        sel = traj.topology.select('resid 0 to '+str(max_length-1))
        traj = traj.atom_slice(sel)
        os.system('mkdir tmp')
        traj.save('tmp/tmp.dcd', force_overwrite=True)
        os.system("calc_native_contact_fraction.pl -i "+cor_file+" -d "+dom_def+" -s '"+sec_def+"' -t tmp/tmp.dcd -m 1 -o ./tmp/ -k '"+act_mask+"' > /dev/null")
        f = open('tmp/qbb_tmp.dat')
        C = f.readlines()
        f.close()
        os.system('rm -rf tmp/')
        C = [C[k].strip().split() for k in range(1,len(C))]
        C = np.array(C, dtype=np.float32)
        Q_ts_0 = C[:,-1].reshape((len(C[:,-1]),1))
        for j in range(rep_per_traj):
            fo = open(po_dir+'analysis/active_pocket_qbb/qbb_traj_'+str(i+1)+'_'+str(j+1)+'.dat', 'r')
            C = fo.readlines()
            fo.close()
            C = [C[k].strip().split() for k in range(1,len(C))]
            C = np.array(C, dtype=np.float32)
            c_start = C.shape[1]-1
            c_end = C.shape[1]
            Q_ts = np.vstack((Q_ts_0, C[:,c_start:c_end].reshape((C.shape[0],c_end-c_start))))
            qbb_list.append(Q_ts)
            print('  %s Traj %d_%d Done'%(mutant_type, i+1, j+1))
    np.save('%s_QBB_act.npy'%mutant_type, qbb_list)
    
def get_entanglement():
    global n_traj, rep_per_traj, co_dir, max_length, cor_file, po_dir, mutant_type
    print(mutant_type+' G:')
    G_number_list = [];
    for i in range(n_traj):
        traj = mdt.load(co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'_ejection.dcd', top=co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'.psf')
        traj += mdt.load(co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'_dissociation.dcd', top=co_dir+'traj/'+str(i+1)+'/rnc_l'+str(max_length)+'.psf')
        sel = traj.topology.select('resid 0 to '+str(max_length-1))
        traj = traj.atom_slice(sel)
        os.system('mkdir tmp')
        traj.save('tmp/tmp.dcd', force_overwrite=True)
        os.system("calc_entanglement_number.pl -i "+cor_file+" -t tmp/tmp.dcd -o ./tmp/ > /dev/null")
        f = open('tmp/G_tmp.dat')
        C = f.readlines()
        f.close()
        os.system('rm -rf tmp/')
        C0 = [C[k].strip().split() for k in range(8,len(C))]
        for CC in C0:
            gmax = np.array(CC[0].split(','), dtype=np.int)
            CC[0] = np.abs(gmax).max()
        C0 = np.array(C0, dtype=np.float32)
        for j in range(rep_per_traj):
            G_list_1 = []
            f = open(po_dir+'analysis/entanglement/G_traj_'+str(i+1)+'_'+str(j+1)+'.dat')
            C = f.readlines()
            f.close()
            C1 = [C[k].strip().split() for k in range(8,len(C))]
            for CC in C1:
                gmax = np.array(CC[0].split(','), dtype=np.int)
                CC[0] = np.abs(gmax).max()
            C1 = np.array(C1, dtype=np.float32)
            G = np.vstack((C0, C1))
            G_number_list.append(G)
            print('  %s Traj %d_%d Done %d'%(mutant_type, i+1, j+1, G_number_list[-1].shape[0]))   
    np.save('%s_Entanglement.npy'%mutant_type, G_number_list)
    
################################## MAIN #######################################
for mutant_type in mutant_type_list:
    co_dir = prefix_dir+'/continuous_synthesis/'+mutant_type+('/1-%d/'%(n_traj))
    po_dir = prefix_dir+'/post_translation/'+mutant_type+('/1-%d/'%(n_traj))
    co_traj_dir = co_dir+'analysis/qbb_full_vs_T/'
    po_traj_dir = po_dir+'analysis/qbb/'

    cor_file = os.popen('ls %s/setup/*_ca.cor'%po_dir).readlines()[0].strip()
    psf_file = os.popen('ls %s/setup/*_ca.psf'%po_dir).readlines()[0].strip()

    sec_def = os.popen('ls %s/setup/secondary_*'%po_dir).readlines()[0].strip()
    dom_def = os.popen('ls %s/setup/domain_*'%po_dir).readlines()[0].strip()

    structure = pmd.load_file(psf_file)
    max_length = len(structure.residues)

    # Get Qbb trajectory
    get_Qbb_act()

    # Get Entanglement trajectory
    get_entanglement()
