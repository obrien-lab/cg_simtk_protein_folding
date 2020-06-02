#!/usr/bin/env python3
import getopt, os, time, sys, random, math, traceback, io, sys
import parmed as pmd
import numpy as np
import mdtraj as mdt

def renumber_structure(new_struct, ref_struct):
	for res in new_struct.residues:
		res.number = ref_struct.residues[res.idx].number
		res.segid = ref_struct.residues[res.idx].segid
	for atm in new_struct.atoms:
		atm.number = atm.idx + 1

####################### MAIN #######################
if len(sys.argv) != 5 and len(sys.argv) != 6:
	print("Usage: visualize_cont_synth.py [aa_pdb] [ribo_type] [traj_id] [start_id] [run_type]")
	sys.exit()

alpha = 4331293

ref_aa_pdb = sys.argv[1]
ribo_type = sys.argv[2]
traj_id = int(sys.argv[3])
start_id = int(sys.argv[4])

if len(sys.argv) == 6:
	run_type = int(sys.argv[5])
else:
	run_type = 0

traj_dir = './traj/'+str(traj_id)+'/'
prefix = 'rnc_l'

if ribo_type == 'ecoli':
	tcl_script = '/'.join(__file__.split('/')[:-1]) + '/render_ecoli_RNC.tcl'
elif ribo_type == 'yeast':
	tcl_script = '/'.join(__file__.split('/')[:-1]) + '/render_yeast_RNC.tcl'
else:
	print("Wrong ribosome type, could only be ecoli or yeast.")
	sys.exit()

chain_length_list = []
for file in os.listdir(traj_dir):
	if file.startswith(prefix) and file.endswith('_stage_3_final.cor'):
		name = os.path.splitext(file)[0]
		words = name.split('l')
		words = words[-2].split('_')
		chain_length_list.append(int(words[0]))

max_chain_length = max(chain_length_list)

tmp_dir = 'visualize_cont_synth_tmp'
os.system('mkdir '+tmp_dir)
os.system('mkdir '+tmp_dir+'/ppm')

pdb_pmd = pmd.load_file(ref_aa_pdb)

if run_type < 2:
	n_frame = 1
	for cl in range(start_id, max_chain_length+1):
		psf_file = traj_dir + prefix + str(cl) + '.psf'
		pdb_pmd[':1-'+str(cl)].write_pdb(tmp_dir+'/tmp_ref.pdb', charmm=True)
		for stage in range(1,4):
			traj_file = traj_dir + prefix + str(cl) + '_stage_' + str(stage) + '_final.cor'
			psf_pmd = pmd.load_file(psf_file)
			psf_pmd.positions = pmd.load_file(traj_file).positions
			psf_pmd[0:cl].write_pdb(tmp_dir+'/tmp_p.pdb', charmm=True)
			psf_pmd[cl:].write_pdb(tmp_dir+'/tmp_r.pdb', charmm=True)
			os.chdir(tmp_dir)
			os.system('backmap.py -i tmp_ref.pdb -c tmp_p.pdb')
			rebuilt_pdbname = 'tmp_p_rebuilt.pdb'
			pdb_1_pmd = pmd.load_file(rebuilt_pdbname)
			pdb_2_pmd = pmd.load_file('tmp_r.pdb')
			new_pdb_pmd = pdb_1_pmd+pdb_2_pmd
			renumber_structure(new_pdb_pmd, psf_pmd)
			new_pdb_pmd.write_pdb('traj_'+str(cl)+'_'+str(stage)+'.pdb', charmm=True, renumber=False)
			os.system('vmd -dispdev none -eofexit -e '+tcl_script+' -args '+str(cl)+' '+str(stage)+
				(' ppm/traj_1_%010d.tga'%n_frame)+' >/dev/null')
			if stage == 1:
				text = 'Length %d'%(cl-1)
			else:
				text = 'Length %d'%cl
			name = 'ppm/traj_1_%010d.tga'%n_frame
			os.system('convert -fill black -pointsize 60 -font helvetica -draw \'text 10,80 "'+text+'"\' '+name+' '+name)
			os.chdir('../')
			print('Done length '+str(cl)+' stage '+str(stage))
			n_frame += 1
	print('Converting gif...')
	os.system('convert -delay 10 '+tmp_dir+'/ppm/traj_1_*.tga -loop 0 traj_synthesis.gif')

if run_type != 1:
	n_frame = 1
	psf_file = traj_dir + prefix + str(max_chain_length) + '.psf'
	psf_pmd = pmd.load_file(psf_file)
	pdb_pmd[':1-'+str(max_chain_length)].write_pdb(tmp_dir+'/tmp_ref.pdb', charmm=True)
	for traj_type in ['ejection', 'dissociation']:
		traj_file = traj_dir + prefix + str(max_chain_length) + '_'+traj_type+'.dcd'
		traj = mdt.load(traj_file, top=psf_file)
		for f in range(traj.n_frames):
			frame = traj.slice(f)
			frame.atom_slice([i for i in range(max_chain_length)]).save(tmp_dir+'/tmp_p.pdb')
			frame.atom_slice([i for i in range(max_chain_length,traj.n_atoms)]).save(tmp_dir+'/tmp_r.pdb')
			os.chdir(tmp_dir)
			os.system('backmap.py -i tmp_ref.pdb -c tmp_p.pdb')
			rebuilt_pdbname = 'tmp_p_rebuilt.pdb'
			pdb_1_pmd = pmd.load_file(rebuilt_pdbname)
			pdb_2_pmd = pmd.load_file('tmp_r.pdb')
			new_pdb_pmd = pdb_1_pmd+pdb_2_pmd
			renumber_structure(new_pdb_pmd, psf_pmd)
			new_pdb_pmd.write_pdb('traj_'+str(max_chain_length)+'_4.pdb', charmm=True, renumber=False)
			os.system('vmd -dispdev none -eofexit -e '+tcl_script+' -args '+str(max_chain_length)+' 4'+
				(' ppm/traj_2_%010d.tga'%n_frame)+' >/dev/null')
			text = traj_type.capitalize() + ' %.3f ms'%(n_frame*0.015*5000*alpha/1e9)
			name = 'ppm/traj_2_%010d.tga'%n_frame
			os.system('convert -fill black -pointsize 60 -font helvetica -draw \'text 10,80 "'+text+'"\' '+name+' '+name)
			os.chdir('../')
			print('Done length '+str(max_chain_length)+' stage 4')
			n_frame += 1
	print('Converting gif...')
	os.system('convert -delay 5 '+tmp_dir+'/ppm/traj_2_*.tga -loop 0 traj_eject.gif')
if run_type == 0:
	os.system('convert -delay 0 traj_synthesis.gif traj_eject.gif -loop 0 traj_total.gif')
	os.system('rm -rf '+tmp_dir)
print('All Done')
