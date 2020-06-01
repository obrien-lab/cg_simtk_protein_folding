#!/usr/bin/env python3
import getopt, os, sys, time, multiprocessing, random, math
import parmed as pmd
import numpy as np

usage = '\nUsage: python gen_ribosome_FF.py\n' \
        '       --pdblist | -p <XXX YYY> input pdb file list\n'\
        '       [-h] Print this information\n\n'


#################################### Data ####################################
protein_cg_name_list = ['SA', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH', 'SI', 'SK', 'SL',
						'SM', 'SN', 'SP', 'SQ', 'SR', 'SS', 'ST', 'SV', 'SW', 'SY']
ribo_cg_name_list = ['P', 'R', 'PU1', 'PU2', 'PY']

# from Dan
#protein_cg_R = [3.753086, 4.078775, 4.231107, 4.414068, 4.740325, 3.358232, 4.532358, 4.607911, 4.740325, 
#				4.607911, 4.607911, 4.141056, 4.480555, 4.887845, 3.861931, 4.186565, 4.359399, 5.047691, 4.809537]

# from Ed
protein_cg_R = [2.862278, 3.030648, 3.142894, 3.367386, 3.535755, 2.525540, 3.423509, 3.423509, 3.535755, 3.423509, 
				3.423509, 3.199017, 3.086771, 3.423509, 3.704125, 2.918401, 3.142894, 3.311263, 3.816371, 3.591879]

#################################### MAIN ####################################
pdb_list = ''

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hp:", ["pdblist="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-p", "--pdblist"):
        pdb_list = arg.strip().split()

if pdb_list == '' or pdb_list == []:
	print(usage)
	sys.exit()

fo = open('NC_results.dat', 'w')
fo.close()
fo = open('L24_results.dat', 'w')
fo.close()
fo = open('RP_results.dat', 'w')
fo.close()

dist_map_NC = [[['-' for i in (ribo_cg_name_list)] for j in protein_cg_name_list] for k in pdb_list]
dist_map_L24 = [[['-' for i in (ribo_cg_name_list)] for j in protein_cg_name_list] for k in pdb_list]
collision_diameter = [['-' for i in (protein_cg_name_list)] for k in pdb_list]
for i in range(len(pdb_list)):
	pdb = pdb_list[i]
	pdb_struct = pmd.load_file(pdb+'_cg/50S_cg.pdb')
	tag = 0
	for res in pdb_struct.residues:
		if res.segid == 'NC':
			tag = 1
			break
	if tag == 0:
		print('Warning: no nascent chain found in pdb %s!'%pdb)
		#sys.exit()

	############ NC -- 23S+5S ###########
	for nc_atm_idx, nc_atm_name in enumerate(protein_cg_name_list):
		for atm_1 in pdb_struct.atoms:
			if atm_1.residue.segid == 'NC' and atm_1.name == nc_atm_name:
				for rr_atm_idx, rr_atm_name in enumerate(ribo_cg_name_list):
					for atm_2 in pdb_struct.atoms:
						if (atm_2.residue.segid == '23S' or atm_2.residue.segid == '5S') and atm_2.name == rr_atm_name:
							if dist_map_NC[i][nc_atm_idx][rr_atm_idx] == '-':
								dist_map_NC[i][nc_atm_idx][rr_atm_idx] = 9999999
							dist = ((atm_1.xx - atm_2.xx)**2 + (atm_1.xy - atm_2.xy)**2 + (atm_1.xz - atm_2.xz)**2)**0.5
							if dist_map_NC[i][nc_atm_idx][rr_atm_idx] > dist:
								dist_map_NC[i][nc_atm_idx][rr_atm_idx] = dist
	fo = open('NC_results.dat', 'a')
	fo.write('PDB: %s\n'%pdb)
	fo.write('%8s '%' ')
	for rn in ribo_cg_name_list:
		fo.write('%8s '%rn)
	fo.write('\n')
	for j, pn in enumerate(protein_cg_name_list):
		fo.write('%8s '%pn)
		for k, rn in enumerate(ribo_cg_name_list):
			if dist_map_NC[i][j][k] == '-':
				fo.write('%8s '%dist_map_NC[i][j][k])
			else:
				fo.write('%8.4f '%dist_map_NC[i][j][k])
		fo.write('\n')
	fo.write('\n')
	fo.close()
	############ L24 -- 23S+5S ###########
	for nc_atm_idx, nc_atm_name in enumerate(protein_cg_name_list):
		for atm_1 in pdb_struct.atoms:
			if atm_1.residue.segid == 'L24' and atm_1.name == nc_atm_name:
				for rr_atm_idx, rr_atm_name in enumerate(ribo_cg_name_list):
					for atm_2 in pdb_struct.atoms:
						if (atm_2.residue.segid == '23S' or atm_2.residue.segid == '5S') and atm_2.name == rr_atm_name:
							if dist_map_L24[i][nc_atm_idx][rr_atm_idx] == '-':
								dist_map_L24[i][nc_atm_idx][rr_atm_idx] = 9999999
							dist = ((atm_1.xx - atm_2.xx)**2 + (atm_1.xy - atm_2.xy)**2 + (atm_1.xz - atm_2.xz)**2)**0.5
							if dist_map_L24[i][nc_atm_idx][rr_atm_idx] > dist:
								dist_map_L24[i][nc_atm_idx][rr_atm_idx] = dist
	fo = open('L24_results.dat', 'a')
	fo.write('PDB: %s\n'%pdb)
	fo.write('%8s '%' ')
	for rn in ribo_cg_name_list:
		fo.write('%8s '%rn)
	fo.write('\n')
	for j, pn in enumerate(protein_cg_name_list):
		fo.write('%8s '%pn)
		for k, rn in enumerate(ribo_cg_name_list):
			if dist_map_L24[i][j][k] == '-':
				fo.write('%8s '%dist_map_L24[i][j][k])
			else:
				fo.write('%8.4f '%dist_map_L24[i][j][k])
		fo.write('\n')
	fo.write('\n')
	fo.close()
	########## get ribosomal protein collision diameter ##########
	for prot_atm_idx, prot_atm_name in enumerate(protein_cg_name_list):
		for atm_1 in pdb_struct.atoms:
			if atm_1.name == prot_atm_name and atm_1.residue.segid != 'NC':
				segid_1 = atm_1.residue.segid
				for atm_2 in pdb_struct.atoms:
					if atm_2.residue.segid != segid_1 and not (atm_2.residue.segid in ['NC', '23S', '5S', 'AtR', 'PtR', 'EtR']):
						if collision_diameter[i][prot_atm_idx] == '-':
							collision_diameter[i][prot_atm_idx] = 9999999
						dist = ((atm_1.xx - atm_2.xx)**2 + (atm_1.xy - atm_2.xy)**2 + (atm_1.xz - atm_2.xz)**2)**0.5
						if collision_diameter[i][prot_atm_idx] > dist:
							collision_diameter[i][prot_atm_idx] = dist
	fo = open('RP_results.dat', 'a')
	fo.write('PDB: %s\n'%pdb)
	for j, pn in enumerate(protein_cg_name_list):
		fo.write('%8s '%pn)
		if collision_diameter[i][j] == '-':
			fo.write('%8s\n'%collision_diameter[i][j])
		else:
			fo.write('%8.4f\n'%collision_diameter[i][j])
	fo.write('\n')
	fo.close()



########## Average ###########
NC_avg = [['-' for i in (ribo_cg_name_list)] for j in protein_cg_name_list]
L24_avg = [['-' for i in (ribo_cg_name_list)] for j in protein_cg_name_list]
collision_diameter_avg = ['-' for j in protein_cg_name_list]
fo = open('NC_results.dat', 'a')
fo.write('D_min_AVERAGE\n')
fo.write('%8s '%' ')
for rn in ribo_cg_name_list:
	fo.write('%8s '%rn)
fo.write('\n')
for j, pn in enumerate(protein_cg_name_list):
	fo.write('%8s '%pn)
	for k, rn in enumerate(ribo_cg_name_list):
		avg = 0
		n_avg = 0
		for i, pdb in enumerate(pdb_list):
			if dist_map_NC[i][j][k] != '-':
				avg += dist_map_NC[i][j][k]
				n_avg += 1
		if n_avg == 0:
			avg = '-'
		else:
			avg = avg / n_avg
			NC_avg[j][k] = avg
		if avg == '-':
			fo.write('%8s '%avg)
		else:
			fo.write('%8.4f '%avg)
	fo.write('\n')
fo.write('\n')
fo.close()

fo = open('L24_results.dat', 'a')
fo.write('D_min_AVERAGE\n')
fo.write('%8s '%' ')
for rn in ribo_cg_name_list:
	fo.write('%8s '%rn)
fo.write('\n')
for j, pn in enumerate(protein_cg_name_list):
	fo.write('%8s '%pn)
	for k, rn in enumerate(ribo_cg_name_list):
		avg = 0
		n_avg = 0
		for i, pdb in enumerate(pdb_list):
			if dist_map_L24[i][j][k] != '-':
				avg += dist_map_L24[i][j][k]
				n_avg += 1
		if n_avg == 0:
			avg = '-'
		else:
			avg = avg / n_avg
			L24_avg[j][k] = avg
		if avg == '-':
			fo.write('%8s '%avg)
		else:
			fo.write('%8.4f '%avg)
	fo.write('\n')
fo.write('\n')
fo.close()

fo = open('RP_results.dat', 'a')
fo.write('COLLISION DIAMETER\n')
for j, pn in enumerate(protein_cg_name_list):
	fo.write('%8s '%pn)
	avg = 0
	n_avg = 0
	for i, pdb in enumerate(pdb_list):
		if collision_diameter[i][j] != '-':
			avg += collision_diameter[i][j]
			n_avg += 1
	if n_avg == 0:
		avg = '-'
	else:
		avg = avg / n_avg
		collision_diameter_avg[j] = avg
	if avg == '-':
		fo.write('%8s '%avg)
	else:
		fo.write('%8.4f '%avg)
	fo.write('\n')
fo.write('\n')

protein_cg_R2 = []
for i in range(len(protein_cg_name_list)):
	protein_cg_R2.append(collision_diameter_avg[i]*2**(1/6)/2)

fo.write('Rmin/2\n')
for j, pn in enumerate(protein_cg_name_list):
	fo.write('%8s '%pn)
	if protein_cg_R2[j] == '-':
		fo.write('%8s '%protein_cg_R2[j])
	else:
		fo.write('%8.4f '%protein_cg_R2[j])
	fo.write('\n')
fo.write('\n')
fo.close()

#protein_cg_R = protein_cg_R2

######## opt FF #########
NC_FF = [['-' for i in range(3)] for j in protein_cg_name_list]
fo = open('NC_results.dat', 'a')
fo.write('FF\n')
fo.write('%10s %10s %10s %10s\n'%(' ', 'P', 'R', 'BR'))
for j, pn in enumerate(protein_cg_name_list):
	fo.write('%10s '%pn)
	for k, rn in enumerate(ribo_cg_name_list[0:2]):
		if NC_avg[j][k] == '-':
			fo.write('%10s '%NC_avg[j][k])
		else:
			R = NC_avg[j][k]*2**(1/6)-protein_cg_R[j]
			fo.write('%10.6f '%R)
			NC_FF[j][k] = R
	avg = 0
	n_avg = 0
	for k in range(2, len(ribo_cg_name_list)):
		if NC_avg[j][k] != '-':
			avg += NC_avg[j][k]
			n_avg += 1
	if n_avg == 0:
		avg = '-'
	else:
		avg = avg / n_avg
	if avg == '-':
		fo.write('%10s '%avg)
	else:
		R = avg*2**(1/6)-protein_cg_R[j]
		fo.write('%10.6f '%R)
		NC_FF[j][2] = R
	fo.write('\n')
fo.write('\n')
fo.close()

L24_FF = [['-' for i in range(3)] for j in protein_cg_name_list]
fo = open('L24_results.dat', 'a')
fo.write('FF\n')
fo.write('%10s %10s %10s %10s\n'%(' ', 'P', 'R', 'BR'))
for j, pn in enumerate(protein_cg_name_list):
	fo.write('%10s '%pn)
	for k, rn in enumerate(ribo_cg_name_list[0:2]):
		if L24_avg[j][k] == '-':
			fo.write('%10s '%L24_avg[j][k])
		else:
			R = L24_avg[j][k]*2**(1/6)-protein_cg_R[j]
			fo.write('%10.6f '%R)
			L24_FF[j][k] = R
	avg = 0
	n_avg = 0
	for k in range(2, len(ribo_cg_name_list)):
		if L24_avg[j][k] != '-':
			avg += L24_avg[j][k]
			n_avg += 1
	if n_avg == 0:
		avg = '-'
	else:
		avg = avg / n_avg
	if avg == '-':
		fo.write('%10s '%avg)
	else:
		R = avg*2**(1/6)-protein_cg_R[j]
		fo.write('%10.6f '%R)
		L24_FF[j][2] = R
	fo.write('\n')
fo.write('\n')
fo.close()

######## Avg opt FF #########
fo = open('NC_results.dat', 'a')
fo.write('AVERAGE FF\n')
str_list=['P', 'R', 'BR']
for k in range(3):
	avg = 0
	n_avg = 0
	for j, pn in enumerate(protein_cg_name_list):
		if NC_FF[j][k] != '-':
			avg += NC_FF[j][k]
			n_avg += 1
	if n_avg == 0:
		avg = '-'
	else:
		avg = avg / n_avg
	if avg == '-':
		fo.write('%-8s 0.0  -0.000132 %11s\n'%(str_list[k], avg))
	else:
		fo.write('%-8s 0.0  -0.000132 %11.6f\n'%(str_list[k], avg))
for k in range(len(protein_cg_name_list)):
	fo.write('%-8s 0.0  -0.000132 %11.6f\n'%(protein_cg_name_list[k], protein_cg_R[k]))
fo.close()

fo = open('L24_results.dat', 'a')
fo.write('AVERAGE FF\n')
fo.write('%10s %10s %10s\n'%('P', 'R', 'BR'))
for k in range(3):
	avg = 0
	n_avg = 0
	for j, pn in enumerate(protein_cg_name_list):
		if L24_FF[j][k] != '-':
			avg += L24_FF[j][k]
			n_avg += 1
	if n_avg == 0:
		avg = '-'
	else:
		avg = avg / n_avg
	if avg == '-':
		fo.write('%10s '%avg)
	else:
		fo.write('%10.6f '%avg)
fo.write('\n')
fo.close()
