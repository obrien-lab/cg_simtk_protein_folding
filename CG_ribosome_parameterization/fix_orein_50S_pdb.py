#!/usr/bin/env python3
import getopt, os, sys, time, multiprocessing, random, math
import parmed as pmd
import numpy as np
import pdbfixer as pfx

usage = '\nUsage: python fix_orein_50S_pdb.py\n' \
        '       --pdbfile | -f <PDB.pdb> input pdb file\n'\
        '       [--help | -h] Print this information\n\n'

########################################################
def rotate_coor(structure, vector1, vector2):
	theta = np.arccos(np.dot(vector1, vector2)/(np.linalg.norm(vector1)*(np.linalg.norm(vector2))))
	normal_vector = np.cross(vector1, vector2)/np.linalg.norm(np.cross(vector1, vector2))
	w = np.cos(theta/2)
	x = normal_vector[0]*np.sin(theta/2)
	y = normal_vector[1]*np.sin(theta/2)
	z = normal_vector[2]*np.sin(theta/2)
	Q = np.mat([[1-2*y*y-2*z*z, 2*x*y-2*z*w, 2*x*z+2*y*w], 
			    [2*x*y+2*z*w, 1-2*x*x-2*z*z, 2*y*z-2*x*w],
			    [2*x*z-2*y*w, 2*y*z+2*x*w, 1-2*x*x-2*y*y]])
	new_pos = pmd.unit.Quantity([], pmd.unit.angstroms)
	for pos in structure.positions:
		P = np.mat(pos.value_in_unit(pmd.unit.angstroms)).T
		pos = Q*P
		pos = pos.A1
		pos = pmd.unit.Quantity(list(pos), pmd.unit.angstroms)
		new_pos.append(pos)
	structure.positions = new_pos
#########################################################

def calc_average_position(atom_cor_list):
	avg_cor = pmd.unit.Quantity((0,0,0), pmd.unit.angstroms)
	for cor in atom_cor_list:
		avg_cor += cor
	avg_cor /= len(atom_cor_list)
	return avg_cor
########################################################

############################# MAIN #######################
pdbfile = ''

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:", ["pdbfile="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-f", "--pdbfile"):
        pdbfile = arg

nonstandard_res = ['2MA', '3AU', '4SU', '5MU', '6MZ', '7MG', 'CM0', 'PSU', 'QUO', 'T6A', 'U8U']
standard_res = ['A', 'U', 'U', 'U', 'A', 'G', 'U', 'U', 'G', 'A', 'U']
nonstandard_map = {}
for i in range(len(nonstandard_res)):
	nonstandard_map[nonstandard_res[i]] = standard_res[i]

fixer = pfx.PDBFixer(filename=pdbfile)
pdb_struct = pmd.load_file(pdbfile)
resid_list = []
segid_list = []
mutate_dict = {}
print('--> Checking nonstandard RNA linker')
for res in pdb_struct.residues:
	resid_list.append(res.number)
	segid_list.append(res.segid)
	if not (res.chain in mutate_dict):
		mutate_dict[res.chain] = []
	if res.name in nonstandard_res:
		mutate_dict[res.chain].append(res.name + '-' + str(res.number) + '-' + nonstandard_map[res.name])
		print('    ' + res.name + '-' + str(res.number) + '-' + nonstandard_map[res.name] + ' in chain ' + res.chain)

for chain in mutate_dict:
	fixer.applyMutations(mutate_dict[chain], chain)
print('    Done')

fixer.findMissingResidues()
fixer.findMissingAtoms()

miss_atm_dict = fixer.missingAtoms
print('--> Missing atoms:')
for res in miss_atm_dict:
	miss_atm_list = miss_atm_dict[res]
	resname = res.name
	resid = res.id
	chain_id = res.chain.id
	out_str = '    Chain %s residue %s %s: '%(chain_id, resname, resid)
	for atm in miss_atm_list:
		atmname = atm.name
		out_str += atmname + ' '
	print(out_str)

print('--> Add missing atoms')
fixer.addMissingAtoms()
print('    Done')

new_pdb_struct = pmd.openmm.load_topology(fixer.topology, xyz=fixer.positions)

print('--> Renumber new pdb')
for res in new_pdb_struct.residues:
	res.number = resid_list[res.idx]
	res.segid = segid_list[res.idx]
for atm in new_pdb_struct.atoms:
	atm.number = atm.idx + 1
print('    Done')

print('--> Translate new pdb as 23S A2602 N6 at origin')
idx_1 = 0
idx_2 = 0
for atm in new_pdb_struct.atoms:
	if atm.residue.segid == '23S' and atm.residue.number == 2602 and atm.name == 'N6':
		idx_1 = atm.idx
	elif atm.residue.segid == 'L24' and atm.residue.number == 51 and atm.name == 'N':
		idx_2 = atm.idx
coor_1 = new_pdb_struct.positions[idx_1]
new_pos = pmd.unit.Quantity([], pmd.unit.angstroms)
for pos in new_pdb_struct.positions:
	pos = pos - coor_1
	new_pos.append(pos)
new_pdb_struct.positions = new_pos
print('    Done')

print('--> Rotate new pdb')
print('    rotate coor to make vector from 23S:2602@N6 to L24:51@N along x-axis')
# rotate coor to make vector from 23S:2602@N6 to L24:51@N along x-axis 
coor_1 = new_pdb_struct.positions[idx_1]
coor_2 = new_pdb_struct.positions[idx_2]
vector1 = np.array(coor_2.value_in_unit(pmd.unit.angstroms) - coor_1.value_in_unit(pmd.unit.angstroms), dtype=np.float64)
vector2 = np.array([1, 0, 0], dtype=np.float64)
rotate_coor(new_pdb_struct, vector1, vector2)

tag_AtR = 0
for res in new_pdb_struct.residues:
	if res.segid == 'AtR':
		tag_AtR = 1
		break

if tag_AtR == 1:
	print('    rotate coor to make yz-component of vector from centroid of AtR:76@ribose_ring'+
		'to centroid of PtR:76@ribose_ring along y-axis')
	# rotate coor to make yz-component of vector from centroid of AtR:76@ribose_ring 
	# to centroid of PtR:76@ribose_ring along y-axis
	ribose_ring_cor_list = []
	for atm in new_pdb_struct.atoms:
		if atm.residue.number == 76 and atm.residue.segid == 'AtR':
			if atm.name == "C1'" or atm.name == "C2'" or atm.name == "C3'" or atm.name == "C4'" or atm.name == "C5'":
				ribose_ring_cor_list.append(new_pdb_struct.positions[atm.idx])
	coor_1 = calc_average_position(ribose_ring_cor_list)
	ribose_ring_cor_list = []
	for atm in new_pdb_struct.atoms:
		if atm.residue.number == 76 and atm.residue.segid == 'PtR':
			if atm.name == "C1'" or atm.name == "C2'" or atm.name == "C3'" or atm.name == "C4'" or atm.name == "C5'":
				ribose_ring_cor_list.append(new_pdb_struct.positions[atm.idx])
	coor_2 = calc_average_position(ribose_ring_cor_list)
	vector1 = np.array(coor_2.value_in_unit(pmd.unit.angstroms) - coor_1.value_in_unit(pmd.unit.angstroms), dtype=np.float64)
	vector1[0] = 0.0
	vector2 = np.array([0, 1, 0], dtype=np.float64)
	rotate_coor(new_pdb_struct, vector1, vector2)

print('    Done')

new_pdb_struct.write_pdb(pdbfile.split('.pdb')[0] + '_model.pdb', renumber=False, charmm=True)
