#!/usr/bin/env python3
import getopt, os, sys, time, multiprocessing, random, math
import parmed as pmd
import numpy as np

usage = '\nUsage: python truncate_ribosome.py\n' \
        '       --psf | -p <ribosome.psf> input psf file\n'\
        '       --cor | -c <ribosome.cor> input cor file\n'\
        '       [--tunnel | -t] <tunnel_def.pdb> tunnel defined by CAVER 3.0\n'\
        '                       If not defined, use x-axis as the tunnel central line.\n'\
        '       [-h] Print this information\n\n'

#################################### MAIN ####################################
psf_file = ''
cor_file = ''
tunnel_file = ''

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hp:c:t:", ["psf=", "cor=", "tunnel="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-p", "--psf"):
        psf_file = arg
    elif opt in ("-c", "--cor"):
        cor_file = arg
    elif opt in ("-t", "--tunnel"):
        tunnel_file = arg

if psf_file == '' or cor_file == '':
	print(usage)
	sys.exit()

print('--> Reading ribosome structure')
structure = pmd.charmm.psf.CharmmPsfFile(psf_file)
cor = pmd.load_file(cor_file)
structure.positions = cor.positions

coor = structure.positions

segid_list = []
resid_list = []
chainid_list = []
mask_array = [0 for i in range(len(structure.atoms))]
print('--> Done.')


## Read tunnel defination ##
if tunnel_file != '':
	print('--> Reading tunnel defination')
	tunnel_struct = pmd.load_file(tunnel_file)
	tunnel_cor = tunnel_struct.positions.value_in_unit(pmd.unit.angstroms)
	print('    Done.')
## END Read tunnel defination ##

print('--> Truncating ribosome')
## Truncate ribosome ##
for res in structure.residues:
	tag = 0
	for atom in res.atoms:
		if tunnel_file != '':
			dx_min = 9999999
			idx_min = 0
			for i in range(len(tunnel_cor)):
				if dx_min > abs(tunnel_cor[i][0]-atom.xx):
					dx_min = abs(tunnel_cor[i][0]-atom.xx)
					idx_min = i
			d = ((atom.xy-tunnel_cor[idx_min][1]) ** 2 + (atom.xz-tunnel_cor[idx_min][2]) ** 2) ** 0.5
		else:
			d = (atom.xy ** 2 + atom.xz ** 2) ** 0.5
		if d <= 30 and atom.xx >= -8 and atom.xx <= 58:
			tag = 1
			break
		elif atom.xx >= 58:
			tag = 1
			break
	if tag == 1:
		for atom in res.atoms:
			mask_array[atom.idx] = 1
		segid_list.append(res.segid)
		resid_list.append(res.number)
		chainid_list.append(res.chain)

new_structure = structure[mask_array]
print('    Done.')
## END Truncate ribosome ##

print('--> Writing outputs')
for res in new_structure.residues:
	res.number = resid_list[res.idx]
	res.segid = segid_list[res.idx]
	res.chain = chainid_list[res.idx]
for atm in new_structure.atoms:
	atm.number = atm.idx + 1

out_name = psf_file.split('.psf')[0]

new_structure.write_psf(out_name+'_truncated.psf')
new_structure.save(out_name+'_truncated.cor', format='charmmcrd', overwrite=True)
print('    Done.')
