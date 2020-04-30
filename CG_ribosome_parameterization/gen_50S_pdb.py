#!/usr/bin/env python3
import getopt, os, sys, time, multiprocessing, random, math
import parmed as pmd
import numpy as np

usage = '\nUsage: python gen_50S_pdb.py\n' \
        '       --ctrlfile | -f <50S.ctrl> Control file for generating 50S pdb\n'\
        '       [--help | -h] Print this information\n\n'\
        ' Example of cntrol file:\n'\
        '  cif_file = 4v9d.cif\n'\
        '  sub_unit = PtR 23S 5S L2 L3 L4 L5 L6 L9 L11 L13 L14 L15 L16 L17 L18 L19 L20 L21 L22 L23 L24 L25 L27 L28 L29 L30 L32 L33 L34 L35 L36\n'\
        '  chain_id = BV DA DB DC DD DE DF DG DH DI DJ DK DL DM DN DO DP DQ DR DS DT DU DV DW DX DY DZ D0 D1 D2 D3 D4\n'\
        '  new_chain_id = A B C D E F G H I J K L M N O P Q R S T U V W X Y Z 0 1 2 3 4 5\n'

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
        if line.startswith('cif_file'):
            words = line.split('=')
            cif_file = words[1].strip()
            continue
        if line.startswith('sub_unit'):
            words = line.split('=')[1].strip().split()
            sub_unit = words
            continue
        if line.startswith('chain_id'):
            words = line.split('=')[1].strip().split()
            chain_id = words
            continue
        if line.startswith('new_chain_id'):
            words = line.split('=')[1].strip().split()
            new_chain_id = words
            continue
finally:
     file_object.close()

if len(sub_unit) != len(chain_id):
	print('Error: sub_unit length %d and chain_id length %d mismatch'%(len(sub_unit), len(chain_id)))
elif len(sub_unit) != len(new_chain_id):
	print('Error: sub_unit length %d and new_chain_id length %d mismatch'%(len(sub_unit), len(new_chain_id)))
elif len(chain_id) != len(new_chain_id):
	print('Error: chain_id length %d and new_chain_id length %d mismatch'%(len(chain_id), len(new_chain_id)))

print('--> Reading cif file %s'%cif_file)
cif = pmd.load_file(cif_file)
print('    Done')
	
sub_unit_struct = {}

df = cif.to_dataframe()

resid_list = []

print('--> Spliting cif file')
for i in range(len(sub_unit)):
	struct = cif[(df.chain == chain_id[i]) & (df.resname != 'HOH') & (df.resname != 'ZN') & (df.resname != 'MG')]
	for res in struct.residues:
		res.chain = new_chain_id[i]
		res.segid = sub_unit[i]
		resid_list.append(res.number)
	sub_unit_struct[sub_unit[i]] = struct
print('    Done')

print('--> Combining 50S subunits')
comp = sub_unit_struct[sub_unit[0]]
for i in range(1, len(sub_unit)):
	comp = comp + sub_unit_struct[sub_unit[i]]

for res in comp.residues:
	res.number = resid_list[res.idx]
for atm in comp.atoms:
	atm.number = atm.idx + 1
print('    Done')

comp.write_pdb(cif_file.split('.cif')[0] + '_50S_tRNA.pdb', renumber=False, charmm=True)
