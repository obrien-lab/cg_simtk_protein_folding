#!/usr/bin/env python3
import getopt, os, sys, time, multiprocessing, random, math
import parmed as pmd
import numpy as np

usage = '\nUsage: python create_cg_ribosome_model.py\n' \
        '       --pdbfile | -f <PDB.pdb> input pdb file\n'\
        '       --cg_protein_segid | -s <"SEGID LIST"> for cg protein segid\n'\
        '       --cg_protein_residue_prefix | -r <"RES NAME LIST"> for cg protein residue name\n'\
        '       --outname | -o <OUT FILE NAME> for output top, psf and pdb\n'\
        '       [-h] Print this information\n\n'

# calculate centroid position of a set of coordinates
def calc_average_position(atom_cor_list):
	avg_cor = pmd.unit.Quantity((0,0,0), pmd.unit.angstroms)
	for cor in atom_cor_list:
		avg_cor += cor
	avg_cor /= len(atom_cor_list)
	return avg_cor
# END calculate centroid position of a set of coordinates

# generate charmm .top
def gen_rtf(struct):
	global out_name, cg_protein_segid_list, atom_list_map, cg_protein_residue_prefix
	fo = open(out_name+'.top', 'w')
	fo.write('* This CHARMM .top file describes CG model of ribosome 50S subunit\n*\n20 1\n')
	# MASS section
	fo.write('! backbone masses\n'+
		'MASS 1   P        95.000000\n'+
		'MASS 2   R        92.000000\n'+
		'MASS 3   BR       64.000000\n'+
		'MASS 4   SG       57.000000\n'+
		'MASS 5   SA       71.000000\n'+
		'MASS 6   SV       99.000000\n'+
		'MASS 7   SL       113.000000\n'+
		'MASS 8   SI       113.000000\n'+
		'MASS 9   SM       131.000000\n'+
		'MASS 10  SF       147.000000\n'+
		'MASS 11  SP       114.000000\n'+
		'MASS 12  SS       87.000000\n'+
		'MASS 13  ST       101.000000\n'+
		'MASS 14  SC       114.000000\n'+
		'MASS 15  SN       114.000000\n'+
		'MASS 16  SQ       128.000000\n'+
		'MASS 17  SY       163.000000\n'+
		'MASS 18  SW       186.000000\n'+
		'MASS 19  SD       114.000000\n'+
		'MASS 20  SE       128.000000\n'+
		'MASS 21  SH       114.000000\n'+
		'MASS 22  SK       128.000000\n'+
		'MASS 23  SR       114.000000\n')
	idx = 23
	for cg_protein_segid in cg_protein_segid_list:
		for atm in atom_list_map[cg_protein_segid]:
			idx += 1
			fo.write('MASS %-3s %-8s %.6f\n'%(str(idx), atm.residue.name, atm.mass))
	fo.write('\n')
	# residue section
	fo.write('RESI CGA        -1.0 ! Adenine\n'+
			 'GROUP\n'+
			 'ATOM P  P       -1.0\n'+
			 'ATOM R  R        0.0\n'+
			 'ATOM PU1 BR      0.0\n'+
			 'ATOM PU2 BR      0.0\n\n')
	fo.write('RESI CGU        -1.0 ! Uracil\n'+
			 'GROUP\n'+
			 'ATOM P  P       -1.0\n'+
			 'ATOM R  R        0.0\n'+
			 'ATOM PY BR       0.0\n\n')
	fo.write('RESI CGG        -1.0 ! Guanine\n'+
			 'GROUP\n'+
			 'ATOM P  P       -1.0\n'+
			 'ATOM R  R        0.0\n'+
			 'ATOM PU1 BR      0.0\n'+
			 'ATOM PU2 BR      0.0\n\n')
	fo.write('RESI CGC        -1.0 ! Cytosine\n'+
			 'GROUP\n'+
			 'ATOM P  P       -1.0\n'+
			 'ATOM R  R        0.0\n'+
			 'ATOM PY BR       0.0\n\n')
	fo.write('RESI SG          0.0 ! GLY alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SG  SG      0.0\n\n')
	fo.write('RESI SA          0.0 ! ALA alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SA  SA      0.0\n\n')
	fo.write('RESI SV          0.0 ! VAL alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SV  SV      0.0\n\n')
	fo.write('RESI SL          0.0 ! LEU alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SL  SL      0.0\n\n')
	fo.write('RESI SI          0.0 ! ILE alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SI  SI      0.0\n\n')
	fo.write('RESI SM          0.0 ! MET alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SM  SM      0.0\n\n')
	fo.write('RESI SF          0.0 ! PHE alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SF  SF      0.0\n\n')
	fo.write('RESI SP          0.0 ! PRO alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SP  SP      0.0\n\n')
	fo.write('RESI SS          0.0 ! SER alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SS  SS      0.0\n\n')
	fo.write('RESI ST          0.0 ! THR alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM ST  ST      0.0\n\n')
	fo.write('RESI SC          0.0 ! CYS alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SC  SC      0.0\n\n')
	fo.write('RESI SN          0.0 ! ASN alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SN  SN      0.0\n\n')
	fo.write('RESI SQ          0.0 ! GLN alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SQ  SQ      0.0\n\n')
	fo.write('RESI SY          0.0 ! TYR alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SY  SY      0.0\n\n')
	fo.write('RESI SW          0.0 ! TRP alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SW  SW      0.0\n\n')
	fo.write('RESI SD         -1.0 ! ASP alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SD  SD     -1.0\n\n')
	fo.write('RESI SE         -1.0 ! GLU alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SE  SE     -1.0\n\n')
	fo.write('RESI SH          0.0 ! HIS alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SH  SH      0.0\n\n')
	fo.write('RESI SK          1.0 ! LYS alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SK  SK      1.0\n\n')
	fo.write('RESI SR          1.0 ! ARG alpha carbon\n'+
			 'GROUP\n'+
			 'ATOM SR  SR      1.0\n\n')
	for cg_protein_segid in cg_protein_segid_list:
		atm_name = cg_protein_residue_prefix[cg_protein_segid_list.index(cg_protein_segid)]
		fo.write('DECL +%s\n'%atm_name)
		fo.write('DECL -%s\n'%atm_name)
		fo.write('DECL #%s\n'%atm_name)
		for atm in atom_list_map[cg_protein_segid]:
			fo.write('RESI %-6s %s\n'%(atm.residue.name, '%.1f'%atm.charge))
			fo.write('GROUP\n')
			fo.write('ATOM %s %-6s %s\n'%(atm.name, atm.residue.name, '%.1f'%atm.charge))
			fo.write('BOND %s +%s\n'%(atm.name, atm.name))
			fo.write('ANGL -%s %s +%s\n'%(atm.name, atm.name, atm.name))
			fo.write('DIHE -%s %s +%s #%s\n\n'%(atm.name, atm.name, atm.name, atm.name))
	# end section
	fo.write('PRES LINK         0.00 ! linkage for IMAGES or for joining segments\n'+
			 '                       ! 1 refers to previous (N terminal)\n'+
			 '                       ! 2 refers to next (C terminal)\n'+
			 '                       ! use in B patch statement\n'+
			 '                       ! follow with AUTOgenerate ANGLes DIHEdrals command\n'+
			 'BOND 1A 2A\n\n'+
			 'END')
	fo.close()
# END generate charmm .top

# defination of data set
AA_resname_list = ["GLY","ALA","VAL","LEU","ILE","MET","PHE","PRO","SER","THR","CYS","ASN",
				   "GLN","TYR","TRP","ASP","GLU","HIS","HSD","LYS","ARG"]
RNA_resname_list = ["A","U","G","C","ADE","URA","GUA","CYT","THY"]
AA_resname_dict = {"ALA" : "A",
				   "CYS" : "C",
				   "ASP" : "D",
				   "GLU" : "E",
				   "PHE" : "F",
				   "GLY" : "G",
				   "HIS" : "H",
				   "HSD" : "H",
				   "HSE" : "H",
				   "HSP" : "H",
				   "ILE" : "I",
				   "LYS" : "K",
				   "LEU" : "L",
				   "MET" : "M",
				   "ASN" : "N",
				   "PRO" : "P",
				   "GLN" : "Q",
				   "ARG" : "R",
				   "SER" : "S",
				   "THR" : "T",
				   "VAL" : "V",
				   "TRP" : "W",
				   "TYR" : "Y"}
charge_dict = {"ALA" : 0,
			   "CYS" : 0,
			   "ASP" : -1,
			   "GLU" : -1,
			   "PHE" : 0,
			   "GLY" : 0,
			   "HIS" : 0,
			   "HSD" : 0,
			   "HSE" : 0,
			   "HSP" : 1,
			   "ILE" : 0,
			   "LYS" : 1,
			   "LEU" : 0,
			   "MET" : 0,
			   "ASN" : 0,
			   "PRO" : 0,
			   "GLN" : 0,
			   "ARG" : 1,
			   "SER" : 0,
			   "THR" : 0,
			   "VAL" : 0,
			   "TRP" : 0,
			   "TYR" : 0,
			   "P" : -1,
			   "R" : 0,
			   "PU1" : 0,
			   "PU2" : 0,
			   "PY" : 0}
AA_mass_dict = {"A" : 71.0,
				"C" : 114.0,
				"D" : 114.0,
				"E" : 128.0,
				"F" : 147.0,
				"G" : 57.0,
				"H" : 114.0,
				"I" : 113.0,
				"K" : 128.0,
				"L" : 113.0,
				"M" : 131.0,
				"N" : 114.0,
				"P" : 114.0,
				"Q" : 128.0,
				"R" : 114.0,
				"S" : 87.0,
				"T" : 101.0,
				"V" : 99.0,
				"W" : 186.0,
				"Y" : 163.0}
ribo_mass_dict = {"P" : 95.0,
				  "R" : 92.0,
				  "PU1" : 64.0,
				  "PU2" : 64.0,
				  "PY" : 64.0,}
# END defination of data set

#################################### MAIN ####################################
aa_pdb = ''
cg_protein_segid = ''
cg_protein_residue_prefix = ''
out_name = ''

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:s:r:o:", ["pdbfile=", "cg_protein_segid=", "cg_protein_residue_prefix=", "outname="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-f", "--pdbfile"):
        aa_pdb = arg
    elif opt in ("-s", "--cg_protein_segid"):
        cg_protein_segid = arg
    elif opt in ("-r", "--cg_protein_residue_prefix"):
        cg_protein_residue_prefix = arg
    elif opt in ("-o", "--outname"):
        out_name = arg

if aa_pdb == '' or out_name == '':
	print(usage)
	sys.exit()

cg_protein_map = {} # store the current resid for a cg protein segid (KEY)
atom_list_map = {} # store atoms for a cg protein segid (KEY)

cg_protein_segid_list = cg_protein_segid.strip().split()
cg_protein_residue_prefix_list = cg_protein_residue_prefix.strip().split()
for si in cg_protein_segid_list:
	cg_protein_map[si] = 0
	atom_list_map[si] = []

aa_struct = pmd.load_file(aa_pdb)
aa_position = aa_struct.positions
cg_struct = pmd.structure.Structure()
cg_position = pmd.unit.Quantity([], pmd.unit.angstroms)

for res in aa_struct.residues:
	if res.segid in cg_protein_segid_list:
		cg_protein_map[res.segid] += 1
		cg_atom_name = cg_protein_residue_prefix[cg_protein_segid_list.index(res.segid)]
		for atm in res.atoms:
			if atm.name == 'CA':
				cg_atom_cor = aa_position[atm.idx]
				cg_position.append(cg_atom_cor)
				cg_res_name = cg_atom_name+str(cg_protein_map[res.segid])
				cg_struct.add_atom(pmd.topologyobjects.Atom(name=cg_atom_name, type=cg_res_name, 
					charge=charge_dict[res.name], mass=AA_mass_dict[AA_resname_dict[res.name]]), 
					cg_res_name, res.number, chain=res.chain, inscode='', segid=res.segid)
				break
	else:
		if res.name in AA_resname_list:
			cg_atom_name = 'S' + AA_resname_dict[res.name]
			for atm in res.atoms:
				if atm.name == 'CA':
					cg_atom_cor = aa_position[atm.idx]
					cg_position.append(cg_atom_cor)
					cg_struct.add_atom(pmd.topologyobjects.Atom(name=cg_atom_name, type=cg_atom_name, 
						charge=charge_dict[res.name], mass=AA_mass_dict[AA_resname_dict[res.name]]), 
						cg_atom_name, res.number, chain=res.chain, inscode='', segid=res.segid)
					break
		elif res.name in RNA_resname_list:
			ribose_ring_cor_list = []
			for atm in res.atoms:
				if atm.name == 'P':
					cg_atom_name = 'P'
					cg_atom_cor = aa_position[atm.idx]
					cg_position.append(cg_atom_cor)
					cg_struct.add_atom(pmd.topologyobjects.Atom(name=cg_atom_name, type='P', 
						charge=charge_dict['P'], mass=ribo_mass_dict['P']), 
						'CG'+res.name, res.number, chain=res.chain, inscode='', segid=res.segid)
				elif atm.name == "C1'" or atm.name == "C2'" or atm.name == "C3'" or atm.name == "C4'" or atm.name == "C5'":
					ribose_ring_cor_list.append(aa_position[atm.idx])
			cg_struct.add_atom(pmd.topologyobjects.Atom(name='R', type='R', 
				charge=charge_dict['R'], mass=ribo_mass_dict['R']), 
				'CG'+res.name, res.number, chain=res.chain, inscode='', segid=res.segid)
			cg_position.append(calc_average_position(ribose_ring_cor_list))
			if res.name == 'A' or res.name == 'G':
				ring_cor_list_1 = []
				ring_cor_list_2 = []
				for atm in res.atoms:
					if atm.name == 'N1' or atm.name == 'N3' or atm.name == 'C2' or atm.name == 'C4' or atm.name == 'C5' or atm.name == 'C6':
						ring_cor_list_1.append(aa_position[atm.idx])
				for atm in res.atoms:
					if atm.name == 'N7' or atm.name == 'N9' or atm.name == 'C4' or atm.name == 'C5' or atm.name == 'C8':
						ring_cor_list_2.append(aa_position[atm.idx])
				cg_struct.add_atom(pmd.topologyobjects.Atom(name='PU1', type='BR', 
					charge=charge_dict['PU1'], mass=ribo_mass_dict['PU1']), 
					'CG'+res.name, res.number, chain=res.chain, inscode='', segid=res.segid)
				cg_position.append(calc_average_position(ring_cor_list_1))
				cg_struct.add_atom(pmd.topologyobjects.Atom(name='PU2', type='BR', 
					charge=charge_dict['PU2'], mass=ribo_mass_dict['PU2']), 
					'CG'+res.name, res.number, chain=res.chain, inscode='', segid=res.segid)
				cg_position.append(calc_average_position(ring_cor_list_2))
			elif res.name == 'U' or res.name == 'C':
				ring_cor_list_1 = []
				for atm in res.atoms:
					if atm.name == 'N1' or atm.name == 'N3' or atm.name == 'C2' or atm.name == 'C4' or atm.name == 'C5' or atm.name == 'C6':
						ring_cor_list_1.append(aa_position[atm.idx])
				cg_struct.add_atom(pmd.topologyobjects.Atom(name='PY', type='BR', 
					charge=charge_dict['PY'], mass=ribo_mass_dict['PY']), 
					'CG'+res.name, res.number, chain=res.chain, inscode='', segid=res.segid)
				cg_position.append(calc_average_position(ring_cor_list_1))

# renumber atoms
cg_struct.positions = cg_position
for atm in cg_struct.atoms:
	atm.number = atm.idx + 1
# END renumber atoms

# generate bonds, angles and dihedrals
bond_list = []
angle_list = []
dihedral_list = []
for atm in cg_struct.atoms:
	if atm.residue.segid in cg_protein_segid_list:
		atom_list_map[atm.residue.segid].append(atm)
for key in atom_list_map:
	for i in range(len(atom_list_map[key])-1):
		bond = pmd.topologyobjects.Bond(atom_list_map[key][i], atom_list_map[key][i+1])
		bond_list.append(bond)
	for i in range(len(atom_list_map[key])-2):
		angle = pmd.topologyobjects.Angle(atom_list_map[key][i], atom_list_map[key][i+1], atom_list_map[key][i+2])
		angle_list.append(angle)
	for i in range(len(atom_list_map[key])-3):
		dihedral = pmd.topologyobjects.Dihedral(atom_list_map[key][i], atom_list_map[key][i+1], atom_list_map[key][i+2], 
			atom_list_map[key][i+3])
		dihedral_list.append(dihedral)
cg_struct.bonds = bond_list
cg_struct.angles = angle_list
cg_struct.dihedrals = dihedral_list
# END generate bonds, angles and dihedrals

cg_struct.write_pdb(out_name+'.pdb', renumber=False, charmm=True)
cg_struct.write_psf(out_name+'.psf')
cg_struct.save(out_name+'.cor', format='charmmcrd', overwrite=True)
gen_rtf(cg_struct)
