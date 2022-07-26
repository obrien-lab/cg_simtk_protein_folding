#!/usr/bin/env python3
import sys, getopt, math, os
import xml.etree.cElementTree as ET
import xml.dom.minidom as MD
import numpy
import parmed as pmd


usage = 'Usage: python parse_cg_prm.py\n' \
		'              --prmfile | -p <CG.prm> Charmm prm file for CG model\n'\
		'              --topfile | -t <CG.top> Charmm top file for CG model\n'

prmfile = ''
topfile = ''

if len(sys.argv) == 1:
	print(usage)
	sys.exit()

try:
	opts, args = getopt.getopt(sys.argv[1:],"hp:t:",["prmfile=", "topfile="])
except getopt.GetoptError:
	print(usage)
	sys.exit()
for opt, arg in opts:
	if opt == '-h':
		print(usage)
		sys.exit()
	elif opt in ("-p", "--prmfile"):
		prmfile = arg
	elif opt in ("-t", "--topfile"):
		topfile = arg

top_file_list = topfile.strip().split()

command = 'pmd.charmm.CharmmParameterSet('
for tf in top_file_list:
	command += '"'+tf + '", '
command += 'prmfile)'

param=eval(command)

name = prmfile.split('.prm')
file_name = name[0]

openmm_param=pmd.openmm.parameters.OpenMMParameterSet.from_parameterset(param)
openmm_param.write(file_name+'_tmp.xml', skip_duplicates=False)
dom = MD.parse(file_name+'_tmp.xml')
root = dom.documentElement
atom_type = root.getElementsByTagName('AtomTypes')
residue = root.getElementsByTagName('Residues')
os.remove(file_name+'_tmp.xml')

root = ET.Element("ForceField")

pf = open(prmfile, 'r')
section = None
node = None
nbxmod = None
ep = None
kc = 138.935485
ld = 1 # 10 Angstrom
atom_type_list = [];
num_atom = 0
dihedral_array = []
acoef_array = None
bcoef_array = None
ccoef_array = None
nb_table = []
nbfix_table = []
try:
	for line in pf:
		line = line.strip()
		if not line:
			# This is a blank line
			continue
		if line.startswith('!'):
			# This is a comment line
			continue
		if line.startswith('ATOM'):
			section = 'ATOM'
			#node = ET.SubElement(root, "AtomTypes")
			continue
		if line.startswith('BOND'):
			section = 'BOND'
			node = ET.SubElement(root, 'HarmonicBondForce')
			continue
		if line.startswith('ANGLE'):
			section = 'ANGLE'
			node = ET.SubElement(root, 'CustomAngleForce', 
				energy='-1/gamma*log(e); e=exp(-gamma*(k_alpha*(theta-theta_alpha)^2+epsilon_alpha))+exp(-gamma*k_betta*(theta-theta_betta)^2)')
			ET.SubElement(node, 'PerAngleParameter', name='k_alpha')
			ET.SubElement(node, 'PerAngleParameter', name='theta_alpha')
			ET.SubElement(node, 'PerAngleParameter', name='k_betta')
			ET.SubElement(node, 'PerAngleParameter', name='theta_betta')
			ET.SubElement(node, 'PerAngleParameter', name='gamma')
			ET.SubElement(node, 'PerAngleParameter', name='epsilon_alpha')
			continue
		if line.startswith('DIHEDRAL'):
			section = 'DIHEDRAL'
			node = ET.SubElement(root, 'PeriodicTorsionForce')
			continue
		if line.startswith('IMPHI'):
			section = 'IMPROPER'
			if len(dihedral_array) != 0:
				proper_node = ET.SubElement(node, 'Proper', type1=dihedral_array[0], type2=dihedral_array[1], 
						type3=dihedral_array[2], type4=dihedral_array[3])
			n0 = 1
			for index in range(4, len(dihedral_array), 3):
				proper_node.set('k'+str(n0), dihedral_array[index])
				proper_node.set('periodicity'+str(n0), dihedral_array[index+1])
				proper_node.set('phase'+str(n0), dihedral_array[index+2])
				n0 += 1
			node = ET.SubElement(root, 'CustomTorsionForce', 
				energy='k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535')
			ET.SubElement(node, 'PerTorsionParameter', name='k')
			ET.SubElement(node, 'PerTorsionParameter', name='theta0')
			continue
		if line.startswith('NONBONDED'):
			section = 'NONBONDED'
			words = line.split()
			nbxmod = int(words[2])
			continue
		if line.startswith('CUTNB'):
			words = line.split()
			ep = float(words[7])
			node = ET.SubElement(root, 'CustomNonbondedForce', 
				energy='ke*charge1*charge2/ep/r*exp(-r/ld)+kv*(a/r^12+b/r^10+c/r^6); '+
				'ke=ke1*ke2; ep=ep1*ep2; ld=ld1*ld2; kv=kv1*kv2; '+
				'a=acoef(index1, index2); b=bcoef(index1, index2); c=ccoef(index1, index2)',
				bondCutoff=str(nbxmod-1))
			ET.SubElement(node, 'PerParticleParameter', name='ke')
			ET.SubElement(node, 'PerParticleParameter', name='kv')
			ET.SubElement(node, 'PerParticleParameter', name='ep')
			ET.SubElement(node, 'PerParticleParameter', name='ld')
			ET.SubElement(node, 'PerParticleParameter', name='charge')
			ET.SubElement(node, 'PerParticleParameter', name='index')
			acoef_array = numpy.zeros((num_atom, num_atom))
			bcoef_array = numpy.zeros((num_atom, num_atom))
			ccoef_array = numpy.zeros((num_atom, num_atom))
			nb_table = [[] for i in atom_type_list]
			continue
		if line.startswith('NBFIX'):
			section = 'NBFIX'
			continue
		# It seems like files? sections? can be terminated with 'END'
		if line.startswith('END'): # should this be case-insensitive?
			section = None
			continue
		# If we have no section, skip
		if section is None: continue
		# Now handle each section specifically
		if section == 'ATOM':
			words = line.split()
			idx = int(words[1])
			name = words[2]
			mass = float(words[3])
			#atom_node = ET.SubElement(node, 'Type', name=name, element='C', mass=str(mass))
			#atom_node.set('class', name)
			num_atom += 1
			atom_type_list.append(name)
		if section == 'BOND':
			words = line.split()
			ET.SubElement(node, 'Bond', type1=words[0], type2=words[1], length=str(float(words[3])/10), k=str(float(words[2])*4.184*100*2))
		if section == 'ANGLE':
			words = line.split()
			ET.SubElement(node, 'Angle', type1=words[0], type2=words[1], type3= words[2],
				k_alpha=str(float(words[3])*4.184), theta_alpha=str(float(words[4])/180*math.pi),
				k_betta=str(float(words[5])*4.184), theta_betta=str(float(words[6])/180*math.pi),
				gamma=str(float(words[7])/4.184), epsilon_alpha=str(float(words[8])*4.184))
		if section == 'DIHEDRAL':
			words = line.split()
			type1 = words[0]
			type2 = words[1]
			type3 = words[2]
			type4 = words[3]
			k = str(float(words[4])*4.184)
			n = words[5]
			phase = str(float(words[6])/180*math.pi)
			if len(dihedral_array) == 0:
				dihedral_array.append(type1)
				dihedral_array.append(type2)
				dihedral_array.append(type3)
				dihedral_array.append(type4)
				dihedral_array.append(k)
				dihedral_array.append(n)
				dihedral_array.append(phase)
			elif (type1 == dihedral_array[0] and type2 == dihedral_array[1] and 
				type3 == dihedral_array[2] and type4 == dihedral_array[3]):
				dihedral_array.append(k)
				dihedral_array.append(n)
				dihedral_array.append(phase)
			else:
				proper_node = ET.SubElement(node, 'Proper', type1=dihedral_array[0], type2=dihedral_array[1], 
					type3=dihedral_array[2], type4=dihedral_array[3])
				n0 = 1
				for index in range(4, len(dihedral_array), 3):
					proper_node.set('k'+str(n0), dihedral_array[index])
					proper_node.set('periodicity'+str(n0), dihedral_array[index+1])
					proper_node.set('phase'+str(n0), dihedral_array[index+2])
					n0 += 1
				dihedral_array = []
				dihedral_array.append(type1)
				dihedral_array.append(type2)
				dihedral_array.append(type3)
				dihedral_array.append(type4)
				dihedral_array.append(k)
				dihedral_array.append(n)
				dihedral_array.append(phase)
		if section == 'IMPROPER':
			# No improper torsion energy term for Ca model
			continue
		if section == 'NONBONDED':
			words = line.split()
			name = words[0]
			epsilon = -float(words[2])*4.184
			R_min_half = float(words[3])/10
			index = atom_type_list.index(name)
			nb_table[index] = [epsilon, R_min_half]
		if section == 'NBFIX':
			words = line.split()
			type1 = words[0]
			type2 = words[1]
			index1 = atom_type_list.index(type1)
			index2 = atom_type_list.index(type2)
			epsilon = -float(words[2])*4.184
			R_min_half = float(words[3])/10
			nbfix_table.append([index1, index2, epsilon, R_min_half])
finally:
	pf.close()

#Build acoef, bcoef, ccoef tables
for index1 in range(num_atom):
	epsilon1 = nb_table[index1][0]
	R_min1 = nb_table[index1][1]
	for index2 in range(num_atom):
		epsilon2 = nb_table[index2][0]
		R_min2 = nb_table[index2][1]
		epsilon = numpy.sqrt(epsilon1 * epsilon2)
		R_min = R_min1 + R_min2
		a = 13 * epsilon * pow(R_min, 12)
		b = -18 * epsilon * pow(R_min, 10)
		c = 4 * epsilon * pow(R_min, 6)
		acoef_array[index1, index2] = a
		bcoef_array[index1, index2] = b
		ccoef_array[index1, index2] = c
for nbfix_list in nbfix_table:
	index1 = nbfix_list[0]
	index2 = nbfix_list[1]
	epsilon = nbfix_list[2]
	R_min = nbfix_list[3]
	a = 13 * epsilon * pow(R_min, 12)
	b = -18 * epsilon * pow(R_min, 10)
	c = 4 * epsilon * pow(R_min, 6)
	acoef_array[index1, index2] = a
	bcoef_array[index1, index2] = b
	ccoef_array[index1, index2] = c
	acoef_array[index2, index1] = a
	bcoef_array[index2, index1] = b
	ccoef_array[index2, index1] = c

#build tabulated function for acoef, bcoef, ccoef
acoef_node = ET.SubElement(node, "Function", name='acoef', type='Discrete2D',
	xsize=str(num_atom), ysize=str(num_atom))
text = ''
for index1 in range(num_atom):
	for index2 in range(num_atom):
		text += str(acoef_array[index1, index2]) + " "
acoef_node.text = text

bcoef_node = ET.SubElement(node, "Function", name='bcoef', type='Discrete2D',
	xsize=str(num_atom), ysize=str(num_atom))
text = ''
for index1 in range(num_atom):
	for index2 in range(num_atom):
		text += str(bcoef_array[index1, index2]) + " "
bcoef_node.text = text

ccoef_node = ET.SubElement(node, "Function", name='ccoef', type='Discrete2D',
	xsize=str(num_atom), ysize=str(num_atom))
text = ''
for index1 in range(num_atom):
	for index2 in range(num_atom):
		text += str(ccoef_array[index1, index2]) + " "
ccoef_node.text = text

#add custom nonbond parameters
ET.SubElement(node, 'UseAttributeFromResidue', name='charge')
for index in range(num_atom):
	name = atom_type_list[index]
	ET.SubElement(node, 'Atom', type=name, index=str(index), ke=str(kc**0.5), ep=str(ep**0.5), ld=str(ld**0.5), kv='1')

dom = MD.parseString(ET.tostring(root))
root = dom.documentElement
root = root.toprettyxml(indent=' ', newl='\n')
dom = MD.parseString(root)
root = dom.documentElement
bond = root.getElementsByTagName('HarmonicBondForce')
root.insertBefore(atom_type[0], bond[0])
if len(residue) > 0:
	root.insertBefore(residue[0], bond[0])

xf = open(file_name+'.xml', 'w')
#dom.writexml(xf, indent='', addindent=' ', newl='\n')
dom.writexml(xf, indent='')
