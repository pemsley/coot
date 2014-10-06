
import sys
import os
import copy
from subprocess import call
from rdkit import Chem
from rdkit.Chem import AllChem

import coot_svn_repo_revision
import pyrogen_swig as pysw
import libpyrogen_boost as pyrogen_boost

from optparse import OptionParser

import tautomer

from jay_util import *

global pyrogen_version
pyrogen_version = "0.0-pre"

global run_mogul
global smiles_dict
run_mogul = True
smiles_dict = False

def make_mogul_ins_file(mogul_ins_file_name, mogul_out_file_name, sdf_file_name):
   f = open(mogul_ins_file_name, 'w')
   if f:
     f.write('mogul molecule file ')
     f.write(sdf_file_name)
     f.write('\n')
     f.write('mogul output   file ')
     f.write(mogul_out_file_name)
     f.write('\n')
     f.write('mogul output distribution all on\n')
     f.write('bond all\n')
     f.write('angle all\n')
#      f.write('torsion all\n')
#      f.write('ring all\n')
     f.write('config output format CSV\n')
     f.write('config output items fragment_type atom_indices query_value nhits mean median sd z-score dmin\n')
     f.write('config search all filter exclude_solvents\n')
     f.write('config output invalid_fragments exclude\n')
     f.close()
   return f


# return True for good, False for bad/not-run
#
def execute_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name):
   f = make_mogul_ins_file(mogul_ins_file_name, mogul_out_file_name, sdf_file_name)
   if f: 
      # print 'now run mogul using ins file %s' % mogul_ins_file_name
      if run_mogul:
          call(['mogul', '-ins', mogul_ins_file_name])
          return True
      else:
          return False
   else:
      return False

def atom_name_from_atomic_number_and_count(element, count):
    name = element
    name += str(count)
    return name

def add_atom_names(mol):
    nz = {}
    atom_names = []
    for atom in mol.GetAtoms():
       try:
          n = atom.GetProp('name')
          atom_names.append(n)
       except KeyError:
          z = atom.GetAtomicNum()
          if z in nz:
             nz[z]  = nz[z] + 1
          else:
             nz[z] = 1;
          ele = atom.GetSymbol().upper()
          name = atom_name_from_atomic_number_and_count(ele, nz[z])
          p_name = pad_atom_name(name, ele)
          atom.SetProp("name", p_name)
          atom_names.append(p_name)
    return atom_names

def convert_to_coot_bond_type(rdkit_type):
    out_type = 'single'
    if (rdkit_type == Chem.rdchem.BondType.SINGLE):
        out_type = 'single'
    if (rdkit_type == Chem.rdchem.BondType.AROMATIC):
        out_type = 'aromatic'
    if (rdkit_type == Chem.rdchem.BondType.DOUBLE):
        out_type = 'double'
    if (rdkit_type == Chem.rdchem.BondType.TRIPLE):
        out_type = 'triple'
    if (rdkit_type == Chem.rdchem.BondType.ONEANDAHALF):
        out_type = 'deloc'
    return out_type

def pad_atom_name(name, element):

    padded = name
    if (len(element) == 1):
        if (len(name) == 2):
            padded = ' ' + name + ' '
        if (len(name) == 3):
            padded = ' ' + name
    if (len(element) == 2):
        if (len(name) == 2):
            padded = name + '  '
        if (len(name) == 3):
            padded = name + ' '
    return padded

def is_smiles_file(file_name):
    bits = file_name.rsplit(".")
    if len(bits) > 1:
        return bits[1] == 'smi'
    else:
        return False

def is_comp_id(comp_id):
   return len(comp_id) == 3

def is_mdl_file(file_name):
    bits = file_name.rsplit(".")
    if (len(bits) < 2):
       return False
    else:
       idx = len(bits) - 1
       if (bits[idx] == 'mol'):
          return True
       else:
          if (bits[idx] == 'mdl'):
             return True
          else:
             return False
        

# return the contents of file_name
def read_file(file_name):
    f = open(file_name)
    return f.read()

        
def make_restraints_for_bond_orders(mol):
    restraints = {}
    bond_list = []
    for bond in mol.GetBonds(): 
        type = bond.GetBondType()
        coot_bond_type = convert_to_coot_bond_type(type)
        at_1 = bond.GetBeginAtom()
        at_2 = bond.GetEndAtom()
        name_1 = at_1.GetProp('name')
        name_2 = at_2.GetProp('name')
        item = [name_1, name_2, coot_bond_type, 1.0, 1.0]
        bond_list.append(item)
    restraints['_chem_comp_bond'] = bond_list
    restraints['_chem_comp'] = [mol.GetProp('comp_id'),
                                mol.GetProp('comp_id'),
                                mol.GetProp('name'),
                                'non-polymer',
                                mol.GetNumAtoms(),
                                mol.GetNumAtoms(),
                                '.']
    return restraints

# match_atom_index can be of type int or a list - otherwise trouble.
#
# Note that atom_type properties can also have been set in hydrogen_transformations():
#
def set_atom_type(match, match_atom_index, mol, atom_type):
    # print "   trying to set ", match_atom_index, " of ", match
    try:
        this_atom = match[match_atom_index]
        # print "   this_atom:", this_atom
        try:
            current_type = mol.GetAtomWithIdx(this_atom).GetProp("atom_type")
            # print "   oops - atom ", mol.GetAtomWithIdx(this_atom).GetProp("name"), " already has type ", current_type
        except KeyError:
            mol.GetAtomWithIdx(this_atom).SetProp("atom_type", atom_type)
            name = mol.GetAtomWithIdx(this_atom).GetProp("name")
	    if False:
		print '    set atom number %s having name %s with type %s ' % (str(this_atom).rjust(2),
									       repr(name), repr(atom_type))
    except TypeError:
        for match_atom in match_atom_index:
            set_atom_type(match, match_atom, mol, atom_type)

def ele_to_smart(v):
    return (v.upper(), '['+v+']', 0)

# those not handled by hand-coding
def smarts_by_element():
   eles = [
      "He", "Li", "Be", "B",  "Ne", "Na", "Mg", "Al", 
      "Ar", "K", "Ca",  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
      "Zn", "Ga", "Ge", "As", "Se",      "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
      "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
      "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
      "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
      "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"]
   return map(ele_to_smart, eles)
   

def set_atom_types(mol):
    smarts_list = [

        # Full coverage for C, H, O.

        # Oxygen
        ('O2',  "[OX2;H0]", 0), # ester, Os between P and C are O2, not OP
        ('OP',  'O~P',   0),
        ('OS',  'O~S',   0),
        ('OB',  'O~B',   0),
        ('OC',  '*C(=O)[OH]', (2,3)), # carboxylic acid
        ('OC',  '*C(=O)O',    (2,3)), # carboxylate, doesn't match deloc bonds
        ('OH1', '[OH1]', 0), # alcohol
        ('O2',  "[oX2;H0]", 0), # ring oxygen
        ('O',   'O=*',   0), # carbonyl oxygen

        # OH2 no examples
        # OHA no examples
        # OHB no examples
        # OHC no examples
        # OC2 no exmampes
        
        # Fallback oxygen
        ('O',   'O',   0),

        # Carbon SP
        ("CSP1", '[H][C]#*',  1), # e.g. in 2GT
        ("CSP",  '[C]#[C]',   (0,1)),
        ("CSP",  '[C]#*',     0),
        
        # Carbon SP2
        ('CR56', 'c12aaaac1aaa2',  (0,5)), # works on indole
        ('CR56', 'c12aaaan1aaa2',  0), # same pattern as (below) N56, but catching first 56 atom
        ('CR66', 'c12aaaac1aaaa2', (0,5)),
        ('CR6',  'c12ccccc1OCO2',  (0,5)),   # mouse, fused atoms in 6-ring not non-arom 5-ring
        ('CR66', 'c12aaaac1AAAA2', (0,5)),   # one 6-ring aromatic, other not. Needed for XXX?
                                             # but makes a fail on 113.
        ('CR6',  'c12caccc1***2',  (0,5)),  # aromatic 6, (non-)aromatic 5, maybe this should be CR56?

        # note CR1  missing - can't find example
        #      CR1H missing - can't find example

        ('CR16', '[cr6;H1]',  0),
        ('CR6',  '[cr6;H0]',  0),
        ('CR15', '[cr5;H1]',  0),
#        ('CR5',  'C1(=O)[C,c][C,c]C=N1', 0), # carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        ('CR5',  '[cr5;H0]',  0),
        ('CR5',  '[CR5;H0]',  0),
        ('C1',   '[CX3;H1]',    0),  # double bond, single bond and one H
        ('C2',   '[CX3;H2]=*',  0),  # double bond, and 2 H
        ('C',    '[CX3;H0;^2]', 0),
        ('C',    '[CX3]=[OX1]', 0),  # carbonyl carbon
        ('C',    '[$([CX2](=C)=C)]',   0), # bonded to 3 things not hydrogen

        # Carbon SP3
        ('CT',   '[CX4H0]', 0), # single bonded to 4 things not hydrogen
        ('CH3',  '[C;H3;^3]',   0), # bonded to something via single bond and 3 Hs
        ('CH2',  '[C;^3;H2]',   0), # standard aliphatic C.
        ('CH1',  '*[C;H1](*)*', 1), # bonded to H and 3 things 

        # sp??? needs sorting 
        ('CH2',  '[CH2]',   0), # bonded to 2 hydrogens
        
        # Carbon fallback
        ('C', '[C,c]', 0),

        # Hydrogen
        ('HCH1', '[H][CH1]',    0),
        ('HCH2', '[H][C;H2^3]', 0),
        ('HCH3', '[H][CH3]',    0),
        ('HNC1', '[H][NX2;H1;^2]', 0), # H of N of N=C ? 
        ('HNC2', '[H][NX3;H2;^2]', 0), # H on a NC2 (NH1 and NH2 of ARG)
        ('HNC3', '[H][NX3;H3;^2]', 0), # guess - no examples
        ('HNT1', '[H][NX4;H1;^3]', 0),
        ('HNT1', '[H][NX3;H1;^3]', 0),
        ('HNT2', '[H][NX3;H2;^3]', 0), # H connected to type NT2
        ('HNT3', '[N^3;H3][H]', 1), # NH3+ 
        ('HNH2', '[H][NH2;^2]', 0), # NH2 is sp2
        ('HNH1', '[H][NX3;H1;^2]',    0),
        ('HCR6', '[H][cr6;H1]', 0),
        ('HCR5', '[H][cr5;H1]', 0), # connected to aromatic ring C with 1 H
        ('HNR5', '[H][nr5;H1]', 0), # connected to aromatic ring C with 1 H
        ('HNR5', '[H][Nr5;H1]', 0), # guess based on above, connected to aromatic N in a 5-ring
        ('HNR6', '[H][nr6;H1]', 0), # connected to aromatic 6-ring N with 1 H
        ('HNR6', '[H][NR6;H1]', 0), # guess based on above

        # HCR missing - no examples (and how is it different to HCR1?)
        ('HCR1', '[H]c',        0),
        ('HNH1', '[H][NH1]',    0),
        ('HOH1', '[H][OH1]',    0),
        ('HOH2', '[H][OH2][H]', (0,2)), # H of HOH - water
        ('H',    '[H]',         0),
        
        # Nitrogen, SP3

        ('NT1', '[NX4;H1;^3]',  0),
        ('NT1', '[NX3;H1;^3]',  0),
        ('NT2', '[NX3;H2;^3]',  0), # different to mon-lib!
        ('NT3', '[NX4;H3;^3]',  0),
        ('NT',  '[NX3;H0;^3]',  0),
        
        # Nitrogen, SP2
        ('NR66', 'c12aaaan1aaaa2', 5), # (second) 66 atom is an N.
        ('NR56', 'c12aaaan1aaa2',  5), # (second) 56 atom is an N.
        ('NR55', 'c12aaan1aaa2',   4), # (second) 55 atom is an N.
        ('NC2',  '[NX3;H2^2]', 0),     # N of sp2 NH2 (as in ARG).
        ('NH2',  '[NX3^2][CX3^2]=[N^2;X3+]', (0,2)), # amidinium (charged)... 
        ('NR15', '[nr5;X3;H1]',    0),
        ('NR5',  '[nr5;X3;H0]',    0),
        ('NR5',  '[NR;X3;H0;^2]',  0), # [NR5;X3;H0;^2] fails on 14C (also at daylight)
        ('NRD5', '[nr5;X2;H0]',    0), # guess from 071
        ('NRD5', 'C1(=O)[C,c][C,c]C=N1', 5), # N bonded to carbonyl C in a (non-percieved?) 5 ring, 0CE (a)
        ('NR16', '[nr6;H1]',    0),
        ('NRD6', 'a:[nr6;X2;H0]:a',  1), # aromatic N with no H, i.e. one double one single
        ('NR6',  '[nr6]',    0),
        ('NC1',  '[NX2;H1;^2]', 0),  # N of N=C ? 
        ('NH1',  '[NX3;H1;^2]', 0),
        ('NH2',  '[NX3;H2;^2]', 0),  # sp2, e.g. ND2 of an ASP
        ('NT',   '*n1~[o]~[o]1', 1), # guess from 16X dioxaziridine (bleugh)
        # (NT needs checking?)
        # NC2 no examples
        # NC3 no examples
        # NPA no examples
        # NPB no examples
        

        # Nitrogen SP1
        ('NS',   '[N^1]', 0),
        # NS1 no examples
        

        # fall-back nitrogen
        ('N',    '[N,n]',      0),  

        # Phosphorus
        ('P',    'P', 0),
        # Cl
        ('CL',   '[Cl]', 0),
        # F
        ('F',    '[F]',  0),
        # Br
        ('BR',    '[Br]',  0),

        # Sulfur
        ('SH1',  '[SH1]', 0),  # SG of CYS
        ('ST',   '[SX4]', 0), # tetrahedral (2 single bonds, 2 double)
        ('S1',   '[S]=*', 0),
        ('S2',   '[SX2,sX2]', 0),
        ('S3',   '[SX3,sX3]', 0),
        ('S',    '[S,s]', 0),

        ('SI1',  '[Si;X4]', 0), # tetragonal Si
        ('SI',   '[Si]',    0)  # Si any other

        ]

    full_list = smarts_list

    for item in smarts_by_element():
       full_list.append(item)

    for smarts_info in full_list:
        atom_type, smarts, match_atom_index = smarts_info
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
	    if False:
		print "SMARTS ", smarts
		print "  ", atom_type, ": ", matches
            for match in matches:
                set_atom_type(match, match_atom_index, mol, atom_type)
        else:
            # print "SMARTS ", smarts, " --- No hits  "
            pass

    # do we return success (everything has a type) or not?
    #
    for atom in mol.GetAtoms():
       try:
          atom_type = atom.GetProp('atom_type')
       except KeyError:
          is_aromatic = atom.GetIsAromatic()
          hybrid      = atom.GetHybridization()
          print "Error:: Missing type for atom \"", atom.GetProp('name'), "\" is_aromatic: ", is_aromatic, " hybridization: ", hybrid
          return False
    # we got to the end, good
    return True

# return True if mogul is not run or mogul exe is in place.
# return False if mogul is expected but not found.
def test_for_mogul():
   if (run_mogul):
      mogol_exe = which('mogul')
      if (mogol_exe == None):
         print "mogul not found in path"
         return False
      else:
         return True
   else:
      return True # OK, really

# this can throw a TypeError
#
def get_smiles_from_comp_id(comp_id):
   global smiles_dict
   if (not smiles_dict):
      read_smiles_tab('smiles.tab')
   return smiles_dict[comp_id]

# return a dictionary or False (if the file does not exist)
# (can this go inside get_smiles_from_comp_id?)
#
def read_smiles_tab(file_name):
    global smiles_dict
    try:
       smiles_dict = {}
       f = open(file_name)
       lines = f.readlines()
       for line in lines:
           bits = line.rstrip().rsplit()
           smiles_dict[bits[0]] = bits[2]
       return True
    except IOError as e:
       smiles_dict = True # we've tested for it
       return False

def get_smiles_from_file(file_name):
    f = open(file_name)
    smi = f.readline()
    return smi

def make_picture(mol, conf_id, comp_id, output_postfix):

   output_file_name = comp_id + "-" + output_postfix + '.png'
   make_picture_to_file(mol, conf_id, output_file_name)
      
def make_picture_to_file(mol, conf_id, output_file_name):

   try:
      from rdkit.Chem import Draw
      import Image
      state = Draw.MolToFile(mol, size=(300,300), fileName=output_file_name, confId=conf_id)
      # print 'INFO:: wrote PNG   "' + output_file_name + '"'

      # img = Draw.MolToImage(mol, fitImage=True, size=(900,900))
      # img2 = img.resize((300, 300), Image.ANTIALIAS)
      # img2.save(output_file_name + "resampled.png")
      
   except ImportError as e:
      print 'ImportError:', e
   except ValueError as e:
      print 'ValueError in make_picture():', e

def make_restraints_from_smiles(smiles_string, comp_id, compound_name, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name, quartet_planes, quartet_hydrogen_planes, use_mmff):

   if not test_for_mogul():
      # return False
      exit(1)
   m = Chem.MolFromSmiles(smiles_string)
   if compound_name:
       m.SetProp('_Name', compound_name)
   return make_restraints(m, comp_id, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name, quartet_planes, quartet_hydrogen_planes, use_mmff)

# return the molecule and return value from make_restraints
# 
def make_restraints_from_mdl(mol_file_name, comp_id, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name, quartet_planes, quartet_hydrogen_planes, use_mmff):

   if (not (test_for_mogul())):
      # return False, False
      exit(1)

   if not os.path.exists(mol_file_name):
      print "No such file:", mol_file_name
      exit(1)

   compound_name = '.'
   m = Chem.MolFromMolFile(mol_file_name)
   return m, make_restraints(m, comp_id, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name,
			     quartet_planes, quartet_hydrogen_planes, use_mmff)


# return a list of (mol, comp_id) pairs for every ligand in the cif
# file.  Often only one of course.
#
def make_restraints_from_mmcif_dict(cif_file_name_in, comp_id, mogul_dir, output_postfix,
				    quartet_planes, quartet_hydrogen_planes, use_mmff):

   if not test_for_mogul():
       return [(None, None)]

   if comp_id == "TRY_ALL_COMP_IDS":
      types = pysw.types_from_mmcif_dictionary(cif_file_name_in)
      l = []
      for type in types:
	    t_mol = make_restraints_from_mmcif_dict_single(cif_file_name_in, type, mogul_dir,
							   output_postfix,
							   quartet_planes,
							   quartet_hydrogen_planes, use_mmff)
	    l.append((t_mol, type))
      return l
   else:
       # just the one
       m = make_restraints_from_mmcif_dict_single(cif_file_name_in, comp_id, mogul_dir, output_postfix,
						  quartet_planes, quartet_hydrogen_planes, use_mmff)
       return [(m, comp_id)]

# return a mol, given a sensible comp_id.
#
# Return None on failure
#
def make_restraints_from_mmcif_dict_single(cif_file_name_in, comp_id, mogul_dir, output_postfix,
					   quartet_planes, quartet_hydrogen_planes, use_mmff):

   # print 'in make_restraints_from_mmcif_dict() comp_id is ', comp_id
   # print 'in make_restraints_from_mmcif_dict() cif_file_name_in is ', cif_file_name_in

   if not test_for_mogul():
       return [(None, None)]

   if comp_id == "TRY_ALL_COMP_IDS":
      types = pysw.types_from_mmcif_dictionary(cif_file_name_in)
      l = []
      for type in types:
	    t_mol = make_restraints_from_mmcif_dict(cif_file_name_in, type, mogul_dir, output_postfix,
						    quartet_planes, quartet_hydrogen_planes)
	    l.append((t_mol, type))
      return l
   else:
      file_name_stub           = comp_id + '-' + output_postfix 
      pdb_out_file_name        = file_name_stub + '.pdb'
      cif_restraints_file_name = file_name_stub + '.cif'
      m = pyrogen_boost.rdkit_mol_chem_comp_pdbx(cif_file_name_in, comp_id)

      if False:  # debugging
	  for atom in m.GetAtoms():
	      try:
		  name   = atom.GetProp('name')
		  chir   = atom.GetProp('_CIPCode')
		  print ' atom', atom, 'name', name, 'chir', chir
	      except KeyError as e:
		  print 'pyrogen.py:: atom', atom, " with name ", name, ' has no _CIPCode property'
		  pass


      # maybe user didn't select the correct comp_id for the given dictionary mmcif
      if m.GetNumAtoms() == 0:
	  print 'No atoms for comp_id', comp_id
	  return False
      else :

	  name = ''
	  try:
	      name = m.GetProp('_Name')
	  except KeyError:
	      print 'caught KeyError in make_restraints_from_mmcif_dict() trying GetProp _Name'

	  return make_restraints(m, comp_id, mogul_dir, file_name_stub,
				 pdb_out_file_name, cif_restraints_file_name,
				 quartet_planes, quartet_hydrogen_planes, use_mmff)



def n_hydrogens(mol):
    n_H = 0
    for atom in mol.GetAtoms():
	if atom.GetAtomicNum() == 1:
	    n_H += 1
    return n_H


# return sane_H_mol
# 
def make_restraints(m, comp_id, mogul_dir, file_name_stub, pdb_out_file_name, mmcif_dict_name,
                    quartet_planes, quartet_hydrogen_planes, use_mmff):

   # pH-dependent protonation or deprotonation
   #
   do_hydrogen_atoms_shift = True

   try:
      compound_name = m.GetProp('_Name');
   except KeyError:
      # this happens all the time when we start from a SMILES, users don't need to see it.
      # print 'caught key error in trying to get _Name in make_restraints() for m'
      compound_name = '.'
   except AttributeError as e:
      # Do we need to see this? Perhaps make_restraints() needs to return a status.
      # print 'AttributeError: problem with molecule in make_restraints()', e, ' on object:', m
      return

   m_H = m
   if n_hydrogens(m) == 0:
       m_H = AllChem.AddHs(m)

   if do_hydrogen_atoms_shift:
      # simple sane pH H-exchanges
      sane_H_mol = pyrogen_boost.hydrogen_transformations(m_H)
      print >>file('sane_H.mol','w+'),Chem.MolToMolBlock(sane_H_mol)
   else:
      sane_H_mol = m_H
 
   AllChem.EmbedMolecule(sane_H_mol)

   if use_mmff:
      AllChem.MMFFOptimizeMolecule(sane_H_mol)
      
      if False:  # debugging output
         ba = pyrogen_boost.mmff_bonds_and_angles(sane_H_mol) # uses _forcefield_ of the molecule
	 n_bonds = ba.bonds_size()
	 if n_bonds > 0:
	    for i_bond in range(n_bonds):
	       bond = ba.get_bond(i_bond)
	       print bond.get_idx_1(), bond.get_idx_2(), bond.get_type(), \
		     bond.get_resting_bond_length(), bond.get_sigma()
         n_angles = ba.angles_size()
	 if n_angles > 0:
	     for i_angle in range(n_angles):
		 angle = ba.get_angle(i_angle)
		 print angle.get_idx_1(), angle.get_idx_2(), angle.get_idx_3(), \
		       angle.get_resting_angle(), angle.get_sigma()
         
   else:
      AllChem.UFFOptimizeMolecule(sane_H_mol)

   atom_names = add_atom_names(sane_H_mol)

   sane_H_mol.SetProp('comp_id', comp_id)
   sane_H_mol.SetProp('name', compound_name)
   all_set = set_atom_types(sane_H_mol)

   if (all_set != True):
      return False
   else:

      sd_local = file_name_stub + ".sdf"
      sdf_file_name       = os.path.join(mogul_dir, file_name_stub + '-mogul.sdf')
      mogul_ins_file_name = os.path.join(mogul_dir, file_name_stub + '-mogul.ins')
      mogul_out_file_name = os.path.join(mogul_dir, file_name_stub + '-mogul.out')
      Chem.AllChem.ComputeGasteigerCharges(sane_H_mol)

      moguled_mol = pyrogen_boost.mogulify(sane_H_mol) # Nitro bond orders (and other things?)
      if not os.path.isdir(mogul_dir):
	  checked_mkdir(mogul_dir)
	  if os.path.isdir(mogul_dir):
	      mb = Chem.MolToMolBlock(moguled_mol)
	      print >> file(sdf_file_name,'w'), mb
      else:
	  mb = Chem.MolToMolBlock(moguled_mol)
	  print >> file(sdf_file_name,'w'), mb

      
      bor = make_restraints_for_bond_orders(sane_H_mol)

      # print out the set types:
      print_atom_props = False
      if print_atom_props:
	  print '--- Atom Props ---'
      for atom in sane_H_mol.GetAtoms():
         charge = atom.GetProp('_GasteigerCharge') # string?
         name   = atom.GetProp('name')
         try:
            atom_type   = atom.GetProp('atom_type')
            is_aromatic = atom.GetIsAromatic()
            hybrid      = atom.GetHybridization()
            f_charge    = float(charge)
	    if print_atom_props:
		print "  atom: %s %s type: %s arom: %s hybrid: %s charge: %6.3f" % (name, atom.GetSymbol(),
										    atom_type.ljust(4),
										    str(is_aromatic).ljust(5),
										    str(hybrid).rjust(3),
										    f_charge)
         except KeyError:
            print "miss", name, atom.GetSymbol(), charge

      #
      replace_with_mmff_b_a_restraints = False
      if use_mmff:
	  replace_with_mmff_b_a_restraints = True

      # execute_mogul() tests if mogul is executable
      #
      mogul_state = execute_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name)
      if mogul_state:

         # Here we need to think about matching to reference
         # dictionary of amino acids (for standard atom names).
         # That function takes a dictionary and a mmdb::Residue.
         # How does that fit in here?
         #
         restraints = pysw.mogul_out_to_mmcif_dict_by_mol(mogul_out_file_name, comp_id,
                                                          compound_name, sane_H_mol, bor,
							  mmcif_dict_name,
                                                          quartet_planes,
							  quartet_hydrogen_planes,
							  replace_with_mmff_b_a_restraints)
         pysw.regularize_and_write_pdb(sane_H_mol, restraints, comp_id, pdb_out_file_name)

      else:

	  restraints = pysw.mmcif_dict_from_mol(comp_id, compound_name, sane_H_mol,
						mmcif_dict_name,
						quartet_planes, quartet_hydrogen_planes,
						replace_with_mmff_b_a_restraints)
	  if restraints == None:
	      print "No restraints"
	      return True # hacked in value
	  
	  pysw.write_pdb_from_mol(sane_H_mol, comp_id, pdb_out_file_name)
	  
      return sane_H_mol


def score_and_print_tautomers(mol, comp_id, output_postfix, do_drawings):

    results = tautomer.enumerate_tautomers(mol)
    for i in range(len(results)):
	m = results[i]
	s = Chem.MolToSmiles(m)
	print "comp_id :", comp_id, ": SMILES", s, 'score:', tautomer.tautomer_score(m)
	if do_drawings:
	    file_name = comp_id + '-tautomer-' + str(i)
	    file_name += '-' + options.output_postfix + '.png'
	    n = m.GetNumConformers()
	    conf_id = 0
	    if n == 0:
		conf_id = AllChem.Compute2DCoords(m)
	    conf = m.GetConformer(conf_id)
	    if conf.Is3D():
		mol_for_drawing = Chem.RemoveHs(m, implicitOnly=False)
		conf2D_id = AllChem.Compute2DCoords(mol_for_drawing)
		make_picture_to_file(mol_for_drawing, conf2D_id, file_name)
	    else:
		make_picture_to_file(m, -1, file_name)
   

if __name__ == "__main__":

    def checked_mkdir(dirname):
        if not os.path.exists(dirname):
	    os.makedirs(dirname)
	else:
	   if os.path.isdir(dirname):
	       pass # this happens most of the time, I imagine
	   else:
	       print 'Stop:: File', dirname, 'exists but is not a directory'

    def smiles_from(smi_raw):
       extension = os.path.splitext(smi_raw)[1]
       smiles_string = ''
       if extension == '.smi' or extension == '.smiles':
          smiles_string = get_smiles_from_file(smi_raw)
       else:
         smiles_string = smi_raw
       return smiles_string

    parser = OptionParser(usage='pyrogen [options] file-or-SMILES'+
                          '\n       if file-or-SMILES has extension ".smi" or ".smiles" ' +
                          'then it is treated as a file')
    parser.add_option("-c", "--mmcif", dest="mmcif_file",
		      help="Make restraints from input mmcif FILE", metavar="FILE")
    parser.add_option("-m", "--mol", dest="sdf_file",
		      help="Make restraints from input sdf/mol FILE", metavar="FILE")
    parser.add_option("-r", "--residue-type", dest="comp_id", default='default',
		      help="Create restraints for this type. Default is LIG")
    parser.add_option("-4", "--quartet-planes", dest="quartet_planes",
		      default=False,
		      help="Use 4-atom plane restraints,\n                    " +
                      "forces --quartet-hydrogens", action="store_true")
    parser.add_option("-H", "--quartet-hydrogens", dest="quartet_hydrogen_planes",
		      default=False,
		      help="Use 4-atom hydrogen plane restraints",
                      action="store_true")
    parser.add_option("-n", "--no-mogul", dest="use_mogul",
		      default=True, action="store_false",
                      help='Don\'t run CSD Mogul to update bond and angle restraints')
    parser.add_option("-N", '--name', dest='compound_name', default=False,
		      help='Compound name')
    parser.add_option("-t", "--tautomers", dest="show_tautomers",
		      default=False, action="store_true",
                      help='Show SMILES for tautomers, don\'t generate restraints')
    parser.add_option("-d", '--directory', dest='mogul_dir',
                      help='Directory into which the tmp files (e.g. for mogul) are written',
                      default='pyrogen-mogul')
    parser.add_option('-o', '--output-postfix', default='pyrogen',
                      dest='output_postfix',
                      help='string to add to output file names, default is "pyrogen"')
    parser.add_option('-p', '--picture', dest='drawing',
                      help='Additionally output a chemical diagram PNG',
                      action='store_true', default=False)
    parser.add_option('-v', '--version', dest='show_version', default=False,
                      action='store_true', help='Print version information')
    parser.add_option('-M', '--testing', dest='use_mmcif', default=False,
                      action='store_true', help='Testing function')
    parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="print less messages")

    (options, args) = parser.parse_args()
    # print 'DEBUG:: options:', options

    if options.show_version:
       print 'pyrogen-' + pyrogen_version, "revision", coot_svn_repo_revision.revision_number()

    comp_id = options.comp_id
    if options.comp_id == 'default':
	comp_id = 'LIG'
    if options.mmcif_file != None:
	if options.comp_id == 'default':
	    comp_id = 'TRY_ALL_COMP_IDS'

    pdb_out_file_name        = comp_id + '-' + options.output_postfix + '.pdb'
    cif_restraints_file_name = comp_id + '-' + options.output_postfix + '.cif'
    file_name_stub           = comp_id + '-' + options.output_postfix

    # this is a bit ugly, perhaps.  this value is inspected inside
    # the following functions
    #
    if options.use_mogul == False:
       run_mogul = False
    if run_mogul:
       if len(options.mogul_dir) > 0:
          if options.mogul_dir[0] == '-':
             print 'Stop:: you probably didn\'t mean that you wanted',options.mogul_dir, 'as your tmp directory.'
             exit(1)
          checked_mkdir(options.mogul_dir)


    if options.show_tautomers:
	mol = False
	if len(args) > 0:
	    smi_raw = args[0]
	    smiles = smiles_from(smi_raw)
	    mol = Chem.MolFromSmiles(smiles)
	else:
	    if options.sdf_file != None:
		mol = Chem.MolFromMolFile(options.sdf_file)
	    else:
		if options.mmcif_file != None:
		    types = pysw.types_from_mmcif_dictionary(options.mmcif_file)
		    print '-- tautomer mode: mmcif file types:', types
		    for type in types:
			mol_local = pyrogen_boost.rdkit_mol_chem_comp_pdbx(options.mmcif_file, type)
			score_and_print_tautomers(mol_local, type, options.output_postfix, options.drawing)

	if mol:
	    score_and_print_tautomers(mol, comp_id, options.output_postfix, options.drawing)

    else:

	if options.mmcif_file:
	    mol_pairs = make_restraints_from_mmcif_dict(options.mmcif_file, comp_id, options.mogul_dir,
							options.output_postfix,
							options.quartet_planes,
							options.quartet_hydrogen_planes,
							options.use_mmcif)

	    # this needs to be in a try block, I suppose, for example if the mmcif file
	    # does not exist.
	    
	    for mol_info in mol_pairs:
		(mol, comp_id) = mol_info
		if not mol:
		    print 'No molecule'
		else:
		    # Happy path
		    if options.drawing:
			# make_picture() by default draws the first conformer in the given molecule.
			# For mol, that is a 3D conformer.  We want to draw a nice 2D diagram
			#
			mol_for_drawing = Chem.RemoveHs(mol, implicitOnly=False)
			conf2D_id = AllChem.Compute2DCoords(mol_for_drawing)
			make_picture(mol_for_drawing, conf2D_id, comp_id, options.output_postfix)

	else:

	    if options.sdf_file != None:
		(mol, results) = make_restraints_from_mdl(options.sdf_file, comp_id,
							  options.mogul_dir, file_name_stub,
							  pdb_out_file_name, cif_restraints_file_name,
							  options.quartet_planes,
							  options.quartet_hydrogen_planes,
							  options.use_mmcif)
		if options.drawing:
		    make_picture(mol, -1, comp_id, options.output_postfix)

	    else:

		if len(args) > 0:
		    smi_raw = args[0]
		    smiles = smiles_from(smi_raw)
		    status = make_restraints_from_smiles(smiles, comp_id, options.compound_name,
							 options.mogul_dir, file_name_stub,
							 pdb_out_file_name,
							 cif_restraints_file_name,
							 options.quartet_planes,
							 options.quartet_hydrogen_planes,
							 options.use_mmcif)
		    if options.drawing:
			mol = Chem.MolFromSmiles(smiles)
			make_picture(mol, -1, comp_id, options.output_postfix)
