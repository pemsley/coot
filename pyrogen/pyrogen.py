# Copyright 2012 by the University of Oxford
# Copyright 2014, 2015 by Medical Research Council
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

"""
pyrogen contains tools for converting molecules to mmcif restraints dictionaries,
utilities for retrival, extraction and depction.
"""

import sys
import os
import copy
from subprocess import call
from rdkit import Chem
from rdkit.Chem import AllChem

import coot_git

# Hello from 20170701, library resolution problems?
# $ otool -L .libs/_pyrogen_swig.so -> @rpath substition... how does that work?
import pyrogen_swig as pysw # relies on above rdkit.Chem import AllChem (but ideally should not)
import pyrogen_boost        # ditto
import atom_types

from optparse import OptionParser

import tautomer

import urllib

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
          state = call(['mogul', '-ins', mogul_ins_file_name])
	  return (state == 0)
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
    r = f.read()
    f.close()
    return r


# return False or a file_name
#
# downloaded file is put in the CCD directory as CCD/x/xyz.cif
#
def get_pdbe_cif_for_comp_id(comp_id):

   try:
      url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif/' + comp_id + '.cif'
      CCD_dir = 'CCD'
      # file_name = "PDBe-" + comp_id + ".cif"
      first_char = comp_id[0]
      try:
         sub_dir = os.path.join(CCD_dir, first_char)
         file_name = os.path.join(sub_dir, comp_id + ".cif")
         if not os.path.isdir(CCD_dir):
            os.mkdir(CCD_dir)
         if not os.path.isdir(sub_dir):
            os.mkdir(sub_dir)
         if os.path.isfile(file_name):
            return file_name
         else:
            url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif/' + comp_id + '.cif'
            status = urllib.urlretrieve(url, file_name)
            print('urllib.urllib returned with status', status)
            return file_name

      except OSError as e:
         print e
         print "Failed: Can't ftp from", url, "and write file", file_name

   except IOError as e:
      print e
      print "Failed: Can't ftp from", url, "and write file", file_name


def MolFromFetchedCode(code):
   f = get_pdbe_cif_for_comp_id(code)
   m = pyrogen_boost.MolFromPDBXr(f, code)
   m.Compute2DCoords()
   return m
        
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


# return True if mogul is not run or mogul exe is in place.
# return False if mogul is expected but not found.
def test_for_mogul():
   if run_mogul:
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

# return a pair, the smiles string and the molecule name (which might be blank)
#
def get_smiles_from_file(file_name):
    if not os.path.exists(file_name):
	return False,False
    else:
	f = open(file_name)
	smi_line = f.readline()
	parts = smi_line.split()
	return parts[0], ' '.join(parts[1:])

def make_picture(mol, conf_id, comp_id, output_postfix):

   output_file_name = comp_id + "-" + output_postfix + '.png'
   make_picture_to_file(mol, conf_id, output_file_name)
      
def make_picture_to_file(mol, conf_id, output_file_name):

   try:
      from rdkit.Chem import Draw
      # The import of Image may change depending on how it was provided.
      # What about pillow? Hmm. Not sure of the details.
      # import Image
      from PIL import Image
      state = Draw.MolToFile(mol, size=(300,300), fileName=output_file_name, confId=conf_id)
      # print 'INFO:: wrote PNG   "' + output_file_name + '"'

      # img = Draw.MolToImage(mol, fitImage=True, size=(900,900))
      # img2 = img.resize((300, 300), Image.ANTIALIAS)
      # img2.save(output_file_name + "resampled.png")
      
   except ImportError as e:
      print 'ImportError:', e
   except ValueError as e:
      print 'ValueError in make_picture():', e

def make_restraints_from_smiles(smiles_string, comp_id, compound_name, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name, quartet_planes, quartet_hydrogen_planes, use_mmff, match_atom_names_to_dict_flag, comp_id_list_for_names_match, dict_file_for_names_match):

   if not test_for_mogul():
      # return False
      exit(1)
   m = Chem.MolFromSmiles(smiles_string)
   if compound_name:
       m.SetProp('_Name', compound_name)
   return make_restraints(m, comp_id, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name, quartet_planes, quartet_hydrogen_planes, use_mmff, match_atom_names_to_dict_flag, comp_id_list_for_names_match, dict_file_for_names_match)

# return the molecule and return value from make_restraints
# 
def make_restraints_from_mdl(mol_file_name, comp_id, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name, quartet_planes, quartet_hydrogen_planes, use_mmff, match_atom_names_to_dict_flag, comp_id_list_for_names_match, dict_files_for_names_match):

   if (not (test_for_mogul())):
      # return False, False
      exit(1)

   if not os.path.exists(mol_file_name):
      print "No such file:", mol_file_name
      exit(1)

   compound_name = '.'
   m = Chem.MolFromMolFile(mol_file_name)
   return m, make_restraints(m, comp_id, mogul_dir, name_stub, pdb_out_file_name, mmcif_dict_name,
			     quartet_planes, quartet_hydrogen_planes, use_mmff,
                             match_atom_names_to_dict_flag, comp_id_list_for_names_match,
                             dict_files_for_names_match)


# return a list of (mol, comp_id) pairs for every ligand in the cif
# file.  Often only one of course.
#
def make_restraints_from_mmcif_dict(cif_file_name_in, comp_id, mogul_dir,
                                    output_dir, output_postfix,
				    quartet_planes, quartet_hydrogen_planes, use_mmff,
                                    pdb_out_file_name, mmcif_restraints_out_file_name):

   if not test_for_mogul():
       return [(None, None)]

   if comp_id == "TRY_ALL_COMP_IDS":
      types = pysw.types_from_mmcif_dictionary(cif_file_name_in)
      l = []
      for r_type in types:

         file_name_stub = r_type + "-" + output_postfix
         if options.output_dir != ".":
            file_name_stub = os.path.join(options.output_dir, file_name_stub)
            
         pdb_out_file_name_local              = file_name_stub + ".pdb"
         mmcif_restraints_out_file_name_local = file_name_stub + ".cif"
         #
         t_mol = make_restraints_from_mmcif_dict_single(cif_file_name_in, r_type, mogul_dir,
                                                        output_postfix,
                                                        quartet_planes,
                                                        quartet_hydrogen_planes, use_mmff,
                                                        pdb_out_file_name_local,
                                                        mmcif_restraints_out_file_name_local)
         l.append((t_mol, r_type))
      return l
   else:
       # just the one
       m = make_restraints_from_mmcif_dict_single(cif_file_name_in, comp_id, mogul_dir, output_postfix,
						  quartet_planes, quartet_hydrogen_planes, use_mmff,
                                                  pdb_out_file_name, mmcif_restraints_out_file_name)
       return [(m, comp_id)]

# return a mol, given a sensible comp_id.
#
# Return None on failure
#
def make_restraints_from_mmcif_dict_single(cif_file_name_in, comp_id, mogul_dir, output_postfix,
					   quartet_planes, quartet_hydrogen_planes, use_mmff,
                                           pdb_out_file_name, mmcif_restraints_out_file_name):

   # print 'in make_restraints_from_mmcif_dict_single() comp_id is ', comp_id
   # print 'in make_restraints_from_mmcif_dict_single() cif_file_name_in is ', cif_file_name_in

   if not test_for_mogul():
       return [(None, None)]
    
   mogul_file_name_stub = comp_id + '-' + output_postfix # file component of files within mogul_dir

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
         print 'caught KeyError in make_restraints_from_mmcif_dict_single() trying GetProp _Name'

      return make_restraints(m, comp_id, mogul_dir, mogul_file_name_stub,
                             pdb_out_file_name, mmcif_restraints_out_file_name,
                             quartet_planes, quartet_hydrogen_planes, use_mmff, False, False, False)


def n_hydrogens(mol):
    n_H = 0
    for atom in mol.GetAtoms():
	if atom.GetAtomicNum() == 1:
	    n_H += 1
    return n_H


# return sane_H_mol
# 
def make_restraints(m, comp_id, mogul_dir, mogul_file_name_stub, pdb_out_file_name, mmcif_dict_name,
                    quartet_planes, quartet_hydrogen_planes, use_mmff,
                    match_atom_names_to_dict_flag,
                    comp_id_list_for_names_match,
                    dict_files_for_names_match):

   # test here (or in calling functions) if m is sane (i.e. is an rdkit molecule)

   if not isinstance(m, Chem.rdchem.Mol):
      print 'ERROR:: not a molecule'
      return False

   n_attempts = 20 * m.GetNumAtoms() # default is 10 * number of atoms.

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
      # print >>file('sane_H.mol','w+'),Chem.MolToMolBlock(sane_H_mol)
   else:
      sane_H_mol = m_H

   # This makes UFF types, which can fail sometimes.
   # conf_id = AllChem.EmbedMolecule(sane_H_mol, AllChem.ETKDG(), maxAttempts=n_attempts)
   conf_id = AllChem.EmbedMolecule(sane_H_mol, AllChem.ETKDG())

   if use_mmff:
      AllChem.MMFFOptimizeMolecule(sane_H_mol, confId=conf_id)
      
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
      AllChem.UFFOptimizeMolecule(sane_H_mol, confId=conf_id)

   # AllChem.UFFOptimizeMolecule(sane_H_mol)

   atom_names = add_atom_names(sane_H_mol) 
   all_set = atom_types.set_atom_types(sane_H_mol)  # has deloc bonds now, potentially

   # debug sane_H_mol
   if True:
      molblock = Chem.MolToMolBlock(sane_H_mol)
      print >> file("sane_H_mol.mol",'w'), molblock


   if (all_set != True):
      return False
   else:

      sane_H_mol.SetProp('comp_id', comp_id)
      sane_H_mol.SetProp('name', compound_name)
      
      sd_local = mogul_file_name_stub + ".sdf"
      sdf_file_name       = os.path.join(mogul_dir, mogul_file_name_stub + '-mogul.sdf')
      mogul_ins_file_name = os.path.join(mogul_dir, mogul_file_name_stub + '-mogul.ins')
      mogul_out_file_name = os.path.join(mogul_dir, mogul_file_name_stub + '-mogul.out')
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
            atom_type   = atom.GetProp('type_energy')
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
							  mmcif_dict_name, # not used
                                                          quartet_planes,
							  quartet_hydrogen_planes,
							  replace_with_mmff_b_a_restraints)


         # match_atom_names_to_dict_flag, comp_id_list_for_names_match, dict_file_for_names_match
         if match_atom_names_to_dict_flag:

             restraints = atom_match_dictionary(restraints, sane_H_mol,
                                                comp_id_list_for_names_match,
                                                dict_files_for_names_match)

         pysw.write_restraints(restraints, mmcif_dict_name)
         pysw.regularize_and_write_pdb(sane_H_mol, restraints, comp_id, pdb_out_file_name)

      else:

	  # mogul failed or was not in the path:
	  
	  if run_mogul == False:

	      # ... but that's OK if we told pyrogen to run without mogul

              # sane_H_mol:
              # print >>file('debug_sane_H.mol','w+'),Chem.MolToMolBlock(sane_H_mol)

	      restraints = pysw.mmcif_dict_from_mol(comp_id, compound_name, sane_H_mol,
						    mmcif_dict_name,
						    quartet_planes, quartet_hydrogen_planes,
						    replace_with_mmff_b_a_restraints)
              
	      if restraints == None:
		  print "No restraints"
		  return True # hacked in value

              if match_atom_names_to_dict_flag:
                  
                  restraints = atom_match_dictionary(restraints, sane_H_mol,
                                                     comp_id_list_for_names_match,
                                                     dict_files_for_names_match)
                  pysw.write_restraints(restraints, mmcif_dict_name)
                  
	      pysw.write_pdb_from_mol(sane_H_mol, comp_id, pdb_out_file_name)

	  else:
	      # ... but not if we wanted to use mogul.
	      # (We get here if there is a licence error for mogul)
	      exit(1)
	  
      return sane_H_mol

def atom_match_dictionary(restraints, sane_H_mol, comp_id_list_for_names_match, dict_files_for_names_match):

    template_comp_ids = ['CYS', 'ASP', 'GLU',        'HIS', 'ILE', 'LYS', 'LEU', 'MET',
                         'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR',
                         'G',   'A',     'C',   'U',    'GLC', 'MAN']

    if isinstance(comp_id_list_for_names_match, basestring):
        template_comp_ids = comp_id_list_for_names_match.split(',')

    template_cif_dict_files_names = []
    if isinstance(dict_files_for_names_match, basestring):
        template_cif_dict_files_names = dict_files_for_names_match.split(',')
        # don't use my set of comp_ids then
        template_comp_ids = []

    success,new_restraints,at_name_list = pysw.match_restraints_to_dictionaries(restraints,
                                                                                template_comp_ids,
                                                                                template_cif_dict_files_names)
    if success:
        for at_name in at_name_list:
            # print at_name
            pass
        n = len(sane_H_mol.GetAtoms())
        if len(restraints['_chem_comp_atom']) == n:
            restraints = new_restraints
            for iat in range(n):
                name = sane_H_mol.GetAtomWithIdx(iat).GetProp('name')
		if name != restraints['_chem_comp_atom'][iat][0]:
		    # print "   changing name from", name, "to", restraints['_chem_comp_atom'][iat][0]
		    sane_H_mol.GetAtomWithIdx(iat).SetProp('name', restraints['_chem_comp_atom'][iat][0]);

    return restraints



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

    def smiles_and_name_from(smi_raw):
       extension = os.path.splitext(smi_raw)[1]
       smiles_string = ''
       name=''
       if extension == '.smi' or extension == '.smiles':
	   if not os.path.exists(smi_raw):
	       print "File not found:", smi_raw
	       exit(1)
	   else: 
	       smiles_string,name = get_smiles_from_file(smi_raw)
       else:
         smiles_string = smi_raw
       return smiles_string,name

    parser = OptionParser(usage='pyrogen [options] file-or-SMILES'+
                          '\n       if file-or-SMILES has extension ".smi" or ".smiles" ' +
                          'then it is treated as a file')
    parser.add_option("-c", "--mmcif", dest="mmcif_file_name",
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
    parser.add_option('-S', '--smiles', dest="show_smiles",
                      default=False, action="store_true", help="Write the SMILES for the input molecule")
    parser.add_option("-t", "--tautomers", dest="show_tautomers",
		      default=False, action="store_true",
                      help='Show SMILES for tautomers, don\'t generate restraints')
    parser.add_option("-T", '--tmp-directory', dest='mogul_dir',
                      help='Directory into which the tmp files (e.g. for mogul) are written',
                      default='pyrogen-mogul')
    parser.add_option("-d", '--directory', dest='output_dir',
                      help='Directory into which the output files (e.g. mmCIF and PDB) are written',
                      default='.')
    parser.add_option('-o', '--output-postfix', default='pyrogen',
                      dest='output_postfix',
                      help='string to add to output file names, default is "pyrogen"')
    parser.add_option('-p', '--picture', dest='drawing',
                      help='Additionally output a chemical diagram PNG',
                      action='store_true', default=False)
    parser.add_option('-v', '--version', dest='show_version', default=False,
                      action='store_true', help='Print version information')
    parser.add_option('-M', '--MMFF', dest='use_mmff', default=False,
                      action='store_true', help='Use MMFF fallbacks for bonds and angles')
    parser.add_option('-a', '--no-match-vs-reference-dictionaries', default=False,
                      action='store_true', dest='no_match_names_flag',
                      help="Don't match atom names vs. dictionary molecules (default False)")
    parser.add_option('-R', '--reference-dictionary-files', dest='dict_files_for_names_match',
                      help='Try to match the atom names of the output molecule '+
                      'to this dictionary in these files (comma-separated list)', default=False)
    parser.add_option('-C', '--reference-dictionary-comp-ids', dest='comp_id_list_for_names_match',
                      help='Try to match the atom names of the output molecule to these comp-ids' +
                      ' (comma-separated list)',
                      default=False)
    parser.add_option('-w', '--wwPDB', default=False, dest="wwPDB", action="store_true",
                      help='Fetch the wwPDB ligand definition and use that')
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="print less messages")

    (options, args) = parser.parse_args()
    # print 'DEBUG:: options:', options
    # print 'DEBUG:: args:', args

    if len(args) == 0:
       print("Usage: pyrogen --help")

    if options.show_version:
       print 'pyrogen-' + pyrogen_version, "revision", coot_git.revision_count()

    comp_id = options.comp_id
    if options.comp_id == 'default':
	comp_id = 'LIG'
    if options.mmcif_file_name != None:
	if options.comp_id == 'default':
	    comp_id = 'TRY_ALL_COMP_IDS'

    file_name_stub           = comp_id + '-' + options.output_postfix
    
    if options.output_dir != ".":
       file_name_stub = os.path.join(options.output_dir, file_name_stub)

    pdb_out_file_name              = file_name_stub + '.pdb'
    mmcif_restraints_out_file_name = file_name_stub + '.cif'

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

    if options.show_tautomers or options.show_smiles:

        # ------------------------ Tautomers and SMILES  ---------------------------------------------

	mol = False
	if len(args) > 0:
	    smi_raw = args[0]
	    smiles,compound_name = smiles_and_name_from(smi_raw)
	    mol = Chem.MolFromSmiles(smiles)
	else:
	    if options.sdf_file != None:
		mol = Chem.MolFromMolFile(options.sdf_file)
	    else:
		if options.mmcif_file_name != None:
		    types = pysw.types_from_mmcif_dictionary(options.mmcif_file_name)
		    print '-- tautomer mode: mmcif file types:', types
		    for type in types:
			mol_local = pyrogen_boost.rdkit_mol_chem_comp_pdbx(options.mmcif_file_name, type)
			score_and_print_tautomers(mol_local, type, options.output_postfix, options.drawing)

	if mol:
           if options.show_tautomers:
              score_and_print_tautomers(mol, comp_id, options.output_postfix, options.drawing)
           if options.show_smiles:
              s = Chem.MolToSmiles(mol);
              print s

    else:

        # ------------------------ dict-build-mode ---------------------------------------------------

        mmcif_file_name = options.mmcif_file_name
        # shall we go get the dictionary?
        if options.wwPDB:
           mmcif_file_name = get_pdbe_cif_for_comp_id(comp_id)
           if os.path.isfile(mmcif_file_name):
              pass # good
           else:
              print "Missing downloaded file for comp-id:",  comp_id
              exit(2)

        # JED mode for hydrogen planes
        #
        quartet_hydrogen_planes = options.quartet_hydrogen_planes
        if options.quartet_planes:
            quartet_hydrogen_planes = True

        match_names_flag = True
        if options.no_match_names_flag:
            match_names_flag = False

	if mmcif_file_name:
	    mol_pairs = make_restraints_from_mmcif_dict(mmcif_file_name,
                                                        comp_id,
                                                        options.mogul_dir,
                                                        options.output_dir,
							options.output_postfix,
							options.quartet_planes,
							quartet_hydrogen_planes,
							options.use_mmff,
                                                        pdb_out_file_name,
                                                        mmcif_restraints_out_file_name)

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
                        conf = mol.GetConformer(conf2D_id)
                        Chem.WedgeMolBonds(mol_for_drawing, conf)
			make_picture(mol_for_drawing, conf2D_id, comp_id, options.output_postfix)

	else:

	    if options.sdf_file != None:
		(mol, results) = make_restraints_from_mdl(options.sdf_file, comp_id,
							  options.mogul_dir, file_name_stub,
							  pdb_out_file_name,
                                                          mmcif_restraints_out_file_name,
							  options.quartet_planes,
							  quartet_hydrogen_planes,
							  options.use_mmff,
                                                          match_names_flag,
                                                          options.comp_id_list_for_names_match,
                                                          options.dict_files_for_names_match)
		if options.drawing:
		    make_picture(mol, -1, comp_id, options.output_postfix)

	    else:

		if len(args) > 0:
		    smi_raw = args[0]
		    smiles,compound_name_from_file = smiles_and_name_from(smi_raw)
                    compound_name=False
                    if len(compound_name_from_file) > 0:
                       compound_name = compound_name_from_file
                    if isinstance(options.compound_name, basestring):
                       compound_name = options.compound_name
		    status = make_restraints_from_smiles(smiles, comp_id, compound_name,
							 options.mogul_dir, file_name_stub,
							 pdb_out_file_name,
							 mmcif_restraints_out_file_name,
							 options.quartet_planes,
							 quartet_hydrogen_planes,
							 options.use_mmff,
                                                         match_names_flag,
                                                         options.comp_id_list_for_names_match,
                                                         options.dict_files_for_names_match)
		    if options.drawing:
			mol = Chem.MolFromSmiles(smiles)
			make_picture(mol, -1, comp_id, options.output_postfix)
