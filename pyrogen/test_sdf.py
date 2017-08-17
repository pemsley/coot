
# coding: utf-8
from __future__ import print_function
import pyrogen
import pyrogen_boost
import atom_types
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Draw
from rdkit.Chem.Draw import IPythonConsole

import pyrogen_swig as pysw

def write_mols(suppl):
    for idx,mol in enumerate(suppl):
        fn = 'pf-'+str(idx)+'.sdf'
        # print('writing {}'.format(fn))
        print(Chem.MolToMolBlock(mol), file=file(fn, 'w'))

def get_amber_types(m, mol_idx, ref_types=None, check_2_chars_only=True, make_dictionaries=False):
    Chem.AllChem.ComputeGasteigerCharges(m)
    names = pyrogen.add_atom_names(m)
    atom_types.set_atom_types(m)
    atom_types.set_parmfrosst_atom_types(m)
    aa_types=[]
    cif_file_name = "test-pyrogen-"+str(mol_idx)+'.cif'
    if make_dictionaries:
       pysw.mmcif_dict_from_mol("XYZ", m.GetProp('_Name'), m, cif_file_name, False, False, True)
    for idx,atom in enumerate(m.GetAtoms()):
        try:
            # in pathological cases, pf_atom_type is not set
            coot_amber_type = atom.GetProp('pf_atom_type')
            # if we were passed ref_types, compare the atom type
            name = atom.GetProp('name')
            ref_type=ref_types[idx]
            if ref_type == "Nstar":
               ref_type = 'N*'
            check_coot_amber_type = coot_amber_type
            if check_2_chars_only:
                check_coot_amber_type = coot_amber_type[:2]
		# C -> Ca or Cb or Cc etc is a special case
		if coot_amber_type == "Ca":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cb":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cc":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cd":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Ce":
		    check_coot_amber_type = "C"
		if coot_amber_type == "Cf":
		    check_coot_amber_type = "C"
            if check_coot_amber_type != ref_type:
                conf = m.GetConformer()
                coords = conf.GetAtomPosition(idx)
                cc = "(set-rotation-centre " + \
                     str(coords.x) + ' ' + str(coords.y) + ' ' + str(coords.z) + ')'
                s = 'failed to match mol_idx {} atom_idx {:>2} name {} coot_type {:>4} ref_type {:>2} at {}'
                print(s.format(mol_idx, idx, name, coot_amber_type, ref_type, cc))
            else:
                s = 'matched         mol_idx {} atom_idx {:>2} name {} coot_type {:>4} ref_type {}'
                print(s.format(mol_idx, idx, name, coot_amber_type, ref_type))

            p=(atom.GetProp('name'), coot_amber_type)
            aa_types.append(p)
        except TypeError:
            pass
	except KeyError:
	    # pf_atom_type was not set on an atom
	    pass
    return aa_types



if __name__ == '__main__':

   top_n = 100

   # write_mols(suppl) # only need do this once

   suppl = Chem.SDMolSupplier('parm/zinc-frag-100.sdf', removeHs=False)
   n_tot = 0
   ref_idx_offset={}
   for i,m in enumerate(suppl):
       n = m.GetNumAtoms()
       ref_idx_offset[i]=(n_tot, n_tot+n)
       n_tot +=n
       print(i,n, n_tot)
   print(ref_idx_offset)

   f = open('parm/zinc_p_f_types.txt')
   ref_types=[line.strip() for line in f.readlines()]
   f.close()

   for sdf_index in range(top_n):
      ref_types_for_mol=ref_types[ref_idx_offset[sdf_index][0]:ref_idx_offset[sdf_index][1]]    
      types = get_amber_types(suppl[sdf_index], sdf_index, ref_types_for_mol, check_2_chars_only=True)

