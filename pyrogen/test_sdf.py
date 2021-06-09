
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
import enum
from enum import Enum

# for dictionaries
import pyrogen_swig as pysw

def write_mols(suppl):
    for idx,mol in enumerate(suppl):
        fn = 'pf-'+str(idx)+'.sdf'
        # print('writing {}'.format(fn))
        print(Chem.MolToMolBlock(mol), file=file(fn, 'w'))

# return True or Less or Greater or Confused_T_1 or Confused_T_2
match_type = Enum('match', 'MATCH ABOVE BELOW CONFUSED_T1A CONFUSED_T1B CONFUSED_T2A CONFUSED_T2B UNKNOWN WRONG_ELEMENT')

# for example CT, CJ returns BELOW because CT (check_coot_amber_type) is below CJ (ref_type)
#
class compare_types:

    def check(self, check_coot_amber_type, ref_type, Cs):
        if check_coot_amber_type == ref_type:
            return match_type.MATCH
        if not check_coot_amber_type in Cs:
            return match_type.UNKNOWN
        rank_c_forward  = self.get_rank(check_coot_amber_type, Cs)
        rank_r_forward  = self.get_rank(ref_type, Cs)
        rank_c_backward = self.get_rank(check_coot_amber_type, Cs, reverse_search=True)
        rank_r_backward = self.get_rank(ref_type, Cs, reverse_search=True)
        name = "####"
        try:
            name = self.atom.GetProp('name')
        except 'green' as e:
            print('no name for {}'.format(self.atom))
        s = "{}: {} should be {} in {} "
        s += "rank_c_forward {} rank_r_forward {} rank_c_backward {} rank_r_backward {}"
        if self.debug_output:
           print(s.format(name, check_coot_amber_type, ref_type, Cs, rank_c_forward, rank_r_forward, rank_c_backward, rank_r_backward))
        if rank_c_forward < rank_r_forward:
            if rank_c_backward < rank_r_backward:
                return match_type.ABOVE
        if rank_c_forward > rank_r_forward:
            if rank_c_backward > rank_r_backward:
                return match_type.BELOW
        # there is some type of confusion
        if rank_c_forward < rank_r_forward:
            if rank_c_backward < rank_r_backward:
                return match_type.CONFUSED_T1A
            else:
                return match_type.CONFUSED_T1B
        if rank_c_forward > rank_r_forward:
            if rank_c_backward < rank_r_backward:
                return match_type.CONFUSED_T2A
            else:
                return match_type.CONFUSED_T1B
        return match_type.UNKNOWN

    def get_rank(self, a_type, types, reverse_search=False):
        if not reverse_search:
            for i,t in enumerate(types):
                if t == a_type:
                    return i
            print('compare_types() error')  # shouldn't happen
            return 'rank-not-found'
        else:
            n = len(types)
            for i in range(n-1,-1,-1):
               if a_type == types[i]:
                  return i
        print("error in get_rank()")
        return -1

    def test_compare_2(self):

       Ns = ['NC', 'N2', 'NB', 'N', 'N2', 'N']
       # Ns = ['NC', 'N2', 'NC']
       for N_1 in Ns:
          r1 = self.get_rank(N_1, Ns)
          print("type {:>2} rank {:>2}".format(N_1, r1))

       N_1 = 'NB'
       N_2 = 'N2'
       r = self.check(N_1, N_1, Ns)
       print("{:>2} and {:>2} {:>2}".format(N_1, N_1, r))
       r = self.check(N_2, N_2, Ns)
       print("{:>2} and {:>2} {:>2}".format(N_2, N_2, r))
       r = self.check(N_1, N_2, Ns)
       print("{:>2} and {:>2} {:>2}".format(N_1, N_2, r))
       r = self.check(N_2, N_1, Ns)
       print("{:>2} and {:>2} {:>2}".format(N_2, N_1, r))

    def compare(self, check_coot_amber_type, ref_type):

       if check_coot_amber_type == ref_type:
           return match_type.MATCH
       ele = ref_type[0];
       if ele == 'C':
           r = self.check(check_coot_amber_type, ref_type, self.Cs)
           return r
       else:
           if ele == 'N':
               r = self.check(check_coot_amber_type, ref_type, self.Ns)
               return r
           else:
               if ele == 'O':
                   r = self.check(check_coot_amber_type, ref_type, self.Os)
                   return r
               else:
                   if ele == 'S':
                       r = self.check(check_coot_amber_type, ref_type, self.Ss)
                       return r
                   else:
                       return match_type.WRONG_ELEMENT # unknown element


    def __init__(self, atom_in, debug_output=False):

       self.debug_output = False
       if debug_output:
           self.debug_output = True
       self.atom = atom_in
       self.Cs = ['CB', 'CR', 'CB', 'CW', 'CC', 'CR', 'CW', 'CC', 'CR', 'C', 'CP', 'C',
                  'CR', 'CB', 'C*', 'CA', 'CM', 'C2', 'CJ', 'CT']
       self.Ns = ['NB', 'NJ', 'NL', 'N3',  'N3', 'NC', 'N2', 'NB', 'N', 'N2', 'N',
                  'NA', 'N3', 'ND', "N2", "N", "N*", "N3", "Nu"]
       self.Os = ['OW', 'OH', 'OS', 'O2', 'OS', 'O2', 'O', 'Ou']
       self.Ss = ['SH', 'S', 'SO', 'SD', 'Su']

       # self.test_compare_2()

def get_parmfrosst_types(m, mol_idx, ref_types=None, check_2_chars_only=True, make_dictionaries=False):
    Chem.AllChem.ComputeGasteigerCharges(m)
    names = pyrogen.add_atom_names(m)
    atom_types.set_atom_types(m)
    atom_types.set_parmfrosst_atom_types(m)
    aa_types=[]

    if make_dictionaries:
       cif_file_name = "test-pyrogen-"+str(mol_idx)+'.cif'
       comp_id = "PF_" + str(mol_idx)
       pysw.mmcif_dict_from_mol(comp_id, m.GetProp('_Name'), m, cif_file_name, False, False, True)

    for idx,atom in enumerate(m.GetAtoms()):
        try:
            # in pathological cases, pf_atom_type is not set
            coot_amber_type = atom.GetProp('pf_atom_type')
            # if we were passed ref_types, compare the atom type
            name = atom.GetProp('name')
            ref_type=ref_types[idx]
            if ref_type == "Nstar":
               ref_type = 'N*'
            if ref_type == "Cstar":
               ref_type = 'C*'
            check_coot_amber_type = coot_amber_type
            if check_2_chars_only:
                check_coot_amber_type = coot_amber_type[:2]
                # C -> Ca or Cb or Cc etc is a special case
                if coot_amber_type in ["Ca", "Cax", "Cb", "Cc", "Cd", "Ce", "Cf", "Cg", "Ch", "Ci", "Cj", "C-"]:
                    check_coot_amber_type = "C"
                if coot_amber_type in ["Sa", "Sax", "Sb", "Sc", "Sd", "Se", "Sf", "Sg", "Sh", "Si", "Sj"]:
                    check_coot_amber_type = "S"
            if check_coot_amber_type != ref_type:
                ct = compare_types(atom)
                match_result = ct.compare(check_coot_amber_type, ref_type)
                conf = m.GetConformer()
                coords = conf.GetAtomPosition(idx)
                cc = "(set-rotation-centre " + \
                     str(coords.x) + ' ' + str(coords.y) + ' ' + str(coords.z) + ')'
                s = 'failed to match mol_idx {:>2} atom_idx {:>2} name {} coot_type {:>10} '
                s += "ref_type {:>2} rank {:>18} at {}"
                # Note: A match.BELOW means that we too specific in the SMILES for the
                #       coot_amber_type. To correct, make the SMILES or general
                #       or add another one that picks up this case.
                #       A match.ABOVE means that we were too general in the SMILES.
                #       This is more problematic.  We need to get rid of these by
                #       making the SMILES more specific
                print(s.format(mol_idx, idx, name, coot_amber_type, ref_type, match_result, cc))
            else:
                s = 'matched         mol_idx {} atom_idx {:>2} name {} coot_type {:>4} ref_type {}'
                print(s.format(mol_idx, idx, name, coot_amber_type, ref_type))

            p=(atom.GetProp('name'), coot_amber_type)
            aa_types.append(p)
        except TypeError as e:
            print('TypeError {}'.format(e))
            pass
        except KeyError as e:
            # pf_atom_type was not set on an atom
            print('KeyError {}'.format(e))
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
       # print(i,n, n_tot)
   # print(ref_idx_offset)

   f = open('parm/zinc_p_f_types.txt')
   ref_types=[line.strip() for line in f.readlines()]
   f.close()

   for sdf_index in range(top_n):
      ref_types_for_mol=ref_types[ref_idx_offset[sdf_index][0]:ref_idx_offset[sdf_index][1]]
      types = get_amber_types(suppl[sdf_index], sdf_index, ref_types_for_mol, check_2_chars_only=True)
#
