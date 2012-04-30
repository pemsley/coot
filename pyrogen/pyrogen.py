
from subprocess import call
from rdkit import Chem
from rdkit.Chem import AllChem
import coot_libs as coot

import jay-util

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
     f.write('torsion all\n')
     f.write('ring all\n')
     f.close()
   return f


def run_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name):
   f = make_mogul_ins_file(mogul_ins_file_name, mogul_out_file_name, sdf_file_name)
   if f: 
      print 'now run mogul using ins file %s' % mogul_ins_file_name
      call(['mogul', '-ins', mogul_ins_file_name])
      return True
   else:
      return False

def atom_name_from_atomic_number_and_count(element, count):
    name = element
    name += str(count)
    return name

def add_atom_names(mol):
    nz = {}
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z in nz:
            nz[z]  = nz[z] + 1
        else:
            nz[z] = 1;
        ele = atom.GetSymbol()
        name = atom_name_from_atomic_number_and_count(ele, nz[z])
        atom.SetProp("name", name)
        
#     for i in nz:
#         print "atomic_number", i, ": count", nz[i]

    

def make_restraints(smiles_string, sdf_file_name, mmcif_dict_name):

   mogol_exe = which('mogul')

   if (mogol_exe == None):
      print "mogul not found in path"
      return False
   else:
      m = Chem.MolFromSmiles(smiles_string)
      AllChem.EmbedMolecule(m)
      AllChem.UFFOptimizeMolecule(m)
      add_atom_names(m)
#       for atom in m.GetAtoms():
#           print atom, atom.GetProp('name')

      mb = Chem.MolToMolBlock(m)
      print >> file(sdf_file_name,'w'), mb
      mogul_ins_file_name = 'mogul.ins'
      mogul_out_file_name = 'mogul.out'
      #nmogul_state = run_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name)
      mogul_state = False
      if mogul_state:
          mogul = coot.mogul(mogul_out_file_name)
          bor = make_restraints_for_bond_orders(m)
          write_restraints(bor, mmcif_dict_name)
          return True
      return mogul_state


smiles_string = 'Cc1ccccc1'
smiles_string = 'Oc1ccccc1'
# smiles_string = 'O(c3cc4c2c1c[n+](ccc1ccc2nc4cc3)CCN5CCC(CC5)CCCC%10CCN(CC[n+]9cc8c7c6cc(OC)ccc6nc7ccc8cc9)CC%10)C'
# smiles_string = 'COc1ccc2[nH]c3ccc4cc[n+](CCN5CCC(CCCC6CCN(CC6)CC[n+]7ccc8ccc9[nH]c%10ccc(OC)cc%10c9c8c7)CC5)cc4c3c2c1'
sdf_file_name = 'test.sdf'
make_restraints(smiles_string, sdf_file_name, 'rest-test.cif')

