
from subprocess import call
from rdkit import Chem
from rdkit.Chem import AllChem

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

def make_restraints(smiles_string, sdf_file_name, mmcif_dict_name):

   m = Chem.MolFromSmiles(smiles_string)
   AllChem.EmbedMolecule(m)
   AllChem.UFFOptimizeMolecule(m)
   mb = Chem.MolToMolBlock(m)
   print >> file(sdf_file_name,'w'), mb
   mogul_ins_file_name = 'mogul.ins'
   mogul_out_file_name = 'mogul.out'
   run_mogul(sdf_file_name, mogul_ins_file_name, mogul_out_file_name)


smiles_string = 'Cc1ccccc1'
smiles_string = 'O(c3cc4c2c1c[n+](ccc1ccc2nc4cc3)CCN5CCC(CC5)CCCC%10CCN(CC[n+]9cc8c7c6cc(OC)ccc6nc7ccc8cc9)CC%10)C'
smiles_string = 'COc1ccc2[nH]c3ccc4cc[n+](CCN5CCC(CCCC6CCN(CC6)CC[n+]7ccc8ccc9[nH]c%10ccc(OC)cc%10c9c8c7)CC5)cc4c3c2c1'
sdf_file_name = 'test.sdf'
make_restraints(smiles_string, sdf_file_name, 'rest-test.cif')

