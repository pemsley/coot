
import gzip
from rdkit import Chem
from rdkit.Chem import Recap
global fragment_dict

fragment_dict = {}

class fragment_mol_info_t():

    # fragment_in is the key
    def __init__(self, fragment_in):
        self.mol = Chem.MolFromSmiles(fragment_in)
        self.n_atoms = self.mol.GetNumAtoms()
        self.orig_mol_list = []
        

    def add_mol(self, mol_in):
        self.orig_mol_list.append(mol_in)

def read_zinc_sdf(gzip_file_name):
    suppl = Chem.SDMolSupplier(gzip_file_name)
    # suppl = Chem.ForwardSDMolSupplier(gzip.open(gzip_file_name))
    for mol in suppl:
        try:
            # print mol.GetNumAtoms()
            hierarch = Recap.RecapDecompose(mol)
            keys = hierarch.GetLeaves().keys()
            mol_info = False
            for key in keys:
                try:
                    fragment_dict[key].add_mol(mol)
                except KeyError:
                    mol_info = fragment_mol_info_t(key)
                    mol_info.add_mol(mol)
                    fragment_dict[key] = mol_info
                    
        except TypeError as e:
            print "outer attribute error", e.what()
            pass

    print "n items in fragment_dict:", len(fragment_dict), " from", len(suppl), "molecules."
    for i in fragment_dict:
        if (fragment_dict[i].n_atoms < 5):
            frag_info = fragment_dict[i]
            mol_smiles = Chem.MolToSmiles(frag_info.orig_mol_list[0])
            print "   ", frag_info.n_atoms, i, "    from   ", mol_smiles

        

f = '/extra/paule/zinc/23_p0.0.sdf.gz'
# f = 'LIG.sdf'
f = 'x2.sdf'

read_zinc_sdf(f)

