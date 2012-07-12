
import gzip
from rdkit import Chem
from rdkit.Chem import Recap
global fragment_dict

fragment_dict = {}

class fragment_mol_info_t:
    mol = False

    def __init__(self, mol_in):
        self.mol = mol_in
        

def read_zinc_sdf(gzip_file_name):
    suppl = Chem.SDMolSupplier(gzip_file_name)
    # suppl = Chem.ForwardSDMolSupplier(gzip.open(gzip_file_name))
    for mol in suppl:
        try:
            # print mol.GetNumAtoms()
            hierarch = Recap.RecapDecompose(mol)
            keys = hierarch.GetLeaves().keys()
            for key in keys:
                try:
                    leaf_mol = hierarch.GetLeaves()[key].mol
                    try:
                        fragment_dict[key].append(mol)
                    except KeyError:
                        new_list = []
                        mol_info = fragment_mol_info_t(mol)
                        new_list.append(mol_info)
                        fragment_dict[key] = new_list
                    
                except KeyError:
                    print "type error in leave from decomposition"

                
        except TypeError as e:
            print "outer attribute error", e.what()
            pass

    print " ================================== result ================================= "
    print "n items in fragment_dict:", len(fragment_dict), " from", len(suppl), "molecules."
    l = fragment_dict.items()
    for i in fragment_dict:
        print "dict item ", i, "from",Chem.MolToSmiles(fragment_dict[i][0].mol)

f = '/extra/paule/zinc/23_p0.0.sdf.gz'
# f = 'LIG.sdf'
f = 'x2.sdf'

read_zinc_sdf(f)

