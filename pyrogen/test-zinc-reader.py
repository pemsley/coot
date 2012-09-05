
import os
import gzip
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem import ChemicalFeatures

global fragment_dict

fragment_dict = {}

class fragment_mol_info_t():

    def add_mol(self, mol_in):
        self.orig_mol_list.append(mol_in)

    # does mol contain any aromatic atoms?
    def get_is_aromatic(self, fragment_mol):
        for atom in fragment_mol.GetAtoms():
            if (atom.GetIsAromatic()):
                return True
        return False

    # fragment_in is the key
    def __init__(self, fragment_in, feature_factory):
        self.mol = Chem.MolFromSmiles(fragment_in)
        self.n_atoms = self.mol.GetNumAtoms()
        self.is_aromatic = self.get_is_aromatic(self.mol)
        self.feats = factory.GetFeaturesForMol(self.mol)
        self.orig_mol_list = []


class fragments_dictionary():
    must_contain_donor=False
    must_contain_acceptor=False
    must_contain_hydrophobe=False
    must_contain_aromatic=False
    max_n_atoms=100

    def passes_queries(self, frag_info):
        if (self.must_contain_donor):
            if (not (self.find_feat(frag_info.feats, "Donor"))):
                return False
        if (self.must_contain_acceptor):
            if (not (self.find_feat(frag_info.feats, "Acceptor"))):
                return False
        if (self.must_contain_hydrophobe):
            if (not (self.find_feat(frag_info.feats, "Hydrophobe"))):
                return False
        if (self.must_contain_aromatic):
            if (not (self.find_feat(frag_info.feats, "Aromatic"))):
                return False
        return True

    def find_feat(self, feats, feat_type):
        found = False
        for f in feats:
            if (f.GetFamily() == feat_type):
                found = True
                break
        return found
        
                
    def __init__(self):
        pass

    def query_dictionary(self, max_n_atoms, contains_donor=False, contains_acceptor=False, contains_hydrophobe=False, contains_aromatic=False):
        print "n items in fragment_dict:", len(fragment_dict)
        self.must_contain_donor      = contains_donor
        self.must_contain_acceptor   = contains_acceptor
        self.must_contain_hydrophobe = contains_hydrophobe
        self.must_contain_aromatic   = contains_aromatic
        self.max_n_atoms = max_n_atoms
        results = []
        
        for i in fragment_dict:
            if (fragment_dict[i].n_atoms <= max_n_atoms):
                frag_info = fragment_dict[i]
                
                if (self.passes_queries(frag_info)):
                    results.append(frag_info)
                else:
                    if (0):
                        print "did not pass query tests"
                        for f in frag_info.feats:
                            print "        ", f.GetFamily(), f.GetType()
        return results


def read_zinc_sdf(suppl, feat_factory):
    # print "Supplier provides", len(suppl), "molecules."
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
                    mol_info = fragment_mol_info_t(key, feat_factory)
                    mol_info.add_mol(mol)
                    fragment_dict[key] = mol_info
                    
        except TypeError as e:
            print "Type Error", e
            pass



if __name__ == "__main__":

    fdef_file_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdef_file_name)

    if factory:

        gzf = '/extra/paule/zinc/23_p0.0.sdf.gz'
        f = 'x3.sdf'
        f = 'x2.sdf'
        # f = 'x1.sdf'
        
        # read_zinc_sdf(Chem.ForwardSDMolSupplier(gzip.open(gzf)), factory)
        read_zinc_sdf(Chem.SDMolSupplier(f), factory)
        fd = fragments_dictionary()
        results = fd.query_dictionary(7, False, True, False, True)
        for frag_info in results:
            mol_smiles = Chem.MolToSmiles(frag_info.orig_mol_list[0])
            zinc_name = ''
            try:
                zinc_name = frag_info.orig_mol_list[0].GetProp('_Name')
            except KeyError:
                pass
            print " ", frag_info.n_atoms, "    from   ", mol_smiles, zinc_name
            for f in frag_info.feats:
                print "        ", f.GetFamily(), f.GetType()
            

