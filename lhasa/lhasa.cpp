#include "lhasa.hpp"
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

std::unique_ptr<RDKit::RWMol> lhasa::rdkit_mol_from_smiles(std::string smiles) {
    std::unique_ptr<RDKit::RWMol> ret(RDKit::SmilesToMol(smiles));
    return ret;
}

std::string lhasa::rdkit_mol_to_smiles(RDKit::ROMol& mol) {
    auto ret = RDKit::MolToSmiles(mol);
    return ret;
}
