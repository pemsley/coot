#include "ligand_builder_utils.hpp"

void coot::ligand_editor::remove_non_polar_hydrogens(RDKit::RWMol* mol) {
    std::vector<RDKit::Atom*> atoms_to_be_removed;

    auto atoms = mol->atoms();
    for(RDKit::Atom* atom: atoms) {
        if(atom->getAtomicNum() == 1) {
            if(atom->getFormalCharge() == 0) {
                atoms_to_be_removed.push_back(atom);
            }
        }
    }

    for(RDKit::Atom* atom: atoms_to_be_removed) {
        mol->removeAtom(atom);
    }
}