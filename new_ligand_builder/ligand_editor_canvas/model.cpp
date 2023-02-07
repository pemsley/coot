#include "model.hpp"
#include <stdexcept>

using namespace coot::ligand_editor_canvas;


void CanvasMolecule::draw(GtkSnapshot* snapshot) const noexcept {
    for(const auto& bond: bonds) {

    }
    for(const auto& atom: atoms) {

    }
}

CanvasMolecule::CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) {
    this->rdkit_molecule = std::move(rdkit_mol);
    this->lower_from_rdkit();
}

void CanvasMolecule::lower_from_rdkit() {
    // 1. Clear what we have
    this->atoms.clear();
    this->bonds.clear();
    // 2. Do the lowering
    throw std::runtime_error("Lowering unimplemented!");
}

