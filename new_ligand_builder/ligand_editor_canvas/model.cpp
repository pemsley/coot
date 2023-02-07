#include "model.hpp"
#include <stdexcept>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2D.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DCairo.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DSVG.h>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/Geometry/point.h>

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

    // 2.1 Get geometry info

    // Maps atom indices to 2D points
    RDGeom::INT_POINT2D_MAP coordinate_map;
    RDDepict::compute2DCoords(*this->rdkit_molecule,&coordinate_map);

    // 2.2 Process atoms and compute bonds
    for(const auto& [atom_idx,plane_point]: coordinate_map) {
        
    }
    
}

