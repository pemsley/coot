#include "model.hpp"
#include "cairo-deprecated.h"
#include "cairo.h"
#include <boost/range/iterator_range_core.hpp>
#include <stdexcept>
#include <set>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2D.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DCairo.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DSVG.h>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/Geometry/point.h>
#include <boost/range/iterator_range.hpp>

using namespace coot::ligand_editor_canvas;


void CanvasMolecule::draw(GtkSnapshot* snapshot, const graphene_rect_t *bounds) const noexcept {
    cairo_t *cr = gtk_snapshot_append_cairo(snapshot, bounds);
    // todo: change
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    for(const auto& bond: bonds) {
        cairo_move_to(cr, bond.first_atom_x, bond.first_atom_y);
        cairo_line_to(cr, bond.second_atom_x, bond.second_atom_y);
        g_debug("TODO: Implement drawing various bond kinds, colors, hightlights etc.");
    }
    for(const auto& atom: atoms) {
        g_debug("TODO: Implement drawing atoms");
    }
    g_object_unref(cr);
}

CanvasMolecule::CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) {
    this->rdkit_molecule = std::move(rdkit_mol);
    this->lower_from_rdkit();
}

CanvasMolecule::BondType CanvasMolecule::bond_type_from_rdkit(RDKit::Bond::BondType rdkit_bond) {
    switch (rdkit_bond) {
        case RDKit::Bond::DOUBLE: {
            return BondType::Single;
        }
        case RDKit::Bond::TRIPLE: {
            return BondType::Triple;
        }
        case RDKit::Bond::UNSPECIFIED:
        case RDKit::Bond::SINGLE:
        case RDKit::Bond::QUADRUPLE:
        case RDKit::Bond::QUINTUPLE:
        case RDKit::Bond::HEXTUPLE:
        case RDKit::Bond::ONEANDAHALF:
        case RDKit::Bond::TWOANDAHALF:
        case RDKit::Bond::THREEANDAHALF:
        case RDKit::Bond::FOURANDAHALF:
        case RDKit::Bond::FIVEANDAHALF:
        case RDKit::Bond::AROMATIC:
        case RDKit::Bond::IONIC:
        case RDKit::Bond::HYDROGEN:
        case RDKit::Bond::THREECENTER:
        case RDKit::Bond::DATIVEONE:
        case RDKit::Bond::DATIVE:
        case RDKit::Bond::DATIVEL:
        case RDKit::Bond::DATIVER:
        case RDKit::Bond::OTHER:
        case RDKit::Bond::ZERO:
        default: {
            g_warning("Unhandled RDKit bond type: %i", rdkit_bond);
            return BondType::Single;
        }
    }
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

    /// Used to avoid duplicating bonds
    std::set<unsigned int> processed_atoms_indices;

    // 2.2 Process atoms and compute bonds
    for(const auto& [atom_idx,plane_point]: coordinate_map) {
        const auto* rdkit_atom = this->rdkit_molecule->getAtomWithIdx(atom_idx);
        auto canvas_atom = CanvasMolecule::Atom();
        // todo: change this
        canvas_atom.color = CanvasMolecule::AtomColor::Black;
        canvas_atom.highlighted = false;
        canvas_atom.idx = atom_idx;
        canvas_atom.symbol = rdkit_atom->getSymbol();
        canvas_atom.x = plane_point.x;
        canvas_atom.y = plane_point.y;

        for(const auto& bond: boost::make_iterator_range(this->rdkit_molecule->getAtomBonds(rdkit_atom))) {
            // Based on `getAtomBonds` documentation.
            // Seems weird but we have to do it that way.
            const auto* bond_ptr = (*this->rdkit_molecule)[bond];
            // We don't want to have duplicate bonds of atoms that we have already processed
            // so we skip them.
            auto first_atom_idx = bond_ptr->getBeginAtomIdx();
            auto second_atom_idx = bond_ptr->getEndAtomIdx();;
            if(processed_atoms_indices.find(first_atom_idx) != processed_atoms_indices.end() 
            || processed_atoms_indices.find(second_atom_idx) != processed_atoms_indices.end()) {
                continue;
            }

            auto canvas_bond = CanvasMolecule::Bond();

            canvas_bond.first_atom_idx = first_atom_idx;
            canvas_bond.first_atom_x = coordinate_map[first_atom_idx].x;
            canvas_bond.first_atom_y = coordinate_map[first_atom_idx].y;
            
            canvas_bond.second_atom_idx = second_atom_idx;
            canvas_bond.second_atom_x = coordinate_map[second_atom_idx].x;
            canvas_bond.second_atom_y = coordinate_map[second_atom_idx].y;

            canvas_bond.highlighted = false;
            // todo: implement processing aromatic bonds
            canvas_bond.type = bond_type_from_rdkit(bond_ptr->getBondType());

            this->bonds.push_back(std::move(canvas_bond));
        }

        this->atoms.push_back(std::move(canvas_atom));
        // Mark the atom as processed
        processed_atoms_indices.insert(atom_idx);
    }
    
}

