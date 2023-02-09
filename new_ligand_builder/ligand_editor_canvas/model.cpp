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
#include <GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/Geometry/point.h>
#include <boost/range/iterator_range.hpp>

using namespace coot::ligand_editor_canvas;


void CanvasMolecule::draw(GtkSnapshot* snapshot, const graphene_rect_t *bounds) const noexcept {
    auto x_offset = bounds->size.width / 2.0;
    auto y_offset = bounds->size.height / 2.0;
    auto scale_factor = 30.f;

    cairo_t *cr = gtk_snapshot_append_cairo(snapshot, bounds);
    // todo: change
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 3.0);
    for(const auto& bond: bonds) {
        cairo_move_to(cr, bond.first_atom_x * scale_factor + x_offset, bond.first_atom_y * scale_factor + y_offset);
        cairo_line_to(cr, bond.second_atom_x * scale_factor + x_offset, bond.second_atom_y * scale_factor + y_offset);
        cairo_stroke(cr);
        g_debug("TODO: Implement drawing various bond kinds, colors, hightlights etc.");
    }
    for(const auto& atom: atoms) {
        g_debug("TODO: Implement drawing atoms");
    }
    cairo_destroy(cr);
}

CanvasMolecule::CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) {
    this->rdkit_molecule = std::move(rdkit_mol);
    this->lower_from_rdkit();
}

CanvasMolecule::BondType CanvasMolecule::bond_type_from_rdkit(RDKit::Bond::BondType rdkit_bond) {
    switch (rdkit_bond) {
        case RDKit::Bond::SINGLE: {
            return BondType::Single;
        }
        case RDKit::Bond::DOUBLE: {
            return BondType::Double;
        }
        case RDKit::Bond::TRIPLE: {
            return BondType::Triple;
        }
        case RDKit::Bond::AROMATIC:{
            g_warning("Todo: Take care of aromatic bonds.");
            return BondType::Single;
        }
        case RDKit::Bond::UNSPECIFIED:
        case RDKit::Bond::QUADRUPLE:
        case RDKit::Bond::QUINTUPLE:
        case RDKit::Bond::HEXTUPLE:
        case RDKit::Bond::ONEANDAHALF:
        case RDKit::Bond::TWOANDAHALF:
        case RDKit::Bond::THREEANDAHALF:
        case RDKit::Bond::FOURANDAHALF:
        case RDKit::Bond::FIVEANDAHALF:
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

    // The following code is based on RDKit documentation.
    // In fact, it's almost a copy-paste.
    // I'm not exactly sure if what I did here is the most optimal solution
    // but it seems to work.
    RDDepict::compute2DCoords(*this->rdkit_molecule,&coordinate_map,true,false);


    RDKit::MatchVectType matchVect;
    if(! RDKit::SubstructMatch( *this->rdkit_molecule , *this->rdkit_molecule , matchVect ) ) {
        throw std::runtime_error("SubstractMatch failed.");
    }
    RDKit::Conformer &conf = this->rdkit_molecule->getConformer();
    for(auto mv: matchVect) {
        RDGeom::Point3D pt3 = conf.getAtomPos( mv.first );
        RDGeom::Point2D pt2( pt3.x , pt3.y );
        coordinate_map[mv.second] = pt2;
    }
    RDDepict::compute2DCoords( *this->rdkit_molecule , &coordinate_map );


    /// Used to avoid duplicating bonds
    std::set<unsigned int> processed_atoms_indices;

    if(coordinate_map.empty()) {
        throw std::runtime_error("RDKit coordinate mapping is empty");
    }
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

