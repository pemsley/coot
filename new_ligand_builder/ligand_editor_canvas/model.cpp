#include "model.hpp"
#include "cairo-deprecated.h"
#include "cairo.h"
#include <boost/range/iterator_range_core.hpp>
#include <optional>
#include <stdexcept>
#include <algorithm>
#include <set>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2D.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DCairo.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DSVG.h>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/Geometry/point.h>
#include <cmath>
#include <boost/range/iterator_range.hpp>

using namespace coot::ligand_editor_canvas;

const float CanvasMolecule::ATOM_HITBOX_RADIUS = 10.f;
const float CanvasMolecule::BASE_SCALE_FACTOR = 30.f;

float CanvasMolecule::get_scale() const noexcept {
    // todo: Add internal canvas scaling
    return BASE_SCALE_FACTOR;
}

CanvasMolecule::MaybeAtomOrBond CanvasMolecule::resolve_click(int x, int y) const noexcept {
    float scale = this->get_scale();
    auto x_offset = this->_x_offset;
    auto y_offset = this->_y_offset;
    // atoms first 
    for(const auto& atom: this->atoms) {
        float atom_x = atom.x * scale + x_offset;
        float atom_y = atom.y * scale + y_offset;
        // Circle equation. Checks if click coords are within atom's radius
        if (ATOM_HITBOX_RADIUS * ATOM_HITBOX_RADIUS >= std::pow(atom_x - x,2.f) + std::pow(atom_y - y,2.f)) {
            return atom;
        }
    }
    // then bonds
    for(const auto& bond: this->bonds) {
        // todo: figure out how to do hitboxes of bonds. Clipper might be very useful.
        g_warning_once("todo: implement hitboxes of bonds.");
    }
    return std::nullopt;
}

void CanvasMolecule::set_offset_from_bounds(const graphene_rect_t *bounds) noexcept {
    this->_x_offset = bounds->size.width / 2.0;
    this->_y_offset = bounds->size.height / 2.0;
}

void CanvasMolecule::draw(GtkSnapshot* snapshot, PangoLayout* pango_layout, const graphene_rect_t *bounds) const noexcept {
    auto scale_factor = this->get_scale();
    auto x_offset = this->_x_offset;
    auto y_offset = this->_y_offset;

    cairo_t *cr = gtk_snapshot_append_cairo(snapshot, bounds);
    

    for(const auto& bond: bonds) {
        if(bond.highlighted) {
            cairo_set_line_width(cr, 5.0);
            cairo_set_source_rgb(cr, 0.0, 1.0, 0.0);
        } else {
            cairo_set_line_width(cr, 3.0);
            cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        }
        cairo_move_to(cr, bond.first_atom_x * scale_factor + x_offset, bond.first_atom_y * scale_factor + y_offset);
        cairo_line_to(cr, bond.second_atom_x * scale_factor + x_offset, bond.second_atom_y * scale_factor + y_offset);
        cairo_stroke(cr);
        g_warning_once("TODO: Implement drawing various bond kinds, colors etc.");
    }

    cairo_set_line_width(cr, 0.5);
    for(const auto& atom: atoms) {
        // Used to make the texts centered where they should be.
        int layout_width, layout_height;

        cairo_move_to(cr, atom.x * scale_factor + x_offset + ATOM_HITBOX_RADIUS, atom.y * scale_factor + y_offset);

        if(atom.highlighted) {
            cairo_set_source_rgb(cr, 0.0, 1.0, 0.0);
            cairo_arc(cr, atom.x * scale_factor + x_offset, atom.y * scale_factor + y_offset,ATOM_HITBOX_RADIUS,0,M_PI * 2.0);
            cairo_stroke_preserve(cr);
            cairo_set_source_rgba(cr, 0.0, 1.0, 0.0, 0.5);
            cairo_fill(cr);
        }

        auto render_text = [&](const std::string& t){
            //todo: color and size
            std::string markup = "<span color=\"navy\" weight=\"bold\" size=\"x-large\">" + t + "</span>";
            pango_layout_set_markup(pango_layout,markup.c_str(),-1);
            pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
            cairo_move_to(cr, atom.x * scale_factor + x_offset - layout_width/2.f, atom.y * scale_factor + y_offset - layout_height/2.f);
            pango_cairo_show_layout(cr, pango_layout);
        };

        g_warning_once("TODO: Correctly implement drawing atoms");
        
        if (atom.symbol == "C") {
            // Ignore drawing Carbon
        } else if (atom.symbol == "H") {
            // Ignore hydrogens for now
        } else {
           render_text(atom.symbol);
        }
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
    std::sort(this->atoms.begin(),this->atoms.end(),[](const auto& lhs, const auto& rhs){
        return lhs.idx < rhs.idx;
    });
    
}

void CanvasMolecule::highlight_atom(int atom_idx) {
    this->atoms.at(atom_idx).highlighted = true;
}

void CanvasMolecule::highlight_bond(int atom_a, int atom_b) {
    auto target = std::find_if(this->bonds.begin(),this->bonds.end(),[=](const auto& bond){
        return bond.first_atom_idx == atom_a && bond.second_atom_idx == atom_b;
    });
    if (target == this->bonds.end()) {
        throw std::runtime_error("Bond doesn't exist");
    }
    target->highlighted = true;
}

void CanvasMolecule::clear_highlights() {
    for(auto& bond: this->bonds) {
        bond.highlighted = false;
    }
    for(auto& atom: this->atoms) {
        atom.highlighted = false;
    }
}