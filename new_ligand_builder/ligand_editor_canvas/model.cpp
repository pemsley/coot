#include "model.hpp"
#include "cairo.h"
#include <optional>
#include <stdexcept>
#include <algorithm>
#include <set>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2D.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DCairo.h>
// // #include <rdkit/GraphMol/MolDraw2D/MolDraw2DSVG.h>
// #include <rdkit/GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MolOps.h>
#include <cmath>
#include <boost/range/iterator_range.hpp>
#include <string>
#include <tuple>
#include <utility>

using namespace coot::ligand_editor_canvas;

const float CanvasMolecule::ATOM_HITBOX_RADIUS = 12.f;
const float CanvasMolecule::BOND_DISTANCE_BOUNDARY = 10.f;
const float CanvasMolecule::BASE_SCALE_FACTOR = 30.f;
const float CanvasMolecule::BOND_LINE_SEPARATION = 0.2f;

float CanvasMolecule::get_scale() const noexcept {
    return BASE_SCALE_FACTOR * this->canvas_scale;
}

void CanvasMolecule::set_canvas_scale(float scale) {
    this->canvas_scale = scale;
}

CanvasMolecule::MaybeAtomOrBond CanvasMolecule::resolve_click(int x, int y) const noexcept {
    float scale = this->get_scale();
    auto x_offset = this->x_canvas_size_adjustment + this->x_canvas_translation;
    auto y_offset = this->y_canvas_size_adjustment + this->y_canvas_translation;
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
        // 1. Find the point lying in the middle of the segment representing the bond
        float first_atom_x = bond.first_atom_x * scale + x_offset;
        float first_atom_y = bond.first_atom_y * scale + y_offset;
        float second_atom_x = bond.second_atom_x * scale + x_offset;
        float second_atom_y = bond.second_atom_y * scale + y_offset;
        float bond_center_x = (first_atom_x + second_atom_x) / 2.0f;
        float bond_center_y = (first_atom_y + second_atom_y) / 2.0f;
        // 2. Compute its' distance to the edges and let it be the radius of a circle
        //    inside of which the click location has to be found (as a pre-condition)
        float bounding_circle_radius_squared = std::pow(first_atom_x - bond_center_x,2.f) + std::pow(first_atom_y - bond_center_y,2.f);
        // If our distance from the bond's center point is greater than the radius
        // then this means that we're outside of the circle of interest and so the current bond can be skipped.
        // Operating on squared values saves us from the expensive square root operation.
        if (bounding_circle_radius_squared < std::pow(bond_center_x - x,2.f) + std::pow(bond_center_y - y,2.f)) {
            continue;
        }
        // 3. Use the formula for distance to a point from a line to determine if
        //    the click location qualifies as matching.

        float distance = std::fabs((second_atom_x - first_atom_x)*(first_atom_y - y) - (first_atom_x - x)*(second_atom_y - first_atom_y)) 
            / std::sqrt(std::pow(first_atom_x - second_atom_x,2.f) + std::pow(first_atom_y - second_atom_y,2.f));
        // We've got a match
        if (BOND_DISTANCE_BOUNDARY >= distance) {
            return bond;
        }
    }
    // nothing matches the click
    return std::nullopt;
}

void CanvasMolecule::set_canvas_size_adjustment_from_bounds(const graphene_rect_t *bounds) noexcept {
    this->x_canvas_size_adjustment = bounds->size.width / 2.0;
    this->y_canvas_size_adjustment = bounds->size.height / 2.0;
}

void CanvasMolecule::apply_canvas_translation(int delta_x, int delta_y) noexcept {
    this->x_canvas_translation += delta_x;
    this->y_canvas_translation += delta_y;
}

std::tuple<float,float,float> CanvasMolecule::atom_color_to_rgb(CanvasMolecule::AtomColor color) noexcept {
    switch (color) {
        case AtomColor::Green:{
            return std::make_tuple(0.0,1.0,0.0);
        }
        case AtomColor::Blue:{
            return std::make_tuple(0.0,0.0,1.0);
        }
        case AtomColor::Red:{
            return std::make_tuple(1.0,0.0,0.0);
        }
        case AtomColor::Black:
        default: {
            return std::make_tuple(0.0,0.0,0.0);
        }
    }
}

std::string CanvasMolecule::atom_color_to_html(CanvasMolecule::AtomColor color) noexcept {
    switch (color) {
        case AtomColor::Green:{
            return "#00FF00";
        }
        case AtomColor::Blue:{
            return "#0000FF";
        }
        case AtomColor::Red:{
            return "#FF0000";
        }
        case AtomColor::Black:
        default: {
            return "#000000";
        }
    }
}

CanvasMolecule::AtomColor CanvasMolecule::atom_color_from_rdkit(const RDKit::Atom * atom) noexcept {
    auto atomic_number = atom->getAtomicNum();
    switch(atomic_number) {
        // Nitrogen
        case 7: {
            return AtomColor::Blue;
        }
        // Oxygen
        case 8: {
            return AtomColor::Red;
        }
        // Chlorine
        case 17: {
            return AtomColor::Green;
        }
        // Carbon
        case 6:
        default: {
            return AtomColor::Black;
        }
    }
}

std::pair<float,float> CanvasMolecule::Bond::get_perpendicular_versor() const noexcept {
    float bond_vector_x = second_atom_x - first_atom_x;
    float bond_vector_y = second_atom_y - first_atom_y;
    float bond_vector_len = std::sqrt(std::pow(bond_vector_x,2.f) + std::pow(bond_vector_y,2.f));
    if (bond_vector_len == 0) {
        return std::make_pair(0.f, 0.f);
    } else {
        return std::make_pair(-bond_vector_y/bond_vector_len,bond_vector_x/bond_vector_len);
    }
}

void CanvasMolecule::draw(GtkSnapshot* snapshot, PangoLayout* pango_layout, const graphene_rect_t *bounds) const noexcept {
    auto scale_factor = this->get_scale();
    auto x_offset = this->x_canvas_size_adjustment + this->x_canvas_translation;
    auto y_offset = this->y_canvas_size_adjustment + this->y_canvas_translation;

    cairo_t *cr = gtk_snapshot_append_cairo(snapshot, bounds);
    

    for(const auto& bond: bonds) {
        if(bond.highlighted) {
            cairo_set_line_width(cr, 5.0);
            cairo_set_source_rgb(cr, 0.0, 1.0, 0.5);
        } else {
            cairo_set_line_width(cr, 3.0);
            cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        }
        cairo_move_to(cr, bond.first_atom_x * scale_factor + x_offset, bond.first_atom_y * scale_factor + y_offset);
        cairo_line_to(cr, bond.second_atom_x * scale_factor + x_offset, bond.second_atom_y * scale_factor + y_offset);
        cairo_stroke(cr);

        auto draw_side_bond_line = [&](bool addOrSub){
            auto [pv_x,pv_y] = bond.get_perpendicular_versor();
            if (!addOrSub) { // change sign of the versor
                pv_x *= -1.f;
                pv_y *= -1.f;
            }
            // Convert the versor to a vector of the desired length
            pv_x *= BOND_LINE_SEPARATION;
            pv_y *= BOND_LINE_SEPARATION;
            cairo_move_to(cr, (bond.first_atom_x + pv_x) * scale_factor + x_offset, (bond.first_atom_y + pv_y) * scale_factor + y_offset);
            cairo_line_to(cr, (bond.second_atom_x + pv_x) * scale_factor + x_offset, (bond.second_atom_y + pv_y) * scale_factor + y_offset);
            cairo_stroke(cr);
        };

        g_warning_once("TODO: Implement drawing different bond types correctly.");
        if(bond.type == BondType::Double) {
            // todo: store the information about what the boolean should be
            draw_side_bond_line(false);
        } else if(bond.type == BondType::Triple) {
            // "to the left"
            draw_side_bond_line(false);
            // "to the right"
            draw_side_bond_line(true);
        }
    }

    cairo_set_line_width(cr, 0.5);
    for(const auto& atom: atoms) {
        // Used to make the texts centered where they should be.
        int layout_width, layout_height;

        auto process_highlight = [&,cr,x_offset,y_offset,scale_factor](){
            if(atom.highlighted) {
                cairo_move_to(cr, atom.x * scale_factor + x_offset + ATOM_HITBOX_RADIUS, atom.y * scale_factor + y_offset);
                cairo_set_source_rgb(cr, 0.0, 1.0, 0.5);
                cairo_arc(cr, atom.x * scale_factor + x_offset, atom.y * scale_factor + y_offset,ATOM_HITBOX_RADIUS,0,M_PI * 2.0);
                cairo_stroke_preserve(cr);
                cairo_set_source_rgba(cr, 0.0, 1.0, 0.5, 0.5);
                cairo_fill(cr);
            }
        };
        
        // an alternative to doing this is to shorten the bonds at the time of lowering
        auto render_white_background = [&,cr,x_offset,y_offset,scale_factor]{
            cairo_move_to(cr, atom.x * scale_factor + x_offset + ATOM_HITBOX_RADIUS, atom.y * scale_factor + y_offset);
            cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
            cairo_arc(cr, atom.x * scale_factor + x_offset, atom.y * scale_factor + y_offset,ATOM_HITBOX_RADIUS,0,M_PI * 2.0);
            cairo_fill(cr);
        };

        auto render_text = [&](const std::string& t, AtomColor color,bool bold = false){
            auto [r,g,b] = atom_color_to_rgb(color);
            std::string color_str = atom_color_to_html(color);
            std::string weight_str = bold ? "bold" : "normal";
            std::string markup = "<span color=\"" + color_str + "\" weight=\"" + weight_str + "\" size=\"x-large\">" + t + "</span>";
            pango_layout_set_markup(pango_layout,markup.c_str(),-1);
            pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
            cairo_move_to(cr, atom.x * scale_factor + x_offset - layout_width/2.f, atom.y * scale_factor + y_offset - layout_height/2.f);
            pango_cairo_show_layout(cr, pango_layout);
        };

        
        if (atom.symbol == "C" || atom.symbol == "H") {
            process_highlight();
            // Ignore drawing Carbon and hydrogen now.
            // This should probably be computed at lowering time in the future.
            g_warning_once("TODO: Implement drawing atoms correctly");
        } else {
            render_white_background();
            process_highlight();
            render_text(atom.symbol,atom.color,atom.highlighted);
        }
    }
    cairo_destroy(cr);
}

CanvasMolecule::CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) {
    this->rdkit_molecule = std::move(rdkit_mol);
    this->lower_from_rdkit();
    this->x_canvas_size_adjustment = 0;
    this->y_canvas_size_adjustment = 0;
    this->x_canvas_translation = 0;
    this->y_canvas_translation = 0;
}

void CanvasMolecule::update_source_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) {
    this->rdkit_molecule = rdkit_mol;
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
            throw std::runtime_error("An aromatic bond remained after kekulization!");
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

CanvasMolecule::BondGeometry CanvasMolecule::bond_geometry_from_rdkit(RDKit::Bond::BondDir) {
    throw std::runtime_error("todo: bond_geometry_from_rdkit");
}

RDKit::Bond::BondType CanvasMolecule::bond_type_to_rdkit(CanvasMolecule::BondType ty) noexcept {
    switch (ty) {
        case BondType::Single:{
            return RDKit::Bond::SINGLE;
        }
        case BondType::Double:{
            return RDKit::Bond::DOUBLE;
        }
        case BondType::Triple:{
            return RDKit::Bond::TRIPLE;
        }
    }
}

void CanvasMolecule::lower_from_rdkit() {
    // 1. Clear what we have
    this->atoms.clear();
    this->bonds.clear();

    // 2. Do the lowering

    // 2.0 Kekulize
    RDKit::MolOps::Kekulize( *this->rdkit_molecule );
    // 2.1 Get geometry info

    // Maps atom indices to 2D points
    RDGeom::INT_POINT2D_MAP coordinate_map;

    // The following code is based on RDKit documentation.
    // In fact, it's almost a copy-paste.
    // I'm not exactly sure if what I did here is the most optimal solution
    // but it seems to work.
    RDDepict::compute2DCoords(*this->rdkit_molecule,nullptr,true,true);


    RDKit::MatchVectType matchVect;
    if(! RDKit::SubstructMatch(*this->rdkit_molecule, *this->rdkit_molecule, matchVect ) ) {
        throw std::runtime_error("SubstractMatch failed.");
    }
    RDKit::Conformer &conf = this->rdkit_molecule->getConformer();
    for(auto mv: matchVect) {
        RDGeom::Point3D pt3 = conf.getAtomPos( mv.first );
        RDGeom::Point2D pt2( pt3.x , pt3.y );
        coordinate_map[mv.second] = pt2;
    }
    // what is going on here?
    // That doesn't seem to change much
    // RDDepict::compute2DCoords( *this->rdkit_molecule, &coordinate_map, true, true);


    /// Used to avoid duplicating bonds
    std::set<unsigned int> processed_atoms_indices;

    if(coordinate_map.empty()) {
        throw std::runtime_error("RDKit coordinate mapping is empty");
    }
    // 2.2 Process atoms and compute bonds
    for(const auto& [atom_idx,plane_point]: coordinate_map) {
        const auto* rdkit_atom = this->rdkit_molecule->getAtomWithIdx(atom_idx);
        auto canvas_atom = CanvasMolecule::Atom();
        canvas_atom.color = atom_color_from_rdkit(rdkit_atom);
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
            canvas_bond.type = bond_type_from_rdkit(bond_ptr->getBondType());
            // todo: add support for geometry
            canvas_bond.geometry = BondGeometry::Flat;

            this->bonds.push_back(std::move(canvas_bond));
        }

        this->atoms.push_back(std::move(canvas_atom));
        // Mark the atom as processed
        processed_atoms_indices.insert(atom_idx);
    }
    std::sort(this->atoms.begin(),this->atoms.end(),[](const auto& lhs, const auto& rhs){
        return lhs.idx < rhs.idx;
    });
    // Reverse kekulization on the original molecule after lowering.
    RDKit::MolOps::sanitizeMol(*this->rdkit_molecule);
    
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