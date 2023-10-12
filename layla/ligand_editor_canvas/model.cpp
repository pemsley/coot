/* layla/ligand_editor_canvas/model.cpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "model.hpp"
#include <exception>
#include <iterator>
#include <memory>
#include <numeric>
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
#include <cctype>
#include <complex>

using namespace coot::ligand_editor_canvas;

const float CanvasMolecule::ATOM_HITBOX_RADIUS = 12.f;
const float CanvasMolecule::BOND_DISTANCE_BOUNDARY = 10.f;
const float CanvasMolecule::BASE_SCALE_FACTOR = 30.f;
const float CanvasMolecule::BOND_LINE_SEPARATION = 0.3f;
const float CanvasMolecule::CENTERED_DOUBLE_BOND_LINE_SEPARATION = 0.2f;
// 10 degrees
const float CanvasMolecule::GEOMETRY_BOND_SPREAD_ANGLE = M_PI/18.f;
const float CanvasMolecule::WAVY_BOND_ARC_LENGTH = 0.25f;
const float CanvasMolecule::GEOMETRY_BOND_DASH_SEPARATION = 0.2f;

const char* coot::ligand_editor_canvas::display_mode_to_string(DisplayMode mode) noexcept {
    switch (mode) {
        case DisplayMode::AtomIndices:{
            return "Atom Indices";
        }
        case DisplayMode::AtomNames:{
            return "Atom Names";
        }
        default:
        case DisplayMode::Standard:{
            return "Standard";
        }
    }
}
std::optional<DisplayMode> coot::ligand_editor_canvas::display_mode_from_string(const char* value_raw) noexcept {
    std::string value(value_raw);
    if(value == "Standard") {
        return DisplayMode::Standard;
    } else if (value == "Atom Indices") {
        return DisplayMode::AtomIndices;
    } else if (value == "Atom Names") {
        return DisplayMode::AtomNames;
    } else {
        return std::nullopt;
    }
}

float CanvasMolecule::get_scale() const noexcept {
    return BASE_SCALE_FACTOR * this->canvas_scale;
}

void CanvasMolecule::set_canvas_scale(float scale) {
    this->canvas_scale = scale;
}

CanvasMolecule::MaybeAtomOrBond CanvasMolecule::resolve_click(int x, int y) const noexcept {
    float scale = this->get_scale();
    auto x_offset = this->x_canvas_translation * scale;
    auto y_offset = this->y_canvas_translation * scale;
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
        float first_atom_x = bond->first_atom_x * scale + x_offset;
        float first_atom_y = bond->first_atom_y * scale + y_offset;
        float second_atom_x = bond->second_atom_x * scale + x_offset;
        float second_atom_y = bond->second_atom_y * scale + y_offset;
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
            return *bond;
        }
    }
    // nothing matches the click
    return std::nullopt;
}

void CanvasMolecule::apply_canvas_translation(int delta_x, int delta_y) noexcept {
    float scale = this->get_scale();
    this->x_canvas_translation += (float) delta_x / scale;
    this->y_canvas_translation += (float) delta_y / scale;
}

std::pair<float,float> CanvasMolecule::get_on_screen_coords(float x, float y) const noexcept {
    float scale = this->get_scale();
    auto x_offset = this->x_canvas_translation * scale;
    auto y_offset = this->y_canvas_translation * scale;
    return std::make_pair(x * scale + x_offset, y * scale + y_offset);
}

std::optional<std::pair<float,float>> CanvasMolecule::get_on_screen_coords_of_atom(unsigned int atom_idx) const noexcept {
    if(this->atoms.size() <= atom_idx) {
        return std::nullopt;
    }
    const Atom& a = this->atoms[atom_idx];
    return this->get_on_screen_coords(a.x, a.y);
}

graphene_rect_t CanvasMolecule::get_on_screen_bounding_rect() const noexcept {
    float scale = this->get_scale();
    auto x_offset = this->x_canvas_translation * scale;
    auto y_offset = this->y_canvas_translation * scale;
    graphene_rect_t ret;
    ret.origin.x = this->bounding_atom_coords.first.x * scale + x_offset;
    ret.origin.y = this->bounding_atom_coords.first.y * scale + y_offset;
    ret.size.width = (this->bounding_atom_coords.second.x - this->bounding_atom_coords.first.x) * scale + x_offset;
    ret.size.height = (this->bounding_atom_coords.second.y - this->bounding_atom_coords.first.y) * scale + y_offset;
    return ret;
}


void CanvasMolecule::perform_flip(FlipMode flip_mode) {
    for(auto& atom: this->cached_atom_coordinate_map.value()) {
        if(flip_mode == FlipMode::Horizontal) {
            atom.second.x *= -1;
        } else {
            atom.second.y *= -1;
        }
    }
}

void CanvasMolecule::rotate_by_angle(double radians) {
    std::complex<double> rotation_mult(std::cos(radians),std::sin(radians));
    for(auto& atom: this->cached_atom_coordinate_map.value()) {
        std::complex<double> atom_cn(atom.second.x,atom.second.y);
        atom_cn *= rotation_mult;
        atom.second.x = atom_cn.real();
        atom.second.y = atom_cn.imag();
    }
}

std::tuple<float,float,float> CanvasMolecule::atom_color_to_rgb(CanvasMolecule::AtomColor color) noexcept {
    switch (color) {
        case AtomColor::Green:{
            return std::make_tuple(0.0,0.75,0.0);
        }
        case AtomColor::Blue:{
            return std::make_tuple(0.0,0.0,1.0);
        }
        case AtomColor::Red:{
            return std::make_tuple(1.0,0.0,0.0);
        }
        case AtomColor::Brown:{
            return std::make_tuple(0.5,0.5,0.0);
        }
        case AtomColor::DarkRed:{
            return std::make_tuple(0.5,0.0,0.0);
        }
        case AtomColor::Orange:{
            return std::make_tuple(1.0,0.5,0.0);
        }
        case AtomColor::DarkBlue:{
            return std::make_tuple(0.0,0.0,0.5);
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
            return "#00C000";
        }
        case AtomColor::Blue:{
            return "#0000FF";
        }
        case AtomColor::Red:{
            return "#FF0000";
        }
        case AtomColor::Brown:{
            return "#808000";
        }
        case AtomColor::DarkRed:{
            return "#800000";
        }
        case AtomColor::Orange:{
            return "#FF8000";
        }
        case AtomColor::DarkBlue:{
            return "#000080";
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
        // Phosphorus
        case 15: {
            return AtomColor::Orange;
        }
        // Sulphur
        case 16: {
            return AtomColor::Brown;
        }
        // Flourine
        case 9:
        // Chlorine
        case 17: {
            return AtomColor::Green;
        }
        // Bromine
        case 35: {
            return AtomColor::DarkRed;
        }
        // Iodine
        case 53: {
            return AtomColor::DarkBlue;
        }
        // Carbon
        case 6:
        default: {
            return AtomColor::Black;
        }
    }
}

CanvasMolecule::Atom::Appendix::Appendix() noexcept 
    :charge(0),
    reversed(false) {
    
}

std::pair<float,float> CanvasMolecule::Bond::get_vector() const noexcept {
    float bond_vector_x = second_atom_x - first_atom_x;
    float bond_vector_y = second_atom_y - first_atom_y;
    return std::make_pair(bond_vector_x,bond_vector_y);
}

std::pair<float,float> CanvasMolecule::Bond::get_versor() const noexcept {
    auto [bond_vector_x, bond_vector_y] = this->get_vector();
    float bond_vector_len = std::sqrt(std::pow(bond_vector_x,2.f) + std::pow(bond_vector_y,2.f));
    if (bond_vector_len == 0) {
        return std::make_pair(0.f, 0.f);
    } else {
        return std::make_pair(bond_vector_x/bond_vector_len,bond_vector_y/bond_vector_len);
    }
}

std::pair<float,float> CanvasMolecule::Bond::get_perpendicular_versor() const noexcept {
    auto [x, y] = this->get_versor();
    return std::make_pair(-y,x);
}

float CanvasMolecule::Bond::get_length() const noexcept {
    auto [bond_vector_x, bond_vector_y] = this->get_vector();
    return std::sqrt(std::pow(bond_vector_x,2.f) + std::pow(bond_vector_y,2.f));
}

void CanvasMolecule::draw(cairo_t* cr, PangoLayout* pango_layout, DisplayMode display_mode) const noexcept {
    auto scale_factor = this->get_scale();
    auto x_offset = scale_factor * this->x_canvas_translation;
    auto y_offset = scale_factor * this->y_canvas_translation;
    
    cairo_set_line_width(cr, 0.5);
    // Used to truncate bonds not to cover atoms
    std::map<unsigned int,graphene_rect_t> atom_idx_to_canvas_rect;
    for(const auto& atom: atoms) {

        auto process_highlight = [&,cr,x_offset,y_offset,scale_factor](){
            if(atom.highlighted) {
                //cairo_move_to(cr, atom.x * scale_factor + x_offset + ATOM_HITBOX_RADIUS, atom.y * scale_factor + y_offset);
                cairo_new_sub_path(cr);
                cairo_set_source_rgb(cr, 0.0, 1.0, 0.5);
                cairo_arc(cr, atom.x * scale_factor + x_offset, atom.y * scale_factor + y_offset,ATOM_HITBOX_RADIUS,0,M_PI * 2.0);
                cairo_stroke_preserve(cr);
                cairo_set_source_rgba(cr, 0.0, 1.0, 0.5, 0.5);
                cairo_fill(cr);
            }
        };

        /// Returns the markup string + info if the appendix is reversed
        auto process_appendix = [&](const std::string& symbol, const std::optional<Atom::Appendix>& appendix) -> std::tuple<std::string,bool> {
            std::string ret = symbol;
            bool reversed = false;
            if(appendix.has_value()) {
                const auto& ap = appendix.value();
                //ret += "<span>";
                std::string ap_root;
                for(auto i = ap.superatoms.begin(); i != ap.superatoms.end(); i++) {
                    if(std::isdigit(*i)) {
                        ap_root += "<sub>";
                        ap_root.push_back(*i);
                        ap_root += "</sub>";
                    } else {
                        ap_root.push_back(*i);
                    }
                }
                if (ap.reversed) {
                    ret = ap_root + ret;
                    reversed = true;
                } else {
                    ret += ap_root;
                }
                //ret += "</span>";
                if(ap.charge != 0) {
                    // The string below begins with 
                    // the invisible U+200B unicode character.
                    // This is a workaround for what's likely 
                    // a bug in pango font rendering engine.
                    // Without it, the superscript is relative 
                    // to the subscript (atom count)
                    // instead of the atom's symbol
                    ret += "​<sup>";
                    unsigned int charge_no_sign = std::abs(ap.charge);
                    if(charge_no_sign > 1) {
                        ret += std::to_string(charge_no_sign);
                    }
                    ret.push_back(ap.charge > 0 ? '+' : '-');
                    ret += "</sup>";
                }
            }
            return std::make_tuple(ret,reversed);
        };

        // Returns a pair of atom index and bonding rect
        auto render_atom = [&](const Atom& atom, DisplayMode render_mode = DisplayMode::Standard) -> std::pair<unsigned int,graphene_rect_t> {
            // pre-process text
            auto [r,g,b] = atom_color_to_rgb(atom.color);
            const std::string color_str = atom_color_to_html(atom.color);
            const std::string weight_str = atom.highlighted ? "bold" : "normal";
            const std::string size_str = render_mode != DisplayMode::AtomIndices ? "x-large" : "medium";
            const std::string markup_beginning = "<span color=\"" + color_str + "\" weight=\"" + weight_str + "\" size=\"" + size_str + "\">";
            const std::string markup_ending = "</span>";

            bool reversed = false;
            // Markup for the whole thing - includes symbol, appendix, index etc.
            std::string markup;
            // Markup for the atom symbol - solely.
            // This allows us to measure the size of the atom's symbol 
            // and then properly center/align the whole text
            std::string markup_no_appendix;

            switch (render_mode) {
                case DisplayMode::AtomIndices: {
                    markup_no_appendix = markup_beginning + atom.symbol + markup_ending;
                    markup = markup_beginning + atom.symbol + ":" + std::to_string(atom.idx) + markup_ending;
                    break;
                }
                case DisplayMode::AtomNames: {
                    if(atom.name.has_value()) {
                        std::string atom_name = atom.name.value();
                        markup_no_appendix = markup_beginning + atom_name + markup_ending;
                        markup = markup_no_appendix;
                        break;
                    } 
                    // break;
                    // We want to fall back to the standard case if the atom has no name.
                }
                default:
                case DisplayMode::Standard: {
                    auto [raw_markup,p_reversed] = process_appendix(atom.symbol,atom.appendix);
                    reversed = p_reversed;
                    markup_no_appendix = markup_beginning + atom.symbol + markup_ending;
                    markup = markup_beginning + raw_markup + markup_ending;
                    break;
                }
            }

            pango_layout_set_markup(pango_layout,markup_no_appendix.c_str(),-1);
            // Used to make the texts centered where they should be (appendix).
            int layout_height_no_ap, layout_width_no_ap;
            // Measure the size of the "main" atom, without "appendix"
            pango_layout_get_pixel_size(pango_layout,&layout_width_no_ap,&layout_height_no_ap);
            pango_layout_set_markup(pango_layout,markup.c_str(),-1);
            // Used to make the texts centered where they should be.
            int layout_width, layout_height;
            // Measure full size of the text
            pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);

            // todo: get rid of this '5' magic number - figure out what's wrong
            int layout_x_offset = reversed ? layout_width - layout_height_no_ap / 2.f + 5 : layout_width_no_ap / 2.f;
            double origin_x = atom.x * scale_factor + x_offset - layout_x_offset;
            double origin_y = atom.y * scale_factor + y_offset - layout_height_no_ap/2.f;

            graphene_rect_t rect;
            rect.origin.x = origin_x;
            // Magic number. This should be removed.
            // Workaround for pango giving us too high layout size.
            const float layout_to_high = 3.f;
            rect.origin.y = origin_y + layout_to_high;
            rect.size.width = layout_width;
            rect.size.height = layout_height - layout_to_high * 2.f;

            // highlight
            process_highlight();
            // text
            cairo_move_to(cr, origin_x, origin_y);
            pango_cairo_show_layout(cr, pango_layout);

            return std::make_pair(atom.idx,rect);
        };

        switch (display_mode) {
            case DisplayMode::AtomIndices: {
                atom_idx_to_canvas_rect.emplace(render_atom(atom,DisplayMode::AtomIndices));
                break;
            }
            case DisplayMode::AtomNames: {
                if(atom.name.has_value()) {
                    atom_idx_to_canvas_rect.emplace(render_atom(atom,DisplayMode::AtomNames));
                    break;
                }
                // We want to fall back to the standard case if the atom has no name
                // break;
            }
            default:
            case DisplayMode::Standard: {
                if(atom.symbol == "C") {
                    if(atom.appendix.has_value()) {
                        atom_idx_to_canvas_rect.emplace(render_atom(atom));
                    } else {
                        process_highlight();
                    }
                } else {
                    atom_idx_to_canvas_rect.emplace(render_atom(atom));
                }
                break;
            }
        }
    }

    for(const auto& bond: bonds) {
        if(bond->highlighted) {
            cairo_set_line_width(cr, 4.0);
            cairo_set_source_rgb(cr, 0.0, 1.0, 0.5);
        } else {
            cairo_set_line_width(cr, 2.0);
            cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        }

        // Returns on-screen bond coordinates
        // cropped not to overlap with atoms' symbols and appendices.
        // Accepts on-screen coordinates of bond atoms.
        auto cropped_bond_coords = [&](const graphene_point_t& first_atom, unsigned int first_atom_idx, const graphene_point_t& second_atom, unsigned int second_atom_idx) -> std::pair<graphene_point_t,graphene_point_t> {

            // We pass in the bond vector so that it always points "away" from the point
            auto crop_line_against_rect = [](const graphene_rect_t& rect, float bond_vec_x, float bond_vec_y, const graphene_point_t& point){
                graphene_point_t ret = point;
                if(rect.origin.x < point.x && rect.origin.x + rect.size.width >= point.x
                && rect.origin.y < point.y && rect.origin.y + rect.size.height >= point.y) { // inside rectangle
                    // We need to use the point of intersection formula
                    // to find the intersection point between the bond and the edges.
                    // We solve against two edges and then pick the solution lying the closest to 'point'.

                    // Line from bond
                    graphene_point_t point2 = point;
                    point2.x += bond_vec_x;
                    point2.y += bond_vec_y;
                    const float bond_a_quot = (point.x - point2.x);
                    const float bond_a = bond_a_quot == 0 ? -point.x : (point.y - point2.y) / bond_a_quot;
                    const float bond_b = -1;
                    const float bond_c = point.y - bond_a * point.x;

                    // Pick line from the relevant vertical edge
                    float vert_egde_c;
                    const float vert_edge_a = 1;
                    const float vert_edge_b = 0;
                    if(bond_vec_x > 0) { // solve against right side
                        vert_egde_c = -(rect.origin.x + rect.size.width);
                    } else { // solve against left side
                        vert_egde_c = -rect.origin.x;
                    }

                    // Pick line from the relevant horizontal edge
                    float horiz_egde_c;
                    const float horiz_edge_a = 0;
                    const float horiz_edge_b = 1;
                    if(bond_vec_y > 0) { // solve against the bottom
                        horiz_egde_c = -(rect.origin.y + rect.size.height);
                    } else { // solve against the top
                        horiz_egde_c = -rect.origin.y;
                    }
                    graphene_point_t vert_edge_solution;
                    vert_edge_solution.x = (bond_b * vert_egde_c - vert_edge_b * bond_c)/(bond_a * vert_edge_b - vert_edge_a * bond_b);
                    vert_edge_solution.y = (vert_edge_a * bond_c - bond_a * vert_egde_c)/(bond_a * vert_edge_b - vert_edge_a * bond_b);
                    graphene_point_t horiz_edge_solution;
                    horiz_edge_solution.x = (bond_b * horiz_egde_c - horiz_edge_b * bond_c)/(bond_a * horiz_edge_b - horiz_edge_a * bond_b);
                    horiz_edge_solution.y = (horiz_edge_a * bond_c - bond_a * horiz_egde_c)/(bond_a * horiz_edge_b - horiz_edge_a * bond_b);

                    // if vert_edge_solution is farther from the point than horiz_edge_solution
                    if(std::pow(vert_edge_solution.x - point.x, 2.f) + std::pow(vert_edge_solution.y - point.y, 2.f) 
                     > std::pow(horiz_edge_solution.x - point.x, 2.f) + std::pow(horiz_edge_solution.y - point.y, 2.f)) {
                        ret = horiz_edge_solution;
                    } else {
                        ret = vert_edge_solution;
                    }
                } 
                return ret;
            };

            graphene_point_t a = first_atom;
            graphene_point_t b = second_atom;

            float bond_vec_x = second_atom.x - first_atom.x;
            float bond_vec_y = second_atom.y - first_atom.y;

            auto first_rect_iter = atom_idx_to_canvas_rect.find(first_atom_idx);
            if(first_rect_iter != atom_idx_to_canvas_rect.end()) {
                a = crop_line_against_rect(first_rect_iter->second, bond_vec_x, bond_vec_y, first_atom);
            }
            auto second_rect_iter = atom_idx_to_canvas_rect.find(second_atom_idx);
            if(second_rect_iter != atom_idx_to_canvas_rect.end()) {
                b = crop_line_against_rect(second_rect_iter->second, -bond_vec_x, -bond_vec_y, second_atom);
            }
            
            return std::make_pair(a,b);
        };

        auto draw_central_bond_line = [&](){
            graphene_point_t first_atom;
            first_atom.x = bond->first_atom_x * scale_factor + x_offset;
            first_atom.y = bond->first_atom_y * scale_factor + y_offset;

            graphene_point_t second_atom;
            second_atom.x = bond->second_atom_x * scale_factor + x_offset;
            second_atom.y = bond->second_atom_y * scale_factor + y_offset;

            auto [first,second] = cropped_bond_coords(first_atom,bond->first_atom_idx,second_atom,bond->second_atom_idx);
            
            cairo_move_to(cr, first.x, first.y);
            cairo_line_to(cr, second.x, second.y);
            cairo_stroke(cr);
        };

        if(bond->geometry != BondGeometry::Flat && bond->type == BondType::Single) {

            auto draw_straight_wedge = [&](bool reversed){
                graphene_point_t origin;
                origin.x = reversed ? bond->first_atom_x : bond->second_atom_x;
                origin.y = reversed ? bond->first_atom_y : bond->second_atom_y;
                auto origin_idx = reversed ? bond->first_atom_idx : bond->second_atom_idx;

                graphene_point_t target;
                target.x = reversed ? bond->second_atom_x : bond->first_atom_x;
                target.y = reversed ? bond->second_atom_y : bond->first_atom_y;
                auto target_idx = reversed ? bond->second_atom_idx : bond->first_atom_idx;

                origin.x *= scale_factor;
                origin.y *= scale_factor;
                target.x *= scale_factor;
                target.y *= scale_factor;

                origin.x += x_offset;
                origin.y += y_offset;
                target.x += x_offset;
                target.y += y_offset;

                auto [origin_cropped,target_cropped] = cropped_bond_coords(origin, origin_idx, target, target_idx);
                
                auto [pv_x,pv_y] = bond->get_perpendicular_versor();
                auto cropped_bond_len = std::sqrt(std::pow(target_cropped.x - origin_cropped.x, 2.f) + std::pow(target_cropped.y - origin_cropped.y, 2.f));
                auto v_x = pv_x * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;
                auto v_y = pv_y * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;

                cairo_new_path(cr);
                cairo_move_to(cr, origin_cropped.x, origin_cropped.y);
                cairo_line_to(cr, target_cropped.x + v_x, target_cropped.y + v_y);
                cairo_stroke_preserve(cr);

                cairo_line_to(cr, target_cropped.x - v_x, target_cropped.y - v_y);
                cairo_stroke_preserve(cr);

                cairo_line_to(cr, origin_cropped.x, origin_cropped.y);
                cairo_stroke_preserve(cr);

                cairo_close_path(cr);
                cairo_fill(cr);
            };
            auto draw_straight_dashed_bond = [&](bool reversed){
                graphene_point_t origin;
                origin.x = reversed ? bond->first_atom_x : bond->second_atom_x;
                origin.y = reversed ? bond->first_atom_y : bond->second_atom_y;
                auto origin_idx = reversed ? bond->first_atom_idx : bond->second_atom_idx;

                graphene_point_t target;
                target.x = reversed ? bond->second_atom_x : bond->first_atom_x;
                target.y = reversed ? bond->second_atom_y : bond->first_atom_y;
                auto target_idx = reversed ? bond->second_atom_idx : bond->first_atom_idx;

                origin.x *= scale_factor;
                origin.y *= scale_factor;
                target.x *= scale_factor;
                target.y *= scale_factor;

                origin.x += x_offset;
                origin.y += y_offset;
                target.x += x_offset;
                target.y += y_offset;

                auto [current,target_cropped] = cropped_bond_coords(origin, origin_idx, target, target_idx);
                
                auto [pv_x,pv_y] = bond->get_perpendicular_versor();
                auto cropped_bond_len = std::sqrt(std::pow(target_cropped.x - current.x, 2.f) + std::pow(target_cropped.y - current.y, 2.f));

                float dashes = cropped_bond_len / (GEOMETRY_BOND_DASH_SEPARATION * scale_factor);
                unsigned int full_dashes = std::floor(dashes);

                float step_x = (target_cropped.x - current.x) / dashes;
                float step_y = (target_cropped.y - current.y) / dashes;

                auto v_x = pv_x * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;
                auto v_y = pv_y * std::sin(GEOMETRY_BOND_SPREAD_ANGLE / 2.f) * cropped_bond_len;
                
                for(unsigned int i = 0; i <= full_dashes; i++) {
                    float spread_multiplier = (float) i / dashes;
                    cairo_move_to(cr, current.x - v_x * spread_multiplier, current.y - v_y * spread_multiplier);
                    cairo_line_to(cr, current.x + v_x * spread_multiplier, current.y + v_y * spread_multiplier);
                    cairo_stroke(cr);
                    current.x += step_x;
                    current.y += step_y;
                }
            };
            switch (bond->geometry) {
                default:
                case BondGeometry::Unspecified:{
                    graphene_point_t first_atom;
                    first_atom.x = bond->first_atom_x * scale_factor + x_offset;
                    first_atom.y = bond->first_atom_y * scale_factor + y_offset;

                    graphene_point_t second_atom;
                    second_atom.x = bond->second_atom_x * scale_factor + x_offset;
                    second_atom.y = bond->second_atom_y * scale_factor + y_offset;

                    auto [first,second] = cropped_bond_coords(first_atom,bond->first_atom_idx,second_atom,bond->second_atom_idx);
                    auto full_vec_x = second.x - first.x;
                    auto full_vec_y = second.y - first.y;

                    const float wave_arc_radius = WAVY_BOND_ARC_LENGTH * scale_factor / 2.f;

                    // The angle at which the bond points
                    float base_angle = std::atan(full_vec_y / full_vec_x);
                    float arcs_count = std::sqrt(std::pow(full_vec_x,2.f) + std::pow(full_vec_y,2.f)) / (WAVY_BOND_ARC_LENGTH * scale_factor);
                    unsigned int rounded_arcs_count = std::floor(arcs_count);
                    float step_x = full_vec_x / arcs_count;
                    float step_y = full_vec_y / arcs_count;
                    float current_x = first.x + step_x / 2.f;
                    float current_y = first.y + step_y / 2.f;
                    // It seems that for positive base_angle, 
                    // 'true' is counter-clockwise
                    // and 'false is clockwise.
                    // For negative base_angle, it's the opposite.
                    bool arc_direction = true;
                    // for debugging
                    // float l_angle_one = 0, l_angle_two = 0;

                    // Two core angles for semi-circles
                    float p1 = base_angle;
                    float p2 = base_angle - M_PI;

                    float angle_one, angle_two;

                    for (unsigned int i = 0; i < rounded_arcs_count; i++) {
                        float next_x = current_x + step_x;
                        float next_y = current_y + step_y;
                        
                        if(arc_direction) {
                            angle_one = p1;
                            angle_two = p2;
                        } else {
                            angle_one = p2;
                            angle_two = p1;
                        }
                        // l_angle_one = angle_one;
                        // l_angle_two = angle_two;
                        cairo_new_sub_path(cr);
                        cairo_arc(cr, current_x, current_y, wave_arc_radius, angle_one, angle_two);
                        cairo_stroke(cr);
                        current_x = next_x;
                        current_y = next_y;
                        arc_direction = !arc_direction;
                    }
                    // Final part of the path. Truncated arc.
                    float partial_arc_proportion = arcs_count - (float) rounded_arcs_count;
                    float arccos_arg = 1.f - (partial_arc_proportion / WAVY_BOND_ARC_LENGTH / 2.f);
                    
                    // This is the angle for the final arc.
                    float theta = std::acos(arccos_arg);

                    // The magic behind 'step_x > 0'
                    // is a bit of a mystery (derived empirically).
                    // There's certain correlation with the sign of base_angle
                    // but that's not the whole story.
                    // What matters is that it allows for differentiating
                    // between various cases, with angles from different quarters.
                    float starting_angle = step_x > 0 ? p2 : p1;
                    if(arc_direction) {
                        if(step_x > 0) {
                            angle_one = starting_angle - theta;
                            angle_two = starting_angle;
                        } else {
                            angle_one = starting_angle;
                            angle_two = starting_angle + theta;
                        }
                    } else {
                        if(step_x > 0) {
                            angle_one = starting_angle;
                            angle_two = starting_angle + theta;
                        } else {
                            angle_one = starting_angle - theta;
                            angle_two = starting_angle;
                        }
                    }

                    // debugging stuff

                    // std::string case_info;
                    // if(theta > M_PI_2) {
                    //     case_info += "T";
                    // } else {
                    //     case_info += "G";
                    // }
                    // if(base_angle > 0) {
                    //     case_info += "A";
                    // } else {
                    //     case_info += "E";
                    // }
                    // if(arc_direction) {
                    //     case_info += "K";
                    // } else {
                    //     case_info += "V";
                    // }
                    // case_info += angle_two - angle_one > 0 ? "O" : "U";
                    // if(step_x > 0) {
                    //     case_info += "Z";
                    // } else {
                    //     case_info += "V";
                    // }
                    // g_debug(
                    //     "theta=%f, base_angle=%f a1=%f, a2=%f a2-a1=%f abs(a2-a1)=%f direction=%s p1=%f, p2=%f case_codename=%s",
                    //     theta / M_PI * 180.f,
                    //     base_angle / M_PI * 180.f,
                    //     angle_one / M_PI * 180.f,
                    //     angle_two / M_PI * 180.f,
                    //     (angle_two - angle_one) / M_PI * 180.f,
                    //     std::fabs(angle_two - angle_one) / M_PI * 180.f,
                    //     arc_direction ? "true" : "false",
                    //     l_angle_one / M_PI * 180.f,
                    //     l_angle_two / M_PI * 180.f,
                    //     case_info.c_str()
                    // );
                    cairo_new_sub_path(cr);
                    cairo_arc(cr, current_x, current_y, wave_arc_radius, angle_one, angle_two);
                    cairo_stroke(cr);
                    break;
                }
                case BondGeometry::WedgeTowardsFirst:{
                    draw_straight_wedge(true);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
                case BondGeometry::WedgeTowardsSecond:{
                    draw_straight_wedge(false);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
                case BondGeometry::DashedTowardsFirst:{
                    draw_straight_dashed_bond(true);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
                case BondGeometry::DashedTowardsSecond:{
                    draw_straight_dashed_bond(false);
                    g_warning_once("todo: rendering bond geometry in rings");
                    break;
                }
            }
        } else {
            auto draw_side_bond_line = [&](
                bool addOrSub, 
                std::optional<float> first_shortening_proportion, 
                std::optional<float> second_shortening_proportion
                ) {

                auto [pv_x,pv_y] = bond->get_perpendicular_versor();
                if (!addOrSub) { // change sign of the versor
                    pv_x *= -1.f;
                    pv_y *= -1.f;
                }
                // Convert the versor to a vector of the desired length
                pv_x *= BOND_LINE_SEPARATION;
                pv_y *= BOND_LINE_SEPARATION;

                auto [bond_vec_x, bond_vec_y] = bond->get_vector();
                auto first_x = bond->first_atom_x;
                auto second_x = bond->second_atom_x;
                auto first_y = bond->first_atom_y;
                auto second_y = bond->second_atom_y;

                if(first_shortening_proportion.has_value()) {
                    first_x += first_shortening_proportion.value() * bond_vec_x;
                    first_y += first_shortening_proportion.value() * bond_vec_y;
                }
                if(second_shortening_proportion.has_value()) {
                    second_x -= second_shortening_proportion.value() * bond_vec_x;
                    second_y -= second_shortening_proportion.value() * bond_vec_y;
                }

                // Points a and b represent the off-center bond before cropping.
                graphene_point_t a;
                a.x = (first_x + pv_x) * scale_factor + x_offset;
                a.y = (first_y + pv_y) * scale_factor + y_offset;
                graphene_point_t b;
                b.x = (second_x + pv_x) * scale_factor + x_offset;
                b.y = (second_y + pv_y) * scale_factor + y_offset;

                // We need to make sure that the off-center bond 
                // after cropping is not going to be longer than the center bond
                graphene_point_t first_atom_centered;
                first_atom_centered.x = first_x * scale_factor + x_offset;
                first_atom_centered.y = first_y * scale_factor + y_offset;

                graphene_point_t second_atom_centered;
                second_atom_centered.x = second_x * scale_factor + x_offset;
                second_atom_centered.y = second_y * scale_factor + y_offset;

                /// Centered bond cropped
                auto [first_c,second_c] = cropped_bond_coords(
                    first_atom_centered,
                    bond->first_atom_idx,
                    second_atom_centered,
                    bond->second_atom_idx
                );
                // Now we offset the center bond after cropping
                first_c.x += pv_x * scale_factor;
                first_c.y += pv_y * scale_factor;
                second_c.x += pv_x * scale_factor;
                second_c.y += pv_y * scale_factor;

                // Points a_cropped and b_cropped represent the off-center bond after cropping.
                auto [a_cropped,b_cropped] = cropped_bond_coords(a,bond->first_atom_idx,b,bond->second_atom_idx);

                // Now we need to make sure that the off-center bond 
                // after cropping is not going to be longer than the center bond
                if(bond_vec_x > 0) {
                    // The beginning is shorter for the centered bond
                    if(first_c.x > a_cropped.x) {
                        a_cropped = first_c;
                    }
                    // The end is shorter for the centered bond
                    if(second_c.x < b_cropped.x) {
                        b_cropped = second_c;
                    }
                } else {
                    // The beginning is shorter for the centered bond
                    if(first_c.x < a_cropped.x) {
                        a_cropped = first_c;
                    }
                    // The end is shorter for the centered bond
                    if(second_c.x > b_cropped.x) {
                        b_cropped = second_c;
                    }
                }
                if(bond_vec_y > 0) {
                    // The beginning is shorter for the centered bond
                    if(first_c.y > a_cropped.y) {
                        a_cropped = first_c;
                    }
                    // The end is shorter for the centered bond
                    if(second_c.y < b_cropped.y) {
                        b_cropped = second_c;
                    }
                } else {
                    // The beginning is shorter for the centered bond
                    if(first_c.y < a_cropped.y) {
                        a_cropped = first_c;
                    }
                    // The end is shorter for the centered bond
                    if(second_c.y > b_cropped.y) {
                        b_cropped = second_c;
                    }
                }

                cairo_move_to(cr, a_cropped.x, a_cropped.y);
                cairo_line_to(cr, b_cropped.x, b_cropped.y);
                cairo_stroke(cr);
            };

            switch(bond->type) {
                case BondType::Double:{
                    DoubleBondDrawingDirection direction = bond->bond_drawing_direction.has_value() ? bond->bond_drawing_direction.value() : DoubleBondDrawingDirection::Primary;
                    bool direction_as_bool = true;

                    switch (direction) { 
                        case DoubleBondDrawingDirection::Secondary:{
                            direction_as_bool = false;
                            // no break here.
                        }
                        case DoubleBondDrawingDirection::Primary:{
                            draw_central_bond_line();
                            draw_side_bond_line(
                                direction_as_bool,
                                bond->first_shortening_proportion,
                                bond->second_shortening_proportion
                            );
                            break;
                        }
                        case DoubleBondDrawingDirection::Centered:{
                            auto [pv_x,pv_y] = bond->get_perpendicular_versor();

                            // Convert the versor to a vector of the desired length
                            pv_x *= CENTERED_DOUBLE_BOND_LINE_SEPARATION / 2.f * scale_factor;
                            pv_y *= CENTERED_DOUBLE_BOND_LINE_SEPARATION / 2.f * scale_factor;

                            graphene_point_t first_atom;
                            first_atom.x = bond->first_atom_x * scale_factor + x_offset;
                            first_atom.y = bond->first_atom_y * scale_factor + y_offset;

                            graphene_point_t second_atom;
                            second_atom.x = bond->second_atom_x * scale_factor + x_offset;
                            second_atom.y = bond->second_atom_y * scale_factor + y_offset;

                            auto [first,second] = cropped_bond_coords(first_atom,bond->first_atom_idx,second_atom,bond->second_atom_idx);

                            cairo_move_to(cr, first.x + pv_x, first.y + pv_y);
                            cairo_line_to(cr, second.x + pv_x, second.y + pv_y);
                            cairo_stroke(cr);

                            cairo_move_to(cr, first.x - pv_x, first.y - pv_y);
                            cairo_line_to(cr, second.x - pv_x, second.y - pv_y);
                            cairo_stroke(cr);
                            break;
                        }
                    }
                    break;
                }
                case BondType::Triple:{
                    draw_central_bond_line();
                    g_warning_once("todo: Triple bonds might need truncating too.");
                    // "to the left"
                    draw_side_bond_line(false,std::nullopt,std::nullopt);
                    // "to the right"
                    draw_side_bond_line(true,std::nullopt,std::nullopt);
                    break;
                }
                default:
                case BondType::Single:{
                    draw_central_bond_line();
                    break;
                }
            }
        }
    }
}

CanvasMolecule::CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) {
    this->rdkit_molecule = std::move(rdkit_mol);
    this->cached_atom_coordinate_map = std::nullopt;
    this->lower_from_rdkit(true);
    this->x_canvas_translation = 0;
    this->y_canvas_translation = 0;
    this->bounding_atom_coords = std::make_pair(RDGeom::Point2D(0,0),RDGeom::Point2D(0,0));
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

CanvasMolecule::BondGeometry CanvasMolecule::bond_geometry_from_rdkit(RDKit::Bond::BondDir dir) noexcept {
    //Not handled here: EITHERDOUBLE
    switch (dir) {
        default:{
            g_warning("Unhandled RDKit bond geometry: %i! Falling back to flat.", dir);
        }
        case RDKit::Bond::NONE:{
            return BondGeometry::Flat;
        }
        case RDKit::Bond::UNKNOWN:{
            return BondGeometry::Unspecified;
        }
        case RDKit::Bond::BEGINWEDGE:{
            return BondGeometry::WedgeTowardsSecond;
        }
        case RDKit::Bond::BEGINDASH:{
            return BondGeometry::DashedTowardsSecond;
        }
        case RDKit::Bond::ENDDOWNRIGHT:{
            // todo: make sure that this makes sense
            return BondGeometry::WedgeTowardsFirst;
        }
        case RDKit::Bond::ENDUPRIGHT:{
            // todo: make sure that this makes sense
            return BondGeometry::DashedTowardsFirst;
        }
    }
}

CanvasMolecule::BondGeometry CanvasMolecule::cycle_bond_geometry(BondGeometry geom) noexcept {
    switch (geom) {
        default: 
        case BondGeometry::Flat: {
            return BondGeometry::WedgeTowardsFirst;
        }
        case BondGeometry::Unspecified: {
            return BondGeometry::Flat;
        }
        case BondGeometry::WedgeTowardsFirst: {
            return BondGeometry::WedgeTowardsSecond;
        }
        case BondGeometry::WedgeTowardsSecond: {
            return BondGeometry::DashedTowardsFirst;
        }
        case BondGeometry::DashedTowardsFirst: {
            return BondGeometry::DashedTowardsSecond;
        }
        case BondGeometry::DashedTowardsSecond: {
            return BondGeometry::Unspecified;
        }
    }
}

RDKit::Bond::BondDir CanvasMolecule::bond_geometry_to_rdkit(BondGeometry geom) noexcept {
    switch (geom) {
        default: 
        case BondGeometry::Flat: {
            return RDKit::Bond::BondDir::NONE;
        }
        case BondGeometry::Unspecified: {
            return RDKit::Bond::BondDir::UNKNOWN;
        }
        case BondGeometry::WedgeTowardsFirst: {
            return RDKit::Bond::BondDir::ENDDOWNRIGHT;
        }
        case BondGeometry::WedgeTowardsSecond: {
            return RDKit::Bond::BondDir::BEGINWEDGE;
        }
        case BondGeometry::DashedTowardsFirst: {
            return RDKit::Bond::BondDir::ENDUPRIGHT;
        }
        case BondGeometry::DashedTowardsSecond: {
            return RDKit::Bond::BondDir::BEGINDASH;
        }
    }
}


RDKit::Bond::BondType CanvasMolecule::bond_type_to_rdkit(CanvasMolecule::BondType ty) noexcept {
    switch (ty) {
        default:
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


RDGeom::INT_POINT2D_MAP CanvasMolecule::compute_molecule_geometry() const {
    // The following code is heavily based on RDKit documentation.

    const RDGeom::INT_POINT2D_MAP* previous_coordinate_map = nullptr;
    
    // Pruned coordinate map will only be instantiated if the last coordinate map
    // contains atoms which are no longer in the molecule (and thus need to be removed).
    std::unique_ptr<RDGeom::INT_POINT2D_MAP> pruned_previous_coordinate_map = nullptr;

    if (this->cached_atom_coordinate_map.has_value()) {
        // g_debug("Computing 2D coords using a reference");
        const RDGeom::INT_POINT2D_MAP& prev_coord_map_ref = this->cached_atom_coordinate_map.value();
        previous_coordinate_map = &prev_coord_map_ref;
        // We need to make sure that each atom in the last_atom_coordinate_map still exists in the molecule.
        // If it doesn't, RDDepict::compute2DCoords() throws an exception. We cannot let this happen.
        for(const auto [atom_idx,_point]: prev_coord_map_ref) {
            try {
                auto* _atom_ptr = this->rdkit_molecule->getAtomWithIdx(atom_idx);
                // All good.
            } catch(...) { // Atom does not exist
                g_info("Atom with id=%u missing in the molecule, but found in 2D coordinate reference. It will be omitted.",atom_idx);
                // If there's any atom which is no longer there, we need to copy the whole coordinate map
                // into pruned_previous_coordinate_map and then remove the atom there.
                if(! pruned_previous_coordinate_map) {
                    // Copy the original coordinate map
                    pruned_previous_coordinate_map = std::make_unique<RDGeom::INT_POINT2D_MAP>(this->cached_atom_coordinate_map.value());
                    // Also we need to update the previous_coordinate_map pointer to refer to our pruned map.
                    previous_coordinate_map = pruned_previous_coordinate_map.get();
                }
                pruned_previous_coordinate_map->erase(pruned_previous_coordinate_map->find(atom_idx));
            }
        }
    } else {
        g_info("Computing fresh 2D coords (without previous reference).");
    }

    try {
        RDDepict::compute2DCoords(*this->rdkit_molecule,previous_coordinate_map,true,true);
    } catch(std::exception& e) {
        throw std::runtime_error(std::string("Failed to compute 2D coords with RDKit! ")+e.what());
    }


    RDKit::MatchVectType matchVect;
    if(!RDKit::SubstructMatch(*this->rdkit_molecule, *this->rdkit_molecule, matchVect)) {
        throw std::runtime_error("SubstractMatch failed.");
    }

    // Maps atom indices to 2D points
    RDGeom::INT_POINT2D_MAP coordinate_map;

    RDKit::Conformer& conf = this->rdkit_molecule->getConformer();
    for(auto mv: matchVect) {
        RDGeom::Point3D pt3 = conf.getAtomPos( mv.first );
        RDGeom::Point2D pt2( pt3.x , pt3.y );
        coordinate_map[mv.second] = pt2;
    }
    // what is going on here?
    // That doesn't seem to change much
    // RDDepict::compute2DCoords( *this->rdkit_molecule, &coordinate_map, true, true);

    if(coordinate_map.empty()) {
        throw std::runtime_error("RDKit coordinate mapping is empty");
    }

    return coordinate_map;
}

void CanvasMolecule::process_alignment_in_rings() {
    const auto& rings = this->rdkit_molecule->getRingInfo();
    // g_debug("Number of rings: %lu", rings->atomRings().size());

    for(const auto& ring: rings->atomRings()) {

        float ring_center_x = 0.f;
        float ring_center_y = 0.f;
        for(int atom_idx: ring) {
            ring_center_x += this->atoms.at(atom_idx).x;
            ring_center_y += this->atoms.at(atom_idx).y;
        }
        ring_center_x /= ring.size();
        ring_center_y /= ring.size();
        
        int i = 0;
        int j = 1;

        // Go over every bond
        while(i!=ring.size()) {
            int atom_one_idx = ring[i];
            int atom_two_idx = ring[j];
            auto bonds_of_atom_one = this->bond_map.find(atom_one_idx);
            if(bonds_of_atom_one == this->bond_map.end()) {
                throw std::runtime_error("Critical internal error: Could not find a bond while processing rings.");
            }
            // Find iterator pointing to the bond
            auto bond = std::find_if(bonds_of_atom_one->second.begin(),bonds_of_atom_one->second.end(),[=](const auto bond){
                return (bond->first_atom_idx == atom_one_idx && bond->second_atom_idx == atom_two_idx) 
                    || (bond->first_atom_idx == atom_two_idx && bond->second_atom_idx == atom_one_idx);
            });
            if(bond == bonds_of_atom_one->second.end()) {
                throw std::runtime_error("Critical internal error: Could not find a bond while processing rings.");
            }
            auto bond_ptr = *bond;
            if(bond_ptr->type == BondType::Double) {
                float x_offset_from_center = (bond_ptr->first_atom_x + bond_ptr->second_atom_x) / 2.f - ring_center_x;
                // negative y on screen is actually "higher" so we need to flip the sign
                float y_offset_from_center = ring_center_y - (bond_ptr->first_atom_y + bond_ptr->second_atom_y) / 2.f;
                bool sign_of_x_offset_from_center = x_offset_from_center > 0.f;
                bool sign_of_y_offset_from_center = y_offset_from_center > 0.f;
                bool x_requirement = bond_ptr->second_atom_x > bond_ptr->first_atom_x == sign_of_y_offset_from_center;
                // negative y on screen is actually "higher" so we need to flip the sign
                bool y_requirement = bond_ptr->second_atom_y <= bond_ptr->first_atom_y != sign_of_x_offset_from_center;
                bool bond_direction = x_requirement && y_requirement;
                // g_debug(
                //     "Bond: %i->%i DeltaX: %f DeltaY: %f CX: %f CY: %f XO: %f SignXO: %i YO: %f SignYO: %i ReqX: %i ReqY: %i DIR: %i",
                //     bond->first_atom_idx,
                //     bond->second_atom_idx,
                //     bond->second_atom_x - bond->first_atom_x,
                //     bond->first_atom_y - bond->second_atom_y,
                //     ring_center_x,ring_center_y,
                //     x_offset_from_center,
                //     sign_of_x_offset_from_center,
                //     y_offset_from_center,
                //     sign_of_y_offset_from_center,
                //     x_requirement,
                //     y_requirement,bond_direction
                // );
                bond_ptr->bond_drawing_direction = bond_direction ? DoubleBondDrawingDirection::Primary : DoubleBondDrawingDirection::Secondary;
            }
            // Lastly, process appendices' alignment relative to ring center
            auto& atom = this->atoms.at(atom_one_idx);
            if(atom.appendix.has_value()) {
                bool should_reverse = atom.x < ring_center_x;
                atom.appendix->reversed = should_reverse;
            }
            i++;
            j++;
            // Process the last bond
            if(j==ring.size()) {
                // Loop j to point to the first atom in the ring
                j = 0;
            }
        }

    }
}

void CanvasMolecule::shorten_double_bonds() {
    typedef std::pair<const Bond*,float> bond_ptr_and_angle;
    for(auto& bond: this->bonds) {
        if(bond->type != BondType::Double) {
            continue;
        }
        if(bond->bond_drawing_direction.has_value()) {
            if(bond->bond_drawing_direction == DoubleBondDrawingDirection::Centered) {
                continue;
            }
        }
        auto find_angle_between_bonds = [&](const Bond* other_bond, bool flip){
            // We can find the angle between two bonds
            // by computing cosinus arcus (reverse cosine)
            // of ( dot product / (length a) * (length b) )
            auto [bond_vec_x,bond_vec_y] = bond->get_vector();
            auto [other_bond_vec_x,other_bond_vec_y] = other_bond->get_vector();
            if(flip) {
                other_bond_vec_x *= -1.f;
                other_bond_vec_y *= -1.f;
            }
            auto dot_product = (bond_vec_x*other_bond_vec_x) + (bond_vec_y*other_bond_vec_y);
            auto bond_length = bond->get_length();
            auto other_bond_length = other_bond->get_length();
            auto result = std::acos(dot_product/(bond_length*other_bond_length));
            return result;
        };
        // 1. Find the adjacent bond(s)
        auto find_adjacent_bonds = [this,&bond,&find_angle_between_bonds]() -> std::pair<std::vector<bond_ptr_and_angle>,std::vector<bond_ptr_and_angle>> {
            // Adjacent bonds touching the first atom
            std::vector<bond_ptr_and_angle> first_bonds;
            // Adjacent bonds touching the second atom
            std::vector<bond_ptr_and_angle> second_bonds;

            auto first_bonds_iter = this->bond_map.find(bond->first_atom_idx);
            if(first_bonds_iter != this->bond_map.end()) {
                // Going over bonds of the first atom in the currently evaluated bond
                for(const auto i: first_bonds_iter->second) {
                    if(i->first_atom_idx == bond->first_atom_idx) {
                        if(i->second_atom_idx == bond->second_atom_idx) {
                            // We're looking at the 'bond' itself. We must skip it.
                            continue;
                        }
                        first_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), false)));
                    } else if(i->second_atom_idx == bond->first_atom_idx) {
                        // i's second atom is bond's first, so we need to flip the sign of the bond vectors
                        // so that we can correctly compute the angle between them.
                        first_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), true)));
                    } else {
                        throw std::runtime_error("Internal error: bond_map is inconsistent!");
                    }
                }
            }
            auto second_bonds_iter = this->bond_map.find(bond->second_atom_idx);
            if(second_bonds_iter != this->bond_map.end()) {
                // Going over bonds of the second atom in the currently evaluated bond
                for(const auto i: second_bonds_iter->second) {
                    if(i->first_atom_idx == bond->second_atom_idx) {
                        // i's first atom is bond's second, so we need to flip the sign of the bond vectors
                        // so that we can correctly compute the angle between them.
                        second_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), true)));
                    } else if(i->second_atom_idx == bond->second_atom_idx) {
                        if(i->first_atom_idx == bond->first_atom_idx) {
                            // We're looking at the 'bond' itself. We must skip it.
                            continue;
                        }
                        second_bonds.push_back(std::make_pair(i.get(), find_angle_between_bonds(i.get(), false)));
                    } else {
                        throw std::runtime_error("Internal error: bond_map is inconsistent!");
                    }
                }
            }
            return std::make_pair(first_bonds,second_bonds);
        };
        auto compute_shortening_proportion = [&](const Bond* other_bond, float angle_between_bonds){
            // 3. Do a little trigonometry to find the length to be shortened
            auto absolute_shortened_length = BOND_LINE_SEPARATION / std::tan(angle_between_bonds/2.f);
            // 4. Find the proportion of the shortening
            auto bond_length = bond->get_length();
            return absolute_shortened_length / bond_length;
        };
        auto [first_bonds, second_bonds] = find_adjacent_bonds();

        auto element_with_smallest_angle_between_bonds = [&](const std::vector<bond_ptr_and_angle>& adjacent_bonds) -> bond_ptr_and_angle {
            auto min = std::min_element(adjacent_bonds.cbegin(),adjacent_bonds.cend(),[&](const auto& lhs, const auto& rhs){
                auto [lbond, langle] = lhs;
                auto [rbond, rangle] = rhs;
                return langle < rangle;
            });
            return std::make_pair(min->first,min->second);
        };
        if(!first_bonds.empty()) {
            auto [adjbond, angle] = element_with_smallest_angle_between_bonds(first_bonds);
            bond->first_shortening_proportion = compute_shortening_proportion(adjbond, angle);
        }
        if(!second_bonds.empty()) {
            auto [adjbond, angle] = element_with_smallest_angle_between_bonds(second_bonds);
            bond->second_shortening_proportion = compute_shortening_proportion(adjbond, angle);
        }
    }
}

void CanvasMolecule::build_internal_molecule_representation(const RDGeom::INT_POINT2D_MAP &coordinate_map) {
    // First, clear what we have

    this->atoms.clear();
    this->bonds.clear();
    this->bond_map.clear();
    this->bounding_atom_coords = std::make_pair(RDGeom::Point2D(0,0),RDGeom::Point2D(0,0));

    /// Used to avoid duplicating bonds
    std::set<unsigned int> processed_atoms_indices;
    
    // 1. Process atoms and compute bonds
    for(const auto& [atom_idx,plane_point]: coordinate_map) {
        const auto* rdkit_atom = this->rdkit_molecule->getAtomWithIdx(atom_idx);
        auto canvas_atom = CanvasMolecule::Atom();
        canvas_atom.color = atom_color_from_rdkit(rdkit_atom);
        if(rdkit_atom->hasProp("name")) {
            std::string atom_name;
            rdkit_atom->getProp("name", atom_name);
            canvas_atom.name = atom_name;
        }
        canvas_atom.highlighted = false;
        canvas_atom.idx = atom_idx;
        canvas_atom.symbol = rdkit_atom->getSymbol();
        canvas_atom.x = plane_point.x;
        canvas_atom.y = plane_point.y;

        if(canvas_atom.x < bounding_atom_coords.first.x) {
            bounding_atom_coords.first.x = canvas_atom.x;
        }
        if(canvas_atom.x > bounding_atom_coords.second.x) {
            bounding_atom_coords.second.x = canvas_atom.x;
        }
        if(canvas_atom.y < bounding_atom_coords.first.y) {
            bounding_atom_coords.first.y = canvas_atom.y;
        }
        if(canvas_atom.y > bounding_atom_coords.second.y) {
            bounding_atom_coords.second.y = canvas_atom.y;
        }

        auto surrounding_hydrogen_count = rdkit_atom->getTotalNumHs(false);
        auto surrounding_non_hydrogen_count = 0;
        auto charge = rdkit_atom->getFormalCharge();
        if (charge != 0) {
            Atom::Appendix ap;
            ap.charge = charge;
            canvas_atom.appendix = ap;
        }

        // Bond pointers to be stored in the `bond_map`
        std::vector<std::shared_ptr<Bond>> bonds_to_be_cached;
        // Used to determine if the 'appendix' should be 'reversed'
        std::optional<float> x_coordinate_of_bonded_atom;

        for(const auto& bond: boost::make_iterator_range(this->rdkit_molecule->getAtomBonds(rdkit_atom))) {
            // Based on `getAtomBonds` documentation.
            // Seems weird but we have to do it that way.
            const auto* bond_ptr = (*this->rdkit_molecule)[bond];
            auto first_atom_idx = bond_ptr->getBeginAtomIdx();
            auto second_atom_idx = bond_ptr->getEndAtomIdx();
            auto the_other_atom_idx = first_atom_idx == atom_idx ? second_atom_idx : first_atom_idx;
            const auto* the_other_atom =  this->rdkit_molecule->getAtomWithIdx(the_other_atom_idx);
            if(the_other_atom->getSymbol() != "H") {
                surrounding_non_hydrogen_count++;
                x_coordinate_of_bonded_atom = coordinate_map.at(the_other_atom_idx).x;
            } 
            // else {
            //     g_warning("Skipping explicit hydrogen bound to atom with idx=%u!",canvas_atom.idx);
            //     continue;
            // }

            // We don't want to have duplicate bonds of atoms that we have already processed
            // so we skip them.
            if(processed_atoms_indices.find(first_atom_idx) != processed_atoms_indices.end() 
            || processed_atoms_indices.find(second_atom_idx) != processed_atoms_indices.end()) {
                continue;
            }

            auto canvas_bond = CanvasMolecule::Bond();

            canvas_bond.first_atom_idx = first_atom_idx;
            canvas_bond.first_atom_x = coordinate_map.at(first_atom_idx).x;
            canvas_bond.first_atom_y = coordinate_map.at(first_atom_idx).y;
            
            canvas_bond.second_atom_idx = second_atom_idx;
            canvas_bond.second_atom_x = coordinate_map.at(second_atom_idx).x;
            canvas_bond.second_atom_y = coordinate_map.at(second_atom_idx).y;

            canvas_bond.highlighted = false;
            canvas_bond.type = bond_type_from_rdkit(bond_ptr->getBondType());
            canvas_bond.geometry = bond_geometry_from_rdkit(bond_ptr->getBondDir());

            auto canvas_bond_ptr = std::make_shared<Bond>(canvas_bond);
            this->bonds.push_back(canvas_bond_ptr);

            bonds_to_be_cached.push_back(canvas_bond_ptr);
            auto cached_bonds_for_other_atom = this->bond_map.find(the_other_atom_idx);
            if(cached_bonds_for_other_atom == this->bond_map.end()) {
                std::vector<std::shared_ptr<Bond>> vec;
                vec.push_back(canvas_bond_ptr);
                this->bond_map.emplace(std::pair(the_other_atom_idx,std::move(vec)));
            } else {
                cached_bonds_for_other_atom->second.push_back(canvas_bond_ptr);
            }
        }

        auto cached_bonds_for_this_atom = this->bond_map.find(atom_idx);
        if(cached_bonds_for_this_atom == this->bond_map.end()) {
            this->bond_map.emplace(std::pair(atom_idx,std::move(bonds_to_be_cached)));
        } else {
            // In the for-loop above,
            // We skip bonds with atoms that have already been processed.
            // Therefore just appending the whole vector 
            // should never result in a duplicate.
            std::move(bonds_to_be_cached.begin(),bonds_to_be_cached.end(),std::back_inserter(cached_bonds_for_this_atom->second));
        }

        bool terminus = surrounding_non_hydrogen_count < 2;
        if(canvas_atom.symbol != "H" && (canvas_atom.symbol != "C" || terminus)) {
            //todo: oxygens I guess?
            if(surrounding_hydrogen_count > 0) {
                Atom::Appendix ap = canvas_atom.appendix.value_or(Atom::Appendix());
                ap.superatoms = "H";
                if(surrounding_hydrogen_count > 1) {
                    ap.superatoms += std::to_string(surrounding_hydrogen_count);
                }
                if(terminus && x_coordinate_of_bonded_atom.has_value()) {
                    float diff = x_coordinate_of_bonded_atom.value() - canvas_atom.x;
                    ap.reversed = diff > 0.2;
                }
                canvas_atom.appendix = ap;
            }
        }
        auto setup_potential_centered_double_bond = [&](const int atom_idx){
            auto bonds_of_this_atom = this->bond_map.find(atom_idx);
            if(bonds_of_this_atom == this->bond_map.end()) {
                return;
            }
            if(!terminus || bonds_of_this_atom->second.empty()) {
                return;
            }
            // This means that we only have one bond
            Bond* bond = bonds_of_this_atom->second.front().get();
            if(bond->type != BondType::Double) {
                return;
            }
            unsigned int the_other_atom_idx = bond->first_atom_idx == atom_idx ? bond->second_atom_idx : bond->first_atom_idx;
            if(this->rdkit_molecule->getAtomWithIdx(the_other_atom_idx)->getAtomicNum() == 6) {
                // This should always be a valid iterator at this point
                auto bonds_of_the_other_atom = this->bond_map.find(the_other_atom_idx);
                if(bonds_of_the_other_atom->second.size() != 3 && bonds_of_the_other_atom->second.size() != 1) {
                    return;
                }
            }
            bond->bond_drawing_direction = DoubleBondDrawingDirection::Centered;
        };
        setup_potential_centered_double_bond(atom_idx);

        this->atoms.push_back(std::move(canvas_atom));
        
        // Mark the atom as processed
        processed_atoms_indices.insert(atom_idx);
    }
    std::sort(this->atoms.begin(),this->atoms.end(),[](const auto& lhs, const auto& rhs){
        return lhs.idx < rhs.idx;
    });
    // Make sure that double bonds are aligned properly
    this->process_alignment_in_rings();
    this->shorten_double_bonds();
}

void CanvasMolecule::lower_from_rdkit(bool sanitize_after) {

    // 2. Do the lowering

    // 2.0 Kekulize
    RDKit::MolOps::Kekulize(*this->rdkit_molecule);

    /// 2.1 Compute geometry
    auto geometry = this->compute_molecule_geometry();

    // 2.2 Build internal repr
    this->build_internal_molecule_representation(geometry);    
    this->cached_atom_coordinate_map = std::move(geometry);

    // 2.3 Reverse kekulization on the original molecule after lowering.
    if (sanitize_after) {
        RDKit::MolOps::sanitizeMol(*this->rdkit_molecule);
    }
}

void CanvasMolecule::highlight_atom(int atom_idx) {
    //g_debug("Highlighted atom with idx=%i",atom_idx);
    this->atoms.at(atom_idx).highlighted = true;
}

void CanvasMolecule::highlight_bond(int atom_a, int atom_b) {
    //g_debug("Highlighted bond between atoms with indices %i and %i",atom_a,atom_b);
    auto bonds_for_atom_a = this->bond_map.find(atom_a);
    if(bonds_for_atom_a == this->bond_map.end()) {
        throw std::runtime_error("Bond doesn't exist");
    }
    auto target = std::find_if(bonds_for_atom_a->second.begin(),bonds_for_atom_a->second.end(),[=](const auto& bond){
        return bond->second_atom_idx == atom_b && bond->first_atom_idx == atom_a;
    });
    if (target == this->bonds.end()) {
        throw std::runtime_error("Bond doesn't exist");
    }
    (*target)->highlighted = true;
}

void CanvasMolecule::clear_highlights() {
    for(auto& bond: this->bonds) {
        bond->highlighted = false;
    }
    for(auto& atom: this->atoms) {
        atom.highlighted = false;
    }
}

void CanvasMolecule::clear_cached_atom_coordinate_map() {
    this->cached_atom_coordinate_map = std::nullopt;
}

void CanvasMolecule::update_cached_atom_coordinate_map_after_atom_removal(unsigned int removed_atom_idx) {
    if (this->cached_atom_coordinate_map.has_value()) {
        auto& coordinate_map = this->cached_atom_coordinate_map.value();
        // There's no point in working with an empty cache.
        // Might as well delete it.
        if (coordinate_map.empty()) {
            this->cached_atom_coordinate_map = std::nullopt;
            return;
        }
        auto to_be_removed = coordinate_map.find(removed_atom_idx);
        if (to_be_removed == coordinate_map.end()) {
            g_warning("Atom to be removed (idx=%u) does not exist in the cached coordinate map!",removed_atom_idx);
            return;
        }
        // We can now remove the atom from the cache.
        coordinate_map.erase(to_be_removed);
        // If that was the only atom in the cache, let's get rid of the cache.
        if (coordinate_map.empty()) {
            this->cached_atom_coordinate_map = std::nullopt;
            return;
        }
        // After the atom has been removed, all of the indices which are
        // greater than its' index have to be decremented.
        //
        // Now, I don't want to mess around with iterator invalidation.
        // So the easiest course of action is to just copy away all 
        // of the affected elements, change them, remove them from the original map
        // and put them back there again.
        std::vector<std::pair<int,RDGeom::Point2D>> altered_elements;

        // This is the highest index AFTER the atom has been removed from the cache.
        // The only reason we get it is that we can make a smart allocation for our vector.
        auto highest_idx = coordinate_map.crbegin()->first;
        // Allocate all the space that we might need
        altered_elements.reserve(highest_idx - removed_atom_idx + 1);
        // Get the affected elements, copy them over to the vector, with decreased indices
        for(auto i = coordinate_map.upper_bound(removed_atom_idx);i != coordinate_map.end(); i++) {
            altered_elements.push_back(std::make_pair(i->first - 1,i->second));
        }
        // Now, remove the original elements from the map
        coordinate_map.erase(coordinate_map.upper_bound(removed_atom_idx),coordinate_map.end());
        // Then put them back again after they've been modified
        for(auto x: std::move(altered_elements)) {
            coordinate_map.emplace(x.first,x.second);
        }
    }
}