/* layla/ligand_editor_canvas/render.cpp
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

#include "render.hpp"
#include "model.hpp"

using namespace coot::ligand_editor_canvas::impl;
using Atom = coot::ligand_editor_canvas::CanvasMolecule::Atom;
using BondGeometry = coot::ligand_editor_canvas::CanvasMolecule::BondGeometry;
using BondType = coot::ligand_editor_canvas::CanvasMolecule::BondType;
using DoubleBondDrawingDirection = coot::ligand_editor_canvas::CanvasMolecule::DoubleBondDrawingDirection;

const float MoleculeRenderContext::CENTERED_DOUBLE_BOND_LINE_SEPARATION = 0.2f;
// 10 degrees
const float MoleculeRenderContext::GEOMETRY_BOND_SPREAD_ANGLE = M_PI/18.f;
const float MoleculeRenderContext::WAVY_BOND_ARC_LENGTH = 0.25f;
const float MoleculeRenderContext::GEOMETRY_BOND_DASH_SEPARATION = 0.2f;

#ifndef __EMSCRIPTEN__
Renderer::Renderer(cairo_t* cr, PangoLayout* pango_layout) {
    this->cr = cr;
    this->pango_layout = pango_layout;
}
#else // __EMSCRIPTEN__ defined
// Lhasa-specific includes/definitions
Renderer::Renderer() {

}
#endif

Renderer::~Renderer() {
    #ifndef __EMSCRIPTEN__
    g_object_unref(this->pango_layout);
    cairo_destroy(this->cr);
    #else // __EMSCRIPTEN__ defined
    // Lhasa-specific includes/definitions
    #endif
}

MoleculeRenderContext::MoleculeRenderContext(const CanvasMolecule& cm, Renderer& ren, DisplayMode mode) 
:canvas_molecule(cm), ren(ren), display_mode(mode) {
    scale_factor = canvas_molecule.get_scale();
    x_offset = scale_factor * canvas_molecule.x_canvas_translation;
    y_offset = scale_factor * canvas_molecule.y_canvas_translation;
}

MoleculeRenderContext::~MoleculeRenderContext() {

}

void MoleculeRenderContext::draw_atoms() {
    #ifndef __EMSCRIPTEN__
    cairo_t* cr = ren.cr;
    PangoLayout* pango_layout = ren.pango_layout;
    cairo_set_line_width(cr, 0.5);
    // Used to truncate bonds not to cover atoms
    std::map<unsigned int,graphene_rect_t> atom_idx_to_canvas_rect;
    for(const auto& atom: canvas_molecule.atoms) {
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
        auto render_atom = [&](const CanvasMolecule::Atom& atom, DisplayMode render_mode = DisplayMode::Standard) -> std::pair<unsigned int,graphene_rect_t> {
            // pre-process text
            auto [r,g,b] = CanvasMolecule::atom_color_to_rgb(atom.color);
            const std::string color_str = CanvasMolecule::atom_color_to_html(atom.color);
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
            process_atom_highlight(atom);
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
                        process_atom_highlight(atom);
                    }
                } else {
                    atom_idx_to_canvas_rect.emplace(render_atom(atom));
                }
                break;
            }
        }
    }
    #else
    #warning TODO: Abstract-away drawing atoms for Lhasa
    #endif
}

void MoleculeRenderContext::process_atom_highlight(const Atom& atom) {
    if(atom.highlighted) {
        #ifndef __EMSCRIPTEN__
        cairo_t* cr = ren.cr;
        PangoLayout* pango_layout = ren.pango_layout;
        //cairo_move_to(cr, atom.x * scale_factor + x_offset + ATOM_HITBOX_RADIUS, atom.y * scale_factor + y_offset);
        cairo_new_sub_path(cr);
        cairo_set_source_rgb(cr, 0.0, 1.0, 0.5);
        cairo_arc(cr, atom.x * scale_factor + x_offset, atom.y * scale_factor + y_offset, CanvasMolecule::ATOM_HITBOX_RADIUS, 0, M_PI * 2.0);
        cairo_stroke_preserve(cr);
        cairo_set_source_rgba(cr, 0.0, 1.0, 0.5, 0.5);
        cairo_fill(cr);
        #else
        #warning TODO: Abstract-away drawing atom highlights for Lhasa
        #endif
    }
}

void MoleculeRenderContext::draw_bonds() {
    #ifndef __EMSCRIPTEN__
    cairo_t* cr = ren.cr;
    PangoLayout* pango_layout = ren.pango_layout;
    for(const auto& bond: canvas_molecule.bonds) {
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
                pv_x *= CanvasMolecule::BOND_LINE_SEPARATION;
                pv_y *= CanvasMolecule::BOND_LINE_SEPARATION;

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
    #else
    #warning TODO: Abstract-away drawing bonds for Lhasa
    #endif
}