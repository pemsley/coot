/* layla/ligand_editor_canvas/render.hpp
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
#ifndef COOT_LIGAND_EDITOR_CANVAS_RENDER_HPP
#define COOT_LIGAND_EDITOR_CANVAS_RENDER_HPP

#ifndef __EMSCRIPTEN__
    #include <gtk/gtk.h>
    #include <pango/pango-layout.h>
#else // Lhasa-specific includes
    #include <graphene.h>
#endif
#include <map>
#include <tuple>
#include "model.hpp"

namespace coot::ligand_editor_canvas::impl {

struct Renderer {
    #ifndef __EMSCRIPTEN__
    cairo_t* cr;
    PangoLayout* pango_layout;
    /// Takes ownership of the pointers
    Renderer(cairo_t*,PangoLayout*);

    #else // __EMSCRIPTEN__ defined
    // Lhasa-specific includes/definitions
    Renderer();
    #endif

    // todo: Common rendering API

    ~Renderer();
};

/// Encapsulates routines and state to render CanvasMolecule
class MoleculeRenderContext {

    const CanvasMolecule& canvas_molecule;
    Renderer& ren;
    DisplayMode display_mode;
    float scale_factor;
    float x_offset;
    float y_offset;
    // Used to truncate bonds not to cover atoms
    std::map<unsigned int,graphene_rect_t> atom_idx_to_canvas_rect;


    void process_atom_highlight(const CanvasMolecule::Atom& atom);
    /// Returns the markup string + info if the appendix is reversed
    std::tuple<std::string, bool> process_appendix(const std::string& symbol, const std::optional<CanvasMolecule::Atom::Appendix>& appendix);
    /// Returns a pair of atom index and bonding rect
    std::pair<unsigned int,graphene_rect_t> render_atom(const CanvasMolecule::Atom& atom, DisplayMode render_mode = DisplayMode::Standard);
    // Returns on-screen bond coordinates
    // cropped not to overlap with atoms' symbols and appendices.
    // Accepts on-screen coordinates of bond atoms.
    std::pair<graphene_point_t,graphene_point_t> cropped_bond_coords(const graphene_point_t& first_atom, unsigned int first_atom_idx, const graphene_point_t& second_atom, unsigned int second_atom_idx);

    void draw_central_bond_line(const CanvasMolecule::Bond& bond);
    void draw_straight_wedge(const CanvasMolecule::Bond& bond, bool reversed);
    void draw_straight_dashed_bond(const CanvasMolecule::Bond& bond, bool reversed);
    void draw_wavy_bond(const CanvasMolecule::Bond& bond);
    void draw_side_bond_line(const CanvasMolecule::Bond& bond, bool addOrSub, std::optional<float> first_shortening_proportion, std::optional<float> second_shortening_proportion);
    void draw_centered_double_bond(const CanvasMolecule::Bond& bond);

    static const float CENTERED_DOUBLE_BOND_LINE_SEPARATION;
    static const float GEOMETRY_BOND_SPREAD_ANGLE;
    static const float WAVY_BOND_ARC_LENGTH;
    static const float GEOMETRY_BOND_DASH_SEPARATION;

    public:
    MoleculeRenderContext(const CanvasMolecule& cm, Renderer& ren, DisplayMode mode);
    ~MoleculeRenderContext();

    void draw_atoms();
    void draw_bonds();
};

} //coot::ligand_editor_canvas::impl

#endif // COOT_LIGAND_EDITOR_CANVAS_RENDER_HPP