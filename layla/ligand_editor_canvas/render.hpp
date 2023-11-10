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

    std::map<unsigned int,graphene_rect_t> atom_idx_to_canvas_rect;


    void process_atom_highlight(const CanvasMolecule::Atom& atom);

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