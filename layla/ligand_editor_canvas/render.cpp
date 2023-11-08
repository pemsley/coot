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
using namespace coot::ligand_editor_canvas::impl;

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