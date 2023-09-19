/* layla/ui.hpp
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

#ifndef LAYLA_UI_HPP
#define LAYLA_UI_HPP
#include <gtk/gtk.h>
#include "ligand_editor_canvas.hpp"
#include "state.hpp"

namespace coot::layla {

// This doesn't need to be public
// void setup_actions(LaylaState* state, GtkApplicationWindow* win, GtkBuilder* builder);

/// Setups the main window using GtkBuilder (created from 'layla.ui').
///
/// Instantiates CootLigandEditor canvas and connects its' signals.
/// Setups window actions.
/// Instantiates the global instance of LaylaState
GtkApplicationWindow* setup_main_window(GtkApplication* app, GtkBuilder* builder);

/// Fetches and loads 'layla.ui'
GtkBuilder* load_gtk_builder();

}

#endif //  LAYLA_UI_HPP