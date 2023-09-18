/* layla/layla_embedded.cpp
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

#include "layla_embedded.hpp"
#include "ligand_builder_ui.hpp"
#include "ligand_builder_state.hpp"
#include "ligand_builder_generators.hpp"

using namespace coot::layla;

GtkApplicationWindow* coot::launch_layla(GtkApplication* app) {
    
    GtkBuilder* builder;

    if(! coot::ligand_editor::global_layla_gtk_builder) {
        builder = load_gtk_builder();
        coot::ligand_editor::global_layla_gtk_builder = builder;
    } else {
        builder = coot::ligand_editor::global_layla_gtk_builder;
    }

    auto *win = coot::layla::setup_main_window(app, builder);
    
    gtk_window_set_hide_on_close(GTK_WINDOW(win), TRUE);

    coot::ligand_editor::global_generator_request_task_cancellable = nullptr;

    g_signal_connect(win, "hide", G_CALLBACK(+[](GtkWidget* self, gpointer user_data){
        g_info("Deinitializing global ligand editor state...");
        delete coot::ligand_editor::global_instance;
    }), nullptr);
    gtk_window_present(GTK_WINDOW(win));
    gtk_application_add_window(app, GTK_WINDOW(win));

    return win;
}

GtkApplicationWindow* coot::launch_layla(GtkApplication* app, std::unique_ptr<RDKit::RWMol>&& mol) {
    auto* win = launch_layla(app);
    // todo: setup
    return win;
}
