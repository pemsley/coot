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
#include "ui.hpp"
#include "state.hpp"
#include "generators.hpp"

using namespace coot::layla;

bool coot::is_layla_initialized() {
    return coot::layla::global_layla_gtk_builder != nullptr;
}

GtkApplicationWindow* coot::initialize_layla(GtkApplication* app) {
    if(coot::layla::global_layla_gtk_builder) {
        g_warning("Layla has already been initialized!");
        auto* ptr = gtk_builder_get_object(coot::layla::global_layla_gtk_builder, "layla_window");
        return GTK_APPLICATION_WINDOW(ptr);
    }

    GtkBuilder* builder = load_gtk_builder();
    coot::layla::global_layla_gtk_builder = builder;
    auto *win = coot::layla::setup_main_window(app, builder);

    gtk_window_set_hide_on_close(GTK_WINDOW(win), TRUE);

    coot::layla::global_generator_request_task_cancellable = nullptr;

    g_signal_connect(win, "hide", G_CALLBACK(+[](GtkWidget* self, gpointer user_data){
        g_info("Resetting global ligand editor state...");
        coot::layla::global_instance->reset();
    }), nullptr);

    gtk_application_add_window(app, GTK_WINDOW(win));

    return win;
}

void coot::deinitialize_layla() {
    if(!is_layla_initialized()) {
        g_error("coot::deinitialize_layla() called before coot::initialize_layla()");
    }
    auto* win = gtk_builder_get_object(coot::layla::global_layla_gtk_builder, "layla_window");
    gtk_window_destroy(GTK_WINDOW(win));
    delete coot::layla::global_instance;
    coot::layla::global_instance = nullptr;
    g_object_unref(coot::layla::global_layla_gtk_builder);
    coot::layla::global_layla_gtk_builder = nullptr;
    coot::layla::global_generator_request_task_cancellable = nullptr;
}

void coot::launch_layla() {
    if(!is_layla_initialized()) {
        g_error("coot::launch_layla() called before coot::initialize_layla()");
    }
    auto *win = gtk_builder_get_object(coot::layla::global_layla_gtk_builder, "layla_window");
    if(gtk_widget_get_visible(GTK_WIDGET(win))) {
        g_warning("Layla window is already visible!");
        return;
    }
    gtk_window_present(GTK_WINDOW(win));
}

void coot::launch_layla(std::shared_ptr<RDKit::RWMol> mol) {
    launch_layla();
    CootLigandEditorCanvas* canvas = coot::layla::global_instance->get_canvas();
#if GLIB_MAJOR_VERSION == 2 && GLIB_MINOR_VERSION >= 74 || GLIB_MAJOR_VERSION > 2
    struct _cb_data_t {
        CootLigandEditorCanvas* canvas;
        std::shared_ptr<RDKit::RWMol> mol;
    };
    auto* cbd = new _cb_data_t;
    cbd->canvas = canvas;
    cbd->mol = std::move(mol);
    g_idle_add_once([](gpointer user_data){
        _cb_data_t* cbd = (_cb_data_t*) user_data;
        coot_ligand_editor_canvas_append_molecule(cbd->canvas, std::move(cbd->mol));
        delete cbd;
    }, cbd);
#else
    std::cout << "WARNING:: Rebuild Coot against Glib >= 2.74. Functionality is broken." << std::endl;
#endif
}
