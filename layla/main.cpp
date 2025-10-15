/* layla/main.cpp
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

// #include "Python.h"
#include <gtk/gtk.h>
#include <thread>
#include "utils/coot-utils.hh"
#include "state.hpp"
#include "generators.hpp"
#include "ui.hpp"
// #include "python_utils.hpp"


int main(int argc, char** argv) {

    using namespace coot::layla;

    // std::thread python_init_thread([argc,argv](){
    //     setup_python_basic(argc, argv);
    //     setup_python_module("coot");
    //     setup_python_module("coot_utils");
    // });

    // python_init_thread.detach();

    gtk_init();

    GtkApplication* app = gtk_application_new("org.pemsley.Layla", G_APPLICATION_NON_UNIQUE);

    GError *error = NULL;
    g_application_register(G_APPLICATION(app), NULL, &error);

    g_signal_connect(app, "activate", G_CALLBACK(+[](GtkApplication* app, gpointer user_data) {

        auto* builder = load_gtk_builder();
        coot::layla::global_layla_gtk_builder = builder;
        
        std::string package_dir = coot::package_data_dir();
        auto *win = coot::layla::setup_main_window(app, builder);

        auto* icon_theme = gtk_icon_theme_get_for_display(gtk_widget_get_display(GTK_WIDGET(win)));
        gboolean dark_mode_flag = FALSE;
        g_object_get(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", &dark_mode_flag, NULL);

        std::string prefix_dir = coot::prefix_dir();
        std::string pixmaps_path_for_icons = coot::util::append_dir_dir(package_dir, "pixmaps");
        if(dark_mode_flag) {
            g_info("Enabling dark theme icons.");
            pixmaps_path_for_icons = coot::util::append_dir_dir(pixmaps_path_for_icons, "dark");
        } else {
            g_info("Not using dark theme icons.");
        }
        std::string prefix_main_path_for_icons = coot::util::append_dir_dir(prefix_dir, "share");
        prefix_main_path_for_icons = coot::util::append_dir_dir(prefix_main_path_for_icons, "icons");
        gtk_icon_theme_add_search_path(icon_theme, pixmaps_path_for_icons.c_str());
        g_info("Appended icon theme search path: %s", pixmaps_path_for_icons.c_str());
        gtk_icon_theme_add_search_path(icon_theme, prefix_main_path_for_icons.c_str());

        coot::layla::global_generator_request_task_cancellable = nullptr;

        gtk_window_present(GTK_WINDOW(win));
        gtk_application_add_window(app, GTK_WINDOW(win));
    }), NULL);

    auto ret = g_application_run(G_APPLICATION(app), 0, 0);
    g_info("Exiting...");
    delete coot::layla::global_instance;
    return ret;
}
