/*
 * src/widget-headers.hh
 *
 * Copyright 2019 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <vector>
#include <string>
#include <utility>
#include <gtk/gtk.h>

std::string pre_directory_file_selection(GtkWidget *sort_button);
void filelist_into_fileselection_clist(GtkWidget *fileselection, const std::vector<std::string> &v);

GtkWidget *wrapped_nothing_bad_dialog(const std::string &label);

std::pair<short int, float> float_from_entry(GtkWidget *entry);
std::pair<short int, int>   int_from_entry(GtkWidget *entry);

GtkWidget *
add_validation_mol_menu_item(int imol, const std::string &name, GtkWidget *menu, GCallback callback);
void create_initial_validation_graph_submenu_generic(GtkWidget *window1,
                                                     const std::string &menu_name,
                                                     const std::string &sub_menu_name);



// To be used to (typically) get the menu item text label from chain
// option menus (rather than the ugly/broken casting of
// GtkPositionType data.  A wrapper to a static graphics_info_t
// function.
std::string menu_item_label(GtkWidget *menu_item);

void
add_map_colour_mol_menu_item(int imol, const std::string &name,
                                  GtkWidget *sub_menu, GCallback callback);

void add_map_scroll_wheel_mol_menu_item(int imol, 
                                        const std::string &name,
                                        GtkWidget *menu, 
                                        GCallback callback);

