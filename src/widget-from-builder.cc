/*
 * src/widget-from-builder.cc
 *
 * Copyright 2021 by Medical Research Council
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
//
#include "widget-from-builder.hh"
#include "graphics-info.h"

GtkWidget *widget_from_builder(const std::string &w_name, GtkBuilder *builder) {
   GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(GTK_BUILDER(builder), w_name.c_str()));
   return  w;
}

GtkWidget *widget_from_builder(const std::string &w_name) {

   GtkWidget *w = graphics_info_t::get_widget_from_builder(w_name);
   return w;
}

GtkWidget *widget_from_preferences_builder(const std::string &w_name) {

   GtkWidget *w = graphics_info_t::get_widget_from_preferences_builder(w_name);
   return w;
}
