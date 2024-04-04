/*
 * src/draw-2.hh
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

#include <string>
#include <gtk/gtk.h>
#include <glm/glm.hpp>

GtkWidget *create_gtkglarea_widget();

// void my_glarea_add_signals_and_events(GtkWidget *glarea);

// void init_central_cube();
// void draw_central_cube(GtkGLArea *glarea);

// glm::vec4 new_unproject(float z);

// glm::vec4 new_unproject(float mouse_x, float mouse_y, float z);

// glm::mat4 get_molecule_mvp();

// void setup_key_bindings();

// map the function already in the key map with name description to the given key
void remap_key(const std::string &description, int);

// void reset_frame_buffers();

// moved into graphics_info_t
// glm::vec4 new_unproject(float z);


glm::vec3 get_camera_up_direction(const glm::mat4 &mouse_quat_mat);

// gboolean
// glarea_tick_func(GtkWidget *widget,
//                  GdkFrameClock *frame_clock,
//                  gpointer data);

