/*
 * src/draw.hh
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

#ifndef DRAW_HH
#define DRAW_HH

void setup_monkey_head();

void draw_monkey_head();

void draw_one_triangle();

void draw_molecular_triangles(GtkWidget *widget);

void setup_for_single_triangle();

void draw_single_triangle();

void gtk3_draw_molecules();
void test_gtk3_adjustment_changed(GtkAdjustment *adj, GtkWidget *window);

struct shader_program_source {
   std::string VertexSource;
   std::string FragmentSource;
};

shader_program_source parse_shader(const std::string &file_name);
unsigned int CreateShader(const std::string &vertex_shader, const std::string &fragment_shader);

// static int programID_global;
// static GLuint VertexArrayID;
// static GLuint molecules_imol_VertexArrayID[10];
// static int location_global;



#endif // DRAW_HH
