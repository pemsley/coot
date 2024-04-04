/*
 * src/LigandViewMesh.hh
 *
 * Copyright 2020 by Medical Research Council
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

#ifndef LIGAND_VIEW_MESH_HH
#define LIGAND_VIEW_MESH_HH

#include <vector>
#include "Shader.hh"

class LigandViewMesh {

   // contains vectors for both lines and triangles. Don't use indexing to draw.   

   enum { VAO_NOT_SET = 99999999 };
   GLuint vao_triangles;
   GLuint vao_text;
   GLuint triangles_buffer_id;
   std::vector<glm::vec2> triangles_vertices;
   // I need a container for text too.
   std::string name;
   void init();
   bool first_time;
   bool draw_this_mesh;
   bool this_mesh_is_closed;
   void setup_buffers();

public:
   LigandViewMesh() { init(); }
   explicit LigandViewMesh(const std::string &n) : name(n) { init(); }
   std::string get_name() const { return name; }
   void import(const std::vector<glm::vec2> &triangle_vertices);
   bool get_draw_status() const { return draw_this_mesh; }
   void draw(Shader *shader_p, float widget_height, float widget_width);
   void close();
   void clear();
};


#endif // LIGAND_VIEW_MESH_HH
