/*
 * api/coot-molecule-blender.cc
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include "coot-molecule.hh"

void
coot::molecule_t::make_mesh_for_bonds_for_blender(const std::string &mode, protein_geometry *geom, bool against_a_dark_background,
                                      float bond_width, float atom_radius_to_bond_width_ratio,
                                      int smoothness_factor) {

   bool show_atoms_as_aniso_flag = true; // pass this
   bool show_aniso_atoms_as_ortep = true; // pass this also
   bool show_aniso_atoms_as_empty = false; // pass this also
   // 20241202-PE start the stopwatch... when will Lucrezia complain?

   // from the swig version:
   // PyList_SetItem(r_py, 0, vertices_py);
   // PyList_SetItem(r_py, 1, tris_py);
   // PyList_SetItem(r_py, 2, face_colours_dict_py);
   // return r_py;
   // where a tri_py is PyList_New(4);
   // with PyList_SetItem(tri_py, 3, PyLong_FromLong(colour_index));

   bool draw_hydrogen_atoms_flag = true;
   bool draw_missing_loops_flag = true;
   float aniso_probability = 0.5f;

   instanced_mesh_t im = get_bonds_mesh_instanced(mode, geom, against_a_dark_background,
                                                  bond_width, atom_radius_to_bond_width_ratio,
                                                  show_atoms_as_aniso_flag,
                                                  aniso_probability,
                                                  show_aniso_atoms_as_ortep,
                                                  show_aniso_atoms_as_empty,
                                                  smoothness_factor, draw_hydrogen_atoms_flag, draw_missing_loops_flag);

   blender_mesh_t bm(im);
   blender_mesh = std::move(bm);

}

void
coot::molecule_t::make_mesh_for_molecular_representation_for_blender(const std::string &cid,
                                                                     const std::string &colour_scheme,
                                                                     const std::string &style,
                                                                     int secondary_structure_usage_flag) {

   simple_mesh_t mesh = get_molecular_representation_mesh(cid, colour_scheme, style, secondary_structure_usage_flag);

   blender_mesh_t bm(mesh);
   blender_mesh = std::move(bm);

}

void
coot::molecule_t::make_mesh_for_map_contours_for_blender(Cartesian position, float contour_level, float radius) {

   bool use_thread_pool = false; // pass this
   ctpl::thread_pool *thread_pool_p = nullptr; // pass this
   clipper::Coord_orth pos_co(position.x(), position.y(), position.z());
   simple_mesh_t sm = get_map_contours_mesh(pos_co, radius, contour_level, use_thread_pool, thread_pool_p);
   blender_mesh_t bm(sm);
   blender_mesh = std::move(bm);
}

void
coot::molecule_t::make_mesh_for_goodsell_style_for_blender(coot::protein_geometry *geom_p,
                                                           float colour_wheel_rotation_step,
                                                           float saturation,
                                                           float goodselliness) {

   simple_mesh_t sm = get_goodsell_style_mesh(geom_p, colour_wheel_rotation_step, saturation, goodselliness);
   blender_mesh_t bm(sm);
   blender_mesh = std::move(bm);
}



void
coot::molecule_t::make_mesh_for_gaussian_surface_for_blender(float sigma, float contour_level, float box_radius, float grid_scale,float b_factor) {

   simple_mesh_t sm = get_gaussian_surface(sigma, contour_level, box_radius, grid_scale, b_factor);
   blender_mesh_t bm(sm);
   blender_mesh = std::move(bm);
}


std::vector<float>
coot::molecule_t::get_vertices_for_blender() const {

   float sf = 1.0;
   std::vector<float> v(blender_mesh.vertices.size() * 3, 0);
   for (unsigned int i=0; i<blender_mesh.vertices.size(); i++) {
      v[i*3  ] = sf * blender_mesh.vertices[i].x;
      v[i*3+1] = sf * blender_mesh.vertices[i].y;
      v[i*3+2] = sf * blender_mesh.vertices[i].z;
   }
   return v;
}


std::vector<int>
coot::molecule_t::get_triangles_for_blender() const {

   std::vector<int> v(blender_mesh.triangles.size() * 4);
   for (unsigned int i=0; i<blender_mesh.triangles.size(); i++) {
      v[i*4  ] = blender_mesh.triangles[i].triangle.point_id[0];
      v[i*4+1] = blender_mesh.triangles[i].triangle.point_id[1];
      v[i*4+2] = blender_mesh.triangles[i].triangle.point_id[2];
      // std::cout << "colour index for tri-i " << blender_mesh.triangles[i].colour_index << std::endl;
      v[i*4+3] = blender_mesh.triangles[i].colour_index;
   }
   return v;
}

std::vector<float>
coot::molecule_t::get_colour_table_for_blender() const {

   int max_colour_index = -1;
   std::vector<float> v;
   std::map<int, glm::vec4>::const_iterator it;
   for (it=blender_mesh.colour_table.begin(); it!=blender_mesh.colour_table.end(); ++it) {
      if (it->first > max_colour_index)
         max_colour_index = it->first;
   }

   if (max_colour_index >= 0) {
      // these colours don't have alpha (at the moment) so that it is easier to import into Blender (I think)
      std::vector<float> vv((max_colour_index+1) * 4, -1);
      v = vv;
      for (int i=0; i<=max_colour_index; i++) {
         it = blender_mesh.colour_table.find(i);
         if (it == blender_mesh.colour_table.end()) {
            // missing colour
            v[i*4  ] = 0.3;
            v[i*4+1] = 0.0;
            v[i*4+2] = 0.3;
            v[i*4+3] = 1.0;
         } else {
            v[i*4  ] = it->second.r;
            v[i*4+1] = it->second.g;
            v[i*4+2] = it->second.b;
            v[i*4+3] = it->second.a;
         }
      }
   }
   return v;
}

