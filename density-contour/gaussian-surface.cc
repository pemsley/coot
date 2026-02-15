/*
 * density-contour/gaussian-surface.cc
 * 
 * Copyright 2023 by Medical Research Council
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#include <functional>
#include <chrono>

#include <clipper/core/nxmap.h>
#include <clipper/core/xmap.h>
#include <clipper/ccp4/ccp4_map_io.h>

#include "gaussian-surface.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "density-contour-triangles.hh"
#include "CIsoSurface.h"
#include "coot-utils/coot-map-utils.hh"

coot::simple_mesh_t
coot::gaussian_surface_t::get_surface() const {

   // std::cout << "debug::get_surface() n-vertices: " << mesh.vertices.size()
   //           << " n-triangles: " << mesh.triangles.size() << std::endl;
   return mesh;
}

void
coot::gaussian_surface_t::using_an_xmap_with_atom_selection(mmdb::Manager *mol, const std::string &cid,
                                                            float sigma_in, float contour_level_in,
                                                            float box_radius_in, float grid_scale,
                                                            float fft_b_factor) {

   std::cout << "using_an_xmap_with_atom_selection() grid_scale " << grid_scale << std::endl;
   int by_chain = 0; // i.e. use the cid for atom selection.
   using_an_xmap_internal(mol, cid, by_chain, sigma_in, contour_level_in, box_radius_in, grid_scale, fft_b_factor);
}

void
coot::gaussian_surface_t::using_an_xmap(mmdb::Manager *mol, const std::string &chain_id,
                                        float sigma_in, float contour_level_in, float box_radius_in, float grid_scale,
                                        float fft_b_factor) {

   std::cout << "using_an_xmap() grid_scale " << grid_scale << std::endl;
   int by_chain = 1;
   using_an_xmap_internal(mol, chain_id, by_chain, sigma_in, contour_level_in, box_radius_in, grid_scale, fft_b_factor);

}

//  [expo_index_to_float, mmdb_to_clipper, place_atom_in_grid]
clipper::Xmap<float>
coot::gaussian_surface_t::make_and_fill_map(mmdb::Manager *mol, int sel_hnd, float gs,
                                            std::pair<clipper::Coord_orth, clipper::Coord_orth> extents,
                                            float sigma, float box_radius) {

   std::cout << "------------------- make_and_fill_map() with gs " << gs << std::endl;

   auto expo_index_to_float = [] (int i) {
      return static_cast<float>(i) * 0.01;
   };

   auto mmdb_to_clipper = [] (mmdb::Atom *at) {
      return clipper::Coord_orth(at->x, at->y, at->z);
   };

   auto place_atom_in_grid = [] (const clipper::Coord_orth &pt,
                                 clipper::Xmap<float> &xmap,
                                 float sigma,
                                 float box_radius) {

      // std::cout << "place_atom_in_grid " << pt.format() << std::endl;
      clipper::Coord_frac centre_f = pt.coord_frac(xmap.cell());
      float box_radius_sqrd = box_radius * box_radius;

      clipper::Coord_frac box0(
                               centre_f.u() - box_radius/xmap.cell().descr().a(),
                               centre_f.v() - box_radius/xmap.cell().descr().b(),
                               centre_f.w() - box_radius/xmap.cell().descr().c() );
      clipper::Coord_frac box1(
                               centre_f.u() + box_radius/xmap.cell().descr().a(),
                               centre_f.v() + box_radius/xmap.cell().descr().b(),
                               centre_f.w() + box_radius/xmap.cell().descr().c() );

      clipper::Grid_map grid( box0.coord_grid(xmap.grid_sampling()),
                              box1.coord_grid(xmap.grid_sampling()));

      clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min() ), iu, iv, iw;
      for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() )  {
         for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
               clipper::Coord_grid c_g = iw.coord();
               clipper::Coord_frac c_f = c_g.coord_frac(xmap.grid_sampling());
               clipper::Coord_orth c_o = c_f.coord_orth(xmap.cell());
               // float plength = clipper::Coord_orth::length(c_o, pt);
               float plength_sqrd = (c_o - pt).lengthsq();
               float z_sqrd = plength_sqrd/(sigma*sigma);
               float x_prime = - z_sqrd;
               if (plength_sqrd < box_radius_sqrd) {
                  // int idx = float_to_expo_index(x_prime);
                  // float v = expo[idx];
                  float atomic_number_scaling = 1.0;
                  float v = atomic_number_scaling * expf(x_prime);
                  xmap[iw] += v;
                  if (false)
                     std::cout << "adding " << v << " to " << c_g.format() << " from "
                               << x_prime << " using z " << std::sqrt(z_sqrd) << " x_prime " << x_prime
                               << std::endl;
               }
            }
         }
      }
   };

   mmdb::PAtom *atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(sel_hnd, atom_selection, n_selected_atoms);

   int nx = (extents.second.x() - extents.first.x()) * gs;
   int ny = (extents.second.y() - extents.first.y()) * gs;
   int nz = (extents.second.z() - extents.first.z()) * gs;

   const clipper::Coord_orth &coords_base = extents.first;
   clipper::Coord_orth dimensions = extents.second - extents.first;

   clipper::Grid grid(nx, ny, nz);
   clipper::Grid_sampling grid_sampling(nx, ny, nz);
   float ninety = clipper::Util::d2rad(90.0);
   clipper::Cell_descr cell_desc(dimensions.x(), dimensions.y(), dimensions.z(), ninety, ninety, ninety);
   clipper::Cell cell(cell_desc);

   if (true) {
      std::cout << "debug:: cell and grid sampling: " << cell.format() << " "
                << grid_sampling.format() << " using gs " << gs << "\n";
   }

   clipper::Spacegroup spacegroup(clipper::Spacegroup::P1);
   clipper::Xmap<float> xmap(spacegroup, cell, grid_sampling);

   unsigned int n_expo_values = 10000;
   std::vector<float> expo(n_expo_values);
   for (unsigned int i=0; i<n_expo_values; i++) {
      float x = expo_index_to_float(i);
      expo[i] = exp(x);
   }

   for (int iat=0; iat<n_selected_atoms; iat++) {
      mmdb:: Atom *atom_p = atom_selection[iat];
      if (! atom_p->isTer()) {
         mmdb::Residue *residue_p = atom_p->GetResidue();
         if (residue_p) {
            std::string res_name_this = residue_p->GetResName();
            if (res_name_this == "HOH") continue;
            clipper::Coord_orth pt = mmdb_to_clipper(atom_p);
            clipper::Coord_orth position_in_grid = pt - coords_base;
            place_atom_in_grid(position_in_grid, xmap, sigma, box_radius); // change xmap
         }
      }
   }
   return xmap;
}


void
coot::gaussian_surface_t::using_an_xmap_internal(mmdb::Manager *mol, const std::string &chain_id, int chain_selection_mode,
                                        float sigma_in, float contour_level_in, float box_radius_in, float grid_scale,
                                        float fft_b_factor) {

   std::cout << "using_an_xmap_internal with grid_scale " << grid_scale << std::endl;

   // these affect the smoothness/resolution of the surface
   //
   // small sigmas mean higher resolution - 0.6 is very high resolution.

   // 20230206-PE use the passed parameters
   float sigma = sigma_in;
   float contour_level = contour_level_in;
   float box_radius = box_radius_in;
   float gs = grid_scale;
   float b_factor = fft_b_factor;

   auto clipper_to_cart = [] (const clipper::Coord_orth &co) {
      return coot::Cartesian(co.x(), co.y(), co.z());
   };

   int sel_hnd = mol->NewSelection(); // d

   if (chain_selection_mode) {
      mol->SelectAtoms (sel_hnd, 0, chain_id.c_str(),
                        mmdb::ANY_RES, // starting resno, an int
                        "*", // any insertion code
                        mmdb::ANY_RES, // ending resno
                        "*", // ending insertion code
                        "!HOH", // any residue name
                        "*",
                        "*", // elements
                        "*"  // alt loc.
                        );
   } else {
      // use the CID - which is the "chain_id" argument!
      mol->Select(sel_hnd, mmdb::STYPE_ATOM, chain_id.c_str(), mmdb::SKEY_OR);
   }

   mmdb::PAtom *atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(sel_hnd, atom_selection, n_selected_atoms);
   if (n_selected_atoms == 0)  return;

   std::pair<clipper::Coord_orth, clipper::Coord_orth> extents = util::extents(mol, sel_hnd);
   const clipper::Coord_orth &coords_base = extents.first;
   coot::Cartesian coords_base_cart = clipper_to_cart(coords_base);
   glm::vec3 coords_base_glm(coords_base.x(), coords_base.y(), coords_base.z());

   if (false)
      std::cout << "debug:: extents: "
                << extents.first.format() << " "
                << extents.second.format() << std::endl;

   // extend the extents
   clipper::Coord_orth delta(5,5,5);
   extents.first  -= delta;
   extents.second += delta;

   std::cout << "calling  make_and_fill_map with grid_scale " << grid_scale << std::endl;
   xmap = make_and_fill_map(mol, sel_hnd, grid_scale, extents, sigma, box_radius);

   if (b_factor > 0.0)
      util::sharpen_blur_map(&xmap, b_factor);

   Cartesian centre(10, 10, 10);
   // std::pair<bool, clipper::Coord_orth> cc = coot::centre_of_molecule(mol);
   // if (cc.first)
      // centre = coot::Cartesian(cc.second.x(), cc.second.y(), cc.second.z());
   clipper::Coord_orth e_midpoint = 0.5 * (extents.first + extents.second);
   centre = clipper_to_cart(e_midpoint);
   centre -= coords_base_cart;
   // std::cout << "--------- Centre of molecule " << centre << " ----------" << std::endl;
   int isample_step = 1;
   int iream_start = 0;
   int n_reams = 1;
   bool is_em_map = true;
   double diag_distance = clipper::Coord_orth::length(extents.first, extents.second);
   float dy_radius = diag_distance / 2.0;

   bool use_vertex_gradients_for_map_normals_flag = true;
   CIsoSurface<float> my_isosurface;
   coot::density_contour_triangles_container_t tri_con =
        my_isosurface.GenerateTriangles_from_Xmap(std::cref(xmap),
                                                  contour_level, dy_radius, centre, isample_step,
                                                  iream_start, n_reams, is_em_map,
                                                  use_vertex_gradients_for_map_normals_flag);

   if (false)
      std::cout << "tri_con points: " << tri_con.points.size()
                << " n-triangles " << tri_con.point_indices.size() << " "
                << std::endl;

   // transfer tri_con to simple_mesh_t
   mesh.vertices.reserve(tri_con.points.size());
   glm::vec4 col(0.5, 0.5, 0.5, 1.0);
   for (size_t i = 0; i < tri_con.points.size(); i++) {
      const auto &p = tri_con.points[i];
      const auto &n = tri_con.normals[i];
      glm::vec3 pos(p.x(), p.y(), p.z());
      glm::vec3 normal(n.x(), n.y(), n.z());
      coot::api::vnc_vertex vertex(pos + coords_base_glm, normal, col);
      mesh.vertices.push_back(vertex);
   }

   mesh.triangles.reserve(tri_con.point_indices.size());
   for (size_t i = 0; i < tri_con.point_indices.size(); i++) {
      const auto &t = tri_con.point_indices[i];
      g_triangle tri(t.pointID[0], t.pointID[1], t.pointID[2]);
      mesh.triangles.push_back(tri);
   }

   mol->DeleteSelection(sel_hnd);

   // 20240911-PE do I need this now?
   // I think not.
   // normals_from_function_gradient(xmap, coords_base_glm); // changes the normal in the verties of the mesh
}

clipper::Xmap<float> coot::gaussian_surface_t::get_xmap() const {

   return xmap;
}

void
coot::gaussian_surface_t::normals_from_function_gradient(const clipper::Xmap<float> &xmap,
                                                         const glm::vec3 &coords_base_glm) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   float delta = 0.92;
   for (unsigned int i=0; i<mesh.vertices.size(); i++) {
      const auto &pos = mesh.vertices[i].pos - coords_base_glm;
      clipper::Coord_orth p_x_1(pos.x - delta, pos.y, pos.z);
      clipper::Coord_orth p_x_2(pos.x + delta, pos.y, pos.z);
      clipper::Coord_orth p_y_1(pos.x, pos.y - delta, pos.z);
      clipper::Coord_orth p_y_2(pos.x, pos.y + delta, pos.z);
      clipper::Coord_orth p_z_1(pos.x, pos.y, pos.z - delta);
      clipper::Coord_orth p_z_2(pos.x, pos.y, pos.z + delta);
      float f_x_1 = util::density_at_point(xmap, p_x_1);
      float f_x_2 = util::density_at_point(xmap, p_x_2);
      float f_y_1 = util::density_at_point(xmap, p_y_1);
      float f_y_2 = util::density_at_point(xmap, p_y_2);
      float f_z_1 = util::density_at_point(xmap, p_z_1);
      float f_z_2 = util::density_at_point(xmap, p_z_2);
      float f = util::density_at_point(xmap, clipper::Coord_orth(pos.x, pos.y, pos.z));
      glm::vec3 grr(f_x_1 - f_x_2, f_y_1 - f_y_2, f_z_1 - f_z_2);
      if (false)
         std::cout << "pos " << pos.x << " " << pos.y << " " << pos.z << " "
                   << "grr " << grr.x << " " << grr.y << " " << grr.z << " "
                   << " f " << f << " "
                   << "fx  " << f_x_1 << " " << f_x_2 << " "
                   << "fy  " << f_y_1 << " " << f_y_2 << " "
                   << "fz  " << f_z_1 << " " << f_z_2 << " "
                   << std::endl;
      if (grr.x != 0.0) {
         if (grr.y != 0.0) {
            if (grr.z != 0.0) {
               glm::vec3 gr = glm::normalize(grr);
               mesh.vertices[i].normal = gr;
            }
         }
      }
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();

   if (false)
      std::cout << "normals_from_function_gradient(): time " << d10 << " ms " << std::endl;
}

void
coot::gaussian_surface_t::using_calc_density(mmdb::Manager *mol) {

   int atom_selection_handle = mol->NewSelection();
   mol->SelectAtoms (atom_selection_handle, 0, "*",
                     mmdb::ANY_RES, // starting resno, an int
                     "*", // any insertion code
                     mmdb::ANY_RES, // ending resno
                     "*", // ending insertion code
                     "*", // any residue name
                     "*",
                     "*", // elements
                     "*"  // alt loc.
                     );


   mmdb::PAtom *atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(atom_selection_handle, atom_selection, n_selected_atoms);
   std::cout << "INFO:: selected " << n_selected_atoms << " atoms" << std::endl;

   clipper::Cell cell(clipper::Cell_descr(100, 100, 100, 90, 90, 90));
   clipper::Spacegroup space_group(clipper::Spacegroup::P1);
   clipper::Grid_sampling grid_sampling(100,100,100);

   clipper::Xmap<float> xmap =
      coot::util::calc_atom_map(mol, atom_selection_handle, cell, space_group, grid_sampling);

   if (true) {
      clipper::CCP4MAPfile outmapfile;
      outmapfile.open_write("atom_calc.map");
      outmapfile.export_xmap(xmap);
      outmapfile.close_write();
   }

   // now how do I contour that?

   float contour_level = 0.5;
   coot::Cartesian centre(50, 50, 50);
   int isample_step = 1;
   int iream_start = 0;
   int n_reams = 100;
   bool is_em_map = true;
   float dy_radius = 50.0;

   bool use_vertex_gradients_for_map_normals_flag = true;
   CIsoSurface<float> my_isosurface;
   coot::density_contour_triangles_container_t tri_con =
        my_isosurface.GenerateTriangles_from_Xmap(std::cref(xmap),
                                                  contour_level, dy_radius, centre, isample_step,
                                                  iream_start, n_reams, is_em_map,
                                                  use_vertex_gradients_for_map_normals_flag);

   std::cout << "tri_con points: " << tri_con.points.size() << " vertices " << tri_con.point_indices.size() << " triangles"
             << std::endl;

   // transfer tri_con to simple_mesh_t
   glm::vec4 col(0.5, 0.5, 0.5, 1.0);
   for (size_t i = 0; i < tri_con.points.size(); i++) {
      const auto &p = tri_con.points[i];
      const auto &n = tri_con.normals[i];
      glm::vec3 pos(p.x(), p.y(), p.z());
      glm::vec3 normal(n.x(), n.y(), n.z());
      coot::api::vnc_vertex vertex(pos, normal, col);
   }


}

// optional args: float sigma=4.4, float contour_level=4.0, float box_radius=5.0, float grid_scale=0.7);
coot::gaussian_surface_t::gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id,
                                             float sigma, float contour_level, float box_radius, float grid_scale,
                                             float b_factor) {

   // using_calc_density(mol);
   using_an_xmap(mol, chain_id, sigma, contour_level, box_radius, grid_scale, b_factor);

}

// optional args: float sigma=4.4, float contour_level=4.0, float box_radius=5.0, float grid_scale=0.7);
coot::gaussian_surface_t::gaussian_surface_t(mmdb::Manager *mol, const std::string &cid,
                                             int chain_vs_cid_mode,
                                             float sigma, float contour_level, float box_radius, float grid_scale,
                                             float b_factor) {

   std::cout << "constructor with grid_scale " << grid_scale << std::endl;

   // using_calc_density(mol);
   using_an_xmap_with_atom_selection(mol, cid, sigma, contour_level, box_radius, grid_scale, b_factor);

}
