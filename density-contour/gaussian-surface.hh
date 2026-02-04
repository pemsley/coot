/*
 * density-contour/gaussian-surface.hh
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

#ifndef GAUSSIAN_SURFACE_HH
#define GAUSSIAN_SURFACE_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>

#include "coot-utils/simple-mesh.hh"

namespace coot {

   class gaussian_surface_t {
      clipper::Xmap<float> xmap;
      simple_mesh_t mesh;
      void using_an_nxmap(mmdb::Manager *mol);
      void using_an_xmap(mmdb::Manager *mol, const std::string &chain_id,
                         float sigma, float contour_level, float box_radius, float grid_scale,
                         float b_factor);
      void using_an_xmap_with_atom_selection(mmdb::Manager *mol, const std::string &cid,
                         float sigma, float contour_level, float box_radius, float grid_scale,
                         float b_factor);
      void using_an_xmap_internal(mmdb::Manager *mol, const std::string &cid,
                                  int chain_mode,
                                  float sigma, float contour_level, float box_radius,
                                  float grid_scale, float b_factor);
      void using_calc_density(mmdb::Manager *mol);
      void normals_from_function_gradient(const clipper::Xmap<float> &xmap,
                                          const glm::vec3 &cb); // changes mesh normals
      clipper::Xmap<float> make_and_fill_map(mmdb::Manager *mol, int sel_hnd, float gs,
                                             std::pair<clipper::Coord_orth, clipper::Coord_orth> extents,
                                             float sigma, float box_radius);
   public:
      // explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id);
      explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id,
                                  float sigma=4.4, float contour_level=4.0, float box_radius=5.0,
                                  float grid_scale=0.7, float fft_b_factor=100.0);
      // dummy_mode is used here to change the args of the constuctor
      explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &atom_selection_cid,
                                  int dummy_mode,
                                  float sigma=4.4, float contour_level=4.0, float box_radius=5.0,
                                  float grid_scale=0.7, float fft_b_factor=100.0);
      simple_mesh_t get_surface() const;
      clipper::Xmap<float> get_xmap() const;
   };

}

#endif // GAUSSIAN_SURFACE_HH
