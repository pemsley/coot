/*
 * density-contour/density-contour-triangles.hh
 * 
 * Copyright 2007 by The University of York
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


#ifndef DENSITY_CONTOUR_TRIANGLES_HH
#define DENSITY_CONTOUR_TRIANGLES_HH

#include <map>

#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include "geometry/residue-and-atom-specs.hh"

struct POINT3DID {
   unsigned int newID;
   float x, y, z;
};

typedef std::map<unsigned int, POINT3DID> ID2POINT3DID;

class TRIANGLE {
public:
   unsigned int pointID[3];
   bool reject_this; // set in post-processing
   clipper::Coord_orth mid_point;
   double back_front_projection_distance;
   float occlusion_factor;
   clipper::Coord_orth normal_for_flat_shading;
   bool operator<(const TRIANGLE &t) const {
      return (back_front_projection_distance < t.back_front_projection_distance);
   }
   TRIANGLE() {
      reject_this = false;
      mid_point = clipper::Coord_orth(0,0,0);
      back_front_projection_distance = 999.9;
      normal_for_flat_shading = clipper::Coord_orth(0,0,1);
      occlusion_factor = 0.0;
   }
};

namespace coot {


   class density_contour_triangles_container_t {
   public:
      // every vertex has a point, a normal and an occlusion_factor
      std::vector<clipper::Coord_orth> points;
      std::vector<clipper::Coord_orth> normals;
      std::vector<float> occlusion_factor;
      std::vector<TRIANGLE> point_indices;
      void depth_sort(const clipper::Coord_orth &back_plane_point,
                      const clipper::Coord_orth &front_plane_point);
      void calculate_normals(); // average normals on shared points
      void calculate_normals_for_vertices(const clipper::Xmap<float> &xmap); // 20240328-PE function gradient at vertex
      void remove_small_triangles();
      bool empty() const { return (points.empty()); }
      void clear() {
         points.clear();
         normals.clear();
         point_indices.clear();
      }
   };

}

#endif // DENSITY_CONTOUR_TRIANGLES_HH
