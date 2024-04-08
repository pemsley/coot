/*
 * ligand/map-point-cluster.hh
 *
 * Copyright 2017 by Medical Research Council
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

#ifndef MAP_POINT_CLUSTER_HH
#define MAP_POINT_CLUSTER_HH

#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>

namespace coot {

   class map_point_cluster {
   public:
      map_point_cluster() : std_dev(clipper::Coord_orth(0,0,0)),
                            eigenvectors_and_centre(clipper::Mat33<double>::identity()) {
         score = 0.0;
      };
      std::vector<clipper::Coord_grid> map_grid;
      float score;
      // clipper::Coord_orth centre;
      // clipper::Mat33<double> eigenvectors;
      clipper::Coord_orth std_dev;
      clipper::RTop_orth eigenvectors_and_centre;
      std::vector<double> eigenvalues;
      bool operator==(const map_point_cluster &mpc) const {
	 return (mpc.map_grid == map_grid);
      }
      double volume(const clipper::Xmap<float> &xmap_ref) const;
   }; 

   bool compare_clusters(const map_point_cluster &a,
			 const map_point_cluster &b);

}

#endif // MAP_POINT_CLUSTER_HH
