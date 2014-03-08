/* build/CalphaBuild.hh
 * 
 * Copyright 2001, 2002, 2003, 2006 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#ifndef HAVE_CALPHABUILD_HH
#define HAVE_CALPHABUILD_HH
#include "clipper/core/coords.h"
#include "clipper/core/xmap.h"
#include "clipper/contrib/skeleton.h"


#ifndef HAVE_VECTOR
#include <vector>
#define HAVE_VECTOR
#endif
#include "angles/AngleInfo.h"

namespace coot {

   // typedef std::pair<float, clipper::Coord_orth> scored_skel_coord;

   class scored_skel_coord {

   public:
      scored_skel_coord(float s, const clipper::Coord_orth &p, const clipper::Coord_grid &g) {
	 score = s;
	 position = p;
	 near_grid_pos = g;
      } 
      clipper::Coord_orth position;
      float score;
      float skel_score;
      clipper::Coord_grid near_grid_pos;
      int skeletonization_level;
      int depth;
   };

   // These are points (clipper::Coord_grids) in a skeleton map and are
   // used to see if target point is "close to" the current point (by
   // being less than some (gridding dependent) critical value number of
   // nodes away).
   // 
   class SkeletonTreeNode { 
      
   public: 
      
      // I get stuck with de-referencing the ix of a neighbour (I mean I
      // don't know how to set it - we have a Coord_grid, not a
      // Map_reference_index).
      // 
      // clipper::Xmap_base::Map_reference_index ix;
      
      // clipper::Coord_grid cg; 
      
      // please be fast enough...
      // 
      // Do not include yourself.
      std::vector<clipper::Coord_grid> neighbs; 

      // itself
      clipper::Coord_grid near_grid_point; 
      
   };


   class CalphaBuild {

      AngleInfo ai;
      float skel_pts_mean_density;
      float skel_pts_mean_std_dev_density;
      std::vector<coot::scored_skel_coord> cluster(const std::vector<coot::scored_skel_coord> &interesting,
						   const clipper::Coord_orth &start) const;
      
      float score_position_by_angles(const std::vector<clipper::Coord_orth> &Previous_ca_positions,
				     const clipper::Coord_orth &this_position) const; 

      std::vector <coot::scored_skel_coord>
      next_ca_internal(const std::vector<clipper::Coord_orth> &Previous_ca_positions,
		       const clipper::Coord_grid &coord_grid_start,
		       short int use_coord_grid_start_flag,
		       float ca_bond_length,
		       const clipper::Xmap<int> &skel,
		       const clipper::Xmap<float> &map,
		       float map_cut_off,
		       const clipper::Xmap<SkeletonTreeNode> &treenodemap,
		       int idepth) const ; // when depth == 0, we are at the end,
                                          // i.e. depth > 0 (currently 1) at the start
                                          // of the search.

      // We get cluster centre, which are the average of skeleton
      // points.  Now these cluster centres are not at 3.8A (they are
      // often closer than that).  We want the point at the end of 3.8
      // (ca_bond_length) of a vector that starts at a and goes
      // through (strictly or points at) b.
      // 
      clipper::Coord_orth extend_to_ca_bond_length_position(const clipper::Coord_orth &a,
							    const clipper::Coord_orth &b,
							    float ca_bond_length) const;


      // return 0.1 for unreachable (in (20?) connections).
      // return 100.0 for connected.
      // 
      float score_by_skel(const clipper::Coord_grid &start,
			  const clipper::Coord_grid &this_grid,
			  short int use_coord_grid_start_flag,
			  const clipper::Xmap<int> &skel,
			  const clipper::Xmap<float> &map,
			  float map_cut_off,
			  const clipper::Xmap<SkeletonTreeNode> &treenodemap) const;

      float score_by_skel_internal(const clipper::Coord_grid &start,
				   const clipper::Coord_grid &target,
				   const clipper::Coord_grid &previous,
				   const clipper::Coord_grid &prev_prev,
				   const clipper::Xmap<int>  &skel_map,
				   const clipper::Xmap<float> &map,
				   float map_cut_off,
				   const clipper::Xmap<SkeletonTreeNode> &treenodemap,
				   int depth) const;

      // changes treenodemap
      // void build_tree_node_map(const clipper::Xmap<int> &skel_map); 

      int start_skel_depth;  // defaults to 10 but can be set by user for
			     // use in a high-res map.
							    
   public:

      CalphaBuild();
      CalphaBuild(int max_skeleton_search_depth);

      // Return a sorted by goodness set of potential ca_positions
      //
      // We now include a test to see if the potiential point is
      // reachable through the skeleton.  In order to do so, we need
      // to have a starting coord_grid that corresponds to the
      // starting coordinate.  That is given to us by the graphics and
      // has been tied to the coord_orth that it selects.
      //
      // Now, we must know of this coord_grid is to be used to find
      // skeleton connectivity - in the case of the first atom to be
      // placed (i.e. we have just generated the baton "here" and that
      // has no associated skeleton coordinate, it definately is not.
      //
      // TMP: we currently modify treenodemap by calling
      // score_by_skel.  This may change in the future, when
      // treenodemap is calculated externally.
      // 
      std::vector <scored_skel_coord>
      next_ca_by_skel(const std::vector<clipper::Coord_orth> &Previous_ca_positions,
		      const clipper::Coord_grid &coord_grid_start,
		      short int use_coord_grid_start_flag,
		      float ca_bond_length,
		      const clipper::Xmap<int> &skel,
		      const clipper::Xmap<float> &map,
		      float map_cut_off,
		      const clipper::Xmap<SkeletonTreeNode> &treenodemap) const;

      // find place in the skeleton that is close to the given point:
      // return a flag for having found it:
      //
      std::pair<short int, clipper::Coord_grid>
      search_for_skeleton_near(const clipper::Coord_orth &co,
			       const clipper::Xmap<int> &skel,
			       const clipper::Xmap<SkeletonTreeNode> &treenodemap) const; 
   };

}

short int
compare_cluster_centres(const coot::scored_skel_coord &a,
			const coot::scored_skel_coord &b);

#endif // HAVE_CALPHABUILD_HH
