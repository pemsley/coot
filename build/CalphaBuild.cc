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

#include <algorithm>  // sort std::sort on gcc 2.96
#include "math.h" // for round
#include "CalphaBuild.hh"

coot::CalphaBuild::CalphaBuild() {

//    std::cout << "Angle/Torsion probability table: " << std::endl;
//    for (float angle=0.0; angle<180.0; angle += 5.0) {
//       for (float torsion=0.0; torsion< 360.0; torsion += 10.0) {
// 	 float p = ai.prob_angle_torsion(angle, torsion);
// 	 std::cout << angle << "  " << torsion << "  " << p << std::endl;
//       }
//    }
//   std::cout  << "--" << std::endl;

   start_skel_depth = 10; 

} 

coot::CalphaBuild::CalphaBuild(int max_skeleton_search_depth) {

   start_skel_depth = max_skeleton_search_depth;

} 


// Here we do the look-ahead:
// 
std::vector <coot::scored_skel_coord>
coot::CalphaBuild::next_ca_by_skel(const std::vector<clipper::Coord_orth> &Previous_ca_positions,
				   const clipper::Coord_grid &coord_grid_start, 
				   short int use_coord_grid_start_flag,
				   float ca_bond_length,
				   const clipper::Xmap<int>  &skel_map,
				   const clipper::Xmap<float> &map,
				   float map_cut_off,
				   const clipper::Xmap<SkeletonTreeNode> &treenodemap) const {


   std::vector <coot::scored_skel_coord> scored_positions = 
      next_ca_internal(Previous_ca_positions, coord_grid_start, use_coord_grid_start_flag,
		       ca_bond_length, skel_map, map, map_cut_off, treenodemap, 1);


   std::vector <coot::scored_skel_coord> ahead_positions;
   for (unsigned int i=0; i<scored_positions.size(); i++) { 
       
      // 3 elements: need 1, 2
      // 2 elements: need 0, 1
      // 1 elements: need -, 0
      std::vector<clipper::Coord_orth> new_ca_positions;
      int psize = Previous_ca_positions.size();  
      if ((psize-2)>=0) 
	 new_ca_positions.push_back(Previous_ca_positions[psize-2]);
      new_ca_positions.push_back(Previous_ca_positions[psize-1]);

      // now add the newly created (trial) position
      new_ca_positions.push_back(scored_positions[i].position);

      ahead_positions = next_ca_internal(new_ca_positions,
					 scored_positions[i].near_grid_pos,
					 1,
					 ca_bond_length, skel_map, map, map_cut_off, treenodemap, 0);
      if (ahead_positions.size() > 0) {
// 	 std::cout << "combining scores : " << scored_positions[i].score
// 		   << " and " << ahead_positions[0].score << " (ahead) " << std::endl;
	 scored_positions[i].score *= ahead_positions[0].score; // it's sorted, best first
      }
   }

   std::sort(scored_positions.begin(), scored_positions.end(), compare_cluster_centres);
   return scored_positions;
}

// Here we get a list of coordinates for the next position, based on 3
// previous positions, the score for that position and an associated
// coord_grid.  There is no recursion going on here.
// 
std::vector <coot::scored_skel_coord>
coot::CalphaBuild::next_ca_internal(const std::vector<clipper::Coord_orth> &Previous_ca_positions,
				    const clipper::Coord_grid &coord_grid_start,
				    short int use_coord_grid_start_flag,
				    float ca_bond_length,
				    const clipper::Xmap<int>  &skel_map,
				    const clipper::Xmap<float> &map,
				    float map_cut_off, 
				    const clipper::Xmap<SkeletonTreeNode> &treenodemap,
				    int deep) const {

   std::vector<coot::scored_skel_coord> ssc;

   // previous atom (i.e. baton rotation centre) to skeleton point distance limits
   //
   float skel_dist_lim_max = 4.3; // A
   float skel_dist_lim_min = 2.5; // A la Oldfield 2003.
   
   int prev_size = Previous_ca_positions.size();

   if (prev_size > 0) {
      
      clipper::Coord_orth centre = Previous_ca_positions.back();
      clipper::Coord_frac centre_f = centre.coord_frac(map.cell());
      float box_radius = 5.0; // (min 3.8), 5 is slightly excessive.
      
      clipper::Coord_frac box0(
			       centre_f.u() - box_radius/map.cell().descr().a(),
			       centre_f.v() - box_radius/map.cell().descr().b(),
			       centre_f.w() - box_radius/map.cell().descr().c() );
      clipper::Coord_frac box1(
			       centre_f.u() + box_radius/map.cell().descr().a(),
			       centre_f.v() + box_radius/map.cell().descr().b(),
			       centre_f.w() + box_radius/map.cell().descr().c() );

      clipper::Grid_map grid( box0.coord_grid(map.grid_sampling()),
			      box1.coord_grid(map.grid_sampling()));

      std::vector<coot::scored_skel_coord> interesting;
      clipper::Xmap_base::Map_reference_coord ix( skel_map, grid.min() ), iu, iv, iw;
      for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() )
	 for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() )
	    for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {

	       if (skel_map[iw] > 0) {

		  if (map[iw] > map_cut_off) { 

		     clipper::Coord_grid c_g = iw.coord();
		     clipper::Coord_frac c_f = c_g.coord_frac(skel_map.grid_sampling());
		     clipper::Coord_orth c_o = c_f.coord_orth(skel_map.cell());
		     double plength = clipper::Coord_orth::length(c_o, centre);
		     if (plength > 2.5 && plength < 4.0) {

			// now check that we are away from the "centre-1" (start atom)
			//
			if (prev_size > 1) {
			
			   double slength =
			      clipper::Coord_orth::length(c_o,
							  Previous_ca_positions[prev_size-2]);

			   if (slength > 3.8) { // equilatral triangle -> 60 degrees,
			                        // we know we must have at least 70 degrees.
			      interesting.push_back(coot::scored_skel_coord(1.0, // dummy value
									    c_o,
									    c_g));
			   }
			} else {
			   // there is no previous atom to check (this is
			   // the first atom probably), so chuck it on
			   // the list.
			   interesting.push_back(coot::scored_skel_coord(1.0, // dummy value
									 c_o,
									 c_g));
			}
		     }
		  }
	       }
	    }

      // ssc = cluster(interesting, centre);
      ssc = interesting;

//       std::cout << "DEBUG:: There are " << ssc.size() << " clusters "
// 		<< " from " << interesting.size() << " interesting points "
// 		<< std::endl;

//       std::cout << "------ ssc scored positions ------" << std::endl;
//       for (int i=0; i<ssc.size(); i++) {
// 	 std::cout << i << " "
// 		   << ssc[i].score << "  "
// 		   << ssc[i].position.format() << "  "
// 		   << ssc[i].near_grid_pos.format() << std::endl;
//       }
//       std::cout << "- " << std::endl;

      for(unsigned int i=0; i<ssc.size(); i++) {
	 float score = score_position_by_angles(Previous_ca_positions,
						ssc[i].position);
	 ssc[i].score = score;
	 float skel_score = score_by_skel(coord_grid_start,
					  ssc[i].near_grid_pos,
					  use_coord_grid_start_flag, skel_map, map, map_cut_off,
					  treenodemap);
// 	 std::cout << "DEBUG: score_by_skel gives " <<  skel_score << std::endl;

	 int depth = int (skel_score/1000.0);
	 int skeleton_lev = int (skel_score - depth * 1000.0 + 0.1);

// 	 std::cout << "DEBUG: got depth " << depth
// 		   << " and skeleton_lev " << skeleton_lev << std::endl;
	 
	 // ssc[i].score *= skel_score;
	 ssc[i].skel_score = skel_score;
	 ssc[i].depth = depth;
	 ssc[i].skeletonization_level = skeleton_lev;

	 // experimental scoring:
	 ssc[i].score *= float(skeleton_lev*skeleton_lev)/100.0;
	 if (depth == 0) {
	    ssc[i].score *= 0.0001;
	 }

	 // score by branch point:

	 if (treenodemap.get_data(ssc[i].near_grid_pos).neighbs.size() > 2) {
	    ssc[i].score *= 10.0;
	 }
	    
	 // score by length

      }

      // extend to baton length (ca_bond_length)
      // 
      for(unsigned int i=0; i<ssc.size(); i++) {
	 ssc[i].position = extend_to_ca_bond_length_position(ssc[i].position, centre, ca_bond_length);
      }
   }

   // sort them 
   std::sort(ssc.begin(), ssc.end(), compare_cluster_centres);

   // Print them out
   // 
//    std::cout << "DEBUG:: scored sites: " << std::endl;
//    for(int i=0; i<ssc.size(); i++) {
//       std::cout << "    " << ssc[i].score << " " << ssc[i].skel_score << "  "
// 		<< ssc[i].position.format() << " "
// 		<< ssc[i].near_grid_pos.format() << std::endl;
//    }
//    std::cout << "---" << std::endl;
   

   return ssc;
}

// We get cluster centre, which are the average of skeleton points.
// Now these cluster centres are not at 3.8A (they are often closer
// than that).  We want the point at the end of 3.8 (ca_bond_length)
// of a vector that starts at a and goes through (strictly or points
// at) b.
// 
clipper::Coord_orth
coot::CalphaBuild::extend_to_ca_bond_length_position(const clipper::Coord_orth &b,
						     const clipper::Coord_orth &a,
						     float ca_bond_length) const {
   clipper::Coord_orth target_dir = b - a;
   return a + ca_bond_length*clipper::Coord_orth(target_dir.unit());
} 


float
coot::CalphaBuild::score_position_by_angles(const std::vector<clipper::Coord_orth> &Previous_ca_positions,
					    const clipper::Coord_orth &this_position) const {

   float score = 0.0001;
   if (Previous_ca_positions.size() == 3) {
      double this_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(this_position,
									  Previous_ca_positions[2],
									  Previous_ca_positions[1]));
      double this_torsion = clipper::Util::rad2d(clipper::Coord_orth::torsion(this_position,
									      Previous_ca_positions[2],
									      Previous_ca_positions[1],
									      Previous_ca_positions[0]));


      if (this_torsion < 0.0)
	 this_torsion += 360.0;
      
      score = ai.prob_angle_torsion(this_angle, this_torsion);
//       std::cout << "DEBUG:: angle: " << this_angle
// 		<< " torsion: " << this_torsion << " score: " << score << std::endl;
   } else {

      if (0) 
	 std::cout << "WARNING: default low geometry score because Previous_ca_positions.size()" 
		   << " is " << Previous_ca_positions.size() << std::endl;

   } 
      
   return score;
   
} 

// Here we (mis)-use the score column to hold the weights
// 
std::vector<coot::scored_skel_coord>
coot::CalphaBuild::cluster(const std::vector<coot::scored_skel_coord> &interesting,
			   const clipper::Coord_orth &starting_position) const {

   // relative to the starting position, fixed before returning
   // 
   std::vector<coot::scored_skel_coord> sum_centres;
   short int ifound;

   for (unsigned int i=0; i<interesting.size(); i++) {

      ifound = 0; 
      clipper::Coord_orth vec = interesting[i].position - starting_position;
      for (unsigned int j=0; j<sum_centres.size(); j++) {

	 clipper::Coord_orth this_centre_average =
	    (1/sum_centres[j].score)*sum_centres[j].position;

	 if (this_centre_average*vec/(sqrt(vec.lengthsq())*sqrt(this_centre_average.lengthsq())) > 0.95) {
	    sum_centres[j].score += 1.0;
	    sum_centres[j].position += vec;
	    ifound = 1;
	    sum_centres[j].near_grid_pos = interesting[i].near_grid_pos;
	 }
      }

      if (ifound == 0) {
	 scored_skel_coord ssc(1.0, vec, interesting[i].near_grid_pos);
	 sum_centres.push_back(ssc);
      }
   }

   for (unsigned int j=0; j<sum_centres.size(); j++) {
      sum_centres[j].position =
	 (1.0/sum_centres[j].score) * sum_centres[j].position + starting_position;
   }

   return sum_centres;
}


short int
compare_cluster_centres(const coot::scored_skel_coord &a,
			const coot::scored_skel_coord &b) {

   return a.score > b.score;
}


// start at start, let's see if we can get to target:
float
coot::CalphaBuild::score_by_skel(const clipper::Coord_grid &start,
				 const clipper::Coord_grid &this_grid,
				 short int use_coord_grid_start_flag,
				 const clipper::Xmap<int>  &skel_map,
				 const clipper::Xmap<float> &map,
				 float map_cut_off,
				 const clipper::Xmap<SkeletonTreeNode> &treenodemap) const {

   if (use_coord_grid_start_flag == 0)
      return 1.0; // can't do skeleton depth tracing


//    clipper::Coord_grid previous;
//    clipper::Coord_grid prev_prev;

   return score_by_skel_internal(start, this_grid, this_grid,
				 this_grid, skel_map, map, map_cut_off, treenodemap, start_skel_depth);
} 


float
coot::CalphaBuild::score_by_skel_internal(const clipper::Coord_grid &start,
					  const clipper::Coord_grid &target,
					  const clipper::Coord_grid &previous,
					  const clipper::Coord_grid &prev_prev,
					  const clipper::Xmap<int>  &skel_map,
					  const clipper::Xmap<float> &map,
					  float map_cut_off,
					  const clipper::Xmap<SkeletonTreeNode> &treenodemap,
					  int depth) const {

   // std::cout << "DEBUG: checking " << start.format() << " depth: " << depth << std::endl;
   
   clipper::Coord_grid cg;
   clipper::Coord_grid new_grid_point;
   float score;

   if (depth == 0) {

      return 0.0; 

   } else { 

      std::vector<clipper::Coord_grid> neighbs =  treenodemap.get_data(start).neighbs;
      int ns = neighbs.size();
//       std::cout << " skel_map.get_data(" << start.format() << ")" << " is "
// 		<< skel_map.get_data(start) << " (start)"  << std::endl;
//       std::cout << "There are " << ns << " neighbs for " << start.format() << std::endl;

      clipper::Coord_grid as_target = treenodemap.get_data(target).near_grid_point;
      
      for (int in=0; in<ns; in++) {
      
	 if (neighbs[in] == as_target) {
// 	    std::cout << "####### found target " << target.format()
// 		      << " at depth: " << depth << " level "
// 		      << skel_map.get_data(as_target) << std::endl;
	    // 	    return 100.0;
	    return depth * 1000.0 + skel_map.get_data(as_target);

	 }
      }

      for (int in=0; in<ns; in++) {
	 new_grid_point = neighbs[in];
// 	 std::cout << " skel_map.get_data(" << new_grid_point.format() << ")" << " is "
// 		   << skel_map.get_data(new_grid_point) << std::endl;
	 if (new_grid_point != previous) {
	    if (new_grid_point != prev_prev) { 
	       if (skel_map.get_data(new_grid_point) > 0) { // is skeleton
		  if (map.get_data(new_grid_point) > map_cut_off) { // matches graphics
		     score = score_by_skel_internal(new_grid_point,
						    as_target,
						    start,
						    previous,
						    skel_map,
						    map,
						    map_cut_off,
						    treenodemap,
						    depth-1);
// 		     if (score > (100.0 - 1.0) ) {
// 			return score;
// 		     }

		     if (score > 0.0)
			return score;
		  }
	       }
	    }
	 }
      }
   }
   return 0.0;
} 

std::pair<short int, clipper::Coord_grid>
coot::CalphaBuild::search_for_skeleton_near(const clipper::Coord_orth &co,
					    const clipper::Xmap<int> &skel,
					    const clipper::Xmap<SkeletonTreeNode> &treenodemap) const {

 
   clipper::Coord_grid grid_dummy(0,0,0); 
   std::pair<short int, clipper::Coord_grid> r(0, grid_dummy); // make the compiler happy
   r.first = 0;

   clipper::Coord_orth centre = co;
   clipper::Coord_frac centre_f = centre.coord_frac(skel.cell());
   float box_radius = 5.0; // (min 3.8), 5 is slightly excessive.
   
   clipper::Coord_frac box0(
			    centre_f.u() - box_radius/skel.cell().descr().a(),
			    centre_f.v() - box_radius/skel.cell().descr().b(),
			    centre_f.w() - box_radius/skel.cell().descr().c() );
   clipper::Coord_frac box1(
			    centre_f.u() + box_radius/skel.cell().descr().a(),
			    centre_f.v() + box_radius/skel.cell().descr().b(),
			    centre_f.w() + box_radius/skel.cell().descr().c() );
   
   clipper::Grid_map grid( box0.coord_grid(skel.grid_sampling()),
			   box1.coord_grid(skel.grid_sampling()));

   double best_length = 999.9;
   clipper::Xmap_base::Map_reference_coord ix( skel, grid.min() ), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() )
	 for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
	    
	    if (skel[iw] > 0) {
	       
	       clipper::Coord_grid c_g = iw.coord();
	       clipper::Coord_frac c_f = c_g.coord_frac(skel.grid_sampling());
	       clipper::Coord_orth c_o = c_f.coord_orth(skel.cell());
	       double plength = clipper::Coord_orth::length(c_o, centre);
	       if (plength < best_length ) {
		  r.first = 1;
		  r.second = c_g;
		  best_length = plength;
	       }
	    }
	 }

   return r;
} 
