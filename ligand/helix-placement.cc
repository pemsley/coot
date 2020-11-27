/* ligand/helix-placement.cc
 * 
 * Copyright 2005 The University of York
 * Author: Paul Emsley, Kevin Cowtan
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

// Portability gubbins
#ifndef _MSC_VER
#include <unistd.h> // for getopt(3)
#else
#define stat _stat
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <algorithm>

#include <clipper/clipper.h>
#include "clipper/ccp4/ccp4_map_io.h"
#include "residue_by_phi_psi.hh"
#include "helix-placement.hh"
#include "utils/coot-utils.hh"
#include "rigid-body.hh"
#include "mini-mol/mini-mol-utils.hh"
#include "db-main/db-strands.hh"

#include <sys/types.h> // for stating
#include <sys/stat.h>

coot::helix_placement_info_t
coot::helix_placement::place_alpha_helix_near(const clipper::Coord_orth &pt,
					      int n_residues,
					      float density_level_for_trim,
					      float b_factor,
					      float map_rmsd) const {

   clipper::Coord_orth ptc = pt;
   for (int i=0; i<10; i++) { 
      clipper::Coord_orth new_pt = move_helix_centre_point_guess(ptc, density_level_for_trim);
//       std::cout << "starting point: " << ptc.format() << std::endl;
//       std::cout << "  ending point: " << new_pt.format() << std::endl;
      // double d = clipper::Coord_orth::length(ptc, new_pt);
      // std::cout << "moved dist " << d << std::endl;
      ptc = new_pt;
      if ( (new_pt - pt).lengthsq() > 9.0) {
	 std::cout << "INFO:: helix placement centre point movement limit reached\n";
	 break;
      }
   }

   // std::cout << "centred position " << ptc.format() << std::endl;

   // So ptc is at the centre of the helix.  What are the eigen
   // vectors of the electron density about this point.

   float acell[6];
   acell[0] = xmap.cell().descr().a();
   acell[1] = xmap.cell().descr().b();
   acell[2] = xmap.cell().descr().c();
   acell[3] = clipper::Util::rad2d(xmap.cell().descr().alpha());
   acell[4] = clipper::Util::rad2d(xmap.cell().descr().beta());
   acell[5] = clipper::Util::rad2d(xmap.cell().descr().gamma());
   // mmmol.set_cell(a);
   std::string spacegroup_str_hm = xmap.spacegroup().symbol_hm();

   // a 20 residue helix is 27A long
   float search_radius = float(n_residues)/2.0;
   coot::eigen_info_t rtop_eigen_info = helix_eigen_system(ptc, search_radius);
   clipper::RTop_orth rtop_eigen = rtop_eigen_info.rtop;
// int best_eigenvalue_index = rtop_eigen_info.best_eigen_value_index;
//    std::cout << "DEBUG:: matrix\n" << rtop_eigen.format()
// 	     << "\n   has determinant "
// 	     << rtop_eigen.rot().det() << std::endl;

   coot::minimol::molecule mr; // part of returned value
   coot::helix_placement_info_t m = get_20_residue_helix_standard_orientation(n_residues,
									      b_factor);
   short int success = m.success;
   std::string failure_message = m.failure_message;

   if (m.mol[0].fragments.size() > 0 ) {

      // now apply rtop_eigen to m!
      //
      // we don't know which way round the helix is at this stage, of course.
      // So we will loop over 2 reorienation matrixes:

      std::vector<clipper::Mat33<double> > flip_rotated_matrices;
      flip_rotated_matrices.push_back(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1));
      flip_rotated_matrices.push_back(clipper::Mat33<double>(-1, 0, 0, 0, 1, 0, 0, 0, -1));

      clipper::Mat33<double> no_rotation    (1, 0,  0, 0, 1,  0, 0, 0, 1);
      clipper::Mat33<double> y_axis_rotation(0, 0, -1, 0, 1,  0, 1, 0, 0);
      clipper::Mat33<double> x_axis_rotation(1, 0,  0, 0, 0, -1, 0, 1, 0);
      clipper::RTop_orth no_rotation_op(no_rotation, clipper::Coord_orth(0,0,0));
      clipper::RTop_orth  y_axis_op(y_axis_rotation, clipper::Coord_orth(0,0,0));
      clipper::RTop_orth  x_axis_op(x_axis_rotation, clipper::Coord_orth(0,0,0));

      std::vector<clipper::RTop_orth> eigen_orientations;
      eigen_orientations.push_back(no_rotation_op);
      eigen_orientations.push_back(y_axis_op);
      eigen_orientations.push_back(x_axis_op);

      float score;
      std::vector<coot::scored_helix_info_t> eigens_scored_helices(3);
      std::vector<std::pair<float, float> > cb_density_info_for_scored_helices(3);

      for (int ieigen=0; ieigen<3; ieigen++) { 
	 std::vector<coot::scored_helix_info_t> this_eigen_scored_helices(2);
	 for (unsigned int irm=0; irm<flip_rotated_matrices.size(); irm++) {
	    coot::minimol::molecule rotated_helix = m.mol[0];
	    // apply flip (180 degrees (upsidedown)) matrices to rotated_helix then:
	    clipper::RTop_orth flip_rtop(flip_rotated_matrices[irm],
					 clipper::Coord_orth(0,0,0));
	    for(unsigned int ifrag=0; ifrag<rotated_helix.fragments.size(); ifrag++)
	       rotated_helix[ifrag].transform(flip_rtop);

	    // now search for best rotation round the helix axis:
	    //
	    for (float theta=0; theta<2*M_PI; theta += 0.2) {

	       coot::minimol::molecule theta_mol = rotated_helix;

	       double sin_t = sin(theta);
	       double cos_t = cos(theta);
	       clipper::Mat33<double> theta_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);
	       clipper::RTop_orth theta_rtop(theta_mat, clipper::Coord_orth(0,0,0));
	       for(unsigned int ifrag=0; ifrag<rotated_helix.fragments.size(); ifrag++)
		  theta_mol[ifrag].transform(theta_rtop);

	       for(unsigned int ifrag=0; ifrag<theta_mol.fragments.size(); ifrag++)
		  theta_mol[ifrag].transform(eigen_orientations[ieigen]);
	 
	       // finally apply eigen to shift to position in map.
	       for(unsigned int ifrag=0; ifrag<theta_mol.fragments.size(); ifrag++)
		  theta_mol[ifrag].transform(rtop_eigen);

	       // 	    std::string fitted_name = "fitted-helix-";
	       // 	    fitted_name += coot::util::int_to_string(irm);
	       // 	    fitted_name += "-";
	       // 	    fitted_name += coot::util::float_to_string(theta*180/M_PI);
	       // 	    fitted_name += ".pdb";
	       // 	    theta_mol.write_file(fitted_name);

	       score = score_helix_position(theta_mol);
	       if (score > this_eigen_scored_helices[irm].score) {
		  this_eigen_scored_helices[irm] = coot::scored_helix_info_t(theta_mol, score);
	       }
	       // std::cout << "score " << irm << " " << theta*180/M_PI << " " << score
	       // << std::endl;
	    }
	 }

	 //    scored_helices[0].mol.write_file("best-pre-fit-0.pdb"); 
	 //    scored_helices[1].mol.write_file("best-pre-fit-1.pdb");

	 // score = score_helix_position(scored_helices[0].mol);
	 // std::cout << "Score 0 pre rigid body fitting: " << score << std::endl;
	 // score = score_helix_position(scored_helices[1].mol);
	 // std::cout << "Score 1 pre rigid body fitting: " << score << std::endl;

	 coot::rigid_body_fit(&this_eigen_scored_helices[0].mol, xmap, map_rmsd);
	 coot::rigid_body_fit(&this_eigen_scored_helices[1].mol, xmap, map_rmsd);

	 score = score_helix_position(this_eigen_scored_helices[0].mol);
	 this_eigen_scored_helices[0].score = score;
	 // std::cout << "Score 0 post rigid body fitting: " << score << std::endl;
	 score = score_helix_position(this_eigen_scored_helices[1].mol);
	 this_eigen_scored_helices[1].score = score;
	 // std::cout << "Score 1 post rigid body fitting: " << score << std::endl;

	 // scored_helices[0].mol.write_file("best-post-fit-0.pdb"); 
	 // scored_helices[1].mol.write_file("best-post-fit-1.pdb");

	 std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> > p0 = 
	    decompose_helix_by_cbeta(this_eigen_scored_helices[0].mol);
	 coot::util::density_stats_info_t ds0 = score_atoms(p0.second);
	 std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> > p1 = 
	    decompose_helix_by_cbeta(this_eigen_scored_helices[1].mol);
	 coot::util::density_stats_info_t ds1 = score_atoms(p1.second);

	 std::pair<float, float> mv0 = ds0.mean_and_variance();
	 std::pair<float, float> mv1 = ds1.mean_and_variance();
	 if (mv0.second >= 0)
	    std::cout << "comparison of CB scores: " << mv0.first << " " << sqrt(mv0.second)
		      << std::endl;
	 if (mv1.second >= 0)
	    std::cout << "comparison of CB scores: " << mv1.first << " " << sqrt(mv1.second)
		      << std::endl;

	 if (mv0.first > mv1.first) {
	    eigens_scored_helices[ieigen] = this_eigen_scored_helices[0];
	    cb_density_info_for_scored_helices[ieigen] = mv0;
	 } else {
	    eigens_scored_helices[ieigen] = this_eigen_scored_helices[1];
	    cb_density_info_for_scored_helices[ieigen] = mv1;
	 }
	 // debug: write out the helices: forwards and back for each
	 // eigen, check the flips.
	 //       for (int i=0; i<2; i++) { 
	 // 	 std::string f = "post-fit";
	 // 	 f += "-eigen:";
	 // 	 f += coot::util::int_to_string(ieigen);
	 // 	 f += "-ori:";
	 // 	 f +=  coot::util::int_to_string(i);
	 // 	 f += ".pdb";
	 // 	 this_eigen_scored_helices[i].mol.write_file(f);
	 //       }
      }

      int best_index = 0;
      if (eigens_scored_helices[1].score > eigens_scored_helices[0].score)
	 best_index = 1;
      if (eigens_scored_helices[2].score > eigens_scored_helices[best_index].score)
	 best_index = 2;
       std::cout << "INFO:: CB density stats for best fit : "
 		<< cb_density_info_for_scored_helices[best_index].first << " "
		 << sqrt(cb_density_info_for_scored_helices[best_index].second) << " "
 		<< cb_density_info_for_scored_helices[best_index].first/sqrt(cb_density_info_for_scored_helices[best_index].second) << std::endl;

      float descrimination_ratio = 1.0; // seems a reasonable cut off, various testing.
      descrimination_ratio = 0.5; 
	 
      if ( (cb_density_info_for_scored_helices[best_index].first  > 0.0) &&
	   (cb_density_info_for_scored_helices[best_index].second > 0.0) &&
	   (cb_density_info_for_scored_helices[best_index].first/sqrt(cb_density_info_for_scored_helices[best_index].second) > descrimination_ratio) ) {
	 mr = eigens_scored_helices[best_index].mol;
	 mr.set_cell(acell);
	 mr.set_spacegroup(spacegroup_str_hm);
	 success = 1;
	 trim_and_grow(&mr, density_level_for_trim, b_factor); // modify mr
	 float final_score = score_helix_position(mr);
	 if (final_score < 0) {
	    success = 0;
	    failure_message = "Sorry, low overall density. Nowhere sensible for a helix.";
	 }
      } else {
	 success = 0;
	 failure_message = "Sorry, low density for CB positions - best fit ignored.";
      }
   } // m.mol.fragments.size() test.
   std::cout << "INFO:: helix addition success status: " << success << "\n";
   return coot::helix_placement_info_t(mr, success, failure_message);
}

// The Engine here is Kevin Cowtan's MR-style search.
//
// min_density_limit is typically 1 sigma
// 
coot::helix_placement_info_t
coot::helix_placement::place_alpha_helix_near_kc_version(const clipper::Coord_orth &pt,
							 int n_residues,
							 float min_density_limit,
							 float high_density_turning_point,
							 float b_factor,
							 float map_rmsd) const {

   clipper::Coord_orth ptc = pt;
   for (int i=0; i<10; i++) {
      float max_density_lim = high_density_turning_point; // min_density_limit is the 
						          // wrong name, I think.
      ptc = move_helix_centre_point_guess(ptc, max_density_lim);
   }

   // OK, let's limit the movement of the centre position to 1.5A.
   //
   double centre_moved_dist = sqrt((pt-ptc).lengthsq());
   double max_allowed_move = 1.0;
   if (centre_moved_dist > max_allowed_move) {
      clipper::Coord_orth v(max_allowed_move * (ptc.x()-pt.x())/centre_moved_dist,
			    max_allowed_move * (ptc.y()-pt.y())/centre_moved_dist,
			    max_allowed_move * (ptc.z()-pt.z())/centre_moved_dist);
      ptc = pt + v;
   }

   if (0) { // debugging initial moving
      std::cout << "=================== called with pt " <<  pt.format() << std::endl;
      std::cout << "=================== ptc now        " << ptc.format() << std::endl;
      std::cout << "=================== difference     " << clipper::Coord_orth(ptc-pt).format()
		<< std::endl;
   }
   

   float acell[6];
   acell[0] = xmap.cell().descr().a();
   acell[1] = xmap.cell().descr().b();
   acell[2] = xmap.cell().descr().c();
   acell[3] = clipper::Util::rad2d(xmap.cell().descr().alpha());
   acell[4] = clipper::Util::rad2d(xmap.cell().descr().beta());
   acell[5] = clipper::Util::rad2d(xmap.cell().descr().gamma());
   std::string spacegroup_str_hm = xmap.spacegroup().symbol_hm();

   std::cout << "DEBUG info:: xmap cell: "
	     << clipper::Util::rad2d(xmap.cell().descr().alpha()) << " " 
	     << clipper::Util::rad2d(xmap.cell().descr().beta()) << " " 
	     << clipper::Util::rad2d(xmap.cell().descr().gamma()) << "\n"; 
   std::cout << "DEBUG info:: xmap HM symbol: " << xmap.spacegroup().symbol_hm() << std::endl;

   std::vector<coot::minimol::molecule> mr(2); // part of returned value
   coot::helix_placement_info_t m =
      get_20_residue_helix_standard_orientation(n_residues, b_factor);
   short int success = m.success;
   std::string failure_message = m.failure_message;
   std::vector<std::pair<float, float> > cbeta_score_holder(2);
   // initializse the score holder to bad (helix rejected)
   cbeta_score_holder[0] = std::pair<float, float> (-1.0, -1.0);
   cbeta_score_holder[1] = std::pair<float, float> (-1.0, -1.0);

   if (m.mol[0].fragments.size() > 0 ) {

      // search model params
      double cyl_len = 6.0;  // half length
      double cyl_rad = 2.8;  // half diameter (radius)
      cyl_len = double(n_residues)/2.0; // PE parameter (for helices)

      // ------------- elided ----------------

      float max_density_limit = high_density_turning_point;
      clipper::RTop_orth ops_resultops =
	 find_best_tube_orientation(ptc, cyl_len, cyl_rad, max_density_limit);

      // ops[resultops] is the matrix which orients the to the helix
      // axis around the origin.
      std::vector<clipper::RTop_orth> ofm =
	 optimize_rotation_fit(m.mol[0], ops_resultops, ptc);

      // There are 2 orientation in ofm correspoding to the up and
      // down orientations of the initial helix.  We have already had
      // a guess at which would be the best, but we need to do a full
      // optimization of both directions before we can decide the best.
      //
      // We will make the "Helix", i.e. the first solution the one
      // with the best cb score.
      //
      // We will not add the second helix if it's cbeta score is not
      // at least 0.5 * the standard deviation of the CB density. More
      // below.

      for (int iofm=0; iofm<2; iofm++) { 
	 // now we need to apply ops[resultops] to the molecule.
	 //
	 if (1) { 
	    mr[iofm] = m.mol[0];

	    clipper::RTop_orth to_point(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), ptc);
	    for(unsigned int ifrag=0; ifrag<mr[iofm].fragments.size(); ifrag++)
	       mr[iofm][ifrag].transform(ofm[iofm]);
	    for(unsigned int ifrag=0; ifrag<mr[iofm].fragments.size(); ifrag++)
	       mr[iofm][ifrag].transform(ops_resultops);
	    for(unsigned int ifrag=0; ifrag<mr[iofm].fragments.size(); ifrag++)
	       mr[iofm][ifrag].transform(to_point);

	    clipper::Coord_orth d( 0.0, 0.0, 2.0 );
	    clipper::Coord_orth pt1(ptc + ops_resultops*d);
	    clipper::Coord_orth pt2(ptc - ops_resultops*d);

	    // debugging orientation:
	    if (0) { 
	     	 std::cout << " Helix axis: \n   (set-rotation-centre " << pt1.x() << " "
	     		   << pt1.y() << " "  << pt1.z() << ")\n";
	     	 std::cout << "     (place-atom-at-pointer)\n";
	     	 std::cout << " Helix axis: \n   (set-rotation-centre " << ptc.x() << " "
	     		   << ptc.y() << " "  << ptc.z() << ")\n";
	     	 std::cout << "     (place-atom-at-pointer)\n";
	     	 std::cout << " Helix axis: \n   (set-rotation-centre " << pt2.x() << " "
	     		   << pt2.y() << " "  << pt2.z() << ")\n";
	     	 std::cout << "     (place-atom-at-pointer)\n";
	    }

	    mr[iofm].set_cell(acell);
	    mr[iofm].set_spacegroup(spacegroup_str_hm);
	    // here set occupancies of CBetas to 0 so that they are not fitted.
	    mr[iofm].set_atom_occ(" CB ", 0.0);
	    coot::rigid_body_fit(&mr[iofm], xmap, map_rmsd);
	    // restore the occupancies of CBetas
	    mr[iofm].set_atom_occ(" CB ", 1.0);
	    success = 1;
	    trim_and_grow(&mr[iofm], min_density_limit, b_factor); // modify mr
	    float final_score = score_helix_position(mr[iofm]);
	    std::cout << "INFO:: Helix position score: " << final_score << "\n";
	    if (final_score < 0) { 
	       success = 0;
	       failure_message = "Low density for fitted helix.  Rejected.";
	    } else { 
	       // pair(other atoms, c-betas)
	       std::pair<std::vector<clipper::Coord_orth>,
		  std::vector<clipper::Coord_orth> > p1 =
		  decompose_helix_by_cbeta(mr[iofm]);
	       coot::util::density_stats_info_t ds = score_atoms(p1.second);
	       std::pair<float, float> mv0 = ds.mean_and_variance();
	       if (mv0.second >= 0) { 
		  std::cout << "INFO:: CB density scores: mean: " << mv0.first
			    << " stddev: " << sqrt(mv0.second) << std::endl;

		  cbeta_score_holder[iofm] = std::pair<float, float> (mv0.first, mv0.second);
		  
		  float descrimination_ratio = 1.0; // seems a reasonable cut off,
		  // various testing.
		  descrimination_ratio = 0.5; // liberal, let questionable
		  // helix through.
		  float mean_over_stddev = mv0.first/sqrt(mv0.second);
		  if (mean_over_stddev > descrimination_ratio) {
		     success = 1;
		  } else {
		     std::cout << "Failed CB density test" << std::endl;
		     failure_message = "Failed CB density test - rejected";
		     success = 1; // debugging
		  } 
	       } else {
		  std::cout      << "WARNING:: Bad fit for CBs - helix rejected\n" ;
		  failure_message = "WARNING:: Bad fit for CBs - helix rejected";
		  std::cout << "INFO:: --------------- density_stats_info_t ds mean and var "
			    << mv0.first << " " << mv0.second << std::endl;
		  std::cout << "INFO:: --------------- n " << ds.n << std::endl;
		  success = 0;
	       }
	    }
	 }
      }
   }
   if (success) {
      // If there was at least one helix that fitted tolerably well,
      // turn the helices in the vector if the score for the second
      // was beter.  Also filter the molecule vector if the poorer
      // molecule has bad stats (the mean cb density must be at least
      // half the standard deviation of its cbeta).  It also must have
      // a score at least 0.4*the better helix score.
      if (cbeta_score_holder[1].first > cbeta_score_holder[0].first) {
	 // second helix was better than first...
	 coot::minimol::molecule mtmp = mr[0];
	 mr[0] = mr[1];
	 mr[1] = mtmp;
	 // also swap the scores, we use them in the next test...
	 std::pair<float, float> tmp = cbeta_score_holder[0];
	 cbeta_score_holder[0] = cbeta_score_holder[1];
	 cbeta_score_holder[1] = tmp;
      }
      // trim out "Reverse" helices with really bad scores...
      if ((cbeta_score_holder[1].first < 0.0) ||
	  (cbeta_score_holder[1].first*2.0 < cbeta_score_holder[1].second) ||
	  (cbeta_score_holder[1].first*2.5 < cbeta_score_holder[0].first)) {
	 mr.resize(1);
      } 
   } 
   return coot::helix_placement_info_t(mr, success, failure_message);
}

// when density is above max_density_limit, start reducing the score
// 
clipper::RTop_orth 
coot::helix_placement::find_best_tube_orientation(clipper::Coord_orth ptc, double cyl_len,
						  double cyl_rad, float high_density_turning_point) const {

   std::cout << "INFO:: high_density_turning_point: " << high_density_turning_point << std::endl;

   // density map score for density d:
   //     below mdl: d
   //
   //     above mdl:
   //               2*mdl - d

   // make the search angles
   std::vector<clipper::RTop_orth> ops;
   int resultops;

   // make a list of rotation ops to try
   float step  =5.0; // test value
   float beta1 =90.0;   // beta search half range for mirror sym
   float alpha1=360.0;  // no gamma search for cyl sym
   // do a uniformly sampled search of orientation space
   for ( float bdeg=step/2; bdeg < beta1; bdeg += step ) {
      float beta = clipper::Util::d2rad(bdeg);
      // float s = anglim/clipper::Util::intf(sin(beta)*anglim/step+1);
      for ( float adeg=step/2; adeg < alpha1; adeg += step ) {
	 float alpha = clipper::Util::d2rad(adeg);
	 clipper::Euler_ccp4 euler( alpha, beta, 0.0 );
	 ops.push_back( clipper::RTop_orth(clipper::Rotation(euler).matrix()) );
      }
   }
   unsigned int ops_size = ops.size();
      
   // make inverse ops
   std::vector<clipper::RTop_orth> opsi( ops.size() );
   for (unsigned int j = 0; j < ops.size(); j++ ) opsi[j] = ops[j].inverse();

   typedef clipper::Xmap<float>::Map_reference_coord MRC;
   MRC i0, iu, iv, iw;
   std::vector<double> scores( ops.size() );
   std::vector<unsigned int> n_cyl_points( ops.size() );
   clipper::Cell cell          = xmap.cell();
   clipper::Grid_sampling grid = xmap.grid_sampling();
   clipper::Grid_range gr( cell, grid, sqrt(cyl_len*cyl_len+cyl_rad*cyl_rad) );
   clipper::Coord_grid g, g0, g1;
   clipper::Coord_orth c0, c1;
   float d; 
   for (unsigned int j = 0; j < ops.size(); j++ ) { 
      scores[j] = 0.0;
      n_cyl_points[j] = 0;
   }

   g = ptc.coord_frac(cell).coord_grid(grid);
   g0 = g + gr.min();
   g1 = g + gr.max();
   i0 = MRC( xmap, g0 );
   for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ) {
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	 for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    c0 = iw.coord_orth() - ptc;
	    for (unsigned int j = 0; j < ops_size; j++ ) {
	       c1 = opsi[j] * c0;
	       if ( fabs(c1.z()) < cyl_len &&
		    c1.x()*c1.x()+c1.y()*c1.y() < cyl_rad*cyl_rad ) {
		  d = xmap[iw];
		  double this_score = (d > high_density_turning_point ? 2*high_density_turning_point-d : d);
		  // std::cout << "this_score: " << this_score << " using d " << d << std::endl;
		  scores[j] += this_score;
		  // scores[j] += d;
		  n_cyl_points[j]++;
	       }
	    }
	 }
   }
   
   if (0) // debug
      for (unsigned int j = 0; j < ops.size(); j++)
	 std::cout << "orientation " << j << " score " << scores[j] << " " << n_cyl_points[j]
		   << " points" << std::endl;
   
   unsigned int jmax = 0;
   for (unsigned int j = 0; j < ops.size(); j++)
      if (scores[j] > scores[jmax])
	 jmax = j;
   resultops = jmax;
   return ops[resultops];
}


// Return the orientation for the best fit and the orientation for the
// best fit with the helix oriented the other way up.  The orientation
// matrices contain the flipping operator too.
// 
std::vector<clipper::RTop_orth>
coot::helix_placement::optimize_rotation_fit(const coot::minimol::molecule &helix,
					     const clipper::RTop_orth &axis_ori,
					     const clipper::Coord_orth &helix_point) const {

   //    helix.write_file("helix-pre-optim.pdb");
   std::vector<clipper::Mat33<double> > flip_rotated_matrices;
   flip_rotated_matrices.push_back(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1));
   flip_rotated_matrices.push_back(clipper::Mat33<double>(-1, 0, 0, 0, 1, 0, 0, 0, -1));
   float score;
   std::vector<float> best_score(2,-9999.9);
   std::vector<clipper::RTop_orth> best_rtop(2);

   for (unsigned int irm=0; irm<flip_rotated_matrices.size(); irm++) {
      coot::minimol::molecule rotated_helix = helix;
      // apply flip (180 degrees (upsidedown)) matrices to rotated_helix then:
      clipper::RTop_orth flip_rtop(flip_rotated_matrices[irm],
				   clipper::Coord_orth(0,0,0));
      for(unsigned int ifrag=0; ifrag<rotated_helix.fragments.size(); ifrag++)
	 rotated_helix[ifrag].transform(flip_rtop);


//       if (1) {
// 	 std::string s("rotated-helix-");
// 	 s += coot::util::int_to_string(irm);
// 	 s += ".pdb";
// 	 rotated_helix.write_file(s);
//       } 
      
      // now search for best rotation round the helix axis:
      //
      int theta_count = 0; // for debugging
      for (float theta=0; theta<2*M_PI; theta += 0.2) {
	 
	 coot::minimol::molecule theta_mol = rotated_helix;
	 
	 double sin_t = sin(theta);
	 double cos_t = cos(theta);
	 clipper::Mat33<double> theta_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);
	 clipper::RTop_orth theta_rtop(theta_mat, clipper::Coord_orth(0,0,0));
	 for(unsigned int ifrag=0; ifrag<rotated_helix.fragments.size(); ifrag++)
	    theta_mol[ifrag].transform(theta_rtop);
	 
	 // finally apply eigen to shift to position in map.
	 clipper::RTop_orth t(axis_ori.rot(), helix_point + axis_ori.trn());
	 for(unsigned int ifrag=0; ifrag<theta_mol.fragments.size(); ifrag++)
	    theta_mol[ifrag].transform(t);
	 score = score_helix_position(theta_mol);
	 if (score > best_score[irm]) { 
	    best_rtop[irm] = clipper::RTop_orth(theta_rtop * flip_rtop);
	    best_score[irm] = score;
	 }

	 if (0) { // debug
	    std::string s("theta-helix-flip-");
	    s += coot::util::int_to_string(irm);
	    s += "-theta-";
	    s += coot::util::int_to_string(theta_count);
	    s += ".pdb";
	    theta_mol.write_file(s, 20.0);
	 }
	 theta_count++;
      }

      if (0) { // debug
	 coot::minimol::molecule t = helix;
	 for(unsigned int ifrag=0; ifrag<t.fragments.size(); ifrag++)
	    t[ifrag].transform(best_rtop[irm]);
	 std::string s = "flipped-";
	 s += coot::util::int_to_string(irm);
	 s += ".pdb";
	 t.write_file(s, 20.0);
      } 
   }
   std::vector<clipper::RTop_orth> flip_best_rtop(2);

   flip_best_rtop[0] = best_rtop[0];
   flip_best_rtop[1] = best_rtop[1];
   if (best_score[1] > best_score[0]) {
      flip_best_rtop[0] = best_rtop[1];
      flip_best_rtop[1] = best_rtop[0];
   } 
   return flip_best_rtop;
}


clipper::Coord_orth
coot::helix_placement::move_helix_centre_point_guess(const clipper::Coord_orth &pt,
						     float density_max) const {

   int i=0;

   // First sample density in sphere around point and find the
   // weighted averge.

   float search_radius = 5.2; // A
   float search_radius_sq = search_radius * search_radius; 

   clipper::Coord_frac cf = pt.coord_frac(xmap.cell());
   clipper::Coord_frac box0(cf.u() - search_radius/xmap.cell().descr().a(),
			    cf.v() - search_radius/xmap.cell().descr().b(),
			    cf.w() - search_radius/xmap.cell().descr().c());

   clipper::Coord_frac box1(cf.u() + search_radius/xmap.cell().descr().a(),
			    cf.v() + search_radius/xmap.cell().descr().b(),
			    cf.w() + search_radius/xmap.cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(xmap.grid_sampling()),
			  box1.coord_grid(xmap.grid_sampling()));

   clipper::Coord_orth sum(0,0,0);
   float weight_sum = 0;
   float weight;
   int n_samples = 0;
   clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min()), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
      float ivcount = 0;
      // std::cout << "iv limits: " << iv.coord().v() << " " << grid.max().v() << std::endl;
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) { 
	 // std::cout << "   iv..." << ivcount << std::endl;
	 ivcount++;
	 for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
	    // sample a sphere
	    if ( (iw.coord().coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()) - pt).lengthsq() < search_radius_sq) {
	       // Get the position of that grid coord and weight it by xmap[iw]

	       float d = xmap[iw];
	       if (d > density_max)
		  d = density_max;

	       weight = d * d;
	       sum += clipper::Coord_orth(iw.coord_orth().x() * weight,
					  iw.coord_orth().y() * weight,
					  iw.coord_orth().z() * weight);
	       weight_sum += weight;
	       n_samples++;
	    }
	 }
      }
   }

   clipper::Coord_orth new_pt = pt;
   if (n_samples > 0) { 
      double one_over = 1/(weight_sum);
      // double one_over = 1/(double(n_samples));

      new_pt = (one_over * sum);

      if (0) { 
	 std::cout << "DEBUG:: starting point: " << pt.format() << std::endl;
	 std::cout << "DEBUG::     mean point: " << new_pt.format() << std::endl;
	 std::cout << "DEBUG:: moved... " << sqrt((pt - new_pt).lengthsq()) << std::endl;
	 std::cout << "DEBUG:: nsamples: " << n_samples << std::endl;
      }
      
      i = 1; // a good thing
   } 

   // std::cout << "DEBUG:: total moved " << sqrt((pt-new_pt).lengthsq()) << std::endl;
   return new_pt;
}

// Return the index of the highest eigen value and the raw rotation matix
// 
coot::eigen_info_t
coot::helix_placement::helix_eigen_system(const clipper::Coord_orth mean_pos,
					  float search_radius) const {

   // set up the covariance matrix 
   
   clipper::Matrix<double> mat(3,3);
   for (int ii=0; ii<3; ii++) 
      for (int jj=0; jj<3; jj++) 
	 mat(ii,jj) = 0.0; 
   
   // float search_radius = 8.0; // A
   float search_radius_sq = search_radius * search_radius; 

   clipper::Coord_frac cf = mean_pos.coord_frac(xmap.cell());
   clipper::Coord_frac box0(cf.u() - search_radius/xmap.cell().descr().a(),
			    cf.v() - search_radius/xmap.cell().descr().b(),
			    cf.w() - search_radius/xmap.cell().descr().c());

   clipper::Coord_frac box1(cf.u() + search_radius/xmap.cell().descr().a(),
			    cf.v() + search_radius/xmap.cell().descr().b(),
			    cf.w() + search_radius/xmap.cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(xmap.grid_sampling()),
			  box1.coord_grid(xmap.grid_sampling()));

   float d;
   int n_samples = 0;
   float sum_weights = 0;
   clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min()), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) { 
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) { 
	 for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
	    // sample a sphere
	    clipper::Coord_orth co = iw.coord_orth();
	    if ( (co - mean_pos).lengthsq() < search_radius_sq) {
	       // Get the position of that grid coord and weight it by xmap[iw]

	       d = xmap[iw];
	       if (d < 0)
		  d = 0;
	       d *= d;

	       mat(0,0) += d * (co.x() - mean_pos.x()) * (co.x() - mean_pos.x()); 
	       mat(0,1) += d * (co.x() - mean_pos.x()) * (co.y() - mean_pos.y()); 
	       mat(0,2) += d * (co.x() - mean_pos.x()) * (co.z() - mean_pos.z()); 
	       mat(1,0) += d * (co.y() - mean_pos.y()) * (co.x() - mean_pos.x()); 
	       mat(1,1) += d * (co.y() - mean_pos.y()) * (co.y() - mean_pos.y()); 
	       mat(1,2) += d * (co.y() - mean_pos.y()) * (co.z() - mean_pos.z()); 
	       mat(2,0) += d * (co.z() - mean_pos.z()) * (co.x() - mean_pos.x()); 
	       mat(2,1) += d * (co.z() - mean_pos.z()) * (co.y() - mean_pos.y()); 
	       mat(2,2) += d * (co.z() - mean_pos.z()) * (co.z() - mean_pos.z());
	       sum_weights += d;
	       n_samples++;
	    }
	 }
      }
   }
   // change mat
   std::vector<double> eigens = mat.eigen(false);

   int best_eigen_index = 0;
   if (eigens[1] > eigens[0])
      best_eigen_index = 1;
   if (eigens[2] > eigens[best_eigen_index])
      best_eigen_index = 2;

//    std::cout << " Biggest eigen is number " << best_eigen_index << std::endl;
//    std::cout << "eigens 0 " <<eigens[0] << " " << sqrt(eigens[0]/float(sum_weights))
// 	     << std::endl;
//    std::cout << "eigens 1 " <<eigens[1] << " " << sqrt(eigens[1]/float(sum_weights))
// 	     << std::endl;
//    std::cout << "eigens 2 " <<eigens[2] << " " << sqrt(eigens[2]/float(sum_weights))
// 	     << std::endl;
      // double eigen_descrimination_factor = eigens[best_eigen_index] * eigens[best_eigen_index] * eigens[best_eigen_index]/(eigens[0] * eigens[1] * eigens[2]);
//    std::cout << "Eigens descrimination factor: " << eigen_descrimination_factor
// 	     << std::endl;
   
//    std::cout << "n samples: " << n_samples << std::endl;
//    std::cout << "sum weight " << sum_weights << std::endl;

    clipper::Coord_orth e0(mat(0,0), mat(1,0), mat(2,0));
    clipper::Coord_orth e1(mat(0,1), mat(1,1), mat(2,1));
    clipper::Coord_orth e2(mat(0,2), mat(1,2), mat(2,2));
    //   clipper::Coord_orth e0(mat(0,0), mat(0,1), mat(0,2));
    //   clipper::Coord_orth e1(mat(1,0), mat(1,1), mat(1,2));
    //   clipper::Coord_orth e2(mat(2,0), mat(2,1), mat(2,2));
   e0 = sqrt(eigens[0]/float(sum_weights)) * e0;
   e1 = sqrt(eigens[1]/float(sum_weights)) * e1;
   e2 = sqrt(eigens[2]/float(sum_weights)) * e2;

   clipper::Coord_orth tip_0_1(mean_pos + e0);
   clipper::Coord_orth tip_0_2(mean_pos - e0);
   clipper::Coord_orth tip_1_1(mean_pos + e1);
   clipper::Coord_orth tip_1_2(mean_pos - e1);
   clipper::Coord_orth tip_2_1(mean_pos + e2);
   clipper::Coord_orth tip_2_2(mean_pos - e2);

   
//    std::cout << "(set-rotation-centre "
// 	     << mean_pos.x() << " "
// 	     << mean_pos.y() << " "
// 	     << mean_pos.z() << ") ; mean-pos" << std::endl;
//    std::cout << "(set-rotation-centre "
// 	     << tip_0_1.x() << " "
// 	     << tip_0_1.y() << " "
// 	     << tip_0_1.z() << ") ; tip x1" << std::endl;
//    std::cout << "(set-rotation-centre "
// 	     << tip_0_2.x() << " "
// 	     << tip_0_2.y() << " "
// 	     << tip_0_2.z() << ") ; tip x2" << std::endl;
//    std::cout << "(set-rotation-centre "
// 	     << tip_1_1.x() << " "
// 	     << tip_1_1.y() << " "
// 	     << tip_1_1.z() << ") ; tip y1" << std::endl;
//    std::cout << "(set-rotation-centre "
// 	     << tip_1_2.x() << " "
// 	     << tip_1_2.y() << " "
// 	     << tip_1_2.z() << ") ; tip y2" << std::endl;
//    std::cout << "(set-rotation-centre "
// 	     << tip_2_1.x() << " "
// 	     << tip_2_1.y() << " "
// 	     << tip_2_1.z() << ") ; tip z1" << std::endl;
//    std::cout << "(set-rotation-centre "
// 	     << tip_2_2.x() << " "
// 	     << tip_2_2.y() << " "
// 	     << tip_2_2.z() << ") ; tip z2" << std::endl;

   clipper::Mat33<double> mat33; 
   for (int j=0; j<3; j++)
      for (int i=0; i<3; i++)
	 mat33(i,j) = mat(i,j);
	    
   clipper::RTop_orth rtop(mat33, mean_pos);
   return coot::eigen_info_t(rtop, best_eigen_index, eigens);
}


coot::helix_placement_info_t
coot::helix_placement::get_20_residue_helix_standard_orientation(int nresidues,
								 float b_factor) const {

   coot::helix_placement_info_t m = get_20_residue_helix(nresidues);
   if (m.mol[0].fragments.size() > 0) {
//       std::cout << "DEBUG:: residue fragment residue range "
// 		<< m.mol[0][0].min_res_no() << " to "
// 		<< m.mol[0][0].max_residue_number()
// 		<< " requested: " << nresidues << std::endl;

      // We need to orient this arbitarily positioned fragment up the
      // z axis: So let's create 20 master points either side of the
      // origin along the z axis and use these as a reference frame to
      // rotate the CAs of the m molecule.
      //
      // We will then apply the transformation to all atoms in the m
      // (well, the 0th fragment of m).
      //
      std::vector<clipper::Coord_orth> z_axis_markers;
      std::vector<clipper::Coord_orth> frag_positions;
      int nres = nresidues; 
      for (int i=0; i<nres; i++) {
	 coot::minimol::atom a = m.mol[0][0][i+1][" CA "];
	 if (a.name != "FAIL") {
	    double z_pos =  (double(i) - double((nres+1)/2.0))*1.5;
	    clipper::Coord_orth pa(0, 0, z_pos);
	    z_axis_markers.push_back(pa);
	    frag_positions.push_back(a.pos);
	 }
      }
      // Note that we get garbage nans for our rtop when we put the
      // z_axis_markers down the x_axis (strangely enough).  So we put
      // the markers down the x axis and then multiply by:
      // ( 0 0 1) 
      // ( 0 1 0)
      // (-1 0 0)
      // No, we don't
      //
      clipper::RTop_orth rtop(frag_positions, z_axis_markers);
//       clipper::Mat33<double> y_axis_rotation(0, 0, -1, 0, 1, 0, 1, 0, 0);
//       clipper::RTop_orth y_axis_op(y_axis_rotation, clipper::Coord_orth(0,0,0));

      // we don't need (use) the x_axis_rotation
//       clipper::Mat33<double> x_axis_rotation(1, 0, 0, 0, 0, -1, 0, 1, 0);
//       clipper::RTop_orth x_axis_op(x_axis_rotation, clipper::Coord_orth(0,0,0));

      std::vector<coot::minimol::atom *> av = m.mol[0].select_atoms_serial();
      for (unsigned int iat=0; iat<av.size(); iat++) {
	 av[iat]->pos = av[iat]->pos.transform(rtop);
	 av[iat]->temperature_factor = b_factor;
      }
   }
   return m;
}


// Return molecule with no fragments on error reading file.  Try
// thowing an error, perhaps.
// 
coot::helix_placement_info_t
coot::helix_placement::get_20_residue_helix(int n_residues) const {

   std::string helix_filename = "theor-helix-70.pdb";
   std::string pkgdatadir = coot::package_data_dir();
   std::string full_path_helix_filename = pkgdatadir;
   full_path_helix_filename += "/";
   full_path_helix_filename += helix_filename;

   coot::minimol::molecule m;
   short int success = 0;
   std::string failure;
   
   struct stat buf;
   int status =  stat(full_path_helix_filename.c_str(), &buf);
   if (status == 0) {
      m.read_file(full_path_helix_filename);
      if (m.fragments.size() > 0) {
	 coot::minimol::fragment f;
	 // fill f with the first 20 residues of m;
	 if (m[0].max_residue_number() > n_residues) { 
	    for (int ires=1; ires<=n_residues; ires++) {
	       try { 
		  f.addresidue(m[0][ires], 0);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "ERROR:: get_20_residue_helix() " << rte.what() << std::endl;
	       }
	    }
	    m[0] = f;
	    success = 1;
	 }
      }
   } else {
      std::cout << "Error finding library theoretical helix pdb." << std::endl;
      failure = "Error finding library theoretical helix pdb.";
      success = 0;
   }
   return coot::helix_placement_info_t(m, success, failure);
}

coot::util::density_stats_info_t
coot::helix_placement::score_atoms(const std::vector<clipper::Coord_orth> &atom_pos) const {

   coot::util::density_stats_info_t dsi;
   float score;
   for (unsigned int i=0; i<atom_pos.size(); i++) {
      score = coot::util::density_at_point(xmap, atom_pos[i]);
      dsi.add(score);
   }
   return dsi;
}

float
coot::helix_placement::score_helix_position(const coot::minimol::molecule &m) const {

   float score = 0.0;
   for (unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
      for (int ires=m[ifrag].min_res_no(); ires<m[ifrag].max_residue_number();
	   ires++) {
	 for (unsigned int iat=0; iat<m[ifrag][ires].atoms.size(); iat++) {
	    score += coot::util::density_at_point(xmap, m[ifrag][ires][iat].pos) *
	       m[ifrag][ires][iat].occupancy;
	 }
      }
   }
   return score; 
}

coot::util::density_stats_info_t
coot::helix_placement::score_residue(const coot::minimol::residue &residue) const {

   coot::util::density_stats_info_t dsi;
   float score;
   for (unsigned int i=0; i<residue.atoms.size(); i++) {
       score = coot::util::density_at_point(xmap, residue.atoms[i].pos);
       dsi.add(score);
   }
   return dsi;
}

// pair(other atoms, c-betas)
std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> >
coot::helix_placement::decompose_helix_by_cbeta(coot::minimol::molecule &m) const {

   std::vector<clipper::Coord_orth> other_atoms;
   std::vector<clipper::Coord_orth> c_betas;
   for (unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
      for (int ires=m[ifrag].min_res_no(); ires<m[ifrag].max_residue_number();
	   ires++) {
	 for (unsigned int iat=0; iat<m[ifrag][ires].atoms.size(); iat++) {
	    if (m[ifrag][ires][iat].name == " CB ") {
	       c_betas.push_back(m[ifrag][ires][iat].pos);
	    } else {
	       other_atoms.push_back(m[ifrag][ires][iat].pos);
	    }
	 }
      }
   }
   return std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> > (other_atoms, c_betas);
} 


void
coot::helix_placement::discrimination_map() const {

   clipper::Xmap_base::Map_reference_index ix;
   clipper::Xmap<float> dmap;
   dmap.init(xmap.spacegroup(),
	     xmap.cell(),
	     xmap.grid_sampling());
   for (ix = xmap.first(); !ix.last(); ix.next() )  { // iterator index.
      coot::eigen_info_t e = helix_eigen_system(ix.coord_orth(), 10.0);
      double eigen_descrimination_factor = e.eigen_values[e.best_eigen_value_index] * e.eigen_values[e.best_eigen_value_index] * e.eigen_values[e.best_eigen_value_index]/(e.eigen_values[0] * e.eigen_values[1] * e.eigen_values[2]);
      dmap[ix] = eigen_descrimination_factor;
   }

   std::string filename("descrimination.map");
   clipper::CCP4MAPfile mapout;
   mapout.open_write(filename);
   mapout.export_xmap(dmap);
   mapout.close_write(); 
}


// modify f,
// return n trim flag, c trim flag
//
std::pair<int, int>
coot::helix_placement::trim_ends(coot::minimol::fragment *f, float min_density_limit) const { 

   // For each end of m, from the end residue (1 or 20) check if the residue is in 
   // positive density (well, above a particular level).  If it is, return, 
   // If it isn't, chop off the residue and repeat test for test residue.

   int ntrim = trim_end(f, 0, min_density_limit); // Nterm;
   int ctrim = trim_end(f, 1, min_density_limit); // Nterm;
    return std::pair<int, int> (ntrim, ctrim);
} 

// return if we did a chop or not
int
coot::helix_placement::trim_end(coot::minimol::fragment *f, short int end_type, float min_density_limit) const { 

   int ichop = 0;
   int dir_sign = 1;
   int start_resno = f->min_res_no();
   int last_resno  = f->max_residue_number();
   if (end_type == 1) { 
      dir_sign = -1;
      start_resno = f->max_residue_number();
      last_resno  = f->min_res_no();
   }

   for (int ires=start_resno; ires!=last_resno; ires+=dir_sign) { 
      coot::util::density_stats_info_t dsi = score_residue((*f)[ires]);
      std::pair<float, float> mv = dsi.mean_and_variance();
      if (dsi.n > 0.5) {  // this residue has atoms...
	   if (mv.first > min_density_limit) {
              break;
           } else {
              // chop out this residue;
              (*f)[ires].atoms.resize(0);
	      ichop = 1;
           }
      } 
   } 
   return ichop;
}


void
coot::helix_placement::build_on_N_end(coot::minimol::fragment *f,
				      float min_density_limit,
				      float b_factor) const { 

   // We only come here if we actually mean to do some building.  The
   // test for whether the end was trimmed must be performed in the
   // calling function.

   int ires_current_last = f->min_res_no();
   coot::minimol::residue current_last = (*f)[ires_current_last];

   bool found_CA = 0;
   bool found_N  = 0;
   bool found_C  = 0;
   clipper::Coord_orth ca;
   clipper::Coord_orth n;
   clipper::Coord_orth c;

//    std::cout << "Trying to build in N terminal residue " << ires_current_last-1
// 	     << std::endl;

   for (unsigned int iat=0; iat<current_last.atoms.size(); iat++) {
      if (current_last[iat].name == " CA ") {
	 found_CA = 1;
	 ca = current_last[iat].pos;
      }
      if (current_last[iat].name == " C  ") {
	 found_C = 1;
	 c = current_last[iat].pos;
      }
      if (current_last[iat].name == " N  ") {
	 found_N = 1;
	 n = current_last[iat].pos;
      }
   }
   if (found_CA && found_C && found_N) {
      coot::minimol::residue r = build_N_terminal_ALA(n, ca, c, ires_current_last-1, b_factor);
      // if that residue is in positive density then add it to f and try to add another res
      coot::util::density_stats_info_t dsi = score_residue(r);
      std::pair<float, float> mv = dsi.mean_and_variance();
      if (mv.first > min_density_limit) {
	 // Add on the C-Beta and add the residue to the fragment
	 std::pair<short int, clipper::Coord_orth> cbeta_info = coot::cbeta_position(r);
	 if (cbeta_info.first)
	    r.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 30.0);
	 try { 
	    f->addresidue(r, 0);
	    std::cout << "INFO build of N terminal residue " << ires_current_last-1
		      << " success" << std::endl;
	    build_on_N_end(f, min_density_limit, b_factor);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "ERROR:: build_on_N_end() " << rte.what() << std::endl;
	 } 
      } else {
	 std::cout << "INFO build of N terminal residue " << ires_current_last-1
		   << " failed (bad density fit)" << std::endl;
      }
   } else {
      std::cout << "Not all needed mainchain atoms found in prev. res" << std::endl;
   }
}

void
coot::helix_placement::build_on_C_end(coot::minimol::fragment *f,
				      float min_density_limit,
				      float b_factor) const { 

   // We only come here if we actually mean to do some building.  The
   // test for whether the end was trimmed must be performed in the
   // calling function.

   int ires_current_last = f->max_residue_number();
   coot::minimol::residue current_last = (*f)[ires_current_last];

//    std::cout << "Trying to build in C terminal residue " << ires_current_last-1
// 	     << std::endl;

   bool found_CA = 0;
   bool found_N  = 0;
   bool found_C  = 0;
   clipper::Coord_orth ca;
   clipper::Coord_orth n;
   clipper::Coord_orth c;

   for (unsigned int iat=0; iat<current_last.atoms.size(); iat++) {
      if (current_last[iat].name == " CA ") {
	 found_CA = 1;
	 ca = current_last[iat].pos;
      }
      if (current_last[iat].name == " C  ") {
	 found_C = 1;
	 c = current_last[iat].pos;
      }
      if (current_last[iat].name == " N  ") {
	 found_N = 1;
	 n = current_last[iat].pos;
      }
   }
   if (found_CA && found_C && found_N) {
      coot::minimol::residue r = build_C_terminal_ALA(n, ca, c, ires_current_last+1, b_factor);
      // if that residue is in positive density then add it to f and try to add another res
      coot::util::density_stats_info_t dsi = score_residue(r);
      std::pair<float, float> mv = dsi.mean_and_variance();
      if (mv.first > min_density_limit) { 
	 // Add on the C-Beta and add the residue to the fragment
	 std::pair<short int, clipper::Coord_orth> cbeta_info = coot::cbeta_position(r);
	 if (cbeta_info.first)
	    r.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 30.0);
	 try { 
	    f->addresidue(r, 0);
	    std::cout << "INFO build of N terminal residue " << ires_current_last-1
		      << " success" << std::endl;
	    build_on_C_end(f, min_density_limit, b_factor);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "ERROR:: build_on_C_end() " << rte.what() << std::endl;
	 } 
      } else {
	 std::cout << "INFO build of C terminal residue " << ires_current_last-1
		   << " failed (bad density fit)" << std::endl;
      }
   } else {
      std::cout << "Not all needed mainchain atoms found in prev. res" << std::endl;
   }
}

coot::minimol::residue
coot::helix_placement::build_N_terminal_ALA(const clipper::Coord_orth &prev_n,
					    const clipper::Coord_orth &prev_ca,
					    const clipper::Coord_orth &prev_c, 
					    int seqno,
					    float b_factor) const { 

   // theortical helix CA-N bond: -57.82 phi
   //                  CA-C bond: -47    psi
   float phi = -57.82; 
   float psi = -47.0;
   coot::minimol::residue r = coot::build_N_terminal_ALA(phi, psi, seqno, prev_n, prev_ca, prev_c,
							 b_factor);

   return r;
}

coot::minimol::residue
coot::helix_placement::build_C_terminal_ALA(const clipper::Coord_orth &prev_n,
					    const clipper::Coord_orth &prev_ca,
					    const clipper::Coord_orth &prev_c, 
					    int seqno,
					    float b_factor) const { 

   // theortical helix CA-N bond: -57.82 phi
   //                  CA-C bond: -47    psi
   float phi = -57.82; 
   float psi = -47.0;
   coot::minimol::residue r = coot::build_C_terminal_ALA(phi, psi, seqno, prev_n, prev_ca, prev_c,
							 b_factor);

   return r;
}


// Tinker with m
void
coot::helix_placement::trim_and_grow(minimol::molecule *m, float min_density_limit,
				     float b_factor) const {

   for (unsigned int ifrag=0; ifrag<m->fragments.size(); ifrag++) { 
      std::pair<int, int> trim_pair = trim_ends(&(*m)[ifrag], min_density_limit);
      if (trim_pair.first == 0) {
	 build_on_N_end(&(*m)[ifrag], min_density_limit, b_factor);
      } else {
	 std::cout << "N terminal of placed helix was trimmed" << std::endl;
      } 
      if (trim_pair.second == 0) {
	 build_on_C_end(&(*m)[ifrag], min_density_limit, b_factor);
      } else {
	 std::cout << "C terminal of placed helix was trimmed" << std::endl;
      }
   }
}


coot::helix_placement_info_t
coot::helix_placement::place_strand(const clipper::Coord_orth &pt, int strand_length,
				    int n_strand_samples, float sigma_level) {

   bool debug = false;
   float map_rmsd = sigma_level / 3; // hack (used in rigid body fitting shift scaling).
   coot::minimol::molecule rm;
   coot::helix_placement_info_t r(rm, 0, "Not Done");

   // What is the orientation of the strand?  Let's define a tube and search it
   //
   double cyl_len = 10.0;
   double cyl_rad = 1.0;

   float dd = 2.5 * sigma_level;
   clipper::RTop_orth best_op = find_best_tube_orientation(pt, cyl_len, cyl_rad, dd);

   clipper::RTop_orth op_plus_trans(best_op.rot(), pt);

   if (debug) {
      std::cout << "DEBUG:: best_op for strand orientation:\n" << best_op.format() << std::endl;
      std::cout << "DEBUG:: best_op with translation :\n" << op_plus_trans.format() << std::endl;
   }

   // get some strands of that length
   //
   coot::db_strands dbs;
   std::vector<coot::minimol::molecule> strands =
      dbs.get_reference_strands(n_strand_samples, strand_length);

   // for each strand, fit as rigid body, both directions and return the score.
   float best_score = -999.0;
   coot::minimol::molecule best_mol;
   for (unsigned int imol=0; imol<strands.size(); imol++) {
      std::cout << "Scoring fragment " << imol+1 << " of " << strands.size() << std::endl;
      coot::scored_helix_info_t info = fit_strand(strands[imol], op_plus_trans, imol, map_rmsd);
      if (info.score > best_score) {
	 std::cout << "   better score: " << info.score << std::endl;
	 best_score = info.score;
	 best_mol = info.mol;
      }
   }
   if (best_score > 0.0) {
      r.mol[0] = best_mol;
      r.success = 1;
      r.failure_message = "success";
      // copy over the cell and symmetry:
      float acell[6];
      acell[0] = xmap.cell().descr().a();
      acell[1] = xmap.cell().descr().b();
      acell[2] = xmap.cell().descr().c();
      acell[3] = clipper::Util::rad2d(xmap.cell().descr().alpha());
      acell[4] = clipper::Util::rad2d(xmap.cell().descr().beta());
      acell[5] = clipper::Util::rad2d(xmap.cell().descr().gamma());
      r.mol[0].set_cell(acell);
      r.mol[0].set_spacegroup(xmap.spacegroup().symbol_hm());
   }
   return r;
}

// imol is passed just for writting out score info.
// 
coot::scored_helix_info_t
coot::helix_placement::fit_strand(const coot::minimol::molecule &mol,
				  const clipper::RTop_orth &rtop,
				  int imol,
				  float map_rmsd) const {

   coot::scored_helix_info_t sir;
   float best_score = -9999.9;

   std::vector<coot::scored_helix_info_t> v =
      find_strand_candidates_by_shift_sampling(mol, rtop);
   std::cout << "Fitting " << v.size() << " shifted frag candidates from "
	     << " candidate fragment number " << imol+1 << std::endl;

   for (unsigned int iv=0; iv<v.size(); iv++) {
      if (v[iv].score < best_score*0.6) {
// 	 std::cout << "Reject - too far " << v[iv].score << " but refined to "
// 		   << score << std::endl;
      } else { 
	 coot::rigid_body_fit(&v[iv].mol, xmap, map_rmsd);
	 float score = score_helix_position(v[iv].mol);
	 if (score > best_score) {
	    best_score = score;
	    std::cout << "Got a better fit in fragment number " << imol+1
		      << " from " << v[iv].score
		      << " to " << score << std::endl;
	    sir = v[iv];
	    sir.score = score;
	 }
      }
   }
   return sir;
}

std::vector<coot::scored_helix_info_t>
coot::helix_placement::find_strand_candidates_by_shift_sampling(const coot::minimol::molecule &mol,
								const clipper::RTop_orth &rtop) const {
   
   // OK hold on to your hats here!

   // Return a vector of the "best-fitting" strands before we do rigid
   // body refinement.  i.e. there's no need to do rigid body
   // refinement on all strands, just the best-fitting ones.  So here
   // we return the best-fitting ones. And in the calling function of
   // this routine, we do the rigid-body fitting of each strand and
   // return the best one.  Return all strands that have a score of
   // more than half the top score.

   float top_cut = 0.5; // above this gets through
   top_cut = 0.8; 

   // rtop transform object near the origin (oriented along z) to
   // match the orientation and position of the place of interest.
   // That transformation happens last.  We will move things around
   // the origin first.
   //
   // We need to do some searching of the density, so we need
   //
   // 1) A forward/backwards transformation (flip)
   // 2) a rotation search around the strand (z) axis
   // 3) a translation search along the strand axis

   float z_shift_step_size = 0.35; // the step size along the strand
				   // axis, to get the zig-zag in the
				   // right place.

   std::vector<coot::scored_helix_info_t> scored_strand_vector_working;
   std::vector<coot::scored_helix_info_t> scored_strand_vector_returned;
   
   coot::scored_helix_info_t scored_strand;
   coot::minimol::molecule mc = mol;

   // flip up/down (reverse direction search)
   std::vector<clipper::Mat33<double> > flip_rotated_matrices;
   flip_rotated_matrices.push_back(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1));
   flip_rotated_matrices.push_back(clipper::Mat33<double>(-1, 0, 0, 0, 1, 0, 0, 0, -1));
   for (unsigned int iflip=0; iflip<flip_rotated_matrices.size(); iflip++) {
//       if (iflip == 0) {
// 	 std::cout << "   testing forward direction..." << std::endl;
//       } else { 
// 	 std::cout << "   testing reverse direction..." << std::endl;
//       }
      
      coot::minimol::molecule flip = mc;
      clipper::RTop_orth flip_rtop(flip_rotated_matrices[iflip],
				   clipper::Coord_orth(0,0,0));
      flip.transform(flip_rtop);

      for (float theta=0; theta<2*M_PI; theta += 0.3) {
	 
	 double sin_t = sin(theta);
	 double cos_t = cos(theta);
	 clipper::Mat33<double> theta_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);
	 clipper::RTop_orth theta_rtop(theta_mat, clipper::Coord_orth(0,0,0));
	 coot::minimol::molecule theta_mol = flip;
	 theta_mol.transform(theta_rtop);

	 for (float z_shift= -2.0; z_shift<2.0; z_shift+=z_shift_step_size) { 
	    clipper::RTop_orth z_shift_rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
					    clipper::Coord_orth(0, 0, z_shift));
	    coot::minimol::molecule z_shift_mol = theta_mol;
	    z_shift_mol.transform(z_shift_rtop);

	    // Now the transformation we first thought of...
	    z_shift_mol.transform(rtop); 
	    float score = score_helix_position(z_shift_mol);

	    // debug
	    if (0) { 
	       std::string filename = "f-";
	       filename += "flip:";
	       filename += coot::util::int_to_string(iflip);
	       filename += "-theta:";
	       filename += coot::util::float_to_string(theta);
	       filename += "-shift:";
	       filename += coot::util::float_to_string(z_shift);
	       filename += ".pdb";
	       z_shift_mol.write_file(filename, 30.0);
	    }

	    scored_strand.score = score;
	    scored_strand.mol = z_shift_mol;
	    scored_strand_vector_working.push_back(scored_strand);

	 }
      }
   }

   // now we have lots of vectors in scored_strand_vector_working. Sort them.
   std::sort(scored_strand_vector_working.begin(),
 	     scored_strand_vector_working.end(),
 	     compare_scored_strands);

   float top_score = 0.0;
   for (unsigned int i=0; i<scored_strand_vector_working.size(); i++) {
      // std::cout << "    " << i << " " << scored_strand_vector_working[i].score << std::endl;
      if (scored_strand_vector_working[i].score > top_score) {
	 top_score = scored_strand_vector_working[i].score;
      }
   }

   if (top_score > 0.0) { 
      for (unsigned int i=0; i<scored_strand_vector_working.size(); i++) {
	 if (scored_strand_vector_working[i].score > top_score*top_cut) {
	    scored_strand_vector_returned.push_back(scored_strand_vector_working[i]);
	 }
      }
   }
   
   return scored_strand_vector_returned;
} 


bool
coot::compare_scored_strands(const scored_helix_info_t &a,
			     const scored_helix_info_t &b) {
   return (a.score > b.score); 
} 
