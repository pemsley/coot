/* ligand/coot-fffear.cc
 * 
 * Copyright 2005 by The University of York
 * Copyright 2006 by The University of York
 * Copyright 2016 by Medical Research Council
 * Author: Kevin Cowtan and Paul Emsley
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */

// Much this class was copied from the exmple code of clipper and
// then fixed with the help of Kevin Cowtan.


// CXX11 version
#include <thread>

#include "clipper/core/map_interp.h"
#include "clipper/core/hkl_compute.h"

#include "coot-map-heavy.hh"
#include "coot-map-utils.hh"
#include "coot-coord-utils.hh"

// for FFfearing
#include "clipper/core/atomsf.h"
#include "clipper/core/rotation.h"
#include "clipper/contrib/fffear.h"
#include "peak-search.hh"


coot::util::fffear_search::fffear_search(mmdb::Manager *mol, int SelectionHandle, const clipper::Xmap<float> &xmap, float angular_resolution, bool translation_search_only) {

   std::pair<clipper::Coord_orth, clipper::Coord_orth> e = extents(mol, SelectionHandle);
   float bx  = e.second.x() - e.first.x();
   float by  = e.second.y() - e.first.y();
   float bz  = e.second.z() - e.first.z();

   float box_size = bx;
   float min_size = bx; 
   if (by > bx)
      box_size = by;
   if (bz > box_size)
      box_size = bz;
   // now similar for min_size:
   if (by < bx)
      min_size = by;
   if (bz < min_size)
      min_size = bz;

   // We have the min box size in an arbitary relative orientation of
   // the molecule to the orothogonal axes, if we wanted to do this
   // correctly we'd orient the molecule along the major, minor axes
   // before we put a box round it.

   // So multiply by a "consevative" ratio.  It might be wrong for a
   // extreme cases of long thin molecule such as errr... P69
   // pertactin.
   
   min_molecule_radius_ = 0.66*min_size/2.0; // set class data for use in
                                                   // peak filtering

   // get mid point of molecule here and take that away from coords of
   // the points that generate the nxmap and the mask.  So mid_point
   // it passed not lower_left.

   box_size += 3.0; // border.  Is this needed?  Yes, I think it is.
                    // Electron density spills over from peripheral
                    // atoms.

   // about the mid point of the molecule.  The nxmap and nxmap_mask
   // constructors constuct a grid about the origin (they do in this
   // case).
   mid_point_ = clipper::Coord_orth((e.first.x() + e.second.x())*0.5,
                                    (e.first.y() + e.second.y())*0.5,
                                    (e.first.z() + e.second.z())*0.5);

   clipper::Coord_orth low_left = e.first;
   std::cout << "Coords extents " << e.first.format() << " " << e.second.format() << "\n";
   std::cout << "Mid point: " << mid_point_.format() << std::endl;
   std::cout << "Box size: " << box_size << "\n";
      
   // clipper::Grid_sampling grid_sampling = xmap.grid_sampling();
   clipper::Resolution reso(2.0);
   clipper::Grid_sampling grid_sampling(xmap.spacegroup(), xmap.cell(), reso);

   // clear results at start
   results.init(clipper::Spacegroup::p1(), xmap.cell(), xmap.grid_sampling()); // fixed grid_sampling.
   std::pair<float, int> p(9999999999.9, -1); // low fffear results are good.
   clipper::Xmap<float>::Map_reference_index ix;
   for ( ix = results.first(); !ix.last(); ix.next() )
      results[ix] = p;

   if (translation_search_only)
      ops.push_back(clipper::RTop_orth(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
                                       clipper::Coord_orth(0,0,0)));
   else
      generate_search_rtops(angular_resolution); // fill class member ops

   // sets of ops indices for each thread
   std::vector<std::vector<unsigned int> > ops_index_set;
   unsigned int n_threads = get_max_number_of_threads();
   if (n_threads > 0) {
      ops_index_set.resize(n_threads);
      int unsigned thread_idx = 0;
      for (std::size_t i=0; i<ops.size(); i++) {
         ops_index_set[thread_idx].push_back(i);
         thread_idx++;
         if (thread_idx == n_threads)
            thread_idx = 0;
      }
   }

   clipper::Grid_range grid_extent(xmap.cell(), grid_sampling, box_size);

   // For each grid point store a best score and the rtop that
   // corresponds to it.
   // 
   // clipper::Xmap<std::pair<float, int> > results; // member data.

   std::cout << "INFO searching " << ops.size() << " orientations\n";
   std::cout << "Searched map grid sampling: " << xmap.grid_sampling().format()
             << "\n        i.e. "
             << xmap.cell().a()/float(xmap.grid_sampling().nu()) << " "
             << xmap.cell().b()/float(xmap.grid_sampling().nv()) << " " 
             << xmap.cell().c()/float(xmap.grid_sampling().nw()) 
             << " A/grid\n";
   std::cout << "Model map grid sampling: " << grid_sampling.format()
             << "\n        i.e. "
             << xmap.cell().a()/float(grid_sampling.nu()) << " "
             << xmap.cell().b()/float(grid_sampling.nv()) << " " 
             << xmap.cell().c()/float(grid_sampling.nw()) << " A/grid\n";
   int icount = 0;

   nxmap.init(xmap.cell(), grid_sampling, grid_extent);
   int npoints1 =  fill_nxmap(mol, SelectionHandle, mid_point_);
   nxmap_mask.init(xmap.cell(), grid_sampling, grid_extent);
   int npoints2 = fill_nxmap_mask(mol, SelectionHandle, mid_point_);
   std::cout << "Initialized search map and mask with grid sampling: "
             << grid_sampling.format() << std::endl;

   if (npoints1 == 0) {
      std::cout << "No point" << std::endl;

   } else {
      
      std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
      float xmap_mean  = mv.first;
      float xmap_stddev = mv.second;

      if (xmap_stddev > 0.0) { 
   
         post_process_nxmap(xmap_mean, xmap_stddev); // make nxmap match these.

         float sum_non_masked_nxmap = 0.0;
         int   n_non_masked = 0;
         clipper::NXmap<float>::Map_reference_index inx;
         for (inx = nxmap_mask.first(); !inx.last(); inx.next()) { 
            if (nxmap_mask[inx] > 0.0) {
               sum_non_masked_nxmap += nxmap[inx];
               n_non_masked++;
            }
         }

         if (n_non_masked == 0) {
            std::cout << "VERY STRANGE:: No non masked points!" << std::endl;
         } else { 
            //          std::cout << "DEBUG:: sum non masked map: " << sum_non_masked_nxmap << " average over "
            //                    << n_non_masked << " points is " << sum_non_masked_nxmap/float(n_non_masked)
            //                    << std::endl;
         }

         std::vector<std::thread> threads;

         for (std::size_t ith=0; ith<n_threads; ith++)
            threads.push_back(std::thread(fffear_search_inner_threaded,
                                          std::cref(xmap), std::cref(nxmap), std::cref(nxmap_mask),
                                          std::cref(ops), ops_index_set.at(ith), std::ref(results)));

         // wait
         for (std::size_t ith=0; ith<n_threads; ith++)
            threads.at(ith).join();
      }
   }
}

// static (of course)
//
// Fill results
void
coot::util::fffear_search::fffear_search_inner_threaded(const clipper::Xmap<float> &xmap,
                                                        const clipper::NXmap<float> &nxmap,
                                                        const clipper::NXmap<float> &nxmap_mask,
                                                        const std::vector<clipper::RTop_orth> &ops,
                                                        const std::vector<unsigned int> &ops_idx_set,
                                                        clipper::Xmap<std::pair<float, int> > &results) {

   unsigned int icount = 0;
   for (unsigned int i=0; i<ops_idx_set.size(); i++) {
      unsigned int iop = ops_idx_set[i];
      clipper::Xmap<float> r1;
      r1.init(clipper::Spacegroup::p1(), xmap.cell(), xmap.grid_sampling()); // fixed to be consistent with above results init.
      clipper::FFFear_fft<float> search(xmap);
      clipper::NX_operator nxop(xmap, nxmap, ops[iop]);
      search(r1, nxmap, nxmap_mask, nxop); // fill r1 with fffear search values at this orientation

      clipper::Xmap<float>::Map_reference_index ix;
      for (ix = r1.first(); !ix.last(); ix.next() ) {
         if (r1[ix] < results[ix].first) {
            results[ix].first = r1[ix];
            results[ix].second = iop;
         }
      }

      icount++;

      std::cout.flush();
      if (icount == 50) {
         std::cout << " " <<100*(float(iop)/float(ops.size())) << "%";
         std::cout.flush();
         icount = 0;
      }
   }
}


int
coot::util::fffear_search::fill_nxmap(mmdb::Manager *mol, int SelectionHandle,
                                      const clipper::Coord_orth &mid_point) {

   // I use as a base clipper's contrib/edcalc.cpp
   // template<class T> bool EDcalc_iso<T>::operator() ( NXmap<T>& nxmap, const Atom_list& atoms ) const
   //

   mmdb::PAtom *atom_selection;
   int n_atoms;
   mol->GetSelIndex(SelectionHandle, atom_selection, n_atoms);

   clipper::ftype radius_ = 2.5;

   clipper::NXmap<float>::Map_reference_index im;
   nxmap = 0.0;

   int n_points = 0;
   clipper::Coord_grid g0, g1;
   g0 = clipper::Coord_map( nxmap.operator_orth_grid().rot() *
                            clipper::Vec3<>(radius_,radius_,radius_) ).coord_grid();
   clipper::Grid_range gd( -g0, g0 );
   clipper::Grid_range box( clipper::Coord_grid(0,0,0),
                            clipper::Coord_grid(nxmap.grid()) - clipper::Coord_grid(1,1,1) );
   clipper::NXmap<float>::Map_reference_coord i0, iu, iv, iw;
   for ( int i = 0; i < n_atoms; i++ )
      if ( atom_selection[i] ) {
         clipper::Coord_orth p(atom_selection[i]->x, atom_selection[i]->y, atom_selection[i]->z);
         p -= mid_point;
         clipper::AtomShapeFn sf( p, std::string(atom_selection[i]->element),
                                  atom_selection[i]->tempFactor,
                                  atom_selection[i]->occupancy);
         g0 = nxmap.coord_map(p).coord_grid() + gd.min();
         g1 = nxmap.coord_map(p).coord_grid() + gd.max();
         i0 = clipper::NXmap<float>::Map_reference_coord( nxmap, g0 );
         float r;
         for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
            for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
               for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
                  if ( box.in_grid( iw.coord() ) ) { 
                     r = sf.rho( iw.coord_orth() );
                     if (clipper::Util::isnan(r)) {
                        // std::cout << "ERROR:: adding nan! " << std::endl;
                        // caused by B factor of 0.00
                     } else {
                        nxmap[iw] += r;
                     } 
                     n_points++;
                  }
      }
   std::cout << "INFO:: Number of non-zero points in atom search map: " << n_points
             << std::endl;
   
   // debugging check for nans:
   clipper::NXmap<float>::Map_reference_index inx;
   float d;
   int n_nan = 0;
   int n_nx_points = 0;
   for (inx = nxmap.first(); !inx.last(); inx.next()) {
      n_nx_points++;
      d = nxmap[inx];
      if (clipper::Util::isnan(d)) {
         n_nan++;
         // If the following line is commented out with compiler flag
         // -O2, then it doesn't get executed (compiler bug, I'd
         // imagine)
         // std::cout << "ERROR:: post fill nan in nxmap " << std::endl;
      }
   }
   if (n_nan > 0) { 
      std::cout << "----:: " << n_nan << " of " << n_nx_points
                << " map points were nans" << std::endl;
      std::cout << "----:: " << n_points << " were set to density values" << std::endl;
   }
   return n_points;
}

int
coot::util::fffear_search::fill_nxmap_mask(mmdb::Manager *mol, int SelectionHandle,
                                           const clipper::Coord_orth &mid_point) {

   // I use as a base clipper's contrib/edcalc.cpp
   // template<class T> bool EDcalc_mask<T>::operator() ( NXmap<T>& nxmap, const Atom_list& atoms ) const

   //
   mmdb::PAtom *atom_selection;
   int n_atoms;
   mol->GetSelIndex(SelectionHandle, atom_selection, n_atoms);

   int n_points = 0;
   clipper::ftype radius_ = 2.5;
   nxmap_mask = 0.0;
   clipper::Coord_grid g0, g1;
   g0 = clipper::Coord_map( nxmap.operator_orth_grid().rot() *
                            clipper::Vec3<>(radius_,radius_,radius_) ).coord_grid();
   clipper::Grid_range gd( -g0, g0 );
   clipper::Grid_range box( clipper::Coord_grid(0,0,0),
                            clipper::Coord_grid(nxmap.grid()) - clipper::Coord_grid(1,1,1) );
   clipper::NXmap<float>::Map_reference_coord i0, iu, iv, iw;
   for ( int i = 0; i < n_atoms; i++ )
      if ( atom_selection[i] ) {
         clipper::Coord_orth xyz(atom_selection[i]->x, atom_selection[i]->y, atom_selection[i]->z);
         xyz -= mid_point;
         g0 = nxmap.coord_map( xyz ).coord_grid() + gd.min();
         g1 = nxmap.coord_map( xyz ).coord_grid() + gd.max();
         i0 = clipper::NXmap<float>::Map_reference_coord( nxmap, g0 );
         for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
            for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
               for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
                  if ( box.in_grid( iw.coord() ) )
                     if ( ( xyz - iw.coord_orth() ).lengthsq() < radius_*radius_ ) { 
                        nxmap_mask[iw] = 1.0;
                        n_points++;
                     }
      }

   std::cout << "INFO:: Number of non-zero points in search mask : " << n_points
             << std::endl;
   return n_points;
}


clipper::RTop_orth
coot::util::fffear_search::get_best_transformation() const {

   clipper::RTop_orth r;


   return r;
}

// Kevin Cowtan writes this function:
//
void
coot::util::fffear_search::generate_search_rtops(float angular_resolution) {

   // Something like this will make you a nice list of RTop_orths

   // make a list of rotation ops to try
   float gamma1=360.0;
   float beta1 =180.0;
   float alpha1=360.0;
   // do a uniformly sampled search of orientation space
   float anglim = clipper::Util::min(alpha1,gamma1);
   float step = angular_resolution;
   for ( float bdeg=step/2; bdeg < beta1; bdeg += step ) {
     float beta = clipper::Util::d2rad(bdeg);
     float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
     float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
     for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
       for ( float thmi=smi/2; thmi < 360.0; thmi += smi ) {
        float alpha = clipper::Util::d2rad(0.5*(thpl+thmi));
        float gamma = clipper::Util::d2rad(0.5*(thpl-thmi));
        clipper::Euler_ccp4 euler( alpha, beta, gamma );
        ops.push_back( clipper::RTop_orth(clipper::Rotation(euler).matrix()) );
       }
   }
}


clipper::Xmap<float>
coot::util::fffear_search::get_results_map() const {


   clipper::Xmap<float> r(results.spacegroup(), results.cell(), results.grid_sampling());
   clipper::Xmap<float>::Map_reference_index ix;
   for (ix = results.first(); !ix.last(); ix.next())
      r[ix] = results[ix].first; 

   return r;
}

std::vector<std::pair<float, clipper::RTop_orth> >
coot::util::fffear_search::scored_orientations() const {

   std::vector<std::pair<float, clipper::RTop_orth> > r;
      
   clipper::Xmap<float> xmap = get_results_map();
   coot::peak_search ps(xmap);
   float n_sigma = 4;

   // peak search should find peaks.

   // returns garbage mri for all but the first for some reason
//    std::vector<std::pair<clipper::Xmap<float>::Map_reference_index, float> > peaks =
//       ps.get_peak_map_indices(xmap, n_sigma);

   std::vector<std::pair<clipper::Coord_grid, float> > peaks =
      ps.get_peak_grid_points(xmap, n_sigma);

   std::cout << "DEBUG: get_peak_grid_points returned " << peaks.size() << " peaks\n";

   if (peaks.size() > 4) {
      for (int i=0; i<4; i++)
         std::cout << "in scored_orientations "
                   << peaks[i].first.format() << " " << peaks[i].second << "\n";
   }

   for (unsigned int i=0; i< peaks.size(); i++) { 
      // std::cout << "Peak: " << peaks[i].second << " at " << peaks[i].first.coord().format()
      // << std::endl;
      // std::cout << "Peak: " << peaks[i].second << std::endl;
      int i_rtop_index = results.get_data(peaks[i].first).second;
      if (i_rtop_index >= 0) { 
         if (i_rtop_index < int(ops.size())) {
            clipper::Coord_orth trn =
               peaks[i].first.coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell());
//             std::cout << "debug:: converted " << peaks[i].first.format()
//                       << " to " << trn.format() << std::endl;
            clipper::RTop_orth rtop_from_map(ops[i_rtop_index]);
            clipper::RTop_orth rtop_with_trans(rtop_from_map.rot(), trn);
            std::pair<float, clipper::RTop_orth> p(peaks[i].second, rtop_with_trans);
            r.push_back(p);
         } else {
            std::cout << "ERROR:: i_rtop_index is " << i_rtop_index << " but peaks.size() is "
                      << peaks.size() << " for peak " << i << std::endl;
         }
      } else {
         std::cout << "ERROR:: this shouldn't happen! " << std::endl;
         std::cout << "  trapped unset rtop at peak with index " << i_rtop_index  << std::endl;
      } 
   }

   // debugging
   if (r.size() > 10) {
      std::cout << "Top 4 orientation matrices: " << std::endl;
      for (int i=0; i<4; i++)
         std::cout << "Peak number " << i << " at " << r[i].first << "\n"
                   << r[i].second.format() << "\n\n";
   }
   return filter_by_distance_to_higher_peak(r);
}


// The map is all negative, we pass the find the least negative peak
// (i.e. it is highest point)
// 
// Also filter on peak height.  Don't let through peaks lower than
// 0.5*max_peak (now that we use positive peaks).
// 
std::vector<std::pair<float, clipper::RTop_orth> >
coot::util::fffear_search::filter_by_distance_to_higher_peak(const std::vector<std::pair<float, clipper::RTop_orth> > &vr) const {

   // add in a bit of wiggle room that accounts for error in position
   // of an otherwise correct high peak
   // 
   double min_dist    = 2 * (0.8 * min_molecule_radius_);
   double min_dist_sq = min_dist * min_dist;
   std::vector<std::pair<float, clipper::RTop_orth> > r;
   std::cout << "INFO:: Maximum plausible inter-peak distance: " << min_dist << "\n";
   

   if (vr.size() > 0) { 
      double max_peak = vr[0].first;
      for (unsigned int i=0; i<vr.size(); i++) {
         if (vr[i].first > 0.5*max_peak) {
            short int too_close_flag = 0;
            for (unsigned int j=0; j<i; j++) {
               clipper::Coord_orth co(vr[i].second.trn() - vr[j].second.trn()); 
               double d2 = clipper::Coord_orth(co).lengthsq(); 
               if (d2 < min_dist_sq) {
                  std::cout << "Filtered peak " << vr[i].second.trn().format()
                            << " by " << vr[j].second.trn().format() << " dist: "
                            << sqrt(d2) << "\n";
                  too_close_flag = 1;
                  break;
               }
            }
            if (! too_close_flag)
               r.push_back(vr[i]);
         }
      }
   }
   std::cout << "There are " << vr.size() << " raw peaks" << std::endl;
   std::cout << "There are " <<  r.size() << " filtered peaks" << std::endl;
   return r;
}

clipper::RTop_orth
coot::util::fffear_search::mid_point_transformation() const {

   clipper::Mat33<double> mat33(1,0,0,0,1,0,0,0,1);
   return clipper::RTop_orth(mat33, -mid_point_);
}


void
coot::util::fffear_search::post_process_nxmap(float xmap_mean, float xmap_stddev) {

   // Let's now process nxmap, so that then mean and variance of the
   // non-masked regions match the mean and the variance of the
   // nxmap.

   //
   // float density_sum = 0;
   // float run_sum;
   // int density_point_count = 0;
   float density_sum_simple = 0.0;
   float sq_density_sum_simple = 0.0;
   int n_points = 0;

   std::cout << "INFO: target map mean: " << xmap_mean
             << " stddev: " << xmap_stddev << std::endl;
   
   clipper::NXmap<float>::Map_reference_index inx;
   float d;
   for (inx = nxmap.first(); !inx.last(); inx.next()) {
      if (nxmap_mask[inx] > 0.0) {
         d = nxmap[inx];
//          if (clipper::Util::isnan(d)) {
//             std::cout << "Ooops nan in nxmap" << std::endl;
//          } 
         density_sum_simple += d;
         sq_density_sum_simple += d*d;
         n_points++;
      }
   }

   if (n_points) {

//       std::cout << "DEBUG:: mean from " << density_sum_simple << " and "
//                 << n_points << std::endl;
      float current_mean = density_sum_simple/float(n_points);
      float var = sq_density_sum_simple/float(n_points) -
         current_mean*current_mean;
      if (var < 0.00000001)
         var = 0.0000001; 
      float current_std = sqrt(var);

      std::cout << "   nxmap initial mean: " << current_mean
                << " stddev: " << current_std << std::endl;

      float scale_factor = xmap_stddev/current_std;
      float add_factor = xmap_mean - current_mean*scale_factor;

      for (inx = nxmap.first(); !inx.last(); inx.next()) {
         nxmap[inx] *= scale_factor;
         nxmap[inx] += add_factor;
      }

      // checking:
      float d;
      float sum_sq = 0.0;
      float sum = 0.0;
      int npoints = 0; 
      for (inx = nxmap.first(); !inx.last(); inx.next()) {
         if (nxmap_mask[inx] > 0.0) {
            d = nxmap[inx];
            sum += d;
            sum_sq += d*d;
            npoints++;
         }
      }
      float mean = sum/float(npoints);
      var = sum_sq/float(npoints);
      std::cout << "   post-process  mean:  " << mean << " stddev: "
                << sqrt(var) << std::endl;
   }
}
