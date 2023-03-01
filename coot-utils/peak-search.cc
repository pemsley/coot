/* coot-utils/peak-search.cc
 * 
 * Copyright 2005, 2006 by The University of York
 * Copyright 2009 by The University of Oxford
 * Copyright 2015 by Medical Research Council
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

#include <string>
#include <queue>
#include <algorithm>

#include "peak-search.hh"
#include "clipper/core/map_interp.h"

coot::peak_search::peak_search(const clipper::Xmap<float> &xmap) { 
   //                                float n_sigma_in) { 

    max_closeness = 2.0; // don't allow "smaller" peaks that are
                         // within max_closeness of a larger one.

   float sum   = 0.0;
   float sumsq = 0.0;
   float n = 0;
   float v;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next() )  { // iterator index.
      v = xmap[ix];
      sum   += v;
      sumsq += v*v;
      n += 1.0;
   }
   
   float mean = sum/n;
   map_rms = sqrt(sumsq/n - mean*mean);
   // n_sigma = n_sigma_in; 
}


std::vector<clipper::Coord_orth> 
coot::peak_search::get_peaks(const clipper::Xmap<float> &xmap,
                             float n_sigma) {

   std::vector<clipper::Coord_orth> r;
   clipper::Xmap<short int> marked_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = marked_map.first(); !ix.last(); ix.next() )  { // iterator index.
      marked_map[ix] = 0;
   }

   peak_search_0(xmap, &marked_map, n_sigma);

   for (ix = marked_map.first(); !ix.last(); ix.next())  { 
      if (marked_map[ix] == 2) {
         // why do I move peak pos to grid?
         r.push_back(move_grid_to_peak(xmap, ix.coord()));
         // r.push_back(ix.coord().coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()));
      }
   }
   return r;
}

std::vector<clipper::Coord_orth> 
coot::peak_search::get_peaks_for_flooding(const clipper::Xmap<float> &xmap,
                                          float n_sigma) {

   std::vector<clipper::Coord_orth> r;
   clipper::Xmap<short int> marked_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = marked_map.first(); !ix.last(); ix.next() )  { // iterator index.
      marked_map[ix] = 0;
   }

   peak_search_for_flooding(xmap, &marked_map, n_sigma);

   for (ix = marked_map.first(); !ix.last(); ix.next())  { 
      if (marked_map[ix] == 2) {
         r.push_back(move_grid_to_peak(xmap, ix.coord()));
      }
   }
   return r;
}


std::vector<std::pair<clipper::Xmap<float>::Map_reference_index, float> >
coot::peak_search::get_peak_map_indices(const clipper::Xmap<float> &xmap,
                                        float n_sigma) const {

   std::vector<std::pair<clipper::Xmap<float>::Map_reference_index, float> > v;
   clipper::Xmap<short int> marked_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = marked_map.first(); !ix.last(); ix.next() )  { // iterator index.
      marked_map[ix] = 0;
   }

   peak_search_0(xmap, &marked_map, n_sigma);

   for (ix = marked_map.first(); !ix.last(); ix.next())  { 
      if (marked_map[ix] == 2) {
         std::cout << "Peak at " << ix.coord().format() << " " << xmap[ix] << std::endl;
         v.push_back(std::pair<clipper::Xmap<float>::Map_reference_index, float> (ix, xmap[ix]));
      }
   }
   std::sort(v.begin(), v.end(), compare_ps_peaks_mri);

   // debuggin
   if (v.size() > 4) {
      for (int i=0; i<4; i++)
         std::cout << v[i].first.coord().format() << " " << v[i].second << " \n";
   }

   return v;
}


std::vector<std::pair<clipper::Coord_grid, float> >
coot::peak_search::get_peak_grid_points(const clipper::Xmap<float> &xmap,
                                        float n_sigma) const {

   std::vector<std::pair<clipper::Coord_grid, float> > v;
   clipper::Xmap<short int> marked_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = marked_map.first(); !ix.last(); ix.next() )  { // iterator index.
      marked_map[ix] = 0;
   }

   peak_search_0(xmap, &marked_map, n_sigma);

   for (ix = marked_map.first(); !ix.last(); ix.next())  { 
      if (marked_map[ix] == 2) {
//          std::cout << "Peak at " << ix.coord().format() << " " << xmap[ix] << std::endl;
          v.push_back(std::pair<clipper::Coord_grid, float> (ix.coord(), xmap[ix]));
      }
   }

   std::sort(v.begin(), v.end(), compare_ps_peaks_cg);
   return v;
}


std::vector<std::pair<clipper::Coord_grid, float> >
coot::peak_search::get_minima_grid_points(const clipper::Xmap<float> &xmap,
                                          float n_sigma) const {

   std::vector<std::pair<clipper::Coord_grid, float> > v;
   clipper::Xmap<short int> marked_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = marked_map.first(); !ix.last(); ix.next() )  { // iterator index.
      marked_map[ix] = 0;
   }

   peak_search_0_minima(xmap, &marked_map);

   for (ix = marked_map.first(); !ix.last(); ix.next())  { 
      if (marked_map[ix] == 2) {
          v.push_back(std::pair<clipper::Coord_grid, float> (ix.coord(), xmap[ix]));
      }
   }

   std::sort(v.begin(), v.end(), compare_ps_peaks_cg);
   std::reverse(v.begin(), v.end());
   return v;
}


// We introduce now (20090901) peak search filtering, i.e. if there
// are a cluster of points above the n_sigma level, then aggregate the
// peaks of the cluster and return only 1 (the highest of the peaks).
//
// (previously we were not doing clustering and that lead to
// several/many peaks in the map that were in the same peak
// (e.g. ligand cluster points) and that was undesirable - just one
// point in that cluster was neeeded.
//
//
// Note that you can mess with marked_map, but end up with 2 at the
// peaks.  We come here with marked_map having all 0s.
// 
void
coot::peak_search::peak_search_0(const clipper::Xmap<float> &xmap,
                                 clipper::Xmap<short int> *marked_map_p,
                                 float n_sigma) const {

   clipper::Xmap_base::Map_reference_index ix;
   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.25, 1.75); // 3x3x3 cube, not centre
   clipper::Coord_grid c_g;
   int is_peak;
   float v;
   float cut_off = map_rms * n_sigma;
   short int IN_CLUSTER = 3; 

   if (false)
      std::cout << "debug:: peak_search_0():: map rms: " << map_rms << ", peak cut-off: "
                << cut_off << "\n";
   
//    std::cout << "There are " << neighb.size() << " neighbours\n";
//    for (int i=0; i<neighb.size(); i++) {
//       std::cout << i << " " << neighb[i].format() << "\n";
//    }

   for (ix = marked_map_p->first(); !ix.last(); ix.next())  { // iterator index.
      if ((*marked_map_p)[ix] == 0) { 
         v = xmap[ix];
         if (v > cut_off) {
            is_peak = 1;
            for (int i=0; i<neighb.size(); i++) {
               c_g = ix.coord() + neighb[i];
               if (v < xmap.get_data(c_g)) {
                  is_peak = 0;
                  break;
               }
            }
            if (is_peak) {
               // Instead of 
               // (*marked_map_p)[ix] = 2;
               //
               // we now search around the ix point and find those
               // that are above cut_off, but in the same cluster and
               // mark them as IN_CLUSTER
               //
               std::queue<clipper::Coord_grid> q;
               q.push(ix.coord());
               clipper::Coord_grid best_highest_peak_grid = ix.coord(); // updates
               float best_highest_peak = xmap[ix];
               while (q.size()) {
                  clipper::Coord_grid c_g_start = q.front();
                  q.pop();
                  for (int i=0; i<neighb.size(); i++) {
                     c_g = c_g_start + neighb[i];
                     if (marked_map_p->get_data(c_g) == 0) { 
                        marked_map_p->set_data(c_g, IN_CLUSTER);
                        float density_level = xmap.get_data(c_g);
                        if (density_level > cut_off) {
                           q.push(c_g);
                           if (density_level > best_highest_peak) {
                              best_highest_peak = density_level;
                              best_highest_peak_grid = c_g;
                           }
                        }
                     }
                  }
               }
               // mark this as the peak then.
               marked_map_p->set_data(best_highest_peak_grid, 2);
            }
         }
      }
   }
}

// Find troughs.
void
coot::peak_search::peak_search_0_negative(const clipper::Xmap<float> &xmap,
                                          clipper::Xmap<short int> *marked_map_p,
                                          float n_sigma) { 

   clipper::Xmap_base::Map_reference_index ix;
   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.25, 1.75); // 3x3x3 cube, not centre
   clipper::Coord_grid c_g;
   bool is_peak;
   float v;
   float cut_off = -map_rms * n_sigma; 

   // find negative peaks
   // 
   for (ix = marked_map_p->first(); !ix.last(); ix.next())  { // iterator index.
      v = xmap[ix];
      if (v < cut_off) {
         is_peak = true;
         for (int i=0; i<neighb.size(); i++) {
            c_g = ix.coord() + neighb[i];
            if (v > xmap.get_data(c_g)) {
               is_peak = false;
               break;
            }
         }
         if (is_peak) 
            (*marked_map_p)[ix] = 2;
      }
   }
}

void
coot::peak_search::peak_search_for_flooding(const clipper::Xmap<float> &xmap,
                                            clipper::Xmap<short int> *marked_map_p,
                                            float n_sigma) const {

   clipper::Xmap_base::Map_reference_index ix;
   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.25, 1.75); // 3x3x3 cube, not centre
   clipper::Coord_grid c_g;
   bool is_peak;
   float v;
   float cut_off = map_rms * n_sigma;
   short int IN_CLUSTER = 3; 
   
   std::cout << "debug:: peak_search_for_flooding():: map rms: " << map_rms << ", peak cut-off: "
             << cut_off << "\n";
   
   for (ix = marked_map_p->first(); !ix.last(); ix.next())  { // iterator index.
      if ((*marked_map_p)[ix] == 0) { 
         v = xmap[ix];
         if (v > cut_off) {
            is_peak = true;
            for (int i=0; i<neighb.size(); i++) {
               c_g = ix.coord() + neighb[i];
               if (v < xmap.get_data(c_g)) {
                  is_peak = false;
                  break;
               }
            }
            if (is_peak) {
               // mark this as the peak then.
               marked_map_p->set_data(ix.coord(), 2);
            }
         }
      }
   }
}
 

// Find troughs, but no cut-off
void
coot::peak_search::peak_search_0_minima(const clipper::Xmap<float> &xmap,
                                        clipper::Xmap<short int> *marked_map_p) const { 

   clipper::Xmap_base::Map_reference_index ix;
   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.25, 1.75); // 3x3x3 cube, not centre
   clipper::Coord_grid c_g;
   bool is_peak;
   float v;

   // find negative peaks
   // 
   for (ix = marked_map_p->first(); !ix.last(); ix.next())  { // iterator index.
      v = xmap[ix];
      is_peak = true; // initially
      for (int i=0; i<neighb.size(); i++) {
         c_g = ix.coord() + neighb[i];
         if (v > xmap.get_data(c_g)) {
            is_peak = false;
            break;
         }
      }
      if (is_peak) 
         (*marked_map_p)[ix] = 2;
   }
}



clipper::Coord_orth
coot::peak_search::move_grid_to_peak(const clipper::Xmap<float> &xmap,
                        const clipper::Coord_grid &c_g) { 

   clipper::Coord_orth pos = c_g.coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell());
   float dv; // density value
   clipper::Grad_frac<float> grad_frac;
   clipper::Grad_map<float>  grad_map;
   // clipper::Curv_map<float> curv_map;
   float shift_len = 1.0; // Angstroms
   int n_cycle = 0; 
   int n_cycle_max = 500;
   double gradient_scale = 0.01;

   while ((n_cycle < n_cycle_max) && (shift_len > 0.001)) { // Angstroms

      clipper::Coord_frac a_cf = pos.coord_frac(xmap.cell());
      clipper::Coord_map  a_cm = a_cf.coord_map(xmap.grid_sampling());
      clipper::Interp_cubic::interp_grad(xmap, a_cm, dv, grad_map);

      grad_frac = grad_map.grad_frac(xmap.grid_sampling());
      clipper::Grad_orth<float> grad_o = grad_frac.grad_orth(xmap.cell());

      clipper::Coord_orth shift(gradient_scale*grad_o.dx(),
                                gradient_scale*grad_o.dy(),
                                gradient_scale*grad_o.dz());

      shift_len = sqrt(shift.lengthsq());
      pos += shift;
      n_cycle++;
   }
   return pos;
}


void
coot::peak_search::mask_map(clipper::Xmap<float> *xmap_p,
                            const std::vector<clipper::Coord_orth> &ps_peaks) const {

   float radius = 1.2; // A
   for (unsigned int i=0; i<ps_peaks.size(); i++) {
      mask_around_coord(xmap_p, ps_peaks[i], radius);
   }
}

void
coot::peak_search::mask_around_coord(clipper::Xmap<float> *xmap_p,
                                     const clipper::Coord_orth &co,
                                     float atom_radius) const {
   
   clipper::Coord_frac cf = co.coord_frac(xmap_p->cell());
   float masked_map_val = 0.0;

   clipper::Coord_frac box0(
                            cf.u() - atom_radius/xmap_p->cell().descr().a(),
                            cf.v() - atom_radius/xmap_p->cell().descr().b(),
                            cf.w() - atom_radius/xmap_p->cell().descr().c());

   clipper::Coord_frac box1(
                            cf.u() + atom_radius/xmap_p->cell().descr().a(),
                            cf.v() + atom_radius/xmap_p->cell().descr().b(),
                            cf.w() + atom_radius/xmap_p->cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(xmap_p->grid_sampling()),
                          box1.coord_grid(xmap_p->grid_sampling()));

   float atom_radius_sq = atom_radius * atom_radius;
   // int nhit = 0;
   // int nmiss = 0;

   clipper::Xmap_base::Map_reference_coord ix( *xmap_p, grid.min() ), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) { 
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) { 
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            if ( (iw.coord().coord_frac(xmap_p->grid_sampling()).coord_orth(xmap_p->cell()) - co).lengthsq() < atom_radius_sq) { 
               (*xmap_p)[iw] = masked_map_val;
            } else {
            }
         }
      }
   }
}


void
coot::peak_search::add_peak_vectors(std::vector<clipper::Coord_orth> *in,
                                    const std::vector<clipper::Coord_orth> &extras) const {

   for(unsigned int i=0; i<extras.size(); i++)
      in->push_back(extras[i]);
}



// filter_peaks_by_closeness() called at the end (and the distance is
// set by set_max_closeness()).
// 
std::vector<std::pair<clipper::Coord_orth, float> >
coot::peak_search::get_peaks(const clipper::Xmap<float> &xmap,
                             float n_sigma,
                             int do_positive_levels_flag,
                             int also_negative_levels_flag) {

   std::vector<std::pair<clipper::Coord_orth, float> > r;

   clipper::Xmap<short int> marked_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = marked_map.first(); !ix.last(); ix.next()) 
      marked_map[ix] = 0;

   peak_search_0(xmap, &marked_map, n_sigma);

   if (do_positive_levels_flag) 
      for (ix = marked_map.first(); !ix.last(); ix.next()) 
         if (marked_map[ix] == 2)
            r.push_back(std::pair<clipper::Coord_orth, float> (move_grid_to_peak(xmap, ix.coord()),
                                                               xmap[ix]));


   // negative levels too?
   if (also_negative_levels_flag) {

      // reset the marked map
      for (ix = marked_map.first(); !ix.last(); ix.next()) 
         marked_map[ix] = 0;
      peak_search_0_negative(xmap, &marked_map, n_sigma);
      for (ix = marked_map.first(); !ix.last(); ix.next())  { 
         if (marked_map[ix] == 2) {
            r.push_back(std::pair<clipper::Coord_orth, float> (ix.coord_orth(), xmap[ix]));
         }
      }
   }
   std::sort(r.begin(), r.end(), compare_ps_peaks);

   std::vector<std::pair<clipper::Coord_orth, float> > rf = filter_peaks_by_closeness(r);
   
   return rf;
}

#include "coot-coord-utils.hh"

std::vector<std::pair<clipper::Coord_orth, float> >
coot::peak_search::get_peaks(const clipper::Xmap<float> &xmap,
                             mmdb::Manager *mol, 
                             float n_sigma,
                             int do_positive_levels_flag,
                             int also_negative_levels_flag,
                             int only_around_protein_flag) {

   auto sample_all_atoms = [mol] () {
                              int imod = 1;
                              mmdb::Model *model_p = mol->GetModel(imod);
                              std::vector<clipper::Coord_orth> coords;
                              if (model_p) {
                                 int n_chains = model_p->GetNumberOfChains();
                                 for (int ichain=0; ichain<n_chains; ichain++) {
                                    mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                    int n_res = chain_p->GetNumberOfResidues();
                                    for (int ires=0; ires<n_res; ires++) {
                                       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                       if (residue_p) {
                                          int n_atoms = residue_p->GetNumberOfAtoms();
                                          for (int iat=0; iat<n_atoms; iat++) {
                                             mmdb::Atom *at = residue_p->GetAtom(iat);
                                             if (! at->isTer()) {
                                                coords.push_back(co(at));
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                              return coords;
                           };

   std::vector<std::pair<clipper::Coord_orth, float> > peaks =
      get_peaks(xmap, n_sigma, do_positive_levels_flag, also_negative_levels_flag);

   std::vector<std::pair<clipper::Coord_orth, float> > r;

   std::vector<clipper::Coord_orth> sampled_protein_coords = make_sample_protein_coords(mol);

   if (only_around_protein_flag) {

      sampled_protein_coords = sample_all_atoms();
      // std::cout << "sampled_protein_coords size is " << sampled_protein_coords.size() << std::endl;

      double max_dist = 4.0; // maybe a bit too much, actually!
      double max_dist_sqrd = max_dist * max_dist;
      unsigned int n_outside = 0;
      for (const auto &peak : peaks) {
         bool within_dist = false;
         for (const auto &atom_pos : sampled_protein_coords) {
            double dd = (atom_pos-peak.first).lengthsq();
            if (dd < max_dist_sqrd) {
               within_dist = true;
               break;
            }
         }
         if (within_dist) {
            r.push_back(peak);
         } else {
            n_outside++;
         }
      }
      // std::cout << "........... n_outside " << n_outside << std::endl;

   } else {

      // The protein could be somewhere away from the origin.  In that
      // case, (indeed, the general case) we should find what the
      // translation operation is to move the centre of the protein
      // coords as closs to the origin it can get (by translation only).
      //
      // We will then move sampled protein coords by this translation
      // vector - so that all tests for distance to new waters are to a
      // protein that is centred round the origin.
      //
      // When we come to actually define the position of the water
      // though, it must apply the reverse translation operator to bring
      // the water point back to where the protein actually is.

      // The unit cell translations to move the centre of the protein as
      // closs to the origin as possible.
      // 
      std::vector<int> iprotein_trans =
         find_protein_to_origin_translations(sampled_protein_coords, xmap);

      //
      if (false)
         std::cout << "DEBUG:: iprotein_trans: "
                   << iprotein_trans[0] << " "
                   << iprotein_trans[1] << " "
                   << iprotein_trans[2] << std::endl;

      // move sampled protein by this translation then:
      //
      clipper::Mat33<double> m_identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
      for (unsigned int i=0; i<sampled_protein_coords.size(); i++) { 
         clipper::Coord_frac cell_shift(iprotein_trans[0],
                                        iprotein_trans[1],
                                        iprotein_trans[2]);
         clipper::RTop_frac rtf(m_identity, cell_shift);
         clipper::RTop_orth orthop = rtf.rtop_orth(xmap.cell());
         sampled_protein_coords[i] = sampled_protein_coords[i].transform(orthop);
      }

      for (unsigned int i=0; i<peaks.size(); i++) {
         clipper::Coord_orth pt = move_point_close_to_protein(peaks[i].first,
                                                              sampled_protein_coords,
                                                              iprotein_trans,
                                                              xmap);
         r.push_back(std::pair<clipper::Coord_orth, float> (pt, peaks[i].second));
      }
   }

   return r;
}

std::vector<std::pair<clipper::Coord_orth, float> >
coot::peak_search::filter_peaks_by_closeness(const std::vector<std::pair<clipper::Coord_orth, float> > &v) const {

   
   std::vector<std::pair<clipper::Coord_orth, float> > nv;

   if (max_closeness > 0) { 
      for (unsigned int ipeak=0; ipeak<v.size(); ipeak++) {
         bool found_a_close_one = 0;
         for (unsigned int jpeak=0; jpeak<nv.size(); jpeak++) {
            double d = clipper::Coord_orth::length(v[ipeak].first, nv[jpeak].first);
            if (d < max_closeness) { 
               found_a_close_one = 1;
               break;
            }
         }
         if (! found_a_close_one)
            nv.push_back(v[ipeak]);
      }
   } else {
      nv = v;
   } 
   // std::cout << "INFO:: peak filtering: npeaks: in: " << v.size() << " out: " << nv.size() << std::endl;
   return nv;
} 


std::vector<clipper::Coord_orth>
coot::peak_search::make_sample_protein_coords(mmdb::Manager *mol) const {

   return make_sample_protein_coords(mol, 5);
}

std::vector<clipper::Coord_orth>
coot::peak_search::make_sample_protein_coords(mmdb::Manager *mol, int every_n) const {

   std::vector<clipper::Coord_orth> r;
   int atom_count = every_n;

   int imod = 1;
      
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::PResidue residue_p;
      mmdb::Atom *at;
      for (int ires=0; ires<nres; ires++) { 
         residue_p = chain_p->GetResidue(ires);
         int n_atoms = residue_p->GetNumberOfAtoms();
         
         for (int iat=0; iat<n_atoms; iat++) {
            if (atom_count == every_n) { 
               at = residue_p->GetAtom(iat);
               r.push_back(clipper::Coord_orth(at->x, at->y, at->z));
               atom_count = 0; 
            }
            atom_count++;
         }
      }
   }
   return r;
}

const std::vector<int>
coot::peak_search::find_protein_to_origin_translations(const std::vector<clipper::Coord_orth> &sampled_protein_coords, const clipper::Xmap<float> &xmap) const {
   
   std::vector<int> r(3, 0);
   clipper::Coord_orth origin(0,0,0);

   // get the protein centre:
   clipper::Coord_orth running(0,0,0);
   for (unsigned int i=0; i<sampled_protein_coords.size(); i++)
      running += sampled_protein_coords[i];

   if (sampled_protein_coords.size() > 0) {
      double frac = 1.0/double(sampled_protein_coords.size());
      clipper::Coord_orth centre(running.x()*frac,
                                 running.y()*frac,
                                 running.z()*frac);

      short int was_shifted = 1;
      double min_distance_to_origin = 99999999999.9;
      clipper::Mat33<double> m_identity(1, 0, 0,
                                        0, 1, 0,
                                        0, 0, 1);
      int cur_x = 0;
      int cur_y = 0;
      int cur_z = 0;
      while (was_shifted) {
         was_shifted = 0;
         // std::cout << "shifted loop " << min_distance_to_origin << std::endl;
         int start_cur_x = cur_x;
         int start_cur_y = cur_y;
         int start_cur_z = cur_z;
         for (int x_shift = -1; x_shift<2; x_shift++) { 
            for (int y_shift = -1; y_shift<2; y_shift++) { 
               for (int z_shift = -1; z_shift<2; z_shift++) {
                  clipper::Coord_frac cell_shift(start_cur_x + x_shift,
                                                 start_cur_y + y_shift,
                                                 start_cur_z + z_shift); 
                  clipper::RTop_frac rtf(m_identity, cell_shift);
                  clipper::RTop_orth orthop = rtf.rtop_orth(xmap.cell());
                  clipper::Coord_orth t_point = centre.transform(orthop);
                  double d = clipper::Coord_orth::length(t_point, origin);
                  // the 0.001 added to ameliorate numerical
                  // instabilility (bizarrely enough).
                  if (d < (min_distance_to_origin - 0.001)) {
//                      std::cout << "hit " << d << " < " << min_distance_to_origin << std::endl;
                     min_distance_to_origin = d;
                     cur_x = start_cur_x + x_shift;
                     cur_y = start_cur_y + y_shift;
                     cur_z = start_cur_z + z_shift;
//                      std::cout << "hit x " << cur_x << " " << start_cur_x << " " << x_shift << std::endl;
//                      std::cout << "hit y " << cur_y << " " << start_cur_y << " " << y_shift << std::endl;
//                      std::cout << "hit z " << cur_z << " " << start_cur_z << " " << z_shift << std::endl;
                     was_shifted = 1;
                     r[0] = cur_x;
                     r[1] = cur_y;
                     r[2] = cur_z;
                  }
               }
            }
         }
      }
      // debug
//       {
//          clipper::Coord_frac cell_shift(r[0], r[1], r[2]);
//          clipper::RTop_frac rtf(m_identity, cell_shift);
//          clipper::RTop_orth orthop = rtf.rtop_orth(xmap.cell());
//          std::cout << "DEBUG:: applying r ("
//                    << r[0] << " " << r[1] << " " << r[2] << ") "
//                    << "to protein centre at "
//                    << centre.format() << " moves to "
//                    << centre.transform(orthop).format() << std::endl;
//          std::cout << "DEBUG:: rtop frac: \n" << rtf.format() << std::endl;
//          std::cout << "DEBUG:: rtop orth: \n" << orthop.format() << std::endl;
//       }
   }
   return r;

}


clipper::Coord_orth
coot::peak_search::move_point_close_to_protein(const clipper::Coord_orth &pt,
                                               const std::vector<clipper::Coord_orth> &sampled_protein_coords,
                                               const std::vector<int> &itrans,
                                               const clipper::Xmap<float> &xmap) const {

   // Recall that sampled_protein_coords have been moved close to the
   // origin, so when we make r, we need to move it back to where the
   // protein really is.  The itrans is how we move the protein close
   // to the origin, so we need to apply the reverse of those shifts
   // when we make r.
   
   int shift_range = 2; // was 1, changed 20091101, Kevin Keating bug (P212121).
   clipper::Coord_orth r = pt; 
   int n = sampled_protein_coords.size();
   if (n > 0) {
      int nsyms = xmap.spacegroup().num_symops();
      clipper::Coord_frac cell_shift;
      double min_dist = 9999999999.9;
      for (int isym=0; isym<nsyms; isym++) {
         for (int x_shift = -shift_range; x_shift<=shift_range; x_shift++) { 
            for (int y_shift = -shift_range; y_shift<=shift_range; y_shift++) { 
               for (int z_shift = -shift_range; z_shift<=shift_range; z_shift++) {
                  cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift);
                  clipper::RTop_orth orthop =
                     clipper::RTop_frac(xmap.spacegroup().symop(isym).rot(),
                                        xmap.spacegroup().symop(isym).trn() +
                                        cell_shift).rtop_orth(xmap.cell());
                  clipper::Coord_orth t_point = pt.transform(orthop);
                  double t_dist = min_dist_to_protein(t_point, sampled_protein_coords);
                  if (t_dist < min_dist) {
                     min_dist = t_dist;
                     r = t_point;
                  }
               }
            }
         }
      }
   }
   // so now r is in the "close to origin/in unit cell" system.  Let's
   // move it to where the protein actually is.
   clipper::Mat33<double> m_identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
   clipper::Coord_frac cell_shift(-itrans[0], -itrans[1], -itrans[2]);
   clipper::RTop_frac rtf(m_identity, cell_shift);
   clipper::RTop_orth orthop = rtf.rtop_orth(xmap.cell());
   r = r.transform(orthop);
   return r;
}
 

double
coot::peak_search::min_dist_to_protein(const clipper::Coord_orth &point,
                                       const std::vector<clipper::Coord_orth> &sampled_protein_coords) const {

   double dist = 9999999.9;
   double this_dist;
   int n = sampled_protein_coords.size();
   if (n > 0) { 
      for (int i=0; i<n; i++) {
         this_dist = clipper::Coord_orth::length(point, sampled_protein_coords[i]);
         if (this_dist < dist)
            dist = this_dist;
      }
   } else {
      dist = 0.0;
   } 
   return dist;
}


// static
bool
coot::peak_search::compare_ps_peaks(const std::pair<clipper::Coord_orth, float> &a,
                                    const std::pair<clipper::Coord_orth, float> &b) {


   float a1 = fabs(a.second);
   float b1 = fabs(b.second);

   return a1 > b1;

}

// static
bool
coot::peak_search::compare_ps_peaks_mri(const std::pair<clipper::Xmap<float>::Map_reference_index, float> &a,
                                        const std::pair<clipper::Xmap<float>::Map_reference_index, float> &b) {

   float a1 = a.second;
   float b1 = b.second;
   return a1 > b1;
}
                                    

// static
bool
coot::peak_search::compare_ps_peaks_cg(const std::pair<clipper::Coord_grid, float> &a,
                                       const std::pair<clipper::Coord_grid, float> &b) {

   float a1 = a.second;
   float b1 = b.second;
   return a1 > b1;
}
                                    

