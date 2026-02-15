/* coot-utils/coot-map-utils.cc
 *
 * Copyright 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009, 2010 by The University of Oxford
 * Copyright 2013, 2014, 2015 by Medical Research Council
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


#include <algorithm> // for sorting.
#include <queue>
#include <fstream>
#include <thread>
#include <iomanip>
#include <filesystem>

#include <gsl/gsl_sf_bessel.h>

#include "clipper/core/coords.h"
#include "clipper/core/map_interp.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/contrib/skeleton.h"
#include <clipper/contrib/edcalc.h>

#include "coot-coord-utils.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "utils/coot-utils.hh"
#include "coot-map-utils.hh"
#include "geometry/main-chain.hh"
#include "exp-fit.hh"

#include "utils/logging.hh"
extern logging logger;

bool
coot::util::map_fill_from_mtz(clipper::Xmap<float> *xmap,
                              std::string mtz_file_name,
                              std::string f_col,
                              std::string phi_col,
                              std::string weight_col,
                              short int use_weights,
                              float sampling_rate) {

   // sampling_rate is optional arg with default value 1.5
   return coot::util::map_fill_from_mtz(xmap, mtz_file_name, f_col, phi_col, weight_col,
                                        use_weights, 0, 0, sampling_rate);
}

bool
coot::util::map_fill_from_mtz(clipper::Xmap<float> *xmap,
                              std::string mtz_file_name,
                              std::string f_col,
                              std::string phi_col,
                              std::string weight_col,
                              short int use_weights,
                              float reso_limit_high,
                              short int use_reso_limit_high,
                              float sampling_rate) {

   // sampling_rate is optional arg with default value 1.5

   auto path_to_file = [] (const std::string &p_col_in) {

                         std::filesystem::path p(p_col_in);
                         std::filesystem::path p_col_path = p.filename();
                         std::string p_col = p_col_path.string();
                         return p_col;
                       };


   // I am not sure that stripping the dataset info is a good thing.
   //
   auto make_import_datanames = [path_to_file] (const std::string &f_col_in,
                                                const std::string &phi_col_in,
                                                const std::string &weight_col_in,
                                                bool use_weights) {

                                  std::pair<std::string, std::string> p("", ""); // return this

                                  std::string      f_col = path_to_file(f_col_in);
                                  std::string    phi_col = path_to_file(phi_col_in);
                                  std::string weight_col = path_to_file(weight_col_in);

                                  std::string no_xtal_dataset_prefix= "/*/*/";
                                  if (use_weights) {
                                    p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " +      f_col + "]";
                                    p.second = no_xtal_dataset_prefix + "[" + phi_col + " " + weight_col + "]";
                                  } else {
                                    p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " + phi_col + "]";
                                  }
                                  return p;
                                };


   if (!file_exists(mtz_file_name))
      return false;

   clipper::HKL_info myhkl;
   clipper::MTZdataset myset;
   clipper::MTZcrystal myxtl;

   // std::cout << "reading mtz file..." << std::endl;
   clipper::CCP4MTZfile mtzin;
   mtzin.open_read( mtz_file_name );       // open new file
   mtzin.import_hkl_info( myhkl );         // read sg, cell, reso, hkls
   clipper::HKL_data< clipper::datatypes::F_sigF<float> >   f_sigf_data(myhkl, myxtl);
   clipper::HKL_data< clipper::datatypes::Phi_fom<float> > phi_fom_data(myhkl, myxtl);
   clipper::HKL_data< clipper::datatypes::F_phi<float> >       fphidata(myhkl, myxtl);

   // If use weights, use both strings, else just use the first
   std::pair<std::string, std::string> datanames = make_import_datanames(f_col, phi_col, weight_col, use_weights);

   if (false) {
      std::cout << ":::::::::::::::::::::: datanames:" << std::endl;
      std::cout << "                      " << datanames.first << std::endl;
      std::cout << "                      " << datanames.second << std::endl;
   }

   if (use_weights) {
      std::string dataname = datanames.first;
      mtzin.import_hkl_data(  f_sigf_data, myset, myxtl, dataname );
      dataname = datanames.second;
      mtzin.import_hkl_data( phi_fom_data, myset, myxtl, dataname );
      mtzin.close_read();
      std::cout << "in map_fill_from_mtz(): We should use the weights: " << weight_col << std::endl;

      fphidata.compute(f_sigf_data, phi_fom_data,
                       clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

   } else {

      clipper::String dataname = datanames.first;
      mtzin.import_hkl_data(fphidata, myset, myxtl, dataname);
      mtzin.close_read();

   }
   // std::cout << "Number of reflections: " << myhkl.num_reflections() << "\n";
   // std::cout << "finding ASU unique map points..." << std::endl;

   clipper::Resolution fft_reso;
   if (use_reso_limit_high) {
      clipper::Resolution user_resolution(reso_limit_high);
      fft_reso = user_resolution;
      coot::util::filter_by_resolution(&fphidata, 99999.0, reso_limit_high);
   } else {
      fft_reso = clipper::Resolution(1.0/sqrt(fphidata.invresolsq_range().max()));
   }

   // std::cout << "FFT Reso..." << fft_reso.invresolsq_limit() << "\n";
   // std::cout << "Sampling rate..." << sampling_rate << "\n";
   logger.log(log_t::INFO, "FFT Resolution:", fft_reso.invresolsq_limit());
   logger.log(log_t::INFO, "Map Sampling Rate:", sampling_rate);
   clipper::Grid_sampling gs(myhkl.spacegroup(), myhkl.cell(), fft_reso, sampling_rate);
   // std::cout << "Grid..." << gs.format() << "\n";
   // std::cout << "Cell..." << myhkl.cell().format() << "\n";
   // std::cout << "Spacegroup..." << myhkl.spacegroup().symbol_hm() << "\n";
   logger.log(log_t::INFO, "Grid:", gs.format());
   logger.log(log_t::INFO, "Cell:", myhkl.cell().format());
   logger.log(log_t::INFO, "Spacegroup:", myhkl.spacegroup().symbol_hm());
   if (gs.nu() == 0) { std::cout << "Bad Grid\n"; return false; }
   if (gs.nv() == 0) { std::cout << "Bad Grid\n"; return false; }
   if (gs.nw() == 0) { std::cout << "Bad Grid\n"; return false; }
   xmap->init( myhkl.spacegroup(), myhkl.cell(), gs);
   // std::cout << "doing fft..." << std::endl;
   xmap->fft_from( fphidata );                  // generate map
   // std::cout << "done fft..." << std::endl;
   return true;
}

// Return a map that is a copy of the given map with interpolation,
// with grid spacing at most 0.5A (by generated by integer scaling
// factor of the input map)
//
clipper::Xmap<float>
coot::util::reinterp_map_fine_gridding(const clipper::Xmap<float> &xmap) {

   float cell_length[3];
   cell_length[0] = xmap.cell().descr().a();
   cell_length[1] = xmap.cell().descr().b();
   cell_length[2] = xmap.cell().descr().c();

   int sampling[3];
   sampling[0] = xmap.grid_sampling().nu();
   sampling[1] = xmap.grid_sampling().nv();
   sampling[2] = xmap.grid_sampling().nw();

   float rate[3];
   for (int i=0; i<3; i++)
      rate[i] = float(sampling[i])/(cell_length[i] + 0.0001);
   float smallest_rate = 100000;
   for (int i=0; i<3; i++)
      if (rate[i] < smallest_rate)
         smallest_rate = rate[i];
   // so now we have smallest_rate;

   // what do we need to multiply smallest_rate by to get a number
   // better than 2?
   int multiplier = 1;
   int count = 1000;
   for (int icount = 0; icount<count; icount++) {
      if (smallest_rate * float(multiplier) > 2.0)
         break;
      multiplier += 1;
   }
   // so now we have multiplier

   if (multiplier == 1) {
      return xmap; // the input map
   } else {

      clipper::Grid_sampling gs_old = xmap.grid_sampling();
      clipper::Grid_sampling gs(sampling[0]*multiplier,
                                sampling[1]*multiplier,
                                sampling[2]*multiplier);

      clipper::Xmap<float> xmap_new(xmap.spacegroup(), xmap.cell(), gs);

      float dv;
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = xmap_new.first(); !ix.last(); ix.next() )  { // iterator index.
         clipper::Coord_grid cg = ix.coord();
         clipper::Coord_frac cf = cg.coord_frac(gs);
         clipper::Coord_grid cg_old = cf.coord_grid(gs_old);
         dv = xmap.interp<clipper::Interp_cubic>(cf);
         xmap_new[ix] = dv;
      }
      return xmap_new;
   }
}

clipper::Xmap<float>
   coot::util::reinterp_map(const clipper::Xmap<float> &xmap_in, float sampling_multiplier) {

   clipper::Grid_sampling gs_old = xmap_in.grid_sampling();
   clipper::Grid_sampling gs(xmap_in.grid_sampling().nu()*sampling_multiplier,
                             xmap_in.grid_sampling().nv()*sampling_multiplier,
                             xmap_in.grid_sampling().nw()*sampling_multiplier);

   clipper::Xmap<float> xmap_new(xmap_in.spacegroup(), xmap_in.cell(), gs);

   float dv;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_new.first(); !ix.last(); ix.next() )  { // iterator index.
      clipper::Coord_grid cg = ix.coord();
      clipper::Coord_frac cf = cg.coord_frac(gs);
      clipper::Coord_grid cg_old = cf.coord_grid(gs_old);
      dv = xmap_in.interp<clipper::Interp_cubic>(cf);
      xmap_new[ix] = dv;
   }
   return xmap_new;
}

// Note density_at_point in molecule-class-info looks as if its
// returning grid values, not interpollated values.  Is that bad?
//
float
coot::util::density_at_point(const clipper::Xmap<float> &xmap,
                             const clipper::Coord_orth &pos) {

   float dv;
   clipper::Coord_frac a_cf = pos.coord_frac(xmap.cell());
   clipper::Coord_map  a_cm = a_cf.coord_map(xmap.grid_sampling());
   clipper::Interp_cubic::interp(xmap, a_cm, dv);
   return dv;
}

float
coot::util::density_at_point_by_cubic_interp(const clipper::NXmap<float> &nxmap,
                                             const clipper::Coord_map &a_cm) {

   float dv;
   clipper::Interp_cubic::interp(nxmap, a_cm, dv);
   return dv;

}


//
float
coot::util::density_at_point_by_linear_interpolation(const clipper::Xmap<float> &xmap,
                                                     const clipper::Coord_orth &pos) {

   float dv;
   clipper::Coord_frac a_cf = pos.coord_frac(xmap.cell());
   clipper::Coord_map  a_cm = a_cf.coord_map(xmap.grid_sampling());
   clipper::Interp_linear::interp(xmap, a_cm, dv);
   return dv;
}

// nearest grid point - faster yet
float
coot::util::density_at_point_by_nearest_grid(const clipper::Xmap<float> &xmap,
                                             const clipper::Coord_orth &co) {

   float dv;
   clipper::Coord_frac a_cf = co.coord_frac(xmap.cell());
   clipper::Coord_map  a_cm = a_cf.coord_map(xmap.grid_sampling());
   clipper::Interp_nearest::interp(xmap, a_cm, dv);
   return dv;
}

// NXmap nearest-grid sampling
float
coot::util::density_at_point_by_nearest_grid(const clipper::NXmap<float> &nxmap,
                                             const clipper::Coord_orth &co) {
   float dv = 0.0f;
   clipper::Coord_map a_cm = nxmap.coord_map(co);
   clipper::Interp_nearest::interp(nxmap, a_cm, dv);
   return dv;
}

// NXmap linear interpolation sampling
float
coot::util::density_at_point_by_linear_interp(const clipper::NXmap<float> &nxmap,
                                              const clipper::Coord_orth &co) {
   float dv = 0.0f;
   clipper::Coord_map a_cm = nxmap.coord_map(co);
   clipper::Interp_linear::interp(nxmap, a_cm, dv);
   return dv;
}


void
coot::util::filter_by_resolution(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata,
                                 const float &reso_low,
                                 const float &reso_high) {

   float inv_low  = 1.0/(reso_low*reso_low);
   float inv_high = 1.0/(reso_high*reso_high);
   int n_data = 0;
   int n_reset = 0;

   for (clipper::HKL_info::HKL_reference_index hri = fphidata->first(); !hri.last(); hri.next()) {
      n_data++;

      if ( hri.invresolsq() > inv_low &&
           hri.invresolsq() < inv_high) {
      } else {
         (*fphidata)[hri].f() = 0.0;
         n_reset++;
      }
   }
}


float
coot::util::density_at_map_point(const clipper::Xmap<float> &xmap,
                                 const clipper::Coord_map   &cm) {

   float dv;
   clipper::Interp_cubic::interp(xmap, cm, dv);
   return dv;
}

clipper::Grad_orth<double>
coot::util::gradient_at_point(const clipper::Xmap<float> &xmap_in,
                             const clipper::Coord_orth &co) {

   clipper::Grad_map<double> grad;
   double dv;

   clipper::Coord_frac af = co.coord_frac(xmap_in.cell());
   clipper::Coord_map  am = af.coord_map(xmap_in.grid_sampling());
   // clipper::Interp_linear::Interp_grad(*xmap_p, am, dv, grad) is not a thing
   clipper::Interp_cubic::interp_grad(xmap_in, am, dv, grad);
   clipper::Grad_frac<double> grad_frac = grad.grad_frac(xmap_in.grid_sampling());
   return grad_frac.grad_orth(xmap_in.cell());
}


coot::util::density_stats_info_t
coot::util::density_around_point(const clipper::Coord_orth &point,
                                 const clipper::Xmap<float> &xmap,
                                 float d) {

   coot::util::density_stats_info_t s;
   // what are the points?
   // +/- d along the x, y and z axes
   // and 2x4 samples at 45 degrees inclination and 45 degree longitudinal rotation.
   //
   std::vector<clipper::Coord_orth> sample_points(14);

   sample_points[0] = clipper::Coord_orth( 0.0,  0.0,  1.0);
   sample_points[1] = clipper::Coord_orth( 0.0,  0.0, -1.0);
   sample_points[2] = clipper::Coord_orth( 0.0,  1.0,  0.0);
   sample_points[3] = clipper::Coord_orth( 0.0, -1.0,  0.0);
   sample_points[4] = clipper::Coord_orth(-1.0,  0.0,  0.0);
   sample_points[5] = clipper::Coord_orth( 1.0,  0.0,  0.0);

   sample_points[6] = clipper::Coord_orth( 0.5,  0.5,  0.7071);
   sample_points[7] = clipper::Coord_orth(-0.5,  0.5,  0.7071);
   sample_points[8] = clipper::Coord_orth(-0.5, -0.5,  0.7071);
   sample_points[9] = clipper::Coord_orth( 0.5, -0.5,  0.7071);

   sample_points[10] = clipper::Coord_orth( 0.5,  0.5, -0.7071);
   sample_points[11] = clipper::Coord_orth(-0.5,  0.5, -0.7071);
   sample_points[12] = clipper::Coord_orth(-0.5, -0.5, -0.7071);
   sample_points[13] = clipper::Coord_orth( 0.5, -0.5, -0.7071);

   float dv;
   for (float scale = 0.2; (scale-0.0001)<=1.0; scale += 0.4) {
      for (int i=0; i<14; i++) {
         dv = density_at_point(xmap, point + d*scale*sample_points[i]);
         s.add(dv, scale); // scale is the weight, which multiplies
                           // additions internally.
      }
   }

   return s;
}


float
coot::util::map_score(mmdb::PPAtom atom_selection,
                      int n_selected_atoms,
                      const clipper::Xmap<float> &xmap,
                      short int with_atomic_weighting) {

   // Thanks Ezra.

   float f = 0.0;
   float f1;

   for (int i=0; i<n_selected_atoms; i++) {
      mmdb::Atom *at = atom_selection[i];
      if (! at->isTer()) {
         f1 = density_at_point(xmap, clipper::Coord_orth(atom_selection[i]->x,
                                                         atom_selection[i]->y,
                                                         atom_selection[i]->z));
         f1 *= atom_selection[i]->occupancy;
         f += f1;
         // std::cout << "debug:: map_score() adding " << atom_spec_t(at) << " f1 " << f1 << std::endl;
      }
   }
   return f;
}

float
coot::util::map_score(std::vector<mmdb::Atom *> atoms,
                      const clipper::Xmap<float> &xmap) {

   float f = 0.0;
   for (unsigned int i=0; i<atoms.size(); i++) {
      if (atoms[i]) {
         float f1 = density_at_point(xmap, co(atoms[i]));
         f1 *= atoms[i]->occupancy;
         f += f1;
      }
   }
   return f;
}


float coot::util::map_score_atom(mmdb::Atom *atom,
                                 const clipper::Xmap<float> &xmap) {

   float f = 0;
   if (atom) {
      f = density_at_point(xmap, clipper::Coord_orth(atom->x, atom->y, atom->z));
   }
   return f;
}

float
coot::util::map_score_by_residue_specs(mmdb::Manager *mol,
                                       const std::vector<residue_spec_t> &res_specs,
                                       const clipper::Xmap<float> &xmap,
                                       bool main_chain_only_flag) {

   float f = 0;
   for (std::size_t i=0; i<res_specs.size(); i++) {
      mmdb::Residue *residue_p = get_residue(res_specs[i], mol);
      if (residue_p) {
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! main_chain_only_flag || is_main_chain_or_cb_p(at)) {
               if (false) // debug
                  std::cout << "map_score_by_residue_specs "
                            << atom_spec_t(at) << " " << map_score_atom(at, xmap) << std::endl;
               f += map_score_atom(at, xmap);
            }
         }
      }
   }
   return f;
}




// A/grid
//
float
coot::util::max_gridding(const clipper::Xmap<float> &xmap) {

   float a_gridding = xmap.cell().a()/xmap.grid_sampling().nu();
   float b_gridding = xmap.cell().b()/xmap.grid_sampling().nv();
   float c_gridding = xmap.cell().c()/xmap.grid_sampling().nw();

   float gridding_max = 0;

   if ( a_gridding > gridding_max)
      gridding_max = a_gridding;
   if ( b_gridding > gridding_max)
      gridding_max = b_gridding;
   if ( c_gridding > gridding_max)
      gridding_max = c_gridding;

   return gridding_max;
}

clipper::RTop_orth
coot::util::make_rtop_orth_from(mmdb::mat44 mat) {

   clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
                                      mat[1][0], mat[1][1], mat[1][2],
                                      mat[2][0], mat[2][1], mat[2][2]);
   clipper::Coord_orth  cco(mat[0][3], mat[1][3], mat[2][3]);
   clipper::RTop_orth rtop(clipper_mat, cco);

   return rtop;
}

clipper::Xmap<float>
coot::util::sharpen_map(const clipper::Xmap<float> &xmap_in, float sharpen_factor) {

   // Does this function work?

   clipper::HKL_info myhkl;
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);

   xmap_in.fft_to(fphis);

   clipper::HKL_info::HKL_reference_index hri;
   for (hri = fphis.first(); !hri.last(); hri.next()) {
      float irs = hri.invresolsq();
      // std::cout << hri.format() << " has reso " << irs << std::endl;
      float fac = exp(-sharpen_factor * irs * 0.25);
      fphis[hri].f() *= fac;
   }

   clipper::Xmap<float> r;
   r.fft_from(fphis);
   return r;
}

clipper::Xmap<float>
coot::util::sharpen_blur_map(const clipper::Xmap<float> &xmap_in, float b_factor) {

   float mg = coot::util::max_gridding(xmap_in);
   clipper::Resolution reso(2.0 * mg);
   clipper::HKL_info myhkl(xmap_in.spacegroup(), xmap_in.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
   clipper::Xmap<float> xmap_out(xmap_in.spacegroup(), xmap_in.cell(), xmap_in.grid_sampling());
   xmap_in.fft_to(fphis);
   clipper::HKL_info::HKL_reference_index hri;

   /*
   // using a map to cache the scale factors is 10 times slower
   std::map<float, float> reso_map;
   std::map<float, float>::const_iterator it;
   //...and inside loop:
   it = reso_map.find(irs);
   if (it != reso_map.end()) {
      fphis[hri].f() *= it->second;
   } else {
      float esf = exp(-b_factor * irs * 0.25);
      reso_map[irs] = esf;
      fphis[hri].f() *= esf;
   }
   */

   int count = 0;
   auto tp_1 = std::chrono::high_resolution_clock::now();
   for (hri = fphis.first(); !hri.last(); hri.next()) {
      float f = fphis[hri].f();
      if (! clipper::Util::is_nan(f)) {
         float irs =  hri.invresolsq();
         fphis[hri].f() *= exp(-b_factor * irs * 0.25);
         count++;
      }
   }
   auto tp_2 = std::chrono::high_resolution_clock::now();
   xmap_out.fft_from(fphis);
   auto tp_3 = std::chrono::high_resolution_clock::now();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   // FFT takes ~50 times more time than the adjust of the Fs.
   // std::cout << "::::::: Timings " << d21 << " " << d32 << " milliseconds"  << std::endl;
   return xmap_out;
}

clipper::Xmap<float>
coot::util::sharpen_blur_map_with_resample(const clipper::Xmap<float> &xmap_in,
                                           float b_factor, float resample_factor) {

   // Normal x-ray maps have a resample value of 1.5
   // Normal EM maps have a resample value of 1.0

   // So to get a "nice" EM map call with (xmap, 0, 1.5)


   if (resample_factor >= 1.0) {
      // normal case
      // resample_factor (say 2) means finer grid and higher resolution
      float mg = coot::util::max_gridding(xmap_in);
      std::cout << "INFO:: Map max gridding " << mg << " A/grid-point" << std::endl;
      clipper::Resolution reso(2.0 * mg); // for map -> data, the resolution is half the gridding
      clipper::HKL_info myhkl(xmap_in.spacegroup(), xmap_in.cell(), reso, true);
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
      clipper::Grid_sampling gs(xmap_in.spacegroup(), xmap_in.cell(), reso, resample_factor);
      clipper::Xmap<float> xmap_out(xmap_in.spacegroup(), xmap_in.cell(), gs);
      xmap_in.fft_to(fphis);
      std::cout << "DEBUG:: n-reflections: input map " << fphis.num_obs() << " reso-limit " << reso.limit() << " A" << std::endl;
      clipper::HKL_info::HKL_reference_index hri;

      auto tp_1 = std::chrono::high_resolution_clock::now();
      for (hri = fphis.first(); !hri.last(); hri.next()) {
         float f = fphis[hri].f();
         if (! clipper::Util::is_nan(f)) {
            float irs = hri.invresolsq();
            fphis[hri].f() *= exp(-b_factor * irs * 0.25);
         }
      }
      auto tp_2 = std::chrono::high_resolution_clock::now();
      xmap_out.fft_from(fphis);
      auto tp_3 = std::chrono::high_resolution_clock::now();
      auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
      auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
      return xmap_out;
   } else {
      return sharpen_blur_map_with_reduced_sampling(xmap_in, b_factor, resample_factor);
   }
}

clipper::Xmap<float>
coot::util::sharpen_blur_map_with_reduced_sampling(const clipper::Xmap<float> &xmap_in,
                                                   float b_factor, float resample_factor) {

      // cut the resolution (by using a more coarsely-sampled map)

      float mg = coot::util::max_gridding(xmap_in);
      std::cout << "INFO:: Map max gridding " << mg << " A/grid-point" << std::endl;

      // if resample_factor is less than 1.0, we want a less high resolution map
      clipper::Resolution reso_in(2.0 * mg);
      clipper::Resolution reso_out(2.0 * mg / resample_factor);
      clipper::Grid_sampling gs_out(xmap_in.spacegroup(), xmap_in.cell(), reso_out, 1.0);
      clipper::Xmap<float> xmap_out(xmap_in.spacegroup(), xmap_in.cell(), gs_out);
      clipper::HKL_info  input_map_hkl_info(xmap_in.spacegroup(),  xmap_in.cell(),  reso_in,  true);
      clipper::HKL_info output_map_hkl_info(xmap_out.spacegroup(), xmap_out.cell(), reso_out, true);
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_in ( input_map_hkl_info);
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_out(output_map_hkl_info);
      xmap_in.fft_to(fphis_in);

      // transfer fphis within the resolution limit from fphis_in to fphis_out

      std::cout << "DEBUG:: n-reflections: input map " << fphis_in.num_obs() << " output map: " << fphis_out.num_obs() << std::endl;
      clipper::HKL_info::HKL_reference_index hri;
      float max_reso = reso_out.invresolsq_limit();
      for (hri = fphis_out.first(); !hri.last(); hri.next()) {
         clipper::HKL hkl = hri.hkl();
         fphis_out[hri] = fphis_in[hkl];
      }
      for (hri = fphis_out.first(); !hri.last(); hri.next()) {
         float irs = hri.invresolsq();
         if (irs < max_reso) {
            fphis_out[hri].f() *= exp(-b_factor * irs * 0.25);
         } else {
            fphis_out[hri].f() = 0.0;
         }
      }
      xmap_out.fft_from(fphis_out);
      return xmap_out;

}



void
coot::util::multi_sharpen_blur_map(const clipper::Xmap<float> &xmap_in,
                                   const std::vector<float> &b_factors,
                                   std::vector<clipper::Xmap<float> > *xmaps_p) {

   float mg = coot::util::max_gridding(xmap_in);
   clipper::Resolution reso(2.0 * mg);
   clipper::HKL_info myhkl(xmap_in.spacegroup(), xmap_in.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
   xmap_in.fft_to(fphis);
   clipper::HKL_info::HKL_reference_index hri;

   for (std::size_t i=0; i<b_factors.size(); i++) {
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_loop = fphis;
      xmaps_p->at(i).init(xmap_in.spacegroup(), xmap_in.cell(), xmap_in.grid_sampling());
      const float &b_factor = b_factors[i];
      for (hri = fphis_loop.first(); !hri.last(); hri.next()) {
         float f = fphis[hri].f();
         if (! clipper::Util::is_nan(f)) {
            float irs =  hri.invresolsq();
            fphis_loop[hri].f() *= exp(-b_factor * irs * 0.25);
         }
      }
      xmaps_p->at(i).fft_from(fphis_loop);
   }
}


// sharpen/blur self
void
coot::util::sharpen_blur_map(clipper::Xmap<float> *xmap_p, float b_factor) {

   clipper::Xmap<float> &xmap(*xmap_p);
   float mg = coot::util::max_gridding(xmap);
   clipper::Resolution reso(2.0 * mg);
   clipper::HKL_info myhkl(xmap.spacegroup(), xmap.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
   clipper::Xmap<float> xmap_out(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   xmap.fft_to(fphis);
   clipper::HKL_info::HKL_reference_index hri;

   int count = 0;
   auto tp_1 = std::chrono::high_resolution_clock::now();
   for (hri = fphis.first(); !hri.last(); hri.next()) {
      float f = fphis[hri].f();
      if (! clipper::Util::is_nan(f)) {
         float irs =  hri.invresolsq();
         fphis[hri].f() *= exp(-b_factor * irs * 0.25);
         count++;
      }
   }
   auto tp_2 = std::chrono::high_resolution_clock::now();
   xmap.fft_from(fphis);
   auto tp_3 = std::chrono::high_resolution_clock::now();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   std::cout << "INFO:: sharpen_blur self: Timings " << d21 << " " << d32 << " milliseconds"  << std::endl;
}



#include <gsl/gsl_fit.h>


float
coot::util::b_factor(const std::vector<coot::amplitude_vs_resolution_point> &fsqrd_data,
                     std::pair<bool, float> reso_low_invresolsq,
                     std::pair<bool, float> reso_high_invresolsq) {

   // we want to chop off data that are zeros (from over-sampled maps (routine, of course)):
   // check if the current data point is 200 times smaller (say) than the
   // previous resolution range.

   std::cout << "debug:: b_factor() fsqrd_data size " << fsqrd_data.size() << std::endl;

   float b = 0.0f;
   std::vector<std::pair<double, double> > data;
   data.reserve(fsqrd_data.size());
   float prev_log = -100.0f;
   for (std::size_t i=0; i<fsqrd_data.size(); i++) {
      const amplitude_vs_resolution_point &d = fsqrd_data[i];
      float reso = d.get_invresolsq();
      float lf = log10(d.get_average_fsqrd());
      if (true)
         std::cout << "debug::raw " << d.count << " " << reso << " " << lf << " "
                   << reso_low_invresolsq.first << " " << reso_low_invresolsq.second << " "
                   << reso_high_invresolsq.first << " " << reso_high_invresolsq.second << std::endl;
      if (d.count > 0) {
         if (!reso_low_invresolsq.first  || reso >= reso_low_invresolsq.second) {
            if (!reso_high_invresolsq.first || reso <= reso_high_invresolsq.second) {
               if (lf > (prev_log-2.3)) {
                  std::pair<double, double> p(reso, lf);
                  data.push_back(p);
                  prev_log = lf;
               } else {
                  // no more data
                  std::cout << "breaking on " << reso << " " << lf << std::endl;
                  break;
               }
            }
         }
      }
   }

   std::cout << "debug:: b_fact(): data size " << data.size() << std::endl;
   if (data.size() > 1) {
      unsigned int n = data.size();
      double *x_p = new double[n];
      double *y_p = new double[n];
      for (std::size_t i=0; i<data.size(); i++) {
         std::cout << "debug::b-factor estimation: adding graph data " << data[i].first << " " << data[i].second << std::endl;
         x_p[i] = data[i].first;
         y_p[i] = data[i].second;
      }
      double cov00, cov01, cov11, sum_sq;
      double c_0, c_1; // c and m
      gsl_fit_linear(x_p, 1, y_p, 1, n, &c_0, &c_1, &cov00, &cov01, &cov11, &sum_sq);
      // b = -8.0 * M_PI * M_PI * c_1;
      b = -0.5 * c_1;
      delete [] x_p;
      delete [] y_p;
   }

   return b;
}


// if n_bins is -1, let the function decide how many bins
//
// actually, we return bins of amplitude squares
//
std::vector<coot::amplitude_vs_resolution_point>
coot::util::amplitude_vs_resolution(const clipper::Xmap<float> &xmap_in,
                                    int n_bins_in) {

   std::vector<coot::amplitude_vs_resolution_point> v;
   int n_bins = n_bins_in;
   if (n_bins_in == -1)
      n_bins = 60;
   v.resize(n_bins);
   if (n_bins < 1) return v;

   float mg = coot::util::max_gridding(xmap_in);
   clipper::Resolution reso(3.0 * mg); // tricky number - crystallographic maps are oversampled
                                       // by at least this: 2.0 x 1.5
   clipper::HKL_info myhkl(xmap_in.spacegroup(), xmap_in.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
   clipper::Xmap<float> xmap_out(xmap_in.spacegroup(), xmap_in.cell(), xmap_in.grid_sampling());
   xmap_in.fft_to(fphis);
   clipper::HKL_info::HKL_reference_index hri;

   float irs_max = 0.0f;

   for (hri = fphis.first(); !hri.last(); hri.next()) {
      float f = fphis[hri].f();
      if (! clipper::Util::is_nan(f)) {
         float irs =  hri.invresolsq();
         if (irs > irs_max)
            irs_max = irs;
      }
   }

   for (hri = fphis.first(); !hri.last(); hri.next()) {
      float f = fphis[hri].f();
      if (! clipper::Util::is_nan(f)) {
         float irs =  hri.invresolsq();
         float res_frac = irs/irs_max;
         int bin = static_cast<int> (n_bins * res_frac);
         if (bin == n_bins) {
            bin = n_bins-1;
         }
         v[bin].add(f*f, irs);
      }
   }

   for (std::size_t i=0; i<v.size(); i++)
      v[i].finish(); // calculate averages

   return v;
}


clipper::Xmap<float>
coot::util::transform_map(const clipper::Xmap<float> &xmap_in,
                          const clipper::Spacegroup &new_space_group,
                          const clipper::Cell &new_cell,
                          const clipper::RTop_orth &rtop,
                          const clipper::Coord_orth &to_pt,
                          float box_size) {

   // we now need to create about_pt: i.e. where the map is pulled *from*
   clipper::Coord_orth about_pt = to_pt.transform(rtop.inverse());

   clipper::Xmap<float> xmap;
   clipper::Grid_sampling new_gs = coot::util::suggested_grid_sampling(xmap_in.grid_sampling(),
                                                                       xmap_in.cell(),
                                                                       new_space_group,
                                                                       new_cell);

   std::cout << "INFO:: creating new map for transformed map with spacegroup: " << new_space_group.symbol_hm()
             << " cell: " << new_cell.format() << " grid-sampling " << new_gs.format()
             << std::endl;

   xmap.init(new_space_group, new_cell, new_gs);

   std::cout << "INFO:: coord info:         to_pt: " << to_pt.format()    << std::endl;
   std::cout << "INFO:: coord info:      about_pt: " << about_pt.format() << std::endl;

   clipper::Grid_sampling grid = xmap.grid_sampling();
   clipper::Grid_range gr(xmap.cell(), xmap.grid_sampling(), box_size);
   clipper::Coord_grid g, g0, g1;
   clipper::RTop_orth rtop_inv = rtop.inverse();
   typedef clipper::Xmap<float>::Map_reference_coord MRC;
   MRC i0, iu, iv, iw;
   g = to_pt.coord_frac(new_cell).coord_grid(new_gs);

   if (false) { // 20250128-PE  have you got the correct cell angles?
      std::cout << "INFO:: new_cell: " << new_cell.format() << std::endl;
      std::cout << "INFO:: new_gs: " << new_gs.format() << std::endl;
      std::cout << "INFO:: to_pt: " << to_pt.format() << std::endl;
      clipper::Coord_frac cf = to_pt.coord_frac(new_cell);
      std::cout << "INFO:: to_pt.coord_frac: " << cf.format() << std::endl;
      std::cout << "INFO:: to_pt.coord_frac.coord_grid: " << to_pt.coord_frac(new_cell).coord_grid(new_gs).format()
                << std::endl;
   }

   std::cout << "DEBUG:: pulling map from point:   " << about_pt.format() << std::endl;
   std::cout << "DEBUG:: creating map about point: " << to_pt.format() << std::endl;

   std::cout << "DEBUG:: grid point g: " << g.format() << std::endl;
   std::cout << "DEBUG:: grid range gr: " << gr.format() << std::endl;
   std::cout << "DEBUG:: grid range gr.min: " << gr.min().format() << std::endl;
   std::cout << "DEBUG:: grid range gr.max: " << gr.max().format() << std::endl;

   std::vector<coot::util::map_ref_triple_t> density_points;

   g0 = g + gr.min();
   g1 = g + gr.max();
   i0 = MRC( xmap, g0 );

   clipper::Coord_orth iw_pos;
   clipper::Coord_orth dpt;
   double d2;

   for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
         for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
            iw_pos = iw.coord().coord_frac(grid).coord_orth(xmap.cell());
            dpt =    iw_pos.transform(rtop_inv);
            d2 = (iw_pos - to_pt).lengthsq();
            density_points.push_back(coot::util::map_ref_triple_t(d2, iw, coot::util::density_at_point(xmap_in, dpt)));
         }
   std::sort(density_points.begin(), density_points.end());
   for (unsigned int ic=0; ic<density_points.size(); ic++)
      xmap[density_points[ic].iw] = density_points[ic].density;

   return xmap;
}


clipper::Grid_sampling
coot::util::suggested_grid_sampling(const clipper::Grid_sampling &orig_sampling,
                                   const clipper::Cell &orig_cell,
                                   const clipper::Spacegroup &new_space_group,
                                   const clipper::Cell &new_cell) {

   float sampling_a = orig_cell.a()/float(orig_sampling.nu());
   float sampling_b = orig_cell.b()/float(orig_sampling.nv());
   float sampling_c = orig_cell.c()/float(orig_sampling.nw());

   float best_sampling = sampling_a;
   if (sampling_b < best_sampling)
      best_sampling = sampling_b;
   if (sampling_c < best_sampling)
      best_sampling = sampling_c;

   clipper::Resolution resolution(best_sampling * 3.0);

   return clipper::Grid_sampling(new_space_group, new_cell, resolution, 2.0);
}



std::pair<float, float>
coot::util::mean_and_variance(const clipper::Xmap<float> &xmap) {

   clipper::Xmap_base::Map_reference_index ix;
   double sum = 0.0;
   double sum_sq = 0.0;
   int npoints = 0;
   float d;
   for (ix = xmap.first(); !ix.last(); ix.next() ) {
      npoints++;
      d = xmap[ix];
      sum += d;
      sum_sq += d*d;
   }
   double mean = 0.0;
   double var = -1.0;

   if (npoints > 0) {
      mean = sum/float(npoints);
      var = sum_sq/float(npoints) - mean*mean;
   }

   return std::pair<float, float> (mean, var);
}


clipper::Xmap<float>
coot::util::laplacian_transform(const clipper::Xmap<float> &xmap_in) {

   clipper::Xmap<float> laplacian = xmap_in;

   clipper::Coord_map pos;
   float val;
   clipper::Grad_map<float> grad;
   clipper::Curv_map<float> curv;

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_in.first(); !ix.last(); ix.next())  {
      // xmap_in.interp_curv(pos, val, grad, curv);
      clipper::Interp_cubic::interp_curv(xmap_in, ix.coord().coord_map(), val, grad, curv);
      val = curv.det();
      laplacian[ix] = -val;
   }

   return laplacian;
}

// Spin the torsioned atom round the rotatable bond and find the
// orientation (in degrees) that is in the highest density.
//
// return two torsions, the first is relative to the current position,
// the second is relative to the first atom (pa1).
// Return first -1111 (less than -1000) on failure
// Atoms in tors are CA (or N), CA, CB, CG.
//
std::pair<float, float>
coot::util::spin_search(const clipper::Xmap<float> &xmap, mmdb::Residue *res, coot::torsion tors) {

   // The plan:
   //
   // 1) find 4 atoms in residue that correspond to the torsion
   //
   // 2) Use rotate point about vector to generate new points
   //
   //      Test density at new points

   float best_ori = -1111.1; //returned thing
   float torsion_relative_to_N = -1111.1;

   std::vector<mmdb::Atom * > match_atoms = tors.matching_atoms(res);

   if (match_atoms.size() != 4) {
      std::cout << "ERROR:: not all atoms for torsion found in residue!" << std::endl;
      std::cout << "        (found " << match_atoms.size() << " atoms.)" << std::endl;
   } else {

      clipper::Coord_orth pa1(match_atoms[0]->x, match_atoms[0]->y, match_atoms[0]->z);
      clipper::Coord_orth pa2(match_atoms[1]->x, match_atoms[1]->y, match_atoms[1]->z);
      clipper::Coord_orth pa3(match_atoms[2]->x, match_atoms[2]->y, match_atoms[2]->z);
      clipper::Coord_orth pa4(match_atoms[3]->x, match_atoms[3]->y, match_atoms[3]->z);

      float best_d = -99999999.9;
      clipper::Coord_orth best_pos;
      for (double theta=0; theta <=360; theta+=3.0) {

         clipper::Coord_orth dir   = pa3 - pa2;
         clipper::Coord_orth pos   = pa4;
         clipper::Coord_orth shift = pa3;
         clipper::Coord_orth co = rotate_around_vector(dir, pos, shift, clipper::Util::d2rad(theta));
         float this_d = density_at_point(xmap, co);
         if (this_d > best_d) {
            best_d = this_d;
            best_ori = theta;
            best_pos = co;
         }
      }
      double ct = clipper::Coord_orth::torsion(pa1, pa2, pa3, best_pos);
      torsion_relative_to_N = clipper::Util::rad2d(ct);
   }


   std::pair<float, float> p(best_ori, torsion_relative_to_N);
   return p;
}

// return a map and its standard deviation.  scale is applied to
// map_in_2 before substraction.
std::pair<clipper::Xmap<float>, float>
coot::util::difference_map(const clipper::Xmap<float> &xmap_in_1,
                           const clipper::Xmap<float> &xmap_in_2,
                           float map_scale) {
   clipper::Xmap<float> r = xmap_in_1;
   float std_dev = 0.2;

    clipper::Xmap_base::Map_reference_index ix;
    for (ix = r.first(); !ix.last(); ix.next())  {
       clipper::Coord_map  cm1 = ix.coord().coord_map();
       clipper::Coord_frac cf1 = cm1.coord_frac(xmap_in_1.grid_sampling());
       clipper::Coord_orth co  = cf1.coord_orth(xmap_in_1.cell());
       clipper::Coord_frac cf2 =  co.coord_frac(xmap_in_2.cell());
       clipper::Coord_map  cm2 = cf2.coord_map(xmap_in_2.grid_sampling());
       float map_2_val;
       clipper::Interp_cubic::interp(xmap_in_2, cm2, map_2_val);
       r[ix] = xmap_in_1[ix] - map_scale * map_2_val;
    }
   return std::pair<clipper::Xmap<float>, float> (r, std_dev);
}

// make a copy of map_in, but in the cell and gridding of reference_map
clipper::Xmap<float>
coot::util::reinterp_map(const clipper::Xmap<float> &xmap_in,
                         const clipper::Xmap<float> &reference_xmap) {

   clipper::Xmap<float> rmap;
   rmap.init(reference_xmap.spacegroup(), reference_xmap.cell(), reference_xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = rmap.first(); !ix.last(); ix.next())  {
      clipper::Coord_map  cm1 = ix.coord().coord_map();
      clipper::Coord_frac cf1 = cm1.coord_frac(rmap.grid_sampling());
      clipper::Coord_orth co  = cf1.coord_orth(rmap.cell());
      clipper::Coord_frac cf2 =  co.coord_frac(xmap_in.cell());
      clipper::Coord_map  cm2 = cf2.coord_map( xmap_in.grid_sampling());
      float map_2_val;
      clipper::Interp_cubic::interp(xmap_in, cm2, map_2_val);
      rmap[ix] = map_2_val;
   }
   return rmap;

}


clipper::Xmap<float>
coot::util::average_map(const std::vector<std::pair<clipper::Xmap<float>, float> > &maps_and_scales_vec) {

   clipper::Xmap<float> rmap;

   if (maps_and_scales_vec.size() > 0) {

      for (unsigned int imap=0; imap<maps_and_scales_vec.size(); imap++)
         std::cout << "INFO:: multiplying map (function index) " << imap << " "
                   << maps_and_scales_vec[imap].first.grid_sampling().format()
                   << " by " << maps_and_scales_vec[imap].second << std::endl;

      // set the first map and scale it.
      rmap = maps_and_scales_vec[0].first;
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = rmap.first(); !ix.last(); ix.next())
         rmap[ix] *= maps_and_scales_vec[0].second;

      for (unsigned int iothers=1; iothers<maps_and_scales_vec.size(); iothers++) {
         for (ix = rmap.first(); !ix.last(); ix.next())  {
            clipper::Coord_map  cm1 = ix.coord().coord_map();
            clipper::Coord_frac cf1 = cm1.coord_frac(rmap.grid_sampling());
            clipper::Coord_orth co  = cf1.coord_orth(rmap.cell());
            clipper::Coord_frac cf2 =  co.coord_frac(maps_and_scales_vec[iothers].first.cell());
            clipper::Coord_map  cm2 = cf2.coord_map(maps_and_scales_vec[iothers].first.grid_sampling());
            float map_2_val;
            clipper::Interp_cubic::interp(maps_and_scales_vec[iothers].first, cm2, map_2_val);
            rmap[ix] += map_2_val * maps_and_scales_vec[iothers].second;
         }
      }
      float scale_sum = 0.0;
      for (unsigned int i=0; i<maps_and_scales_vec.size(); i++) {
         scale_sum += maps_and_scales_vec[i].second;
      }
      if (scale_sum != 0.0) {
         float sf = 1.0/scale_sum;
         for (ix = rmap.first(); !ix.last(); ix.next())  {
            rmap[ix] *= sf;
         }
      }
   }
   return rmap;
}

// like above, but modify the map, don't return a new one. Also this presumes that the maps haave the same gridding
void
coot::util::regen_weighted_map(clipper::Xmap<float> *xmap_in,
                              const std::vector<std::pair<clipper::Xmap<float> *, float> > &maps_and_scales_vec) {

   if (!maps_and_scales_vec.empty()) {
      for (unsigned int imap=0; imap<maps_and_scales_vec.size(); imap++) {
         const clipper::Xmap<float> &xmap = *maps_and_scales_vec[imap].first;
         const float &sf                  = maps_and_scales_vec[imap].second;
         clipper::Xmap_base::Map_reference_index ix;
         for (ix = xmap_in->first(); !ix.last(); ix.next()) {
            if (imap == 0) {
               (*xmap_in)[ix] = xmap[ix] *sf;
            } else {
               (*xmap_in)[ix] += xmap[ix] *sf;
            }
         }
      }
   }
}


// like above, but variance, scales are ignored.
//
clipper::Xmap<float>
coot::util::variance_map(const std::vector<std::pair<clipper::Xmap<float>, float> > &maps_and_scales_vec) {

   clipper::Xmap<float> rmap;
   clipper::Xmap<float> sum_map;
   clipper::Xmap<float> sum_squares_map;
   if (maps_and_scales_vec.size()) {
      // set the first map
      int n_maps = maps_and_scales_vec.size();
      const clipper::Xmap<float> &first_map = maps_and_scales_vec[0].first;
      sum_map = first_map;
      sum_squares_map.init(first_map.spacegroup(), first_map.cell(), first_map.grid_sampling());
      rmap.init(           first_map.spacegroup(), first_map.cell(), first_map.grid_sampling());
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = sum_map.first(); !ix.last(); ix.next())  {
         sum_squares_map[ix] = sum_map[ix] * sum_map[ix];
      }
      for (unsigned int iothers=1; iothers<maps_and_scales_vec.size(); iothers++) {
         for (ix = sum_map.first(); !ix.last(); ix.next())  {
            clipper::Coord_map  cm1 = ix.coord().coord_map();
            clipper::Coord_frac cf1 = cm1.coord_frac(sum_map.grid_sampling());
            clipper::Coord_orth co  = cf1.coord_orth(sum_map.cell());
            clipper::Coord_frac cf2 =  co.coord_frac(maps_and_scales_vec[iothers].first.cell());
            clipper::Coord_map  cm2 = cf2.coord_map(maps_and_scales_vec[iothers].first.grid_sampling());
            float map_2_val;
            clipper::Interp_cubic::interp(maps_and_scales_vec[iothers].first, cm2, map_2_val);
            float v = map_2_val;
            sum_map[ix] += v;
            sum_squares_map[ix] += v * v;
         }
      }
      float f = 1.0/float(n_maps);
      for (ix = rmap.first(); !ix.last(); ix.next()) {
         float mean = sum_map[ix] * f;
         rmap[ix] = sum_squares_map[ix] * f  - mean*mean;
      }
   }
   return rmap;
}



// Delete atoms that don't have this alt conf (or "" alt conf).
void
coot::util::backrub_residue_triple_t::trim_this_residue_atoms() {

   std::vector<std::string> vec;
   trim_residue_atoms_generic(this_residue, vec, 0);
}

void
coot::util::backrub_residue_triple_t::trim_residue_atoms_generic(mmdb::Residue *residue_p,
                                                                 std::vector<std::string> keep_atom_vector,
                                                                 bool use_keep_atom_vector) {

   // Hmmm... in other places, I don't go round the houses (make the
   // vector, then delete).  I just DeleteAtom(), so it may not be
   // necessary here.

   if (residue_p) {
      std::vector<int> delete_atom_index_vec;
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
         std::string atom_name(residue_atoms[i]->name);
         std::string atom_alt_conf(residue_atoms[i]->altLoc);

         bool delete_this_atom_flag = 1;
         if (use_keep_atom_vector) {
            for (unsigned int ikeep=0; ikeep<keep_atom_vector.size(); ikeep++) {
               if (atom_name == keep_atom_vector[ikeep]) {
                  delete_this_atom_flag = 0;
                  break;
               }
            }
         } else {
            // i.e. is the this_residue, keep all of these atoms
            delete_this_atom_flag = 0;
         }
         if (delete_this_atom_flag) {
            delete_atom_index_vec.push_back(i);
         } else {
            if ((atom_alt_conf != alt_conf) &&
                (atom_alt_conf != "")) {
               delete_atom_index_vec.push_back(i);
            }
         }
      }

      if (delete_atom_index_vec.size() > 0) {
         for (unsigned int i=0; i<delete_atom_index_vec.size(); i++)
            residue_p->DeleteAtom(i);
         residue_p->TrimAtomTable();
      }
   }
}



// As above, and also delete all atoms that are not C or O
void
coot::util::backrub_residue_triple_t::trim_prev_residue_atoms() {

   std::vector<std::string> vec;
   vec.push_back(" C  ");
   vec.push_back(" O  ");
   trim_residue_atoms_generic(this_residue, vec, 0);
}

   // As trim_this_residue_atoms, and also delete all atoms that are not N.
void
coot::util::backrub_residue_triple_t::trim_next_residue_atoms() {

   std::vector<std::string> vec;
   vec.push_back(" N  ");
   vec.push_back(" H  ");
   trim_residue_atoms_generic(this_residue, vec, 0);
}



// as in the verb, not the noun., return the number of segments (0 is
// also useful segment).
//
std::pair<int, clipper::Xmap<int> >
coot::util::segment_map::segment_emsley_flood(const clipper::Xmap<float> &xmap, float low_level) {

   clipper::Xmap<int> xmap_int(xmap.spacegroup(),
                               xmap.cell(),
                               xmap.grid_sampling());

   int UNASSIGNED = -1;
   int TOO_LOW    = -2;

   // how many points are there in an xmap are there above the
   // (user-defined) noise level (low_level)?
   //
   long n_points = 0;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next() )  {
      if (xmap[ix] < low_level) {
         xmap_int[ix] = TOO_LOW;
      } else {
         xmap_int[ix] = UNASSIGNED;
         n_points++;
      }
   }

   std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float> > density_values(n_points);
   long int i=0;
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      if (xmap[ix] >= low_level) {
         density_values[i] = std::pair<clipper::Xmap_base::Map_reference_index, float> (ix, xmap[ix]);
         i++;
      }
   }

   std::sort(density_values.begin(), density_values.end(), compare_density_values_map_refs);

   int i_segment_index = 0;
   clipper::Coord_grid c_g;

   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.5, 3.1); // 3x3x3 cube, not centre
   for (i=0; i<n_points; i++) {

      // Flood down points around xmap[density_values[i].first] (if
      // they are not already in the neigbourhood of a peak).
      //
      //
      if (xmap_int[density_values[i].first] == UNASSIGNED) {

         // OK, so this point starts a new segment

         int q_set = 0; // debugging counter
         float v = xmap[density_values[i].first]; // all points in the segment should be
                                                  // lower than (or equal to) this density
                                                  // at the  "centre" of the segment blob.
         xmap_int[density_values[i].first] = i_segment_index;
         std::queue<clipper::Coord_grid> q;
         for (int i_n=0; i_n<neighb.size(); i_n++) {
            c_g = density_values[i].first.coord() + neighb[i_n];
            if (xmap_int.get_data(c_g) == UNASSIGNED) {
               if (xmap.get_data(c_g) <= v) {
                  xmap_int.set_data(c_g, i_segment_index);
                  q.push(c_g);
                  q_set++;
               }
            }
         }

         while (q.size()) {
            clipper::Coord_grid c_g_start = q.front(); // new local centre
            v = xmap.get_data(c_g_start);
            q.pop();
            for (int i_n=0; i_n<neighb.size(); i_n++) {
               c_g = c_g_start + neighb[i_n];
               if (xmap_int.get_data(c_g) == UNASSIGNED) {
                  if (xmap.get_data(c_g) <= v) {

                     // But only set this neighbour if the steepest ascent from this
                     // neighbour is c_g. (otherwise we run into problems of crossing
                     // boundary at the border between segments when we take the steepest
                     // path up from them)

                     clipper::Coord_grid steepest_neighb(0,0,0);
                     float biggest = -1.0;
                     for (int i_nn=0; i_nn<neighb.size(); i_nn++) {
                        float tv = xmap.get_data(c_g + neighb[i_nn]);
                        if (tv > biggest) {
                           biggest = tv;
                           steepest_neighb = neighb[i_nn];
                        }
                     }

                     if ((c_g + steepest_neighb) == c_g_start) {

                        xmap_int.set_data(c_g, i_segment_index);
                        q.push(c_g);
                        q_set++;
                     }
                  }
               }
            }
         }
         i_segment_index++;
      }
   }
   int n_segments = i_segment_index;
   return std::pair<int, clipper::Xmap<int> > (n_segments, xmap_int);
}



// as in the verb, not the noun., return the number of segments (0 is
// also useful segment).
//
std::pair<int, clipper::Xmap<int> >
coot::util::segment_map::segment(const clipper::Xmap<float> &xmap, float low_level) {

   clipper::Xmap<int> xmap_int(xmap.spacegroup(),
                               xmap.cell(),
                               xmap.grid_sampling());

   int UNASSIGNED = -1;
   int TOO_LOW    = -2;

   // how many points are there in an xmap are there above the
   // (user-defined) noise level (low_level)?
   //
   long n_points = 0;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next() )  {
      if (xmap[ix] < low_level) {
         xmap_int[ix] = TOO_LOW;
      } else {
         xmap_int[ix] = UNASSIGNED;
         n_points++;
      }
   }

   std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float> > density_values(n_points);
   long int i=0;
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      if (xmap[ix] >= low_level) {
         density_values[i] = std::pair<clipper::Xmap_base::Map_reference_index, float> (ix, xmap[ix]);
         i++;
      }
   }

   std::sort(density_values.begin(), density_values.end(), compare_density_values_map_refs);

   int i_segment_index = 0;
   clipper::Coord_grid c_g;
   std::map<int, int> segment_id_counter_map;

   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.5, 3.1); // 3x3x3 cube, not centre
   for (i=0; i<n_points; i++) {

      if (xmap_int[density_values[i].first] == UNASSIGNED) {

         // Does this grid point have a neighbour which is already in
         // a segment?  If so, which one?  If there are many
         // neighbours (this is a watershed point) then we need to make
         // a note of all of them and choose later

         std::map<int, std::vector<clipper::Coord_grid> > neighbour_count_map; // (as in the stl map)
         for (int i_n=0; i_n<neighb.size(); i_n++) {
            c_g = density_values[i].first.coord() + neighb[i_n];
            int segment_index = xmap_int.get_data(c_g);
            if (segment_index >= 0) {
               neighbour_count_map[segment_index].push_back(c_g);
            }
         }

         if (neighbour_count_map.size()) {

            // OK, so add this grid point to a previously existing
            // segment.  Now, which one?
            //
            if (neighbour_count_map.size() == 1) {
               xmap_int[density_values[i].first] = neighbour_count_map.begin()->first;
               segment_id_counter_map[neighbour_count_map.begin()->first]++;
            } else {

               // The complicated case, we have to decide which is the
               // biggest of the neighbouring segments and add this
               // grid point to that segment.

               // int big_seg_id = find_biggest_segment(neighbour_count_map, segment_id_counter_map);
               int big_seg_id = find_smallest_segment(neighbour_count_map, segment_id_counter_map);
               if (0)
                  std::cout << "watershed point " << density_values[i].first.coord().format()
                            << " has " << neighbour_count_map.size() << " neighboring segments and "
                            << " is of segment " << big_seg_id << "\n";
               xmap_int[density_values[i].first] = big_seg_id;
               segment_id_counter_map[big_seg_id]++;

            }

         } else {

            // this didn't have any neighbours so we start a new segment:
            //
            xmap_int[density_values[i].first] = i_segment_index;
            segment_id_counter_map[i_segment_index]++;
            i_segment_index++; // for next round
         }
      }
   }
   resegment_watershed_points(&xmap_int, xmap); // Pintilie et al. didn't mention that we need this.

   int n_segments = i_segment_index;
   return std::pair<int, clipper::Xmap<int> > (n_segments, xmap_int);
}



// sorting function used by segmentation
//
// static
bool
coot::util::segment_map::compare_density_values_map_refs(const std::pair<clipper::Xmap_base::Map_reference_index, float> &v1,
                                            const std::pair<clipper::Xmap_base::Map_reference_index, float> &v2) {
   return (v2.second < v1.second);
}

// Pintilie et al. didn't mention that we need this.
void
coot::util::segment_map::resegment_watershed_points(clipper::Xmap<int> *xmap_int_in,
                                                    const clipper::Xmap<float> &xmap) const {

   clipper::Xmap<int> &xmap_int = *xmap_int_in;
   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.5, 3.1); // 3x3x3 cube, not centre

   clipper::Xmap_base::Map_reference_index ix;
   int is;
   int ns;
   clipper::Coord_grid c_g;

   for (ix = xmap_int.first(); !ix.last(); ix.next()) {
      is = xmap_int[ix];
      if (is >= 0) {
         std::map<int, int> segment_id_map;
         for (int i_n=0; i_n<neighb.size(); i_n++) {
            c_g = ix.coord() + neighb[i_n];
            ns = xmap_int.get_data(c_g);
            if (ns >= 0)
               segment_id_map[ns]++;
         }
         if (segment_id_map.size() > 1) {
            // OK, we have a watershed point, we need to (potentially) change the segment to the
            // one which has the steepest gradient from this point.
            float v = xmap[ix];
            float best_vn = -1;
            clipper::Coord_grid best_n;
            for (int i_n=0; i_n<neighb.size(); i_n++) {
               c_g = ix.coord() + neighb[i_n];
               float vn = xmap.get_data(c_g);
               if (v > best_vn) {
                  if (vn > v) {
                     best_vn = vn;
                     best_n = neighb[i_n];
                  }
               }
            }
            // OK, best_vn was set
            if (best_vn > -0.9) {
               int i_seg_neighb = xmap_int.get_data(ix.coord() + best_n);
               if (xmap_int[ix] != i_seg_neighb) {
                  if (0)
                     std::cout << "Resegmenting " << ix.coord().format() << " from " << xmap_int[ix] << " to "
                               << i_seg_neighb << "\n";
                  xmap_int[ix] = i_seg_neighb;
               }
            }
         }
      }
   }

}


// "scale-space" segmentation (i.e. merge segments by progressive blurring)
//
std::pair<int, clipper::Xmap<int> >
coot::util::segment_map::segment(const clipper::Xmap<float> &xmap_in,
                                 float low_level,
                                 float gaussian_sigma, // per round
                                 int n_rounds) {

   std::cout << "DEBUG:: start of segment with low_level " << low_level
             << " gaussian_sigma " << gaussian_sigma
             << " n_rounds " << n_rounds
             << std::endl;

   // This algorithm is critically dependent on the gradient around
   // points on the borders between segments (the watershed regions).
   //
   // If we
   //
   // 1) segment then input map
   //
   // 2) do a high-res fft of map, then re-fft to get a new map and
   //    then segment that
   //
   // then the segmentation is very different.
   //
   // Well, so be it.
   //
   // So, the first map that we segment should use exactly the same
   // sfs throughout, the only difference being a blurring gaussian
   // function which scales the f values.
   //
   // the f and phis for the map are fphis.  Each round they get
   // manipulated and they are fphis_for_loop, a local copy
   //

   int n_segments = 0;


   clipper::Xmap<std::pair<bool, int> > segmented;
   segmented.init(xmap_in.spacegroup(),
                  xmap_in.cell(),
                  xmap_in.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;

   float mg = coot::util::max_gridding(xmap_in);
   // std::cout << "We got max gridding " << mg << std::endl;
   clipper::Resolution reso(mg); // nice high res.
   clipper::HKL_info myhkl(xmap_in.spacegroup(), xmap_in.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> >       fphis(myhkl);
   clipper::Xmap<float> xmap_first(xmap_in.spacegroup(), xmap_in.cell(), xmap_in.grid_sampling());
   xmap_in.fft_to(fphis);
   xmap_first.fft_from(fphis);
   std::pair<int, clipper::Xmap<int> > segmented_first = segment(xmap_first, low_level);
   for (ix = segmented_first.second.first(); !ix.last(); ix.next())
      segmented[ix] = std::pair<bool, int> (0 , segmented_first.second[ix]);

   // debugging
   std::cout << " ===================== initial segmentation: ===================\n";
   std::map<int, int> debug_segment_id_map;
   for (ix=segmented.first(); !ix.last(); ix.next()) {
      debug_segment_id_map[segmented[ix].second]++;
   }
   std::map<int, int>::const_iterator dit;
   for (dit=debug_segment_id_map.begin(); dit!=debug_segment_id_map.end(); ++dit)
      std::cout << "  segment " << dit->first << " has " << dit->second << " grid points" << std::endl;
   std::cout << " ========================================\n";


   for (int i_round=0; i_round<n_rounds; i_round++) {

      std::cout << "INFO Segmentation: Round " << i_round << " of " << n_rounds << std::endl;

      // float round_factor = 1.0/float(i_round + 1);
      float round_factor = (i_round + 1);
      clipper::HKL_data< clipper::datatypes::F_phi<float> >  fphis_for_loop = fphis;

      clipper::HKL_info::HKL_reference_index hri;
      for (hri = fphis_for_loop.first(); !hri.last(); hri.next()) {
         if (! clipper::Util::is_nan(fphis_for_loop[hri].f())) {
            float irs =  hri.invresolsq();
            float scale = expf(- round_factor * irs * irs/(2 * gaussian_sigma * gaussian_sigma));
            scale = 1;
            fphis_for_loop[hri].f() *= scale;
            if (0)
               std::cout << " new f: " << hri.hkl().format() << " " <<  fphis_for_loop[hri].f()
                         << " reso= " << irs
                         << "   and scale " << scale << "\n";
         }
      }

      clipper::Xmap<float> xmap_new(xmap_in.spacegroup(),
                                    xmap_in.cell(),
                                    xmap_in.grid_sampling());

      xmap_new.fft_from(fphis_for_loop);

      for (ix = xmap_new.first(); !ix.last(); ix.next()) {
         if (segmented[ix].second != TOO_LOW) {

            // assert
            if (segmented[ix].second == UNASSIGNED)
               std::cout << "  segmentation failure, UNASSIGNED value"
                         << ix.coord().format() << std::endl;

            // path includes the starting point ix.coord();
            std::vector<clipper::Coord_grid> path = path_to_peak(ix.coord(), xmap_new);
            if (path.size() > 1) {
               int segmented_start = segmented[ix].second;
               int segmented_local_peak = segmented.get_data(path.back()).second;

               if (segmented_local_peak != segmented_start) {

                  if (0) {
                     std::cout << " ============= yeay! map point "
                               << ix.coord().format() << " " << segmented_start
                               << " with path\n";
                     for (unsigned int i=0; i<path.size(); i++) {
                        std::cout << "    " << path[i].format() << " segment: "
                                  << segmented.get_data(path[i]).second << " value: "
                                  << xmap_new.get_data(path[i])
                                  << std::endl;
                     }
                     std::cout << " ending in segment: " << segmented_local_peak
                               << std::endl;
                  }

                  for (unsigned int ipath=0; ipath<path.size(); ipath++) {
                     if (segmented.get_data(path[ipath]).second != segmented_local_peak)
                        flood_fill_segmented_map(&segmented, xmap_in, path[ipath],
                                                 segmented_start, segmented_local_peak);
                  }
               }
            }
         }
      }
   }

   clipper::Xmap<int> clean_segmented;
   clean_segmented.init(segmented.spacegroup(),
                        segmented.cell(),
                        segmented.grid_sampling());
   clean_segmented = UNASSIGNED;
   std::map<int, int> segment_id_map; // as in the STL map, not clipper map

   // We want to eliminate gaps in the segmentation number: i.e. first
   // segment should be 0 and increment from there.
   std::map<int, int> segment_id_reindexer; // from -> to

   int seg_id;
   int new_index;
   std::map<int, int>::const_iterator it;
   int segment_counter = 0;

   for (ix=segmented.first(); !ix.last(); ix.next()) {
      seg_id = segmented[ix].second;
      if (seg_id == TOO_LOW) {
         clean_segmented[ix] = TOO_LOW;
      } else {
         if (seg_id == UNASSIGNED) {
            clean_segmented[ix] = UNASSIGNED;
         } else {

            it = segment_id_reindexer.find(seg_id);
            if (it != segment_id_reindexer.end()) {
               // it was found
               new_index = it->second;
            } else {
               segment_id_reindexer[seg_id] = segment_counter;
               new_index = segment_counter;
               segment_counter++;
            }
            clean_segmented[ix] = new_index;
         }
      }
   }
   n_segments = segment_counter+1;

   // Now we want to renumber the segments so that the segments with
   // the most points have the lowest segment numbers.

   for (ix=clean_segmented.first(); !ix.last(); ix.next())
      segment_id_map[clean_segmented[ix]]++;

   if (0) {
      std::cout << "======= pre-sort segments analysis =========" << std::endl;
      for (it=segment_id_map.begin(); it!=segment_id_map.end(); it++)
         std::cout << "  segment " << it->first << " has " << it->second << " grid points" << std::endl;
      std::cout << "===================================" << std::endl;
   }

   // create a vector for sorting.
   std::vector<std::pair<int, int> > segment_vec(segment_id_map.size());
   for (unsigned int icount=0; icount<segment_vec.size(); icount++) {
      segment_vec[icount] = std::pair<int, int> (segment_id_map[icount], icount);
   }
   std::sort(segment_vec.begin(), segment_vec.end(), sort_segment_vec);

   for (unsigned int i=0; i<segment_vec.size(); i++) {
      std::cout << "segment index:" << i << " with " << segment_vec[i].first << " grid points and "
                << "old index " << segment_vec[i].second << std::endl;
   }
   std::cout << "===================================" << std::endl;

   // create a map, fill it from the above sorted vector and give it
   // the sorted segment index.

   std::map<int, int> new_index_map;
   for (unsigned int i=0; i<segment_vec.size(); i++)
      new_index_map[segment_vec[i].second] = i;

   for (ix=clean_segmented.first(); !ix.last(); ix.next())
      if (clean_segmented[ix] != UNASSIGNED)
         if (clean_segmented[ix] != TOO_LOW)
            clean_segmented[ix] = new_index_map[clean_segmented[ix]];

   // ------------------------ post sort analysis -------------------------

   segment_id_map.clear();
   for (ix=clean_segmented.first(); !ix.last(); ix.next())
      segment_id_map[clean_segmented[ix]]++;

   std::cout << "======= post-sort segments analysis =========" << std::endl;
   for (it=segment_id_map.begin(); it!=segment_id_map.end(); it++)
      std::cout << "  " << it->first << " has " << it->second << " grid points" << std::endl;
   std::cout << "===================================" << std::endl;

   return std::pair<int, clipper::Xmap<int> > (n_segments, clean_segmented);


}


int
coot::util::segment_map::find_biggest_segment(const std::map<int, std::vector<clipper::Coord_grid> > &segment_id_map, const std::map<int, int> &segment_id_counter_map) const {

   // look through all of the segments in the segment id map and

   int seg_id_biggest = UNASSIGNED;
   int n_gp_in_biggest_segment = 0;
   std::map<int, std::vector<clipper::Coord_grid> >::const_iterator it;
   for (it=segment_id_map.begin(); it!=segment_id_map.end(); it++) {
      std::map<int, int>::const_iterator iti = segment_id_counter_map.find(it->first);
      if (iti != segment_id_counter_map.end()) {
         if (iti->second > n_gp_in_biggest_segment) {
            n_gp_in_biggest_segment = iti->second;
            seg_id_biggest = it->first;
         }
      }
   }

   return seg_id_biggest;
}

int
coot::util::segment_map::find_smallest_segment(const std::map<int, std::vector<clipper::Coord_grid> > &segment_id_map, const std::map<int, int> &segment_id_counter_map) const {

   // look through all of the segments in the segment id map and

   int seg_id_smallest = UNASSIGNED;
   int n_gp_in_smallest_segment = 65500;
   std::map<int, std::vector<clipper::Coord_grid> >::const_iterator it;
   for (it=segment_id_map.begin(); it!=segment_id_map.end(); it++) {
      std::map<int, int>::const_iterator iti = segment_id_counter_map.find(it->first);
      if (iti != segment_id_counter_map.end()) {
         if (iti->second < n_gp_in_smallest_segment) {
            n_gp_in_smallest_segment = iti->second;
            seg_id_smallest = it->first;
         }
      }
   }

   return seg_id_smallest;
}


// static
bool
coot::util::segment_map::sort_segment_vec(const std::pair<int, int> &a,
                                          const std::pair<int, int> &b) {
   return (b.first < a.first);
}


// change values in segmented_map
//
void
coot::util::segment_map::flood_fill_segmented_map(clipper::Xmap<std::pair<bool, int> > *segmented_map,
                                                  const clipper::Xmap<float> &xmap,
                                                  const clipper::Coord_grid &seed_point,
                                                  int from_val, int to_val) {

   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.5, 3.1); // 3x3x3 cube, not centre
   std::queue<clipper::Coord_grid> q;
   clipper::Coord_grid c_g;

   q.push(seed_point);
   int n_converted = 0;

   while (q.size()) {
      clipper::Coord_grid c_g_start = q.front(); // new local centre
      q.pop();
      for (int i_n=0; i_n<neighb.size(); i_n++) {
         c_g = c_g_start + neighb[i_n];
         if (segmented_map->get_data(c_g).second == from_val) {
            segmented_map->set_data(c_g, std::pair<bool,int> (1, to_val));
            if (0)
               std::cout << " converted " << c_g.format() << " from "
                         << from_val << " to " << to_val << std::endl;
            n_converted++; //debugging
            q.push(c_g);
         }
      }
   }

   if (0)
      std::cout << "=== flood_fill_segmented_map: converted "
                << n_converted << " grid points from " << from_val
                << " to " << to_val << std::endl;

   // debugging
   if (n_converted == 0) {
      std::cout << "diagnose 0 conversions: " << seed_point.format() << " "
                << segmented_map->get_data(seed_point).second << " with neighbours: " << std::endl;
      for (int i_n=0; i_n<neighb.size(); i_n++) {
         std::cout << "diagnose 0 conversions:    "
                   << i_n << " " << neighb[i_n].format() << " "
                   << segmented_map->get_data(seed_point + neighb[i_n]).second
                   << std::endl;
      }
   }

}

// Return a vector grid points as we trace the route to the
// local peak from start_point.
//
std::vector<clipper::Coord_grid>
coot::util::segment_map::path_to_peak(const clipper::Coord_grid &start_point,
                                      const clipper::Xmap<float> &xmap) {

   std::vector<clipper::Coord_grid> v;
   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.5, 3.1); // 3x3x3 cube, not centre

   clipper::Coord_grid current_point = start_point;
   v.push_back(start_point);
   bool at_peak = 0;

   while (! at_peak) {
      float current_val = xmap.get_data(current_point);

      // what is the biggest value of the neighbours greater than current_value?
      //
      bool found_bigger = 0;
      clipper::Coord_grid biggest_neighbour;
      float biggest = -1;
      for (int i_n=0; i_n<neighb.size(); i_n++) {
         float v = xmap.get_data(current_point + neighb[i_n]);
         if (0)
            std::cout << "current: " << current_point.format() << " " << current_val
                      << "    neighb " << i_n << " " << neighb[i_n].format() << " "
                      << v << " and biggest " << biggest << std::endl;
         if (v > biggest) {
            if (v > current_val) {
               found_bigger = 1;
               biggest = v;
               biggest_neighbour = neighb[i_n];
            }
         }
      }

      if (! found_bigger) {
         at_peak = 1;
      } else {
         current_point += biggest_neighbour;
         v.push_back(current_point);
      }

   }
   return v;
}


// pass a negative atom_selection_handle to build an atom map for the whole molecule
//
clipper::Xmap<float>
coot::util::calc_atom_map(mmdb::Manager *mol,
                          int atom_selection_handle,
                          const clipper::Cell &cell,
                          const clipper::Spacegroup &space_group,
                          const clipper::Grid_sampling &sampling) {

   clipper::Xmap<float> xmap;
   xmap.init(space_group, cell, sampling);

   std::vector<clipper::Atom> l;

   mmdb::PPAtom sel_atoms = 0;
   int n_atoms;
   mol->GetSelIndex(atom_selection_handle, sel_atoms, n_atoms);


   // 20130710 for some reason this gives (memory?) errors when making clipper::Atom atom
   // (coords are nan and ele is scrambled text)
   //
//    for (unsigned int i=0; i<n_atoms; i++) {
//       mmdb::Atom mmdb_atom;
//       mmdb_atom.Copy(sel_atoms[i]);
//       clipper::MMDBAtom clipper_mmdb_at(mmdb_atom);
//       clipper::Atom atom(clipper_mmdb_at);
//       std::cout << "copied mmdb_atom " << clipper_mmdb_at << " to clipper atom "
//                 << atom.coord_orth().format() << " " << atom.element()
//                 << std::endl;
//       l.push_back(atom);
//    }

   for (int iat=0; iat<n_atoms; iat++) {
      float rescale_b_u = 1/(8*M_PI*M_PI);
      mmdb::Atom *at = sel_atoms[iat];
      clipper::Coord_orth pt(at->x, at->y, at->z);
      std::string ele(at->element);
      clipper::Atom cat;
      cat.set_element(ele);
      cat.set_coord_orth(pt);
      cat.set_u_iso(at->tempFactor * rescale_b_u);
      cat.set_occupancy(at->occupancy);
      l.push_back(cat);
   }

   try {

      clipper::Atom_list al(l);
      if (0) {
         std::cout << "======================= al size(): " << al.size() << std::endl;
         std::cout << "in calc_atom_map() here are some atoms" << std::endl;
         for (unsigned int iat=0; iat<10; iat++)
            std::cout << "           "
                      << al[iat].coord_orth().format() << " u: "
                      << al[iat].u_iso() // u should be b/(8*pi*pi)
                      << " ele: "
                      << al[iat].element() << std::endl;
      }
      clipper::EDcalc_iso<float> e;
      e(xmap, al);
   }
   catch (const clipper::Message_generic &e) {
      std::cout << "ERROR:: some sort of clipper map problem" << std::endl;
      std::cout << e.text() << std::endl;
   }
   return xmap;
}


// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-exclusing CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
// 10: all-atom with b-factor-dependent radius
float
coot::util::map_to_model_correlation(mmdb::Manager *mol,
                                     const std::vector<residue_spec_t> &specs,
                                     const std::vector<residue_spec_t> &specs_for_masking_neighbs,
                                     unsigned short int atom_mask_mode,
                                     float atom_radius,
                                     const clipper::Xmap<float> &reference_map) {

   map_stats_t map_stats = coot::SIMPLE;
   density_correlation_stats_info_t dcs =
      map_to_model_correlation_stats(mol, specs, specs_for_masking_neighbs,
                                     atom_mask_mode, atom_radius, reference_map, map_stats);
   return dcs.correlation();

}


clipper::Xmap<float>
coot::util::mask_map(const clipper::Xmap<float> &xmap_in,
                    const std::vector<mmdb::Residue *> &neighb_residues) {

   clipper::Xmap<float> masked_map = xmap_in;

   for (unsigned int ir=0; ir<neighb_residues.size(); ir++) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      mmdb::Residue *residue_p = neighb_residues[ir];
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         clipper::Coord_orth pt = co(at);
         // float atom_radius = refmac_atom_radius(at);
         float atom_radius = 1.4;

         clipper::Coord_frac cf = pt.coord_frac(masked_map.cell());
         clipper::Coord_frac box0(
                                  cf.u() - atom_radius/masked_map.cell().descr().a(),
                                  cf.v() - atom_radius/masked_map.cell().descr().b(),
                                  cf.w() - atom_radius/masked_map.cell().descr().c());

         clipper::Coord_frac box1(
                                  cf.u() + atom_radius/masked_map.cell().descr().a(),
                                  cf.v() + atom_radius/masked_map.cell().descr().b(),
                                  cf.w() + atom_radius/masked_map.cell().descr().c());

         clipper::Grid_map grid(box0.coord_grid(masked_map.grid_sampling()),
                                box1.coord_grid(masked_map.grid_sampling()));

         float atom_radius_sq = atom_radius * atom_radius;

         clipper::Xmap_base::Map_reference_coord ix(masked_map, grid.min() ), iu, iv, iw;
         for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
            for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
               for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                  if ( (iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()) - pt).lengthsq() < atom_radius_sq) {
                     masked_map[iw] = -10;
                     if (false)
                        std::cout << "cutting into map at point "
                                  << iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()).format()
                                  << " for neighb atom at: " << pt.format() << " "
                                  << (iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()) - pt).lengthsq()
                                  << std::endl;
                  }
               }
            }
         }
      }
   }

   if (false) {
      clipper::CCP4MAPfile mapout;
      mapout.open_write("masked.map");
      mapout.export_xmap(masked_map);
      mapout.close_write();
   }

   return masked_map;
}



// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-exclusing CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
// 10: all-atom with b-factor-dependent radius
coot::util::density_correlation_stats_info_t
coot::util::map_to_model_correlation_stats(mmdb::Manager *mol,
                                           const std::vector<residue_spec_t> &specs,
                                           const std::vector<residue_spec_t> &specs_for_masking_neighbs,
                                           unsigned short int atom_mask_mode,
                                           float atom_radius_in,
                                           const clipper::Xmap<float> &reference_map,
                                           coot::map_stats_t map_stats_flag) {

   int ATOM_MASK_MAINCHAIN = 1;
   int ATOM_MASK_NOT_MAINCHAIN = 2;
   int ATOM_MASK_NOT_MAINCHAIN_OR_CB = 3;
   int ATOM_MASK_ALL_ATOM_B_FACTOR = 10;

   density_correlation_stats_info_t stats; // returned variable

   float atom_radius = atom_radius_in;
   float ret_val = -2;
   int SelHnd = mol->NewSelection(); // d
   bool debug = false;

   if (debug) {
      std::cout << "INFO:: map_to_model_correlation:: there are " << specs.size()
                << " residues " << std::endl;
      for (unsigned int ilocal=0; ilocal<specs.size(); ilocal++)
         std::cout << "   " << specs[ilocal] << std::endl;
   }

   for (unsigned int ilocal=0; ilocal<specs.size(); ilocal++) {

      std::string res_name_selection  = "*";
      std::string atom_name_selection = "*";

      if (atom_mask_mode != 0) { // main chain for standard amino acids
         mmdb::Residue *res = get_residue(specs[ilocal], mol);

         if (res) {
            std::string residue_name(res->GetResName());
            if (is_standard_residue_name(residue_name)) {

               // PDBv3 FIXME
               //
               if (atom_mask_mode == ATOM_MASK_MAINCHAIN)
                  atom_name_selection = " N  , H  , HA , CA , C  , O  ";
               if (atom_mask_mode == ATOM_MASK_NOT_MAINCHAIN)
                  atom_name_selection = "!( N  , H  , HA , CA , C  , O  )";
               if (atom_mask_mode == ATOM_MASK_NOT_MAINCHAIN_OR_CB)
                  atom_name_selection = "!( N  , H  , HA , CA , C  , O  , CB )";
            } else {
               if (atom_mask_mode == 4)
                  atom_name_selection = "%%%%%%"; // nothing (perhaps use "")
               if (atom_mask_mode == 5)
                  atom_name_selection = "%%%%%%"; // nothing
            }
         }
      }

      if (debug) {
         mmdb::Residue *res = get_residue(specs[ilocal], mol);
         std::cout << "debug during atom selection trying residue "
                   << specs[ilocal] << " got residue " << res
                   << std::endl;
      }

      mol->SelectAtoms(SelHnd, 1,
         specs[ilocal].chain_id.c_str(),
         specs[ilocal].res_no,
         specs[ilocal].ins_code.c_str(),
         specs[ilocal].res_no,
         specs[ilocal].ins_code.c_str(),
         res_name_selection.c_str(),
         atom_name_selection.c_str(),
         "*", // elements
         "*", // alt loc.
         mmdb::SKEY_OR
      );
   }

   std::vector<mmdb::Residue *> neighb_residues;
   for (unsigned int inb=0; inb<specs_for_masking_neighbs.size(); inb++) {
      mmdb::Residue *r = get_residue(specs_for_masking_neighbs[inb], mol);
      if (r)
         neighb_residues.push_back(r);
   }

   clipper::Xmap<float> calc_map =
   coot::util::calc_atom_map(mol, SelHnd,
                             reference_map.cell(),
                             reference_map.spacegroup(),
                             reference_map.grid_sampling());

      if (debug)
         std::cout << "DEBUG:: map_to_model_correlation_stats() calc_map null status: "
                   << calc_map.is_null() << std::endl;

      if (! (calc_map.is_null())) {
         clipper::Xmap<short int> masked_map(reference_map.spacegroup(),
         reference_map.cell(),
         reference_map.grid_sampling());
         clipper::Xmap_base::Map_reference_index ix;
         for (ix = masked_map.first(); !ix.last(); ix.next())
            masked_map[ix] = 0;
         mmdb::Atom **atom_selection = 0;
         int n_atoms;
         mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

         // debugging atom selection
         if (debug) {
            std::cout << "debug:: selected " << n_atoms << " atoms " << std::endl;
            for (int iat=0; iat<n_atoms; iat++)
               std::cout << "    " << iat << ": " << atom_selection[iat]->name << " "
                         << atom_spec_t(atom_selection[iat]) << std::endl;
         }

         if (n_atoms == 0)
            return stats;

         // atom selection grid:
         //
         //
         std::pair<clipper::Coord_orth, clipper::Coord_orth> selection_extents = util::extents(mol, specs);

         if (debug)
            std::cout << "INFO:: mol residue set extents: "
                      << selection_extents.first.format() << " to "
                      << selection_extents.second.format() << std::endl;

         // double border = 4.1;
         double border = 3.1; // border is used to create selection_grid.
         // the grid that we check at the end is +/-3 A of the X,Y,Z extends of the ligand

         bool good_frac_coords = false;  // force loop evaluation for the first time

         clipper::Coord_orth ex_pt_1_co;
         clipper::Coord_orth ex_pt_2_co;
         clipper::Coord_frac ex_pt_1_fc;
         clipper::Coord_frac ex_pt_2_fc;
         unsigned int n_rounds = 0;
         while (! good_frac_coords) {
            n_rounds++;
            ex_pt_1_co = clipper::Coord_orth(selection_extents.first.x()-border,
            selection_extents.first.y()-border,
            selection_extents.first.z()-border);
            ex_pt_2_co = clipper::Coord_orth(selection_extents.second.x()+border,
            selection_extents.second.y()+border,
            selection_extents.second.z()+border);
            ex_pt_1_fc = ex_pt_1_co.coord_frac(reference_map.cell());
            ex_pt_2_fc = ex_pt_2_co.coord_frac(reference_map.cell());
            if (false) { // debug
               std::cout << "INFO:: Selection grid construction, ex_pt_1_co: " << ex_pt_1_co.format() << std::endl;
               std::cout << "INFO:: Selection grid construction, ex_pt_2_co: " << ex_pt_2_co.format() << std::endl;
               std::cout << "INFO:: Selection grid construction, ex_pt_1_fc: " << ex_pt_1_fc.format() << std::endl;
               std::cout << "INFO:: Selection grid construction, ex_pt_2_fc: " << ex_pt_2_fc.format() << std::endl;
               std::cout << "using cell " << reference_map.cell().descr().format() << std::endl;
               clipper::Mat33<double> mat = reference_map.cell().matrix_frac();
               std::cout << "mat: \n" << mat.format() << std::endl;
            }
            good_frac_coords = true;
            for (int i=0; i<3; i++) {
               // Yes, this still can happen, but we can rescue (see later)
               // std::cout << "comparing: " << ex_pt_2_fc[i] << " should be greater than "
               // << ex_pt_1_fc[i] << std::endl;
               if (ex_pt_2_fc[i] < ex_pt_1_fc[i]) {
                  // std::cout << "opps false " << std::endl;
                  good_frac_coords = false;
                  border += 1.4;
               }
            }
            // Opps. We didn't converge - but we don't want to stay here.
            if (n_rounds >= 10)
            good_frac_coords = true;
         }

         // If needed, try to rescue using a different method - that works for P1. Maybe this should
         // replace the method above?
         //
         if (n_rounds == 10) {
            std::pair<clipper::Coord_frac, clipper::Coord_frac> nfc =
            find_struct_fragment_coord_fracs_v2(selection_extents, reference_map.cell());
            ex_pt_1_fc = nfc.first;
            ex_pt_2_fc = nfc.second;
         }

         clipper::Grid_map selection_grid(ex_pt_1_fc.coord_grid(reference_map.grid_sampling()),
         ex_pt_2_fc.coord_grid(reference_map.grid_sampling()));
         if (debug) {
            std::cout << "INFO:: Selection grid construction, ex_pt_1_co: " << ex_pt_1_co.format() << std::endl;
            std::cout << "INFO:: Selection grid construction, ex_pt_2_co: " << ex_pt_2_co.format() << std::endl;
            std::cout << "INFO:: Selection grid construction, ex_pt_1_fc: " << ex_pt_1_fc.format() << std::endl;
            std::cout << "INFO:: Selection grid construction, ex_pt_2_fc: " << ex_pt_2_fc.format() << std::endl;
            std::cout << "INFO:: Selection grid: " << selection_grid.format() << std::endl;
         }

         for (int iat=0; iat<n_atoms; iat++) {
            clipper::Coord_orth co(atom_selection[iat]->x,
               atom_selection[iat]->y,
               atom_selection[iat]->z);

               if (atom_mask_mode == ATOM_MASK_ALL_ATOM_B_FACTOR)
               atom_radius = refmac_atom_radius(atom_selection[iat]);

               clipper::Coord_frac cf = co.coord_frac(masked_map.cell());
               clipper::Coord_frac box0(
                  cf.u() - atom_radius/masked_map.cell().descr().a(),
                  cf.v() - atom_radius/masked_map.cell().descr().b(),
                  cf.w() - atom_radius/masked_map.cell().descr().c());

                  clipper::Coord_frac box1(
                     cf.u() + atom_radius/masked_map.cell().descr().a(),
                     cf.v() + atom_radius/masked_map.cell().descr().b(),
                     cf.w() + atom_radius/masked_map.cell().descr().c());

                     clipper::Grid_map grid(box0.coord_grid(masked_map.grid_sampling()),
                     box1.coord_grid(masked_map.grid_sampling()));

                     float atom_radius_sq = atom_radius * atom_radius;

                     clipper::Xmap_base::Map_reference_coord ix(masked_map, grid.min() ), iu, iv, iw;

                     if (debug) {
                        std::cout << "INFO:: masking iat " << iat << " box grid: " << grid.format() << std::endl;
                        std::cout << "INFO:: masking iu range: " << iu.coord().u() << " to " << grid.max().u() << std::endl;
                     }

                     for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
                        for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
                           for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                              if ( (iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()) - co).lengthsq() < atom_radius_sq) {
                                 if (0)
                                 std::cout << "masked point at "
                                 << iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()).format()
                                 << " centre point: " << co.format() << " "
                                 << (iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()) - co).lengthsq()
                                 << std::endl;
                                 masked_map[iw] = 1;
                              }
                           }
                        }
                     }
                  }

                  // masked map is 1 over the atoms of the selected residues (specs).
                  //
                  // Now we need to (potentially) cut into that near the atoms of neighb_residues.
                  //
                  for (unsigned int ir=0; ir<neighb_residues.size(); ir++) {
                     mmdb::PPAtom residue_atoms = 0;
                     int n_residue_atoms;
                     mmdb::Residue *residue_p = neighb_residues[ir];
                     residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                     for (int iat=0; iat<n_residue_atoms; iat++) {
                        mmdb::Atom *at = residue_atoms[iat];
                        clipper::Coord_orth pt = co(at);

                        if (atom_mask_mode == ATOM_MASK_ALL_ATOM_B_FACTOR)
                        atom_radius = refmac_atom_radius(at);
                        else
                        atom_radius = atom_radius_in;

                        clipper::Coord_frac cf = pt.coord_frac(masked_map.cell());
                        clipper::Coord_frac box0(
                           cf.u() - atom_radius/masked_map.cell().descr().a(),
                           cf.v() - atom_radius/masked_map.cell().descr().b(),
                           cf.w() - atom_radius/masked_map.cell().descr().c());

                           clipper::Coord_frac box1(
                              cf.u() + atom_radius/masked_map.cell().descr().a(),
                              cf.v() + atom_radius/masked_map.cell().descr().b(),
                              cf.w() + atom_radius/masked_map.cell().descr().c());

                              clipper::Grid_map grid(box0.coord_grid(masked_map.grid_sampling()),
                              box1.coord_grid(masked_map.grid_sampling()));

                              float atom_radius_sq = atom_radius * atom_radius;

                              clipper::Xmap_base::Map_reference_coord ix(masked_map, grid.min() ), iu, iv, iw;
                              for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
                                 for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
                                    for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                                       if ( (iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()) - pt).lengthsq() < atom_radius_sq) {
                                          if (masked_map[iw] == 1) {
                                             masked_map[iw] = 0;
                                             if (0)
                                             std::cout << "cutting into mask at point "
                                             << iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()).format()
                                             << " for neighb atom at: " << pt.format() << " "
                                             << (iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()) - pt).lengthsq()
                                             << std::endl;
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }


                        double sum_x  = 0;
                        double sum_y  = 0;
                        double sum_x_sqd  = 0;
                        double sum_y_sqd  = 0;
                        double sum_xy = 0;
                        double y;
                        double x;
                        int n = 0;

                        std::vector<double> map_samples; // for KS-test of flatness of difference map
                        bool debug_grid_points = true;

                        //       std::ofstream gp;
                        //       if (debug_grid_points) {
                        //          gp = std::ifstream("xy-points.tab");
                        //       }

                        // scan the selection grid
                        //
                        clipper::Xmap_base::Map_reference_coord iix(masked_map, selection_grid.min() ), iu, iv, iw;

                        if (debug) {
                           std::cout << "INFO:: correl " << " box grid: " << selection_grid.format() << std::endl;
                           std::cout << "INFO:: correl iu range: " << iix.coord().format() << " to " << selection_grid.max().u() << std::endl;
                        }
                        for (iu = iix; iu.coord().u() <= selection_grid.max().u(); iu.next_u() ) {
                           for (iv = iu; iv.coord().v() <= selection_grid.max().v(); iv.next_v() ) {
                              for (iw = iv; iw.coord().w() <= selection_grid.max().w(); iw.next_w() ) {
                                 if (masked_map[iw]) {
                                    x = calc_map[iw];
                                    if (! clipper::Util::is_nan(x)) {
                                       y = reference_map[iw];
                                       if (! clipper::Util::is_nan(y)) {
                                          if (0)
                                          std::cout << "xy-pair: " << x << " " << y << " "
                                          << iw.coord().u() << " "
                                          << iw.coord().v() << " "
                                          << iw.coord().w() << "\n";
                                          sum_x  += x;
                                          sum_y  += y;
                                          sum_xy += x * y;
                                          sum_x_sqd += x*x;
                                          sum_y_sqd += y*y;
                                          if (map_stats_flag == WITH_KOLMOGOROV_SMIRNOV_DIFFERENCE_MAP_TEST)
                                          map_samples.push_back(y);
                                          n++;
                                       } else {
                                          // std::cout << "null reference map data point at " << iw.coord().format() << std::endl;
                                       }
                                    } else {
                                       // std::cout << "null calc map data point at " << iw.coord().format() << std::endl;
                                    }
                                 }
                              }
                           }
                        }

                        //       if (debug_grid_points)
                        //          gp.close();

                        if (false) {
                           // just checking that the maps are what we expect them to be...
                           clipper::CCP4MAPfile mapout;
                           mapout.open_write("calc.map");
                           mapout.export_xmap(calc_map);
                           mapout.close_write();

                           clipper::CCP4MAPfile mapout_mask;
                           mapout_mask.open_write("masked.map");
                           mapout_mask.export_xmap(masked_map);
                           mapout_mask.close_write();
                        }

                        stats = density_correlation_stats_info_t(double(n), sum_xy,
                                                                 sum_x_sqd, sum_y_sqd,
                                                                 sum_x, sum_y);

                        if (map_stats_flag == WITH_KOLMOGOROV_SMIRNOV_DIFFERENCE_MAP_TEST)
                        stats.density_values = map_samples;

                        double top = double(n) * sum_xy - sum_x * sum_y;
                        double b_1 = double(n) * sum_x_sqd - sum_x * sum_x;
                        double b_2 = double(n) * sum_y_sqd - sum_y * sum_y;
                        if (debug) {
                           std::cout << ".... n is " << n << std::endl;
                           std::cout << ".... sum_xy is " << sum_xy << std::endl;
                           std::cout << ".... sum_x is " << sum_x << std::endl;
                           std::cout << ".... sum_y is " << sum_y << std::endl;
                           std::cout << ".... top is " << top << std::endl;
                           std::cout << ".... b_1 is " << b_1 << std::endl;
                           std::cout << ".... b_2 is " << b_2 << std::endl;
                        }

                        if (b_1 < 0) b_1 = 0;
                        if (b_2 < 0) b_2 = 0;

                        double c = top/(sqrt(b_1) * sqrt(b_2));
                        if (debug)
                        std::cout << "INFO:: map vs model correlation: "
                        << c << " vs " << stats.correlation() << std::endl;
                        ret_val = c;
                     }
                     mol->DeleteSelection(SelHnd);

   return stats;
}

// helper
std::pair<clipper::Coord_frac, clipper::Coord_frac>
coot::util::find_struct_fragment_coord_fracs_v2(const std::pair<clipper::Coord_orth, clipper::Coord_orth> &selection_extents,
                                                const clipper::Cell &cell) {

   clipper::Coord_orth sum = selection_extents.first + selection_extents.second;
   const clipper::Coord_orth &pt_1 = selection_extents.first;
   const clipper::Coord_orth &pt_2 = selection_extents.second;
   clipper::Coord_orth mid_pt(0.5 * sum.x(), 0.5 * sum.y(), 0.5 * sum.z());
   double l = 0.7 * clipper::Coord_orth::length(pt_2, pt_1);

   // we want to find the min and max fractional coords of the corner points of a box centred on mid_pt
   //
   std::vector<clipper::Coord_orth> v;
   v.push_back(clipper::Coord_orth(pt_1.x(), pt_1.y(), pt_1.z()));
   v.push_back(clipper::Coord_orth(pt_1.x(), pt_1.y(), pt_2.z()));
   v.push_back(clipper::Coord_orth(pt_1.x(), pt_2.y(), pt_1.z()));
   v.push_back(clipper::Coord_orth(pt_1.x(), pt_2.y(), pt_2.z()));
   v.push_back(clipper::Coord_orth(pt_2.x(), pt_1.y(), pt_1.z()));
   v.push_back(clipper::Coord_orth(pt_2.x(), pt_1.y(), pt_2.z()));
   v.push_back(clipper::Coord_orth(pt_2.x(), pt_2.y(), pt_1.z()));
   v.push_back(clipper::Coord_orth(pt_2.x(), pt_2.y(), pt_2.z()));

   clipper::Coord_frac cf_1( 99, 99,  99);
   clipper::Coord_frac cf_2(-99, -99, -99);
   for (std::size_t i=0; i<v.size(); i++) {
      const clipper::Coord_orth &pt = v[i];
      clipper::Coord_frac cf = pt.coord_frac(cell);
      if (cf.u() < cf_1.u()) cf_1 = clipper::Coord_frac(cf.u(),  cf_1.v(), cf_1.w());
      if (cf.v() < cf_1.v()) cf_1 = clipper::Coord_frac(cf_1.u(), cf.v(),  cf_1.w());
      if (cf.w() < cf_1.w()) cf_1 = clipper::Coord_frac(cf_1.u(), cf_1.v(),  cf.w());
      if (cf.u() > cf_2.u()) cf_2 = clipper::Coord_frac(cf.u(),  cf_2.v(), cf_2.w());
      if (cf.v() > cf_2.v()) cf_2 = clipper::Coord_frac(cf_2.u(), cf.v(),  cf_2.w());
      if (cf.w() > cf_2.w()) cf_2 = clipper::Coord_frac(cf_2.u(), cf_2.v(),  cf.w());
   }

   if (false) {
      std::cout << "......... rescue " << std::endl;
      std::cout << "     " << cf_1.format() << std::endl;
      std::cout << "     " << cf_2.format() << std::endl;
   }

   return std::pair<clipper::Coord_frac, clipper::Coord_frac> (cf_1, cf_2);

}


// the first of the pair contains the correlation for the given residue spec.
//
//! \brief atom-mask-mode is as follows:
// 0: all-atoms
// 1: main-chain atoms if is standard amino-acid, else all atoms
// 2: side-chain atoms if is standard amino-acid, else all atoms
// 3: side-chain atoms-excluding CB if is standard amino-acid, else all atoms
// 4: main-chain atoms if is standard amino-acid, else nothing
// 5: side-chain atoms if is standard amino-acid, else nothing
//
// reformat that:
//
// if (standard-amino-acid):
//    0: all atoms
//    1: main-chain
//    2: side-chain
//    3: side-chain excluding CB
//    4: main-chain atoms
//    5: side-chain atoms
// else
//    0: all atoms
//    1: all atoms
//    2: all atoms
//    3: all atoms
//    4: nothing
//    5: nothing
std::vector<std::pair<coot::residue_spec_t, float> >
coot::util::map_to_model_correlation_per_residue(mmdb::Manager *mol,
                                                 const std::vector<coot::residue_spec_t> &specs,
                                                 unsigned short int atom_mask_mode,
                                                 float atom_radius, // for masking
                                                 const clipper::Xmap<float> &reference_map) {

   if (false)
      std::cout << "DEBUG:: --------------- map_to_model_correlation_per_residue() "
                << specs.size() << std::endl;

   std::vector<std::pair<residue_spec_t, float> > v;
   int SelHnd = mol->NewSelection(); // d

   for (unsigned int ispec=0; ispec<specs.size(); ispec++) {

      std::string res_name_selection  = "*";
      std::string atom_name_selection = "*";

      if (atom_mask_mode != 0) { // main chain for standard amino acids
         mmdb::Residue *res = get_residue(specs[ispec], mol);
         if (res) {
            std::string residue_name(res->GetResName());
            if (is_standard_residue_name(residue_name)) {

               // PDBv3 FIXME
               //
               if (atom_mask_mode == 1 || atom_mask_mode == 4)
                  atom_name_selection = " N  , H  , HA , CA , C  , O  ";
               if (atom_mask_mode == 2 || atom_mask_mode == 5)
                  atom_name_selection = "!( N  , H  , HA , CA , C  , O  )";
               if (atom_mask_mode == 3)
                  atom_name_selection = "!( N  , H  , HA , CA , C  , O  , CB )";
            } else {
               if (atom_mask_mode == 4)
                  atom_name_selection = "%%%%%%"; // nothing (perhaps use "")
               if (atom_mask_mode == 5)
                  atom_name_selection = "%%%%%%"; // nothing
            }
         }
      }

      mol->SelectAtoms(SelHnd, 1,
                       specs[ispec].chain_id.c_str(),
                       specs[ispec].res_no,
                       specs[ispec].ins_code.c_str(),
                       specs[ispec].res_no,
                       specs[ispec].ins_code.c_str(),
                       res_name_selection.c_str(),
                       atom_name_selection.c_str(),
                       "*", // elements
                       "*", // alt loc.
                       mmdb::SKEY_OR
                       );
      if (0) { // debugging selection
         mmdb::PPAtom atom_selection = 0;
         int n_atoms;
         mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

         std::cout << "selected n_atoms " << n_atoms << " where specs[0] is " << specs[0]
                   << " for mask mode " << atom_mask_mode << std::endl;
         for (int iat=0; iat<n_atoms; iat++)
            std::cout << "    " << iat << " " << atom_spec_t(atom_selection[iat]) << std::endl;
      }
   }



   clipper::Xmap<float> calc_map =
      coot::util::calc_atom_map(mol, SelHnd,
                                reference_map.cell(),
                                reference_map.spacegroup(),
                                reference_map.grid_sampling());

   if (! (calc_map.is_null())) {
      // fill with null residue specs
      clipper::Xmap<residue_spec_t> contributor_map(reference_map.spacegroup(),
                                                    reference_map.cell(),
                                                    reference_map.grid_sampling());
      clipper::Xmap_base::Map_reference_index ix;
      mmdb::PPAtom atom_selection = 0;
      int n_atoms;
      mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

      residue_spec_t many_contributors; // we don't want map that is contributed to by
                                        // many (2 or more) residues.
      int MANY_CONTRIBUTORS = 200; // magic value
      many_contributors.int_user_data = MANY_CONTRIBUTORS;


      for (int iat=0; iat<n_atoms; iat++) {
         residue_spec_t res_spec(atom_selection[iat]->GetResidue());
         clipper::Coord_orth co(atom_selection[iat]->x,
                                atom_selection[iat]->y,
                                atom_selection[iat]->z);
         clipper::Coord_frac cf = co.coord_frac(reference_map.cell());
         clipper::Coord_frac box0(
                                  cf.u() - atom_radius/reference_map.cell().descr().a(),
                                  cf.v() - atom_radius/reference_map.cell().descr().b(),
                                  cf.w() - atom_radius/reference_map.cell().descr().c());

         clipper::Coord_frac box1(
                                  cf.u() + atom_radius/reference_map.cell().descr().a(),
                                  cf.v() + atom_radius/reference_map.cell().descr().b(),
                                  cf.w() + atom_radius/reference_map.cell().descr().c());

         clipper::Grid_map grid(box0.coord_grid(reference_map.grid_sampling()),
                                box1.coord_grid(reference_map.grid_sampling()));

         float atom_radius_sq = atom_radius * atom_radius;
         clipper::Xmap_base::Map_reference_coord ix(reference_map, grid.min() ), iu, iv, iw;
         for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
            for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
               for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                  if ( (iw.coord().coord_frac(reference_map.grid_sampling()).coord_orth(reference_map.cell()) - co).lengthsq() < atom_radius_sq) {
                     if (0)
                        std::cout << "specs: masked point at "
//                                   << iw.coord().coord_frac(reference_map.grid_sampling()).coord_orth(reference_map.cell()).format()
//                                   << " centre point: " << co.format() << " "
//                                   << (iw.coord().coord_frac(reference_map.grid_sampling()).coord_orth(reference_map.cell()) - co).lengthsq()
                                  << iw.coord().format()
                                  << std::endl;
                     if (contributor_map[iw].int_user_data != MANY_CONTRIBUTORS) {
                        // so it was either set once or not at all (so far).
                        if (contributor_map[iw].unset_p()) {
                           // was a default residue spec (an as-yet untouched grid point).
                           contributor_map[iw] = res_spec;
                        } else {
                           if (contributor_map[iw] == res_spec) {
                              // OK, cool we are in the same residue
                           } else {
                              // Oops. A bang.  Mask this point out
                              // std::cout << "   4 a bang for grid point " << iw.coord().format() << " "
                              // << contributor_map[iw] << " vs " << res_spec << std::endl;
                              contributor_map[iw] = many_contributors;
                           }
                        }
                     } else {
                        // std::cout << iw.coord().format() << " was already a many_contributor" << std::endl;
                     }
                  }
               }
            }
         }
      }

      std::map<residue_spec_t, map_stats_holder_helper_t> map_stats;
      for (ix = contributor_map.first(); !ix.last(); ix.next()) {
         if (! contributor_map[ix].unset_p()) {
            if (contributor_map[ix].int_user_data != MANY_CONTRIBUTORS) {
               // std::cout << "xy-pair " << calc_map[ix] << " " << reference_map[ix] << "\n";
               map_stats[contributor_map[ix]].add_xy(calc_map[ix], reference_map[ix]);
            }
         }
      }

      std::map<residue_spec_t, map_stats_holder_helper_t>::const_iterator it;
      for (it=map_stats.begin(); it!=map_stats.end(); ++it) {

         if (false)
            std::cout << "   " << it->first << " n: " << it->second.n << " spec user-data: "
                      << it->first.int_user_data << std::endl;

         if (it->second.n > 1) {
            double top = double(it->second.n) * it->second.sum_xy        - it->second.sum_x * it->second.sum_y;
            double b_1 = double(it->second.n) * it->second.sum_x_squared - it->second.sum_x * it->second.sum_x;
            double b_2 = double(it->second.n) * it->second.sum_y_squared - it->second.sum_y * it->second.sum_y;
            if (b_1 < 0) b_1 = 0;
            if (b_2 < 0) b_2 = 0;
            double c = top/(sqrt(b_1) * sqrt(b_2));
            std::pair<residue_spec_t, float> p(it->first, c);
            v.push_back(p);
         }
      }
   }
   mol->DeleteSelection(SelHnd);
   return v;
}

std::map<coot::residue_spec_t, coot::util::density_stats_info_t>
coot::util::map_to_model_correlation_stats_per_residue(mmdb::Manager *mol,
                                                       const std::vector<residue_spec_t> &specs,
                                                       unsigned short int atom_mask_mode,
                                                       float atom_radius, // for masking
                                                       const clipper::Xmap<float> &xmap) {

   std::map<residue_spec_t, density_stats_info_t> res_map;

   int SelHnd = mol->NewSelection(); // d

   for (unsigned int ispec=0; ispec<specs.size(); ispec++) {

      std::string res_name_selection  = "*";
      std::string atom_name_selection = "*";

      if (atom_mask_mode != 0) { // main chain for standard amino acids
         mmdb::Residue *res = get_residue(specs[ispec], mol);
         if (res) {
            std::string residue_name(res->GetResName());
            if (is_standard_residue_name(residue_name)) {

               // PDBv3 FIXME
               //
               if (atom_mask_mode == 1 || atom_mask_mode == 4)
               atom_name_selection = " N  , H  , HA , CA , C  , O  ";
               if (atom_mask_mode == 2 || atom_mask_mode == 5)
               atom_name_selection = "!( N  , H  , HA , CA , C  , O  )";
               if (atom_mask_mode == 3)
               atom_name_selection = "!( N  , H  , HA , CA , C  , O  , CB )";
            } else {
               if (atom_mask_mode == 4)
               atom_name_selection = "%%%%%%"; // nothing (perhaps use "")
               if (atom_mask_mode == 5)
               atom_name_selection = "%%%%%%"; // nothing
            }
         }
      }

      mol->SelectAtoms(SelHnd, 1,
                       specs[ispec].chain_id.c_str(),
                       specs[ispec].res_no,
                       specs[ispec].ins_code.c_str(),
                       specs[ispec].res_no,
                       specs[ispec].ins_code.c_str(),
                       res_name_selection.c_str(),
                       atom_name_selection.c_str(),
                       "*", // elements
                       "*", // alt loc.
                       mmdb::SKEY_OR
                       );

      if (false) { // debugging selection
         mmdb::PPAtom atom_selection = 0;
         int n_atoms;
         mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

         std::cout << "selected n_atoms " << n_atoms << " where specs[0] is " << specs[0]
                   << " for mask mode " << atom_mask_mode << std::endl;
         for (int iat=0; iat<n_atoms; iat++)
            std::cout << "    " << iat << " " << atom_spec_t(atom_selection[iat]) << std::endl;
      }
   }

   clipper::Xmap<float> calc_map =
      coot::util::calc_atom_map(mol, SelHnd, xmap.cell(), xmap.spacegroup(), xmap.grid_sampling());

   if (! (calc_map.is_null())) {
      // fill with null residue specs
      clipper::Xmap<residue_spec_t> contributor_map(xmap.spacegroup(),
                                                    xmap.cell(),
                                                    xmap.grid_sampling());
      clipper::Xmap_base::Map_reference_index ix;
      mmdb::PPAtom atom_selection = 0;
      int n_atoms;
      mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

      residue_spec_t many_contributors; // we don't want map that is contributed to by
                                        // many (2 or more) residues.
      int MANY_CONTRIBUTORS = 200; // magic value
      many_contributors.int_user_data = MANY_CONTRIBUTORS;


      for (int iat=0; iat<n_atoms; iat++) {
         residue_spec_t res_spec(atom_selection[iat]->GetResidue());
         clipper::Coord_orth co(atom_selection[iat]->x,
            atom_selection[iat]->y,
            atom_selection[iat]->z);
            clipper::Coord_frac cf = co.coord_frac(xmap.cell());
            clipper::Coord_frac box0(
               cf.u() - atom_radius/xmap.cell().descr().a(),
               cf.v() - atom_radius/xmap.cell().descr().b(),
               cf.w() - atom_radius/xmap.cell().descr().c());

         clipper::Coord_frac box1(
            cf.u() + atom_radius/xmap.cell().descr().a(),
            cf.v() + atom_radius/xmap.cell().descr().b(),
            cf.w() + atom_radius/xmap.cell().descr().c());

         clipper::Grid_map grid(box0.coord_grid(xmap.grid_sampling()),
                                box1.coord_grid(xmap.grid_sampling()));

         float atom_radius_sq = atom_radius * atom_radius;
         clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min()), iu, iv, iw;
         for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
            for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
               for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                  if ( (iw.coord().coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()) - co).lengthsq() < atom_radius_sq) {
                     if (false)
                        std::cout << "specs: masked point at "
                        //                                   << iw.coord().coord_frac(reference_map.grid_sampling()).coord_orth(reference_map.cell()).format()
                        //                                   << " centre point: " << co.format() << " "
                        //                                   << (iw.coord().coord_frac(reference_map.grid_sampling()).coord_orth(reference_map.cell()) - co).lengthsq()
                                  << iw.coord().format()
                                  << std::endl;
                     if (contributor_map[iw].int_user_data != MANY_CONTRIBUTORS) {
                        // so it was either set once or not at all (so far).
                        if (contributor_map[iw].unset_p()) {
                           // was a default residue spec (an as-yet untouched grid point).
                           contributor_map[iw] = res_spec;
                        } else {
                           if (contributor_map[iw] == res_spec) {
                              // OK, cool we are in the same residue
                           } else {
                              // Oops. A bang.  Mask this point out
                              // std::cout << "   4 a bang for grid point " << iw.coord().format() << " "
                              // << contributor_map[iw] << " vs " << res_spec << std::endl;
                              contributor_map[iw] = many_contributors;
                           }
                        }
                     } else {
                        // std::cout << iw.coord().format() << " was already a many_contributor" << std::endl;
                     }
                  }
               }
            }
         }
      }

      for (ix = contributor_map.first(); !ix.last(); ix.next()) {
         if (! contributor_map[ix].unset_p()) {
            if (contributor_map[ix].int_user_data != MANY_CONTRIBUTORS) {
               const residue_spec_t &res_spec = contributor_map[ix];
               res_map[res_spec].add(xmap[ix]);
            }
         }
      }
   }
   mol->DeleteSelection(SelHnd);

   return res_map;

}

#include "utils/split-indices.hh"

std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
          std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >
coot::util::map_to_model_correlation_stats_per_residue_run(mmdb::Manager *mol,
                                                           const std::string &chain_id,
                                                           const clipper::Xmap<float> &xmap,
                                                           unsigned int n_residues_per_blob,
                                                           bool exclude_CON,
                                                           float atom_mask_radius,
                                                           float NOC_mask_radius) {

   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> res_map_all_atom;
   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> res_map_side_chain;;

   // can be 5s also
   class residue_run_t {
   public:
      unsigned int idx_mid;
      unsigned int n_residues_per_blob;
      residue_run_t() {}
      explicit residue_run_t(unsigned int n_residues_per_blob_in) : n_residues_per_blob(n_residues_per_blob_in) {
         idx_mid = n_residues_per_blob/2;
      }
      explicit residue_run_t(const std::vector<mmdb::Residue *> &rr_in) : residues(rr_in) {
         n_residues_per_blob = rr_in.size();
         idx_mid = residues.size()/2;
      }
      std::vector<mmdb::Residue *> residues;
      void print() const {
         for (unsigned int i=0; i<residues.size(); i++) {
            std::cout << " " << residue_spec_t(residues[i]);
         }
         std::cout << std::endl;
      }
      void add_residue(mmdb::Residue *r) { residues.push_back(r); }
      mmdb::Residue *residue_mid() const {
         if (residues.size() >= n_residues_per_blob) {
            return residues[idx_mid];
         } else {
            std::cout << "ERROR:: indexing residues_mid " << residues.size() << std::endl;
            return 0;
         }
      }
      void add(const std::vector<mmdb::Residue *> &rv) {
         residues.insert(residues.begin(), rv.cbegin(), rv.cend());
      }
   };

   // Thi should be in utils, I guess.
   auto is_het_residue = [] (mmdb::Residue *residue_p) {
                            bool status = false;
                            mmdb::Atom **residue_atoms = 0;
                            int n_residue_atoms = 0;
                            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                            if (n_residue_atoms > 0) {
                               mmdb::Atom *at = residue_atoms[0];
                               if (! at->isTer())
                                  if (at->Het)
                                     status = true;
                            }
                            return status;
                         };

   // triples and 5s
   auto make_residue_runs = [mol, chain_id, n_residues_per_blob, is_het_residue] () {

      std::map<residue_spec_t, residue_run_t> residue_run_map;
      int n_models = mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int i_chain=0; i_chain<n_chains; i_chain++) {
               mmdb::Chain *chain_p = model_p->GetChain(i_chain);
               std::string this_chain_id(chain_p->GetChainID());
               if (this_chain_id == chain_id) {
                  int n_residues = chain_p->GetNumberOfResidues();
                  int idx_last_residue = n_residues - static_cast<int>(n_residues_per_blob);
                  std::vector<mmdb::Residue *> residue_vec;
                  for (int i_res=0; i_res<idx_last_residue; i_res++) {
                     for (int i_run=0; i_run<static_cast<int>(n_residues_per_blob); i_run++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(i_res+i_run);
                        // std::cout << "pushing back " << residue_spec_t(residue_p) << " to residue_vec" << std::endl;
                        std::string res_name(residue_p->GetResName());
                        if (res_name != "HOH" &&
                            res_name != "A"  && res_name != "G"  && res_name != "U"  && res_name != "C" &&
                            res_name != "DA" && res_name != "DG" && res_name != "DT" && res_name != "DC") {
                           if (! is_het_residue(residue_p)) {
                              residue_vec.push_back(residue_p);
                           }
                        }
                     }
                     if (residue_vec.size() >= n_residues_per_blob) {
                        if (! residue_vec.empty()) {
                           residue_spec_t spec(chain_p->GetResidue(i_res)->GetChainID(),
                                               chain_p->GetResidue(i_res)->GetSeqNum(),
                                               chain_p->GetResidue(i_res)->GetInsCode());
                           // if it's not in the map, we need to created it with a constructor that sets
                           // the mid point
                           std::map<residue_spec_t, residue_run_t>::const_iterator it;
                           it = residue_run_map.find(spec);
                           if (it == residue_run_map.end())
                              residue_run_map[spec] = residue_run_t(n_residues_per_blob);
                           residue_run_map[spec].add(residue_vec);
                           residue_vec.clear();
                        }
                     }
                  }
               }
            }
         }
      }
      if (true)
         std::cout << "debug:: make_residue_runs() residue_run_map size " << residue_run_map.size() << std::endl;
      std::vector<residue_run_t> residue_runs;
      std::map<residue_spec_t, residue_run_t>::const_iterator it;
      for (it=residue_run_map.begin(); it!=residue_run_map.end(); ++it) {
         residue_spec_t spec = it->first;
         const residue_run_t &rr(it->second);
         // rr.print();
         residue_runs.push_back(rr);
      }

      if (true)
         std::cout << "debug:: make_residue_runs() returns " << residue_runs.size() << " residue runs" << std::endl;
      return residue_runs;
   };

   auto get_residue_run_stats = [NOC_mask_radius] (const residue_run_t &residue_run,
                                                   const clipper::Xmap<float> &xmap,
                                                   const clipper::Xmap<float> &xmap_calc,
                                                   const clipper::Xmap<std::set<mmdb::Residue *> > &contributor_map,
                                                   float atom_radius,
                                                   bool exclude_NOC) {

                                   // if exclude_NOC is true, then fill dcs_side_chain also. Return the pair.

                                   util::density_correlation_stats_info_t dcs_all_atom;
                                   util::density_correlation_stats_info_t dcs_side_chain;
                                   std::vector<clipper::Coord_orth> atom_positions_to_avoid;
                                   if (exclude_NOC) {
                                      for (unsigned int ires=0; ires<residue_run.residues.size(); ires++) {
                                         mmdb::Residue *residue_p = residue_run.residues[ires];
                                         mmdb::Atom **residue_atoms = 0;
                                         int n_residue_atoms = 0;
                                         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                                         for (int iat=0; iat<n_residue_atoms; iat++) {
                                            mmdb::Atom *at = residue_atoms[iat];
                                            if (! at->isTer()) {
                                               std::string atom_name(at->GetAtomName());
                                               if (atom_name == " N  ") atom_positions_to_avoid.push_back(coot::co(at));
                                               if (atom_name == " C  ") atom_positions_to_avoid.push_back(coot::co(at));
                                               if (atom_name == " O  ") atom_positions_to_avoid.push_back(coot::co(at));
                                               if (atom_name == " H  ") atom_positions_to_avoid.push_back(coot::co(at));
                                            }
                                         }
                                      }
                                   }

                                   // std::cout << "rrrrrrrrrrrrrrrrrrrrrrrrr here with " << atom_positions_to_avoid.size()
                                   // << " atom positions to avoid "<< std::endl;

                                   for (unsigned int ires=0; ires<residue_run.residues.size(); ires++) {
                                      mmdb::Residue *residue_p = residue_run.residues[ires];
                                      mmdb::Atom **residue_atoms = 0;
                                      int n_residue_atoms = 0;
                                      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                                      for (int iat=0; iat<n_residue_atoms; iat++) {
                                         mmdb::Atom *at = residue_atoms[iat];
                                         if (! at->isTer()) {

                                            clipper::Coord_orth co(at->x, at->y, at->z);
                                            clipper::Coord_frac cf = co.coord_frac(contributor_map.cell());
                                            clipper::Coord_frac box0(cf.u() - atom_radius/contributor_map.cell().descr().a(),
                                                                     cf.v() - atom_radius/contributor_map.cell().descr().b(),
                                                                     cf.w() - atom_radius/contributor_map.cell().descr().c());
                                            clipper::Coord_frac box1(cf.u() + atom_radius/contributor_map.cell().descr().a(),
                                                                     cf.v() + atom_radius/contributor_map.cell().descr().b(),
                                                                     cf.w() + atom_radius/contributor_map.cell().descr().c());
                                            clipper::Grid_map grid(box0.coord_grid(contributor_map.grid_sampling()),
                                                                   box1.coord_grid(contributor_map.grid_sampling()));
                                            float atom_radius_sq = atom_radius * atom_radius;
                                            clipper::Xmap_base::Map_reference_coord ix(contributor_map, grid.min()), iu, iv, iw;
                                            for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
                                               for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
                                                  for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                                                     if ( (iw.coord().coord_frac(contributor_map.grid_sampling()).coord_orth(contributor_map.cell()) - co).lengthsq() < atom_radius_sq) {
                                                        // only count this grid point if it is not contributed to by residues other than those
                                                        // in the residue_run
                                                        bool found_other_residue_flag = false; // other than those in the residue_run, I mean
                                                        std::set<mmdb::Residue *>::const_iterator it;
                                                        for (it=contributor_map[iw].begin(); it!=contributor_map[iw].end(); ++it) {
                                                           const mmdb::Residue *contributor_residue = *it;
                                                           // is contributor_residue in the residue run?
                                                           if (std::find(residue_run.residues.begin(), residue_run.residues.end(), contributor_residue) == residue_run.residues.end()) {
                                                              // if it isn't, we can't count this grid point
                                                              found_other_residue_flag = true;
                                                              break;
                                                           }
                                                        }
                                                        if (!found_other_residue_flag) {
                                                           const float &xf = xmap[iw];
                                                           const float &yf = xmap_calc[iw];
                                                           if (! clipper::Util::is_nan(yf)) {
                                                              double x = static_cast<double>(xf);
                                                              double y = static_cast<double>(yf);
                                                              if (exclude_NOC) {
                                                                 bool found_a_close_NOC = false;
                                                                 const double dd_inside = NOC_mask_radius * NOC_mask_radius;
                                                                 for (const auto &NOC_pos : atom_positions_to_avoid) {
                                                                    double d = (co-NOC_pos).lengthsq();
                                                                    if (d < dd_inside) {
                                                                       found_a_close_NOC = true;
                                                                       break;
                                                                    }
                                                                 }
                                                                 if (found_a_close_NOC)
                                                                    dcs_all_atom.add(x,y); // just main-chain in this case (confusing name).
                                                                 else
                                                                    dcs_side_chain.add(x,y);
                                                              } else {
                                                                 dcs_all_atom.add(x,y);
                                                              }
                                                           } else {
                                                              if (false)
                                                                 std::cout << "Null calc map value " << atom_spec_t(at) << " " << iw.coord().format() << std::endl;
                                                           }
                                                        }
                                                     }
                                                  }
                                               }
                                            }
                                         }
                                      }
                                   }
                                   return std::pair<util::density_correlation_stats_info_t, util::density_correlation_stats_info_t>(dcs_all_atom, dcs_side_chain);
                                };

   auto fill_contributor_map = [] (clipper::Xmap<std::set<mmdb::Residue *> > &contributor_map, // reference
                                   mmdb::Manager *mol, int SelHnd, float atom_radius) {
                                  float atom_radius_sq = atom_radius * atom_radius;
                                  mmdb::PPAtom atom_selection = 0;
                                  int n_atoms;
                                  mol->GetSelIndex(SelHnd, atom_selection, n_atoms);
                                  for (int iat=0; iat<n_atoms; iat++) {
                                     mmdb::Atom *at = atom_selection[iat];
                                     mmdb::Residue *residue_p = at->residue;
                                     residue_spec_t res_spec(at->GetResidue());
                                     clipper::Coord_orth co(at->x, at->y, at->z);
                                     clipper::Coord_frac cf = co.coord_frac(contributor_map.cell());
                                     clipper::Coord_frac box0(cf.u() - atom_radius/contributor_map.cell().descr().a(),
                                                              cf.v() - atom_radius/contributor_map.cell().descr().b(),
                                                              cf.w() - atom_radius/contributor_map.cell().descr().c());

                                     clipper::Coord_frac box1(cf.u() + atom_radius/contributor_map.cell().descr().a(),
                                                              cf.v() + atom_radius/contributor_map.cell().descr().b(),
                                                              cf.w() + atom_radius/contributor_map.cell().descr().c());

                                     clipper::Grid_map grid(box0.coord_grid(contributor_map.grid_sampling()),
                                                            box1.coord_grid(contributor_map.grid_sampling()));
                                     clipper::Xmap_base::Map_reference_coord ix(contributor_map, grid.min()), iu, iv, iw;
                                     unsigned int i_count_for_this_atom = 0;
                                     for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
                                        for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
                                           for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                                              if ( (iw.coord().coord_frac(contributor_map.grid_sampling()).coord_orth(contributor_map.cell()) - co).lengthsq() < atom_radius_sq) {
                                                 contributor_map[iw].insert(residue_p);
                                                 i_count_for_this_atom ++;
                                                 // std::cout << "adding contributor " << residue_spec_t(residue_p) << " for " << iw.coord().format() << std::endl;
                                              }
                                           }
                                        }
                                     }
                                     // std::cout << "i_count_for_this_atom " << atom_spec_t(at) << " " << i_count_for_this_atom << std::endl;
                                  }
                               };

   auto multi_get_residue_range_stats = [get_residue_run_stats, exclude_CON] (unsigned int idx_start, unsigned int idx_end,
                                                                              const std::vector<residue_run_t> &residue_runs,
                                                                              const clipper::Xmap<float> &xmap,
                                                                              const clipper::Xmap<float> &calc_map,
                                                                              const clipper::Xmap<std::set<mmdb::Residue *> > &contributor_map,
                                                                              float atom_radius,
                                                                              std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> &res_map_all_atom_l,
                                                                              std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> &res_map_side_chain_l) {

                                           for (unsigned int i=idx_start; i<idx_end; i++) {
                                              const residue_run_t &residue_run = residue_runs[i];
                                              std::pair<util::density_correlation_stats_info_t, util::density_correlation_stats_info_t> stats_all_atom_and_side_chains =
                                                 get_residue_run_stats(residue_run, xmap, calc_map, contributor_map, atom_radius, exclude_CON);
                                              const util::density_correlation_stats_info_t &stats_all_atom    = stats_all_atom_and_side_chains.first;
                                              const util::density_correlation_stats_info_t &stats_side_chain  = stats_all_atom_and_side_chains.second;
                                              res_map_all_atom_l[residue_spec_t(residue_run.residue_mid())]   = stats_all_atom;
                                              res_map_side_chain_l[residue_spec_t(residue_run.residue_mid())] = stats_side_chain;
                                           }
                                        };


   std::vector<residue_run_t> residue_runs = make_residue_runs();

   if (false)
      for (const auto &rr : residue_runs)
         rr.print();

   int SelHnd = mol->NewSelection(); // d
   mol->SelectAtoms(SelHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
   clipper::Xmap<float> calc_map = coot::util::calc_atom_map(mol, SelHnd, xmap.cell(), xmap.spacegroup(), xmap.grid_sampling());

   if (false) {
      clipper::CCP4MAPfile mapout;
      mapout.open_write("calc-map.map");
      mapout.export_xmap(calc_map);
      mapout.close_write();
   }

   if (calc_map.is_null()) {
      std::cout << "OOPS! calc_map is null" << std::endl;
   } else {
      clipper::Xmap<std::set<mmdb::Residue *> > contributor_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
      fill_contributor_map(std::ref(contributor_map), mol, SelHnd, atom_mask_radius);

#if 0 // single thread
      for (unsigned int i=0; i<residue_runs.size(); i++) {
         const residue_run_t &residue_run = residue_runs[i];
         util::density_correlation_stats_info_t stats = get_residue_run_stats(residue_run, xmap, calc_map, contributor_map, atom_radius);
         if (false)
            std::cout << "debug:: residue-run stats for mid-res " << residue_spec_t(residue_run.residue_mid()) << " n: "
                      << stats.n << std::endl;
         if (stats.n > 1)
            res_map[residue_spec_t(residue_run.residue_mid())] = stats;
      }
#endif

#if 1
      if (! residue_runs.empty()) {
         unsigned int n_threads = coot::get_max_number_of_threads();
         std::vector<std::thread> threads;
         std::vector<std::pair<unsigned int, unsigned int> > ranges = atom_index_ranges(residue_runs.size(), n_threads);

         std::vector<std::map<residue_spec_t, util::density_correlation_stats_info_t> > stats_map_all_atom_vec(n_threads);   // local maps
         std::vector<std::map<residue_spec_t, util::density_correlation_stats_info_t> > stats_map_side_chain_vec(n_threads); // local maps
         for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {
            std::map<coot::residue_spec_t, util::density_correlation_stats_info_t> &res_map_all_atom_l   = stats_map_all_atom_vec[i_thread];
            std::map<coot::residue_spec_t, util::density_correlation_stats_info_t> &res_map_side_chain_l = stats_map_side_chain_vec[i_thread];
            if (i_thread < ranges.size())
               threads.push_back(std::thread(multi_get_residue_range_stats, ranges[i_thread].first, ranges[i_thread].second,
                                             std::cref(residue_runs), std::cref(xmap), std::cref(calc_map), std::cref(contributor_map),
                                             atom_mask_radius, std::ref(res_map_all_atom_l), std::ref(res_map_side_chain_l)));
            else {
               logger.log(log_t::ERROR, "bad thread index", i_thread, "vs", ranges.size());
            }
         }

         for (unsigned int i_thread=0; i_thread<threads.size(); i_thread++)
            threads[i_thread].join();

         // now merge the maps
         for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {
            std::map<residue_spec_t, util::density_correlation_stats_info_t> &res_map_l = stats_map_all_atom_vec[i_thread];
            std::map<residue_spec_t, util::density_correlation_stats_info_t>::const_iterator it;
            for (it=res_map_l.begin(); it!=res_map_l.end(); ++it)
               res_map_all_atom[it->first] = it->second; // all the keys are different. But are they for multple models
            res_map_l = stats_map_side_chain_vec[i_thread];
            for (it=res_map_l.begin(); it!=res_map_l.end(); ++it)
               res_map_side_chain[it->first] = it->second; // all the keys are different. But are they for multple models
         }
      }
#endif
   }

   mol->DeleteSelection(SelHnd);
   return std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
                    std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> >(res_map_all_atom, res_map_side_chain);
}





// should this be here, or is it heavy?
//
// Input maps is a difference map, return a vector of doubles for a plot
//
std::vector<std::pair<double, double> >
coot::util::qq_plot_for_map_over_model(mmdb::Manager *mol,
                                       const std::vector<coot::residue_spec_t> &specs,
                                       const std::vector<coot::residue_spec_t> &nb_residues,
                                       int atom_mask_mode,
                                       const clipper::Xmap<float> &xmap) {

   // We need to get the grid points of the input (difference) map.
   // To do so, this is like the correlation operation.

   // First, identify parts of the map which correspond to the atom
   // selection (that are not covered by the neighbouring reisdues.

   std::vector<mmdb::Residue *> neighb_residues;
   for (unsigned int i=0; i<nb_residues.size(); i++) {
      mmdb::Residue *r = get_residue(nb_residues[i], mol);
      if (r)
         neighb_residues.push_back(r);
   }

   // We must delete the selection!
   //
   int SelHnd = specs_to_atom_selection(specs, mol, atom_mask_mode); // d
   mmdb::PPAtom sel_atoms = 0;
   int n_atoms;
   mol->GetSelIndex(SelHnd, sel_atoms, n_atoms);

   clipper::Xmap<short int> mask(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap<short int>::Map_reference_index inx;
   for (inx = mask.first(); !inx.last(); inx.next())
      mask[inx] = 0;

   int n_points_masked = 0;
   for (int iat=0; iat<n_atoms; iat++) {
      mmdb::Atom *at = sel_atoms[iat];
      clipper::Coord_orth c_o = co(at);
      float radius = 1.5 + at->tempFactor*1.5/80.0; // should be some function of tempFactor;
      float radius_sq = radius * radius;

      clipper::Coord_frac cf = c_o.coord_frac(xmap.cell());
      clipper::Coord_frac box0(cf.u() - radius/xmap.cell().descr().a(),
                               cf.v() - radius/xmap.cell().descr().b(),
                               cf.w() - radius/xmap.cell().descr().c());

      clipper::Coord_frac box1(cf.u() + radius/xmap.cell().descr().a(),
                               cf.v() + radius/xmap.cell().descr().b(),
                               cf.w() + radius/xmap.cell().descr().c());

      clipper::Grid_map grid(box0.coord_grid(xmap.grid_sampling()),
                             box1.coord_grid(xmap.grid_sampling()));

      clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min()), iu, iv, iw;
      for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
         for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
               // sample a sphere
               if ( (iw.coord().coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()) - c_o).lengthsq() < radius_sq) {
                  mask[iw] = 1;
                  n_points_masked++;
               }
            }
         }
      }
   }

   mol->DeleteSelection(SelHnd);

   std::vector<double> map_points;
   for (inx = mask.first(); !inx.last(); inx.next()) {
      if (mask[inx]) {
         float v = xmap[inx];
         map_points.push_back(v);
      }
   }

   std::cout << "map_points.size(): " << map_points.size() << " n_points_masked "
             << n_points_masked << std::endl;

   qq_plot_t qq(map_points);
   return qq.qq_norm();
}


std::vector<float>
coot::util::density_map_points_in_sphere(clipper::Coord_orth pt, float search_radius, const clipper::Xmap<float> &xmap) {

   std::vector<float> v;

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

   clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min()), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            // sample a sphere
            if ( (iw.coord().coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()) - pt).lengthsq() < search_radius_sq) {
               v.push_back(xmap[iw]);
            }
         }
      }
   }
   return v;
}


coot::util::map_fragment_info_t::map_fragment_info_t(const clipper::Xmap<float> &ip_xmap,
                                                     const clipper::Coord_orth &centre,
                                                     float radius, bool centre_at_origin) {

   if (centre_at_origin)
      init_making_map_centred_at_origin(ip_xmap, centre, radius);
   else
      init(ip_xmap, centre, radius);
}


void
coot::util::map_fragment_info_t::init(const clipper::Xmap<float> &ip_xmap,
                                      const clipper::Coord_orth &centre,
                                      float radius) {

   clipper::Cell          ip_xmap_cell = ip_xmap.cell();
   clipper::Grid_sampling ip_xmap_grid_sampling = ip_xmap.grid_sampling();

   clipper::Grid_range gr0(ip_xmap_cell, ip_xmap_grid_sampling, radius);
   clipper::Grid_range gr1(gr0.min() + centre.coord_frac(ip_xmap_cell).coord_grid(ip_xmap_grid_sampling),
                           gr0.max() + centre.coord_frac(ip_xmap_cell).coord_grid(ip_xmap_grid_sampling));

   int new_x_u = 2*gr0.max().u();
   int new_x_v = 2*gr0.max().v();
   int new_x_w = 2*gr0.max().w();

   clipper::Coord_grid new_xmap_grid(new_x_u, new_x_v, new_x_w);
   clipper::Coord_grid new_xmap_grid_max(new_x_u-1, new_x_v-1, new_x_w-1);
   clipper::Coord_grid new_xmap_origin(0,0,0);
   clipper::Grid_range new_xmap_grid_range(new_xmap_origin, new_xmap_grid_max);
   clipper::Grid_sampling new_xmap_grid_sampling(new_x_u, new_x_v, new_x_w);
   std::cout << "--------------- new_xmap grid_sampling init with " << new_x_u << " " << new_x_v << " " << new_x_w
             << std::endl;

   std::cout << "--------------- gr0.min() " << gr0.min().format() << std::endl;
   std::cout << "--------------- gr0.max() " << gr0.max().format() << std::endl;
   std::cout << "--------------- new_xmap_grid_sampling "   << new_xmap_grid_sampling.format() << std::endl;
   std::cout << "--------------- new_xmap grid range min: " << new_xmap_grid_range.min().format() << std::endl;
   std::cout << "--------------- new_xmap grid range max: " << new_xmap_grid_range.max().format() << std::endl;

   std::cout << "DEBUG:: input map cells/A: "
             << ip_xmap_grid_sampling.nu()/ip_xmap.cell().descr().a() << " "
             << ip_xmap_grid_sampling.nv()/ip_xmap.cell().descr().b() << " "
             << ip_xmap_grid_sampling.nw()/ip_xmap.cell().descr().c() << " "
             << std::endl;

   double new_x_alpha = M_PI_2;
   double new_x_beta  = M_PI_2;
   double new_x_gamma = M_PI_2;
   double new_x_a = ip_xmap.cell().descr().a() * double(new_xmap_grid_sampling.nu())/double(ip_xmap_grid_sampling.nu());
   double new_x_b = ip_xmap.cell().descr().b() * double(new_xmap_grid_sampling.nv())/double(ip_xmap_grid_sampling.nv());
   double new_x_c = ip_xmap.cell().descr().c() * double(new_xmap_grid_sampling.nw())/double(ip_xmap_grid_sampling.nw());
   clipper::Cell_descr new_xmap_cell_descr(new_x_a, new_x_b, new_x_c, new_x_alpha, new_x_beta, new_x_gamma);
   clipper::Cell new_xmap_cell(new_xmap_cell_descr);

   // init new_xmap
   // clipper::Xmap<float> new_xmap(clipper::Spacegroup::p1(), new_xmap_cell, new_xmap_grid_sampling);
   xmap.init(clipper::Spacegroup::p1(), new_xmap_cell, new_xmap_grid_sampling);

   clipper::Xmap<float>::Map_reference_coord ix(ip_xmap);
   clipper::Coord_orth centre_radius(centre.x() - radius,
                                     centre.y() - radius,
                                     centre.z() - radius);
   clipper::Coord_orth new_xmap_centre(radius, radius, radius);
   offset = ip_xmap.coord_map(centre_radius).coord_grid();

   std::cout << "--------------- xmap offset to centre        " << centre_radius.format() << std::endl;
   std::cout << "--------------- xmap offset to centre (grid) " << offset.format() << std::endl;

   clipper::Xmap<float>::Map_reference_index inx;
   double limited_radius = radius * 0.92;
   for (inx = xmap.first(); !inx.last(); inx.next()) {
      clipper::Coord_orth p = inx.coord().coord_frac(new_xmap_grid_sampling).coord_orth(new_xmap_cell);
      double d_to_c_sq = clipper::Coord_orth(p-new_xmap_centre).lengthsq();
      if (d_to_c_sq > limited_radius*limited_radius) {
         xmap[inx] = 0.0;
      } else {
         // happy path

         // ix indexes the input xmap
         //
         // 20140214: I'm retiring the grid based method of sampling the input map (in part because
         // now the output cell is orthogonal).
         //
         // clipper::Coord_grid gp = p.coord_frac(ip_xmap_cell).coord_grid(ip_xmap_grid_sampling);
         // ix.set_coord(gp + offset);

         float dv = density_at_point(ip_xmap, p+centre);

         // make a function that is
         // y=1   around x=0 and
         // y=0.5 around x=0.8 or so.
         // y=0   around x=1
         //
         double x = sqrt(d_to_c_sq)/limited_radius;
         double gompertz_a = 0.14;
         double gompertz_b = 0.1;
         double gompertz_c = 3;
         double gompertz_scale = 1 - (-gompertz_a*1.1 + gompertz_a * exp (gompertz_b * exp(gompertz_c * x)));

         xmap[inx] = dv * gompertz_scale;
         if (0)
            std::cout << " inx " << inx.coord().format() << " " << d_to_c_sq << "  " << p.format() << " "
                      << centre.format() << " vs " << limited_radius*limited_radius << " " << gompertz_scale
                      << std::endl;
      }
   }
}

#include "utils/radix.hh"


clipper::Grid_map
coot::util::map_fragment_info_t::make_grid_map(const clipper::Xmap<float> &input_xmap,
                                               const clipper::Coord_orth &centre_from_input) const {

   clipper::Coord_frac centre_from_input_f    = centre_from_input.coord_frac(input_xmap.cell());
   clipper::Coord_map  centre_from_input_m    = centre_from_input_f.coord_map(input_xmap.grid_sampling());
   clipper::Coord_grid centre_from_input_g    = centre_from_input_m.coord_grid();
   clipper::Coord_map  centre_at_grid_point_m = centre_from_input_g.coord_map();
   clipper::Coord_frac centre_at_grid_point_f = centre_at_grid_point_m.coord_frac(input_xmap.grid_sampling());

   clipper::Coord_frac box0(centre_at_grid_point_f.u() - box_radius_a_internal/input_xmap.cell().descr().a(),
                            centre_at_grid_point_f.v() - box_radius_b_internal/input_xmap.cell().descr().b(),
                            centre_at_grid_point_f.w() - box_radius_c_internal/input_xmap.cell().descr().c());
   clipper::Coord_frac box1(centre_at_grid_point_f.u() + box_radius_a_internal/input_xmap.cell().descr().a(),
                            centre_at_grid_point_f.v() + box_radius_b_internal/input_xmap.cell().descr().b(),
                            centre_at_grid_point_f.w() + box_radius_c_internal/input_xmap.cell().descr().c());

   clipper::Grid_sampling input_gs = input_xmap.grid_sampling();
   clipper::Grid_map grid(box0.coord_grid(input_gs), box1.coord_grid(input_gs));

   return grid;
}


void
coot::util::map_fragment_info_t::simple_origin_shift(const clipper::Xmap<float> &input_xmap,
                                                     const clipper::Coord_orth &centre_from_input,
                                                     float box_radius) {

   // -----------------------------------------------------------------------------------------
   // -----------------------------------------------------------------------------------------
   //                                  currently broken
   // -----------------------------------------------------------------------------------------
   // -----------------------------------------------------------------------------------------
   // and difficult to fix.

   // I had intended to use this for interactive local map sharpening.
   // i.e. move a small map fragment to the origin,
   //      sharpen/blur that map
   //      transfer those grid coordinates back to where they came from
   // but the function crashes. Currently there is a test function (now disabled) in c-interface-test.cc
   // To try this again, you should disable the current constructor (or create a new one that does nothing).


   // the centre of the box should be a grid point, not the user value
   clipper::Coord_frac centre_from_input_f = centre_from_input.coord_frac(input_xmap.cell());
   clipper::Coord_map centre_from_input_m = centre_from_input_f.coord_map(input_xmap.grid_sampling());
   clipper::Coord_grid centre_from_input_g = centre_from_input_m.coord_grid();
   clipper::Coord_map centre_at_grid_point_m = centre_from_input_g.coord_map();
   clipper::Coord_frac centre_at_grid_point_f = centre_at_grid_point_m.coord_frac(input_xmap.grid_sampling());
   clipper::Coord_orth centre_at_grid_point = centre_at_grid_point_f.coord_orth(input_xmap.cell());

   std::cout << "################## simple_origin_shift() " << centre_from_input_g.format() << "\n";
   std::cout << "################## simple_origin_shift() " << centre_from_input.format() << " moved to "
             << centre_at_grid_point.format() << std::endl;

   clipper::Coord_orth centre = centre_at_grid_point;

   clipper::Coord_frac centre_f = centre.coord_frac(input_xmap.cell());
   clipper::Coord_frac input_box0(centre_f.u() - box_radius/input_xmap.cell().descr().a(),
                                  centre_f.v() - box_radius/input_xmap.cell().descr().b(),
                                  centre_f.w() - box_radius/input_xmap.cell().descr().c());
   clipper::Coord_frac input_box1(centre_f.u() + box_radius/input_xmap.cell().descr().a(),
                                  centre_f.v() + box_radius/input_xmap.cell().descr().b(),
                                  centre_f.w() + box_radius/input_xmap.cell().descr().c());
   clipper::Grid_sampling input_gs = input_xmap.grid_sampling();
   clipper::Grid_map input_grid(input_box0.coord_grid(input_gs), input_box1.coord_grid(input_gs));

   double alpha = M_PI_2;
   double beta  = M_PI_2;
   double gamma = M_PI_2;
   clipper::Cell_descr cd(2.0 * box_radius, 2.0 * box_radius, 2.0 * box_radius, alpha, beta, gamma);
   clipper::Cell_descr cd_input = input_xmap.cell();
   int nnx = static_cast<int>(static_cast<float>(input_gs.nu()) * 2.0 * box_radius / cd_input.a());
   int nny = static_cast<int>(static_cast<float>(input_gs.nv()) * 2.0 * box_radius / cd_input.b());
   int nnz = static_cast<int>(static_cast<float>(input_gs.nw()) * 2.0 * box_radius / cd_input.c());

   int nnx_radix = coot::suggest_radix(nnx);
   int nny_radix = coot::suggest_radix(nny);
   int nnz_radix = coot::suggest_radix(nnz);

   std::cout << "debug orig  nnx nny nnz " << nnx << " " << nny << " " << nnz << std::endl;
   std::cout << "debug radix nnx nny nnz " << nnx_radix << " " << nny_radix << " " << nnz_radix << std::endl;
   double small = 1e-5;
   clipper::Grid_sampling gs_new(nnx_radix, nny_radix, nnz_radix);
   double new_box_a = box_radius * static_cast<double>(nnx_radix)/static_cast<double>(nnx) + small;
   double new_box_b = box_radius * static_cast<double>(nny_radix)/static_cast<double>(nny) + small;
   double new_box_c = box_radius * static_cast<double>(nnz_radix)/static_cast<double>(nnz) + small;
   double new_a = 2.0 * new_box_a;
   double new_b = 2.0 * new_box_b;
   double new_c = 2.0 * new_box_c;
   box_radius_a_internal = new_box_a;
   box_radius_b_internal = new_box_b;
   box_radius_c_internal = new_box_c;
   clipper::Cell_descr cd_new(new_a, new_b, new_c, alpha, beta, gamma);
   clipper::Cell new_xmap_cell(cd_new);
   std::cout << "old grid sampling " << input_gs.format() << std::endl;
   std::cout << "new grid sampling " << gs_new.format() << std::endl;
   std::cout << "old cell " << input_xmap.cell().format() << std::endl;
   std::cout << "new cell " << new_xmap_cell.format() << std::endl;
   xmap.init(clipper::Spacegroup::p1(), new_xmap_cell, gs_new);
   clipper::Xmap<float>::Map_reference_index inx;
   for (inx = xmap.first(); !inx.last(); inx.next()) xmap[inx]= 0.0;

   // offset brings the bottom left hand side of the grid in the input map to the origin
   clipper::Coord_orth offset_o = centre - clipper::Coord_orth(new_box_a, new_box_b, new_box_c);
   clipper::Coord_frac offset_f = offset_o.coord_frac(input_xmap.cell());
   clipper::Coord_map offset_m = input_xmap.coord_map(offset_o);
   offset = offset_m.coord_grid();
   std::cout << "offset o " << offset_o.format() << std::endl;
   std::cout << "offset f " << offset_f.format() << std::endl;
   std::cout << "offset m " << offset.format() << std::endl;

   clipper::Coord_frac box0(centre_at_grid_point_f.u() - new_box_a/input_xmap.cell().descr().a(),
                            centre_at_grid_point_f.v() - new_box_b/input_xmap.cell().descr().b(),
                            centre_at_grid_point_f.w() - new_box_c/input_xmap.cell().descr().c());
   clipper::Coord_frac box1(centre_at_grid_point_f.u() + new_box_a/input_xmap.cell().descr().a(),
                            centre_at_grid_point_f.v() + new_box_b/input_xmap.cell().descr().b(),
                            centre_at_grid_point_f.w() + new_box_c/input_xmap.cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(input_gs), box1.coord_grid(input_gs));
   std::cout << "DEUBG:: grid in reference map: " << grid.format() << std::endl;

   clipper::Xmap_base::Map_reference_coord ix(input_xmap);
   for (int u = grid.min().u(); u < grid.max().u(); u++) {
      for (int v = grid.min().v(); v < grid.max().v(); v++) {
         for (int w = grid.min().w(); w < grid.max().w(); w++) {
            ix.set_coord(clipper::Coord_grid(u, v, w)); // don't copy this. using set_coord() is slow
            float f = input_xmap[ix];
            clipper::Coord_grid cg = ix.coord() - offset;
            xmap.set_data(cg, f);
            if (false)
              std::cout << "set xmap: from " << ix.coord().format() << " " << cg.format() << " " << f << std::endl;
         }
      }
   }

   if (false) {
      clipper::Xmap_base::Map_reference_coord ix(input_xmap);
      for (int u = grid.min().u(); u < grid.max().u(); u++) {
         for (int v = grid.min().v(); v < grid.max().v(); v++) {
            for (int w = grid.min().w(); w < grid.max().w(); w++) {
               ix.set_coord(clipper::Coord_grid(u, v, w));
               clipper::Coord_grid cg = ix.coord() - offset;
               float f1 = input_xmap[ix];
               float f2 = xmap.get_data(cg);
               std::cout << "DEBUG:: test-get from " << ix.coord().format() << " to " << cg.format() << " "
                         << f1 << " " << f2 << std::endl;
            }
         }
      }
   }

}


void
coot::util::map_fragment_info_t::init_making_map_centred_at_origin(const clipper::Xmap<float> &ip_xmap,
                                                                   const clipper::Coord_orth &centre,
                                                                   float radius) {

   std::cout << "------------------------ centre returned map at origin ---------------------------"
             << std::endl;

   clipper::Grid_sampling ip_xmap_grid_sampling = ip_xmap.grid_sampling();
   double gpa = ip_xmap.grid_sampling().nu()/ip_xmap.cell().a();
   int n = int(radius*gpa)+1;

   // this needs to take into account that only a fraction of the cell is sampled. Otherwise
   // we supersample the map. Calculate this after new_xmap_a,b,c have been calculated.
   // (i.e. we have the sampe number of grid points in a much smaller box)
   clipper::Grid_sampling new_xmap_grid_sampling = ip_xmap.grid_sampling();

   double border = 10;
   double new_xmap_alpha = M_PI_2;
   double new_xmap_beta  = M_PI_2;
   double new_xmap_gamma = M_PI_2;
   double new_xmap_a = 2 * radius * double(new_xmap_grid_sampling.nu())/double(ip_xmap_grid_sampling.nu()) + border;
   double new_xmap_b = 2 * radius * double(new_xmap_grid_sampling.nv())/double(ip_xmap_grid_sampling.nv()) + border;
   double new_xmap_c = 2 * radius * double(new_xmap_grid_sampling.nw())/double(ip_xmap_grid_sampling.nw()) + border;
   clipper::Cell_descr new_cell_descr(new_xmap_a, new_xmap_b, new_xmap_c,
                                      new_xmap_alpha, new_xmap_beta, new_xmap_gamma);
   clipper::Cell new_xmap_cell(new_cell_descr);

   clipper::Coord_grid grid_min(-n,-n,-n);
   clipper::Coord_grid grid_max(n,n,n);
   clipper::Grid_range gr(grid_min, grid_max);

   std::cout << "--------------- ip centre:  " << centre.format() << std::endl;
   std::cout << "--------------- ip radius   " << radius << std::endl;
   std::cout << "--------------- new_xmap map gpa  " << gpa << std::endl;
   std::cout << "--------------- new_xmap cell     " << new_xmap_cell.format() << std::endl;
   std::cout << "--------------- new_xmap map n    " << n << std::endl;
   std::cout << "--------------- new_xmap map grid min " << grid_min.format() << std::endl;
   std::cout << "--------------- new_xmap map grid max " << grid_max.format() << std::endl;
   std::cout << "--------------- new_xmap map gr       " << gr.format() << std::endl;
   std::cout << "--------------- input xmap sampling   " << ip_xmap_grid_sampling.format() << std::endl;
   std::cout << "---------------   new xmap sampling   " << new_xmap_grid_sampling.format() << std::endl;


   // why am I using a new_xmap here? Just use xmap and then I don't have to copy it when it's done.
   //
   clipper::Xmap<float> new_xmap;
   clipper::Xmap<float>::Map_reference_index inx;
   new_xmap.init(clipper::Spacegroup::p1(), new_xmap_cell, new_xmap_grid_sampling);
   for (inx = new_xmap.first(); !inx.last(); inx.next()) new_xmap[inx]= 0;

   int u, v;
   clipper::Xmap_base::Map_reference_coord ix(new_xmap);
   for (u = gr.min().u(); u <= gr.max().u(); u++) {
      for (v = gr.min().v(); v <= gr.max().v(); v++) {
         for (ix.set_coord(clipper::Coord_grid(u,v,gr.min().w())); ix.coord().w() <= gr.max().w(); ix.next_w()) {

            // " " << ix.coord().coord_frac(new_xmap_grid_sampling).format()
            // << std::endl;

            // << " " << p.format() << std::endl;

             clipper::Coord_orth p = ix.coord().coord_frac(new_xmap_grid_sampling).coord_orth(new_xmap_cell);
            if (p.lengthsq() < radius*radius) {

               float dv = density_at_point(ip_xmap, p+centre);
               // std::cout << ix.coord().format() << " " << p.format() << std::endl;
               new_xmap[ix] = dv;
            }
         }
      }
   }
   xmap = new_xmap;
}

// transfer xmap (small) into xmap_p (big)
void
coot::util::map_fragment_info_t::unshift(clipper::Xmap<float> *xmap_p,
                                         const clipper::Coord_orth &centre_from_input) {

   // -----------------------------------------------------------------------------------------
   // -----------------------------------------------------------------------------------------
   //                                  currently broken
   // -----------------------------------------------------------------------------------------
   // -----------------------------------------------------------------------------------------
   // and difficult to fix.

   // do the opposite simple_origin_shift

   std::cout << "----------------------- unshift -------------- "<< std::endl;

   // first try, let's match how we transfer it in simple_origin_shift
   clipper::Xmap<float> &input_xmap(*xmap_p);

   clipper::Grid_map grid = make_grid_map(*xmap_p, centre_from_input);

   clipper::Xmap_base::Map_reference_coord ix(*xmap_p);
   for (int u = grid.min().u(); u < grid.max().u(); u++) {
      for (int v = grid.min().v(); v < grid.max().v(); v++) {
         for (int w = grid.min().w(); w < grid.max().w(); w++) {
            ix.set_coord(clipper::Coord_grid(u, v, w)); // don't copy this. using set_coord() is slow
            std::cout << " A " << ix.coord().format() << std::endl;
            clipper::Coord_grid cg = ix.coord() - offset;
            std::cout << " B " << cg.format() << std::endl;
            float f = xmap.get_data(cg);
            std::cout << " C " << f << std::endl;
            std::cout << " put from " << cg.format() << " to " << ix.coord().format() << " " << f << std::endl;
            input_xmap[ix] = f;
            if (false)
               std::cout << " put " << u << " " << v << " " << w << " "
                         << ix.coord().format() << " " << cg.format() << " " << f << std::endl;
            std::cout << "done " << u << " " << v << " " << w << std::endl;
         }
      }
   }
}


bool
coot::util::is_EM_map(const clipper::Xmap<float> &xmap) {

   // 20250524-PE this is not a good function - it doesn't check the labels.
   // The labels are only accessible when the map is read.
   // So asking is question of an xmap is not a good idea.
   // The status needs to be decided at map-read time and caried along
   // beside the map.

   bool is_em = false;
   if (xmap.spacegroup().num_symops() == 1) { // P1
      if (((xmap.cell().descr().alpha() - M_PI/2) <  0.0001) &&
          ((xmap.cell().descr().alpha() - M_PI/2) > -0.0001) &&
          ((xmap.cell().descr().beta()  - M_PI/2) > -0.0001) &&
          ((xmap.cell().descr().beta()  - M_PI/2) <  0.0001) &&
          ((xmap.cell().descr().gamma() - M_PI/2) > -0.0001) &&
          ((xmap.cell().descr().gamma() - M_PI/2) <  0.0001)) {
         is_em = true;
      }
   }
   return is_em;
}

#include "analysis/stats.hh"

// static
bool
coot::util::soi_variance::mri_var_pair_sorter(const std::pair<clipper::Xmap_base::Map_reference_index, float> &p1,
                                              const std::pair<clipper::Xmap_base::Map_reference_index, float> &p2) {
   return (p1.second < p2.second);
}

void
coot::util::soi_variance::proc(float solvent_content_frac) {

   // do I want a coord_grid?
   typedef std::pair<clipper::Xmap_base::Map_reference_index, float>  dd;
   std::vector<dd> data(200000);
   clipper::Xmap<float> variance_xmap = make_variance_map();
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = variance_xmap.first(); !ix.last(); ix.next()) {
      float fv = variance_xmap[ix];
      std::pair<clipper::Xmap_base::Map_reference_index, float> p(ix, fv);
      data.push_back(p);
   }

   std::cout << "INFO:: sorting variances " << std::endl;
   // high variance has high rank index
   std::sort(data.begin(), data.end(), mri_var_pair_sorter);
   std::cout << "INFO:: done sorting " << std::endl;

   clipper::Xmap<unsigned int> variance_rank_xmap;
   variance_rank_xmap.init(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   for (std::size_t i=0; i<data.size(); i++) {
      variance_rank_xmap[data[i].first] = i;
   }
   std::cout << "INFO:: done variance map " << std::endl;

   if (false) { // debug
      for (ix = variance_rank_xmap.first(); !ix.last(); ix.next()) {
         unsigned int iv = variance_rank_xmap[ix];
         clipper::Coord_grid cg = ix.coord();
         std::cout << "   " << cg.format() << " " << iv << std::endl;
      }
      // output rank map -> high variance over the gaps between atoms in the protein region
      clipper::CCP4MAPfile mapout;
      mapout.open_write("var-rank.map");
      mapout.export_xmap(variance_rank_xmap);
      mapout.close_write();

      // data for table
      std::ofstream f("var-ranks.tab");
      for (std::size_t i=0; i<data.size(); i++)
         f << "   " << i << " " << data[i].second << "\n";
      f.close();
   }

   std::size_t size = data.size();
   float size_f(size);
   const clipper::Xmap<float> &pt = protein_treatment_map();
   const clipper::Xmap<float> &st = solvent_treatment_map();

   if (solvent_content_frac > 0.75) solvent_content_frac = 0.75;
   if (solvent_content_frac < 0.25) solvent_content_frac = 0.25;

   // How about a cubic spline using point (solvent_content_frac, 0.5)

   clipper::Xmap<float> soi_xmap; // treatments applied
   // The top 25% is definitely protein, the bottom 25% is definitely solvent
   //
   soi_xmap.init(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      const unsigned int &idx = variance_rank_xmap[ix];
      float variance_rank_frac = static_cast<float>(idx)/size_f; // raw
      float fr = variance_rank_frac; // adjust this according to solvent content

      if (variance_rank_frac < 0.25) {
         fr = 0.0;
      } else {
         if (variance_rank_frac > 0.75) {
            fr = 1.0;
         } else {
            if (variance_rank_frac < solvent_content_frac) {
               float x_prime = fr - 0.25;
               float m = 0.5/(solvent_content_frac - 0.25);
               fr = m * x_prime;
            } else {
               float x_prime = fr - solvent_content_frac;
               float m = 2.0/(3.0 - 4.0 * solvent_content_frac);
               fr = m * x_prime + 0.5;
            }
         }
      }
      float other_frac = 1.0 - fr;
      soi_xmap[ix] = variance_rank_frac * pt[ix] + other_frac * st[ix];
   }

   if (true) { // debug
      clipper::CCP4MAPfile mapout;
      mapout.open_write("soi.map");
      mapout.export_xmap(soi_xmap);
      mapout.close_write();
   }
}

clipper::Xmap<float>
coot::util::soi_variance::make_variance_map() const {

   std::cout << "INFO:: making variance map" << std::endl;

   clipper::Xmap<float> var_map = xmap;

   const clipper::Cell &cell = xmap.cell();
   const clipper::Grid_sampling &gs = xmap.grid_sampling();

   // I don't much like this method, the grid spacing may be too course
   // to get a good sampling of points that are 2.42 A away

   float soi_dist  = 2.42;
   float soi_delta = 0.014; // for high-res test map this is a good starting value.
                            // We can be more clever and set this initial value
                            // dependent on grid_sampling/cell.
   std::vector<clipper::Coord_grid> soi_gps;

   bool n_in_sphere_is_ok = false; // to get started
   //
   while (! n_in_sphere_is_ok) {
      double ddmin = (soi_dist - soi_delta) * (soi_dist - soi_delta);
      double ddmax = (soi_dist + soi_delta) * (soi_dist + soi_delta);

      float radius = 2.5; // a bit more than 2.42
      clipper::Grid_range gr0(cell, gs, radius);
      for(clipper::Coord_grid cg=gr0.min(); cg!=gr0.max(); cg=cg.next(gr0)) {
         // std::cout << "   test " << cg.format() << std::endl;
         clipper::Coord_frac cf = cg.coord_frac(gs);
         clipper::Coord_orth co = cf.coord_orth(cell);
         double ll = co.lengthsq();
         if ((ll < ddmax) && (ll > ddmin)) {
            soi_gps.push_back(cg);
         }
      }

      if (soi_gps.size() < 50) {
         if (soi_delta < 0.6) {  // we will never match an echidnahedron
            soi_delta *= 1.8;
         } else {
            n_in_sphere_is_ok = true; // hacky escape for hideous low res maps
         }
      } else {
         n_in_sphere_is_ok = true;
      }
   }

   std::cout << "INFO:: Found " << soi_gps.size() << " SOI grid points " << std::endl;

   if (false) { // debuging SOI
      std::ofstream f("soi.tab");
      for (std::size_t i=0; i<soi_gps.size(); i++)
         f << i << " " << soi_gps[i].format() << "\n";
      f.close();
   }

   if (true) {
      auto tp_1 = std::chrono::high_resolution_clock::now();

      clipper::Xmap_base::Map_reference_index ix;
      unsigned int n_grids = 0;
      // is there a better way?
      for (ix = xmap.first(); !ix.last(); ix.next())
         ++n_grids;

      unsigned int n_threads = 4;
      std::vector<std::vector<clipper::Xmap_base::Map_reference_index> > grid_indices(n_threads);
      int r_reserve_size = std::lround(static_cast<float>(n_grids)/static_cast<float>(n_threads)) + 2;
      for (std::size_t n=0; n<n_threads; n++)
         grid_indices[n].reserve(r_reserve_size);
      unsigned int i_thread = 0; // insert to vector for this thread
      for (ix=xmap.first(); !ix.last(); ix.next()) {
         grid_indices[i_thread].push_back(ix);
         ++i_thread;
         if (i_thread==n_threads) i_thread=0;
      }

      auto tp_2 = std::chrono::high_resolution_clock::now();
      std::vector<std::thread> threads;
      for (std::size_t n=0; n<n_threads; n++) {
         threads.push_back(std::thread(soi_variance::apply_variance_values,
                                       std::ref(var_map),
                                       std::cref(xmap),
                                       std::cref(soi_gps),
                                       std::cref(grid_indices[n])));
      }

      auto tp_3 = std::chrono::high_resolution_clock::now();

      for (std::size_t ithread=0; ithread<threads.size(); ithread++)
         threads[ithread].join();

      auto tp_4 = std::chrono::high_resolution_clock::now();
      auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
      auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
      auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
      std::cout << "Timings:: make_variance_map(): grid spliting: grid-points  " << std::setw(4) << d21 << " milliseconds" << std::endl;
      std::cout << "Timings:: make_variance_map(): grid spliting: make-threads " << std::setw(4) << d32 << " milliseconds" << std::endl;
      std::cout << "Timings:: make_variance_map(): grid spliting: wait threads " << std::setw(4) << d43 << " milliseconds" << std::endl;
   }

   return var_map;
}

// static
void
coot::util::soi_variance::apply_variance_values(clipper::Xmap<float> &variance_map, // modify this
                                                const clipper::Xmap<float> &xmap,
                                                const std::vector<clipper::Coord_grid> &soi_gps,
                                                const std::vector<clipper::Xmap_base::Map_reference_index> &grid_indices) {

   for (std::size_t i=0; i<grid_indices.size(); i++) {
      const clipper::Xmap_base::Map_reference_index &ix = grid_indices[i];
      clipper::Coord_grid cg = ix.coord();
      std::vector<double> data(soi_gps.size());
      for (std::size_t j=0; j<soi_gps.size(); j++) {
         clipper::Coord_grid cg_soi_gp = cg + soi_gps[j];
         float fv = xmap.get_data(cg_soi_gp);
         data[j] = fv;
      }
      stats::single s(data); // Note: header only (we can't link libcoot-analysis from here)
      variance_map[ix] = s.variance();
   }
}


clipper::Xmap<float>
coot::util::soi_variance::solvent_treatment_map() const {

   clipper::Xmap<float> treated = xmap;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = treated.first(); !ix.last(); ix.next()) {
      treated[ix] = -treated[ix];
   }
   return treated;
}

clipper::Xmap<float>
coot::util::soi_variance::protein_treatment_map() const {

   clipper::Xmap<float> treated = xmap;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = treated.first(); !ix.last(); ix.next()) {
      float fv = treated[ix];
      if (fv < 0) {
         treated[ix] = 0.0;
      }
   }
   return treated;
}


#include "utils/dodec.hh"

coot::util::map_molecule_centre_info_t
coot::util::map_molecule_centre(const clipper::Xmap<float> &xmap) {

   // how about also returning a useful contour level?
   map_molecule_centre_info_t mmci;
   std::vector<clipper::Coord_orth> centres;

   std::pair<bool,clipper::Coord_orth> p(false, clipper::Coord_orth(0,0,0));
   clipper::Cell cell = xmap.cell();

   clipper::Coord_orth cell_centre(cell.descr().a() * 0.5,
                                   cell.descr().b() * 0.5,
                                   cell.descr().c() * 0.5);
   clipper::Coord_orth half_z(0, 0, 0.5 * cell.descr().c());
   centres.push_back(clipper::Coord_orth(0,0,0));
   centres.push_back(cell_centre);
   centres.push_back(half_z);

   clipper::Coord_orth best_centre = cell_centre;
   map_molecule_centre_info_t best_mmci;
   float best_score = 0.0;

   for(unsigned int i=0; i<centres.size(); i++) {
      map_molecule_centre_info_t r = map_molecule_recentre_from_position(xmap, centres[i]);
      if (r.success) {
         if (false)
            std::cout << "Starting at " << centres[i].format() << " gives score "
                      << r.success << " " << r.sum_of_densities << " "
                      << r.updated_centre.format() << std::endl;
         if (r.sum_of_densities > best_score) {
            best_mmci = r;
            best_score = r.sum_of_densities;
            best_centre = r.updated_centre;
         }
      } else {
         // std::cout << "Starting at " << centres[i].format() << " fail" << std::endl;
      }
   }

   if (best_score > 0.0) {
      mmci = best_mmci;
      std::cout << "INFO:: suggested centre " << best_centre.format() << std::endl;
      std::cout << "INFO:: suggested contour level " << best_mmci.suggested_contour_level
                << std::endl;
   }

   float sr = cell.descr().c() * 0.22;
   mmci.suggested_radius = sr;
   return mmci;
}

coot::util::map_molecule_centre_info_t
coot::util::map_molecule_recentre_from_position(const clipper::Xmap<float> &xmap,
                                                const clipper::Coord_orth &initial_centre) {

   map_molecule_centre_info_t mmci;
   clipper::Coord_orth current_centre = initial_centre;
   pentakis_dodec pdodec;
   double normalizer = 1.0/sqrt(3.0);
   std::vector<clipper::Coord_orth> directions;
   directions = pdodec.d.coords(); // 20 vertices
   for (unsigned int i=0; i<pdodec.pyrimid_vertices.size(); i++)
      directions.push_back(clipper::Coord_orth(pdodec.pyrimid_vertices[i] * normalizer));

   if (false) {
      for (unsigned int i=0; i<directions.size(); i++) {
         std::cout << i << " " << directions[i].format() << std::endl;
      }
   }

   unsigned int n_rounds = 20;
   for (unsigned int iround=0; iround<n_rounds; iround++) {
      clipper::Coord_orth sum(0,0,0);
      float sum_density = 0.0;
      double max_radius = 0.25 * xmap.cell().descr().a();
      unsigned int n_radii = 3;
      double radius_step = max_radius/static_cast<double>(n_radii);
      unsigned int point_count = 0;
      for (unsigned int ir=0; ir<n_radii; ir++) {
         double r = static_cast<double>(ir+1) * radius_step;
         for (unsigned int i=0; i<directions.size(); i++) {
            const clipper::Coord_orth &dir = directions[i];
            clipper::Coord_orth offset(dir * r);
            clipper::Coord_orth co(current_centre + offset);
            float d = density_at_point_by_nearest_grid(xmap, co);
            sum += clipper::Coord_orth(d * offset);
            if (d > 0.0) {
               point_count++;
               sum_density += d;
            }
         }
      }
      if (false)
         std::cout << "current_centre " << current_centre.format() << " " << sum.format()
                   << " " << sum_density << std::endl;
      if (point_count > 0) {
         float oopc = 1.0/static_cast<float>(point_count);
         float average_positive_density = sum_density * oopc;
         clipper::Coord_orth centroid(oopc * sum);
         current_centre += (0.1/average_positive_density) * centroid;
         mmci.updated_centre = current_centre;
         mmci.success = true;
         mmci.sum_of_densities = sum_density;
         mmci.suggested_contour_level = 2.0 * average_positive_density;
      } else {
         break;
      }
   }

   return mmci;
}

clipper::Xmap<float>
coot::util::power_scale(const clipper::Xmap<float> &xmap_1, const clipper::Xmap<float> &xmap_2) {

   float mg = coot::util::max_gridding(xmap_1); // A/grid
   clipper::Resolution reso(2.0 * mg); // Angstroms
   std::cout << "# making data info 1" << std::endl;
   clipper::HKL_info hkl_info_1(xmap_1.spacegroup(), xmap_1.cell(), reso, true);
   std::cout << "# making data info 2" << std::endl;
   clipper::HKL_info hkl_info_2(xmap_2.spacegroup(), xmap_2.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_1(hkl_info_1);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_2(hkl_info_2);
   std::cout << "# starting Fouriers" << std::endl;
   xmap_1.fft_to(fphis_1);
   std::cout << "# done map-1" << std::endl;
   xmap_2.fft_to(fphis_2);
   std::cout << "# done map-2" << std::endl;
   std::map<unsigned int, double> fref2_sums;
   std::map<unsigned int, double> fother2_sums;
   std::map<unsigned int, unsigned int> f2_counts;
   clipper::HKL_info::HKL_reference_index hri;
   for (hri = fphis_1.first(); !hri.last(); hri.next()) {
      float fr = fphis_1[hri].f();
      if (! clipper::Util::is_nan(fr)) {
         int h = hri.hkl().h();
         int k = hri.hkl().k();
         int l = hri.hkl().l();
         double rr = std::sqrt(h * h + k * k + l * l);
         unsigned int bin_index = static_cast<int>(rr);
         std::map<unsigned int, double>::const_iterator it = fref2_sums.find(bin_index);
         if (it == fref2_sums.end()) {
            fref2_sums[bin_index]   = fr * fr;
         } else {
            fref2_sums[bin_index]   += fr * fr;
         }
      }
   }
   for (hri = fphis_2.first(); !hri.last(); hri.next()) {
      float fo = fphis_2[hri].f();
      if (! clipper::Util::is_nan(fo)) {
         int h = hri.hkl().h();
         int k = hri.hkl().k();
         int l = hri.hkl().l();
         double rr = std::sqrt(h * h + k * k + l * l);
         unsigned int bin_index = static_cast<int>(rr);
         std::map<unsigned int, double>::const_iterator it = fother2_sums.find(bin_index);
         if (it == fother2_sums.end()) {
            fother2_sums[bin_index] = fo * fo;
            f2_counts[bin_index] = 1;
         } else {
            fother2_sums[bin_index] += fo * fo;
            f2_counts[bin_index] += 1;
         }
      }
   }
   for (hri = fphis_2.first(); !hri.last(); hri.next()) {
      float fo = fphis_2[hri].f();
      if (! clipper::Util::is_nan(fo)) {
         int h = hri.hkl().h();
         int k = hri.hkl().k();
         int l = hri.hkl().l();
         double rr = std::sqrt(h * h + k * k + l * l);
         unsigned int bin_index = static_cast<int>(rr);
         double sf = fref2_sums[bin_index] / fother2_sums[bin_index];
         fphis_2[hri].f() *= sf;
      }
   }

   std::map<unsigned int, double>::const_iterator it;
   for (it=fref2_sums.begin(); it!=fref2_sums.end(); ++it) {
      const auto &key(it->first);
      const auto &fref2_sum = it->second;
      const auto &fother = fother2_sums[key];
      float scale = fref2_sum / fother;
      std::cout << "compare-sums: " << key << "  " << fref2_sum << " " << fother << " scale " << scale << std::endl;
   }

   clipper::Xmap<float> xmap_scaled = xmap_2;
   xmap_scaled.fft_from(fphis_2);
   return xmap_scaled;
}



std::vector<std::pair<clipper::Resolution, double> >
coot::util::fsc(const clipper::Xmap<float> &xmap_1, const clipper::Xmap<float> &xmap_2) {

   std::vector<std::pair<clipper::Resolution, double> > v;
   std::cout << "# starting FSC" << std::endl;

   int n_bins = 100; //  pass this?
   float mg = coot::util::max_gridding(xmap_1); // A/grid
   clipper::Resolution reso(2.0 * mg); // Angstroms
   std::cout << "# making data info 1" << std::endl;
   clipper::HKL_info hkl_info_1(xmap_1.spacegroup(), xmap_1.cell(), reso, true);
   std::cout << "# making data info 2" << std::endl;
   clipper::HKL_info hkl_info_2(xmap_2.spacegroup(), xmap_2.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_1(hkl_info_1);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_2(hkl_info_2);
   std::cout << "# starting Fouriers" << std::endl;
   xmap_1.fft_to(fphis_1);
   std::cout << "# done map-1" << std::endl;
   xmap_2.fft_to(fphis_2);
   std::cout << "# done map-2" << std::endl;
   clipper::HKL_info::HKL_reference_index hri;
   std::vector<double> ff_sum(n_bins, 0.0);
   std::vector<double> f1_sqrd_sum(n_bins, 0.0);
   std::vector<double> f2_sqrd_sum(n_bins, 0.0);
   std::vector<unsigned int> counts(n_bins, 0);
   double max_reso = 0.0;
   for (hri = fphis_1.first(); !hri.last(); hri.next()) {
      float f = fphis_1[hri].f();
      if (! clipper::Util::is_nan(f)) {
         float irs = hri.invresolsq();
         float ir = sqrt(irs);
         if (ir > max_reso)
            max_reso = ir;
      }
   }
   for (hri = fphis_1.first(); !hri.last(); hri.next()) {
      float f_1 = 1000.0 * fphis_1[hri].f();
      try {
         if (! clipper::Util::is_nan(f_1)) {
            float f_2 = 1000.0 * fphis_2[hri].f();
            float irs = hri.invresolsq();
            float ir = sqrtf(irs);
            int bin_no = static_cast<int> (static_cast<double>(n_bins) * ir/max_reso);
            if (bin_no == n_bins) bin_no = n_bins - 1; // catch the reflection at the edge

            float A_f_1 = f_1 * cosf(fphis_1[hri].phi());
            float B_f_1 = f_1 * sinf(fphis_1[hri].phi());
            float A_f_2 = f_2 * cosf(fphis_2[hri].phi());
            float B_f_2 = f_2 * sinf(fphis_2[hri].phi());

            double prod = A_f_1 * A_f_2 + B_f_1 * B_f_2;
            ff_sum[bin_no] += prod;
            if (false)
               std::cout << hri.hkl().format() << " "
                         << f_1 << " " << fphis_1[hri].phi() << " "
                         << f_2 << " " << fphis_2[hri].phi() << " prod " << prod << std::endl;
            f1_sqrd_sum[bin_no] += f_1 * f_1;
            f2_sqrd_sum[bin_no] += f_2 * f_2;
            counts[bin_no]++;
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << rte.what() << std::endl;
      }
   }

   for(int i=0; i<n_bins; i++) {
      double fsc = ff_sum[i]/sqrt(f1_sqrd_sum[i] * f2_sqrd_sum[i]);
      double ir = (static_cast<double>(i) + 0.5) * max_reso/static_cast<double>(n_bins);
      std::cout << i << " " << ir << " "
                << counts[i] << " " << ff_sum[i]/static_cast<double>(counts[i]) << " "
                << f1_sqrd_sum[i]/static_cast<double>(counts[i]) << " "
                << f2_sqrd_sum[i]/static_cast<double>(counts[i]) << "    "
                << fsc << std::endl;
      std::pair<clipper::Resolution, double> p(clipper::Resolution(1.0/ir), fsc);
      v.push_back(p);
   }

   return v;
}

void
coot::util::compare_structure_factors(const clipper::Xmap<float> &xmap_1, const clipper::Xmap<float> &xmap_2) {

   float mg = coot::util::max_gridding(xmap_1); // A/grid
   clipper::Resolution reso(2.0 * mg); // Angstroms
   std::cout << "# making data info 1" << std::endl;
   clipper::HKL_info hkl_info_1(xmap_1.spacegroup(), xmap_1.cell(), reso, true);
   std::cout << "# making data info 2" << std::endl;
   clipper::HKL_info hkl_info_2(xmap_2.spacegroup(), xmap_2.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_1(hkl_info_1);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_2(hkl_info_2);

   auto fft_to_fphi = [] (const clipper::Xmap<float> &xmap,
                          clipper::HKL_data< clipper::datatypes::F_phi<float> > &fphis) { // reference
                         xmap.fft_to(fphis);
                      };

   std::cout << "# starting FFTs" << std::endl;
   std::thread thread_1(fft_to_fphi, std::cref(xmap_1), std::ref(fphis_1));
   std::thread thread_2(fft_to_fphi, std::cref(xmap_2), std::ref(fphis_2));
   thread_1.join();
   thread_2.join();

   std::cout << "# FFTs done" << std::endl;

   std::string file_name="compare-sfs.table";
   std::ofstream f(file_name);
   unsigned int count = 0;
   clipper::HKL_info::HKL_reference_index hri;
   for (hri = fphis_1.first(); !hri.last(); hri.next()) {
      if (count%1000==0) {
         float f_1   = fphis_1[hri].f();
         float phi_1 = fphis_1[hri].phi();
         try {
            if (! clipper::Util::is_nan(f_1)) {
               float f_2   = fphis_2[hri].f();
               if (! clipper::Util::is_nan(f_2)) {
                  float phi_2 = fphis_2[hri].phi();
                  float irs = hri.invresolsq();
                  float ir = sqrtf(irs);
                  clipper::HKL hkl = hri.hkl();
                  f << std::setw(9) << irs << " " << std::setw(9) << 1.0/ir << "   "
                    << std::setw(4) << hkl.h() << " " << std::setw(4) << hkl.k() << " " << std::setw(4) << hkl.l() << " "
                    << "Fs: " << f_1 << " " << f_2
                    << " phases: " << std::setw(9) << phi_1 << " " << std::setw(9) << phi_2 << "\n";
               }
            }
         }
         catch (const std::runtime_error &rte) {
            std::cout << rte.what() << std::endl;
         }
      }
      count++;
   }
   f.close();
}



void
coot::util::flip_hand(clipper::Xmap<float> *xmap_p) {

   std::vector<std::pair<clipper::Resolution, double> > v;
   clipper::Xmap<float> &xmap(*xmap_p);

   float mg = coot::util::max_gridding(xmap); // A/grid
   clipper::Resolution reso(2.0 * mg); // Angstroms
   clipper::HKL_info hkl_info(xmap.spacegroup(), xmap.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(hkl_info);
   xmap.fft_to(fphis);
   clipper::HKL_info::HKL_reference_index hri;
   for (hri = fphis.first(); !hri.last(); hri.next()) {
      fphis[hri].phi() = -fphis[hri].phi();
   }
  xmap.fft_from(fphis);

}



#include "coot-least-squares.hh"
#include <unordered_map>

// return the map of the density gradients
//
clipper::Xmap<float>
coot::util::zero_dose_extrapolation(const std::vector<std::pair<clipper::Xmap<float> *, float> > &xmaps,
                                    const clipper::Xmap<float> &xmap_mask) {

   auto all_maps_have_the_same_grid = [xmaps] () {
                                         return true; // teehee
                                      };

   auto debug_the_mask_map = [] (const clipper::Xmap<float> &xmap_mask) {
                                     clipper::Xmap_base::Map_reference_index ix;
                                     for (ix = xmap_mask.first(); !ix.last(); ix.next()) {
                                        clipper::Coord_grid cg = ix.coord();
                                        float dv = xmap_mask[ix];
                                        std::cout << "   " << cg.format() << " " << dv << std::endl;
                                     }
                                  };

   // returns a modified version of the input xmap
   auto apply_mask_to_map = [&xmap_mask] (const clipper::Xmap<float> &xmap) {
      clipper::Xmap<float> r = xmap;
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = xmap.first(); !ix.last(); ix.next() )  { // iterator index.
          r[ix] *= xmap_mask[ix];
      }
      return r;
   };

   class exponential_fit {
      double linear_m;
      double linear_c;
   public:
      explicit exponential_fit(const std::vector<std::pair<double, double> > &A_data) {
         std::vector<std::pair<double, double> > data(A_data.size());
         for (unsigned int i=0; i<A_data.size(); i++)
            if (A_data[i].second > 0.0) {
               double y = std::log(A_data[i].second);
               // std::cout << "adding data for linear fit " << A_data[i].first << " " << y << std::endl;
               data.push_back(std::make_pair(A_data[i].first, y));
            }

         least_squares_fit lsqf(data);
         linear_m = lsqf.m();
         linear_c = lsqf.c();
      }
      float at(const double &x) const {
         double y_hat = linear_m * x + linear_c;
         return exp(y_hat);
      }
   };

   class coot_hkl {
   public:
      int h,k,l;
      float invresolsq;
      coot_hkl(const clipper::HKL &hkl, const clipper::Cell &cell) {
         h = hkl.h();
         k = hkl.k();
         l = hkl.l();
         invresolsq = hkl.invresolsq(cell);
      }
      bool operator<(const coot_hkl &chkl) const { // for a map to work we need this function
         if (chkl.invresolsq < invresolsq) return true;  // check the sign
         if (chkl.invresolsq > invresolsq) return false;
         if (chkl.h > h) return true;
         if (chkl.h < h) return false;
         if (chkl.k > k) return true;
         if (chkl.k < k) return false;
         if (chkl.l > l) return true;
         if (chkl.l < l) return false;
         return false;
      }
      bool operator==(const coot_hkl &chkl) const {
         if (chkl.h == h)
            if (chkl.k == k)
               if (chkl.l == l)
                  return true;
         return false;
      }
   };

   std::cout << "DEBUG:: starting zero dose extrapolation" << std::endl;
   clipper::Xmap<float> xmap_result(*xmaps[0].first);

   const clipper::Xmap<float> &xmap_0 = *xmaps[0].first;
   clipper::Cell cell = xmap_0.cell();
   float mg = coot::util::max_gridding(xmap_0);
   clipper::Resolution reso(2.0 * mg);
   clipper::HKL_info myhkl_0(xmap_0.spacegroup(), xmap_0.cell(), reso, true);
   clipper::HKL_info myhkl_result_map(xmap_0.spacegroup(), xmap_0.cell(), reso, true);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis_for_result_map(myhkl_result_map);

   std::map<std::string, std::vector<double> > f_data_accum;
   std::map<coot_hkl, std::vector<std::pair<float, float> > > f_data_hkl_accum;
   std::map<std::string, float> f_resolution;

   // resize f_data_accum vectors so that I can fft and stuff the data into them in threads
   unsigned int n_maps = xmaps.size();
   clipper::HKL_info::HKL_reference_index hri;
   unsigned int count = 0;
   for (hri=fphis_for_result_map.first(); !hri.last(); hri.next()) {
      coot_hkl chkl(hri.hkl(), cell);
      f_data_hkl_accum[chkl].resize(n_maps);
      // std::cout << "resize refl " << chkl.h << " " << chkl.k << " " << chkl.l << " vec size "
      // << f_data_hkl_accum[chkl].size() << std::endl;
      count += 1;
   }
   std::cout << "debug:: resized " << count << " accumulation-reflections" << std::endl;

   auto accum_f_phi_data = [xmaps, &xmap_0, cell, reso, apply_mask_to_map] (const std::pair<unsigned int, unsigned int> main_index_range_pair,
                                                std::map<coot_hkl, std::vector<std::pair<float, float> > > *f_data_hkl_accum_p // results go here
                                                ) {

                              for (unsigned int imap=main_index_range_pair.first; imap<main_index_range_pair.second; imap++) {
                                 const clipper::Xmap<float> &xmap_in = *xmaps[imap].first;
                                 clipper::Xmap<float> xmap = apply_mask_to_map(xmap_in);
                                 clipper::HKL_info myhkl(xmap_0.spacegroup(), xmap_0.cell(), reso, true);
                                 clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
                                 xmap.fft_to(fphis);
                                 std::cout << "Accumulating f-data into coot-hkl-based map from map " << imap << std::endl;
                                 clipper::HKL_info::HKL_reference_index hri;
                                 for (hri = fphis.first(); !hri.last(); hri.next()) {
                                    coot_hkl chkl(hri.hkl(), cell);
                                    (*f_data_hkl_accum_p)[chkl][imap] = std::make_pair(fphis[hri].f(), fphis[hri].phi());
                                 }
                              }
                           };

#if 0 // single threaded
   // now we can do this block in parallel - lots of cache-sharing though :-/
   for (unsigned int imap=0; imap<xmaps.size(); imap++) {
      std::cout << "Try 2: FFT for map " << imap << std::endl;
      const clipper::Xmap<float> &xmap = *xmaps[imap].first;
      clipper::HKL_info myhkl(xmap_0.spacegroup(), xmap_0.cell(), reso, true);
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
      xmap.fft_to(fphis);
      std::cout << "Accumulating f-data into coot-hkl-based map from map " << imap << std::endl;
      for (hri = fphis.first(); !hri.last(); hri.next()) {
         coot_hkl chkl(hri.hkl(), cell);
         f_data_hkl_accum[chkl][imap] = std::make_pair(fphis[hri].f(), fphis[hri].phi());
         // std::cout << "refl " << chkl.h << " " << chkl.k << " " << chkl.l << " vec size "
         // << f_data_hkl_accum[chkl].size() << std::endl;
      }
   }
#endif

#if 1 // multi-threaded
   unsigned int n_threads = 8;
   std::vector<std::thread> threads;
   std::vector<std::pair<unsigned int, unsigned int> > ranges = atom_index_ranges(n_maps, n_threads);
   for (unsigned int irange=0; irange<ranges.size(); irange++)
      std::cout << "debug range " << irange << " :  " << ranges[irange].first << " " << ranges[irange].second << std::endl;

   for (unsigned int irange=0; irange<ranges.size(); irange++)
      threads.push_back(std::thread(accum_f_phi_data, ranges[irange], &f_data_hkl_accum));
   for (unsigned int ithread=0; ithread<threads.size(); ithread++)
      threads[ithread].join();
#endif


   std::cout << "Generating extrapolated structure factors" << std::endl;
   float f_max = 0.0;
   float av_devi_min = 999999999999.9;
   count = 0; // reset
   for (hri=fphis_for_result_map.first(); !hri.last(); hri.next()) {
      coot_hkl refl(hri.hkl(), cell);
      std::map<coot_hkl, std::vector<std::pair<float, float> > >::const_iterator it_hkl = f_data_hkl_accum.find(refl); // slow?
      if (it_hkl == f_data_hkl_accum.end()) {
         std::cout << "ERROR:: zde: this should not happen" << std::endl;
      } else {
         const std::vector<std::pair<float, float> > &f_phi_data = it_hkl->second;
         std::vector<std::pair<double, double> >   f_data(n_maps);
         std::vector<std::pair<double, double> > phi_data(n_maps);
         for (unsigned int ii=0; ii<n_maps; ii++) {
            f_data[ii]   = std::make_pair(ii+1, f_phi_data[ii].first);
            phi_data[ii] = std::make_pair(ii+1, f_phi_data[ii].second);
         }
         // exponential_fit exp_fit(f_data);
         coot::exponential_fit_with_offset exp_fit(f_data);
         float f = exp_fit.at(0.0);
         if (f > f_max) {
            coot_hkl chkl(hri.hkl(), cell);
            float f_frame_1 = f_data_hkl_accum[chkl][0].first;
            std::cout << "f_max updated " << f << " f-Frame-1 " << f_frame_1 << " for " << hri.hkl().format() << " reso " << refl.invresolsq << " a b c "
                      << exp_fit.a << " " <<  exp_fit.b << " " << exp_fit.c << " f_data_hkl_accum for-hkl size " << f_data_hkl_accum[chkl].size() << "\n";
            f_max = f;
            for (unsigned int ii=0; ii<f_data_hkl_accum[chkl].size(); ii++) {
               std::cout << "   debug:: " << ii << " " << f_data_hkl_accum[chkl][ii].first << std::endl;
            }
         }

         double av_devi = exp_fit.average_deviation(f_data);
         // std::cout << "av_devi " << av_devi << std::endl;
         if (av_devi < av_devi_min) {
            coot_hkl chkl(hri.hkl(), cell);
            float f_frame_1 = f_data_hkl_accum[chkl][0].first;
            std::cout << "av_devi_min updated " << av_devi << " f-Frame-1 " << f_frame_1 << " for " << hri.hkl().format() << " reso " << refl.invresolsq << " a b c "
                      << exp_fit.a << " " <<  exp_fit.b << " " << exp_fit.c << " f_data_hkl_accum for-hkl size " << f_data_hkl_accum[chkl].size() << "\n";
            f_max = f;
            for (unsigned int ii=0; ii<f_data.size(); ii++) {
               double x = f_data[ii].first;
               double y = f_data[ii].second;
               std::cout << "   debug:: " << ii << " " << x << " " << y << " " << exp_fit.at(x) << std::endl;
            }
            av_devi_min = av_devi;
         }
         // least_squares_fit lsq_f_fit(f_data);
         // f = lsq_f_fit.at(0);
         if (refl.invresolsq > 0.96) f = 0.0;
         // f = f_data[0].second;
         least_squares_fit lsq_phi_fit(phi_data);
         float phi_base = lsq_phi_fit.at(0.0);
         std::vector<std::pair<double, double> > zoned_data = phi_data;
         for (unsigned int ii=0; ii<n_maps; ii++) {
            if ((phi_data[ii].second - phi_base) >  M_PI) zoned_data[ii].second -= 2.0 * M_PI;
            if ((phi_data[ii].second - phi_base) < -M_PI) zoned_data[ii].second += 2.0 * M_PI;
         }
         least_squares_fit lsq_phi_zoned_fit(zoned_data);
         float phi = lsq_phi_zoned_fit.at(0.0);
         // phi = phi_data[0].second;
         fphis_for_result_map[hri].f()   = f;
         fphis_for_result_map[hri].phi() = phi;
         if (count%10000==0) {
            std::cout << "phase-compare reso " << refl.invresolsq << " first-phi " << phi_data[0].second << " extrap: " << phi << "\n";
         }
      }
      count++;
      if (count%100000==0)
         std::cout << "Extapolated " << count << " sfs" << std::endl;
   }

   std::cout << "INFO:: Generating final map..." << std::endl;
   xmap_result.fft_from(fphis_for_result_map);
   std::cout << "INFO:: ZDE done." << std::endl;
   return xmap_result;

}


// return the map of the density gradients
//
clipper::Xmap<float>
coot::util::analyse_map_point_density_change(const std::vector<std::pair<clipper::Xmap<float> *, float> > &xmaps,
                                             const clipper::Xmap<float> &xmap_for_mask) {

   // the caller guarantees that the map files are in the correct order in the vector

   auto all_maps_have_the_same_grid = [xmaps] () {
                                         return true;
                                      };

   auto make_linear_fit_map = [&xmaps] () {

                                 clipper::Xmap<float> linear_fit_map; // return this
                                 clipper::Xmap_base::Map_reference_index ix;
                                 clipper::Xmap<float> &first_map(*xmaps[0].first);

                                 clipper::Xmap<std::vector<float> > store;

                                 store.init(first_map.spacegroup(), first_map.cell(), first_map.grid_sampling());
                                 linear_fit_map.init(first_map.spacegroup(), first_map.cell(), first_map.grid_sampling());
                                 unsigned int n_maps = xmaps.size();
                                 std::cout << "INFO:: resizing store..." << std::endl;
                                 for (ix=first_map.first(); !ix.last(); ix.next())
                                    store[ix].resize(n_maps);
                                 std::cout << "INFO adding data to store..." << std::endl;
                                 for (unsigned int imap=0; imap<n_maps; imap++)
                                    for (ix=first_map.first(); !ix.last(); ix.next())
                                       store[ix][imap] = (*xmaps[imap].first)[ix];

                                 unsigned int count = 0; // for debugging
                                 for (ix=store.first(); !ix.last(); ix.next()) {
                                    const std::vector<float> &vec = store[ix];
                                    std::vector<std::pair<double, double> > data;
                                    for (unsigned int i=0; i<vec.size(); i++) {
                                       float f = store[ix][i];
                                       const float &rmsd = xmaps[i].second;
                                       float r = f/rmsd;
                                       data.push_back(std::make_pair(static_cast<double>(i), static_cast<double>(r)));
                                    }
                                    coot::least_squares_fit lsqf(data);
                                    linear_fit_map[ix] = lsqf.m();

                                    if (false) { // debugging - show me the data for some points within a radius
                                                 // of the centre of the map
                                       if (count%5000 == 0) {
                                          clipper::Cell cell = first_map.cell();
                                          float a = 0.5 * cell.a(); float b = 0.5 * cell.b(); float c = 0.5 * cell.c();
                                          clipper::Coord_orth cell_centre(a,b,c);
                                          clipper::Coord_grid cg = ix.coord();
                                          clipper::Coord_frac cf = cg.coord_frac(first_map.grid_sampling());
                                          clipper::Coord_orth co = cf.coord_orth(first_map.cell());
                                          const double dist_min_sqrd = 37 * 37;
                                          const double dist_max_sqrd = 90 * 90;
                                          double dd = (co-cell_centre).lengthsq();
                                          float f = first_map[ix];
                                          if (f > 0.2) {
                                             if (dd < dist_max_sqrd) {
                                                if (dd > dist_min_sqrd) {
                                                   std::string fn = "analyse-map-" + std::to_string(count) + ".table";
                                                   std::ofstream fout(fn.c_str());
                                                   for (unsigned int i=0; i<vec.size(); i++) {
                                                      fout << count << " " << i << "  " << vec[i] << "\n";
                                                   }
                                                   fout.close();
                                                }
                                             }
                                          }
                                       }
                                    }
                                    count++;
                                 }
                                 return linear_fit_map;
                              };


   clipper::Xmap<float> empty_map; // initially empty

   if (xmaps.size() < 2) return empty_map;
   if (! all_maps_have_the_same_grid()) return empty_map;

   clipper::Xmap<float> resulting_map = make_linear_fit_map();
   return resulting_map;

}

// user checks that the vector is not empty before calling this function.
clipper::Xmap<float>
coot::util::real_space_zero_dose_extrapolation(const std::vector<clipper::Xmap<float> *> &xmaps,
                                                const clipper::Xmap<float> &xmap_mask) {

   if (xmaps.empty())
      throw(std::runtime_error("empty xmaps"));
   clipper::Xmap<float> xmap_result(*xmaps[0]);
   clipper::Xmap_base::Map_reference_index ix;
   unsigned int n_masked = 0;
   unsigned int n_non_masked = 0;
   for (ix = xmap_result.first(); !ix.last(); ix.next() ) {
      if (xmap_mask[ix] > 0.2) {
         n_non_masked++;
         std::vector<std::pair<double, double> > f_values(xmaps.size());
         for (unsigned int imap=0; imap<xmaps.size(); imap++) {
            const auto &xmap = *xmaps[imap];
            float f = xmap[ix];
            f_values.push_back(std::make_pair(static_cast<double>(imap), static_cast<double>(f)));
         }
         least_squares_fit lsq(f_values);
         xmap_result[ix] = lsq.at(0);
      } else {
         n_masked++;
         xmap_result[ix] = 0.0f;
      }
   }
   std::cout << "masked counts " << n_masked << " " << n_non_masked << std::endl;
   unsigned int n_sum = n_masked + n_non_masked;
   float f1 = static_cast<float>(    n_masked)/static_cast<float>(n_sum);
   float f2 = static_cast<float>(n_non_masked)/static_cast<float>(n_sum);
   std::cout << "masked counts " << 100.0 * f1 << "% "  << 100.0 * f2 << "%" << std::endl;
   return xmap_result;
}



// negative becomes positive and positive becomes negative.
// Apply an offset so that most of the map is above zero.
//
void
coot::util::reverse_map(clipper::Xmap<float> *xmap_p) {

   // 20240413-PE do I really want to add a base?
   clipper::Xmap<float> &xmap(*xmap_p);
   std::pair<float, float> mv = mean_and_variance(xmap);
   float base = mv.first - 2.5f * mv.second;

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next() ) {
      float f = xmap[ix];
      xmap[ix] = -f - base;
   }
}


// this function should be passed a mask too.
// Otherwise we would spend a lot of time calculating
// distances of pixel in the map that noone cares about.
// This is currently prevented by testing for 0.0 in the input map
// but using a mask would be better.

std::vector<std::pair<std::string, clipper::Xmap<float> > >
coot::util::partition_map_by_chain(const clipper::Xmap<float> &xmap, mmdb::Manager *mol,
                                   std::string *state_string_p) {

   std::vector<std::pair<std::string, clipper::Xmap<float> > > v;

   // this is only called for amino acids
   auto make_side_chain_centre = [] (mmdb::Residue *residue_p) {
      bool status = false;
      clipper::Coord_orth centre;
      int count = 0;
      mmdb::Atom **residue_atoms = nullptr;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for(int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            if (atom_name == " N  " || atom_name == " H  " || atom_name == " C  " ||
               atom_name == " CA " || atom_name == " O  " || atom_name == " HA ")
               continue;
            centre += coot::co(at);
            count += 1;
         }
      }
      if (count > 0) {
         double sf = 1.0/static_cast<double>(count);
         centre = clipper::Coord_orth(centre.x() * sf, centre.y() * sf, centre.z() * sf);
         status = true;
      }
      return std::make_pair(status, centre);
   };

   auto make_reference_points_for_chains = [make_side_chain_centre] (mmdb::Manager *mol,
                                                const clipper::Cell &cell,
                                                const clipper::Grid_sampling &gs) {
      // chain-id and CA positions as grid points
      std::map<std::string, std::vector<clipper::Coord_grid> > coordinates_grid_points_map;
      for (int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               std::string chain_id(chain_p->GetChainID());
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at->isTer()) {
                           std::string name(at->GetAtomName());
                           if (name == " CA " || name == " P  " || name == " C4'" ||
                               name == " N9 " || name == " N1 " || name == " C4 ") {
                              clipper::Coord_orth co = coot::co(at);
                              clipper::Coord_frac cf = co.coord_frac(cell);
                              clipper::Coord_grid cg = cf.coord_grid(gs);
                              coordinates_grid_points_map[chain_id].push_back(cg);
                           }
                        }
                     }
                     if (residue_p->isAminoacid()) {
                        std::pair<bool, clipper::Coord_orth> side_chain_centre =
                           make_side_chain_centre(residue_p);
                        if (side_chain_centre.first) {
                           clipper::Coord_frac cf = side_chain_centre.second.coord_frac(cell);
                           clipper::Coord_grid cg = cf.coord_grid(gs);
                           coordinates_grid_points_map[chain_id].push_back(cg);
                        }
                     }
                  }
               }
            }
         }
      }
      return coordinates_grid_points_map;
   };

   auto split_biggest_chain = [] (std::map<std::string, std::vector<clipper::Coord_grid> > *rp_p) {
      std::map<std::string, std::vector<clipper::Coord_grid> >::const_iterator it;
      unsigned int n_biggest = 0;
      std::string biggest_chain;
      for (it=rp_p->begin(); it!=rp_p->end(); ++it) {
         unsigned int n = it->second.size();
         if (n > n_biggest) {
            const std::string chain_id = it->first;
            n_biggest = n;
            biggest_chain = chain_id;
         }
      }
      if (n_biggest > 0) {
         std::vector<clipper::Coord_grid> &v = (*rp_p)[biggest_chain];
         size_t half_size = n_biggest / 2;
         if (n_biggest % 2 != 0) half_size += 1;
         std::vector<clipper::Coord_grid> new_v;
         for (unsigned int i=half_size; i<v.size(); i++)
            new_v.push_back(v[i]);
         v.resize(half_size);
         std::string new_chain_id = biggest_chain + "+";
         (*rp_p)[new_chain_id] = new_v;
      }
   };

   std::cout << "Making reference points for chains" << std::endl;
   if (state_string_p) *state_string_p = "Making reference points for chains";
   clipper::Cell cell = xmap.cell();
   clipper::Spacegroup sg = xmap.spacegroup();
   clipper::Grid_sampling gs = xmap.grid_sampling();
   std::map<std::string, std::vector<clipper::Coord_grid> > rp =
      make_reference_points_for_chains(mol, cell, gs);

   split_biggest_chain(&rp);

   std::vector<std::string> chain_ids;
   for(const auto &item : rp)
      chain_ids.push_back(item.first);
   clipper::Xmap<std::map<std::string, int> > distance_map;
   distance_map.init(sg, cell, gs);

   std::cout << "INFO:: Filling distance map with initial values" << std::endl;
   if (state_string_p) *state_string_p = "Filling distance map with initial values";
   std::map<std::string, int> starting_distance_map;
   for (const auto &item : chain_ids)
      starting_distance_map[item] = 999999;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = distance_map.first(); !ix.last(); ix.next())  {
      clipper::Coord_grid cg = ix.coord();
      float f = xmap.get_data(cg);
      if (f != 0.0)
         distance_map[ix] = starting_distance_map;
   }

   auto manhattan_check = +[] (const std::string &chain_id,
                               const std::vector<clipper::Coord_grid> &reference_points,
                               clipper::Xmap<std::map<std::string, int> > *distance_map_p,
                               std::string *info_string_p) {

      if (info_string_p) *info_string_p = "Distance check for " + chain_id;
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = distance_map_p->first(); !ix.last(); ix.next()) {
         std::map<std::string, int>::const_iterator it;
         if ((*distance_map_p)[ix].empty()) continue;
         auto &current_dist = (*distance_map_p)[ix][chain_id];
         clipper::Coord_grid cg = ix.coord();
         for (unsigned int i=0; i<reference_points.size(); i++) {
            const auto &ref_grid_pt = reference_points[i];
            int d_u = ref_grid_pt.u() - cg.u();
            int d_v = ref_grid_pt.v() - cg.v();
            int d_w = ref_grid_pt.w() - cg.w();
            int manhattan_dist = abs(d_u) + abs(d_v) + abs(d_w);
            if (manhattan_dist < current_dist) {
               current_dist = manhattan_dist;
            }
         }
      }
      std::cout << "INFO:: manhattan map check done " << chain_id << std::endl;
   };

   std::vector<std::thread> threads;
   for (const auto &chain_id : chain_ids) {
      // this block can be multi-threaded, I think
      const std::vector<clipper::Coord_grid> &reference_points = rp[chain_id];
      // modify distance map
      // manhattan_check(chain_id, reference_points, &distance_map);
      threads.push_back(std::thread(manhattan_check, chain_id, reference_points, &distance_map,
                                    state_string_p));
   }

   if (state_string_p) *state_string_p = "Joining threads...";
   std::cout << "INFO:: joining threads" << std::endl;
   for (auto &thread : threads) thread.join();

   // now extract each of the maps for each chain
   std::cout << "INFO:: now constructing the map for each chain" << std::endl;
   if (state_string_p) *state_string_p = "Constructing the map for each chain";

   clipper::Xmap<std::string> chain_map;
   chain_map.init(sg, cell, gs);
   for (ix = distance_map.first(); !ix.last(); ix.next()) {
      const auto &dist_map = distance_map[ix];
      std::string best_chain_id;
      int dist_best = 999999;
      std::map<std::string, int>::const_iterator it;
      for (it=dist_map.begin(); it!=dist_map.end(); ++it) {
         if (it->second < dist_best) {
            dist_best = it->second;
            best_chain_id = it->first;
            // std::cout << " best-chain " << ix.coord().format() << " " << best_chain_id << std::endl;
         }
      }
      // chaid-ids ending in "+" really are part of another chain.
      if (best_chain_id.size() == 2)
         if (best_chain_id[1] == '+')
            best_chain_id = best_chain_id[0];
      chain_map[ix] = best_chain_id;
   }

   for(const auto &chain_id : chain_ids) {
      // chaid-ids ending in "+" really are part of another chain.
      if (chain_id.size() == 2)
         if (chain_id[1] == '+')
            continue;
      std::cout << "INFO:: constructing map for chain " << chain_id << std::endl;
      if (state_string_p) *state_string_p = "Constructing map for chain " + chain_id;
      clipper::Xmap<float> map_for_chain;
      map_for_chain.init(sg, cell, gs);
      for (ix = chain_map.first(); !ix.last(); ix.next()) {
         const auto &chain_for_this_grid_point = chain_map[ix];
         clipper::Coord_grid cg = ix.coord();
         if (chain_for_this_grid_point == chain_id) {
            float f = xmap.get_data(cg);
            map_for_chain.set_data(cg, f);
         } else {
            map_for_chain.set_data(cg, 0.0f);
         }
      }
      v.push_back(std::make_pair(chain_id, map_for_chain));
   }

   return v;
}


clipper::Xmap<float>
coot::util::make_map_mask(const clipper::Spacegroup &space_group,
                          const clipper::Cell &cell,
                          const clipper::Grid_sampling &grid_sampling,
                          mmdb::Manager *mol,
                          int atom_selection_handle,
                          float radius,
                          float smooth) {

   float outer_space_f = -1.1f;
   clipper::Xmap<float> xmap(space_group, cell, grid_sampling);
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next())
      xmap[ix] = outer_space_f;

   mmdb::Atom **selected_atoms = 0;
   int n_atoms = 0;
   mol->GetSelIndex(atom_selection_handle, selected_atoms, n_atoms);
   for (int iat=0; iat<n_atoms; iat++) {
      mmdb::Atom *at = selected_atoms[iat];
      if (! at->isTer()) {

         float atom_radius_sq = radius * radius;
         clipper::Coord_orth pt = co(at);
         clipper::Coord_frac cf = pt.coord_frac(cell);
         clipper::Coord_frac box0(
                                  cf.u() - radius/cell.descr().a(),
                                  cf.v() - radius/cell.descr().b(),
                                  cf.w() - radius/cell.descr().c());
         clipper::Coord_frac box1(
                                  cf.u() + radius/cell.descr().a(),
                                  cf.v() + radius/cell.descr().b(),
                                  cf.w() + radius/cell.descr().c());
         clipper::Grid_map grid(box0.coord_grid(grid_sampling),
                                box1.coord_grid(grid_sampling));

         if (smooth > 0.0) {
            clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min() ), iu, iv, iw;
            for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
               for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
                  for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                     clipper::Coord_orth gp_orth = iw.coord().coord_frac(grid_sampling).coord_orth(cell);
                     float dd = (gp_orth - pt).lengthsq();
                     const float &current_v = xmap[iw];
                     if (current_v < 1.0) {
                        if (dd < atom_radius_sq) {
                           float d = sqrtf(dd);
                           float v_delta = 1.0 - d/radius;
                           xmap[iw] += v_delta;
                           if (xmap[iw] > 1.0) xmap[iw] = 1.0;
                        }
                     }
                  }
               }
            }

         } else {
            // no smoothing
            clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min() ), iu, iv, iw;
            for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
               for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
                  for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                     float dd = (iw.coord().coord_frac(grid_sampling).coord_orth(cell) - pt).lengthsq();
                     if (dd < atom_radius_sq) {
                        xmap[iw] = 1.0;
                     }
                  }
               }
            }
         }
      }
   }

   // remove voids, the starting points of which are the masking values between 1.0 and 0.0

   bool continue_removing_voids = true;
   int count = 0;
   float cut_off = 1.0f;

   clipper::Xmap<short int> considered_map(space_group, cell, grid_sampling);
   clipper::Xmap<short int>  outer_surface(space_group, cell, grid_sampling);
   for (ix = considered_map.first(); !ix.last(); ix.next()) considered_map[ix] = 0;
   for (ix = outer_surface.first();  !ix.last(); ix.next())  outer_surface[ix] = 0;

   clipper::Skeleton_basic::Neighbours neighb(xmap);

   while (continue_removing_voids) {

      if (false)
         std::cout << "void removal round " << count << std::endl;
      continue_removing_voids = false;

      std::queue<clipper::Coord_grid> q;
      for (ix = xmap.first(); !ix.last(); ix.next()) {
         const float &xmix = xmap[ix];
         if (xmix < cut_off) {
            if (xmix > 0.0f) {
               q.push(ix.coord());
            }
         }

         if (! q.empty()) {
            if (considered_map.get_data(ix.coord()) == 1) {
            } else {
               considered_map[ix] = 1;
               std::vector<clipper::Coord_grid> cluster_grid_points;
               while (! q.empty()) {
                  clipper::Coord_grid c_g_start = q.front();
                  q.pop();
                  cluster_grid_points.push_back(c_g_start);
                  for (int i=0; i<neighb.size(); i++) {
                     clipper::Coord_grid c_g = c_g_start + neighb[i];
                     if (considered_map.get_data(c_g) == 0) {
                        float f = xmap.get_data(c_g);
                        if (f < cut_off) {
                           considered_map.set_data(c_g, 1);
                           cluster_grid_points.push_back(c_g);
                           if (f > 0.0f)
                              q.push(c_g);
                        }
                     }
                  }
               }
               // std::cout << "cluster size " << cluster_grid_points.size() << " for " << ix.coord().format() << std::endl;

               if (false) { // debug
                  if (cluster_grid_points.size() > 100) {
                     for (unsigned int ii=0; ii<cluster_grid_points.size(); ii++) {
                        const auto &gp = cluster_grid_points[ii];
                        std::cout << " cluster-grid-point " << gp.u() << " " << gp.v() << " " << gp.w() << " "
                                  << xmap.get_data(gp) << std::endl;
                     }
                  }
               }

               if (cluster_grid_points.size() < 100) { // need testing
                  for (unsigned int ii=0; ii<cluster_grid_points.size(); ii++) {
                     const auto &gp = cluster_grid_points[ii];
                     xmap.set_data(gp, 1.0f);
                  }
                  continue_removing_voids = true;
               } else {
                  for (unsigned int ii=0; ii<cluster_grid_points.size(); ii++) {
                     const auto &gp = cluster_grid_points[ii];
                     outer_surface.set_data(gp, 1);
                  }
               }
            }
         }
      }
      count += 1;
   }

   // what's left is voids with values of actually zero

   continue_removing_voids = true;
   while (continue_removing_voids) {

      continue_removing_voids = false;

      for (ix = xmap.first(); !ix.last(); ix.next()) {
         const float &xmix = xmap[ix];
         if (xmix == 0.0f) {
            clipper::Coord_grid c_g_start = ix.coord();
            if (considered_map.get_data(c_g_start) == 1) {
            } else {
               bool has_non_zero_neighb = false;
               for (int i=0; i<neighb.size(); i++) {
                  clipper::Coord_grid c_g = c_g_start + neighb[i];
                  if (xmap.get_data(c_g) > 0.0) {
                     if (outer_surface.get_data(c_g) != 1) {
                        has_non_zero_neighb = true;
                        break;
                     }
                     if (has_non_zero_neighb) {
                        xmap.set_data(c_g, 1.0);
                        continue_removing_voids = true;
                     }
                  }
               }
            }
         }
      }
   }

   return xmap;
}


#include <string.h>
#include "fib-sphere.hh"

int
coot::util::split_residue_using_map(mmdb::Residue *residue_p,
                                    mmdb::Manager *mol,
                                    const clipper::Xmap<float> &xmap) {

   auto direction_of_a_to_b_alt_conf = [] (mmdb::Residue *residue_p) {

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::map<std::string, std::pair<mmdb:: Atom *, mmdb:: Atom *> > atom_name_map;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            // let's handle the perverse case where B alt confs come before the A alt confs
            std::string atom_name(at->name);
            std::map<std::string, std::pair<mmdb:: Atom *, mmdb:: Atom *> >::iterator it;
            it = atom_name_map.find(atom_name);
            std::string alt_loc = at->altLoc;
            if (it == atom_name_map.end()) {
               if (alt_loc == "A") atom_name_map[atom_name] = std::make_pair(at, nullptr);
               if (alt_loc == "B") atom_name_map[atom_name] = std::make_pair(nullptr, at);
            } else {
               if (alt_loc == "A") it->second.first  = at; // shouldn't happen
               if (alt_loc == "B") it->second.second = at; // normal case
            }
         }
      }

      if (false)
         std::cout << ":::::::: debug:: split_residue_using_map() here with residue "
                  << coot::residue_spec_t(residue_p)
                  << " " << residue_p->GetResName() <<  " with atom_name_map size "
                   << atom_name_map.size() << std::endl;

      bool status = false;
      clipper::Coord_orth a_b_uv(0,0,0);
      clipper::Coord_orth sum_diff(0,0,0);
      if (! atom_name_map.empty()) {
         std::map<std::string, std::pair<mmdb:: Atom *, mmdb:: Atom *> >::iterator it;
         for (it=atom_name_map.begin(); it!=atom_name_map.end(); ++it) {
            const auto &atom_name(it->first);
            mmdb:: Atom *at_1 = it->second.first;
            mmdb:: Atom *at_2 = it->second.second;
            if (at_1 && at_2) {
               status = true;
               clipper::Coord_orth pos_1 = co(at_1);
               clipper::Coord_orth pos_2 = co(at_2);
               clipper::Coord_orth diff = pos_2 - pos_1;
               sum_diff += diff;
            }
         }
      }
      if (status) {
         a_b_uv = clipper::Coord_orth(sum_diff.unit());
      }
      return std::pair<bool, clipper::Coord_orth> (status, a_b_uv);
   };

   auto split_residue = [] (mmdb::Residue *residue_p, const clipper::Coord_orth &v) {

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      clipper::Coord_orth h = 0.5 * v;
      std::vector<mmdb::Atom *> atoms_to_be_added;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         mmdb::Atom *at_copy = new mmdb::Atom;
         at_copy->Copy(at);
         strncpy(at->altLoc,      "A", 2);
         strncpy(at_copy->altLoc, "B", 2);
         at->x      += h.x();      at->y += h.y();      at->z += h.z();
         at_copy->x -= h.x(); at_copy->y -= h.y(); at_copy->z -= h.z();
         atoms_to_be_added.push_back(at_copy);
      }
      for(mmdb::Atom *at_copy : atoms_to_be_added)
        residue_p->AddAtom(at_copy);
   };

   int status = 0;

   // 0.6, 0.75, 0.9
   std::vector<clipper::Coord_orth> fs = fibonacci_sphere(50);

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   std::vector<std::pair<clipper::Coord_orth, float> > atom_map_vector;
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         clipper::Coord_orth at_pos = co(at);
         float d_at = density_at_point(xmap, at_pos);
         clipper::Coord_orth vector_sum(0,0,0);
         float w_sum = 0.0;
         for (const auto &sf : {0.6, 0.75, 0.9}) {
            for (unsigned int i=0; i<fs.size(); i++) {
               if (fs[i].z() < 0.0) continue;
               clipper::Coord_orth pt_1 = at_pos - sf * fs[i];
               clipper::Coord_orth pt_2 = at_pos + sf * fs[i];
               float d_1 = density_at_point(xmap, pt_1);
               float d_2 = density_at_point(xmap, pt_2);
               float w = d_1 + d_2 - 2.0f * d_at;
               if (w > 0.0) {
                  clipper::Coord_orth delta = (pt_2 - pt_1);
                  vector_sum += w * delta;
                  w_sum += w;
               }
            }
         }
         atom_map_vector.push_back(std::make_pair(vector_sum, w_sum));
      }
   }

   clipper::Coord_orth vector_sum(0,0,0);
   float w_sum = 0.0;
   for (unsigned int i=0; i<atom_map_vector.size(); i++) {
      const float &w = atom_map_vector[i].second;
      vector_sum += w * atom_map_vector[i].first;
      w_sum += w;
   }

   double one_over_w_sum = 1.0 / w_sum;
   clipper::Coord_orth map_vector = one_over_w_sum * vector_sum;
   clipper::Coord_orth vec = clipper::Coord_orth(map_vector.unit());

   mmdb::Residue *residue_prev = previous_residue(residue_p);
   mmdb::Residue *residue_next = next_residue(residue_p);

   clipper::Coord_orth neighbour_vec(0,0,0);
   bool has_altconfed_neighbour = false;
   if (residue_prev) {
      auto d = direction_of_a_to_b_alt_conf(residue_prev);
      if (d.first) {
         has_altconfed_neighbour = true;
         neighbour_vec += d.second;
      }
   }
   if (residue_next) {
      auto d = direction_of_a_to_b_alt_conf(residue_next);
      if (d.first) {
         has_altconfed_neighbour = true;
         neighbour_vec += d.second;
      }
   }
   if (has_altconfed_neighbour) {
      double dp_1 = clipper::Coord_orth::dot(vec, neighbour_vec);
      double dp_2 = clipper::Coord_orth::dot(vec, -neighbour_vec);

      // std::cout << "dp_1 " << dp_1 << " dp_2 " << dp_2 << std::endl;
      if (dp_1 > 0.0)
         vec = - vec;
   }

   split_residue(residue_p, vec);

   return status;

}

// pt_ref must not be co-linear with pt_1-pt_2
//
//   pt_1 ------------- pt_2
//       /
//      /
//     /
//   pt_ref
//
std::vector<std::vector<float> >
coot::util::get_density_on_cylinder(const clipper::Coord_orth &pt_1, const clipper::Coord_orth &pt_2,
                                    const clipper::Coord_orth &pt_ref, const clipper::Xmap<float> &xmap,
                                    double radius, unsigned int n_length, unsigned int n_ring) {

   bool debug = false;
   std::vector<std::vector<float> > v;
   v.reserve(n_length);
   clipper::Coord_orth v1 = pt_2 - pt_1;
   clipper::Coord_orth v0 = pt_1 - pt_ref;

   clipper::Coord_orth v_up_uv(clipper::Coord_orth::cross(clipper::Coord_orth(v1.unit()),
                                                          clipper::Coord_orth(v0.unit())));
   clipper::Coord_orth v_tube_uv(v1.unit());
   clipper::Coord_orth v_perp_uv(clipper::Coord_orth::cross(v_tube_uv, v_up_uv));

   clipper::Coord_orth v_step = 1.0/static_cast<double>(n_length) * v1;

   for (unsigned int i=0; i<=n_length; i++) {
      clipper::Coord_orth tube_mid = pt_1 + static_cast<double>(i) * v_step;
      clipper::Coord_orth circle_start = tube_mid + radius * v_perp_uv;
      std::vector<float> d_ring(n_ring, 0.0f);
      for (unsigned int j=0; j<n_ring; j++) {
         // rotate circle_start around the vector v1:
         double angle = static_cast<double>(j)/static_cast<double>(n_ring) * 2.0 * M_PI;
         // args: direction, position, origin_shift, angle
         clipper::Coord_orth rotated_pt =
            coot::util::rotate_around_vector(v1, circle_start, pt_1, angle);
         float d = coot::util::density_at_point(xmap, rotated_pt);
         if (debug)
            std::cout << "rotated_pt: " << i << " " << j << " "
                      << rotated_pt.x() << " " << rotated_pt.y() << " " << rotated_pt.z()
                      << " density: " << d << std::endl;
         d_ring[j] = d;
      }
      v.push_back(d_ring);
   }
   return v;

}
