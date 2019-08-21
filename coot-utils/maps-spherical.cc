/* ligand/maps-spherical.cc
 * 
 * Copyright 2014 by The Medical Research Council
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */


#include <gsl/gsl_sf_bessel.h>

#include <clipper/mmdb/clipper_mmdb.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/clipper-contrib.h> // sfc

#include "utils/coot-utils.hh"
#include "coot-map-utils.hh"
#include "emma.hh"
#include "peak-search.hh"
#include "xmap-stats.hh" // needed?

void
coot::util::emma::sfs_from_boxed_molecule(mmdb::Manager *mol_orig, float border) {

   mmdb::Manager *mol = new mmdb::Manager;
   mol->Copy(mol_orig, mmdb::MMDBFCM_All);
   mmdb::PPAtom atom_selection = 0;
   int n_selected_atoms;
   // now do selection
   int SelHnd = mol->NewSelection(); // d
   mol->SelectAtoms(SelHnd, 1, "*",
		    mmdb::ANY_RES, "*",
		    mmdb::ANY_RES, "*",
		    "*", "*", 
		    "*", // elements
		    "*", // alt loc.
		    mmdb::SKEY_NEW);
   mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);

   std::pair<bool, clipper::Coord_orth> centre = centre_of_molecule(mol);
   if (centre.first) {

      // move the coordinates so that the middle of the molecule is at the origin.
      // shift(mol, -centre.second);
      std::pair<clipper::Coord_orth, clipper::Coord_orth> e = extents(mol, SelHnd);
      double x_range = e.second.x() - e.first.x();
      double y_range = e.second.y() - e.first.y();
      double z_range = e.second.z() - e.first.z();

      std::cout << "DEBUG:: molecule  centre after recentering: "
		<< centre.first << " " << centre.second.format() << std::endl;

      double nr = clipper::Util::d2rad(90);
      clipper::Cell_descr cell_descr(x_range + 2*border,
				     y_range + 2*border,
				     z_range + 2*border, nr, nr, nr);
      cell = clipper::Cell(cell_descr);
      spacegroup = clipper::Spacegroup::p1();
      reso = clipper::Resolution(3.0);

      // calculate structure factors
      hkl_info = clipper::HKL_info(spacegroup, cell, reso);
      hkl_info.generate_hkl_list();
      std::cout << "P1-sfs: num_reflections: " << hkl_info.num_reflections() << std::endl;

      // with bulking
      // clipper::SFcalc_obs_bulk<float> sfcb;
      // sfcb( fc, fo, atoms );
      // bulkfrc = sfcb.bulk_frac();
      // bulkscl = sfcb.bulk_scale();

      // debug
      //
      std::cout << "P1-sfs: cell " << cell.format() << std::endl;
      std::cout << "P1-sfs: resolution limit " << reso.limit() << std::endl;

      clipper::MMDBAtom_list atoms(atom_selection, n_selected_atoms);
      std::cout << "P1-sfs: n_selected_atoms: " << n_selected_atoms << std::endl;

      fc_from_model = clipper::HKL_data<clipper::data32::F_phi>(hkl_info, cell);

      clipper::SFcalc_aniso_fft<float> sfc;
      sfc(fc_from_model, atoms);
      std::cout << "P1-sfs: done sfs calculation " << std::endl;

      if (0) { // debug
	 clipper::HKL_info::HKL_reference_index hri;
	 for (hri = fc_from_model.first(); !hri.last(); hri.next()) {
	    std::cout << "fc " << hri.hkl().format() << " "
		      << fc_from_model[hri].f()   << " "
		      << fc_from_model[hri].phi() << "\n";
	 }
      }
   }


   // Write out the map from the calculated sfs:
   //
   {
      clipper::Grid_sampling gs(spacegroup, cell, reso);
      clipper::Xmap<float> fcalc_map(spacegroup, cell, gs);
      fcalc_map.fft_from(fc_from_model);
      clipper::CCP4MAPfile mapout;
      mapout.open_write("fc_from_model.map");
      mapout.export_xmap(fcalc_map);
      mapout.close_write();
      mol->WritePDBASCII("shifted_model.pdb");
   }
      
   
   mol->DeleteSelection(SelHnd);
   delete mol;
}


void
coot::util::emma::overlap(const clipper::Xmap<float> &xmap) const {

   // sfs from map
   clipper::Resolution reso(3);
   clipper::HKL_info hkl_info(xmap.spacegroup(), xmap.cell(), reso);
   hkl_info.generate_hkl_list();
   clipper::HKL_data< clipper::datatypes::F_phi<double> > map_fphidata(hkl_info);
   xmap.fft_to(map_fphidata, clipper::Xmap_base::Normal);

   clipper::Range<double>    fc_range = fc_from_model.invresolsq_range();
   clipper::Range<double> map_f_range = map_fphidata.invresolsq_range();
   std::cout << "fc from model resolution ranges " << fc_range.min() << " " << fc_range.max() << std::endl;
   std::cout << "fc from model resolution ranges " << 1/sqrt(fc_range.min()) << " " << 1/sqrt(fc_range.max())
	     << std::endl;
   std::cout << "SFs from map  resolution ranges " << map_f_range.min() << " " << map_f_range.max() << std::endl;
   std::cout << "SFs from map  resolution ranges " << 1/sqrt(map_f_range.min()) << " " << 1/sqrt(map_f_range.max())
	     << std::endl;


   // N is number of points for Gauss-Legendre Integration.
   // Run over model reflection list, to make T(k), k = 1,N
   // 
   int N = 16;
   float f = 1/float(N);
   gauss_legendre_t gl; // 1-indexed
   float radius_max = 15; // max integration radius
   float radius_min =  0; // min integration radius
   float radius_range = radius_max - radius_min;
   float radius_mid = (radius_max + radius_min)*0.5;
   // store T(k) 
   std::complex<double> zero_c(0.0, 0.0);
   std::vector<std::complex<double> > T(N+1.0, zero_c);
   for (float r_k_i=1; r_k_i<=N; r_k_i++) {
      float w_k = gl.weight  (r_k_i);
      float r_k = gl.abscissa(r_k_i);
      float x_gl = radius_range*gl.abscissa(r_k_i)*0.5 + radius_mid;
      
      clipper::HKL_info::HKL_reference_index hri;
      for (hri = fc_from_model.first(); !hri.last(); hri.next()) {
	 
	 // r and x_gl are in inversely-related spaces: when
	 // fc_range.max() is 0.0623598, the max value of r is
	 // 0.0623598 (=1/(4*4))
	 // clipper::data32::F_phi friedel = fc_from_model[hri];
	 // friedel.friedel();
	 float r = sqrt(hri.invresolsq());
	 float x = 2*M_PI*x_gl*r;
	 float y = gsl_sf_bessel_J0(x);
	 if (0) 
	    std::cout << r_k_i <<  " for " << hri.hkl().format() << " y: " << y << std::endl;
	 
	 std::complex<float> t = fc_from_model[hri];
	 std::complex<double> yt(y*t.real(), y*t.imag());
	 std::complex<double> yt_friedel(y*t.real(), -y*t.imag());
	 if (0) 
	    std::cout << r_k_i <<  " for " << hri.hkl().format() << " x: " << x << " r: " << r
		      << " adding " << yt << std::endl;
	 T[r_k_i] += yt;          // sum now (then multiply after looping over reflection list)
	 T[r_k_i] += yt_friedel; 
      }
      float func_r_k = 1;
      T[r_k_i] *= w_k * func_r_k * r_k *r_k;
   }

   std::cout << "---- T(k) table ----- " << std::endl;
   for (float r_k_i=1; r_k_i<=N; r_k_i++) {
      std::cout << "   rki " << r_k_i << "  val: " << T[r_k_i] << std::endl;
   }

   clipper::HKL_data< clipper::datatypes::F_phi<double> > A_data = map_fphidata; 
   
   // run over the "map fragment" reflection list
   // 
   clipper::HKL_info::HKL_reference_index hri;
   for (hri = map_fphidata.first(); !hri.last(); hri.next()) {
      std::complex<double> sum(0,0);
      for (float r_k_i=1; r_k_i<=N; r_k_i++) {

	 float x_gl = radius_range*gl.abscissa(r_k_i)*0.5 + radius_mid;
	 float r = sqrt(hri.invresolsq());
	 float x = 2*M_PI*x_gl*r;
	 float y = gsl_sf_bessel_J0(x);
	 std::complex<double> prod(T[r_k_i].real() * y, T[r_k_i].imag() * y);
	 sum += prod;
      }
      A_data[hri] = std::complex<double>(map_fphidata[hri]) * sum;
      if (0) 
	 std::cout << "   A_data: " << hri.hkl().format()
		   << " f: " << A_data[hri].f() << " phi: " << A_data[hri].phi()
		   << " sum_abs: " << std::abs(sum) << " sum_arg: " << std::arg(sum)
		   << "\n";
   }

   // scale down A_data to "sensible" numbers:
   // for (hri = map_fphidata.first(); !hri.last(); hri.next()) A_data[hri].f() *= 0.0000001;


   // compare phases:
   if (0)
      for (hri = map_fphidata.first(); !hri.last(); hri.next())
	 std::cout << "  A_phase vs map phase " << hri.hkl().format() << " "
		   << clipper::Util::rad2d(map_fphidata[hri].phi()) << " "
		   << clipper::Util::rad2d(      A_data[hri].phi()) << std::endl;
   


   clipper::Xmap<float> A_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());

   // debug
   { 
      clipper::HKL_data< clipper::datatypes::F_phi<double> > A_inv_data = A_data;
      for (hri = map_fphidata.first(); !hri.last(); hri.next()) A_inv_data[hri].phi() = -A_data[hri].phi();
      clipper::Xmap<float> A_inv_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
      A_inv_map.fft_from(A_inv_data);
      clipper::CCP4MAPfile mapout;
      mapout.open_write("A_inv.map");
      mapout.export_xmap(A_inv_map);
      mapout.close_write();
   }
   
   A_map.fft_from(A_data);
   
   // mean_and_variance<float> mv =
   float v = map_density_distribution(A_map, 40, false).variance;
   float rmsd = sqrt(v);
   peak_search ps(A_map);
   float n_sigma = 3;
   std::vector<std::pair<clipper::Coord_grid, float> > peaks = ps.get_peak_grid_points(A_map, n_sigma);
   std::cout << "==== peaks ==== " << std::endl;
   for (unsigned int ipeak=0; ipeak<peaks.size(); ipeak++)
      std::cout << ipeak << " "
		<< peaks[ipeak].second << " "
		<< peaks[ipeak].first.format() << "  "
		<< peaks[ipeak].second/rmsd
		<< std::endl;

   std::vector<clipper::Coord_orth> vp;
   for (unsigned int ipeak=0; ipeak<peaks.size(); ipeak++)
      vp.push_back(peaks[ipeak].first.coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()));
   mmdb::Manager *mol = create_mmdbmanager_from_points(vp);
   

   clipper::CCP4MAPfile mapout;
   mapout.open_write("A.map");
   mapout.export_xmap(A_map);
   mapout.close_write();
   
   mapout.open_write("map_fragment.map");
   mapout.export_xmap(xmap);
   mapout.close_write();
   
}

void
coot::util::emma::overlap_simple(const clipper::Xmap<float> &xmap) const {

   // sfs from map
   clipper::Resolution reso(3.0);
   clipper::HKL_info hkl_info(xmap.spacegroup(), xmap.cell(), reso);
   hkl_info.generate_hkl_list();
   clipper::HKL_data< clipper::datatypes::F_phi<double> > map_fphidata(hkl_info);
   xmap.fft_to(map_fphidata, clipper::Xmap_base::Normal);

   clipper::Range<double>    fc_range = fc_from_model.invresolsq_range();
   clipper::Range<double> map_f_range = map_fphidata.invresolsq_range();
   std::cout << "fc from model resolution ranges " << fc_range.min() << " " << fc_range.max() << std::endl;
   std::cout << "fc from model resolution ranges " << 1/sqrt(fc_range.min()) << " " << 1/sqrt(fc_range.max())
	     << std::endl;
   std::cout << "SFs from map  resolution ranges " << map_f_range.min() << " " << map_f_range.max() << std::endl;
   std::cout << "SFs from map  resolution ranges " << 1/sqrt(map_f_range.min()) << " " << 1/sqrt(map_f_range.max())
	     << std::endl;

   clipper::HKL_data< clipper::datatypes::F_phi<double> > A_data = map_fphidata; 
   
   double radius_of_integration = 15; // needs sqrt(1/r)
   clipper::HKL_info::HKL_reference_index hri_model;
   clipper::HKL_info::HKL_reference_index hri_map;
   for (hri_map = map_fphidata.first(); !hri_map.last(); hri_map.next()) {

      std::complex<double> sum_F_Sii(0,0);
      
      for (hri_model = fc_from_model.first(); !hri_model.last(); hri_model.next()) {
	 // integrate function from 0 to a in a_step
	 double a_step = 1;
	 double Sii = 0.0; // inner_integration_sum
	 for (double va=0.5; va<=radius_of_integration;va+=a_step) {
	    Sii += f(hri_model, hri_map, va) * a_step; // a_step is strip width (for strip area)
	 }
	 std::complex<float> tf = fc_from_model[hri_model];  // transfer
	 std::complex<double> f_model(tf.real(), tf.imag()); // transfer
	 std::complex<double> f_model_fridel(tf.real(), -tf.imag());
	 sum_F_Sii += f_model        * std::complex<double>(Sii, 0);
	 sum_F_Sii += f_model_fridel * std::complex<double>(Sii, 0);
      }
      std::complex<double> t = std::complex<double>(A_data[hri_map]) * sum_F_Sii;
      A_data[hri_map] = t;
   }


   // ------------------------- output the A map and analyse -----------------------------
   
   clipper::Xmap<float> A_map(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   A_map.fft_from(A_data);
   // mean_and_variance<float> mv =
   float v = map_density_distribution(A_map, 40, false).variance;
   float rmsd = sqrt(v);
   peak_search ps(A_map);
   float n_sigma = 3;
   std::vector<std::pair<clipper::Coord_grid, float> > peaks = ps.get_peak_grid_points(A_map, n_sigma);
   std::cout << "==== peaks ==== " << std::endl;
   for (unsigned int ipeak=0; ipeak<peaks.size(); ipeak++)
      std::cout << ipeak << " "
		<< peaks[ipeak].second << " "
		<< peaks[ipeak].first.format() << "  "
		<< peaks[ipeak].second/rmsd
		<< std::endl;

   std::vector<clipper::Coord_orth> vp;
   for (unsigned int ipeak=0; ipeak<peaks.size(); ipeak++)
      vp.push_back(peaks[ipeak].first.coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()));
   mmdb::Manager *mol = create_mmdbmanager_from_points(vp);
   

   clipper::CCP4MAPfile mapout;
   mapout.open_write("A.map");
   mapout.export_xmap(A_map);
   mapout.close_write();
   
}

double
coot::util::emma::f(const clipper::HKL_info::HKL_reference_index &hri_model,
		    const clipper::HKL_info::HKL_reference_index &hri_map,
		    double va) const{

   double x_a = 2*M_PI*va*sqrt(hri_model.invresolsq());
   double x_b = 2*M_PI*va*sqrt(hri_map.invresolsq());
   double j0_a = gsl_sf_bessel_J0(x_a);
   double j0_b = gsl_sf_bessel_J0(x_b);
   double d = j0_a * j0_b * va * va;
   return d;
} 


void
coot::util::emma::test() const {

   std::cout << "--------------------- start test -------------" << std::endl;
   std::cout << "--------------------- done test -------------" << std::endl;
}


// put atoms on grid, calc sfs then average by invresolsq
// Hmm.. This is the wrong way I think
std::vector<std::pair<double, double> >
coot::util::emma::spherically_averaged_FT() {

   std::vector<std::pair<double, double> > vp;
   int n_bins = 10;

   return vp;
}

std::vector<std::pair<double, double> >
coot::util::spherically_averaged_molecule(const atom_selection_container_t &asc,
					  float angstroms_per_bin) {

   std::vector<std::pair<double, double> > vp;

   std::pair<clipper::Coord_orth, clipper::Coord_orth> e = extents(asc.mol);
   std::pair<bool, clipper::Coord_orth> cc = centre_of_molecule(asc.mol);
   if (!cc.first) return vp;
   clipper::Coord_orth c = cc.second;

   double diag_len_sqrd = (e.second - e.first).lengthsq();
   double diag_len = sqrt(diag_len_sqrd);
   double radius_max = 0.5 * diag_len;

   int n_bins = static_cast<int>(radius_max/angstroms_per_bin) + 1;
   vp.resize(n_bins);
   for (int ibin=0; ibin<n_bins; ibin++) {
      vp[ibin].first  = (static_cast<float>(ibin) + 0.5) * angstroms_per_bin;
      vp[ibin].second = 0.0;
   }


   for (int iat=0; iat<asc.n_selected_atoms; iat++) {
      mmdb::Atom *at = asc.atom_selection[iat];
      clipper::Coord_orth co = coot::co(at);
      float dist = std::sqrt((co - c).lengthsq());
      int bin_id = dist / angstroms_per_bin;
      if (bin_id >= n_bins) {
	 std::cout << "bin error! " << std::endl;
      } else {
	 vp[bin_id].second += 1.0;
      }
   }

   return vp;
}

std::vector<coot::util::phitheta>
coot::util::make_phi_thetas(unsigned int n_pts) {

   std::vector<phitheta> v;

   // ideally I should push these around so that they are equidistant
   //
   double recip = 1.0/static_cast<double> (RAND_MAX);
   for (std::size_t i=0; i<n_pts; i++) {
      double theta = 2 * M_PI * random() * recip; // longitude
      double phi = acos(2.0*random()*recip-1.0);  // latitude
      v.push_back(std::pair<double, double>(phi, theta));
   }

   return v;
}


float
coot::util::average_of_sample_map_at_sphere_points(clipper::Coord_orth &centre,
						   float radius,
						   const std::vector<coot::util::phitheta> &phi_thetas,
						   clipper::Xmap<float> &xmap) {
   float r = 0.0;

   double sum = 0.0;
   for (std::size_t i=0; i<phi_thetas.size(); i++) {
      const double &phi   = phi_thetas[i].first;
      const double &theta = phi_thetas[i].second;
      clipper::Coord_orth pt(radius * cos(theta) * sin(phi),
			     radius * sin(theta) * sin(phi),
			     radius * cos(phi));
      //std::cout << "phi " << phi << " theta " << theta  << " pt: " << pt.format() << std::endl;
      pt += centre;
      sum += density_at_point_by_linear_interpolation(xmap, pt);
   }
   r = sum/static_cast<double>(phi_thetas.size());

   return r;

}
