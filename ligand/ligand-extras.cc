/* ligand/ligand-extras.cc 
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
 * Copyright 2014 by Medical Research Council
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
 * 02110-1301, USA.
 */

// This file was initially created as the USA started its (2003) attack on Iraq.
//
 
#include <stdio.h> // for snprintf

#include <fstream>
#include <list>

#include "clipper/core/xmap.h"
#include "clipper/core/map_utils.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_interp.h"

#include "utils/coot-utils.hh"
#include "ligand.hh" // has mmdb-manager because mmdb::PPAtom is part
                     // of mask_map interface.

#include <mmdb2/mmdb_coormngr.h> // for GetMassCenter

#include "coot-utils/xmap-stats.hh"
#include "coot-utils/peak-search.hh"


// uses index to cluster array
// c.f. version that is passed a coot::map_point_cluster
short int
coot::ligand::cluster_is_possible_water(int i) const {

   // very crude.  Fixme.
   //
   float vol = xmap_pristine.cell().volume();
   float ngrid =
      xmap_pristine.grid_sampling().nu() *
      xmap_pristine.grid_sampling().nv() *
      xmap_pristine.grid_sampling().nw();

   float grid_vol = vol/ngrid;
   float n_grid_lim = water_molecule_volume/grid_vol; // 15 is 4/3 Pi r^3: with r=1.53

//    std::cout << "n_grid_lim is " << n_grid_lim
// 	     << " map_grid size: " << cluster[i].map_grid.size() << std::endl;
      
   if (cluster[i].map_grid.size() < n_grid_lim) {
      return 1;
   } else {
//       std::cout << "Big cluster has " << cluster[i].map_grid.size()
// 		<< " grid points with n_grid_lim " << n_grid_lim << " single grid volume "
// 		<< grid_vol << "A^3" << std::endl;
      return 0;
   }
}

// uses coot::map_point_cluster
// c.f. version that is passed an index to cluster array
short int
coot::ligand::cluster_is_possible_water(const coot::map_point_cluster &mpc) const {

   float vol = xmap_pristine.cell().volume();
   float ngrid =
      xmap_pristine.grid_sampling().nu() *
      xmap_pristine.grid_sampling().nv() *
      xmap_pristine.grid_sampling().nw();

   float grid_vol = vol/ngrid;
   float n_grid_lim = water_molecule_volume/grid_vol; // 15 is 4/3 Pi r^3: with r=1.53

   if (mpc.map_grid.size() < n_grid_lim) {
      return 1;
   } else {
      return 0;
   }
}

// Things above this volume limit are too bit to be waters.  What is that limit?
// 
double
coot::ligand::possible_water_volume_limit() const {

   return water_molecule_volume;
} 


// In the case of waters, the search map is the masked map.
//
clipper::Coord_orth
coot::ligand::move_atom_to_peak(const clipper::Coord_orth &a,
				const clipper::Xmap<float> &search_map) const {

   clipper::Coord_orth pos = a;
   float dv; // density value
   clipper::Grad_frac<float> grad_frac;
   clipper::Grad_map<float>  grad_map;
   clipper::Curv_map<float> curv_map;
   float shift_len = 1.0; // Angstroms
   int n_cycle = 0; 
   int n_cycle_max = 500; 

   while ((n_cycle < n_cycle_max) && (shift_len > 0.001)) { // Angstroms

      clipper::Coord_frac a_cf = pos.coord_frac(search_map.cell());
      //      std::cout << "getting grad_map:  at " << pos.format() << std::endl;
      // std::cout << "getting grad_map:  at " << a_cf.format() << std::endl;
      clipper::Coord_map  a_cm = a_cf.coord_map(search_map.grid_sampling());
      clipper::Interp_cubic::interp_grad(search_map, a_cm, dv, grad_map);

      // Can't get this or something like it to compile...  Ask Kevin.
      // search_map.interp_curv<clipper::Interp_cubic>(a_cm, dv, grad_map, curv_map);

      grad_frac = grad_map.grad_frac(search_map.grid_sampling());
      clipper::Grad_orth<float> grad_o = grad_frac.grad_orth(search_map.cell());
      //
      clipper::Coord_orth shift(gradient_scale*0.8*grad_o.dx(),
				gradient_scale*0.8*grad_o.dy(),
				gradient_scale*0.8*grad_o.dz());

      shift_len = sqrt(shift.lengthsq());

// debugging:
//       std::cout << n_cycle << "   " << shift_len << "  " << shift.x() << "  "
//	<< shift.y() << "  " << shift.z() << std::endl;


      pos += shift;
      n_cycle++;
   }

   if (n_cycle == n_cycle_max) {
      std::cout << "WARNING:: refinement failure" <<  std::endl;
      std::cout << "          start pos: " << a.format()   << std::endl;
      std::cout << "          final pos: " << pos.format() << std::endl << std::endl;
   }
   return pos;

}

#include "geometry/residue-and-atom-specs.hh"

float
coot::ligand::density_at_point(const clipper::Coord_orth &atom_pos,
			       const clipper::Xmap<float> &search_map) const {

   clipper::Coord_frac atom_pos_frc = atom_pos.coord_frac(search_map.cell());
   float dv = search_map.interp<clipper::Interp_cubic>(atom_pos_frc);
   return dv;
}

#include "coot-utils/fib-sphere.hh"
#include "analysis/stats.hh"

std::pair<float, float>
coot::ligand::mean_and_variance_where_the_atoms_are(mmdb::Manager *mol) const {

   std::pair<float, float> r(0,0);
   unsigned int n_test_points = 100;
   std::vector<clipper::Coord_orth> test_points;

   unsigned int n_molecule_atoms = 0;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      const clipper::Xmap<float> &xmap = xmap_pristine;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) {
	    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    std::string rn(residue_p->GetResName());
	    if (rn != "HOH") {
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       for (int iat=0; iat<n_atoms; iat++) {
		  mmdb::Atom *at = residue_p->GetAtom(iat);
		  if (! at->isTer()) {
		     std::string ele = at->element;
		     if (ele != " H")
			n_molecule_atoms++;
		  }
	       }
	    }
	 }
      }

      if (n_molecule_atoms > n_test_points) {
	 float rmi = 1.0f/float(RAND_MAX);
	 float crit_val = static_cast<float>(n_test_points)/static_cast<float>(n_molecule_atoms);
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       std::string rn(residue_p->GetResName());
	       if (rn != "HOH") {
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  for (int iat=0; iat<n_atoms; iat++) {
		     mmdb::Atom *at = residue_p->GetAtom(iat);
		     if (! at->isTer()) {
			std::string ele = at->element;
			if (ele != " H") {
			   float f = coot::util::random() * rmi;
			   if (f < crit_val) {
			      clipper::Coord_orth c(at->x, at->y, at->z);
			      test_points.push_back(c);
			   }
			}
		     }
		  }
	       }
	    }
	 }
      } else {
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       std::string rn(residue_p->GetResName());
	       if (rn != "HOH") {
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  for (int iat=0; iat<n_atoms; iat++) {
		     mmdb::Atom *at = residue_p->GetAtom(iat);
		     if (! at->isTer()) {
			std::string ele = at->element;
			if (ele != " H") {
			   clipper::Coord_orth c(at->x, at->y, at->z);
			   test_points.push_back(c);
			}
		     }
		  }
	       }
	    }
	 }
      }

      if (! test_points.empty()) {
	 coot::stats::single s;
	 for (std::size_t i=0; i<test_points.size(); i++)
	    s.add(density_at_point(test_points[i], xmap));
	 float m = s.mean();
	 float sd = sqrt(s.variance());
	 return std::pair<float, float>(m,sd);
      }
   }
   return r;
}


coot::ligand::spherical_density_score_t
coot::ligand::spherical_density_score(const clipper::Coord_orth &a,
				      float mean_density_of_other_atoms) const {

   double step = 0.4;

   const clipper::Xmap<float> &search_map = xmap_pristine;

   std::vector<int> n_samples = { 0, 30, 80, 150 };
   int n_total = 30 + 80 + 150;

   std::vector<float> var(4,0);

   float sum_for_scaling = 0.0;

   for (int istep=1; istep<=3; istep++) {
      std::vector<clipper::Coord_orth> sphere_points = fibonacci_sphere(n_samples[istep]);
      double sum_sq = 0;
      double sum = 0;

      for (int i=0; i<n_samples[istep]; i++) {
         clipper::Coord_orth pos = a + step * sphere_points[i];
         float dv = density_at_point(pos, search_map);
         sum    += dv;
         sum_sq += dv * dv;
         sum_for_scaling += dv;
      }

      float mean = sum/static_cast<float>(n_samples[istep]);
      var[istep] = sum_sq/static_cast<float>(n_samples[istep]) - mean * mean;

   }
   float overall_mean = sum_for_scaling/static_cast<float>(n_total);
   float non_spherical = 0.0;
   for (int istep=1; istep<=3; istep++)
      non_spherical += 0.333 * sqrt(var[istep])/mean_density_of_other_atoms;
   // if (non_spherical < 0.0) non_spherical = 0.0;
   // if (non_spherical > 1.0) non_spherical = 1.0;
   float sphericalness = 1.0 - non_spherical;

   float dp = density_at_point(a, search_map);
   return spherical_density_score_t(dp, non_spherical);

}

short int
coot::ligand::has_sphericalish_density(const clipper::Coord_orth &a,
				       const clipper::Xmap<float> &search_map) const {

   // var_limit = 0.07;

   short int iret = 0;
   double step = 0.3;
   clipper::Coord_orth pos;
   double dv[6];
   double var[3];
   double peak_height = density_at_point(a, search_map);
   for (int i=1; i<4; i++) {
      double i_f = double(i);
      pos = clipper::Coord_orth(a.x() + i_f*step, a.y(), a.z());
      dv[0] = density_at_point(pos, search_map);
      pos = clipper::Coord_orth(a.x() - i_f*step, a.y(), a.z());
      dv[1] = density_at_point(pos, search_map);

      pos = clipper::Coord_orth(a.x(), a.y() + i_f*step, a.z());
      dv[2] = density_at_point(pos, search_map);
      pos = clipper::Coord_orth(a.x(), a.y() - i_f*step, a.z());
      dv[3] = density_at_point(pos, search_map);

      pos = clipper::Coord_orth(a.x(), a.y(), a.z() + i_f*step);
      dv[4] = density_at_point(pos, search_map);
      pos = clipper::Coord_orth(a.x(), a.y(), a.z() - i_f*step);
      dv[5] = density_at_point(pos, search_map);

      double sum_sq = 0;
      double sum = 0;
      for (int j=0; j<6; j++) {
         sum    += dv[j];
         sum_sq += dv[j]*dv[j];
      }
      var[i-1] = sum_sq/6.0 - sum*sum/36.0;
   }

   double total_var = var[0] + var[1] + var[2];
//    std::cout << "variance test vars: " << var[0] << " " << var[1]
// 	     << " " << var[2] << std::endl;
//    std::cout << "variance test : (peak height " << peak_height << ") "
// 	     << total_var/(peak_height*peak_height)
// 	      << " vs. " << var_limit << std::endl;
   if (total_var/(peak_height*peak_height) < var_limit) {
      iret = 1;
   } else {
      iret = 0;
   }

   return iret;
}

// Is close to protein, H-bonding distance to an O or N atom?
//
// Return OK_GOLDILOCKS // OK, neither too close nor too far from protein
//        TOO_CLOSE     // too close to a protein atom
//        TOO_FAR       // too far from protein atom
//        WATER_STATUS_UNKNOWN // initial value
//
short int
coot::ligand::water_pos_is_chemically_sensible(clipper::Coord_orth water_centre) const {

   short int iret = coot::ligand::WATER_STATUS_UNKNOWN;
   float len;

   for (unsigned int ifrag=0; ifrag<protein_atoms.fragments.size(); ifrag++) {
      for (int ires=protein_atoms.fragments[ifrag].min_res_no();
	   ires<=protein_atoms.fragments[ifrag].max_residue_number();
	   ires++) {
	 for (unsigned int iatom=0; iatom<protein_atoms.fragments[ifrag][ires].atoms.size(); iatom++) {
	    
	    if (protein_atoms[ifrag][ires][iatom].element == " N" ||
		protein_atoms[ifrag][ires][iatom].element == " O") {
	       
	       if (protein_atoms[ifrag][ires].name != "WAT" &&
		   protein_atoms[ifrag][ires].name != "HOH") {
		  
		  len = clipper::Coord_orth::length(protein_atoms[ifrag][ires][iatom].pos, water_centre);
// 		  std::cout << protein_atoms[ifrag][ires][iatom].pos.format() << " "
// 			    << len << " " << water_to_protein_distance_lim_min
// 			    << " " << water_to_protein_distance_lim_max<< std::endl;
		  if (len < water_to_protein_distance_lim_min) {
		     iret = coot::ligand::TOO_CLOSE;
		     break;
		  } else {
		     if (len < water_to_protein_distance_lim_max) {
// 			std::cout << "found H-bond atom " << ires << " "
// 				  << protein_atoms[ifrag][ires][iatom].name 
// 				  << " at " << protein_atoms[ifrag][ires][iatom].pos.format() << std::endl;
			iret = coot::ligand::OK_GOLDILOCKS;
		     }
		  }
	       }
	    }
	 }
	 if (iret == coot::ligand::TOO_CLOSE)
	    break;
      }
      if (iret == coot::ligand::TOO_CLOSE)
	 break;
   }
   return iret;
}


// Return OK_GOLDILOCKS // OK, neither too close nor too far from protein
//        TOO_CLOSE     // too close to a protein atom
//        TOO_FAR       // too far from protein atom
//        WATER_STATUS_UNKNOWN // initial value
// 
short int
coot::ligand::water_pos_is_chemically_sensible(const clipper::Coord_orth &water_centre,
					       const std::vector<clipper::Coord_orth> &extra_sites) const {

   // this allows for a water to be accepted because it is
   // OK_GOLDILOCKS to a water, even though it is too far from protein
   // atoms.
   // 
   short int iret = water_pos_is_chemically_sensible(water_centre);
   if (iret == coot::ligand::WATER_STATUS_UNKNOWN || iret == coot::ligand::TOO_FAR) {
      double len;
      double min_length = 9999.9;
      for (unsigned int i=0; i<extra_sites.size(); i++) {
	 len = clipper::Coord_orth::length(water_centre, extra_sites[i]);
	 if (len < min_length)
	    min_length = len;
      }
      if (min_length < water_to_protein_distance_lim_max) {
	 if (min_length > water_to_protein_distance_lim_min) {
	    iret = OK_GOLDILOCKS;
	 }
      }
   }

   // This rejects water_centre if it is too close to anythin in extra_sites
   if (iret == OK_GOLDILOCKS) {
      double len;
      double min_length = 9999.9;
      for (unsigned int i=0; i<extra_sites.size(); i++) {
	 len = clipper::Coord_orth::length(water_centre, extra_sites[i]);
	 if (len < water_to_protein_distance_lim_min) {
	    iret = coot::ligand::TOO_CLOSE;
	 }
      }
   }
   return iret;
}


void
coot::ligand::write_waters(const std::vector<clipper::Coord_orth> &water_list,
			   const std::string &file_name) const {

   short int separate_residues = 1;
   std::cout << "writing "
	     << water_list.size() << " water atoms to ligand-waters.pdb"
	     << std::endl;
   std::string chain_id = protein_atoms.unused_chain_id("W"); // pass the prefered chain id.
   minimol::molecule mol(water_list, "HOH", " O  ", chain_id);
   mol.write_file(file_name, default_b_factor); 
} 

void
coot::ligand::water_fit(float sigma_cutoff, int n_cycle) {

   clipper::Coord_orth new_centre;
   std::vector<clipper::Coord_orth> water_list;
   std::vector<clipper::Coord_orth> this_round_water_list;
   std::vector<clipper::Coord_orth> raw_water_list;

   short int found_waters_prev_round_flag = 1;

   if (xmap_masked_stats.first == 0) { 
      clipper::Map_stats stats(xmap_cluster);
      xmap_masked_stats.first = 1;
      xmap_masked_stats.second.first  = stats.mean();
      xmap_masked_stats.second.second = stats.std_dev();
   }

   water_list = water_fit_internal(sigma_cutoff, n_cycle);
   
   std::cout << "INFO:: found " << water_list.size()
	     << " waters in water fitting"
      
	     << std::endl;
   std::cout.flush();
   std::string ch = protein_atoms.unused_chain_id("W");
   coot::minimol::molecule mol(water_list, "HOH", " O  ", ch);
   
   mol.set_cell(xmap_cluster.cell());
   std::string spg(xmap_cluster.spacegroup().descr().symbol_hm());
   mol.set_spacegroup(spg);
   water_molecule = mol;
}

std::vector<clipper::Coord_orth> 
coot::ligand::water_fit_internal(float sigma_cutoff, int n_cycle) {

   std::vector<clipper::Coord_orth> water_list;
   std::vector<clipper::Coord_orth> this_round_water_list;
   std::vector<clipper::Coord_orth> raw_water_list;
   std::vector<std::pair<clipper::Coord_orth, double> > blobs;

   // xmap_cluster is fine at start
   // output_map(xmap_cluster, "xmap_cluster_start_water_fit.map");

   if (xmap_masked_stats.first == 0) {
      std::cout << "ERROR: xmap_masked_stats not set" << std::endl;
   } else {
//       std::cout << "DEBUG:: at start of water_fit_internal: map stats: sigma: "
// 		<< xmap_masked_stats.second.second << std::endl;
      
      std::vector <clipper::Coord_orth> sampled_protein_coords =
	 make_sample_protein_coords();
      
      float z_cutoff = sigma_cutoff;
      // std::cout << "DEBUG:: sigma_cutoff in round is " << z_cutoff
      // << std::endl;
      n_clusters = 0;
      cluster.clear();
      find_clusters_internal(z_cutoff, sampled_protein_coords); // fill cluster
      std::cout << "-------------------------------------------------------------"
		<< std::endl;
      // std::cout << "DEBUG:: found " << cluster.size()
      // << " clusters at " << z_cutoff << "z cut." << std::endl;

      std::list<coot::map_point_cluster> cluster_list;
      // convert from a vector to a list:
      for (unsigned int ic=0; ic<cluster.size(); ic++) {
	 cluster_list.push_back(cluster[ic]);
      }
      
	 
      for(int iround = 0; iround < n_cycle; iround++) {
	 
// 	 // useful debugging
// 	 std::string mapfilename = "xmap_cluster_start_water_fit-";
// 	 mapfilename += coot::util::int_to_string(iround);
// 	 mapfilename += ".map";
// 	 output_map(xmap_cluster, mapfilename);
	 
	 std::list<coot::map_point_cluster>::iterator it;
// 	 std::vector<std::list<coot::map_point_cluster>::const_iterator> iterator_remove_list;
	 std::vector<std::list<coot::map_point_cluster>::iterator> iterator_remove_list;
	 //	 std::cout << "DEBUG:: round " << iround << " cluster list size: "
	 // << cluster_list.size() << "\n";
	 for(it=cluster_list.begin(); it!=cluster_list.end(); it++) {
	    
	    clipper::Coord_orth cl_centre(it->eigenvectors_and_centre.trn());
	    if ((do_cluster_size_check_flag && cluster_is_possible_water(*it))
		|| do_cluster_size_check_flag == 0) {
	       
	       clipper::Coord_orth new_centre = move_atom_to_peak(cl_centre, xmap_cluster);
	       short int chem_sensible = water_pos_is_chemically_sensible(new_centre, water_list);
	       if ((do_chemically_sensible_test_flag &&
		    (chem_sensible == coot::ligand::OK_GOLDILOCKS))
		   || (do_chemically_sensible_test_flag == 0)) {
		  water_list.push_back(new_centre);
		  iterator_remove_list.push_back(it);
// 	       } else {
// 		  std::cout << "INFO:: site at " << new_centre.format()
// 			    << " is not chemically sensible\n";
		  // but it may be next round.
	       } 
	       
	    } else { 
	       if (iround == 0) {
		  std::cout << "INFO:: cluster at " << cl_centre.format()
			    << " is too big to be water\n";
		  double cl_vol = it->volume(xmap_pristine);
		  std::pair<clipper::Coord_orth, double> p(cl_centre, cl_vol);
		  blobs.push_back(p);
		  iterator_remove_list.push_back(it);
	       }
	    }
	 }
	 for (unsigned int irl=0; irl<iterator_remove_list.size(); irl++)
	    // cluster_list.remove(*(iterator_remove_list[irl]));
	    cluster_list.erase(iterator_remove_list[irl]);
      }
   }
   if (write_raw_waters)
      write_waters(raw_water_list, "raw-ligand-waters-peaksearch-results.pdb");

   // user can get to these via big_blobs() member function
   keep_blobs = blobs;
   return water_list; 
}
		
void
coot::ligand::flood() {

   clipper::Coord_orth new_centre;
   std::vector<clipper::Coord_orth> water_list;
   std::vector<clipper::Coord_orth> this_round_water_list;
   std::vector<clipper::Coord_orth> raw_water_list;

   std::cout << "positioning waters in " << cluster.size() 
	     << " density clusters " << std::endl;
   
   std::vector<std::pair<clipper::Coord_orth, double> > blobs; // is this used in this function!?
   short int found_waters_prev_round_flag = 1;
   float density_cut_off = 0.2;

   std::vector <clipper::Coord_orth> sampled_protein_coords =
      make_sample_protein_coords();
      
   if (xmap_masked_stats.first == 0) { 
      clipper::Map_stats stats(xmap_cluster);
      xmap_masked_stats.first = 1;
      xmap_masked_stats.second.first  = stats.mean();
      xmap_masked_stats.second.second = stats.std_dev();
   }

   if (xmap_masked_stats.first != 1) { 
      std::cout << "PROGRAMMER ERROR! find_clusters/find_clusters_water_flood first!\n";
      // exit(1); // maybe. // no exit() from libraries
      return;
   }

   int n_cycle = 20;
   for(int iround = 0; iround < n_cycle; iround++) {

      if (1) { 

	 // at start (iround = 0) high cut-off 
	 // at end   (iround = (n_cycle -1) low cut-ff
	 // 

	 double frac = float(iround+1)/float(n_cycle);
	 float z_cut_off = (1.0 - frac) * 4.0 + 1.0;
	 density_cut_off = z_cut_off*xmap_masked_stats.second.second; // n sigma, that is.

	 n_clusters = 0;
	 cluster.clear();
	 find_clusters_water_flood(z_cut_off, sampled_protein_coords);
   
	 this_round_water_list.clear();

	 std::cout << "INFO:: Density cut-off in round " << iround << " is "
		   << density_cut_off << std::endl;
	 for (unsigned int iclust=0; iclust<cluster.size(); iclust++) {
	    
	    // std::cout << "DEBUG:: cluster: " << i << " of " << cluster.size() << std::endl;

	    if (cluster[iclust].map_grid.size() > 0) { 
	       
	       clipper::Coord_orth cl_centre =
		  clipper::Coord_orth(cluster[iclust].eigenvectors_and_centre.trn());
	       raw_water_list.push_back(cl_centre);
	       
	       
	       if (density_at_point(cl_centre, xmap_masked) > density_cut_off) { 

		  this_round_water_list.push_back(cl_centre);
	       } else {
// 		  std::cout << "   Round " << iround << " reject cluster " << iclust
// 			    << " because " << density_at_point(cl_centre, xmap_masked)
// 			    << " < " << density_cut_off << std::endl;
	       } 
	    }
	 }

	 // The total waters are in water_list:
	 
	 //
	 std::cout << "Round " << iround << ": Adding " << this_round_water_list.size() 
		   << " new waters to list which contains " 
		   << water_list.size() << " waters already\n";
	 if (this_round_water_list.size() == 0)
	    found_waters_prev_round_flag = 0; // set flag for no more rounds then...
	 // 
	 for (unsigned int iw=0; iw<this_round_water_list.size(); iw++) { 
	    water_list.push_back(this_round_water_list[iw]);
	 } 

	 if ( (iround < (n_cycle-1)) && (this_round_water_list.size() > 0) ) {
	 
	    // mask new waters
	    for (unsigned int iw=0; iw<this_round_water_list.size(); iw++) { 
	       // change xmap_masked:
	       mask_around_coord(this_round_water_list[iw], map_atom_mask_radius, &xmap_masked);
	    }
	 }
      }
   }
   
   std::cout << "found " << water_list.size() << " waters in water fitting"
	     << std::endl;
   std::cout.flush();
   if (write_solutions)
      write_waters(water_list, "ligand-waters.pdb");
   // std::cout << "constructing a minimol molecule for waters " << std::endl;
   std::string ch = protein_atoms.unused_chain_id("W");
   // std::cout << "DEBUG:: Water fit: unused chain: " << ch << std::endl;
   coot::minimol::molecule mol(water_list, "HOH", " O  ", ch);

   // set space group and cell of mol here.
   // 
//    float cell[6];
//    cell[0] = xmap_masked.cell().a();
//    cell[1] = xmap_masked.cell().b();
//    cell[2] = xmap_masked.cell().c();
//    cell[3] = clipper::Util::rad2d(xmap_masked.cell().alpha());
//    cell[4] = clipper::Util::rad2d(xmap_masked.cell().beta());
//    cell[5] = clipper::Util::rad2d(xmap_masked.cell().gamma());
//    mol.set_cell(cell);


   mol.set_cell(xmap_cluster.cell());
   std::string spg(xmap_cluster.spacegroup().descr().symbol_hm());
   mol.set_spacegroup(spg);

   water_molecule = mol; 

// mol.write_file("coords.pdb");

   if (write_raw_waters)
      write_waters(raw_water_list, "raw-ligand-waters-peaksearch-results.pdb");

   // user can get to these via big_blobs() member function
   keep_blobs = blobs;
}

// debugging 
#include "clipper/ccp4/ccp4_map_io.h"

// find peaks that are in density greater than n_sigma
void
coot::ligand::flood2(float n_sigma) {

   bool ignore_pseudo_zeros = true; // because cryo-EM maps - maybe pass this?

   bool debug = false;

   int n_rounds = 20/map_atom_mask_radius; // 1.4 default.

   if (false) { // debugging hack
      n_rounds = 1;
      clipper::CCP4MAPfile mapout;
      mapout.open_write("flood2-masked.map");
      mapout.export_xmap(xmap_masked);
      mapout.close_write(); 
   }

   std::vector<clipper::Coord_orth> water_list;

   mean_and_variance<float> mv_start = map_density_distribution(xmap_masked, 40, false, ignore_pseudo_zeros);

   int n_added_waters=0;
   for (int iround=0; iround<n_rounds; iround++) {

      mean_and_variance<float> mv_this = map_density_distribution(xmap_masked, 40, false, ignore_pseudo_zeros);
      float n_sigma_crit = n_sigma * sqrt(mv_start.variance/mv_this.variance);

      if (debug)
	 std::cout << "flood2() iround = " << iround << " n_sigma_crit "
		   << n_sigma_crit << std::endl;

      coot::peak_search ps(xmap_masked);
      ps.set_max_closeness(0);
      std::vector<clipper::Coord_orth> peaks = ps.get_peaks_for_flooding(xmap_masked, n_sigma_crit);
      
      if (debug) {
	 for (unsigned int ipeak=0; ipeak<peaks.size(); ipeak++) {
	    float d = density_at_point(peaks[ipeak], xmap_masked);
	    std::cout << "iround = " << iround << " peak " << ipeak << " "
		      << peaks[ipeak].format() << " " << d << std::endl;
	 }
      }
      
      if (peaks.size() == 0) {
	 std::cout << "INFO:: No extra peaks: breaking on round "
		   << iround << " of " << n_rounds << std::endl;
	 break;
      }
      // mask new waters
      for (unsigned int iw=0; iw<peaks.size(); iw++) {
	 if (! close_to_another(peaks[iw], water_list, map_atom_mask_radius)) {
	    // change xmap_masked:
	    mask_around_coord(peaks[iw], map_atom_mask_radius, &xmap_masked);
	    water_list.push_back(peaks[iw]);
	    n_added_waters++;
	 }
      }
   }

   // move these waters to around protein:
   // 
   std::vector <clipper::Coord_orth> sampled_protein_coords = make_sample_protein_coords();
   std::cout << "DEBUG:: in flood2() sampled_protein_coords() size is "
	     << sampled_protein_coords.size() << std::endl;

   std::vector<clipper::Coord_orth> moved_waters;   
   if (sampled_protein_coords.size())
      moved_waters = move_waters_close_to_protein(water_list, sampled_protein_coords);
   else 
      moved_waters = water_list; // unomved waters

   std::cout << "INFO:: added " << n_added_waters << " waters to molecule\n";
   std::string ch = protein_atoms.unused_chain_id("W");
   // coot::minimol::molecule mol(water_list, "DUM", " DUM", ch);
   coot::minimol::molecule mol(moved_waters, "DUM", " DUM", ch);
   mol.set_cell(xmap_masked.cell());
   std::string spg(xmap_masked.spacegroup().descr().symbol_hm());
   mol.set_spacegroup(spg);
   water_molecule = mol;
}

std::vector<clipper::Coord_orth>
coot::ligand::move_waters_close_to_protein(const std::vector<clipper::Coord_orth> &water_list,
					   const std::vector<clipper::Coord_orth> &sampled_protein_coords) const {


   std::vector<clipper::Coord_orth> sample_set = sampled_protein_coords; // gets added to

   // We add the moved waters to sampled_protein_coords so that the
   // waters of this cluster are move likely to be in the same
   // symmetry as the water above it in the list that it might
   // otherwise have been.
   
   std::vector<clipper::Coord_orth> moved_waters = water_list;
   for (unsigned int i=0; i<water_list.size(); i++) { 
      moved_waters[i] = move_water_close_to_protein(water_list[i], sample_set);
      sample_set.push_back(moved_waters[i]);
   }
   return moved_waters;
}

clipper::Coord_orth
coot::ligand::move_water_close_to_protein(const clipper::Coord_orth &water_pos,
					  const std::vector<clipper::Coord_orth> &sampled_protein_coords) const {

   double best_dist = 999999999999999999999999999999.;
   clipper::Coord_orth moved_pos = water_pos;
   
   int n_sampled = sampled_protein_coords.size();
   if (n_sampled > 0) { 
      int n = xmap_pristine.spacegroup().num_symops();
      clipper::Coord_frac cell_shift; 
      clipper::Coord_orth t_point; 
      for (int isym=0; isym<n; isym++) {
	 for (int x_shift = -1; x_shift<2; x_shift++) { 
	    for (int y_shift = -1; y_shift<2; y_shift++) { 
	       for (int z_shift = -1; z_shift<2; z_shift++) {
		  cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift); 
		  clipper::RTop_orth orthop =
		     clipper::RTop_frac(xmap_pristine.spacegroup().symop(isym).rot(),
					xmap_pristine.spacegroup().symop(isym).trn() +
					cell_shift).rtop_orth(xmap_pristine.cell());
		  t_point = water_pos.transform(orthop); 
		  double t_dist = min_dist_to_protein(t_point, sampled_protein_coords);
		  if (t_dist < best_dist) {
		     moved_pos = t_point;
		     best_dist = t_dist;
		  }
	       }
	    }
	 }
      }
   }
   return moved_pos;
} 


bool
coot::ligand::close_to_another(const clipper::Coord_orth &p1,
			       const std::vector<clipper::Coord_orth> &ref,
			       const double &d_crit) const {

   // fast distance check
   bool is_close = 0; 
   for (unsigned int i=0; i<ref.size(); i++) {
      double d1 = ref[i].x() - p1.x();
      if (d1 < d_crit) {
	 double d2 = ref[i].y() - p1.y();
	 if (d2 < d_crit) {
	    double d3 = ref[i].z() - p1.z();
	    if (d3 < d_crit) {
	       if ((d1*d1 + d2*d2 +d3*d3) < d_crit*d_crit) { 
		  is_close = 1;
		  break;
	       }
	    }
	 }
      } 
   }
   return is_close;
}


// void
// coot::ligand::flood_and_fit(float n_sigma_in) { 

//    // now new peak_search code:
//    coot::peak_search ps(xmap_pristine, n_sigma_in);
//    std::vector<clipper::Coord_orth> ps_peaks = ps.get_peaks(xmap_pristine);
//    std::string ch = protein_atoms.unused_chain_id("W");
//    coot::minimol::molecule ps_mol(ps_peaks, 0, " OW1", ch);
//    std::string spg(xmap_masked.spacegroup().descr().symbol_hm());
//    ps_mol.set_spacegroup(spg);
//    ps_mol.set_cell(xmap_pristine.cell());
//    ps_mol.write_file("peaksearch-peaks.pdb");
//    coot::high_res hr(ps_mol);
//    hr.output_pdb("coords.pdb");

// } 


// Xmap_plus_bits
// coot::ligand::masked_map_plus_bits() const {

//    Xmap_plus_bits a;
//    clipper::Map_stats stats(xmap);

//    a.xmap = xmap; // masked map
//    a.mean = stats.mean();
//    a.std_dev = stats.std_dev();
//    a.name = "map masked by protein";

//    return a;

// } 

const clipper::Xmap<float> &
coot::ligand::masked_map() const {
   return xmap_masked;
}


coot::minimol::molecule
coot::ligand::get_solution(unsigned int isolution, unsigned int iclust) const {

   coot::minimol::molecule empty;
   unsigned int n_final_ligands = final_ligand.size();
   if (iclust < n_final_ligands) {
      if (isolution < final_ligand[iclust].size())
	 return final_ligand[iclust][isolution].first;
   } else {
      std::cout << "Error in get_solution: iclust is " << iclust
		<< " but final_size is " << final_ligand.size()
		<< " with inital ligand size "
		<< int(initial_ligand.size()) << std::endl;
   }
   return empty;
}
