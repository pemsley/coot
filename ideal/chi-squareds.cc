/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
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

// #define ANALYSE_REFINEMENT_TIMING

#include <string.h> // for strcmp

#ifdef ANALYSE_REFINEMENT_TIMING
#include <sys/time.h>
#endif // ANALYSE_REFINEMENT_TIMING

// we don't want to compile anything if we don't have gsl

#include <fstream>
#include <algorithm> // for sort
#include <stdexcept>
#include <iomanip>

#ifdef HAVE_CXX_THREAD
#include <thread>
#include <chrono>
#endif // HAVE_CXX_THREAD

#include "utils/split-indices.hh"
#include "geometry/mol-utils.hh"
#include "geometry/main-chain.hh"
#include "simple-restraint.hh"

//
#include "coot-utils/coot-coord-extras.hh"  // is_nucleotide_by_dict

std::vector<coot::refinement_lights_info_t>
coot::restraints_container_t::chi_squareds(std::string title, const gsl_vector *v, bool print_table_flag) const {

   bool print_summary = print_table_flag;

   if (!v) {
      std::cout << "ERROR:: oops null v in chi_squareds()" << std::endl;
   }
   if (verbose_geometry_reporting == QUIET) print_summary = false;
   
   std::vector<refinement_lights_info_t> lights_vec;
   int n_bond_restraints = 0; 
   int n_angle_restraints = 0; 
   int n_torsion_restraints = 0; 
   int n_plane_restraints = 0;
   int n_parallel_plane_restraints = 0;
   int n_improper_dihedral_restraints = 0;
   int n_non_bonded_restraints = 0;
   int n_chiral_volumes = 0;
   int n_rama_restraints = 0;
   int n_start_pos_restraints = 0;
   int n_target_pos_restraints = 0;
   int n_geman_mcclure_distance = 0;
   int n_trans_peptide_restraints = 0;

   double bond_distortion = 0; 
   double gm_distortion = 0; 
   double angle_distortion = 0; 
   double torsion_distortion = 0; 
   double plane_distortion = 0; 
   double parallel_planes_distortion = 0; 
   double non_bonded_distortion = 0;
   double chiral_vol_distortion = 0;
   double rama_distortion = 0;
   double start_pos_distortion = 0;
   double target_pos_distortion = 0;
   double trans_peptide_distortion = 0;
   double improper_dihedral_distortion = 0;

   // const be gone :-) (I only do this because we are interfacing with a
   // GSL function. Ideally params should be const void * for most of it's usages.
   //
   void *params = const_cast<void *>(reinterpret_cast<const void *>(this));
   std::pair<int, double> dist_max_bonds(0,0);
   std::pair<int, double> dist_max_angles(0,0);
   std::pair<int, double> dist_max_planes(0,0);
   std::pair<int, double> dist_max_nbc(0,0);
   std::map<std::string, refinement_lights_info_t::the_worst_t> baddies;
   std::map<std::string, refinement_lights_info_t::the_worst_t>::iterator baddies_iterator;

   for (int i=0; i<size(); i++) {
      {
	 const simple_restraint &restraint = restraints_vec[i];
	 if (restraints_usage_flag & BONDS_MASK) {
	    if (restraint.restraint_type == BOND_RESTRAINT) {
	       n_bond_restraints++;
	       // 	    bond_distortion += distortion_score_bond(restraint, v);
	       double dist = distortion_score_bond(restraint, v);
	       bond_distortion += dist;
	       if (dist > dist_max_bonds.second) {
		  dist_max_bonds.first = i;
		  dist_max_bonds.second = dist;
	       }
	       baddies["Bonds"].update_if_worse(dist, i);
	    }
	 }

	 if (restraints_usage_flag & GEMAN_MCCLURE_DISTANCE_MASK) {
	    if (restraint.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	       n_geman_mcclure_distance++;
	       double d = distortion_score_geman_mcclure_distance(restraint, v, geman_mcclure_alpha);
	       gm_distortion += d;
	       baddies["GemanMcClure"].update_if_worse(d, i);
	    }
	 }

	 if (restraints_usage_flag & ANGLES_MASK) { // 2: angles
	    if (restraint.restraint_type == coot::ANGLE_RESTRAINT) {
	       n_angle_restraints++;
	       double dist = coot::distortion_score_angle(restraint, v);
	       angle_distortion += dist;
	       if (dist > dist_max_angles.second) {
		  dist_max_angles.first = i;
		  dist_max_angles.second = dist;
	       }
	       baddies["Angles"].update_if_worse(dist, i);
	    }
	 }

         if (restraints_usage_flag & IMPROPER_DIHEDRALS_MASK) {
            if (restraint.restraint_type == coot::IMPROPER_DIHEDRAL_RESTRAINT) {
               n_improper_dihedral_restraints++;
               double dist = coot::distortion_score_improper_dihedral(restraint, v);
               improper_dihedral_distortion += dist;
               baddies["ImproperDihedrals"].update_if_worse(dist, i);
            }
         }

	 if (restraints_usage_flag & TORSIONS_MASK) { // 4: torsions
	    if (restraint.restraint_type == coot::TORSION_RESTRAINT) {
	       try {
		  double dist = coot::distortion_score_torsion(i, restraint, v);
		  torsion_distortion += dist;
		  n_torsion_restraints++;
		  baddies["Torsions"].update_if_worse(dist, i);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING:: caught runtime_error torsion " << rte.what() << std::endl;
	       }
	    }
	 }

	 if (restraints_usage_flag & TRANS_PEPTIDE_MASK) {
	    if (restraint.restraint_type == TRANS_PEPTIDE_RESTRAINT) {
	       try {
		  double dist = distortion_score_trans_peptide(i, restraint, v);
		  trans_peptide_distortion += dist;
		  n_trans_peptide_restraints++;
		  baddies["Trans_peptide"].update_if_worse(dist, i);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING:: caught runtime_error trans-pep " << rte.what() << std::endl;
	       }
	    }
	 }

	 if (restraints_usage_flag & PLANES_MASK) { // 8: planes
	    if (restraint.restraint_type == coot::PLANE_RESTRAINT) {
	       n_plane_restraints++;
	       double dist = coot::distortion_score_plane(restraint, v);
	       plane_distortion += dist;
	       if (dist > dist_max_planes.second) {
		  dist_max_planes.first = i;
		  dist_max_planes.second = dist;
	       }
	       baddies["Planes"].update_if_worse(dist, i);
	       if (false) {  // debugging plane restraints.
		  std::cout << " plane distortion " << i << " " 
			    << coot::distortion_score_plane(restraint, v) << " " 
			    << restraint;
		  for (unsigned int jj = 0; jj<restraint.plane_atom_index.size(); jj+=3) { 
		     std::cout << "\n                                ";
		     unsigned int idx = restraint.plane_atom_index[jj].first;
		     std::cout << idx << " " << coot::atom_spec_t(atom[idx]);
		     if ((jj+1) < restraint.plane_atom_index.size()) { 
			unsigned int idx_1 = restraint.plane_atom_index[jj+1].first;
			std::cout << " " << idx_1 << " " << coot::atom_spec_t(atom[idx_1]);
		     }
		     if ((jj+2) < restraint.plane_atom_index.size()) { 
			unsigned int idx_2 = restraint.plane_atom_index[jj+1].first;
			std::cout << " " << idx_2 << " " << coot::atom_spec_t(atom[idx_2]);
		     }
		  }
		  std::cout << std::endl;
	       }
	    }
	 }

         if (restraints_usage_flag & PARALLEL_PLANES_MASK) {
            if (restraint.restraint_type == coot::PARALLEL_PLANES_RESTRAINT) {
               n_parallel_plane_restraints++;
               double dist = coot::distortion_score_parallel_planes(restraint, v);
               parallel_planes_distortion += dist;
               baddies["Parallel Planes"].update_if_worse(dist, i);
               if (true) {
                  std::cout << "parallel plane " << i << " " << restraint << " " << dist << std::endl;
               }
            }
         }

	 if (restraints_usage_flag & coot::NON_BONDED_MASK) { 
	    if ( restraint.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) { 
	       n_non_bonded_restraints++;
	       double dist = coot::distortion_score_non_bonded_contact(restraint, lennard_jones_epsilon, v);
	       non_bonded_distortion += dist;
	       if (dist > dist_max_nbc.second) {
		  dist_max_nbc.first = i;
		  dist_max_nbc.second = dist;
	       }
	       baddies["NonBonded"].update_if_worse(dist, i);
	    }
	 }

	 if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) { 
	    if ( restraint.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
	       n_chiral_volumes++;
	       double dist = coot::distortion_score_chiral_volume(restraint, v);
	       chiral_vol_distortion += dist;
	       baddies["Chirals"].update_if_worse(dist, i);
	    }
	 }

	 if (restraints_usage_flag & coot::RAMA_PLOT_MASK) {
	    if (restraint.restraint_type == coot::RAMACHANDRAN_RESTRAINT) {
	       n_rama_restraints++;
	       if (rama_type == restraints_container_t::RAMA_TYPE_ZO) {
		  double dd = distortion_score_rama( restraint, v, ZO_Rama(), get_rama_plot_weight());
		  rama_distortion += dd;
		  baddies["Rama"].update_if_worse(dd, i);

#if 0	// needs indexing fixup	  
		  if (false) { // debugging rama baddie update
		     baddies_iterator = baddies.find("Rama");
		     if (baddies_iterator != baddies.end()) {
			const refinement_lights_info_t::the_worst_t &w = baddies_iterator->second;
			const simple_restraint &baddie_restraint = restraints_vec[w.restraints_index];
			std::cout << "Running rama worst baddie: w.restraints_index " << w.restraints_index
				  << " w.value " << w.value
				  << " distortion " << baddie_restraint.format(atom, w.value)
				  << std::endl;
		     }
		  }
#endif		  

	       } else {

                  // type is RAMA_TYPE_LOGRAMA

                  double w = get_rama_plot_weight();
		  double dd = distortion_score_rama(restraint, v, lograma, w);
		  rama_distortion += dd;
		  baddies["Rama"].update_if_worse(dd, i);
	       }
	       if (false) {
                  double w = get_rama_plot_weight();
		  double d1 = distortion_score_rama(restraint, v, LogRama(), w);
		  double d2 = coot::distortion_score_rama(restraint, v, ZO_Rama(), get_rama_plot_weight());
		  std::cout << "distortion-comparision logramas " << d1 << " zo " << d2 << std::endl;
	       }
	    }
	 }

	 if (restraint.restraint_type == coot::TARGET_POS_RESTRAINT) {
	    n_target_pos_restraints++;
	    double dist = coot::distortion_score_target_pos(restraint, log_cosh_target_distance_scale_factor, v);
	    target_pos_distortion += dist;
	    baddies["Target_pos"].update_if_worse(dist, i);
	 }

	 if (restraint.restraint_type == coot::START_POS_RESTRAINT) {
	    n_start_pos_restraints++;
	    double dist = distortion_score_start_pos(restraint, params, v);
	    start_pos_distortion += dist;
	    baddies["StartPositions"].update_if_worse(dist, i);
	 }
      }
   }

   std::string r = "";

   r += title;
   r += "\n";
   std::setprecision(3);
   if (print_summary)
      std::cout << "    " << title << std::endl;
   if (n_bond_restraints == 0) {
      if (print_summary)
	 std::cout << "bonds:      N/A " << std::endl;
   } else {
      double bd = bond_distortion/double(n_bond_restraints);
      double sbd = 0.0;
      if (bd > 0)
	 sbd = sqrt(bd);
      if (print_summary)
	 std::cout << "bonds:      " << sbd << std::endl;
      r += "   bonds:  ";
      r += coot::util::float_to_string_using_dec_pl(sbd, 3);
      r += "\n";
      std::string s = "Bonds:  ";
      s += coot::util::float_to_string_using_dec_pl(sbd, 3);
      coot::refinement_lights_info_t rl("Bonds", s, sbd);
      baddies_iterator = baddies.find("Bonds");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_angle_restraints == 0) {
      if (print_summary)
	 std::cout << "angles:     N/A " << std::endl;
   } else {
      double ad = angle_distortion/double(n_angle_restraints);
      double sad = 0.0;
      if (ad > 0.0)
	 sad = sqrt(ad);
      if (print_summary)
	 std::cout << "angles:     " << sad << std::endl;
      r += "   angles: ";
      r += coot::util::float_to_string_using_dec_pl(sad, 3);
      r += "\n";
      std::string s = "Angles: ";
      s += coot::util::float_to_string_using_dec_pl(sad, 3);
      coot::refinement_lights_info_t rl("Angles", s, sad);
      baddies_iterator = baddies.find("Angles");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_improper_dihedral_restraints == 0) {
      if (print_summary)
	 std::cout << "improper-dihedrals: N/A " << std::endl;
   } else {
      double idd = improper_dihedral_distortion/static_cast<double>(n_improper_dihedral_restraints);
      if (print_summary)
         std::cout << "improper-dihedrals: " << idd << std::endl;
      r += "   improper-dihedrals: ";
      r += util::float_to_string_using_dec_pl(idd, 3);
      r += " ";
      r += util::int_to_string(n_improper_dihedral_restraints);
      r += "\n";
      // add worst baddie handling here
   }
   if (n_torsion_restraints == 0) {
      if (print_summary)
	 std::cout << "torsions:   N/A " << std::endl;
   } else {
      double td = torsion_distortion/double(n_torsion_restraints);
      double std = 0.0;
      if (td > 0.0)
	 std = sqrt(td);
      if (print_summary)
	 std::cout << "torsions:   " << std << std::endl;
      r += "   torsions: ";
      r += coot::util::float_to_string_using_dec_pl(std, 3);
      r += "\n";
      std::string s = "Torsions: ";
      s += coot::util::float_to_string_using_dec_pl(std, 3);
      coot::refinement_lights_info_t rl("Torsions", s, std);
      baddies_iterator = baddies.find("Torsions");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_trans_peptide_restraints == 0) {
      if (print_summary)
	 std::cout << "trans-peptide: N/A" << std::endl;
   } else {
      double td = trans_peptide_distortion/double(n_trans_peptide_restraints);
      if (print_summary)
	 std::cout << "trans-peptide: " << td << " (non-sqrt)" << std::endl;
      r += "   trans-peptide: ";
      r += coot::util::float_to_string_using_dec_pl(td, 3);
      r += "\n";
      std::string s = "Trans_peptide: ";
      s += coot::util::float_to_string_using_dec_pl(td, 3);
      coot::refinement_lights_info_t rl("Trans_peptide", s, td);
      baddies_iterator = baddies.find("Trans_peptide");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
	 
   }
   if (n_plane_restraints == 0) {
      if (print_summary)
	 std::cout << "planes:     N/A " << std::endl;
   } else {
      double pd = plane_distortion/static_cast<double>(n_plane_restraints);
      double spd = 0.0;
      if (pd > 0.0)
	 spd = sqrt(pd);
      if (print_summary)
	 std::cout << "planes:     " << spd << " from " << n_plane_restraints << " restraints " << std::endl;
      r += "   planes: ";
      r += coot::util::float_to_string_using_dec_pl(spd, 3);
      r += "\n";
      std::string s = "Planes: ";
      s += coot::util::float_to_string_using_dec_pl(spd, 3);
      coot::refinement_lights_info_t rl("Planes", s, spd);
      baddies_iterator = baddies.find("Planes");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_parallel_plane_restraints == 0) {
      if (print_summary)
	 std::cout << "parallel planes:     N/A " << std::endl;
   } else {
      double ppd = parallel_planes_distortion/static_cast<double>(n_parallel_plane_restraints);
      double sppd = 0.0;
      if (ppd > 0.0)
	 sppd = sqrt(ppd);
      if (print_summary)
	 std::cout << "parallel planes: " << sppd << " from " << n_parallel_plane_restraints
                   << " restraints " << std::endl;
      r += "   parallel planes: ";
      r += coot::util::float_to_string_using_dec_pl(sppd, 3);
      r += "\n";
      std::string s = "Parallel Planes: ";
      s += coot::util::float_to_string_using_dec_pl(sppd, 3);
      coot::refinement_lights_info_t rl("Parallel Planes", s, sppd);
      baddies_iterator = baddies.find("Parallel Planes");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_non_bonded_restraints == 0) {
      if (print_summary)
	 std::cout << "non-bonded: N/A " << std::endl;
   } else {
      double nbd = non_bonded_distortion/double(n_non_bonded_restraints);
      double snbd = 0.0;
      if (nbd > 0.0)
	 snbd = sqrt(nbd);
      if (print_summary)
	 std::cout << "non-bonded: " << nbd << std::endl;
      r += "   non-bonded: ";
      r += coot::util::float_to_string_using_dec_pl(snbd, 3);
      r += "\n";
      std::string s = "Non-bonded: ";
      s += coot::util::float_to_string_using_dec_pl(snbd, 3);
      coot::refinement_lights_info_t rl("Non-bonded", s, snbd);
      baddies_iterator = baddies.find("NonBonded");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_chiral_volumes == 0) { 
      if (print_summary)
	 std::cout << "chiral vol: N/A " << std::endl;
   } else {
      double cd = chiral_vol_distortion/double(n_chiral_volumes);
      double scd = 0.0;
      if (cd > 0.0)
	 scd = sqrt(cd);
      if (print_summary)
	 std::cout << "chiral vol: " << scd << std::endl;
      r += "   chirals: ";
      r += coot::util::float_to_string_using_dec_pl(scd, 3);
      std::string s = "Chirals: ";
      s += coot::util::float_to_string_using_dec_pl(scd, 3);
      coot::refinement_lights_info_t rl("Chirals", s, scd);
      baddies_iterator = baddies.find("Chirals");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_rama_restraints == 0) { 
      if (print_summary)
	 std::cout << "rama plot:  N/A " << std::endl;
   } else {
      double rd = rama_distortion/double(n_rama_restraints);

      if (print_summary)
	 std::cout << "rama plot:  " << rd << " " << n_rama_restraints << std::endl;

      r += "   rama plot: ";
      r += util::float_to_string_using_dec_pl(rd, 3);
      std::string s = "Rama Plot: ";
      s += util::float_to_string_using_dec_pl(rd, 3);
      refinement_lights_info_t rli("Rama", s, rd);
      baddies_iterator = baddies.find("Rama");
      if (baddies_iterator != baddies.end()) {
	 rli.worst_baddie = baddies_iterator->second;
	 const simple_restraint &baddie_restraint = restraints_vec[rli.worst_baddie.restraints_index];
	 if (print_summary)
	    std::cout << "rama worst baddie: index " << rli.worst_baddie.restraints_index
		      << " distortion " << baddie_restraint.format(atom, rli.worst_baddie.value)
		      << std::endl;
      }
      if (rama_type == RAMA_TYPE_ZO)
	 rli.rama_type = RAMA_TYPE_ZO;
      
      lights_vec.push_back(rli);
   }
   if (n_start_pos_restraints == 0) {
      if (print_summary)
	 std::cout << "start_pos:  N/A " << std::endl;
   } else {
      double spd = start_pos_distortion/double(n_start_pos_restraints);
      double sspd = 0.0;
      if (spd > 0.0)
	 sspd = sqrt(spd);
      if (print_summary)
	 std::cout << "start_pos:  " << sspd << std::endl;
      r += "startpos:  ";
      r += util::float_to_string_using_dec_pl(sspd, 3);
      r += "\n";
      std::string s = "Start pos: ";
      s += util::float_to_string_using_dec_pl(sspd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Start_pos", s, sspd));
   }

   if (n_target_pos_restraints == 0) {
      if (print_summary)
	 std::cout << "TargetPos:  N/A " << std::endl;
   } else {
      double tpd = target_pos_distortion/double(n_target_pos_restraints);
      if (print_summary) 
	 std::cout << "target_pos: " << tpd  << " (non-sqrt)" << std::endl;
      r += "targetpos: ";
      r += util::float_to_string_using_dec_pl(tpd, 3);
      r += "\n";
      std::string s = "Target pos:";
      s += util::float_to_string_using_dec_pl(tpd, 3);
      baddies_iterator = baddies.find("Target_pos");
      coot::refinement_lights_info_t rl("Target_pos", s, tpd);
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }
   if (n_geman_mcclure_distance == 0) {
      if (print_summary)
	 std::cout << "GemanMcCl:  N/A " << std::endl;
   } else {
      double spd = gm_distortion/double(n_geman_mcclure_distance);
      double sspd = 0.0;
      if (spd > 0.0)
	 sspd = sqrt(spd);
      if (print_summary)
	 std::cout << "GemanMcCl:  " << sspd << " from " << n_geman_mcclure_distance << " distances"
		   << std::endl;
      r += "GemanMcCl:  ";
      r += util::float_to_string_using_dec_pl(sspd, 3);
      r += "\n";
      std::string s = "GemanMcCl: ";
      s += util::float_to_string_using_dec_pl(sspd, 3);
      coot::refinement_lights_info_t rl("GemanMcCl", s, sspd);
      baddies_iterator = baddies.find("GemanMcClure");
      if (baddies_iterator != baddies.end())
	 rl.worst_baddie = baddies_iterator->second;
      lights_vec.push_back(rl);
   }

   // more about baddies:
   if (false) {
      std::cout << std::endl;
      for (std::size_t i=0; i<lights_vec.size(); i++) {
	 const refinement_lights_info_t &rl = lights_vec[i];
	 if (rl.worst_baddie.is_set) {
	    const simple_restraint &baddie_restraint = restraints_vec[rl.worst_baddie.restraints_index];
	    std::cout << " worst baddie of type " << std::setw(13) << rl.name << " "
		      << rl.worst_baddie.value << " "
		      << std::setw(4) << rl.worst_baddie.restraints_index << " "
		      << std::setprecision(8)
		      << baddie_restraint.format(atom, rl.worst_baddie.value) << std::endl;
	 } else {
	    std::cout << "worst baddie not set " << rl.name << std::endl;
	 }
      }
   }
   return lights_vec;
} 

