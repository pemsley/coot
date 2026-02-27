/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

// #define ANALYSE_REFINEMENT_TIMING

#include <string.h> // for strcmp

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

// #include "mmdb.h" // for printing of mmdb::Atom pointers as info not raw
                     // pointers.  Removed. Too much (linking issues in)
                     // Makefile pain.

#include "compat/coot-sysdep.h"

#include "utils/logging.hh"
extern logging logger;


// We need to fill restraints_vec (which is a vector of
// simple_restraint) using the coordinates () and the dictionary of
// restraints, protein_geometry geom.
//
// The plan is to get a list of residues, and for each of those
// residues, look in geom for a residue of that type and if found,
// fill restraints_vec appropriately.
int
coot::restraints_container_t::make_restraints(int imol,
					      const coot::protein_geometry &geom,
					      coot::restraint_usage_Flags flags_in, 
					      bool do_residue_internal_torsions,
					      bool do_trans_peptide_restraints,
					      float rama_plot_target_weight,
					      bool do_rama_plot_restraints,
					      bool do_auto_helix_restraints,
					      bool do_auto_strand_restraints,
					      bool do_auto_h_bond_restraints,
					      coot::pseudo_restraint_bond_type sec_struct_pseudo_bonds,
					      bool do_link_restraints,
					      bool do_flank_restraints) {


#if 1

   make_restraints_ng(imol, geom, flags_in, do_residue_internal_torsions, do_trans_peptide_restraints,
		      rama_plot_target_weight, do_rama_plot_restraints,
		      do_auto_helix_restraints,
                      do_auto_strand_restraints,
                      do_auto_h_bond_restraints,
		      sec_struct_pseudo_bonds, do_link_restraints, do_flank_restraints);

   return size();
#endif

   // if a peptide is trans, add a restraint to penalize non-trans configuration
   // (currently a torsion restraint on peptide w of 180)
   //

   if (false)
      std::cout << "debug:: make_restraints() called with flags " << flags_in << std::endl;

   // debugging SRS inclusion.
   if (false) {

      std::cout << "------- debug:: here in make_restraints() do_trans_peptide_restraints is "
		<< do_trans_peptide_restraints << std::endl;
      std::cout << "----- make restraints() called with geom of size : " << geom.size() << std::endl;
      std::cout << "    geom ref pointer " << &geom << std::endl;
   }

   restraints_usage_flag = flags_in; // also set in minimize() and geometric_distortions()

   if (false) { // debugging, by forcing the restraints type

      restraints_usage_flag = BONDS_ANGLES_AND_CHIRALS; // looks good

      restraints_usage_flag = BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS; // looks good!

      restraints_usage_flag = BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_TRANS_PEPTIDE_RESTRAINTS;

      // trans peptide dfs look fine, perhaps the presence of trans peptide restraints
      // make other restraints go bad?

      // perhaps it's the chirals?
      // Hmm! seems plausible
      // restraints_usage_flag = BONDS_ANGLES_PLANES_NON_BONDED_AND_TRANS_PEPTIDE_RESTRAINTS;
      //
      restraints_usage_flag = BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_TRANS_PEPTIDE_RESTRAINTS;

      // restraints_usage_flag = GEMAN_MCCLURE_DISTANCE_RESTRAINTS;
      // restraints_usage_flag = NO_GEOMETRY_RESTRAINTS;

      // restraints_usage_flag = BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      // restraints_usage_flag = BONDS_ANGLES_TORSIONS_NON_BONDED_CHIRALS_AND_TRANS_PEPTIDE_RESTRAINTS;

      restraints_usage_flag = TRANS_PEPTIDE_RESTRAINT;

      restraints_usage_flag = JUST_RAMAS;
   }

   if (n_atoms > 0) {

      mark_OXT(geom);
      make_monomer_restraints(imol, geom, do_residue_internal_torsions);

      bool do_link_restraints_internal = true;
      bool do_flank_restraints_internal = true;

      if (! from_residue_vector) {
	 if (istart_res == iend_res)
	    do_link_restraints_internal = false;
	 if (! istart_minus_flag && !iend_plus_flag)
	    do_flank_restraints_internal = false;
      }

      if (! do_link_restraints)
	 do_link_restraints_internal = false;
      if (! do_flank_restraints)
	 do_flank_restraints_internal = false;

      rama_plot_weight = rama_plot_target_weight;

      // sets bonded_pairs_container (note that this doesn't make bonded pairs for
      // residues that are not in the given set of residues) i.e. no bonds
      // between a CYS and a CYS that is in the residue set.
      if (do_link_restraints_internal) {
	 make_link_restraints(geom, do_rama_plot_restraints, do_trans_peptide_restraints);
      }

      if (false)
	 std::cout << "after make_link_restraints() bonded_pairs_container has size "
		   << bonded_pairs_container.size() << std::endl;

      // don't do torsions, ramas maybe.
      coot::bonded_pair_container_t bpc;

      if (do_flank_restraints_internal)
	 bpc = make_flanking_atoms_restraints(geom,
					      do_rama_plot_restraints,
					      do_trans_peptide_restraints);
      int iret_prev = size();

      if (sec_struct_pseudo_bonds == coot::HELIX_PSEUDO_BONDS) {
	 make_helix_pseudo_bond_restraints();
      } 
      if (sec_struct_pseudo_bonds == coot::STRAND_PSEUDO_BONDS) {
	 make_strand_pseudo_bond_restraints();
      }

      if (do_auto_helix_restraints)
	 make_helix_pseudo_bond_restraints_from_res_vec_auto();

      if (restraints_usage_flag & coot::NON_BONDED_MASK) {
	 if ((iret_prev > 0) || are_all_one_atom_residues) {
	    reduced_angle_info_container_t ai(restraints_vec);
	    int n_nbcr = make_non_bonded_contact_restraints(imol, bpc, ai, geom);
	    if (verbose_geometry_reporting != QUIET)
	       // std::cout << "INFO:: make_restraints(): made " << n_nbcr
	       //          << " non-bonded restraints\n";
	       logger.log(log_t::INFO, "make_restraints(): made", n_nbcr, "non-bonded restraints");
	 }
      }
      make_restraint_types_index_limits();
   }
   if (false)
      std::cout << "returning from make_restraints() and restraints_usage_flag is "
	        << restraints_usage_flag << std::endl;

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

   make_df_restraints_indices();
   make_distortion_electron_density_ranges();
#endif //  HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

   // debug
   // info();
   return size();
}


void
coot::restraints_container_t::make_restraint_types_index_limits() {

   // 20170628 why not restraints_limits_target_pos also?
   //          needs investigation

   unsigned int unset = 9999999;
   restraints_limits_bonds = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_angles = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_torsions = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_chirals = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_planes =  std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_parallel_planes =  std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_non_bonded_contacts = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_geman_mclure = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_start_pos = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_trans_peptide = std::pair<unsigned int, unsigned int> (unset,0);

   for (unsigned int i=0; i<restraints_vec.size(); i++) { // needs dual indexing
      {
	 const simple_restraint &restraint = restraints_vec[i];
	 if (restraint.restraint_type == coot::BOND_RESTRAINT) {
	    if (restraints_limits_bonds.first == unset)
	       restraints_limits_bonds.first = i;
	    if (i > restraints_limits_bonds.second)
	       restraints_limits_bonds.second = i;
	 }
	 if (restraint.restraint_type == coot::ANGLE_RESTRAINT) {
	    if (restraints_limits_angles.first == unset)
	       restraints_limits_angles.first = i;
	    if (i > restraints_limits_angles.second)
	       restraints_limits_angles.second = i;
	 }
	 if (restraint.restraint_type == coot::TORSION_RESTRAINT) {
	    if (restraints_limits_torsions.first == unset)
	       restraints_limits_torsions.first = i;
	    if (i > restraints_limits_torsions.second)
	       restraints_limits_torsions.second = i;
	 }
	 if (restraint.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	    if (restraints_limits_chirals.first == unset)
	       restraints_limits_chirals.first = i;
	    if (i > restraints_limits_chirals.second)
	       restraints_limits_chirals.second = i;
	 }
	 if (restraint.restraint_type == coot::PLANE_RESTRAINT) {
	    if (restraints_limits_planes.first == unset)
	       restraints_limits_planes.first = i;
	    if (i > restraints_limits_planes.second)
	       restraints_limits_planes.second = i;
	 }
	 if (restraint.restraint_type == coot::PARALLEL_PLANES_RESTRAINT) {
	    if (restraints_limits_parallel_planes.first == unset)
	       restraints_limits_parallel_planes.first = i;
	    if (i > restraints_limits_parallel_planes.second)
	       restraints_limits_parallel_planes.second = i;
	 }
	 if (restraint.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    if (restraints_limits_non_bonded_contacts.first == unset)
	       restraints_limits_non_bonded_contacts.first = i;
	    if (i > restraints_limits_non_bonded_contacts.second)
	       restraints_limits_non_bonded_contacts.second = i;
	 }
	 if (restraint.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    if (restraints_limits_geman_mclure.first == unset)
	       restraints_limits_geman_mclure.first = i;
	    if (i > restraints_limits_geman_mclure.second)
	       restraints_limits_geman_mclure.second = i;
	 }
	 if (restraint.restraint_type == coot::TRANS_PEPTIDE_RESTRAINT) {
	    if (restraints_limits_trans_peptide.first == unset) {
	       restraints_limits_trans_peptide.first = i;
	    }
	    if (i > restraints_limits_trans_peptide.second)
	       restraints_limits_trans_peptide.second = i;
	 }
	 if (restraint.restraint_type == coot::START_POS_RESTRAINT) {
	    if (restraints_limits_start_pos.first == unset)
	       restraints_limits_start_pos.first = i;
	    if (i > restraints_limits_start_pos.second)
	       restraints_limits_start_pos.second = i;
	 }
      }
   }

   // now check for unsets
   if (restraints_limits_bonds.first  == unset)   restraints_limits_bonds.first = 0;
   if (restraints_limits_angles.first == unset)   restraints_limits_angles.first = 0;
   if (restraints_limits_torsions.first == unset) restraints_limits_torsions.first = 0;
   if (restraints_limits_chirals.first == unset)  restraints_limits_chirals.first = 0;
   if (restraints_limits_planes.first == unset)   restraints_limits_planes.first = 0;
   if (restraints_limits_parallel_planes.first == unset) restraints_limits_parallel_planes.first = 0;
   if (restraints_limits_non_bonded_contacts.first == unset) restraints_limits_non_bonded_contacts.first = 0;
   if (restraints_limits_geman_mclure.first == unset) restraints_limits_geman_mclure.first = 0;
   if (restraints_limits_start_pos.first == unset) restraints_limits_start_pos.first = 0;

   if (false) {
      std::cout << "restraints limits bonds "
	        << restraints_limits_bonds.first << " " << restraints_limits_bonds.second << std::endl;
      std::cout << "restraints limits angles "
	        << restraints_limits_angles.first << " " << restraints_limits_angles.second << std::endl;
      std::cout << "restraints limits torsions "
	        << restraints_limits_torsions.first << " " << restraints_limits_torsions.second << std::endl;
      std::cout << "restraints limits chirals "
	        << restraints_limits_chirals.first << " " << restraints_limits_chirals.second << std::endl;
      std::cout << "restraints limits planes "
	        << restraints_limits_planes.first << " " << restraints_limits_planes.second << std::endl;
      std::cout << "restraints limits nbc "
	        << restraints_limits_non_bonded_contacts.first << " " << restraints_limits_non_bonded_contacts.second
	        << std::endl;
   }

}


void
coot::restraints_container_t::make_helix_pseudo_bond_restraints() {

   // somewhat hacky
   if (from_residue_vector) {
      make_helix_pseudo_bond_restraints_from_res_vec();
      return;
   }


   // This method of making pseudo bonds relies on the residue range
   // being continuous in sequence number (seqNum) and no insertion
   // codes messing up the number scheme.  If these are not the case
   // then we will make bonds between the wrongs atoms (of the wrong
   // residues)... a more sophisticated algorithm would be needed.
   // 
   // The method ignores residues with alt confs.
   //
   float pseudo_bond_esd = 0.04; // just a guess
   int selHnd = mol->NewSelection();
   int nSelResidues;
   mmdb::PPResidue SelResidue;
   mmdb::PPAtom res_1_atoms = NULL;
   mmdb::PPAtom res_2_atoms = NULL;
   int n_res_1_atoms;
   int n_res_2_atoms;
   int index1 = -1; 
   int index2 = -1; 
   mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"",  // altLocs
		mmdb::SKEY_NEW // selection key
		);
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 0) {
      for (int i=4; i<nSelResidues; i++) {
         // nN -> (n-4)O 2.91
         // nN -> (n-3)O 3.18
         // nO -> (n+3)N 3.18  // symmetric.  No need to specify both forwards
         // nO -> (n+4)N 2.91  // and backwards directions.
	 SelResidue[i]->GetAtomTable(res_1_atoms, n_res_1_atoms);
	 for (int iat1=0; iat1<n_res_1_atoms; iat1++) {
	    std::string at_1_name(res_1_atoms[iat1]->name);

	    if (at_1_name == " N  ") {
	       mmdb::Residue *contact_res = SelResidue[i-4];
	       if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 4)) {
		  contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " O  ") {
			res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			add(BOND_RESTRAINT, index1, index2, fixed_flags,
			    2.91, pseudo_bond_esd, 1.2);
			std::cout << "Helix Bond restraint (" << res_1_atoms[iat1]->name << " "
				  << res_1_atoms[iat1]->GetSeqNum() << ") to ("
				  << res_2_atoms[iat2]->name << " "
				  << res_2_atoms[iat2]->GetSeqNum() << ") 2.91" << std::endl;
		     }
		  }
	       }

	       contact_res = SelResidue[i-3];
	       if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 3)) {
		  contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " O  ") {
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			add(BOND_RESTRAINT, index1, index2, fixed_flags,
			    3.18, pseudo_bond_esd, 1.2);
			std::cout << "Helix Bond restraint (" << res_1_atoms[iat1]->name << " "
				  << res_1_atoms[iat1]->GetSeqNum() << ") to ("
				  << res_2_atoms[iat2]->name << " "
				  << res_2_atoms[iat2]->GetSeqNum() << ") 3.18" << std::endl;
		     }
		  }
	       }
	       
	    } 
	 } 
      }
   }
   mol->DeleteSelection(selHnd);
}


#include "coot-utils/helix-like.hh"

void
coot::restraints_container_t::make_helix_pseudo_bond_restraints_from_res_vec_auto() {

   bool console_output_for_restraints_generation_timings = false;

   float pseudo_bond_esd = 0.1; // 0.05 was too tight (why?)
                                // 0.20 was too loose

   unsigned int n_helical_restraints = 0;

   auto tp_0 = std::chrono::high_resolution_clock::now();

   // somehow sometimes the residues of residue_vec can be null here

   std::vector<mmdb::Residue *> sorted_residues;
   sorted_residues.reserve(residues_vec.size());
   for (unsigned int i=0; i<residues_vec.size(); i++)
      if (residues_vec[i].second)
         sorted_residues.push_back(residues_vec[i].second);
   std::sort(sorted_residues.begin(), sorted_residues.end(), util::residues_sort_function);

   if (false)
      for (unsigned int i=0; i<sorted_residues.size(); i++)
         std::cout << "  sorted residue " << residue_spec_t(sorted_residues[i]) << std::endl;

   for (unsigned int i=0; i<sorted_residues.size(); i++) {

      if ((i+3) >= sorted_residues.size()) continue;

      if (sorted_residues[i]->GetChain() != sorted_residues[i+3]->GetChain()) continue;

      // test that these residues are about helical before adding helical restraints
      std::vector<mmdb::Residue *> test_helical_residues;
      // fill test_helical_residues with 4 residues in order.

      // unfortunately(?) test_helical_residues vector is "disconneted"
      // from the atom selection (below) for the restraints, so we do a
      // similar test for the chain there also.
      //
      mmdb::Residue *residue_0 = sorted_residues[i];
      for (unsigned int iir=0; iir<4; iir++) {
         if ((i+iir) < sorted_residues.size()) {
	    mmdb::Residue *residue_p = sorted_residues[i+iir];
	    // only add residues that are in the same chain as the first (0th) residue
	    if (residue_0) {
	       if (residue_p->GetChain() == residue_0->GetChain()) {
	          test_helical_residues.push_back(residue_p);
	       }
	    }
         }
      }

      // std::cout << "calling compare_to_helix() with test_helical_residues size "
      // << test_helical_residues.size() << std::endl;

      bool sane_residue_numbers = false;
      if (test_helical_residues.size() == 4)
         if ((test_helical_residues[0]->GetSeqNum()+3) == test_helical_residues[3]->GetSeqNum())
            sane_residue_numbers = true;

      // maybe we *do* want compare_to_helix() to be run on 5 residues?
      // i -> i+4 is the convention for alpha helical H-bonds, after all.
      //
      if (test_helical_residues.size() == 5)
         if ((test_helical_residues[0]->GetSeqNum()+4) == test_helical_residues[4]->GetSeqNum())
            sane_residue_numbers = true;

      helical_results_t hr = compare_to_helix(test_helical_residues); // tests for 4 residues

      if (false) // useful for debugging helical restraints
         std::cout << "DEBUG:: helix_result " << hr.is_alpha_helix_like << " "
                   << hr.sum_delta << " "
                   << residue_spec_t(sorted_residues[i]) << std::endl;

      if (hr.is_alpha_helix_like && sane_residue_numbers) {

         int index_1 = -1; // O
         int index_2 = -1; // N (n+4)
         int index_3 = -1; // N (n+3)
         mmdb::Atom *at_1 = 0;
         mmdb::Atom *at_2 = 0;
         mmdb::Atom *at_3 = 0;
         mmdb::Residue *residue_p = sorted_residues[i];
         mmdb::Atom **residue_atoms_1 = 0;
         mmdb::Atom **residue_atoms_2 = 0;
         mmdb::Atom **residue_atoms_3 = 0;
         int n_residue_atoms_1;
         int n_residue_atoms_2;
         int n_residue_atoms_3;
	 bool do_i_plus_4_also = true;
	 bool do_i_plus_3 = true;
	 if ((i+4) >= sorted_residues.size())
	    do_i_plus_4_also = false;

	 if (do_i_plus_4_also) {
	    if (sorted_residues[i+4]->GetChain() != residue_0->GetChain())
	       do_i_plus_4_also = false;
	 }
	 if (do_i_plus_4_also) {
	    if (sorted_residues[i+4]->GetSeqNum() != (residue_0->GetSeqNum() + 4))
	       do_i_plus_4_also = false;
	 }
	 if (sorted_residues[i+3]->GetChain() != residue_0->GetChain())
	    do_i_plus_3 = false;
	 if (do_i_plus_3) {
	    if (sorted_residues[i+3]->GetSeqNum() != (residue_0->GetSeqNum() + 3))
	       do_i_plus_3 = false;
	 }
         sorted_residues[i  ]->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
         sorted_residues[i+3]->GetAtomTable(residue_atoms_3, n_residue_atoms_3);
	 if (do_i_plus_4_also) {
	    sorted_residues[i+4]->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
	 }
         for (int iat=0; iat<n_residue_atoms_1; iat++) {
            std::string atom_name_1 = residue_atoms_1[iat]->GetAtomName();
            if (atom_name_1 == " O  ") {
               at_1 = residue_atoms_1[iat];
            }
         }
	 if (do_i_plus_3) {
	    for (int iat=0; iat<n_residue_atoms_3; iat++) {
	       std::string atom_name_3 = residue_atoms_3[iat]->GetAtomName();
	       if (atom_name_3 == " N  ") {
               at_3 = residue_atoms_3[iat];
	       }
	    }
	 }
	 if (do_i_plus_4_also) {
	    for (int iat=0; iat<n_residue_atoms_2; iat++) {
	       std::string atom_name_2 = residue_atoms_2[iat]->GetAtomName();
	       if (atom_name_2 == " N  ") {
		  at_2 = residue_atoms_2[iat];
	       }
	    }
	 }
         if (at_1 && at_2 && at_3) {
            at_1->GetUDData(udd_atom_index_handle, index_1);
            at_2->GetUDData(udd_atom_index_handle, index_2);
            at_3->GetUDData(udd_atom_index_handle, index_3);
            std::vector<bool> fixed_flags_1 = make_fixed_flags(index_1, index_2);
            std::vector<bool> fixed_flags_2 = make_fixed_flags(index_1, index_3);
            double ideal_dist_i_4 = 2.91;
            double ideal_dist_i_3 = 3.18;
            add(BOND_RESTRAINT, index_1, index_2, fixed_flags_1, ideal_dist_i_4, pseudo_bond_esd, 1.2);
            add(BOND_RESTRAINT, index_1, index_3, fixed_flags_2, ideal_dist_i_3, pseudo_bond_esd, 1.2);

            if (verbose_geometry_reporting != QUIET) {
               // std::cout << "INFO:: Alpha Helix Bond restraint ("
               //          << at_1->name << " " << at_1->GetSeqNum() << ") to ("
               //          << at_3->name << " " << at_3->GetSeqNum() << ") " << ideal_dist_i_3 << std::endl;
               logger.log(log_t::INFO, "Alpha Helix Bond restraint (" + std::string(at_1->name) + " " + std::to_string(at_1->GetSeqNum()) +
                          ") to (" + std::string(at_3->name) + " " + std::to_string(at_3->GetSeqNum()) + ") " + std::to_string(ideal_dist_i_3));
               // std::cout << "INFO:: Alpha Helix Bond restraint ("
               //          << at_1->name << " " << at_1->GetSeqNum() << ") to ("
               //          << at_2->name << " " << at_2->GetSeqNum() << ") " << ideal_dist_i_4 << std::endl;
               logger.log(log_t::INFO, "Alpha Helix Bond restraint (" + std::string(at_1->name) + " " + std::to_string(at_1->GetSeqNum()) +
                          ") to (" + std::string(at_2->name) + " " + std::to_string(at_2->GetSeqNum()) + ") " + std::to_string(ideal_dist_i_4));
            }
            n_helical_restraints += 2;
         } else {
	    if (at_1 && at_3) {
	       at_1->GetUDData(udd_atom_index_handle, index_1);
	       at_3->GetUDData(udd_atom_index_handle, index_3);
	       std::vector<bool> fixed_flags_2 = make_fixed_flags(index_1, index_3);
	       double ideal_dist_i_3 = 3.18;
	       add(BOND_RESTRAINT, index_1, index_3, fixed_flags_2, ideal_dist_i_3, pseudo_bond_esd, 1.2);

               if (verbose_geometry_reporting != QUIET) {
                  // std::cout << "INFO:: Alpha Helix Bond restraint ("
                  //          << at_1->name << " " << at_1->GetSeqNum() << ") to ("
                  //          << at_3->name << " " << at_3->GetSeqNum() << ") " << ideal_dist_i_3 << std::endl;
                  logger.log(log_t::INFO, "Alpha Helix Bond restraint (" + std::string(at_1->name) + " " + std::to_string(at_1->GetSeqNum()) +
                             ") to (" + std::string(at_3->name) + " " + std::to_string(at_3->GetSeqNum()) + ") " + std::to_string(ideal_dist_i_3));
               }
               n_helical_restraints += 1;
	    }
	 }
      }

   }

   if (verbose_geometry_reporting) {
      // std::cout << "INFO:: added " << n_helical_restraints << " helical restraints" << std::endl;
      logger.log(log_t::INFO, "added", n_helical_restraints, "helical restraints");
   }

   if (console_output_for_restraints_generation_timings) {
      auto tp_1 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      // std::cout << "INFO:: Timing for auto-helix " << d10 << " microseconds" << std::endl;
      logger.log(log_t::INFO, "Timing for auto-helix " + std::to_string(d10) + " microseconds");
   }

}


// restraint_addition_mode can be AUTO_HELIX - restrain anything that looks like a helix (alpha currently)
// or EVERYTHING_HELICAL - add helix restrains to residue with the same chain id and in a residue range
// that matches a H-bonded residue pair of a helix.
void
coot::restraints_container_t::make_helix_pseudo_bond_restraints_from_res_vec() {

   // Note to self: restraint_addition_mode is no longer used/checked in this function -
   // it can be removed from the call (was an experiment)

   // this doesn't do the right thing if there are insertion codes. Maybe I could check for that
   // here and jump out at the start if so.

   float pseudo_bond_esd = 0.02; // 0.035 seems reasonable.
                                 // 0.035 20191022-PE seems a bit weak now.

   // this double loop might be hideous for many hundreds of residues
   //
   for (std::size_t ir=0; ir<residues_vec.size(); ir++) {
      for (std::size_t jr=0; jr<residues_vec.size(); jr++) {
	 if (residues_vec[ir].second->GetChain() == residues_vec[jr].second->GetChain()) {
	    // check that at least one of them is not fixed
	    if (residues_vec[ir].first == false || residues_vec[jr].first == false) {
	       bool jr_is_upstream = false;
	       bool jr_is_downstream = false; // further along the chain (higher chain id)
	       int res_no_delta = residues_vec[jr].second->GetSeqNum() - residues_vec[ir].second->GetSeqNum();
	       if (res_no_delta == 3)
		  jr_is_downstream = true;
	       if (res_no_delta == 4)
		  jr_is_downstream = true;
	       if (res_no_delta == -3)
		  jr_is_upstream = true;
	       if (res_no_delta == -4)
		  jr_is_upstream = true;

	       // actually, we only need to check downstream because
	       // the double loop with catch the reverse direction

	       if (jr_is_downstream) {

		  mmdb::Atom **residue_atoms_1 = 0;
		  int n_residue_atoms_1;
		  residues_vec[ir].second->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
		  for (int iat=0; iat<n_residue_atoms_1; iat++) {
		     mmdb::Atom *at_1 = residue_atoms_1[iat];
		     std::string atom_name_1 = at_1->GetAtomName();
		     if (atom_name_1 == " O  ") {
			mmdb::Atom **residue_atoms_2 = 0;
			int n_residue_atoms_2;
			residues_vec[jr].second->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
			for (int jat=0; jat<n_residue_atoms_2; jat++) {
			   mmdb::Atom *at_2 = residue_atoms_2[jat];
			   std::string atom_name_2 = at_2->GetAtomName();
			   if (atom_name_2 == " N  ") {
			      std::string alt_conf_1 = at_1->altLoc;
			      std::string alt_conf_2 = at_2->altLoc;
			      if (alt_conf_1 == alt_conf_2) {

				 int index_1 = -1;
				 int index_2 = -1;
				 at_1->GetUDData(udd_atom_index_handle, index_1);
				 at_2->GetUDData(udd_atom_index_handle, index_2);
				 std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);
				 double ideal_dist = 2.919;
				 if (res_no_delta == 3)
				    ideal_dist = 3.181;
				 add(BOND_RESTRAINT, index_1, index_2, fixed_flags, ideal_dist, pseudo_bond_esd, 1.2);
				 std::cout << "Helix Bond restraint ("
					   << at_1->name << " " << at_1->GetSeqNum() << ") to ("
					   << at_2->name << " " << at_2->GetSeqNum() << ") " << ideal_dist << std::endl;
			      }
			   }
			   if (atom_name_2 == " O  ") {
			      std::string alt_conf_1 = at_1->altLoc;
			      std::string alt_conf_2 = at_2->altLoc;
			      if (alt_conf_1 == alt_conf_2) {

				 int index_1 = -1;
				 int index_2 = -1;
				 at_1->GetUDData(udd_atom_index_handle, index_1);
				 at_2->GetUDData(udd_atom_index_handle, index_2);
				 std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);
				 double ideal_dist = 6.16;
				 if (res_no_delta == 3)
				    ideal_dist = 4.92;
                                 double O_O_pseudo_bond_esd = 0.07; // guess
				 add(BOND_RESTRAINT, index_1, index_2, fixed_flags, ideal_dist, O_O_pseudo_bond_esd, 1.2);
				 std::cout << "Helix Bond restraint ("
					   << at_1->name << " " << at_1->GetSeqNum() << ") to ("
					   << at_2->name << " " << at_2->GetSeqNum() << ") " << ideal_dist << std::endl;
			      }
                           }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }


}


void
coot::restraints_container_t::make_strand_pseudo_bond_restraints() { 

   // This method of making pseudo bonds relies on the residue range
   // being continuous in sequence number (seqNum) and no insertion
   // codes messing up the number scheme.  If these are not the case
   // then we will make bonds between the wrongs atoms (of the wrong
   // residues)... a more sophisticated algorithm would be needed.
   // 
   // The method ignores residues with alt confs.
   //
   float pseudo_bond_esd = 0.08; // just a guess
   int selHnd = mol->NewSelection();
   int nSelResidues;
   mmdb::PPResidue SelResidue;
   mmdb::PPAtom res_1_atoms = NULL;
   mmdb::PPAtom res_2_atoms = NULL;
   mmdb::PPAtom res_3_atoms = NULL;
   int n_res_1_atoms;
   int n_res_2_atoms;
   int n_res_3_atoms;
   int index1 = -1; 
   int index2 = -1; 
   int index3 = -1; 
   mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"",  // altLocs
		mmdb::SKEY_NEW // selection key
		);
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 0) {
      for (int i=1; i<nSelResidues; i++) {
         // nO -> (n-1)O 4.64  // symmetric.  No need to specify both forwards
	 // Angle nO-(n+1)O-(n+2)O: 
	 SelResidue[i]->GetAtomTable(res_1_atoms, n_res_1_atoms);
	 if (res_1_atoms) { 
	    for (int iat1=0; iat1<n_res_1_atoms; iat1++) {
	       std::string at_1_name(res_1_atoms[iat1]->name);
	       // O Pseudo bonds and angles
	       if (at_1_name == " O  ") {
		  mmdb::Residue *contact_res = SelResidue[i-1];
		  if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 1)) {
		     contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		     if (res_2_atoms) { 
			for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
			   std::string at_2_name(res_2_atoms[iat2]->name);
			   if (at_2_name == " O  ") {
			      std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			      res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			      res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			      add(BOND_RESTRAINT, index1, index2, fixed_flags,
				  4.64, pseudo_bond_esd, 1.2);
			      std::cout << "Strand Bond restraint ("
					<< res_1_atoms[iat1]->name << " "
					<< res_1_atoms[iat1]->GetSeqNum() << ") to ("
					<< res_2_atoms[iat2]->name << " "
					<< res_2_atoms[iat2]->GetSeqNum() << ") 4.64" << std::endl;

			      // now the pseudo angle
			      if (i<(nSelResidues-1)) { 
				 mmdb::Residue *contact_res_2 = SelResidue[i+1];
				 if (SelResidue[i]->GetSeqNum() == (contact_res_2->GetSeqNum() - 1)) {
				    contact_res_2->GetAtomTable(res_3_atoms, n_res_3_atoms);
				    for (int iat3=0; iat3<n_res_3_atoms; iat3++) {
				       std::string at_3_name(res_3_atoms[iat3]->name);
				       if (at_3_name == " O  ") {
					  std::vector<bool> fixed_flag =
					     make_fixed_flags(index2, index1, index3);
					  res_3_atoms[iat3]->GetUDData(udd_atom_index_handle, index3);
					  // 98.0 degrees
					  add(ANGLE_RESTRAINT, index2, index1, index3,
					      fixed_flag, 98.0, 0.5, false);
					  std::cout << "Strand Angle restraint ("
						    << res_1_atoms[iat1]->name << " "
						    << res_1_atoms[iat1]->GetSeqNum() << ") to ("
						    << res_2_atoms[iat2]->name << " "
						    << res_2_atoms[iat2]->GetSeqNum()
						    << ") to ("
						    << res_3_atoms[iat3]->name << " "
						    << res_3_atoms[iat3]->GetSeqNum()
						    << ") 98.0 " << std::endl;
					  break;
				       }
				    }
				 }
			      }
			      break;
			   }
			}
		     }
		  }
	       } // end of O

	       // Now make a CA-CA-CA pseudo angle of 120 degrees
	       if (at_1_name == " CA ") {
		  mmdb::Residue *contact_res_2 = SelResidue[i-1];
		  if (SelResidue[i]->GetSeqNum() == (contact_res_2->GetSeqNum() + 1)) {
		     contact_res_2->GetAtomTable(res_2_atoms, n_res_2_atoms);
		     if (res_2_atoms) { 
			for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
			   std::string at_2_name(res_2_atoms[iat2]->name);
			   if (at_2_name == " CA ") {
			      if (i<(nSelResidues-1)) {
				 mmdb::Residue *contact_res_3 = SelResidue[i+1];
				 if (SelResidue[i]->GetSeqNum() == (contact_res_3->GetSeqNum() - 1)) {
				    contact_res_3->GetAtomTable(res_3_atoms, n_res_3_atoms);
				    for (int iat3=0; iat3<n_res_3_atoms; iat3++) {
				       std::string at_3_name(res_3_atoms[iat3]->name);
				       if (at_3_name == " CA ") {
					  std::vector<bool> fixed_flag =
					     make_fixed_flags(index1, index2, index3);
					  res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
					  res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
					  res_3_atoms[iat3]->GetUDData(udd_atom_index_handle, index3);
					  add(ANGLE_RESTRAINT, index2, index1, index3,
					      fixed_flag, 120.0, 0.5, false);
					  std::cout << "Strand Angle restraint ("
						    << res_1_atoms[iat1]->name << " "
						    << res_1_atoms[iat1]->GetSeqNum() << ") to ("
						    << res_2_atoms[iat2]->name << " "
						    << res_2_atoms[iat2]->GetSeqNum()
						    << ") to ("
						    << res_3_atoms[iat3]->name << " "
						    << res_3_atoms[iat3]->GetSeqNum()
						    << ") 120.0 " << std::endl;
					  break;
				       }
				    }
				 }
			      }
			      break;
			   }
			}
		     }
		  }
	       }
	    }
	 } 
      }
   }
   mol->DeleteSelection(selHnd);
}

#include "coot-utils/coot-h-bonds.hh"

void
coot::restraints_container_t::make_h_bond_restraints_from_res_vec_auto(const coot::protein_geometry &geom, int imol) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   int SelHnd = mol->NewSelection(); // d
   h_bonds hbs;

   auto tp_1 = std::chrono::high_resolution_clock::now();
   for(unsigned int i=0; i<residues_vec.size(); i++) {
      mmdb::Residue *r = residues_vec[i].second;
      residue_spec_t(r).select_atoms(mol, SelHnd, mmdb::SKEY_OR);
   }
   auto tp_2 = std::chrono::high_resolution_clock::now();
   std::vector<h_bond> v = hbs.get(SelHnd, SelHnd, mol, geom, imol);
   auto tp_3 = std::chrono::high_resolution_clock::now();

   unsigned int n_bonds = 0;

   for (unsigned int i=0; i<v.size(); i++) {
      const h_bond &hb = v[i];
      if (hb.donor) {
         if (hb.acceptor) {
            clipper::Coord_orth b(co(hb.donor) - co(hb.acceptor));

            int index_1 = -1, index_2 = -1;
            int udd_get_data_status_1 = hb.donor->GetUDData(   udd_atom_index_handle, index_1);
            int udd_get_data_status_2 = hb.acceptor->GetUDData(udd_atom_index_handle, index_2);

            if (udd_get_data_status_1 == mmdb::UDDATA_Ok &&
                udd_get_data_status_2 == mmdb::UDDATA_Ok) {

               double bl = sqrt(b.lengthsq()); // bond length
               std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);

               bool as_bond = true;
               if (as_bond) {
                  add_h_bond(BOND_RESTRAINT, index_1, index_2, fixed_flags, bl, 0.1);
               } else {
                  add_geman_mcclure_distance(GEMAN_MCCLURE_DISTANCE_RESTRAINT, index_1, index_2, fixed_flags, bl, 0.1);
               }
               n_bonds++;
            }
         }
      }
   }

   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   std::cout << "------------------- timing: " << d10 << " " << d21 << " " << d32
             <<  " milliseconds to find " << v.size() << " H-bonds " << std::endl;
   std::cout << "DEBUG:: made " << n_bonds << " hydrogen bonds " << std::endl;

   mol->DeleteSelection(SelHnd);

}


int
coot::restraints_container_t::make_monomer_restraints(int imol,
						      const coot::protein_geometry &geom,
						      short int do_residue_internal_torsions) {

   // std::cout << "------------------------ in make_monomer_restraints() "
   // << do_residue_internal_torsions << std::endl;

   if (from_residue_vector)
      return make_monomer_restraints_from_res_vec(imol, geom, do_residue_internal_torsions);
   else
      return make_monomer_restraints_by_linear(imol, geom, do_residue_internal_torsions);

}

int
coot::restraints_container_t::make_monomer_restraints_by_linear(int imol,
								const coot::protein_geometry &geom,
								bool do_residue_internal_torsions) {

   // note: mini-rsr uses only the residue vector method
   
   int iret = 0;
   
   int selHnd = mol->NewSelection();
   int nSelResidues;
   restraint_counts_t sum;

   mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",  // altLocs
		mmdb::SKEY_NEW // selection key
		);
   SelResidue_active = NULL;
   mol->GetSelIndex (selHnd, SelResidue_active, nSelResidues);
//    std::cout << "INFO:: GetSelIndex returned " << nSelResidues
// 	     << " residues (monomer restraints) " << std::endl;
   // save the (new (7Nov2003)) class variables (used in non_bonded
   // stuff) that keep the "active" (as opposed to "flanking") residues:
   nSelResidues_active = nSelResidues;
   // std::cout << "------------------------ in make_monomer_restraints_by_linear() nSelResidues "
   // << nSelResidues << std::endl;
   if (nSelResidues > 0) { 
      for (int i=0; i<nSelResidues; i++) {
	 if (SelResidue_active[i]) {
	    // std::cout << "------- calling make_monomer_restraints_by_residue() " << std::endl;
	    restraint_counts_t local =
	       make_monomer_restraints_by_residue(imol, SelResidue_active[i], geom,
						  do_residue_internal_torsions);
	    sum += local;
	 }
      }
   } else {
      std::cout << "get_monomer_restraints: There were no residues selected!? "
		<< std::endl;
   }
   // mol->DeleteSelection(selHnd); // -> makes crash! 6-feb-2004.
   //
   // This is because SelResidue_active is used elsewhere.
   // 
   // This should go into the destructor, I guess.

   sum.report(do_residue_internal_torsions);
   if (verbose_geometry_reporting != QUIET) {
      // std::cout << "INFO:: by_linear() created " << size() << " restraints" << std::endl;
      logger.log(log_t::INFO, "by_linear() created", size(), "restraints");
      std::cout << std::endl;
   }
   return iret; // return 1 on success.  Hmm... how is this set? (and subsequently used?)
}

int
coot::restraints_container_t::make_monomer_restraints_from_res_vec(int imol,
								   const coot::protein_geometry &geom,
								   bool do_residue_internal_torsions) {

   int iret = 0;

   restraint_counts_t sum;

   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
      if (residues_vec[ir].second) {
         restraint_counts_t local =
	    make_monomer_restraints_by_residue(imol, residues_vec[ir].second, geom,
					       do_residue_internal_torsions);
         sum += local;
      } else {
         std::cout << "ERROR:: in make_monomer_restraints_from_res_vec() null residue "
                   << ir << " of " << residues_vec.size() << std::endl;
      }
   }

   if (verbose_geometry_reporting != QUIET) {
      // std::cout << "INFO:: make_monomer_restraints_from_res_vec() of size " << residues_vec.size()
      // << " created " << size() << " monomer restraints " << std::endl;
      logger.log(log_t::INFO, logging::function_name_t("make_monomer_restraints_from_res_vec"),
                 {std::string("of size"), residues_vec.size(), std::string("created"),
                  size(), "monomer restraints"});
      sum.report(do_residue_internal_torsions);
   }
   return iret;
}



coot::restraints_container_t::restraint_counts_t
coot::restraints_container_t::make_monomer_restraints_by_residue(int imol, mmdb::Residue *residue_p,
								 const protein_geometry &geom,
								 bool do_residue_internal_torsions) {

   if (false)
      std::cout << "----------- make_monomer_restraints_by_residue() "
                << residue_p->GetChain() << " "
                << residue_p->GetSeqNum()<< " : "
                << residue_p->GetResName() << std::endl;

   restraint_counts_t local;

   if (! residue_p) {
      std::cout << "ERROR in make_monomer_restraints_by_residue() null residue" << std::endl;
      return local;
   }

   std::string pdb_resname(residue_p->name);
   if (pdb_resname == "UNK") pdb_resname = "ALA";

   if (false)
      std::cout << "--------------- make_monomer_restraints_by_residue() called "
                << residue_spec_t(residue_p) << " with " << residue_p->GetNumberOfAtoms() << " atoms"
                <<  " and using type :" << pdb_resname << ": and imol "
                << imol << " do_residue_internal_torsions: "
                << do_residue_internal_torsions << std::endl;

   // idr: index dictionary residue
   int idr = geom.get_monomer_restraints_index(pdb_resname, imol, false);

   if (false) {
      std::cout << "debug:: in make_monomer_restraints_by_residue() here with idr " << idr << std::endl;
   }

   if (idr == -1) {

      std::cout << "ERROR:: failed to get restraints index for monomer " << pdb_resname << std::endl;

   } else {

      // if (geom[idr].comp_id == pdb_resname) {
      // old style comp_id usage
      // if (dictionary_name_matches_coords_resname(geom[idr].comp_id,pdb_resname)) {

      // OK, we need the 3 letter code for carbohydrates, the
      // comp_id for nucleotides:
      //
      // comp_id 3-letter-code name group
      // Ar         A        'Adenosine                    ' RNA                33  22 .
      // GAL-b-D    GAL      'beta_D_galactose             ' D-pyranose         24  12 .


      // now get a list of atoms in that residue
      // (SelResidue[i]) and compare them to the atoms in
      // geom[idr].bond_restraint[ib].

      int i_no_res_atoms = 0;
      mmdb::PPAtom res_selection = NULL;
      residue_p->GetAtomTable(res_selection, i_no_res_atoms);

      if (i_no_res_atoms > 0) {

	 if (util::is_standard_amino_acid_name(pdb_resname))
	    local += add_N_terminal_residue_bonds_and_angles_to_hydrogens(residue_p);

	 if (restraints_usage_flag & BONDS_MASK)
	    local.n_bond_restraints += add_bonds(idr, res_selection, i_no_res_atoms,
						 residue_p, geom);

	 if (restraints_usage_flag & ANGLES_MASK)
	    local.n_angle_restraints += add_angles(idr, res_selection, i_no_res_atoms,
						   residue_p, geom);

	 if (restraints_usage_flag & TORSIONS_MASK) {
	    if (do_residue_internal_torsions) {
	       std::string residue_type = residue_p->GetResName();
	       if (residue_type != "PRO")
		  local.n_torsion_restr += add_torsions(idr, res_selection, i_no_res_atoms,
							residue_p, geom, torsion_restraints_weight);
	    }
	 }

	 if (restraints_usage_flag & PLANES_MASK)
	    local.n_plane_restraints += add_planes(idr, res_selection, i_no_res_atoms,
						   residue_p, geom);

         if (restraints_usage_flag & IMPROPER_DIHEDRALS_MASK) {
            // dictionaries need to be converted for this to happen.
	    int n = add_planes_as_improper_dihedrals(idr, res_selection, i_no_res_atoms, residue_p, geom);
            if (false)
               std::cout << "debug:: in make_monomer_restraints_by_residue() made "
                         << n << " improper_dihedrals" << std::endl;
            local.n_improper_dihedral_restr += n;
	 }


	 if (restraints_usage_flag & CHIRAL_VOLUME_MASK) {
	    local.n_chiral_restr += add_chirals(idr, res_selection, i_no_res_atoms, residue_p, geom);
	 }

	 restraint_counts_t mod_counts = apply_mods(idr, res_selection, i_no_res_atoms, residue_p, geom);
	 // now combine mod_counts with local
      }
   }

   // local.report(false);
   return local;
}


int
coot::restraints_container_t::add_bonds(int idr, mmdb::PPAtom res_selection,
					int i_no_res_atoms,
					mmdb::PResidue SelRes,
					const coot::protein_geometry &geom) {

   int n_bond_restr = 0;
   int index1, index2;
   bool debug = false;

   if (debug)
      std::cout << "in add_bonds() for " << residue_spec_t(SelRes) << std::endl;

   const dictionary_residue_restraints_t &dict = geom[idr].second;

   if (debug) {
      std::cout << "debug:: dictionary index idr " << idr << std::endl;
      std::cout << "debug:: idr indexes dictionary with name " << dict.residue_info.comp_id << " "
                << dict.residue_info.three_letter_code << " " << dict.residue_info.name << std::endl;
   }

   for (unsigned int ib=0; ib<dict.bond_restraint.size(); ib++) {
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name1(res_selection[iat]->name);

	 if (debug)
	    std::cout << "comparing first (pdb) :" << pdb_atom_name1
		      << ": with (dict) :"
		      << dict.bond_restraint[ib].atom_id_1_4c()
		      << ":" << std::endl;

	 if (pdb_atom_name1 == dict.bond_restraint[ib].atom_id_1_4c()) {
	    for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

	       std::string pdb_atom_name2(res_selection[iat2]->name);

	       if (debug)
		  std::cout << "comparing second (pdb) :" << pdb_atom_name2
			    << ": with (dict) :"
			    << dict.bond_restraint[ib].atom_id_2_4c()
			    << ":" << std::endl;

	       if (pdb_atom_name2 == dict.bond_restraint[ib].atom_id_2_4c()) {

		  // check that the alt confs aren't different
		  std::string alt_1(res_selection[iat ]->altLoc);
		  std::string alt_2(res_selection[iat2]->altLoc);
		  if (alt_1 == "" || alt_2 == "" || alt_1 == alt_2) {

		     if (debug) {
			std::cout << "atom match 1 " << pdb_atom_name1;
			std::cout << " atom match 2 " << pdb_atom_name2
				  << std::endl;
		     }

		     // now we need the indices of
		     // pdb_atom_name1 and
		     // pdb_atom_name2 in asc.atom_selection:

		     //  		  int index1_old = get_asc_index(pdb_atom_name1,
		     //  					     SelRes->seqNum,
		     //  					     SelRes->GetChainID());
		     //  		  int index2_old = get_asc_index(pdb_atom_name2,
		     //  					     SelRes->seqNum,
		     //  					     SelRes->GetChainID());

		     int udd_get_data_status_1 = res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
		     int udd_get_data_status_2 = res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);

		     // set the UDD flag for this residue being bonded/angle with 
		     // the other

		     if (udd_get_data_status_1 == mmdb::UDDATA_Ok &&
			 udd_get_data_status_2 == mmdb::UDDATA_Ok) { 
		  
			bonded_atom_indices[index1].insert(index2);
			bonded_atom_indices[index2].insert(index1);

			// this needs to be fixed for fixed atom (rather
			// than just knowing that these are not flanking
			// atoms).
			// 
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);

			if (debug) {
			   std::string altconf_1("\"");
			   altconf_1 += atom[index1]->altLoc;
			   altconf_1 += "\"";
			   std::string altconf_2("\"");
			   altconf_2 += atom[index2]->altLoc;
			   altconf_2 += "\"";
			   std::cout << "creating (monomer) bond restraint, idr " << idr
				     << " with fixed flags "
				     << fixed_flags[0] << " " << fixed_flags[1] << " "
				     << atom[index1]->GetSeqNum() << " "
				     << "\"" << atom[index1]->name << "\" "
				     << std::setw(3) << altconf_1 << " to "
				     << atom[index2]->GetSeqNum() << " "
				     << "\"" << atom[index2]->name << "\" "
				     << std::setw(3) << altconf_2 << " "
				     << "bond restraint index " << n_bond_restr << "\n";
			}
			try {
			   add(BOND_RESTRAINT, index1, index2,
			       fixed_flags,
			       dict.bond_restraint[ib].value_dist(),
			       dict.bond_restraint[ib].value_esd(),
			       1.2);  // junk value

			   n_bond_restr++;

			   // now cache the parent energy type: for looking up the type
			   // of the H atom so that we can adjust the bond type (target
			   // distance) when making non-bonded contacts
			   //
			   if (is_hydrogen(atom[index1])) {
			      mmdb::Atom *H_at = atom[index1];
			      mmdb::Atom *parent_at = atom[index2];
			      std::string atom_name(parent_at->name);
			      std::string te = dict.type_energy(atom_name);
			      hb_t hbt = geom.get_h_bond_type(te);
			      H_atom_parent_energy_type_atom_map[H_at] = hbt;
                              // std::cout << "settings H_atom_parent route-A H_atom " << atom_spec_t(H_at) << " "
                              //           << atom_spec_t(parent_at) << " type " << hbt << std::endl;
			   }
			   if (is_hydrogen(atom[index2])) {
			      mmdb::Atom *H_at = atom[index2];
			      mmdb::Atom *parent_at = atom[index1];
			      std::string atom_name(parent_at->name);
			      std::string te = dict.type_energy(atom_name);
			      hb_t hbt = geom.get_h_bond_type(te);
			      H_atom_parent_energy_type_atom_map[H_at] = hbt;
                              // std::cout << "settings H_atom_parent route-B H_atom " << atom_spec_t(H_at) << " "
                              //           << atom_spec_t(parent_at) << " type " << hbt << std::endl;

                              // It seems to me that NR5 in the dictionary for HIS and TRP is wrong. It should
                              // be a donor.
                              // So kludge that in here
                              if (dict.residue_info.comp_id == "HIS") {
                                 mmdb::Atom *parent_at = atom[index1];
                                 std::string parent_atom_name(parent_at->name);
                                 if (parent_atom_name == "ND1 ") // PDBv3 fixme
                                    H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                                 if (parent_atom_name == "NE2 ") // PDBv3 fixme
                                    H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                              }
                              if (dict.residue_info.comp_id == "TRP") {
                                 mmdb::Atom *parent_at = atom[index1];
                                 std::string parent_atom_name(parent_at->name);
                                 if (parent_atom_name == "NE1 ") // PDBv3 fixme
                                    H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                              }
                           }
			}

			catch (const std::runtime_error &rte) {

			   // do nothing, it's not really an error if the dictionary
			   // doesn't have target geometry (the bonding description came
			   // from a Chemical Component Dictionary entry for example).
			   std::cout << "trapped a runtime_error on adding bond restraint "
				     << " no target. " << rte.what() << std::endl;
			}
		     } else {
			std::cout << "ERROR:: Caught Enrico Stura bug.  How did it happen?" << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_bond_restr;
}

int
coot::restraints_container_t::add_angles(int idr, mmdb::PPAtom res_selection,
					 int i_no_res_atoms,
					 mmdb::PResidue SelRes,
					 const coot::protein_geometry &geom) {

   int n_angle_restr = 0;
   int index1, index2, index3;

   std::vector<std::string> string_atom_names(i_no_res_atoms);
   for (int iat=0; iat<i_no_res_atoms; iat++)
      string_atom_names[iat] = res_selection[iat]->name;

//    std::cout << "There are " << geom[idr].angle_restraint.size()
// 	     << " angle restraints for this residue type" << std::endl; 

   for (unsigned int ib=0; ib<geom[idr].second.angle_restraint.size(); ib++) {
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 const std::string &pdb_atom_name1 = string_atom_names[iat];

//  	 std::cout << "angle:  comparing :" << pdb_atom_name1 << ": with :"
//  		   << geom[idr].angle_restraint[ib].atom_id_1_4c()
//  		   << ":" << std::endl;
	 
	 if (pdb_atom_name1 == geom[idr].second.angle_restraint[ib].atom_id_1_4c()) {
	    for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

	       const std::string &pdb_atom_name2 = string_atom_names[iat2];
	       if (pdb_atom_name2 == geom[idr].second.angle_restraint[ib].atom_id_2_4c()) {
				    
// 		  std::cout << "angle: atom match 1 " << pdb_atom_name1;
// 		  std::cout << " atom match 2 " << pdb_atom_name2
// 			    << std::endl;

		  for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		     
                     const std::string &pdb_atom_name3 = string_atom_names[iat3];
		     if (pdb_atom_name3 == geom[idr].second.angle_restraint[ib].atom_id_3_4c()) {

			std::string alt_1(res_selection[iat ]->altLoc);
			std::string alt_2(res_selection[iat2]->altLoc);
			std::string alt_3(res_selection[iat3]->altLoc);

			if (((alt_1 == alt_2) && (alt_1 == alt_3)) ||
			    ((alt_1 == ""   ) && (alt_2 == alt_3)) ||
			    ((alt_2 == ""   ) && (alt_1 == alt_3)) ||
			    ((alt_3 == ""   ) && (alt_1 == alt_2)))
			   {
			
			   res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
			   res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
			   res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);

			   // std::cout << "add_angles: " << index1_old << " " << index1 << std::endl;
			   // std::cout << "add_angles: " << index2_old << " " << index2 << std::endl;
			   // std::cout << "add_angles: " << index3_old << " " << index3 << std::endl;
			
			   // set the UDD flag for this residue being bonded/angle with 
			   // the other
			
			   bonded_atom_indices[index1].insert(index3);
			   bonded_atom_indices[index3].insert(index1);
		  
			   // this needs to be fixed for fixed atom (rather
			   // than just knowing that these are not flanking
			   // atoms).
			   // 
			   std::vector<bool> fixed_flag = make_fixed_flags(index1, index2, index3);
			   bool is_single_H_atom_angle_restraint = false;
			   unsigned int nH = 0;
			   if (is_hydrogen(res_selection[iat ])) nH++;
			   if (is_hydrogen(res_selection[iat3])) nH++;
			   if (nH == 1) is_single_H_atom_angle_restraint = true;

			   add(ANGLE_RESTRAINT, index1, index2, index3, fixed_flag,
			       geom[idr].second.angle_restraint[ib].angle(),
			       geom[idr].second.angle_restraint[ib].esd(),
			       is_single_H_atom_angle_restraint);
			   n_angle_restr++;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_angle_restr;
}


bool
coot::restraints_container_t::add_torsion_internal(const coot::dict_torsion_restraint_t &torsion_restraint,
                                                   mmdb::PPAtom res_selection, int i_no_res_atoms,
                                                   const double &torsion_restraint_weight) {

   bool status = false;

   { // cut and paste

      // Joel Bard fix: Don't add torsion restraints for torsion that
      // have either s.d. or period 0

      if (torsion_restraint.periodicity() > 0) { // we had this test most inner
	 if (torsion_restraint.esd() > 0.000001) { // new test

	    // now find the atoms
	    for (int iat=0; iat<i_no_res_atoms; iat++) {
	       std::string pdb_atom_name1(res_selection[iat]->name);

	       if (pdb_atom_name1 == torsion_restraint.atom_id_1_4c()) {
		  for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

		     std::string pdb_atom_name2(res_selection[iat2]->name);
		     if (pdb_atom_name2 == torsion_restraint.atom_id_2_4c()) {

			for (int iat3=0; iat3<i_no_res_atoms; iat3++) {

			   std::string pdb_atom_name3(res_selection[iat3]->name);
			   if (pdb_atom_name3 == torsion_restraint.atom_id_3_4c()) {

			      for (int iat4=0; iat4<i_no_res_atoms; iat4++) {

				 std::string pdb_atom_name4(res_selection[iat4]->name);
				 if (pdb_atom_name4 == torsion_restraint.atom_id_4_4c()) {

				    // now we need the indices of
				    // pdb_atom_name1 and
				    // pdb_atom_name2 in asc.atom_selection:

                                    // kill off weird dictionary torsions here
                                    if (pdb_atom_name1 == " O  ")
                                       if (pdb_atom_name2 == " C  ")
                                          if (pdb_atom_name3 == " CA ")
                                             continue;
                                    if (pdb_atom_name1 == " CB ")
                                       if (pdb_atom_name2 == " CA  ")
                                          if (pdb_atom_name3 == " N ")
                                             if (pdb_atom_name4 == " H ")
                                                continue;

				    int index1;
				    int index2;
				    int index3;
				    int index4;

				    res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
				    res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
				    res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);
				    res_selection[iat4]->GetUDData(udd_atom_index_handle, index4);

				    double torsion_angle = torsion_restraint.angle();
				    if (torsion_angle < 0)
				       torsion_angle += 360;
				    if (torsion_angle > 360)
				       torsion_angle -= 360;

                                    std::string alt_conf_1(res_selection[iat]->altLoc);
                                    std::string alt_conf_2(res_selection[iat2]->altLoc);
                                    std::string alt_conf_3(res_selection[iat3]->altLoc);
                                    std::string alt_conf_4(res_selection[iat4]->altLoc);

                                    bool alt_confs_match = false;
                                    if (alt_conf_1 == "" || alt_conf_1 == alt_conf_2)
                                       if (alt_conf_2 == "" || alt_conf_2 == alt_conf_3)
                                          if (alt_conf_3 == "" || alt_conf_3 == alt_conf_4 || alt_conf_4 == "")
                                             alt_confs_match = true;

                                    if (alt_confs_match) {
				       std::vector<bool> fixed_flags = make_fixed_flags(index1, index2, index3, index4);
				       add(TORSION_RESTRAINT, index1, index2, index3, index4,
					   fixed_flags,
					   torsion_angle,
					   torsion_restraint.esd(),
					   torsion_restraint_weight,
					   torsion_restraint.periodicity());

				       if (false) {
                                          get_print_lock();
				          std::cout << "debug:: Adding monomer torsion restraint: "
						    << index1 << " "
						    << index2 << " "
						    << index3 << " "
						    << index4 << " torsion "
						    << torsion_restraint.angle() << " esd "
						    << torsion_restraint.esd() << " period "
						    << torsion_restraint.periodicity() << " "
                                                    << " weight " << torsion_restraint_weight
						    << std::endl;
                                          release_print_lock();
                                       }
				       status = true;
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return status;
}

int
coot::restraints_container_t::add_torsions(int idr, mmdb::PPAtom res_selection,
					   int i_no_res_atoms,
					   mmdb::PResidue SelRes,
					   const coot::protein_geometry &geom,
                                           const double &torsion_restraints_weight) {

   int n_torsion_restr = 0;
   const std::vector<dict_torsion_restraint_t> &torsion_restraints = geom[idr].second.torsion_restraint;

   for (unsigned int ib=0; ib<torsion_restraints.size(); ib++) {
      const dict_torsion_restraint_t &torsion_restraint = torsion_restraints[ib];
      bool status = add_torsion_internal(torsion_restraint, res_selection, i_no_res_atoms, torsion_restraints_weight);
      if (status)
         n_torsion_restr++;
   }

   return n_torsion_restr;
}


int
coot::restraints_container_t::add_chirals(int idr, mmdb::PPAtom res_selection,
					  int i_no_res_atoms,
					  mmdb::PResidue SelRes,
					  const coot::protein_geometry &geom) { 

   int n_chiral_restr = 0;
   int index1, index2, index3, indexc;
   
   //   std::cout << "DEBUG:: trying to add chirals for this residue..." << std::endl;
   
   std::vector<std::string> string_atom_names(i_no_res_atoms);
   for (int iat=0; iat<i_no_res_atoms; iat++)
      string_atom_names[iat] = res_selection[iat]->name;

   for (unsigned int ic=0; ic<geom[idr].second.chiral_restraint.size(); ic++) {
      // for now, let's just reject restraints that are a "both",
      // better would be to check the geometry and refine to the one
      // that is closest.

      const dict_chiral_restraint_t &dcr = geom[idr].second.chiral_restraint[ic];

      if (!geom[idr].second.chiral_restraint[ic].is_a_both_restraint()) { 
	 for (int iat1=0; iat1<i_no_res_atoms; iat1++) {
	    const std::string &pdb_atom_name1 = string_atom_names[iat1];
	    if (pdb_atom_name1 == geom[idr].second.chiral_restraint[ic].atom_id_1_4c()) {
	       
	       for (int iat2=0; iat2<i_no_res_atoms; iat2++) {
                  const std::string &pdb_atom_name2 = string_atom_names[iat2];
		  if (pdb_atom_name2 == geom[idr].second.chiral_restraint[ic].atom_id_2_4c()) {
		     
		     for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
                        const std::string &pdb_atom_name3 = string_atom_names[iat3];
			if (pdb_atom_name3 == geom[idr].second.chiral_restraint[ic].atom_id_3_4c()) {
			   
			   for (int iatc=0; iatc<i_no_res_atoms; iatc++) {
                              const std::string &pdb_atom_namec = string_atom_names[iatc];
			      if (pdb_atom_namec == geom[idr].second.chiral_restraint[ic].atom_id_c_4c()) {
				 
//   			      std::cout << "DEBUG:: adding chiral number " << ic << " for " 
//   					<< res_selection[iatc]->GetSeqNum() << " "
//   					<< res_selection[0]->GetResName()
//   					<< pdb_atom_namec << " bonds to "
//   					<< pdb_atom_name1 << " "
//   					<< pdb_atom_name2 << " "
//   					<< pdb_atom_name3 << " "
//   					<< std::endl;

				 std::string alt_conf_c = res_selection[iatc]->altLoc;
				 std::string alt_conf_1 = res_selection[iat1]->altLoc;
				 std::string alt_conf_2 = res_selection[iat2]->altLoc;
				 std::string alt_conf_3 = res_selection[iat3]->altLoc;

				 if (((alt_conf_1 == alt_conf_c) || (alt_conf_1 == "")) &&
				     ((alt_conf_2 == alt_conf_c) || (alt_conf_2 == "")) && 
				     ((alt_conf_3 == alt_conf_c) || (alt_conf_3 == ""))) { 

				    res_selection[iat1]->GetUDData(udd_atom_index_handle, index1);
				    res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
				    res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);
				    res_selection[iatc]->GetUDData(udd_atom_index_handle, indexc);


				    // does this chiral centre have exactly one hydrogen?
				    // If so, set chiral_hydrogen_index to the atom index.
				    // If not set to -1.
				    //

				    // int chiral_hydrogen_index_old = get_chiral_hydrogen_index(indexc, index1, index2, index3);

				    int chiral_hydrogen_index = get_chiral_hydrogen_index(indexc, index1, index2, index3, dcr);

				    if (fabs(geom[idr].second.chiral_restraint[ic].target_volume()) < 1000.0 &&
					fabs(geom[idr].second.chiral_restraint[ic].target_volume()) > 0.00001) {

				       if (false) // debug
					  std::cout << "   Adding chiral restraint for "
						    << res_selection[iatc]->name
						    << " " << res_selection[iatc]->GetSeqNum() <<  " "
						    << res_selection[iatc]->GetChainID()
						    << " with target volume "
						    << geom[idr].second.chiral_restraint[ic].target_volume()
						    << " with volume sigma "
						    << geom[idr].second.chiral_restraint[ic].volume_sigma()
						    << " with volume sign "
						    << geom[idr].second.chiral_restraint[ic].volume_sign
						    << " idr index: " << idr << " ic index: " << ic
						    << " chiral_hydrogen_index: " << chiral_hydrogen_index
						    << std::endl;

				       std::vector<bool> fixed_flags =
					  make_fixed_flags(indexc, index1, index2, index3);
				       simple_restraint sr(CHIRAL_VOLUME_RESTRAINT, indexc,
							   index1, index2, index3,
							   geom[idr].second.chiral_restraint[ic].volume_sign,
							   geom[idr].second.chiral_restraint[ic].target_volume(),
							   geom[idr].second.chiral_restraint[ic].volume_sigma(),
							   fixed_flags, chiral_hydrogen_index);
				       restraints_vec.push_back(sr); // push_back_restraint()
				       n_chiral_restr++;
				    } else {
				       std::cout << "WARNING:: Reject chiral restraint for "
						 << res_selection[iatc]->name
						 << " " << res_selection[iatc]->GetSeqNum() <<  " "
						 << res_selection[iatc]->GetChainID()
						 << " with target volume "
						 << geom[idr].second.chiral_restraint[ic].target_volume()
						 << " with volume sigma "
						 << geom[idr].second.chiral_restraint[ic].volume_sigma()
						 << " with volume sign "
						 << geom[idr].second.chiral_restraint[ic].volume_sign
						 << " idr index: " << idr << " ic index: " << ic
						 << " chiral_hydrogen_index: " << chiral_hydrogen_index
						 << std::endl;
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_chiral_restr;
}

int
coot::restraints_container_t::get_chiral_hydrogen_index(int indexc, int index_1, int index_2, int index_3,
							const coot::dict_chiral_restraint_t &dcr) const {

   int idx_chiral = -1;
   int n_hydrogen = 0;
   if (is_hydrogen(atom[index_1])) { n_hydrogen++; idx_chiral = index_1; }
   if (is_hydrogen(atom[index_2])) { n_hydrogen++; idx_chiral = index_2; }
   if (is_hydrogen(atom[index_3])) { n_hydrogen++; idx_chiral = index_3; }

   if (n_hydrogen == 1) {
      return idx_chiral;
   } else {
      return -1;
   }
}


// ------------------ old, slow -----------------------------
// is there a single hydrogen connected to this chiral centre?
// If so, return the index, if not return -1
// 
int
coot::restraints_container_t::get_chiral_hydrogen_index(int indexc, int index1, int index2, int index3) const {

   int r = -1;
   int n_H = 0;
   int H_atom_index = -1;

   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & coot::BONDS_MASK) {
	 if ( (*this)[i].restraint_type == coot::BOND_RESTRAINT) {
	    mmdb::Atom *at_1 = atom[(*this)[i].atom_index_1]; 
	    mmdb::Atom *at_2 = atom[(*this)[i].atom_index_2];
	    if ((*this)[i].atom_index_1 == indexc) {
	       if (is_hydrogen(at_2)) {
		  H_atom_index = (*this)[i].atom_index_2;
		  n_H++;
	       } 
	    }
	    if ((*this)[i].atom_index_2 == indexc) {
	       if (is_hydrogen(at_1)) {
		  H_atom_index = (*this)[i].atom_index_1;
		  n_H++;
	       }
	    }
	 }
      }
   }
   if (n_H == 1)
      return H_atom_index;
   else
      return -1;
} 

// Creates any number of simple_restraints for this monomer and adds
// them to restraints_vec.
// 
// idr provides the index of the comp_id (e.g. "ALA") match in geom.
// 
int
coot::restraints_container_t::add_planes(int idr, mmdb::PPAtom res_selection,
					 int i_no_res_atoms,
					 mmdb::PResidue SelRes,
					 const coot::protein_geometry &geom) {
   if (false)
      std::cout << "debug:: in add_planes(): with convert_plane_restraints_to_improper_dihedral_restraints_flag "
	        << convert_plane_restraints_to_improper_dihedral_restraints_flag << std::endl;
   if (! convert_plane_restraints_to_improper_dihedral_restraints_flag) {
      int n_added = add_planes_multiatom_eigen(idr, res_selection, i_no_res_atoms, SelRes, geom);
      // std::cout << "debug:: n_added (multiatom-eigen) " << n_added << std::endl;
      return n_added;
   } else {
      int n_added = add_planes_as_improper_dihedrals(idr, res_selection, i_no_res_atoms, SelRes, geom);
      // std::cout << "debug:: n_added (improper_dihedrals) " << n_added << std::endl;
      return n_added;
   }
}

int
coot::restraints_container_t::add_planes_multiatom_eigen(int idr, mmdb::PPAtom res_selection,
							 int i_no_res_atoms,
							 mmdb::PResidue SelRes,
							 const coot::protein_geometry &geom) {

   bool debug = false;

   if (debug)
      std::cout << "There are " << geom[idr].second.plane_restraint.size()
		<< " dictionary plane restraints for " << SelRes->seqNum << " type: "
		<< geom[idr].second.residue_info.comp_id << std::endl;

   int n_plane_restr = 0;
   // either altconfs are all the same,
   // or they are different, in which case, add the atoms with blank altconfs to
   // (only) each of the non-blank ones
   //
   std::vector<std::string> altconfs = util::get_residue_alt_confs(SelRes);
   bool all_altconfs_the_same = true;
   if (altconfs.size() > 1)
      all_altconfs_the_same = false;

   for (unsigned int ip=0; ip<geom[idr].second.plane_restraint.size(); ip++) {
      std::map<std::string, std::vector <std::pair<int, double> > > idx_and_sigmas;
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name(res_selection[iat]->name);
	 std::string alt_conf(res_selection[iat]->altLoc);
	 for (int irest_at=0; irest_at<geom[idr].second.plane_restraint[ip].n_atoms(); irest_at++) {
	    if (pdb_atom_name == geom[idr].second.plane_restraint[ip].atom_id(irest_at)) {
	       // is this slow?
// 	       int idx = get_asc_index(res_selection[iat]->name,
// 				       res_selection[iat]->altLoc,
// 				       SelRes->seqNum,
// 				       SelRes->GetInsCode(),
// 				       SelRes->GetChainID());

	       int idx = get_asc_index(res_selection[iat]);

	       if (idx >= 0) {
		  double sigma = geom[idr].second.plane_restraint[ip].dist_esd(irest_at);
		  if (sigma > 0) {
		     std::pair<int, double> idx_sigma_pair(idx, sigma);
		     if (alt_conf.empty()) {
			if (all_altconfs_the_same) {
			   idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
			} else {
			   idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
			   for (unsigned int ialtconf=0; ialtconf<altconfs.size(); ialtconf++) {
			      if (! altconfs[ialtconf].empty()) {
				 idx_and_sigmas[altconfs[ialtconf]].push_back(idx_sigma_pair);
			      }
			   }
			}
		     } else {
			idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
		     }
		  }
	       }
	    }
	 }
      }

      std::map<std::string, std::vector <std::pair<int, double> > >::const_iterator it;
      for (it=idx_and_sigmas.begin(); it != idx_and_sigmas.end(); it++) {
	 if (it->second.size() > 3 ) {
	    std::vector<int> pos(it->second.size());
	    for (unsigned int i=0; i<it->second.size(); i++) pos[i] = it->second[i].first;
	    std::vector<bool> fixed_flags = make_fixed_flags(pos);
	    add_plane(it->second, fixed_flags);
	    n_plane_restr++;
	 }
      }

   }
   return n_plane_restr; 
}

int
coot::restraints_container_t::add_planes_as_improper_dihedrals(int idr, mmdb::PPAtom res_selection,
                                           int i_no_res_atoms,
                                           mmdb::PResidue SelRes,
                                           const protein_geometry &geom) {
   int n_impropers = 0;
   int index1, index2, index3, index4;

   std::vector<std::string> string_atom_names(i_no_res_atoms);
   for (int iat=0; iat<i_no_res_atoms; iat++)
      string_atom_names[iat] = res_selection[iat]->name;

   for (unsigned int ic=0; ic<geom[idr].second.improper_dihedral_restraint.size(); ic++) {

      const dict_improper_dihedral_restraint_t &dict_restraint = geom[idr].second.improper_dihedral_restraint[ic];
      
      if (true) {
         for (int iat1=0; iat1<i_no_res_atoms; iat1++) {
            const std::string &pdb_atom_name1 = string_atom_names[iat1];
            if (pdb_atom_name1 == dict_restraint.atom_id_1_4c()) {

               for (int iat2=0; iat2<i_no_res_atoms; iat2++) {
                  const std::string &pdb_atom_name2 = string_atom_names[iat2];
                  if (pdb_atom_name2 == dict_restraint.atom_id_2_4c()) {

                     for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
                        const std::string &pdb_atom_name3 = string_atom_names[iat3];
                        if (pdb_atom_name3 == dict_restraint.atom_id_3_4c()) {

                           for (int iat4=0; iat4<i_no_res_atoms; iat4++) {
                              const std::string &pdb_atom_name4 = string_atom_names[iat4];
                              if (pdb_atom_name4 == dict_restraint.atom_id_4_4c()) {
                                 
				 std::string alt_conf_1 = res_selection[iat1]->altLoc;
				 std::string alt_conf_2 = res_selection[iat2]->altLoc;
				 std::string alt_conf_3 = res_selection[iat3]->altLoc;
				 std::string alt_conf_4 = res_selection[iat4]->altLoc;

                                 if (((alt_conf_1 == alt_conf_4) || (alt_conf_1 == "")) &&
                                     ((alt_conf_2 == alt_conf_4) || (alt_conf_2 == "")) &&
                                     ((alt_conf_3 == alt_conf_4) || (alt_conf_3 == ""))) {

				    res_selection[iat1]->GetUDData(udd_atom_index_handle, index1);
				    res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
				    res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);
				    res_selection[iat4]->GetUDData(udd_atom_index_handle, index4);

				    std::vector<bool> fixed_flags =
				    make_fixed_flags(index4, index1, index2, index3);
				    float sigma = dict_restraint.sigma;
                                    simple_restraint sr(IMPROPER_DIHEDRAL_RESTRAINT,
                                                        index1, index2, index3, index4,
                                                        sigma, fixed_flags);
                                    restraints_vec.push_back(sr); // push_back_restraint()
				    n_impropers++;
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_impropers;

}



void
coot::restraints_container_t::convert_plane_restraints_to_improper_dihedral_restraints(const std::vector<std::pair<int, double> > &atom_index_and_sigma,
										       const std::vector<bool> &fixed_atom_flags) {

   if (atom_index_and_sigma.size() == 4) {

#if 0
      double sigma = 0.01;
      simple_restraint sr(IMPROPER_DIHEDRAL_RESTRAINT,
                          atom_and_index_sigma[3].first, atom_and_index_sigma[0].first, atom_and_index_sigma[1].first, atom_and_index_sigma[2].first,
                          sigma, fixed_flags);
      restraints_vec.push_back(sr);
#endif

   } else {

      // is 5 (planar peptide restraints from link restraints). We shouldn't get here by any other
      // (say monomer dictionary) means.

      if (atom_index_and_sigma.size() == 5) {

      } else {

         std::cout << "ERROR:: in convert_plane_restraints_to_improper_dihedral_restraints() "
                   << atom_index_and_sigma.size() << std::endl;
      }
   }

}



int
coot::restraints_container_t::add_rama(const coot::rama_triple_t &rt, const coot::protein_geometry &geom) {

   return add_rama(rt.link_type, rt.r_1, rt.r_2, rt.r_3,
		   rt.fixed_1, rt.fixed_2, rt.fixed_3,
		   geom);
}


// make RAMACHANDRAN_RESTRAINTs, not TORSION_RESTRAINTs these days.
int
coot::restraints_container_t::add_rama(std::string link_type,
				       mmdb::PResidue prev_res,
				       mmdb::PResidue this_res,
				       mmdb::PResidue post_res,
				       bool is_fixed_first,
				       bool is_fixed_second,
				       bool is_fixed_third,
				       const coot::protein_geometry &geom) {

   // Old notes:
   // TRANS    psi      1 N      1 CA     1 C      2 N
   // TRANS    phi      1 C      2 N      2 CA     2 C
   // TRANS    omega    1 CA     1 C      2 N      2 CA
   //
   // New assignements:
   // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
   // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C)
   //
   // So Rama_atoms in this order:
   //   0       1        2      3         4
   // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)


   //std::cout << "DEBUG:: --------- :: Adding RAMA phi_psi_restraints_type for " << this_res
   //     << std::endl;

   int n_rama = 0;

   mmdb::PPAtom prev_sel;
   mmdb::PPAtom this_sel;
   mmdb::PPAtom post_sel;
   int n_first_res_atoms, n_second_res_atoms, n_third_res_atoms;

   prev_res->GetAtomTable(prev_sel,  n_first_res_atoms);
   this_res->GetAtomTable(this_sel, n_second_res_atoms);
   post_res->GetAtomTable(post_sel,  n_third_res_atoms);

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
      // throw 
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
      // throw
   }
   if (n_third_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
      // throw
   }

   std::vector<bool> fixed_flag(5, 0);
   if (is_fixed_first) {
      fixed_flag[0] = 1;
   }
   if (is_fixed_second) {
      fixed_flag[1] = 1;
      fixed_flag[2] = 1;
      fixed_flag[3] = 1;
   }
   if (is_fixed_third) {
      fixed_flag[4] = 1;
   }
   std::vector<mmdb::Atom *> rama_atoms(5);
   for (int ir=0; ir<5; ir++)
      rama_atoms[ir] = 0;
   
   for (int i=0; i<n_first_res_atoms; i++) {
      std::string atom_name(prev_sel[i]->name);
      if (atom_name == " C  ")
	 rama_atoms[0] = prev_sel[i];
   }
   for (int i=0; i<n_second_res_atoms; i++) {
      std::string atom_name(this_sel[i]->name);
      if (atom_name == " N  ")
	 rama_atoms[1] = this_sel[i];
      if (atom_name == " CA ")
	 rama_atoms[2] = this_sel[i];
      if (atom_name == " C  ")
	 rama_atoms[3] = this_sel[i];
   }
   for (int i=0; i<n_third_res_atoms; i++) {
      std::string atom_name(post_sel[i]->name);
      if (atom_name == " N  ")
	 rama_atoms[4] = post_sel[i];
   }

   if (rama_atoms[0] && rama_atoms[1] && rama_atoms[2] && 
       rama_atoms[3] && rama_atoms[4]) {

      std::vector<int> atom_indices(5, -1);
      for (int i=0; i<5; i++) {
// 	 atom_indices[i] = get_asc_index(rama_atoms[i]->name,
// 					 rama_atoms[i]->altLoc,
// 					 rama_atoms[i]->residue->seqNum,
// 					 rama_atoms[i]->GetInsCode(),
// 					 rama_atoms[i]->GetChainID());
	 atom_indices[i] = get_asc_index(rama_atoms[i]);
      }

      if ( (atom_indices[0] != -1) && (atom_indices[1] != -1) && (atom_indices[2] != -1) && 
	   (atom_indices[3] != -1) && (atom_indices[4] != -1)) { 


	 std::string zort = zo_rama.get_residue_type(this_res->GetResName(),
						     post_res->GetResName());

	 if (true)
	    // std::cout << "INFO:: Adding Ramachandran restraint "
	    //           << "type " << std::setw(6) << zort << " for " << residue_spec_t(this_res)
	    //           << " " << this_res->GetResName() << " "
	    //           << "fixed: "
	    //           << fixed_flag[0] << " " << fixed_flag[1] << " "
	    //           << fixed_flag[2] << " " << fixed_flag[3] << " "
	    //           << fixed_flag[4]
	    //           << std::endl;
	    logger.log(log_t::INFO, "Adding Ramachandran restraint type " + zort + " for " +
                       residue_spec_t(this_res).format() + " " + std::string(this_res->GetResName()) +
                       " fixed: " + std::to_string(fixed_flag[0]) + " " + std::to_string(fixed_flag[1]) + " " +
                       std::to_string(fixed_flag[2]) + " " + std::to_string(fixed_flag[3]) + " " + std::to_string(fixed_flag[4]));

	 add(RAMACHANDRAN_RESTRAINT,
	     zort,
	     atom_indices[0], atom_indices[1], atom_indices[2],
	     atom_indices[3], atom_indices[4], fixed_flag);
	 n_rama++;
      }
   }
   // std::cout << "returning..." << n_rama << std::endl;
   return n_rama; 
}


// return in millisecs
// double 
// coot::restraints_container_t::time_diff(const timeval &current, const timeval &start) const {

//    double d = current.tv_sec - start.tv_sec;
//    d *= 1000.0;
//    d += double(current.tv_usec - start.tv_usec)/1000.0;
//    return d;


coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_conventional(int selHnd,
							   const coot::protein_geometry &geom) const {

   float dist_crit = 3.0; // if atoms in different residues are closer
			  // than this, then they are considered
			  // bonded (potentially).
   
   // First add "linear" links
   // 
   coot::bonded_pair_container_t c = bonded_residues_by_linear(selHnd, geom);

   // Now, are there any other links?
   //
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 1) {
      for (int ii=0; ii<nSelResidues; ii++) {
	 for (int jj=0; jj<nSelResidues; jj++) {
	    if (jj>ii) {
	       if (! c.linked_already_p(SelResidue[ii], SelResidue[jj])) {
		  std::pair<bool, float> d = closest_approach(SelResidue[ii], SelResidue[jj]);
		  if (d.first) {
		     if (d.second < dist_crit) {
                        if (false)
                           std::cout << "####################### find_link_type_compli() called from " << __FUNCTION__  << "()"
                                     << std::endl;
			std::pair<std::string, bool> l =
			   find_link_type_complicado(SelResidue[ii], SelResidue[jj], geom);
			if (l.first != "") {

			   // Eeek!  Fill me? [for 0.7]
			   
			} 
		     } 
		  }
	       }
	    }
	 }
      }
   }
   return c;
}
   
void
coot::restraints_container_t::construct_non_bonded_contact_list_conventional() {

   // So first, I need a method/list to determine what is not-bonded
   // to what.
   // 

//    std::cout << "bonded list:" << std::endl;
//    std::cout << "--------------------------------------------------\n";
//    for (int i=0; i<bonded_atom_indices.size(); i++) { 
//       std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
//       for (int j=0; j<bonded_atom_indices[i].size(); j++) { 
// 	 std::cout << bonded_atom_indices[i][j] << " ";
//       } 
//       std::cout << std::endl;
//    } 
//    std::cout << "--------------------------------------------------\n";

  // Now we need to know which indices into the mmdb::PPAtom atoms are in
  // the moving set (rather than the flanking atoms).
  // 
  std::vector<std::vector<int> > non_bonded_atom_indices;
  //
  short int was_bonded_flag;
  // Note: bonded_atom_indices is sized to n_atoms in init_shared_post().
  non_bonded_atom_indices.resize(bonded_atom_indices.size());

  // Set up some mmdb::PPAtom things needed in the loop:
  mmdb::PPAtom res_selection_local;
  int n_res_atoms;
  mmdb::PPAtom res_selection_local_inner;
  int n_res_atoms_inner;
  int atom_index, atom_index_inner;
  int ierr;

  // iar = i active residue, nSelResidues_active is class variable

  // std::cout << "INFO:: There are " << nSelResidues_active << " active residues\n";
  for (int iar=0; iar<nSelResidues_active; iar++) { 

     SelResidue_active[iar]->GetAtomTable(res_selection_local, n_res_atoms);
     // std::cout << "There are " << n_res_atoms << " active atoms in this active residue\n";

     for (int iat=0; iat<n_res_atoms; iat++) { 

	// set atom_index
	ierr = res_selection_local[iat]->GetUDData(udd_atom_index_handle, atom_index);
	if (ierr != mmdb::UDDATA_Ok) { 
	   std::cout << "ERROR:: in getting UDDATA res_selection_local, ierr=" 
		     << ierr << " "
		     << res_selection_local[iat]->GetSeqNum() << " " 
		     << res_selection_local[iat]->GetAtomName() << " \n";
	}
	
	bool matched_oxt = false;
	if (have_oxt_flag) {
	   if (std::string(res_selection_local[iat]->name) == " OXT") {  // PDBv3 FIXME
	      matched_oxt = true;
	   } else { 
	      matched_oxt = false;
	   }
	}

	if (! matched_oxt) { 

	   // For each of the bonds of atom with atom_index index
	   // we need to check if bonded_atom_indices[atom_index][j]
	   // matches any atom index of the active atoms:

	   for (int jar=0; jar<nSelResidues_active; jar++) { 
	   
	      SelResidue_active[jar]->GetAtomTable(res_selection_local_inner, 
						   n_res_atoms_inner);
	   
	      for (int jat=0; jat<n_res_atoms_inner; jat++) { 
	      
		 // set atom_index_inner
		 ierr = res_selection_local_inner[jat]->GetUDData(udd_atom_index_handle,
                                                                  atom_index_inner);

		 if (atom_index == atom_index_inner) { 
		    // std::cout << "skipping same index " << std::endl;
		 } else {
		    
// 		    std::cout << "DEBUG:: checking bond pair " << atom_index << " " 
// 			      << atom_index_inner << " " 
// 			      << atom[atom_index]->name << " " << atom[atom_index]->GetSeqNum() << "    " 
// 			      << atom[atom_index_inner]->name << " " << atom[atom_index_inner]->GetSeqNum()
//          		      << std::endl;
	      
		    was_bonded_flag = 0;

		    // PDBv3 FIXME
		    if (have_oxt_flag) 
		       if (! strcmp(res_selection_local_inner[jat]->name, " OXT")) // matched
			  matched_oxt = true;

		    if (! matched_oxt) { 

		       std::set<int>::const_iterator it;
		       for (it=bonded_atom_indices[atom_index].begin();
			    it!=bonded_atom_indices[atom_index].end(); ++it) {
			  if (*it == atom_index_inner) {
			     was_bonded_flag = true;
			     break;
			  }
		       }
		       
		       if (was_bonded_flag == 0) {
			  non_bonded_atom_indices[atom_index].push_back(atom_index_inner);
		       }
		    }
		 }
	      }
	   }
	}
     }
  }

  if (false) {
     std::cout << "--------------------------------------------------\n";
     std::cout << "   conventional non-bonded list (unfiltered by distance):" << std::endl;
     std::cout << "--------------------------------------------------\n";
     for (unsigned int i=0; i<non_bonded_atom_indices.size(); i++) { 
	std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
	for (unsigned int j=0; j<non_bonded_atom_indices[i].size(); j++) { 
	   std::cout << non_bonded_atom_indices[i][j] << " ";
	} 
	std::cout << std::endl;
     }
     std::cout << "--------------------------------------------------\n";
  }

  filter_non_bonded_by_distance(non_bonded_atom_indices, 8.0);
}

// This function is called by make_non_bonded_contact_restraints()
//
void
coot::restraints_container_t::construct_non_bonded_contact_list_by_res_vec(const coot::bonded_pair_container_t &bpc,
									   const coot::protein_geometry &geom) {

#ifdef HAVE_CXX_THREAD
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();
#endif

   // How frequently does this function get called? - needs optimizing

   //  on a whole chain:
   //  8 -> 2.9 s
   // 11 -> 3.1 s
   //
   const double dist_crit = 8.0;
   const double dist_crit_sqrd = dist_crit * dist_crit;

   filtered_non_bonded_atom_indices.resize(bonded_atom_indices.size());

   if (false) { // debug
      std::cout << "DEBUG:: construct_non_bonded_contact_list_by_res_vec ::::::::::::::::" << std::endl;
      for (unsigned int i=0; i<bpc.size(); i++)
	 std::cout << "   "
		   << coot::residue_spec_t(bpc[i].res_1) << " "
		   << coot::residue_spec_t(bpc[i].res_2) << " "
		   << bpc[i].is_fixed_first << " " 
		   << bpc[i].is_fixed_second << " " 
		   << std::endl;

      std::cout << "--------------- debug:: bonded_atom_indices size "
		<< bonded_atom_indices.size() << std::endl;
      std::cout << "--------------- debug:: n_atoms " << n_atoms << std::endl;

      std::cout << "Bonded atom indices:" << std::endl;
      for (unsigned int i=0; i<bonded_atom_indices.size(); i++) {
	 std::cout << "  " << i << " " << atom_spec_t(atom[i]) << " " << bonded_atom_indices[i].size()
		   << " |  ";
	 // interate through a set now.
	 // for (unsigned int j=0; j<bonded_atom_indices[i].size(); j++)
	 // std::cout << " " << bonded_atom_indices[i][j];
	 std::cout << "\n";
      }
   }

   // Yes, this adds symmetry to filtered_non_bonded_atom_indices, i.e.
   // filtered_non_bonded_atom_indices[0] contains 1
   // filtered_non_bonded_atom_indices[1] contains 0
   // 
   for (unsigned int i=0; i<bonded_atom_indices.size(); i++) {

      // This is a hack.  It removes OXT from all NBCs.  Not the Right Way
      // 
      bool matched_oxt = false;
      if (have_oxt_flag) {
	 if (std::string(atom[i]->name) == " OXT") {  // PDBv3 FIXME
	    matched_oxt = true;
	 }
      }

      for (unsigned int j=0; j<bonded_atom_indices.size(); j++) {

	 if (i != j) {

	    if (have_oxt_flag) {
	       if (std::string(atom[j]->name) == " OXT") {  // PDBv3 FIXME
		  matched_oxt = true;
	       }
	    }

	    if (false)
	       std::cout << "moving->moving: here with atoms "
			 << atom_spec_t(atom[i]) <<  " " << atom_spec_t(atom[j])
			 << " have_oxt_flag: " << have_oxt_flag
			 << " matched_oxt: " << matched_oxt << std::endl;

	    if (! matched_oxt) {

	       // In this section, we don't want NCBs within or to fixed
	       // residues (including the flanking residues), so if both
	       // atoms are in residues that are not in residue_vec, then
	       // we don't add a NCB for that atom pair.

	       // bonded_atom_indices contains indices of atoms that
	       // are angle-related (not just directly bonded)
	       // 
	       //if (is_member_p(bonded_atom_indices[i], j)) {
	       if (bonded_atom_indices[i].find(j) != bonded_atom_indices[i].end()) {

		  // debug bonded atoms
		  
	       } else {

		  // atom j is not bonded to atom i, is it close? (i.e. within dist_crit?)

		  // Go faster than this...
		  // clipper::Coord_orth pt1(atom[i]->x, atom[i]->y, atom[i]->z);
		  // clipper::Coord_orth pt2(atom[j]->x, atom[j]->y, atom[j]->z);
		  // double d = clipper::Coord_orth::length(pt1, pt2);

		  double xd(atom[i]->x - atom[j]->x);
		  double yd(atom[i]->y - atom[j]->y);
		  double zd(atom[i]->z - atom[j]->z);
		  double d_sqrd = xd*xd + yd*yd + zd*zd;
		  if (d_sqrd < dist_crit_sqrd) {
		     mmdb::Residue *r1 = atom[i]->residue;
		     mmdb::Residue *r2 = atom[j]->residue;

		     std::string alt_conf_1 = atom[i]->altLoc;
		     std::string alt_conf_2 = atom[j]->altLoc;

 		     if ((alt_conf_1 == alt_conf_2) ||
 			 (alt_conf_1.length() == 0) ||
 			 (alt_conf_2.length() == 0)) {

			if (is_a_moving_residue_p(r1) && is_a_moving_residue_p(r2)) {
			   filtered_non_bonded_atom_indices[i].push_back(j);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (false) { // debug - how many bonded atoms are we talking about here?
      int n = 0;
      for (int iat=0; iat<n_atoms; iat++) {
	 n += bonded_atom_indices[iat].size();
      }
      std::cout << "DEBUG:: " << n << " bonded atom pairs to check " << std::endl;
   }

   // now add NBC restraints between atoms that are moving and atoms
   // of the neighbour residues.
   // 
   for (int iat=0; iat<n_atoms; iat++) {
      if (bonded_atom_indices[iat].size()) { 
	 mmdb::Residue *bonded_atom_residue = atom[iat]->residue;
	 for (int jat=0; jat<n_atoms; jat++) {

	    if (iat != jat) {

	       mmdb::Residue *other_atom_residue = atom[jat]->residue;
	       if (bonded_atom_residue != other_atom_residue) {

		  if (is_a_moving_residue_p(bonded_atom_residue) &&
		      ! is_a_moving_residue_p(other_atom_residue)) {

		     bool matched_oxt = false;
		     if (have_oxt_flag) {
			if (std::string(atom[jat]->name) == " OXT") {  // PDBv3 FIXME
			   matched_oxt = true;
			} else { 
			   matched_oxt = false;
			}
		     }

		     if (false)
			std::cout << "moving->non-moving: here with atom "
				  << atom_spec_t(atom[jat]) << " have_oxt_flag: "
				  << have_oxt_flag << " matched_oxt: " << matched_oxt
				  << std::endl;

		     if (! matched_oxt) {

			bonded_pair_match_info_t mi = 
			   bpc.match_info(bonded_atom_residue, other_atom_residue);

			if (! mi.state) {

			   // Simple part, the residues were not bonded to each other.

			   // if (! is_member_p(bonded_atom_indices[iat], jat)) {
			   if (bonded_atom_indices[iat].find(jat) != bonded_atom_indices[iat].end()) {

			      // atom j is not bonded to atom i, is it close? (i.e. within dist_crit?)
			      clipper::Coord_orth pt1(atom[iat]->x, atom[iat]->y, atom[iat]->z);
			      clipper::Coord_orth pt2(atom[jat]->x, atom[jat]->y, atom[jat]->z);
			      double d = clipper::Coord_orth::length(pt1, pt2);
			      if (d < dist_crit) {
				 if (false)
				    std::cout << " ///////////////////////////// NBC   here "
					      << iat << " " << jat << " "
					      << coot::atom_spec_t(atom[iat]) << " "
					      << coot::atom_spec_t(atom[jat]) << " "
					      << std::endl;
				 filtered_non_bonded_atom_indices[iat].push_back(jat);
			      }
			   }
			   
			} else {

			   // they were bonded to each other.
			   
			   // add to filtered_non_bonded_atom_indices (which is a class variable)
			
			   if (mi.swap_needed)
			      construct_nbc_for_moving_non_moving_bonded(jat, iat, mi.link_type, geom);
			   else
			      construct_nbc_for_moving_non_moving_bonded(iat, jat, mi.link_type, geom);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (false) {
      std::cout << "--------------------------------------------------\n";
      std::cout << "  res-vec non-bonded list:" << std::endl;
      std::cout << "--------------------------------------------------\n";
      for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
	 std::cout << i << "  " << atom_spec_t(atom[i]) << " ";
	 for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) { 
	    std::cout << filtered_non_bonded_atom_indices[i][j] << " "
	       //  << atom_spec_t(atom[filtered_non_bonded_atom_indices[i][j]])
		      << " ";
	    if (j%20==0)
	       if (j > 0)
		  if (j != (filtered_non_bonded_atom_indices[i].size()-1))
		     std::cout << "\n          ";
	 }
	 std::cout << std::endl;
      }
      std::cout << "--------------------------------------------------\n";
   }

#ifdef HAVE_CXX_THREAD

   end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
   // std::cout << "INFO:: nbc computation " << "elapsed time: " << elapsed_seconds.count() << "s\n";
   logger.log(log_t::INFO, "nbc computation elapsed time:", elapsed_seconds.count(), "s");

#endif // HAVE_CXX_THREAD

}


// Add non-bonded contacts for atoms that are in residues that are
// bonded to each other.  Atom iat is in a moving residue and atom jat
// is in a residue that is not moving.
// 
// Add to filtered_non_bonded_atom_indices (which is a class variable).
// 
void
coot::restraints_container_t::construct_nbc_for_moving_non_moving_bonded(unsigned int iat, unsigned int jat,
									 const std::string &link_type,
									 const coot::protein_geometry &geom) {

   // We dont know if res_1 is in the moving residue or not.  We do
   // know that res_1 and res_2 are in the correct order for the given
   // link_type link.
   // 
   mmdb::Residue *res_1 = atom[iat]->residue;
   mmdb::Residue *res_2 = atom[jat]->residue;

   dictionary_residue_link_restraints_t link = geom.link(link_type);
   // std::cout << "link: " << link.link_id << " " << link.link_bond_restraint.size() << std::endl;
   if (! link.empty()) {
      std::string atom_name_1 = atom[iat]->name;
      std::string atom_name_2 = atom[jat]->name;
      bool add_it = true;
      for (unsigned int i=0; i<link.link_bond_restraint.size(); i++) {
	 if (atom_name_1 == link.link_bond_restraint[i].atom_id_1_4c() && 
	     atom_name_2 == link.link_bond_restraint[i].atom_id_2_4c()) {
	    add_it = false;
	    break;
	 }
      }
      for (unsigned int i=0; i<link.link_angle_restraint.size(); i++) { 
	 if (atom_name_1 == link.link_angle_restraint[i].atom_id_1_4c() && 
	     atom_name_2 == link.link_angle_restraint[i].atom_id_3_4c()) {
	    add_it = false;
	    break;
	 }
      }
      for (unsigned int i=0; i<link.link_torsion_restraint.size(); i++) { 
	 if (atom_name_1 == link.link_torsion_restraint[i].atom_id_1_4c() && 
	     atom_name_2 == link.link_torsion_restraint[i].atom_id_4_4c()) {
	    add_it = false;
	    break;
	 }
      }
      if (add_it) {

	 filtered_non_bonded_atom_indices[iat].push_back(jat);
	 
	 if (0) { // debug.
	    clipper::Coord_orth pt1(atom[iat]->x, atom[iat]->y, atom[iat]->z);
	    clipper::Coord_orth pt2(atom[jat]->x, atom[jat]->y, atom[jat]->z);
	    double d = sqrt((pt1-pt2).lengthsq());
		     
	    std::cout << "moving-non-moving: adding filtered non-bonded atom indices: " 
		      << atom_spec_t(atom[iat]) << " to  " 
		      << atom_spec_t(atom[jat]) << " dist: " << d
		      << std::endl;
	 }
	 
      } else {

	 if (0) 
	    std::cout << "moving-non-moving: REJECT filtered non-bonded atom indices: " 
		      << atom_spec_t(atom[iat]) << " to  " 
		      << atom_spec_t(atom[jat]) 
		      << std::endl;
      } 
   }
}

void 
coot::restraints_container_t::filter_non_bonded_by_distance(const std::vector<std::vector<int> > &non_bonded_atom_indices, double dist) { 

   filtered_non_bonded_atom_indices.resize(non_bonded_atom_indices.size());

   mmdb::Atom *atom_1;
   mmdb::Atom *atom_2;
   double dist2;
   double dist_lim2 = dist*dist;
   int i_at_ind;
   
   for (unsigned int i=0; i<non_bonded_atom_indices.size(); i++) { 
      for (unsigned int j=0; j<non_bonded_atom_indices[i].size(); j++) {
	 
	 atom_1 = atom[i];
	 atom_2 = atom[non_bonded_atom_indices[i][j]];

// 	 dist2 = clipper::Coord_orth::lengthsq(clipper::Coord_orth(atom_1->x, atom_1->y, atom_1->z), 
// 					       clipper::Coord_orth(atom_2->x, atom_2->y, atom_2->z));

	 dist2 = (clipper::Coord_orth(atom_1->x, atom_1->y, atom_1->z) -
		  clipper::Coord_orth(atom_2->x, atom_2->y, atom_2->z)).lengthsq();

	 if (dist2 < dist_lim2) { 
// 	    std::cout << "accepting non-bonded contact between " << atom_1->GetSeqNum() 
// 		      << " " << atom_1->name << " and " << atom_2->GetSeqNum() 
// 		      << " " << atom_2->name  << "\n";
	    atom_2->GetUDData(udd_atom_index_handle, i_at_ind); // sets i_at_ind.
	    filtered_non_bonded_atom_indices[i].push_back(i_at_ind);
	 } else { 
// 	    std::cout << "          reject non-bonded contact between " << atom_1->GetSeqNum() 
// 		      << " " << atom_1->name  << " and " << atom_2->GetSeqNum() 
// 		      << " " << atom_2->name << " rejected by distance\n";
	 } 
      }
   }
}
