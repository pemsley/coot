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
#ifdef HAVE_GSL


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

zo::rama_table_set coot::restraints_container_t::zo_rama;



coot::restraints_container_t::~restraints_container_t() {
   if (from_residue_vector) {
      if (atom) {
	 // this is constructed manually.

	 // Oh we can't do this here because we copy the
	 // restraints in simple_refine_residues() and that
	 // shallow copies the atom pointer - the original
	 // restriants go out of scope and call this destructor.
	 //
	 // We need a new way to get rid of atom - c.f. the
	 // linear/conventional way?
	 //
	 // delete [] atom;
	 // atom = NULL;
      }
   } else {
      // member data item mmdb::PPAtom atom is constructed by an
      // mmdb SelectAtoms()/GetSelIndex() (which includes
      // flanking atoms).
      // 20081207: don't do this here now - because the
      // memory/selection is deleted again in
      // clear_up_moving_atoms(). It *should* be done here of
      // course, but we'll save that for the future.
      //
      //if (atom) {
      // mol->DeleteSelection(SelHnd_atom);
      // atom = NULL;
      // }
   }
}

// iend_res is inclusive, so that 17,17 selects just residue 17.
//   have_disulfide_residues: other residues are included in the
//				residues_mol for disphide restraints.
// 
coot::restraints_container_t::restraints_container_t(int istart_res_in, int iend_res_in,
						     bool have_flanking_residue_at_start,
						     bool have_flanking_residue_at_end,
						     short int have_disulfide_residues,
						     const std::string &altloc,
						     const std::string &chain_id,
						     mmdb::Manager *mol_in, 
						     const std::vector<coot::atom_spec_t> &fixed_atom_specs,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {

   init(true);
   are_all_one_atom_residues = false;
   init_from_mol(istart_res_in, iend_res_in, 
		 have_flanking_residue_at_start, 
		 have_flanking_residue_at_end,
		 have_disulfide_residues,
		 altloc,
		 chain_id, mol_in, fixed_atom_specs);
}

// Used in omega distortion graph
// 
coot::restraints_container_t::restraints_container_t(atom_selection_container_t asc_in,
						     const std::string &chain_id,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {
   init(true);
   mol = asc_in.mol;
   are_all_one_atom_residues = false;

   istart_res = 999999;
   iend_res = -9999999;

   mmdb::PResidue *SelResidues = NULL;
   int nSelResidues;

   // -------- Find the max and min res no -----------------------------
   int selHnd = mol->NewSelection();
   mol->Select(selHnd, mmdb::STYPE_RESIDUE, 1,
	       chain_id.c_str(),
	       mmdb::ANY_RES, "*",
	       mmdb::ANY_RES, "*",
	       "*",  // residue name
	       "*",  // Residue must contain this atom name?
	       "*",  // Residue must contain this Element?
	       "*",  // altLocs
	       mmdb::SKEY_NEW // selection key
	       );
   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

   int resno;
   for (int ires=0; ires<nSelResidues; ires++) {
      resno = SelResidues[ires]->GetSeqNum();
      if (resno < istart_res)
	 istart_res = resno;
      if (resno > iend_res)
	 iend_res = resno;
   }
   mol->DeleteSelection(selHnd);
   // 
   // -------- Found the max and min res no -----------------------------

   // -------------------------------------------------------------------
   // Set class variables atom to the selection that includes the
   // chain (which is not the same as the input atom selection)
   //
   int SelHnd = mol->NewSelection();
   atom = NULL;
   mol->SelectAtoms(SelHnd, 0,
		    chain_id.c_str(),
		    mmdb::ANY_RES, // starting resno, an int
		    "*", // any insertion code
		    mmdb::ANY_RES, // ending resno
		    "*", // ending insertion code
		    "*", // any residue name
		    "*", // atom name
		    "*", // elements
		    "*"  // alt loc.
		    );
   mol->GetSelIndex(SelHnd, atom, n_atoms);

   // -------------------------------------------------------------------

   initial_position_params_vec.resize(3*n_atoms);
   for (int i=0; i<n_atoms; i++) {
      initial_position_params_vec[3*i  ] = atom[i]->x; 
      initial_position_params_vec[3*i+1] = atom[i]->y; 
      initial_position_params_vec[3*i+2] = atom[i]->z;
      // std::cout << "    " << i << "  " << coot::atom_spec_t(atom[i]) << "\n";
   }
}

coot::restraints_container_t::restraints_container_t(mmdb::PResidue *SelResidues, int nSelResidues,
						     const std::string &chain_id,
						     mmdb::Manager *mol_in,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {

   init(true);
   are_all_one_atom_residues = false;

   std::vector<coot::atom_spec_t> fixed_atoms_dummy;
   int istart_res_l = 999999;
   int iend_res_l = -9999999;
   int resno;
   
   for (int i=0; i<nSelResidues; i++) { 
      resno = SelResidues[i]->seqNum;
      if (resno < istart_res_l)
	 istart_res_l = resno;
      if (resno > iend_res_l)
	 iend_res_l = resno;
   }
   
   short int have_flanking_residue_at_start = 0;
   short int have_flanking_residue_at_end = 0;
   short int have_disulfide_residues = 0;
   const char *chn = chain_id.c_str();

   // std::cout << "DEBUG:  ==== istart_res iend_res " << istart_res << " "
   // << iend_res << std::endl; 

   init_from_mol(istart_res_l, iend_res_l,
		 have_flanking_residue_at_start,
		 have_flanking_residue_at_end,
		 have_disulfide_residues, 
		 std::string(""), chn, mol_in, fixed_atoms_dummy);

}

coot::restraints_container_t::restraints_container_t(int istart_res_in, int iend_res_in,
						     short int have_flanking_residue_at_start,
						     short int have_flanking_residue_at_end,
						     short int have_disulfide_residues,
						     const std::string &altloc,
						     const std::string &chain_id,
						     mmdb::Manager *mol_in,
						     const std::vector<coot::atom_spec_t> &fixed_atom_specs,
						     const clipper::Xmap<float> *map_p_in,
						     float map_weight_in) : xmap_p(map_p_in) {

   init(true);
   init_from_mol(istart_res_in, iend_res_in, 		 
		 have_flanking_residue_at_start, 
		 have_flanking_residue_at_end,
		 have_disulfide_residues,
		 altloc,
		 chain_id, mol_in, fixed_atom_specs);
   are_all_one_atom_residues = false;
   map_weight = map_weight_in;
   include_map_terms_flag = 1;

}

bool
coot::residue_sorter(const std::pair<bool, mmdb::Residue *> &r1,
		     const std::pair<bool, mmdb::Residue *> &r2) {

   std::string chain_id_1 = r1.second->GetChainID();
   std::string chain_id_2 = r2.second->GetChainID();
   if (chain_id_1 < chain_id_2) {
      return true;
   } else {
      if (chain_id_1 > chain_id_2) {
	 return false;
      } else {
	 if (r1.second->index < r2.second->index) {
	    return true;
	 } else {
	    if (r1.second->index > r2.second->index) {
	       return false;
	    } else {
	       if (r1.second->GetSeqNum() < r2.second->GetSeqNum()) {
		  return true;
	       } else {
		  if (r1.second->GetSeqNum() > r2.second->GetSeqNum()) {
		     return false;
		  } else {
		     std::string ins_code_1 = r1.second->GetInsCode();
		     std::string ins_code_2 = r2.second->GetInsCode();
		     if (ins_code_1 < ins_code_2) {
			return true;
		     } else {
			if (ins_code_1 > ins_code_2) {
			   return false;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return false;
}


// 20081106 construct from a vector of residues, each of which
// has a flag attached that denotes whether or not it is a fixed
// residue (it would be set, for example in the case of flanking
// residues).
coot::restraints_container_t::restraints_container_t(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						     const std::vector<mmdb::Link> &links,
						     const coot::protein_geometry &geom,
						     mmdb::Manager *mol_in,
						     const std::vector<atom_spec_t> &fixed_atom_specs,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {

   istart_minus_flag = false; // used in make_flanking_atoms_rama_restraints
   iend_plus_flag = false;

   init(true);
   from_residue_vector = 1;
   are_all_one_atom_residues = false;

   std::vector<std::pair<bool,mmdb::Residue *> > residues_local;
   residues_local.reserve(residues.size());

   for(unsigned int i=0; i<residues.size(); i++)
      if (residues[i].second)
         residues_local.push_back(residues[i]);

   // now sort those residues so that polymer linking is easy
   // sorting function should return false at the end - because sorting function is called
   // with the same argument for left and right hand side, I think (std::sort testing for sanity?)
   std::sort(residues_local.begin(), residues_local.end(), residue_sorter);

   if (false)
      for (std::size_t i=0; i<residues_local.size(); i++)
	 std::cout << "    restraints_container_t() constructor: " << residue_spec_t(residues_local[i].second)
		   << " has index " << residues_local[i].second->index << std::endl;

   residues_vec = residues_local;
   init_from_residue_vec(residues_local, geom, mol_in, fixed_atom_specs);


}

coot::restraints_container_t::restraints_container_t(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						     const coot::protein_geometry &geom,
						     mmdb::Manager *mol_in,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {

   init(true);
   from_residue_vector = 1;
   are_all_one_atom_residues = false;

   std::vector<atom_spec_t> fixed_atom_specs;
   std::vector<std::pair<bool,mmdb::Residue *> > residues_local;
   residues_local.reserve(residues.size());

   for(unsigned int i=0; i<residues.size(); i++)
      if (residues[i].second)
         residues_local.push_back(residues[i]);

   residues_vec = residues_local;

   if (false) {
      std::cout << "debug:: in restraints_container_t() constructor with local residue size " << residues_local.size() << std::endl;
      for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
         if (residues_vec[ir].second) {
           std::cout << "INFO:: starting init_from_residue_vec() residue " << residues_vec[ir].second << std::endl;
         } else {
           std::cout << "ERROR:: starting init_from_residue_vec() NUll residue " << ir << std::endl;
         }
      }
   }

   init_from_residue_vec(residues_local, geom, mol_in, fixed_atom_specs);

}



// What are the rules for dealing with alt conf in flanking residues?
// 
// First we want to try to select only atom that have the same alt
// conf, if that fails to select atoms in flanking residues (and there
// are atoms in flanking residues (which we know from
// have_flanking_residue_at_end/start), then we should try an mmdb construction [""|"A"]
// (i.e. blank or "A").  If that fails, give up - it's a badly formed pdb file (it 
// seems to me).
// 
void
coot::restraints_container_t::init_from_mol(int istart_res_in, int iend_res_in,
					    bool have_flanking_residue_at_start,
					    bool have_flanking_residue_at_end,
					    short int have_disulfide_residues,
					    const std::string &altloc,
					    const std::string &chain_id,
					    mmdb::Manager *mol_in,
					    const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   init_shared_pre(mol_in);

   istart_res = istart_res_in;
   iend_res   = iend_res_in;
   chain_id_save = chain_id;

   // internal flags that mirror passed have_flanking_residue_at_* variables
   istart_minus_flag = have_flanking_residue_at_start;
   iend_plus_flag    = have_flanking_residue_at_end;

   int iselection_start_res = istart_res;
   int iselection_end_res   = iend_res;
   // std::cout << "start res range: " << istart_res << " " << iend_res << " " << chain_id << "\n";

   // Are the flanking atoms available in mol_in?
   // mol_in was constructed outside of this class, so the mol_in constructing
   // routine knows if they were there or not.
   // 
   if (have_flanking_residue_at_start) iselection_start_res--;
   if (have_flanking_residue_at_end)   iselection_end_res++;

   SelHnd_atom = mol->NewSelection();
   mol->SelectAtoms(SelHnd_atom,
		    0,
		    chain_id.c_str(),
		    iselection_start_res, "*",
		    iselection_end_res,   "*",
		    "*", // rnames
		    "*", // anames
		    "*", // elements
		    "*"  // altLocs 
		    );

   // set the mmdb::PPAtom atom (class variable) and n_atoms:
   // 
   mol->GetSelIndex(SelHnd_atom, atom, n_atoms);

   if (false) { // debugging;
      std::cout << "debug:: in init_from_mol() here are the " << fixed_atom_indices.size()
		<< " fixed_atom indices: \n";
      std::set<int>::const_iterator it;
      for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++)
	 std::cout << " " << *it;
      std::cout << "\n";

      for (int iat=0; iat<n_atoms; iat++)
	 std::cout << "   " << iat << "  "  << coot::atom_spec_t(atom[iat]) << "  with altloc :"
		   << altloc << ":" << std::endl;
   }

   bool debug = false;
   if (debug) {
      std::cout << "DEBUG:: Selecting residues in chain \"" << chain_id << "\" gives "
		<< n_atoms << " atoms " << std::endl;
      for (int iat=0; iat<n_atoms; iat++) {
	 std::cout << "   " << iat << " " << atom[iat]->name << " "  << atom[iat]->GetSeqNum()
		   << " " << atom[iat]->GetChainID() << std::endl;
      }
   }

   // debugging, you need to uncomment mmdb.h at the top and add coords to the linking
   // atom_selection_container_t tmp_res_asc = make_asc(mol_in);
   //    std::cout << "There are " << tmp_res_asc.n_selected_atoms
   //              << " atoms in tmp_res_asc\n";
//    for (int kk=0; kk<tmp_res_asc.n_selected_atoms; kk++) { 
//       std::cout << "In simple rest " << kk << " "
//                 << tmp_res_asc.atom_selection[kk] << "\n";
//    } 

//    std::cout << "INFO::" << n_atoms << " atoms selected from molecule for refinement" ;
//    std::cout << " (this includes fixed and flanking atoms)." << std::endl;

   if (n_atoms == 0) { 
      std::cout << "ERROR:: atom selection disaster:" << std::endl;
      std::cout << "   This should not happen" << std::endl;
      std::cout << "   residue range: " << iselection_start_res << " " 
		<< iselection_end_res << " chain-id \"" << chain_id << "\" " 
		<< "flanking flags: " << have_flanking_residue_at_start 
		<< " " << have_flanking_residue_at_end << std::endl;
   }

   init_shared_post(fixed_atom_specs); // clears fixed_atom_indices

   add_fixed_atoms_from_flanking_residues(have_flanking_residue_at_start,
					  have_flanking_residue_at_end,
					  iselection_start_res, iselection_end_res);

}

void
coot::restraints_container_t::init_shared_pre(mmdb::Manager *mol_in) {

   needs_reset = false;
   verbose_geometry_reporting = NORMAL;
   do_numerical_gradients_flag = false;
   have_oxt_flag = false; // set in mark_OXT()
   dist_crit_for_bonded_pairs = 3.0;
   // the smaller the alpha, the more like least squares
   geman_mcclure_alpha = 0.2; // Is this a good value? Talk to Rob.
   mol = mol_in;
   lennard_jones_epsilon = 0.1;
   cryo_em_mode = true;
   n_times_called = 0;
   n_small_cycles_accumulator = 0;
   m_s = 0;
   x = 0;
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   n_threads = 0;
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   log_cosh_target_distance_scale_factor = 3000.0;
}

void
coot::restraints_container_t::init_shared_post(const std::vector<atom_spec_t> &fixed_atom_specs) {


   bonded_atom_indices.resize(n_atoms);

   initial_position_params_vec.resize(3*n_atoms); 
   for (int i=0; i<n_atoms; i++) {
      initial_position_params_vec[3*i  ] = atom[i]->x; 
      initial_position_params_vec[3*i+1] = atom[i]->y; 
      initial_position_params_vec[3*i+2] = atom[i]->z; 
   }

   // Set the UDD have_bond_or_angle to initally all "not".  They get
   // set to "have" (1) in make_restraints (and functions thereof).
   // 
   // udd_handle becomes member data so that it can be used in
   // make_restraints() without passing it back (this function is part
   // of a constructor, don't forget).
   //
   // 20131213-PE: I dont see the point of udd_bond_angle.
   // 
   if (mol) { 
      udd_bond_angle = mol->RegisterUDInteger (mmdb::UDR_ATOM, "bond or angle");
      if (udd_bond_angle < 0) { 
	 std::cout << "ERROR:: can't make udd_handle in init_from_mol\n";
      } else { 
	 for (int i=0; i<n_atoms; i++) {
	    atom[i]->PutUDData(udd_bond_angle,0);
	 }
      }
   }

   // Set the UDD of the indices in the atom array (i.e. the thing
   // that get_asc_index returns)
   // 
   if (mol) {
      udd_atom_index_handle = mol->RegisterUDInteger ( mmdb::UDR_ATOM, "atom_array_index");
      if (udd_atom_index_handle < 0) {
	 std::cout << "ERROR:: can't make udd_handle in init_from_mol\n";
      } else {
	 for (int i=0; i<n_atoms; i++) {
	    atom[i]->PutUDData(udd_atom_index_handle,i);
	    // std::cout << "init_shared_post() atom " << atom_spec_t(atom[i])
	    // << " gets udd_atom_index_handle value " << i << std::endl;
	 }
      }
   }

   use_map_gradient_for_atom.resize(n_atoms, false);
   if (! from_residue_vector) {
      // convential way
      for (int i=0; i<n_atoms; i++) {
	 if (atom[i]->residue->seqNum >= istart_res &&
	     atom[i]->residue->seqNum <= iend_res) {
	    if (! is_hydrogen(atom[i]))
	       use_map_gradient_for_atom[i] = true;
	 } else {
	    use_map_gradient_for_atom[i] = false;
	 }
      }
   } else {
      // blank out the non moving atoms (i.e. flanking residues)
      for (int i=0; i<n_atoms; i++) {
	 mmdb::Residue *res_p = atom[i]->residue;
	 if (is_a_moving_residue_p(res_p)) {
	    if (! is_hydrogen(atom[i]))
	       use_map_gradient_for_atom[i] = true;
	 } else {
	    // std::cout << "blanking out density for atom " << i << std::endl;
	    use_map_gradient_for_atom[i] = false;
	 }
      }
   }

   // z weights:
   //
   atom_z_occ_weight.resize(n_atoms);
   std::vector<std::pair<std::string, int> > atom_list = coot::util::atomic_number_atom_list();
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      if (! at->isTer()) {
	 double z = coot::util::atomic_number(at->element, atom_list);
	 double weight = 1.0;
	 double occupancy = atom[i]->occupancy;
	 if (occupancy > 1.0) occupancy = 1.0;
	 if (cryo_em_mode) {
	    // is-side-chain? would be a better test
	    if (! is_main_chain_or_cb_p(at)) {
		  // std::cout << "downweighting atom " << coot::atom_spec_t(atom[i]) << std::endl;
		  weight = 0.2;
	    }
	    std::string at_name = atom[i]->name;
	    if (at_name == " O  ") {
	       weight = 0.4;
	    }
	 }

	 if (z < 0.0) {
	    std::cout << "WARNING:: init_shared_post() atom " << i << " " << atom_spec_t(atom[i])
		      << " Unknown element \"" << atom[i]->element << "\"" << std::endl;
	    z = 6.0; // as for carbon
	 }
	 atom_z_occ_weight[i] = weight * z * occupancy;
      }
   }

   // the fixed atoms:   
   // 
   assign_fixed_atom_indices(fixed_atom_specs); // convert from std::vector<atom_spec_t>
   				                // to std::vector<int> fixed_atom_indices;

   // blank out those atoms from seeing electron density map gradients

   std::set<int>::const_iterator it;
   for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++)
      use_map_gradient_for_atom[*it] = false;

   if (verbose_geometry_reporting == VERBOSE)
      for (int i=0; i<n_atoms; i++)
	 std::cout << atom[i]->name << " " << atom[i]->residue->seqNum << " "
		   << use_map_gradient_for_atom[i] << std::endl;

}

// uses fixed_atom_indices
void
coot::restraints_container_t::set_fixed_during_refinement_udd() {

   int uddHnd = mol->RegisterUDInteger(mmdb::UDR_ATOM , "FixedDuringRefinement");
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      // std::cout << "  setting fixed udd flag on atom " << atom_spec_t(at) << std::endl;
      // if (std::find(fixed_atom_indices.begin(), fixed_atom_indices.end(), i) == fixed_atom_indices.end())
      if (fixed_atom_indices.find(i) == fixed_atom_indices.end())
	 at->PutUDData(uddHnd, 0);
      else
	 at->PutUDData(uddHnd, 1);
   }

}


void
coot::restraints_container_t::init_from_residue_vec(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						    const coot::protein_geometry &geom,
						    mmdb::Manager *mol_in,
						    const std::vector<atom_spec_t> &fixed_atom_specs) {

   // This function is called from the constructor.
   // make_restraints() is called after this function by the user of this class.
   
   init_shared_pre(mol_in);

   if (residues.size() > 2000000) {
      std::cout << "ERROR:: in init_from_residue_vec() - memory error " << residues.size() << std::endl;
      return;
   }
   residues_vec.resize(residues.size());

   if (false) {
      for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
         if (residues_vec[ir].second) {
           std::cout << "INFO:: starting init_from_residue_vec() residue " << residues_vec[ir].second << std::endl;
         } else {
           std::cout << "ERROR:: starting init_from_residue_vec() NUll residue " << ir << std::endl;
         }
      }
   }

   // residues_vec = residues;
   for (std::size_t i=0; i<residues.size(); i++)
      if (residues[i].first == false) // i.e. is moving
	 residues_vec_moving_set.insert(residues[i].second);


   // Need to set class members mmdb::PPAtom atom and int n_atoms.
   // ...
   // 20090620: or do we?

   // debug:
   bool debug = false;
   if (debug) {
      for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
	 mmdb::PAtom *res_atom_selection = NULL;
	 int n_res_atoms;
	 residues_vec[ir].second->GetAtomTable(res_atom_selection, n_res_atoms);
	 std::cout << "debug:: =============== in init_from_residue_vec() residue "
		   << ir << " of " << residues_vec.size() << " has : "
		   << n_res_atoms << " atom " << std::endl;
	 std::cout << "debug:: =============== in init_from_residue_vec() residue "
		   << ir << " of " << residues_vec.size() << " " << residue_spec_t(residues_vec[ir].second)
		   << std::endl;
	 if (false)
	    for (int iat=0; iat<n_res_atoms; iat++) {
	       mmdb::Atom *at =  res_atom_selection[iat];
	       std::cout << "DEBUG:: in init_from_residue_vec: atom "
			 << iat << " of " << n_res_atoms << " \"" 
			 << at->name << "\" \"" << at->altLoc << "\" " 
			 << at->GetSeqNum() << " \"" << at->GetInsCode() << "\" \"" 
			 << at->GetChainID() << "\"" << std::endl;
	 }
      }
   }

   // what about adding the flanking residues?  How does the atom
   // indexing of that work when (say) adding a bond?

   if (false)
      std::cout << "debug::info in init_from_residue_vec() calling bonded_flanking_residues_by_residue_vector() "
		<< std::endl;

   float dist_crit = 5.3; // 20170924-PE was 3.0 but this made a horrible link in a tight turn
                          // (which I suspect is not uncommon) crazy-neighbour-refine-519.pdb
                          // for EMDB 6224.
                          // 520 was bonded to 522 in a neighb (3-residue) refine on 519.
                          // This function is called by init but not make_restraints.
                          // init doesn't set bonded_pairs_container (make_restraints does that).


   // fill neighbour_set from rnr (excluding residues of residues_vec):
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnr = residues_near_residues(residues_vec, mol, dist_crit);
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > neighbour_set;
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it_map;

   for(it_map=rnr.begin(); it_map!=rnr.end(); it_map++) {
      mmdb::Residue *r = it_map->first;
      const std::set<mmdb::Residue *> &s = it_map->second;
      std::set<mmdb::Residue *>::const_iterator it_set;
      for (it_set=s.begin(); it_set!=s.end(); it_set++) {
	 bool found = false;
	 for (std::size_t i=0; i<residues_vec.size(); i++) {
	    if (*it_set == residues_vec[i].second) {
	       found = true;
	       break;
	    }
	 }
	 if (! found) {
	    neighbour_set[r].insert(*it_set);
	 }
      }
   }


   bonded_pair_container_t bpc = bonded_flanking_residues_by_residue_vector(neighbour_set, geom);

   // internal variable non_bonded_neighbour_residues is set by this
   // function:
   set_non_bonded_neighbour_residues_by_residue_vector(neighbour_set, bpc, geom);

   if (false) { // debug

      std::cout << "############## neighbour set: " << std::endl;
      for(it_map=neighbour_set.begin(); it_map!=neighbour_set.end(); it_map++) {
	 const std::set<mmdb::Residue *> &s = it_map->second;
	 std::cout << "::: " << residue_spec_t(it_map->first) << std::endl;
	 std::set<mmdb::Residue *>::const_iterator it_set;
	 for (it_set=s.begin(); it_set!=s.end(); it_set++) {
	    std::cout << "      " << residue_spec_t(*it_set) << std::endl;
	 }
      }

      std::cout << "debug:: in init_from_residue_vec() here with these "
		<< non_bonded_neighbour_residues.size() << " non-bonded neighbours"
		<< std::endl;
      for (std::size_t jj=0; jj<non_bonded_neighbour_residues.size(); jj++) {
	 std::cout << "    " << residue_spec_t(non_bonded_neighbour_residues[jj]) << std::endl;
      }
   }
   

   // std::cout << "   DEBUG:: made " << bpc.size() << " bonded flanking pairs " << std::endl;

   // passed and flanking
   // 
   std::vector<mmdb::Residue *> all_residues;
   std::vector<mmdb::Residue *>::const_iterator it;
   for (unsigned int i=0; i<residues.size(); i++)
      all_residues.push_back(residues[i].second);

   // we don't need to calculate the NBC for all atoms - just these ones:
   n_atoms_limit_for_nbc = 0;
   for (unsigned int i=0; i<residues.size(); i++)
      n_atoms_limit_for_nbc += residues[i].second->GetNumberOfAtoms();


   // Include only the fixed residues, because they are the flankers,
   // the other residues are the ones in the passed residues vector.
   // We don't have members of bonded_pair_container_t that are both
   // fixed.
   //
   // 20151128 only include the residues once (test that they are not there first)
   //
   int n_bonded_flankers_in_total = 0; // debug/info counter
   for (unsigned int i=0; i<bpc.size(); i++) {
      if (bpc[i].is_fixed_first) {
	 it = std::find(all_residues.begin(),
			all_residues.end(),
			bpc[i].res_1);
	 if (it == all_residues.end()) { 
	    all_residues.push_back(bpc[i].res_1);
	    n_bonded_flankers_in_total++;
	 }
      } 
      if (bpc[i].is_fixed_second) {
	 it = std::find(all_residues.begin(),
			all_residues.end(),
			bpc[i].res_2);
	 if (it == all_residues.end()) {
	    all_residues.push_back(bpc[i].res_2);
	    n_bonded_flankers_in_total++;
	 }
      }
   }

   // Finally add the neighbour residues that are not bonded:
   for (unsigned int ires=0; ires<non_bonded_neighbour_residues.size(); ires++) {

      it = std::find(all_residues.begin(),
		     all_residues.end(),
		     non_bonded_neighbour_residues[ires]);
      if (it == all_residues.end())
	 all_residues.push_back(non_bonded_neighbour_residues[ires]);
   }

   if (0) {
      std::cout << "   DEBUG:: There are " << residues.size() << " passed residues and "
		<< all_residues.size() << " residues total (including flankers)"
		<< " with " << non_bonded_neighbour_residues.size()
		<< " non-bonded neighbours" << std::endl;
      std::cout << "These are the non-bonded neighbour residues " << std::endl;
      for (unsigned int i=0; i<non_bonded_neighbour_residues.size(); i++) { 
	 std::cout << "   " << residue_spec_t(non_bonded_neighbour_residues[i]) << std::endl;
      }
   }

   n_atoms = 0;

   // usualy this is reset in the following loop
   are_all_one_atom_residues = true;  // all refining residue (in residues) that is - do
                                      // not consider flanking residues here.
   // Now, is everything to be refined a water or MG or some such?  If
   // not (as is most often the case) are_all_one_atom_residues is
   // false.
   //
   // Something of a kludge because this will fail for residues with
   // alt-conf waters (but perhaps not, because before we get here,
   // we have done an atom selection so that there is only one alt
   // conf in each residue).
   // 
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      if (residues[ires].second->GetNumberOfAtoms() > 1) {
	 are_all_one_atom_residues = false;
	 break;
      }
   }

   for (unsigned int i=0; i<all_residues.size(); i++)
      n_atoms += all_residues[i]->GetNumberOfAtoms();
   atom = new mmdb::PAtom[n_atoms];
   int atom_index = 0;
   for (unsigned int i=0; i<all_residues.size(); i++) {
      mmdb::PPAtom residue_atoms = 0;
      int n_res_atoms;
      all_residues[i]->GetAtomTable(residue_atoms, n_res_atoms);
      for (int iat=0; iat<n_res_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 atom[atom_index] = at;
	 atom_index++;
      }
   }
   atom_is_metal.resize(n_atoms, false);
   for (int iat=0; iat<n_atoms; iat++) {
      atom_is_metal[iat] = geom.atom_is_metal(atom[iat]);
   }
   atom_is_hydrogen.resize(n_atoms, false);
   for (int iat=0; iat<n_atoms; iat++) {
      atom_is_hydrogen[iat] = is_hydrogen(atom[iat]);
   }

   // fill fixed_neighbours_set:
   //
   // to do this quicker, lets make a map of the residues in residues_vec that
   // are fixed - I could use residues_vec_moving_set above.
   //
   fixed_neighbours_set.clear();
   std::set<mmdb::Residue *> fixed_residue_set;
   for (std::size_t i=0; i<residues_vec.size(); i++)
      if (residues_vec[i].first)
	 fixed_residue_set.insert(residues_vec[i].second);

   if (false) { // debug
      for(it_map=rnr.begin(); it_map!=rnr.end(); it_map++) {
	 mmdb::Residue *r = it_map->first;
	 std::cout << "###### debugging rnr: Residue " << coot::residue_spec_t(r) << std::endl;
	 const std::set<mmdb::Residue *> &s = it_map->second;
	 std::set<mmdb::Residue *>::const_iterator it_set;
	 for (it_set=s.begin(); it_set!=s.end(); it_set++) {
	    mmdb::Residue *residue_neighb = *it_set;
	    std::cout << "###### debugging rnr:    Neighb: " << coot::residue_spec_t(residue_neighb)
		      << " " << residue_neighb << std::endl;
	 }
      }
   }

   
   for(it_map=rnr.begin(); it_map!=rnr.end(); it_map++) {
      mmdb::Residue *r = it_map->first;
      const std::set<mmdb::Residue *> &s = it_map->second;
      std::set<mmdb::Residue *>::const_iterator it_set;
      for (it_set=s.begin(); it_set!=s.end(); it_set++) {
	 mmdb::Residue *residue_neighb = *it_set;
	 // if residue_neigh is fixed and r is not then add residue_neigh as a neighbour of r
	 if (residues_vec_moving_set.find(residue_neighb) == residues_vec_moving_set.end()) {
	    if (residues_vec_moving_set.find(r) != residues_vec_moving_set.end()) {
	       fixed_neighbours_set[r].insert(residue_neighb);
	    }
	 }
      }
   }

   init_shared_post(fixed_atom_specs); // use n_atoms, fills fixed_atom_indices

   if (false) {
      std::cout << "---- after init_shared_post(): here are the "<< fixed_atom_indices.size()
		<< " fixed atoms " << std::endl;
      std::set<int>::const_iterator it;
      for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++)
	 std::cout << "    " << atom_spec_t(atom[*it]) << std::endl;
   }

   add_fixed_atoms_from_flanking_residues(bpc);
   add_fixed_atoms_from_non_bonded_neighbours();

   set_fixed_during_refinement_udd(); // uses fixed_atom_indices

   if (debug) {
      std::cout << "DEBUG:: init_from_residue_vec() Selecting residues gives "
		<< n_atoms << " atoms " << std::endl;
      for (int iat=0; iat<n_atoms; iat++) {
	 bool fixed_flag = false;
	 if (std::find(fixed_atom_indices.begin(),
		       fixed_atom_indices.end(), iat) != fixed_atom_indices.end())
	    fixed_flag = true;
	 std::cout << "   " << std::setw(3) << iat << " " << atom[iat]->name << " "
		   << atom[iat]->GetSeqNum() << " " << atom[iat]->GetChainID() << " "
		   << atom[iat]->GetResName() << " fixed: " << fixed_flag << std::endl;
      }
   }
}



void
coot::restraints_container_t::assign_fixed_atom_indices(const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   fixed_atom_indices.clear();

   for (unsigned int i=0; i<fixed_atom_specs.size(); i++) {
      for (int iat=0; iat<n_atoms; iat++) {
	 if (fixed_atom_specs[i].matches_spec(atom[iat])) {
	    fixed_atom_indices.insert(iat);
	 }
      }
   }
   //    std::cout << "Found indices for " << fixed_atom_indices.size()
   // << " fixed atoms" << std::endl;
}

// return false if any of the atoms are fixed
bool
coot::restraints_container_t::none_are_fixed_p(const std::vector<bool> &fixed_atom_flags) const {

   bool r = true;
   for (unsigned int i=0; i<fixed_atom_flags.size(); i++) {
      if (fixed_atom_flags[i]) {
	 r = false;
	 break;
      }
   }
   return r;
}



void
coot::restraints_container_t::debug_atoms() const {

   // use fixed_atom_indices to display if the atoms are fixed or not.

   std::cout << "---- " << n_atoms << " atoms" << std::endl;
   for (int iat=0; iat<n_atoms; iat++) {
      bool is_fixed = false;
      // use fixed_atom_indices.find()
      std::set<int>::const_iterator it;
      for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++) {
	 if (*it == iat) {
	    is_fixed = true;
	    break;
	 }
      }
      std::cout << iat << " " << atom_spec_t(atom[iat]) << "  "
		<< std::right << std::setw(10) << atom[iat]->x << " "
		<< std::right << std::setw(10) << atom[iat]->y << " "
		<< std::right << std::setw(10) << atom[iat]->z
		<< " fixed: " << is_fixed << std::endl;
   }
}

void
coot::restraints_container_t::debug_sets() const {

   std::cout << "-------------------- in debug_sets() residues_vec: " << std::endl;
   for (std::size_t i=0; i<residues_vec.size(); i++)
      std::cout << "   " << residues_vec[i].first << " " << residue_spec_t(residues_vec[i].second)
		<< std::endl;

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
   for (it=fixed_neighbours_set.begin(); it!=fixed_neighbours_set.end(); it++) {
      std::cout << "   Moving residue " << residue_spec_t(it->first) << std::endl;
      const std::set<mmdb::Residue *> &s = it->second;
      std::set<mmdb::Residue *>::const_iterator its;
      for (its=s.begin(); its!=s.end(); its++) {
	 std::cout << "     fixed neigb: " << residue_spec_t(*its) << std::endl;
      }
   }

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 std::cout << "   Chain " << chain_p->GetChainID() << std::endl;
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) {
	    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    std::cout << "      " << residue_spec_t(residue_p) << std::endl;
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    for (int iat=0; iat<n_atoms; iat++) {
	       mmdb::Atom *at = residue_p->GetAtom(iat);
	       std::cout << "   " << atom_spec_t(at) << std::endl;
	    }
	 }
      }
   }

}


void
coot::restraints_container_t::pre_sanitize_as_needed(std::vector<refinement_lights_info_t> lights) {


   bool do_pre_sanitize = false;
   for (std::size_t i=0; i<lights.size(); i++) {
      if (lights[i].worst_baddie.value > 100.0) {
	 do_pre_sanitize = true;
	 break;
      }
   }

   // do_pre_sanitize = false; // hack for debugging.
   int iter = 0;
   if (do_pre_sanitize) {
      if (verbose_geometry_reporting != QUIET)
	 std::cout << "debug:: :::: pre-sanitizing" << std::endl;
      int nsteps_max = 30; // wow! in the basic-no-progress test only 1 round is needed!
      int status;
      int restraints_usage_flag_save = restraints_usage_flag;
      restraints_usage_flag = BONDS_ANGLES_CHIRALS_AND_NON_BONDED;
      do
	 {
	    iter++;
	    status = gsl_multimin_fdfminimizer_iterate(m_s);
	    if (status) {
	       break;
	    }

	    double grad_lim = sqrt(size()) * 0.15;
	    if (grad_lim < 0.3)
	       grad_lim = 0.3;
	    status = gsl_multimin_test_gradient (m_s->gradient, grad_lim);

	    if (status == GSL_SUCCESS) {
	       if (verbose_geometry_reporting != QUIET) {
		  std::cout << "Pre-Sanitize Minimum found (iteration number " << iter << ") at ";
		  std::cout << m_s->f << "\n";
	       }
	    }
	 
	    if (status == GSL_ENOPROG)
	       std::cout << "pre-sanitize (No Progress)\n";
	 }
      while ((status == GSL_CONTINUE) && (iter < nsteps_max));
      gsl_vector_set_zero(m_s->dx);
      // set the 'starting' atom positions to be the result of pre-sanitization
      for (int i=0; i<3*n_atoms; i++)
	 gsl_vector_set(x, i, gsl_vector_get(m_s->x, i));
      restraints_usage_flag = restraints_usage_flag_save;
      gsl_multimin_fdfminimizer_set(m_s, &multimin_func, x, m_initial_step_size, m_tolerance);
      if (verbose_geometry_reporting != QUIET)
	 std::cout << "debug:: :::: pre-sanitization complete" << std::endl;
   }
}

unsigned int
coot::restraints_container_t::n_atom_pull_restraints() const {

   unsigned int n=0;
   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      const simple_restraint &rest = restraints_vec[i];
      if (rest.restraint_type == TARGET_POS_RESTRAINT)
	 n++;
   }

   return n;
}


// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
// 
coot::refinement_results_t
coot::restraints_container_t::minimize(restraint_usage_Flags usage_flags) {

   short int print_chi_sq_flag = 1;
   refinement_results_t rr = minimize(usage_flags, 1000, print_chi_sq_flag);
   // std::cout << "minimize() returns " << rr.progress << std::endl;
   return rr;

}


#include <gsl/gsl_blas.h> // for debugging norm of gradient

// this does a reset/setup
//
void
coot::restraints_container_t::setup_minimize() {

   if (m_s)
      gsl_multimin_fdfminimizer_free(m_s);
   if (x)
      gsl_vector_free(x);

   // T = gsl_multimin_fdfminimizer_conjugate_fr; // not as good as pr
   // T = gsl_multimin_fdfminimizer_steepest_descent; // pathetic
   // T = gsl_multimin_fminimizer_nmsimplex; // you can't just drop this in,
                                             // because simplex is a
                                             // non-gradient method.
   // T = gsl_multimin_fdfminimizer_vector_bfgs;

   //   This is the Polak-Ribiere conjugate gradient algorithm. It is
   //   similar to the Fletcher-Reeves method, differing only in the
   //   choice of the coefficient \beta. Both methods work well when
   //   the evaluation point is close enough to the minimum of the
   //   objective function that it is well approximated by a quadratic
   //   hypersurface.
   //
   const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;

   // setup before the first minimize (n_times_called == 1)

   setup_gsl_vector_variables();  //initial positions
   setup_multimin_func(); // provide functions for f, df, fdf

   m_s = gsl_multimin_fdfminimizer_alloc(T, n_variables());

   m_initial_step_size = 0.5 * gsl_blas_dnrm2(x); // how about just 0.1?

   // std::cout << "debug:: setup_minimize() with step_size " << m_initial_step_size << std::endl;

   // this does a lot of work
   gsl_multimin_fdfminimizer_set(m_s, &multimin_func, x, m_initial_step_size, m_tolerance);

   m_grad_lim = pow(size(), 0.7) * 0.03;
   if (m_grad_lim < 0.3) m_grad_lim = 0.3;

   needs_reset = false;

}
 
// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
//
coot::refinement_results_t
coot::restraints_container_t::minimize(restraint_usage_Flags usage_flags,
				       int nsteps_max,
				       short int print_initial_chi_sq_flag) {

   std::cout << "------------ ::minimize() basic " << restraints_vec.size() << std::endl;
   n_times_called++;
   if (n_times_called == 1 || needs_reset)
      setup_minimize();

   refinement_results_t rr = minimize_inner(usage_flags, nsteps_max, print_initial_chi_sq_flag);

   return rr;
}

// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
//
coot::refinement_results_t
coot::restraints_container_t::minimize(int imol, restraint_usage_Flags usage_flags, int nsteps_max,
				       short int print_initial_chi_sq_flag,
				       const coot::protein_geometry &geom) {

   // std::cout << "------------ ::minimize() " << restraints_vec.size() << std::endl;

   unsigned int n_steps_per_relcalc_nbcs = 3000000;

   n_times_called++;
   n_small_cycles_accumulator += n_times_called * nsteps_max;

   if (n_times_called == 1 || needs_reset)
      setup_minimize();

#if 0

   // This doesn't work - and on reflection as to why that is, it *cannot* work.
   // make_non_bonded_contact_restraints_ng() works using the residues of the
   // input moolecule - and they are not going to change as the atoms in the atom
   // vector move around.  We need a function that takes the atoms in atom and
   // finds neighbours in mol - which is extremely not how refinement works in Coot,
   // because we find the residues to refine before refinement.
   //
   // It will take quite a rewrite to fix this - pass the moving residues and mol
   // and let a function in this class work out what the "residues near residues"
   // are. This is quite doable, but not for now - updating NBCs is not as
   // important as OpenGLv3 (shader) graphics.

   if (n_small_cycles_accumulator >= n_steps_per_relcalc_nbcs) {
      auto tp_0 = std::chrono::high_resolution_clock::now();
      make_non_bonded_contact_restraints_ng(imol, geom);
      auto tp_1 = std::chrono::high_resolution_clock::now();
      setup_minimize();
      auto tp_2 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
      std::cout << "minimize() nbc updates " << d10 << " " << d21 << " milliseconds" << std::endl;
      n_small_cycles_accumulator = 0;
   }
#endif

   refinement_results_t rr = minimize_inner(usage_flags, nsteps_max, print_initial_chi_sq_flag);

   return rr;
}

// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
//
coot::refinement_results_t
coot::restraints_container_t::minimize_inner(restraint_usage_Flags usage_flags, 
					     int nsteps_max,
					     short int print_initial_chi_sq_flag) {

   // Was just checking that minimize() was not being called multiple times concurrently
   // (it wasn't)
   // std::cout << "DEBUG:: incrementing n_refiners_refining to " << n_refiners_refining+1 << std::endl;
   // n_refiners_refining++;

   // check that we have restraints before we start to minimize:
   if (restraints_vec.size() == 0) {
      if (restraints_usage_flag != NO_GEOMETRY_RESTRAINTS) {
	 std::cout << "SPECIFICATION ERROR:  There are no restraints. ";
	 std::cout << "No minimization will happen" << std::endl;
	 return coot::refinement_results_t(0, 0, "No Restraints!");
      }
   }

   restraints_usage_flag = usage_flags;

   if (do_numerical_gradients_flag) {
      std::cout << "debug:: minimize_inner called with usage_flags " << usage_flags << std::endl;
      debug_atoms();
   }

   // BOND + density fail: BONDS regularize works
   // restraints_usage_flag = BONDS_AND_ANGLES;

   // We get ~1ms/residue with bond and angle terms and no density terms.

   std::vector<refinement_lights_info_t> lights; // = chi_squareds("--------", m_s->x, false);

   int iter = 0; 
   int status = GSL_SUCCESS; // some start value for the compiler not to complain about
                             // construction of a refinement_results_t
   std::vector<coot::refinement_lights_info_t> lights_vec;
   bool done_final_chi_squares = false;

   // Note to self, after an atom pull, should I be using
   // gsl_multimin_fdfminimizer_restart() to set to the current position as a starting point?

   typedef struct
   {
      int iter;
      double step;
      double max_step;
      double tol;
      gsl_vector *x1;
      gsl_vector *dx1;
      gsl_vector *x2;
      double pnorm;
      gsl_vector *p;
      double g0norm;
      gsl_vector *g0;
   }
   conjugate_pr_state_t;

   do
      {
	 iter++;

	 // we should not update the atom pull restraints while the refinement is running.
	 // we shouldn't refine when the atom pull restraints are being updated.

#ifdef HAVE_CXX_THREAD
	 bool unlocked = false;
	 while (! restraints_lock.compare_exchange_weak(unlocked, true) && !unlocked) {
	    std::this_thread::sleep_for(std::chrono::microseconds(10));
	    unlocked = false;
	 }
#endif

         if (m_s == 0) {
            std::cout << "ERROR:: !! m_s has disappeared! " << std::endl;
            break;
         } else {
	    status = gsl_multimin_fdfminimizer_iterate(m_s);
         }

#ifdef HAVE_CXX_THREAD
	 restraints_lock = false; // unlock
#endif
	 // this might be useful for debugging rama restraints

	 conjugate_pr_state_t *state = (conjugate_pr_state_t *) m_s->state;
	 double pnorm = state->pnorm;
	 double g0norm = state->g0norm;
	 //
	 if (false)
	    std::cout << "iter: " << iter << " f " << m_s->f << " " << gsl_multimin_fdfminimizer_minimum(m_s)
		      << " pnorm " << pnorm << " g0norm " << g0norm << std::endl;

	 if (status != GSL_SUCCESS) {
	    std::cout << "Unexpected error from gsl_multimin_fdfminimizer_iterate at iter " << iter << std::endl;
	    if (status == GSL_ENOPROG) {
	       std::cout << "Error:: in gsl_multimin_fdfminimizer_iterate() result was GSL_ENOPROG" << std::endl; 
	       if (true)
		  std::cout << "Error:: iter: " << iter << " f " << m_s->f << " "
			    << gsl_multimin_fdfminimizer_minimum(m_s)
			    << " pnorm " << pnorm << " g0norm " << g0norm << "\n";

	       // write out gradients here - with numerical gradients for comparison

	       lights_vec = chi_squareds("Final Estimated RMS Z Scores (ENOPROG)", m_s->x);
	       done_final_chi_squares = true;
	       refinement_lights_info_t::the_worst_t worst_of_all = find_the_worst(lights_vec);
	       if (worst_of_all.is_set) {
		  const simple_restraint &baddie_restraint = restraints_vec[worst_of_all.restraints_index];
		  std::cout << "INFO:: Most dissatisfied restraint (refine no-progress): "
			    << baddie_restraint.format(atom, worst_of_all.value) << std::endl;
	       } else {
		  std::cout << "INFO:: somehow the worst restraint was not set (no-progress)"
			    << std::endl;
	       }
	    }
	    if (status == GSL_ENOPROG) {
               // debugging/analysis
               std::cout << "----------------------- FAIL, ENOPROG --------------- " << std::endl;
               gsl_vector *non_const_v = const_cast<gsl_vector *> (m_s->x); // because there we use gls_vector_set()
               void *params = static_cast<void *>(this);
               // useful - but not for everyone
               // numerical_gradients(non_const_v, params, m_s->gradient, "failed-gradients.tab");
            }
	    break;
	 }

         // std::cout << "Debug:: before gsl_multimin_test_gradient, status is " << status << std::endl;

         if (status == GSL_SUCCESS || status == GSL_CONTINUE) // probably just GSL_SUCCESS is what I want
            status = gsl_multimin_test_gradient(m_s->gradient, m_grad_lim);

         // std::cout << "Debug:: after gsl_multimin_test_gradient, status is " << status << std::endl;

	 if (status == GSL_SUCCESS) {
	    if (verbose_geometry_reporting != QUIET) { 
	       std::cout << "Minimum found (iteration number " << iter << ") at ";
	       std::cout << m_s->f << "\n";
	    }
	 }

	 if (status == GSL_SUCCESS) {
	    std::string title = "Final Estimated RMS Z Scores:";
	    std::vector<coot::refinement_lights_info_t> results = chi_squareds(title, m_s->x);
	    lights_vec = results;
	    done_final_chi_squares = true;
	 }

	 if (verbose_geometry_reporting == VERBOSE)
	    std::cout << "iteration number " << iter << " " << m_s->f << std::endl;

      }
   while ((status == GSL_CONTINUE) && (iter < nsteps_max));

   // std::cout << "Debug:: post loop status is " << status << std::endl;

   if (! done_final_chi_squares) {
      if (status != GSL_CONTINUE) {
	 lights = chi_squareds("Final Estimated RMS Z Scores:", m_s->x);
	 refinement_lights_info_t::the_worst_t worst_of_all = find_the_worst(lights);
	 if (worst_of_all.is_set) {
	    const simple_restraint &baddie_restraint = restraints_vec[worst_of_all.restraints_index];
	    std::cout << "Most dissatisfied restraint (post-refine): "
		      << baddie_restraint.format(atom, worst_of_all.value) << std::endl;
	 }
      }
   }

   if (false) // this is probably useful in some cases.
      if (iter == nsteps_max)
      std::cout << "Hit nsteps_max " << nsteps_max << " " << m_s->f << std::endl;

   update_atoms(m_s->x); // do OXT here

   // (we don't get here unless restraints were found)
   coot::refinement_results_t rr(1, status, lights_vec);

   if (status != GSL_CONTINUE) {

#ifdef HAVE_CXX_THREAD

      // protection so that clearing of the vectors doesn't coincide with geometric_distortions()
      // evaluation

      bool unlocked = false;
      while (! restraints_lock.compare_exchange_weak(unlocked, true) && !unlocked) {
	 std::this_thread::sleep_for(std::chrono::microseconds(10));
	 unlocked = false;
	 }
#endif

      std::cout << "DEBUG:: ---- free/delete/reset m_s and x" << std::endl; // works fine
      gsl_multimin_fdfminimizer_free(m_s);
      gsl_vector_free(x);
      m_s = 0;
      x = 0;
      needs_reset = true;
#ifdef HAVE_CXX_THREAD
      restraints_lock = false; // unlock
#endif
   }

   // the bottom line from the timing test is the only thing that matters
   // is the time spend in the core minimization iterations

   // std::cout << "-------------- Finally returning from minimize_inner() with rr with status "
   //           << rr.progress << std::endl;
   n_refiners_refining--;
   return rr;
}

// this should be a static member function of refinement_lights_info_t
//
coot::refinement_lights_info_t::the_worst_t
coot::restraints_container_t::find_the_worst(const std::vector<refinement_lights_info_t> &lights) const {

   bool debug = false; // well, it's "give me extra baddie info" rather than "debug"
   refinement_lights_info_t::the_worst_t worst_of_all;
   for (std::size_t ii=0; ii<lights.size(); ii++) {
      if (lights[ii].worst_baddie.is_set)
	 worst_of_all.update_if_worse(lights[ii].worst_baddie);
   }

   if (debug) {
      for (std::size_t ii=0; ii<lights.size(); ii++) {
	 std::cout << "Refinements Lights "
		   << lights[ii].name  << " "
		   << lights[ii].label << " "
		   << lights[ii].value;
	 if (lights[ii].worst_baddie.is_set) {
	    std::cout << " worst baddie "
		      << lights[ii].worst_baddie.value << " index: "
		      << lights[ii].worst_baddie.restraints_index;
	 }
	 std::cout << "\n";
      }
   }
   return worst_of_all;
}


std::ostream &
coot::operator<<(std::ostream &s, const simple_restraint &r) {

   s << "{restraint: ";
   if (r.restraint_type == coot::BOND_RESTRAINT)
      s << "Bond   ";
   if (r.restraint_type == coot::TARGET_POS_RESTRAINT)
      s << "Target-Pos ";
   if (r.restraint_type == coot::ANGLE_RESTRAINT)
      s << "Angle  ";
   if (r.restraint_type == coot::TORSION_RESTRAINT)
      s << "Torsion";
   if (r.restraint_type == coot::PLANE_RESTRAINT)
      s << "Plane  ";
   if (r.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT)
      s << "NBC    ";
   if (r.restraint_type == coot::TRANS_PEPTIDE_RESTRAINT)
      s << "Trans-Pep ";
   if (r.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
      s << "Chiral ";
      s << r.atom_index_centre;
   }
   if (r.restraint_type == coot::RAMACHANDRAN_RESTRAINT)
      s << "Rama   ";
   s << "}";
   return s;
}

// this is more useful:
std::string
coot::simple_restraint::format(mmdb::PAtom *atoms_vec, double distortion) const {

   std::string s;
   if (restraint_type == BOND_RESTRAINT) {
      s = "Bond ";
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      atom_spec_t spec_2(atoms_vec[atom_index_2]);
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      if (distortion >= 0) {
	 s += "  z = ";
	 s += util::float_to_string_using_dec_pl(sqrt(distortion), 2);
      }
   }
   if (restraint_type == ANGLE_RESTRAINT) {
      s = "Angle ";
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      atom_spec_t spec_2(atoms_vec[atom_index_2]);
      atom_spec_t spec_3(atoms_vec[atom_index_3]);
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      if (distortion >= 0) {
	 s += "  z = ";
	 s += util::remove_whitespace(util::float_to_string_using_dec_pl(sqrt(distortion), 2));
      }
   }
   if (restraint_type == TORSION_RESTRAINT) {
      s = "Torsion ";
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      atom_spec_t spec_2(atoms_vec[atom_index_2]);
      atom_spec_t spec_3(atoms_vec[atom_index_3]);
      atom_spec_t spec_4(atoms_vec[atom_index_4]);
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      s += spec_4.label();
      s += " ";
      if (distortion >= 0) {
	 s += "  z = ";
	 s += util::remove_whitespace(util::float_to_string_using_dec_pl(sqrt(distortion), 2));
      }
   }
   if (restraint_type == TRANS_PEPTIDE_RESTRAINT) {
      s = "Trans_peptide ";
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      atom_spec_t spec_2(atoms_vec[atom_index_2]);
      atom_spec_t spec_3(atoms_vec[atom_index_3]);
      atom_spec_t spec_4(atoms_vec[atom_index_4]);
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      s += spec_4.label();
      s += " ";
      if (distortion >= 0) {
	 s += "  score: ";
	 s += util::remove_whitespace(util::float_to_string_using_dec_pl(distortion, 2));
	 s += " (non-sqrt)";
      }
   }
   if (restraint_type == PLANE_RESTRAINT) {
      s = "Plane ";
      for (std::size_t j=0; j<plane_atom_index.size(); j++) {
	 atom_spec_t spec(atoms_vec[plane_atom_index[j].first]);
	 s += spec.label();
	 s += " ";
      }
      if (distortion >= 0) {
	 s += "  z = ";
	 s += util::remove_whitespace(util::float_to_string_using_dec_pl(sqrt(distortion), 2));
      }
   }
   if (restraint_type == NON_BONDED_CONTACT_RESTRAINT) {
      s = "Non-Bonded-Contact ";
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      atom_spec_t spec_2(atoms_vec[atom_index_2]);
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      if (distortion >= 0) {
	 s += "  z = ";
	 s += util::remove_whitespace(util::float_to_string_using_dec_pl(sqrt(distortion), 2));
      }
   }
   if (restraint_type == CHIRAL_VOLUME_RESTRAINT) {
      s = "Chiral ";
      atom_spec_t spec_c(atoms_vec[atom_index_centre]);
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      atom_spec_t spec_2(atoms_vec[atom_index_2]);
      atom_spec_t spec_3(atoms_vec[atom_index_3]);
      s += spec_c.label();
      s += " neigbhs: ";
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      if (distortion >= 0) {
	 s += "  z = ";
	 s += util::remove_whitespace(util::float_to_string_using_dec_pl(sqrt(distortion), 2));
      }
   }
   if (restraint_type == RAMACHANDRAN_RESTRAINT) {
      s = "Rama ";
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      atom_spec_t spec_2(atoms_vec[atom_index_2]);
      atom_spec_t spec_3(atoms_vec[atom_index_3]);
      atom_spec_t spec_4(atoms_vec[atom_index_4]);
      atom_spec_t spec_5(atoms_vec[atom_index_5]);
      s += spec_1.label();
      s += " ";
      s += spec_2.label();
      s += " ";
      s += spec_3.label();
      s += " ";
      s += spec_4.label();
      s += " ";
      s += spec_5.label();
      s += " ";
      s += util::remove_whitespace(util::float_to_string_using_dec_pl(distortion, 2));
   }

   if (restraint_type == TARGET_POS_RESTRAINT) {
      s = "Target_pos ";
      atom_spec_t spec_1(atoms_vec[atom_index_1]);
      s += spec_1.label();
      s += " ";
      s += util::remove_whitespace(util::float_to_string_using_dec_pl(distortion, 2));
   }

   return s;
}

std::string
coot::simple_restraint::type() const {

   std::string s;
   if (restraint_type == coot::BOND_RESTRAINT)
      s = "Bond";
   if (restraint_type == coot::ANGLE_RESTRAINT)
      s = "Angle";
   if (restraint_type == coot::TORSION_RESTRAINT)
      s = "Torsion";
   if (restraint_type == coot::PLANE_RESTRAINT)
      s = "Plane";
   if (restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT)
      s = "NBC";
   if (restraint_type == coot::CHIRAL_VOLUME_RESTRAINT)
      s = "Chiral";
   if (restraint_type == coot::RAMACHANDRAN_RESTRAINT)
      s = "Rama";
   return s;
   
}


void
coot::restraints_container_t::adjust_variables(const atom_selection_container_t &asc) { 

}

// 
double
starting_structure_diff_score(const gsl_vector *v, void *params) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   double d;
   double dist = 0; 

   // if (v->size != restraints->init_positions_size() ) {

   for (int i=0; i<restraints->init_positions_size(); i++) { 
      d = restraints->initial_position(i) - gsl_vector_get(v, i);
      dist += 0.01*d*d;
   }
   std::cout << "starting_structure_diff_score: " << dist << std::endl; 
   return dist;
}


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
   double non_bonded_distortion = 0;
   double chiral_vol_distortion = 0;
   double rama_distortion = 0;
   double start_pos_distortion = 0;
   double target_pos_distortion = 0;
   double trans_peptide_distortion = 0;

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

   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      if (restraints_usage_flag & BONDS_MASK) {
	 if (restraints_vec[i].restraint_type == BOND_RESTRAINT) {
	    n_bond_restraints++;
	    // 	    bond_distortion += distortion_score_bond(restraints_vec[i], v);
	    double dist = distortion_score_bond(restraints_vec[i], v);
	    const simple_restraint &rest = restraints_vec[i];
	    bond_distortion += dist;
	    if (dist > dist_max_bonds.second) {
	       dist_max_bonds.first = i;
	       dist_max_bonds.second = dist;
	    }
	    baddies["Bonds"].update_if_worse(dist, i);
	 }
      }

      if (restraints_usage_flag & GEMAN_MCCLURE_DISTANCE_MASK) {
	 if (restraints_vec[i].restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    n_geman_mcclure_distance++;
	    double d = distortion_score_geman_mcclure_distance(restraints_vec[i], v, geman_mcclure_alpha);
	    gm_distortion += d;
	    baddies["GemanMcClure"].update_if_worse(d, i);
	 }
      }

      if (restraints_usage_flag & ANGLES_MASK) { // 2: angles
	 if (restraints_vec[i].restraint_type == coot::ANGLE_RESTRAINT) {
	    n_angle_restraints++;
	    double dist = coot::distortion_score_angle(restraints_vec[i], v);
	    angle_distortion += dist;
	    if (dist > dist_max_angles.second) {
	       dist_max_angles.first = i;
	       dist_max_angles.second = dist;
	    }
	    baddies["Angles"].update_if_worse(dist, i);
	 }
      }

      if (restraints_usage_flag & TORSIONS_MASK) { // 4: torsions
	 if (restraints_vec[i].restraint_type == coot::TORSION_RESTRAINT) {
	    try {
	       double dist = coot::distortion_score_torsion(i, restraints_vec[i], v);
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
	 if (restraints_vec[i].restraint_type == TRANS_PEPTIDE_RESTRAINT) {
	    try {
	       double dist = distortion_score_trans_peptide(i, restraints_vec[i], v);
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
	 if (restraints_vec[i].restraint_type == coot::PLANE_RESTRAINT) {
	    n_plane_restraints++;
	    double dist = coot::distortion_score_plane(restraints_vec[i], v);
	    plane_distortion += dist;
	    if (dist > dist_max_planes.second) {
	       dist_max_planes.first = i;
	       dist_max_planes.second = dist;
	    }
	    baddies["Planes"].update_if_worse(dist, i);
            if (false) {  // debugging plane restraints.
                std::cout << " plane distortion " << i << " " 
                          << coot::distortion_score_plane(restraints_vec[i], v) << " " 
                          << restraints_vec[i];
                for (unsigned int jj = 0; jj<restraints_vec[i].plane_atom_index.size(); jj+=3) { 
                    std::cout << "\n                                ";
		    unsigned int idx = restraints_vec[i].plane_atom_index[jj].first;
                    std::cout << idx << " " << coot::atom_spec_t(atom[idx]);
                    if ((jj+1) < restraints_vec[i].plane_atom_index.size()) { 
		       unsigned int idx_1 = restraints_vec[i].plane_atom_index[jj+1].first;
                       std::cout << " " << idx_1 << " " << coot::atom_spec_t(atom[idx_1]);
                    }
                    if ((jj+2) < restraints_vec[i].plane_atom_index.size()) { 
		       unsigned int idx_2 = restraints_vec[i].plane_atom_index[jj+1].first;
                       std::cout << " " << idx_2 << " " << coot::atom_spec_t(atom[idx_2]);
                    }
                }
                std::cout << std::endl;
            }
	 }
      }

      if (restraints_usage_flag & coot::NON_BONDED_MASK) { 
	 if ( restraints_vec[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) { 
	    n_non_bonded_restraints++;
	    double dist = coot::distortion_score_non_bonded_contact(restraints_vec[i], lennard_jones_epsilon, v);
	    non_bonded_distortion += dist;
	    if (dist > dist_max_nbc.second) {
	       dist_max_nbc.first = i;
	       dist_max_nbc.second = dist;
	    }
	    baddies["NonBonded"].update_if_worse(dist, i);
	 }
      }

      if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) { 
  	 if ( restraints_vec[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
  	    n_chiral_volumes++;
	    double dist = coot::distortion_score_chiral_volume(restraints_vec[i], v);
  	    chiral_vol_distortion += dist;
	    baddies["Chirals"].update_if_worse(dist, i);
  	 }
      }

      if (restraints_usage_flag & coot::RAMA_PLOT_MASK) {
  	 if (restraints_vec[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) {
  	    n_rama_restraints++;
	    if (rama_type == restraints_container_t::RAMA_TYPE_ZO) {
	       double dd = distortion_score_rama( restraints_vec[i], v, ZO_Rama(), get_rama_plot_weight());
	       // std::cout << "Here with index " << i << " rama distortion score " << dd << std::endl;
	       rama_distortion += dd;
	       baddies["Rama"].update_if_worse(dd, i);

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

	    } else {
	       double dd = distortion_score_rama(restraints_vec[i], v, lograma);
	       rama_distortion += dd;
	       baddies["Rama"].update_if_worse(dd, i);
	    }
	    if (false) {
	       double d1 = distortion_score_rama(restraints_vec[i], v, LogRama());
	       double d2 = coot::distortion_score_rama(restraints_vec[i], v, ZO_Rama(), get_rama_plot_weight());
	       std::cout << "distortion-comparision logramas " << d1 << " zo " << d2 << std::endl;
	    }
	 }
      }

      if (restraints_vec[i].restraint_type == coot::TARGET_POS_RESTRAINT) {
	 n_target_pos_restraints++;
	 double dist = coot::distortion_score_target_pos(restraints_vec[i],
							 log_cosh_target_distance_scale_factor, v);
         target_pos_distortion += dist;
	 baddies["Target_pos"].update_if_worse(dist, i);
      }

      if (restraints_vec[i].restraint_type == coot::START_POS_RESTRAINT) {
         n_start_pos_restraints++;
	 double dist = distortion_score_start_pos(restraints_vec[i], params, v);
         start_pos_distortion += dist;
	 baddies["StartPositions"].update_if_worse(dist, i);
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
      double pd = plane_distortion/double(n_plane_restraints);
      double spd = 0.0;
      if (pd > 0.0)
	 spd = sqrt(pd);
      if (print_summary)
	 std::cout << "planes:     " << spd << std::endl;
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


// public
coot::model_bond_deltas
coot::restraints_container_t::resolve_bonds() {

   setup_gsl_vector_variables();
   return resolve_bonds(x);

}

coot::model_bond_deltas
coot::restraints_container_t::resolve_bonds(const gsl_vector *v) const {

   model_bond_deltas resultant;

   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & BONDS_MASK) {
	 const simple_restraint &rest = restraints_vec[i];
	 if (rest.restraint_type == BOND_RESTRAINT) {

	    int idx = 3*(rest.atom_index_1);
	    clipper::Coord_orth a1(gsl_vector_get(v,idx),
				   gsl_vector_get(v,idx+1),
				   gsl_vector_get(v,idx+2));
	    idx = 3*(rest.atom_index_2);
	    clipper::Coord_orth a2(gsl_vector_get(v,idx),
				   gsl_vector_get(v,idx+1),
				   gsl_vector_get(v,idx+2));
	    double ideal = rest.target_value;
	    double bl = clipper::Coord_orth::length(a1,a2);
	    double delta = bl - ideal;
	    clipper::Coord_orth b_uv((a2-a1).unit());
	    clipper::Coord_orth b_uv_abs(std::fabs(b_uv.x()),
					 std::fabs(b_uv.y()),
					 std::fabs(b_uv.z()));
	    clipper::Coord_orth frag(b_uv_abs * delta);
	    resultant.add(clipper::Coord_orth(b_uv_abs * delta));
	 }
      }
   }
   return resultant;
}


// Ah, but (c.f. distortion) we want to return a low value for a good
// fit and a high one for a bad.
double
coot::electron_density_score(const gsl_vector *v, void *params) {

   coot::restraints_container_t *restraints_p = static_cast<restraints_container_t *> (params);
   if (restraints_p->include_map_terms() == 1) {
      return electron_density_score_from_restraints(v, restraints_p);
   } else {
      return 0.0;
   }
}


// Ah, but (c.f. distortion) we want to return a low value for a good
// fit and a high one for a bad.
double
coot::electron_density_score_from_restraints_simple(const gsl_vector *v,
					     coot::restraints_container_t *restraints_p) {

#ifdef HAVE_CXX_THREAD
   auto tp_1 = std::chrono::high_resolution_clock::now();
#endif
   // We weight and sum to get the score and negate.
   // 
   double score = 0;

   if (restraints_p->include_map_terms() == 1) { 

      unsigned int n_atoms = restraints_p->get_n_atoms();
      for (unsigned int iat=0; iat<n_atoms; iat++) {
	 bool use_it = false;
	 if (restraints_p->use_map_gradient_for_atom[iat]) {

	    int idx = 3 * iat;
	    clipper::Coord_orth ao(gsl_vector_get(v,idx),
				   gsl_vector_get(v,idx+1),
				   gsl_vector_get(v,idx+2));

	    score += restraints_p->Map_weight() *
	       restraints_p->atom_z_occ_weight[iat] *
	       restraints_p->electron_density_score_at_point(ao);
	 }
      }
   }
#ifdef HAVE_CXX_THREAD
   auto tp_2 = std::chrono::high_resolution_clock::now();
   auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
   // std::cout << "info:: f electron_density: " << d21 << " microseconds\n";
#endif // HAVE_CXX_THREAD

   return -score;
}


double
coot::electron_density_score_from_restraints(const gsl_vector *v,
					     coot::restraints_container_t *restraints_p) {

   double score = 0.0;
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

   std::vector<std::pair<unsigned int, unsigned int> > &ranges = restraints_p->m_atom_index_ranges;

   // they are no longer "threads" - they are restraints sets - to prevent
   // "thread starvation"
   std::atomic<unsigned int> done_count_for_restraints_sets(0);

   double results[1024]; // we will always have less than 1024 threads

   if (restraints_p->thread_pool_p) {
      for(unsigned int i=0; i<ranges.size(); i++) {
         results[i] = 0.0;
	 restraints_p->thread_pool_p->push(electron_density_score_from_restraints_using_atom_index_range,
					   v, std::cref(ranges[i]), restraints_p, &results[i],
					   std::ref(done_count_for_restraints_sets));
      }
      while (done_count_for_restraints_sets < ranges.size()) {
	 std::this_thread::sleep_for(std::chrono::microseconds(1));
      }

      // consolidate
      for(unsigned int i=0; i<ranges.size(); i++)
	 score += results[i];
   } else {
      // null thread pool. restraints_container_t was created without a call to
      // set the thread_pool. Happens in crankshafting currently.
      electron_density_score_from_restraints_using_atom_index_range(0, v, ranges[0], restraints_p, &score,
								    done_count_for_restraints_sets);
   }

#else
   std::cout << __FUNCTION__ << " no thread pool" << std::endl;
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   return score;
}

#ifdef HAVE_CXX_THREAD
// atom_index_range works "as expected"
// so given atom_index_range of 0,10 we start at the first value (0) and check that the
// current value is less than the atom_index_range.second (10):
// ie. density values for atom indices 0 to 9 inclusive are added.
//
void
coot::electron_density_score_from_restraints_using_atom_index_range(int thread_index,
                                             const gsl_vector *v,
					     const std::pair<unsigned int, unsigned int> &atom_index_range,
					     coot::restraints_container_t *restraints_p,
					     double *result,
					     std::atomic<unsigned int> &done_count_for_threads) {

   auto tp_1 = std::chrono::high_resolution_clock::now();

   // We weight and sum to get the score and negate.
   //
   double score = 0;

   if (restraints_p->include_map_terms() == 1) {

      if (false)
	 std::cout << "debug:: " << __FUNCTION__ << "() range: "
		   << atom_index_range.first << " " << atom_index_range.second << std::endl;

      for (unsigned int iat=atom_index_range.first; iat<atom_index_range.second; iat++) {

	 // we get a crash around here. Protect from wrong index into v vector
	 if (iat < restraints_p->get_n_atoms()) {
	    if (restraints_p->use_map_gradient_for_atom[iat]) {

	       int idx = 3 * iat;
	       clipper::Coord_orth ao(gsl_vector_get(v,idx),
				      gsl_vector_get(v,idx+1),
				      gsl_vector_get(v,idx+2));

	       score += restraints_p->Map_weight() *
		  restraints_p->atom_z_occ_weight[iat] *
		  restraints_p->electron_density_score_at_point(ao);
	    }
	 } else {
	    std::cout << "ERROR:: electron_density_score_from_restraints_using_atom_index_range "
		      << " caught bad atom index " << iat << " " << restraints_p->get_n_atoms()
		      << std::endl;
	 }
      }
   }

   auto tp_2 = std::chrono::high_resolution_clock::now();
   auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
   // std::cout << "info:: f electron_density: " << d21 << " microseconds\n";
   // return -score;

   *result = -score;
   done_count_for_threads++; // atomic
}
#endif




// Note that the gradient for the electron density is opposite to that
// of the gradient for the geometry (consider a short bond on the edge
// of a peak - in that case the geometry gradient will be negative as
// the bond is lengthened and the electron density gradient will be
// positive).
//
// So we want to change that positive gradient for a low score when
// the atoms coinside with the density - hence the contributions that
// we add are negated.
//
void coot::my_df_electron_density(const gsl_vector *v,
				  void *params,
				  gsl_vector *df) {

   // first extract the object from params
   //
   coot::restraints_container_t *restraints_p = static_cast<restraints_container_t *> (params);
   if (restraints_p->include_map_terms() == 1) {
      my_df_electron_density_single(v, restraints_p, df, 0, v->size/3);
   }
}



#ifdef HAVE_CXX_THREAD


// restraints are modified by atomic done_count_for_threads changing.
//
// multi-threaded inner
//
void coot::my_df_electron_density_threaded_single(int thread_idx, const gsl_vector *v,
						  coot::restraints_container_t *restraints,
						  gsl_vector *df,
						  int atom_idx_start, int atom_idx_end,
						  std::atomic<unsigned int> &done_count_for_threads) {

   for (int iat=atom_idx_start; iat<atom_idx_end; ++iat) {
      if (restraints->use_map_gradient_for_atom[iat]) {

	 int idx = 3 * iat;
	 clipper::Coord_orth ao(gsl_vector_get(v,idx), 
				gsl_vector_get(v,idx+1), 
				gsl_vector_get(v,idx+2));
	    
	 clipper::Grad_orth<double> grad_orth = restraints->electron_density_gradient_at_point(ao);
	 float zs = restraints->Map_weight() * restraints->atom_z_occ_weight[iat];

	 if (0) { 
	    std::cout << "electron density df: adding "
		      <<  - zs * grad_orth.dx() << " "
		      <<  - zs * grad_orth.dy() << " "
		      <<  - zs * grad_orth.dz() << " to "
		      <<  gsl_vector_get(df, idx  ) << " "
		      <<  gsl_vector_get(df, idx+1) << " "
		      <<  gsl_vector_get(df, idx+2) << "\n";
	 }
	    
	 // 	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
	 // 	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
	 // 	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());


         // we no longer do this sort of locking

	 // use atomic lock to access derivs of atom idx
	 // unsigned int unlocked = 0;
	 // while (! restraints->gsl_vector_atom_pos_deriv_locks.get()[idx].compare_exchange_weak(unlocked, 1)) {
	 //    std::this_thread::sleep_for(std::chrono::nanoseconds(10));
	 //    unlocked = 0;
	 // }

	 *gsl_vector_ptr(df, idx  ) -= zs * grad_orth.dx();
	 *gsl_vector_ptr(df, idx+1) -= zs * grad_orth.dy();
	 *gsl_vector_ptr(df, idx+2) -= zs * grad_orth.dz();
	 // restraints->gsl_vector_atom_pos_deriv_locks.get()[idx] = 0; // unlock
      }
   }

   // int sleepy_time = 10000 + int(2000*float(coot::util::random())/float(RAND_MAX));
   // std::this_thread::sleep_for(std::chrono::microseconds(sleepy_time));
   ++done_count_for_threads; // atomic
}
#endif	 


//
void coot::my_df_electron_density_single(const gsl_vector *v,
					 coot::restraints_container_t *restraints,
					 gsl_vector *df,
					 int atom_idx_start, int atom_idx_end) {

   //    std::cout << "debug:: in my_df_electron_density_single() " << atom_idx_start << " " << atom_idx_end
   // << std::endl;

   for (int iat=atom_idx_start; iat<atom_idx_end; ++iat) {
      if (restraints->use_map_gradient_for_atom[iat]) {

	 int idx = 3 * iat;
	 clipper::Coord_orth ao(gsl_vector_get(v,idx), 
				gsl_vector_get(v,idx+1), 
				gsl_vector_get(v,idx+2));
	    
	 clipper::Grad_orth<double> grad_orth = restraints->electron_density_gradient_at_point(ao);
	 float zs = restraints->Map_weight() * restraints->atom_z_occ_weight[iat];

	 if (0) { 
	    std::cout << "electron density df: adding "
		      <<  - zs * grad_orth.dx() << " "
		      <<  - zs * grad_orth.dy() << " "
		      <<  - zs * grad_orth.dz() << " to "
		      <<  gsl_vector_get(df, idx  ) << " "
		      <<  gsl_vector_get(df, idx+1) << " "
		      <<  gsl_vector_get(df, idx+2) << "\n";
	 }

	 // 	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
	 // 	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
	 // 	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());

	 *gsl_vector_ptr(df, idx  ) -= zs * grad_orth.dx();
	 *gsl_vector_ptr(df, idx+1) -= zs * grad_orth.dy();
	 *gsl_vector_ptr(df, idx+2) -= zs * grad_orth.dz();
      }
   }
}

void coot::my_df_electron_density_old (gsl_vector *v, 
				       void *params, 
				       gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   if (restraints->include_map_terms() == 1) { 

      double new_S_minu, new_S_plus, tmp, val; 

      std::cout << "density_gradients" << std::endl; 
      for (unsigned int i=0; i<v->size; i++) { 
      
	 tmp = gsl_vector_get(v, i); 
	 gsl_vector_set(v, i, tmp+0.01); 
	 new_S_plus = coot::electron_density_score(v, params); 
	 gsl_vector_set(v, i, tmp-0.01); 
	 new_S_minu = coot::electron_density_score(v, params);
	 // new_S_minu = 2*tmp - new_S_plus; 

	 // restore the initial value: 
	 gsl_vector_set(v, i, tmp);

	 val = (new_S_plus - new_S_minu)/(2*0.01); 
	 std::cout << "density gradient: " << i << " " << val << std::endl;

	 // add this density term to the gradient
	 gsl_vector_set(df, i, gsl_vector_get(df, i) + val);
      } 
   }
}

// Compute both f and df together.
void coot::my_fdf(const gsl_vector *x, void *params, 
		  double *f, gsl_vector *df) { 

   // 20170423 these can be done in parallel? ... check the timings at least.
   *f = coot::distortion_score(x, params); 
    coot::my_df(x, params, df);
}


unsigned int
coot::restraints_container_t::test_function(const coot::protein_geometry &geom) {
   std::cout << "----- test_function() with geom of size : " << geom.size() << std::endl;
   std::cout << "    geom ref pointer " << &geom << std::endl;
   return geom.size();
} 

unsigned int
coot::restraints_container_t::const_test_function(const coot::protein_geometry &geom) const {
   std::cout << "----- const_test_function() with geom of size : " << geom.size() << std::endl;
   std::cout << "    geom ref pointer " << &geom << std::endl;
   return geom.size();
} 

void
coot::restraints_container_t::set_dist_crit_for_bonded_pairs(float dist) {
   dist_crit_for_bonded_pairs = dist;
}

int
coot::restraints_container_t::something_like_make_restraints(
                                              // int imol
					      // bool do_residue_internal_torsions,
					      // bool do_trans_peptide_restraints,
					      // float rama_plot_target_weight,
					      // bool do_rama_plot_restraints,
					      // bool do_auto_helix_restraints,
					      // bool do_auto_strand_restraints,
					      // coot::pseudo_restraint_bond_type sec_struct_pseudo_bonds,
					      // bool do_link_restraints,
					      // bool do_flank_restraints
					      ) {

   return df_by_thread_results_size();

}


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
					      coot::pseudo_restraint_bond_type sec_struct_pseudo_bonds,
					      bool do_link_restraints,
					      bool do_flank_restraints) {


#if 1
   make_restraints_ng(imol, geom, flags_in, do_residue_internal_torsions, do_trans_peptide_restraints,
		      rama_plot_target_weight, do_rama_plot_restraints,
		      do_auto_helix_restraints, do_auto_strand_restraints,
		      sec_struct_pseudo_bonds, do_link_restraints, do_flank_restraints);

   return restraints_vec.size();
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

   if (n_atoms) {

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
      int iret_prev = restraints_vec.size();

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
	       std::cout << "INFO:: make_restraints(): made " << n_nbcr
			 << " non-bonded restraints\n";
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
   return restraints_vec.size();
}

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
void coot::restraints_container_t::make_distortion_electron_density_ranges() {

   // std::vector<std::pair<unsigned int, unsigned int> > ranges =
   // atom_index_ranges(n_atoms, restraints_p->n_threads);

   unsigned int nt = n_threads;
   if (nt == 0) nt = 1;
   m_atom_index_ranges = atom_index_ranges(n_atoms, nt);

}
#endif



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

   for (unsigned int i=0; i<restraints_vec.size(); i++) {
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



// This only marks the first OXT atom we find that has all its
// reference atoms (which is reasonable (I hope)).
// 
void 
coot::restraints_container_t::mark_OXT(const coot::protein_geometry &geom) {

   std::string oxt(" OXT");
   for (int i=0; i<n_atoms; i++) { 
      if (std::string(atom[i]->name) == oxt) {

	 mmdb::Residue *residue = atom[i]->residue;
	 mmdb::Atom *res_atom = NULL;
	 
	 std::string res_name = residue->GetResName();
	 if (coot::util::is_standard_residue_name(res_name)) {
	    // add it if it has not been added before.
	    if (std::find(residues_with_OXTs.begin(),
			  residues_with_OXTs.end(),
			  residue) == residues_with_OXTs.end()) {
	       residues_with_OXTs.push_back(residue);
	       have_oxt_flag = true; // added 20160612 for Miguel's flying OXT bug.
	    }
	 }
      }
   }
}

bool
coot::restraints_container_t::fixed_check(int index_1) const {

   bool r = false;
   if (fixed_atom_indices.find(index_1) != fixed_atom_indices.end())
      r = true;
   return r;
} 

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2) const {

   std::vector<bool> r(2, false);
//    for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
//       if (index1 == fixed_atom_indices[ifixed])
// 	 r[0] = true;
//       if (index2 == fixed_atom_indices[ifixed])
// 	 r[1] = true;
//    }
   if (fixed_atom_indices.find(index1) != fixed_atom_indices.end())
      r[0] = true;
   if (fixed_atom_indices.find(index2) != fixed_atom_indices.end())
      r[1] = true;

   return r;
}

std::vector<bool>
coot::restraints_container_t::make_non_bonded_fixed_flags(int index1, int index2) const {

   std::vector<bool> r(2,false);
   bool set_0 = 0;
   bool set_1 = 0;

   if (fixed_atom_indices.find(index1) != fixed_atom_indices.end()) {
      r[0] = true;
      set_0 = true;
   }
   if (fixed_atom_indices.find(index2) != fixed_atom_indices.end()) {
      r[1] = true;
      set_1 = true;
   }

   if (set_0 && set_1) {
      return r;  // yay, fast.
   }

   if (! set_0) {
      mmdb::Residue *res = atom[index1]->residue;
      if (std::find(non_bonded_neighbour_residues.begin(),
		    non_bonded_neighbour_residues.end(),
		    res) != non_bonded_neighbour_residues.end())
	 r[0] = 1; // if we found the residue in non_bonded_neighbour_residues
                   // then that atom of that residue is fixed
   }
   if (! set_1) {
      mmdb::Residue *res = atom[index2]->residue;
      if (std::find(non_bonded_neighbour_residues.begin(),
		    non_bonded_neighbour_residues.end(),
		    res) != non_bonded_neighbour_residues.end())
	 r[1] = 1; 
   }
   return r;
} 


std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2, int index3) const {

   std::vector<bool> r(3,false);
   if (fixed_atom_indices.find(index1) != fixed_atom_indices.end()) r[0] = true;
   if (fixed_atom_indices.find(index2) != fixed_atom_indices.end()) r[1] = true;
   if (fixed_atom_indices.find(index3) != fixed_atom_indices.end()) r[2] = true;
   return r;
} 

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2, int index3, int index4) const {

   std::vector<bool> r(4,0);
   if (fixed_atom_indices.find(index1) != fixed_atom_indices.end()) r[0] = true;
   if (fixed_atom_indices.find(index2) != fixed_atom_indices.end()) r[1] = true;
   if (fixed_atom_indices.find(index3) != fixed_atom_indices.end()) r[2] = true;
   if (fixed_atom_indices.find(index4) != fixed_atom_indices.end()) r[3] = true;
   return r;
}

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(const std::vector<int> &indices) const {

   std::vector<bool> r(indices.size(), 0);
   for (unsigned int i_index=0; i_index<indices.size(); i_index++) {
      if (fixed_atom_indices.find(indices[i_index]) != fixed_atom_indices.end())
	 r[i_index] = true;
   }
   return r;
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

   float pseudo_bond_esd = 0.05;
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

            std::cout << "INFO:: Alpha Helix Bond restraint ("
               << at_1->name << " " << at_1->GetSeqNum() << ") to ("
               << at_3->name << " " << at_3->GetSeqNum() << ") " << ideal_dist_i_3 << std::endl;
            std::cout << "INFO:: Alpha Helix Bond restraint ("
               << at_1->name << " " << at_1->GetSeqNum() << ") to ("
               << at_2->name << " " << at_2->GetSeqNum() << ") " << ideal_dist_i_4 << std::endl;
            n_helical_restraints += 2;
         } else {
	    if (at_1 && at_3) {
	       at_1->GetUDData(udd_atom_index_handle, index_1);
	       at_3->GetUDData(udd_atom_index_handle, index_3);
	       std::vector<bool> fixed_flags_2 = make_fixed_flags(index_1, index_3);
	       double ideal_dist_i_3 = 3.18;
	       add(BOND_RESTRAINT, index_1, index_3, fixed_flags_2, ideal_dist_i_3, pseudo_bond_esd, 1.2);

	       std::cout << "INFO:: Alpha Helix Bond restraint ("
			 << at_1->name << " " << at_1->GetSeqNum() << ") to ("
			 << at_3->name << " " << at_3->GetSeqNum() << ") " << ideal_dist_i_3 << std::endl;
	       n_helical_restraints += 1;
	    }
	 }
      }

   }
   std::cout << "INFO:: added " << n_helical_restraints << " helical restraints" << std::endl;

   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
   std::cout << "INFO:: Timing for auto-helix " << d10 << " microseconds" << std::endl;

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

   float pseudo_bond_esd = 0.035; // seems reasonable.

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
      std::cout << "created " << restraints_vec.size() << " restraints" << std::endl;
      std::cout << std::endl;
   }
   return iret; // return 1 on success.  Hmm... how is this set? (and subsequently used?)
}

int
coot::restraints_container_t::make_monomer_restraints_from_res_vec(int imol,
								   const coot::protein_geometry &geom,
								   bool do_residue_internal_torsions) {

   bool print_summary = true;
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
      std::cout << "INFO:: created " << restraints_vec.size() << " restraints" << std::endl;
      std::cout << std::endl;
      if (print_summary)
	 sum.report(do_residue_internal_torsions);
   }
   return iret;
}



coot::restraints_container_t::restraint_counts_t
coot::restraints_container_t::make_monomer_restraints_by_residue(int imol, mmdb::Residue *residue_p,
								 const protein_geometry &geom,
								 bool do_residue_internal_torsions) {

   restraint_counts_t local;

   if (! residue_p) {
      std::cout << "ERROR in make_monomer_restraints_by_residue() null residue" << std::endl;
      return local;
   }

   int i_no_res_atoms;
   mmdb::PPAtom res_selection = NULL;
   std::string pdb_resname(residue_p->name);
   if (pdb_resname == "UNK") pdb_resname = "ALA";

   if (false)
      std::cout << "--------------- make_monomer_restraints_by_residue() called "
                << residue_spec_t(residue_p) << " with " << residue_p->GetNumberOfAtoms() << " atoms "
                <<  " and using type :" << pdb_resname << ": and imol "
                << imol << " do_residue_internal_torsions: "
                << do_residue_internal_torsions << std::endl;

   // idr: index dictionary residue
   int idr = geom.get_monomer_restraints_index(pdb_resname, imol, false);
   if (idr >= 0) {

      const dictionary_residue_restraints_t &dict = geom[idr].second;

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
							residue_p, geom);
	    }
	 }

	 if (restraints_usage_flag & PLANES_MASK)
	    local.n_plane_restraints += add_planes(idr, res_selection, i_no_res_atoms,
						   residue_p, geom);

	 if (restraints_usage_flag & CHIRAL_VOLUME_MASK) {
	    local.n_chiral_restr += add_chirals(idr, res_selection, i_no_res_atoms, 
						residue_p, geom);
	 }

	 restraint_counts_t mod_counts =
	    apply_mods(idr, res_selection, i_no_res_atoms, residue_p, geom);
	 // now combine mod_counts with local
      }
   }

   // local.report(false);
   return local;
}


// return in millisecs
// double 
// coot::restraints_container_t::time_diff(const timeval &current, const timeval &start) const {

//    double d = current.tv_sec - start.tv_sec;
//    d *= 1000.0;
//    d += double(current.tv_usec - start.tv_usec)/1000.0;
//    return d;


// Need to add test that residues are linked with trans
//
std::vector<coot::rama_triple_t>
coot::restraints_container_t::make_rama_triples(int SelResHnd,
						const coot::protein_geometry &geom) const {
   std::vector<coot::rama_triple_t> v;
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);

   for (int i=0; i<(nSelResidues-2); i++) {
      if (SelResidue[i] && SelResidue[i+1] && SelResidue[i+2]) {
	 std::string link_type = "TRANS";
	 rama_triple_t t(SelResidue[i], SelResidue[i+1], SelResidue[i+2], link_type);
	 v.push_back(t);
      }
   }
   return v;
}


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

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_by_linear(int SelResHnd,
							const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t c;
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);
   if (nSelResidues > 1) {
   
      std::string link_type("TRANS");
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage

      for (int i=0; i<(nSelResidues-1); i++) {
	 if (SelResidue[i] && SelResidue[i+1]) {


	    // There are a couple of ways we allow links.  First, that
	    // the residue numbers are consecutive.
	    //
	    // If there is a gap in residue numbering, the C and N
	    // have to be within 3A.
	    
	    // Are these residues neighbours?  We can add some sort of
	    // ins code test here in the future.
	    if ( abs(SelResidue[i]->GetSeqNum() - SelResidue[i+1]->GetSeqNum()) <= 1) { 
	       link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
	       // std::cout << "DEBUG ------------ in bonded_residues_by_linear() link_type is :"
	       // << link_type << ":" << std::endl;
	       if (link_type != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					whole_first_residue_is_fixed,
					whole_second_residue_is_fixed, link_type);
		  c.try_add(p);
	       }
	    }

	    // distance check this one.... it could be the opposite of
	    // an insertion code, or simply a gap - and we don't want
	    // to make a bond for a gap.
	    // 
	    if (abs(SelResidue[i]->index - SelResidue[i+1]->index) <= 1) {
	       // link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
	       std::pair<std::string, bool> link_info =
		  find_link_type_complicado(SelResidue[i], SelResidue[i+1], geom);
	       if (false)
		  std::cout << "DEBUG:: ---------- in bonded_residues_by_linear() link_info is :"
			    << link_info.first << " " << link_info.second << ":" << std::endl;
	       
	       if (link_info.first != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  if (link_info.second == 0) { 
		     coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_info.first);
		     c.try_add(p);
		  } else {
		     // order switch
		     coot::bonded_pair_t p(SelResidue[i+1], SelResidue[i],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_info.first);
		     c.try_add(p);
		  }
	       }
	    }
	 }
      }
   }
   return c;
}

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_from_res_vec(const coot::protein_geometry &geom) const {

   bool debug = false; // Are your residues in the same chain?  If not filter() will not bond them.

   coot::bonded_pair_container_t bpc;

   if (verbose_geometry_reporting == VERBOSE)
      debug = true;

   int nres = residues_vec.size();

   if (debug) {
      std::cout << "debug:: ------------------- bonded_residues_from_res_vec() residues_vec.size() "
		<< residues_vec.size() << std::endl;
      for (unsigned int i=0; i<residues_vec.size(); i++) {
	 std::cout << "   fixed: " << residues_vec[i].first << " spec: "
		   << residue_spec_t(residues_vec[i].second) << std::endl;
      }
      for (unsigned int ii=0; ii<residues_vec.size(); ii++) {
	 mmdb::Residue *res_f = residues_vec[ii].second;
	 for (unsigned int jj=ii+1; jj<residues_vec.size(); jj++) {
	    mmdb::Residue *res_s = residues_vec[jj].second;

	    std::cout << "debug:: ------------ test here with res_f and res_s "
		      << residue_spec_t(res_f) << " " << residue_spec_t(res_s) << std::endl;

	 }
      }
   }

   for (unsigned int ii=0; ii<residues_vec.size(); ii++) {
      mmdb::Residue *res_f = residues_vec[ii].second;
      for (unsigned int jj=ii+1; jj<residues_vec.size(); jj++) {
	 mmdb::Residue *res_s = residues_vec[jj].second;

	 if (debug)
	    std::cout << "debug:: ----- in bonded_resdues_from_res_vec " << residue_spec_t(res_f) << " "
		      << residue_spec_t(res_s) << "\n";

	 if (res_f == res_s) continue;

	 // 20180911 I now no longer want to evaluate closest approach here.
	 //
	 // std::pair<bool, float> d = closest_approach(res_f, res_s);
	 // Linking should be resolved by find_link_type_complicado(), not
	 // here by distance between residues.

	 std::pair<std::string, bool> l = find_link_type_complicado(res_f, res_s, geom);

	 // too verbose?
	 if (false)
	    std::cout << "   INFO:: find_link_type_complicado() for: "
		      << coot::residue_spec_t(res_f) << " " << coot::residue_spec_t(res_s)
		      << " returns link_type -> \"" << l.first << "\"" << std::endl;

	 std::string link_type = l.first;
	 if (!link_type.empty()) {

	    bool whole_first_residue_is_fixed = 0;
	    bool whole_second_residue_is_fixed = 0;
	    bool order_switch_flag = l.second;

	    if (!order_switch_flag) {
	       coot::bonded_pair_t p(res_f, res_s,
				     whole_first_residue_is_fixed,
				     whole_second_residue_is_fixed, link_type);
	       bool previously_added_flag = bpc.try_add(p);
	    } else {
	       coot::bonded_pair_t p(res_s, res_f,
				     whole_first_residue_is_fixed,
				     whole_second_residue_is_fixed,
				     link_type);
	       bool previously_added_flag = bpc.try_add(p);
	       // std::cout << "              previously_added_flag " << previously_added_flag
	       // << std::endl;
	    }

	    // if the link type is a straight-forward TRANS link of 2 residue next to each other
	    // in residue numbers and serial numbers, then we don't need to find any other type
	    // of link for this residue (so break out of the inner for-loop).
	    bool was_straight_forward_trans_link = false;
	    int resno_1 = res_f->GetSeqNum();
	    int resno_2 = res_s->GetSeqNum();
	    int ser_num_1 = res_f->index;
	    int ser_num_2 = res_s->index;
	    if (resno_2 == (resno_1 + 1)) {
	       if (ser_num_2 == (ser_num_1 + 1)) {
		  std::string rn_1 = res_f->GetResName();
		  if (rn_1 != "ASN" && rn_1 != "CYS" && rn_1 != "SER" && rn_1 != "TYR") {
		     if (link_type == "TRANS" || link_type == "PTRANS")
			was_straight_forward_trans_link = true;
		  }
	       }
	    }
	    if (was_straight_forward_trans_link) {
	       // std::cout << "------------ was straight_forward TRANS link! - breaking"  << std::endl;
	       break;
	    }

	 } else {
	    if (debug)
	       std::cout << "DEBUG:: find_link_type_complicado() blank result: "
			 << "link_type find_link_type_complicado() for "
			 << coot::residue_spec_t(res_f) << " " << coot::residue_spec_t(res_s)
			 << " returns \"" << l.first << "\" " << l.second << std::endl;
	 }
      }
   }

   bpc.filter(); // removes 1-3 bond items and if 1-2 and 1-3 bonds exist

   // std::cout << "---------------- done bonded_residues_from_res_vec()" << std::endl;

   return bpc;
}

// a pair, first is if C and N are close
//       using enum peptide_order_info_t { IS_PEPTIDE=1, IS_NOT_PEPTIDE=0, UNKNOWN=-1 }
//
// and second if and order switch is needed to make it so.
std::pair<coot::restraints_container_t::peptide_order_info_t, bool>
coot::restraints_container_t::peptide_C_and_N_are_in_order_p(mmdb::Residue *r1, mmdb::Residue *r2) const {

   // If the residues are next to each other in serial and residue number then it's a peptide, no
   // matter how far apart they are.
   //
   // If that is not the case, then sometimes we don't know because this might be a residues pair
   // with an insertion code - and in that case, a different check should be used.

   bool debug = false;
   if (r1->chain == r2->chain) {
      int serial_delta = r2->index - r1->index;
      if (debug)
	 std::cout << "   serial_delta " << serial_delta << std::endl;
      if ((serial_delta == -1) || (serial_delta == 1)) {
	 // ok to proceed
      } else {
	 if (debug)
	    std::cout << "   ------ peptide_C_and_N_are_in_order_p path : A0 - "
		      << "same chain not sequencial" << std::endl;
	 return std::pair<peptide_order_info_t, bool> (IS_NOT_PEPTIDE, false);
      }

      if (serial_delta == 1) {
	 if (debug)
	    std::cout << "   ------ peptide_C_and_N_are_in_order_p path A" << std::endl;
	 std::string ins_code_1 = r1->GetInsCode();
	 std::string ins_code_2 = r2->GetInsCode();
	 int res_no_delta = r2->GetSeqNum() - r1->GetSeqNum();
	 if (ins_code_1 == "") {
	    if (ins_code_2 == "") {
	       if (res_no_delta == 1 || res_no_delta == -1) {
		  return std::pair<peptide_order_info_t, bool> (IS_PEPTIDE, false);
	       }
	    }
	 }
	 if (debug)
	    std::cout << "   ------ peptide_C_and_N_are_in_order_p path A-unk" << std::endl;
	 return std::pair<peptide_order_info_t, bool> (UNKNOWN, false);

      } else {
	 if (debug)
	    std::cout << "   ------ peptide_C_and_N_are_in_order_p path B" << std::endl;

	 std::string ins_code_1 = r1->GetInsCode();
	 std::string ins_code_2 = r2->GetInsCode();
	 int res_no_delta = r2->GetSeqNum() - r1->GetSeqNum();
	 if (ins_code_1 == "") {
	    if (ins_code_2 == "") {
	       if (res_no_delta == 1 || res_no_delta == -1) {
		  return std::pair<peptide_order_info_t, bool> (IS_PEPTIDE, true);
	       }
	    }
	 }
	 if (debug)
	    std::cout << "   ------ peptide_C_and_N_are_in_order_p path B-unk" << std::endl;
	 return std::pair<peptide_order_info_t, bool> (UNKNOWN, true);
      }

   } else {
      // we are considering a link between a residue in the mol and a residue
      // of the neighbouring residues vectors (which are not residues in the mol(!))
      // i.e. the don't have the same indexing (residue serial indexing) scheme.

      // we can't make a decision. We need to be able to tell the caller that - so
      // that the caller can choose to bond the residues by distance (and residue number
      // and insertion code)

      return std::pair<peptide_order_info_t, bool> (UNKNOWN, false);
   }
}


// a pair, first is if C and N are close and second if and order
// switch is needed to make it so.
std::pair<bool, bool>
coot::restraints_container_t::peptide_C_and_N_are_close_p(mmdb::Residue *r1, mmdb::Residue *r2) const {

   // needs PDBv3 fixup.
   float dist_crit = 2.8; // 20150714: this used to be 2.0.  That is too long, I think.
                          //           Try 2.8.
                          // 
                          // 2.0 A for a peptide link - so that we
			  // don't find unintentional peptides - which
			  // would be a disaster.

   std::string C_atom_name = " C  ";
   std::string N_atom_name = " N  ";
   
   mmdb::Atom *at_c_1 = NULL;
   mmdb::Atom *at_n_1 = NULL;
   mmdb::Atom *at_c_2 = NULL;
   mmdb::Atom *at_n_2 = NULL;

   mmdb::PPAtom residue_atoms_1 = NULL;
   mmdb::PPAtom residue_atoms_2 = NULL;
   int n_residue_atoms_1;
   int n_residue_atoms_2;
   r1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   r2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

   for (int iat=0; iat<n_residue_atoms_1; iat++) {
      std::string atom_name(residue_atoms_1[iat]->name);
      if (atom_name == C_atom_name) {
	 at_c_1 = residue_atoms_1[iat];
      } 
      if (atom_name == N_atom_name) {
	 at_n_1 = residue_atoms_1[iat];
      } 
   }

   for (int iat=0; iat<n_residue_atoms_2; iat++) {
      std::string atom_name(residue_atoms_2[iat]->name);
      if (atom_name == C_atom_name) {
	 at_c_2 = residue_atoms_2[iat];
      } 
      if (atom_name == N_atom_name) {
	 at_n_2 = residue_atoms_2[iat];
      } 
   }

   if (at_c_1 && at_n_2) {
      clipper::Coord_orth c1(at_c_1->x, at_c_1->y, at_c_1->z);
      clipper::Coord_orth n2(at_n_2->x, at_n_2->y, at_n_2->z);
      float d = clipper::Coord_orth::length(c1, n2);
      // std::cout << "   C1->N2 dist " << d << std::endl;
      if (d < dist_crit)
	 return std::pair<bool, bool> (1,0);
   } 

   if (at_n_1 && at_c_2) {
      clipper::Coord_orth n1(at_n_1->x, at_n_1->y, at_n_1->z);
      clipper::Coord_orth c2(at_c_2->x, at_c_2->y, at_c_2->z);
      float d = clipper::Coord_orth::length(n1, c2);
      // std::cout << "   N1->C2 dist " << d << std::endl;
      if (d < dist_crit)
	 return std::pair<bool, bool> (1,1);
   } 

   return std::pair<bool, bool> (0, 0);

}



std::pair<bool,float>
coot::restraints_container_t::closest_approach(mmdb::Residue *r1, mmdb::Residue *r2) const {

   return coot::closest_approach(mol, r1, r2);
} 


// 20180224 New-style: Post Weizmann 
//
// find residues in the neighbourhood that are not in the refining set
// and are not already marked as bonded flankers.
//
// set the class variable non_bonded_neighbour_residues
void
coot::restraints_container_t::set_non_bonded_neighbour_residues_by_residue_vector(const std::map<mmdb::Residue *, std::set<mmdb::Residue *> > &neighbour_set,
										  const coot::bonded_pair_container_t &bonded_flanking_pairs, const coot::protein_geometry &geom) {

   // non_bonded_neighbour_residues becomes this:
   //
   std::vector<mmdb::Residue *> nbr; // non-bonded residues 
   float dist_crit = 3.0;

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it_map;

   // don't iterate like this:
   // for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
   // std::vector<mmdb::Residue *> neighbours =
   // coot::residues_near_residue(residues_vec[ir].second, mol, dist_crit);

   for(it_map=neighbour_set.begin(); it_map!=neighbour_set.end(); it_map++) {

      const std::set<mmdb::Residue *> &neighbours = it_map->second;
      std::set<mmdb::Residue *>::const_iterator it_set;

      for (it_set=neighbours.begin(); it_set!=neighbours.end(); it_set++) {
	 mmdb::Residue *test_res = *it_set;
	 if (std::find(nbr.begin(), nbr.end(), test_res) == nbr.end()) {
	    // not already there...
	    bool found = false;

	    if (false) // debug
	       std::cout << ".... about to compare " << residue_spec_t(test_res) << " to "
			 << residues_vec.size() << " refining residues " << std::endl;
	    for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	       if (test_res == residues_vec[ires].second) {
		  found = true;
		  break;
	       }
	    }

	    if (! found) {
	       // OK, so this neighbour was not in the passed set of
	       // moving residues (and not already in nbr)... it can
	       // be a flanking residue then...

	       // check that it is not a bonded flanking residue...
	       for (unsigned int iflank=0; iflank<bonded_flanking_pairs.size(); iflank++) { 
		  if (bonded_flanking_pairs[iflank].res_1 == test_res) {
		     found = 1;
		     // std::cout << "      oops bonded flanking residue res1 " << std::endl;
		     break;
		  } 
		  if (bonded_flanking_pairs[iflank].res_2 == test_res) {
		     found = 1;
		     // std::cout << "   oops bonded flanking residue res2 " << std::endl;
		     break;
		  }
	       }

	       if (! found) {
		  // std::cout << ".... adding non-bonded neighbour " << residue_spec_t(test_res) << std::endl;
		  nbr.push_back(test_res);
	       }
	    }
	 }
      }
   }
   non_bonded_neighbour_residues = nbr;
}

// 20180224 pre-Weizmann
//
// find residues in the neighbourhood that are not in the refining set
// and are not already marked as bonded flankers.
// 
void
coot::restraints_container_t::set_non_bonded_neighbour_residues_by_residue_vector(const coot::bonded_pair_container_t &bonded_flanking_pairs, const coot::protein_geometry &geom) {

   std::vector<mmdb::Residue *> nbr; // non-bonded residues 
   float dist_crit = 3.0;

   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
      std::vector<mmdb::Residue *> neighbours =
	 coot::residues_near_residue(residues_vec[ir].second, mol, dist_crit);
      for (unsigned int ineighb=0; ineighb<neighbours.size(); ineighb++) {
	 mmdb::Residue *test_res = neighbours[ineighb];
	 if (std::find(nbr.begin(), nbr.end(), test_res) == nbr.end()) {
	    // not already there...
	    bool found = 0;

	    if (0)
	       std::cout << ".... about to compare " << residue_spec_t(test_res) << " to "
			 << residues_vec.size() << " refining residues " << std::endl;
	    for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	       if (test_res == residues_vec[ires].second) {
		  found = 1;
		  break;
	       }
	    }

	    if (! found) {
	       // OK, so this neighbour was not in the passed set of
	       // moving residues (and not already in nbr)... it can
	       // be a flanking residue then...

	       // check that it is not a bonded flanking residue...
	       for (unsigned int iflank=0; iflank<bonded_flanking_pairs.size(); iflank++) { 
		  if (bonded_flanking_pairs[iflank].res_1 == test_res) {
		     found = 1;
		     // std::cout << "      oops bonded flanking residue res1 " << std::endl;
		     break;
		  } 
		  if (bonded_flanking_pairs[iflank].res_2 == test_res) {
		     found = 1;
		     // std::cout << "   oops bonded flanking residue res2 " << std::endl;
		     break;
		  }
	       }

	       if (! found) {
		  // std::cout << ".... adding non-bonded neighbour " << residue_spec_t(test_res) << std::endl;
		  nbr.push_back(test_res);
	       } 
	    }
	 }
      }
   }
   non_bonded_neighbour_residues = nbr;
} 

bool
coot::restraints_container_t::H_parent_atom_is_donor(mmdb::Atom *at) {

   bool state = false;
   std::map<mmdb::Atom *, hb_t>::const_iterator it;
   it = H_atom_parent_energy_type_atom_map.find(at);
   if (it != H_atom_parent_energy_type_atom_map.end()) {
      // found it
      const hb_t &hbt = it->second;
      if (hbt == HB_DONOR || hbt == HB_BOTH)
	 state = true;
   } else {
      // not found
   }
   return state;
}

int 
coot::restraints_container_t::make_non_bonded_contact_restraints(int imol, const coot::bonded_pair_container_t &bpc,
								 const coot::protein_geometry &geom) {

   // is this function used any more?
   //
   coot::restraints_container_t::reduced_angle_info_container_t ai(restraints_vec);
   ai.write_angles_map("angles_map.tab");
   return make_non_bonded_contact_restraints(imol, bpc, ai, geom);

} 

#include "coot-utils/contacts-by-bricks.hh"

// Atoms that are not involved in bonds or angles, but are in the
// residue selection should be at least 2.7A away from each other.
// 
// Here are my anti-bumping notes:
//
//
//     Anti-bumping restraints in regularization:
//
//     Considering totally screwed-up geometry: We should add a strong
//     repulsion for atoms that are not bonded so that they go away from
//     each other.  
//
//     Something like a triangle function between 0->2A and 0 beyond that.
//
//     Each atom has to check the distance to each other atom in the
//     selection: if they are not bonded, get a repulsion score for that
//     distance.
//
//     Derivative of that should be not too tricky, similar to bonds, but
//     not the same.
//
//     Instead of 500, use 10*matrix?  Doesn't really matter, I think.
//
//     Instead of using 2.0 as the critical distance, let's instead use
//     d_crit:
//
//     Infact,
//
//     f = 1000-d*(1000/d_crit)   for d<d_crit
//     df/dd = -1000/d_crit
//     df/dx = df/dd dd/dx
//           = -1000/d_crit
//
//     It's like bonds:
//     d = sqrt[ (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 ]
//     => dd/dx = (xi-xj)/d
//
//     So df/dx = -1000/d_crit * (xi-xj)/d
//
//     Need to keep a list of repulsing atom pairs so that we don't have
//     to calculate them once each for distortion_score and derivates..?
//
// Note that if residue-2 is not moving then it will not have angle restraints.  If
// it doesn't have angle resraints then the is_1_4_related test will fail.
// e.g (if n-1 is fixed residue): C(n-1)-N(n)-Ca(n)-C(n) or C(n-1)-N(n)-Ca(n)-CB(n)
// will not be seen as 1-4 related. So that's where strange_exception comes in.
//
int
coot::restraints_container_t::make_non_bonded_contact_restraints(int imol, const coot::bonded_pair_container_t &bpc,
								 const coot::restraints_container_t::reduced_angle_info_container_t &ai,
								 const coot::protein_geometry &geom) {


#ifdef HAVE_CXX_THREAD

   std::cout << "------------------- timing" << std::endl;
   std::set<unsigned int> fixed_atom_flags_set; // fill this properly!
   auto tp_0 = std::chrono::high_resolution_clock::now();
   contacts_by_bricks cb(atom, n_atoms, fixed_atom_flags_set);
   auto tp_1 = std::chrono::high_resolution_clock::now();
   std::vector<std::set<unsigned int> > vcontacts;
   cb.find_the_contacts(&vcontacts);
   auto tp_2 = std::chrono::high_resolution_clock::now();
   cb.find_the_contacts(&vcontacts);
   auto tp_3 = std::chrono::high_resolution_clock::now();

   unsigned int n_nbc = 0;
   for (std::size_t ii=0; ii<vcontacts.size(); ii++) {
      n_nbc += vcontacts.at(ii).size();
   }
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   std::cout << "------------------- timing: " << d10 << " " << d21 << " " << d32
	     << " milliseconds for " << n_nbc << " nbcs " << std::endl;

#endif // HAVE_CXX_THREAD

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map_cache;
   construct_non_bonded_contact_list(bpc, geom);

   // so now filtered_non_bonded_atom_indices is filled.
   // but it is not necessarily symmetric - so we can't do a j > 1 test (yet).
   // 
   // Write a debug/test for symmetry of filtered_non_bonded_atom_indices.
   symmetry_non_bonded_contacts(0);

   // We need to find if atom pairs are in the same ring.
   // We do that by finding the restraints of each residue and putting them in the map.
   // To make things faster in cases where the restraints look-up fails, we add a flag 
   // to the value of the map which let's us know that we have searched this dictionary 
   // type before.
   std::map<mmdb::Residue *, std::pair<bool, dictionary_residue_restraints_t> > restraints_map;

   if (false) {
      std::cout << "--------- make_non_bonded_contact_restraints() the atom array: " << std::endl;
      for (int iat=0; iat<n_atoms; iat++)
	 std::cout << "------- " << iat << " " << atom_spec_t(atom[iat]) << std::endl;
   }

   // Thinking of setting this to true? is the (link) angle in the dictionary? Is one of the
   // residues non-moving? (see above notes).
   if (false)
      ai.write_angles_map("angles-map.tab");
   
   if (false) {
      std::cout << "--------------------------------------------------\n";
      std::cout << "   non-bonded list:" << std::endl;
      std::cout << "--------------------------------------------------\n";
      for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
	 std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
	 for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) { 
	    std::cout << filtered_non_bonded_atom_indices[i][j] << " ";
	 } 
	 std::cout << std::endl;
      } 
      std::cout << "--------------------------------------------------\n";
   }

   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
      mmdb::Atom *at = atom[i];
      std::string res_type = at->GetResName();
      std::map<mmdb::Residue *, std::pair<bool, dictionary_residue_restraints_t> >::const_iterator it;
      it = restraints_map.find(at->residue);
      if (it == restraints_map.end()) {
	 // have_restraints_for() is faster?
	 std::pair<bool, dictionary_residue_restraints_t> p = geom.get_monomer_restraints(res_type, imol);
	 // p.first is false if this is not a filled dictionary
	 restraints_map[at->residue] = p;
      }
   }

   // cache the energy types:
   std::map<mmdb::Atom *, std::string> energy_type_cache;
   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) {
      mmdb::Atom *at = atom[i];
      energy_type_cache[at] = get_type_energy(imol, at, geom);
   }

   int n_nbc_r = 0;
   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
      for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) {

	 mmdb::Atom *at_1 = atom[i];
	 mmdb::Atom *at_2 = atom[filtered_non_bonded_atom_indices[i][j]];

	 std::vector<bool> fixed_atom_flags =
	    make_non_bonded_fixed_flags(i, filtered_non_bonded_atom_indices[i][j]);

	 timeval start_time;
	 timeval current_time;
	 double d;
         if (at_1 && at_2) {

	    std::string type_1 = energy_type_cache[at_1];
	    std::string type_2 = energy_type_cache[at_2];

	    bool add_it = true;

	    // no H-H bumps in the same residue
	    //
	    // [20131212: Why not?  I suppose that there was a reason,
	    // it is not clear to me what it is now].  This needs to
	    // be investigated/fixed.
	    //
	    if (at_2->residue == at_1->residue)
	       if (is_hydrogen(at_1))
	          if (is_hydrogen(at_2))
		     add_it = false;

   	    if (filtered_non_bonded_atom_indices[i][j] < int(i))
  	       add_it = false;

	    int res_no_1 = at_1->GetSeqNum();
	    int res_no_2 = at_2->GetSeqNum();

	    if (add_it) { 

	       // Don't make a bump between the CD of a PRO at residue(n) and the atoms of n-1
	    
	       std::string res_name_1 = at_1->GetResName();
	       std::string res_name_2 = at_2->GetResName();

	       if (false)
		  std::cout << "DEBUG:: here with " << atom_spec_t(at_1) << " " << atom_spec_t(at_2)
			    << " res_names " << res_name_1 << " " << res_name_2 << " "
			    << at_1->GetAtomName() << " " << at_2->GetAtomName()
			    << std::endl;

	       if (res_name_1 == "PRO" || res_name_1 == "HYP") {
		  int res_no_pro   = res_no_1;
		  int res_no_other = res_no_2;
		  if (res_no_pro == (res_no_other + 1)) {
		     std::string atom_name = at_1->name;
		     if (atom_name == " CD ") {  // PDBv3 FIXME
			add_it = false;
		     }
		  }
	       }
	       if (res_name_2 == "PRO" || res_name_2 == "HYP") {
		  int res_no_pro   = res_no_2;
		  int res_no_other = res_no_1;
		  if (res_no_pro == (res_no_other + 1)) {
		     std::string atom_name = at_2->name;
		     if (atom_name == " CD ") {  // PDBv3 FIXME
			add_it = false;
		     } 
		  }
	       }

	       // hack to remove C1-OD1 NBC on N-linked glycosylation
	       //
	       if (res_name_1 == "ASN" || res_name_2 == "NAG") {
		  std::string atom_name_1(at_1->name);
		  std::string atom_name_2(at_2->name);
		  if (atom_name_1 == " OD1")
		     if (atom_name_2 == " C1 ")
			add_it = false;
		  if (atom_name_1 == "HD21")
		     if (atom_name_2 == " C1 ")
			add_it = false;
	       }

	       if (res_name_1 == "NAG" || res_name_2 == "ASN") {
		  std::string atom_name_1(at_1->name);
		  std::string atom_name_2(at_2->name);
		  if (atom_name_1 == " C1 ")
		     if (atom_name_2 == " OD1")
			add_it = false;
		  if (atom_name_1 == " C1 ")
		     if (atom_name_2 == "HD21")
			add_it = false;
	       }
	    }

	    // -------------- OK add_it was set -----
	    
	    if (add_it) {

	       double dist_min = 3.4;

	       bool in_same_ring_flag    = true;
	       bool in_same_residue_flag = true;
	       
	       if (at_2->residue != at_1->residue) {
		  in_same_ring_flag    = false;
		  in_same_residue_flag = false;
	       }
	       
	       if (in_same_ring_flag) {
		  std::string atom_name_1 = at_1->GetAtomName();
		  std::string atom_name_2 = at_2->GetAtomName();

		  // in_same_ring_flag = restraints_map[at_2->residue].second.in_same_ring(atom_name_1,
		  //                                                                       atom_name_2);

		  in_same_ring_flag = is_in_same_ring(imol, at_2->residue,
						      residue_ring_map_cache,
						      atom_name_1, atom_name_2, geom);
	       }
	       
	       // this doesn't check 1-4 over a moving->non-moving peptide link (see comment above function)
	       // because the non-moving atom doesn't have angle restraints.
	       //
	       bool is_1_4_related = ai.is_1_4(i, filtered_non_bonded_atom_indices[i][j]);

	       if (false)
		  std::cout << "here C with at_1 " << atom_spec_t(at_1) << " at_2 " << atom_spec_t(at_2)
			    << " is_1_4_related " << is_1_4_related << std::endl;

	       if (is_1_4_related) {
		  dist_min = 2.64; // was 2.7 but c.f. guanine ring distances
		  if (is_hydrogen(at_1))
		      dist_min -= 0.7;
		  if (is_hydrogen(at_2))
		      dist_min -= 0.7;
	       } else {

		  std::pair<bool, double> nbc_dist = geom.get_nbc_dist(type_1, type_2,
								       in_same_residue_flag,
								       in_same_ring_flag);

		  if (nbc_dist.first) {

		     // In a helix O(n) is close to C(n+1), we should allow it.
		     // 
		     bool is_O_C_1_5_related = check_for_O_C_1_5_relation(at_1, at_2);

		     if (is_O_C_1_5_related) {
			dist_min = 2.84;
		     } else {

			// Perhaps we don't have angle restraints to both atoms because one
			// of the atoms is fixed (and thus miss that these have a 1-4 relationship).
			// e.g. O(n) [moving] -> CA(n+1) [fixed]
			// 
			// (this test will fail on insertion codes)
			//

			bool strange_exception = false;
			int rn_diff = abs(res_no_2 - res_no_1);
			if (rn_diff == 1) {
			   std::string atom_name_1 = at_1->GetAtomName();
			   std::string atom_name_2 = at_2->GetAtomName();
			   if (fixed_atom_flags.size()) {
			      if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
				 if (atom_name_1 == " O  ")
				    if (atom_name_2 == " CA ") 
				       strange_exception = true;
				 if (atom_name_1 == " CA ")
				    if (atom_name_2 == " O  ")
				       strange_exception = true;
				 if (atom_name_1 == " N  ")
				    if (atom_name_2 == " CB ")
				       strange_exception = true;
				 if (atom_name_1 == " CB ")
				    if (atom_name_2 == " N  ")
				       strange_exception = true;
				 if (atom_name_1 == " C  ")
				    if (atom_name_2 == " CB ")
				       strange_exception = true;
			      }
			   }
			   if (strange_exception)
			      dist_min = 2.7;

			   // Strange that these are not marked as 1-4 related.  Fix here...
			   // HA-CA-N-C can be down to ~2.4A.
			   // HA-CA-C-N can be down to ~2.41A.
			   if (res_no_2 > res_no_1) {
			      if (atom_name_1 == " C  ") {
				 if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			      if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
				 if (atom_name_2 == " N  ") {
				    strange_exception = true;
				    dist_min = 2.41;
				 }
			      }
			      if (atom_name_1 == " N  ") {
				 if (atom_name_2 == " H  ") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			   } else {
			      if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
				 if (atom_name_2 == " C  ") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			      if (atom_name_1 == " N  ") {
				 if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
				    strange_exception = true;
				    dist_min = 2.41;
				 }
			      }
			      if (atom_name_2 == " N  ") {
				 if (atom_name_1 == " H  ") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			   }
			}
			if (rn_diff == 2) { 
			   if (fixed_atom_flags.size()) {
			      if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
				 std::string atom_name_1 = at_1->GetAtomName();
				 std::string atom_name_2 = at_2->GetAtomName();
				 if (atom_name_1 == " C  ")
				    if (atom_name_2 == " N  ")
				       strange_exception = true;
				 if (atom_name_1 == " N  ")
				    if (atom_name_2 == " C  ")
				       strange_exception = true; // 3.1 would be enough

				 if (strange_exception)
				    dist_min = 2.7;
			      }
			   }
			}

			if (! strange_exception)
			   dist_min = nbc_dist.second;
		     }
		  } else {
		     // short/standard value
		     dist_min = 2.8;
		  }
	       }

	       if (false) { // debug.
	          clipper::Coord_orth pt1(atom[i]->x, atom[i]->y, atom[i]->z);
	          clipper::Coord_orth pt2(at_2->x,    at_2->y,    at_2->z);
	          double dd = sqrt((pt1-pt2).lengthsq());

		  std::pair<bool, double> nbc_dist = geom.get_nbc_dist(type_1, type_2,
								       in_same_residue_flag,
								       in_same_ring_flag);

	          std::cout << "adding non-bonded contact restraint index " 
			    << i << " to index " << filtered_non_bonded_atom_indices[i][j]
			    << " "
			    << atom_spec_t(atom[i]) << " to "
			    << atom_spec_t(atom[filtered_non_bonded_atom_indices[i][j]])
			    << "  types: " << type_1 <<  " " << type_2 <<  " fixed: "
			    << fixed_atom_flags[0] << " " << fixed_atom_flags[1] << "   current: " << dd
			    << " dist_min: " << dist_min << " using nbc dist " << nbc_dist.second
			    << "\n";
	       }

	       bool is_H_non_bonded_contact = false;

	       if (is_hydrogen(at_1)) {
		  is_H_non_bonded_contact = true;
		  if (H_parent_atom_is_donor(at_1))
		     if (is_acceptor(type_2, geom))
			dist_min -= 0.7;
	       }
	       if (is_hydrogen(at_2)) {
		  is_H_non_bonded_contact = true;
		  if (H_parent_atom_is_donor(at_2))
		     if (is_acceptor(type_1, geom))
			dist_min -= 0.7;
	       }


	       simple_restraint::nbc_function_t nbcf = simple_restraint::LENNARD_JONES;
	       // simple_restraint::nbc_function_t nbcf = simple_restraint::HARMONIC;
	       simple_restraint r(NON_BONDED_CONTACT_RESTRAINT,
				  nbcf,
				  i, filtered_non_bonded_atom_indices[i][j],
				  type_1, type_2, is_H_non_bonded_contact,
				  fixed_atom_flags, dist_min);

	       if (true)
		  std::cout << "Adding NBC " << i << " " << filtered_non_bonded_atom_indices[i][j]
			    << " " << type_1 << " " << type_2 << " " 
			    << is_H_non_bonded_contact << " "
			    << fixed_atom_flags[0] << " " << fixed_atom_flags[1] << " "
			    << dist_min <<  "\n";

	       r.n_atoms_from_all_restraints = n_atoms; // for debugging crash in non-bonded contact
	                                                // restraints
	       r.restraints_index = restraints_vec.size(); // likewise
	       restraints_vec.push_back(r);

	       n_nbc_r++;
	    }
	 }
      }
   }
   return n_nbc_r;
}

bool
coot::restraints_container_t::is_acceptor(const std::string &energy_type,
					  const coot::protein_geometry &geom) const {

   // get_energy_lib_atom() returns a blank atom on failure to look up energy_type
   energy_lib_atom ela = geom.get_energy_lib_atom(energy_type);
   bool acceptor_flag = ((ela.hb_type == HB_ACCEPTOR) || (ela.hb_type == HB_BOTH));
   
   return acceptor_flag;
}

// the bool in the residue_ring_map_cache is a flag that means "I've
// tried before to look this residue up and failed".
//
// static
bool
coot::restraints_container_t::is_in_same_ring(int imol, mmdb::Residue *residue_p,
					      std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > &residue_ring_map_cache,
					      const std::string &atom_name_1,
					      const std::string &atom_name_2,
					      const coot::protein_geometry &geom) {

   // maybe lists are slow.
   bool r = false;

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map;
   std::list<std::string> r1;
   std::list<std::string> r2;
   std::list<std::string> r3;
   std::list<std::string> r4;

   // HIS
   r1.push_back(" CG ");
   r1.push_back(" CD2");
   r1.push_back(" ND1");
   r1.push_back(" CE1");
   r1.push_back(" NE2");

   // PHE/TYR
   r2.push_back(" CG ");
   r2.push_back(" CD1");
   r2.push_back(" CD2");
   r2.push_back(" CE1");
   r2.push_back(" CE2");
   r2.push_back(" CZ ");

   // TRP
   r3.push_back(" CG ");
   r3.push_back(" CD1");
   r3.push_back(" CD2");
   r3.push_back(" CE2");
   r3.push_back(" NE1");
   
   r4.push_back(" CD2");
   r4.push_back(" CE2");
   r4.push_back(" CE3");
   r4.push_back(" CZ2");
   r4.push_back(" CZ3");
   r4.push_back(" CH2");

   if (residue_ring_map_cache.size() == 0) {
      r1.sort();
      r2.sort();
      r3.sort();
      r4.sort();
      residue_ring_map["HIS"].second.push_back(r1);
      residue_ring_map["PHE"].second.push_back(r2);
      residue_ring_map["TYR"].second.push_back(r2);
      residue_ring_map["TRP"].second.push_back(r3);
      residue_ring_map["TRP"].second.push_back(r4);
      residue_ring_map["HIS"].first = false;
      residue_ring_map["PHE"].first = false;
      residue_ring_map["TYR"].first = false;
      residue_ring_map["TRP"].first = false;
   }

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > >::const_iterator it;
   std::string res_name = residue_p->GetResName();

   it = residue_ring_map_cache.find(res_name);
   if (it != residue_ring_map_cache.end()) {

      if (it->second.first == 0) { // not looked up before and failed 
	 for (unsigned int i=0; i<it->second.second.size(); i++) {
	    std::list<std::string>::const_iterator it_1 = std::find(it->second.second[i].begin(), it->second.second[i].end(), atom_name_1);
	    std::list<std::string>::const_iterator it_2 = std::find(it->second.second[i].begin(), it->second.second[i].end(), atom_name_2);
	    if (it_1 != it->second.second[i].end()) {
	       if (it_2 != it->second.second[i].end()) {
		  r = true;
		  break;
	       }
	    }
	 }
      } else {
	 // We tried to look it up before and failed
      }
   } else {

      // add it then
      std::pair<bool, dictionary_residue_restraints_t> rest =
	 geom.get_monomer_restraints(res_name, imol);
      if (rest.first) {
	 std::vector<std::vector<std::string> > ri = rest.second.get_ligand_ring_list();
	 residue_ring_map_cache[res_name].first = false; // not looked up before and failed	 
	 for (unsigned int ii=0; ii<ri.size(); ii++) {
	    std::list<std::string> l;
	    for (unsigned int jj=0; jj<ri[ii].size(); jj++)
	       l.push_back(ri[ii][jj]);
	    l.sort();
	    residue_ring_map_cache[res_name].second.push_back(l);
	 }

	 std::vector<std::list<std::string> > &vl = residue_ring_map_cache[res_name].second;

	 for (unsigned int ii=0; ii<vl.size(); ii++) {
	    std::list<std::string>::const_iterator it_1 = std::find(vl[ii].begin(), vl[ii].end(), atom_name_1);
	    std::list<std::string>::const_iterator it_2 = std::find(vl[ii].begin(), vl[ii].end(), atom_name_2);
	    if (it_1 != vl[ii].end()) {
	       if (it_2 != vl[ii].end()) {
		  r = true;
		  break;
	       }
	    }
	 }
      } else {
	 // OK, the lookup failed
	 std::vector<std::list<std::string> > fv;
	 std::pair<bool, std::vector<std::list<std::string> > > failed_data(true, fv);
	 residue_ring_map_cache[res_name] = failed_data;
      }
   }
   return r;
}


bool
coot::restraints_container_t::check_for_1_4_relation(int idx_1, int idx_2,
						     const reduced_angle_info_container_t &ai) const {

   bool is_1_4 = false;
   is_1_4 = ai.is_1_4(idx_1, idx_2);
   // std::cout << "debug:: check_for_1_4_relation(ai) " << idx_1 << " " << idx_2 << " is " << is_1_4
   // << std::endl;
   return is_1_4;
}

coot::restraints_container_t::reduced_angle_info_container_t::reduced_angle_info_container_t(const std::vector<coot::simple_restraint> &r) {

   init(r);

}

void
coot::restraints_container_t::reduced_angle_info_container_t::init(const std::vector<coot::simple_restraint> &r) {

   // this map is constructed correctly.  If you are here it's because
   // you expect an angle restraint that is not there.
   // 
   for (unsigned int ii=0; ii<r.size(); ii++) {
      if (r[ii].restraint_type == ANGLE_RESTRAINT) {
	 std::pair<int, int> p_1(r[ii].atom_index_2, r[ii].atom_index_3);
	 std::pair<int, int> p_2(r[ii].atom_index_2, r[ii].atom_index_1);
	 angles[r[ii].atom_index_1].push_back(p_1);
	 angles[r[ii].atom_index_3].push_back(p_2);
      }
   }
}


void
coot::restraints_container_t::reduced_angle_info_container_t::write_angles_map(const std::string &file_name) const {

   std::ofstream f(file_name.c_str());
   if (f) {
      std::map<int, std::vector<std::pair<int, int> > >::const_iterator it;
      for (it=angles.begin(); it!=angles.end(); it++) {
	 const std::vector<std::pair<int, int> > &v = it->second;
	 for (unsigned int i=0; i<v.size(); i++) {
	    f << "key: ";
	    f << it->first;
	    f << " value ";
	    f << " " << v[i].first <<  " " << v[i].second << "\n";
	 }
      }
      f.close();
   }

} 

bool
coot::restraints_container_t::reduced_angle_info_container_t::is_1_4(int indx_1, int indx_2) const {

   // this function can be const because we don't use [] operator on the angles map.
   
   bool f = false;

   std::map<int, std::vector<std::pair<int, int> > >::const_iterator it_1, it_2;
   it_1 = angles.find(indx_1);
   if (it_1 != angles.end()) {
      const std::vector<std::pair<int, int> > &v = it_1->second;
      for (unsigned int ii=0; ii<v.size(); ii++) {
	 
	 // what are the angles that have atom_mid as atom_1?  We can ask this because angles
	 // go into this object both way rounds: A-B-C, C-B-A.
	 
	 int idx_mid = v[ii].first;

	 it_2 = angles.find(idx_mid);
	 if (it_2 != angles.end()) {
	    const std::vector<std::pair<int, int> > &v_2 = it_2->second;
	    // are any of these indx_2?
	    for (unsigned int jj=0; jj<v_2.size(); jj++) { 
	       if (v_2[jj].second == indx_2) {
		  f = true;
		  break;
	       }
	    }
	 }

	 if (f)
	    break;

      }
   }

   return f;
} 

bool
coot::restraints_container_t::check_for_1_4_relation(int idx_1, int idx_2) const {

   bool is_1_4 = false;

   for (unsigned int ii=0; ii<restraints_vec.size(); ii++) { 
      if (restraints_vec[ii].restraint_type == coot::ANGLE_RESTRAINT) {

	 if (idx_1 == restraints_vec[ii].atom_index_1 ||
	     idx_1 == restraints_vec[ii].atom_index_3 ||
	     idx_2 == restraints_vec[ii].atom_index_1 ||
	     idx_2 == restraints_vec[ii].atom_index_3) { 

	    for (unsigned int jj=ii; jj<restraints_vec.size(); jj++) {
	       if (jj != ii) { 
		  if (restraints_vec[jj].restraint_type == coot::ANGLE_RESTRAINT) {

		     if (idx_2 == restraints_vec[jj].atom_index_1 ||
			 idx_2 == restraints_vec[jj].atom_index_3 ||
			 idx_1 == restraints_vec[jj].atom_index_1 ||
			 idx_1 == restraints_vec[jj].atom_index_3) {

			if (false)
			   std::cout << "check_for_1_4_relation() indices "
				     << idx_1 << " " << idx_2
				     << " examining angle restraint pair "
				     << restraints_vec[ii].atom_index_1 << " "
				     << restraints_vec[ii].atom_index_2 << " "
				     << restraints_vec[ii].atom_index_3 << " and "
				     << restraints_vec[jj].atom_index_1 << " "
				     << restraints_vec[jj].atom_index_2 << " "
				     << restraints_vec[jj].atom_index_3 << std::endl;

			if ((restraints_vec[ii].atom_index_2 == restraints_vec[jj].atom_index_1) ||
			    (restraints_vec[ii].atom_index_2 == restraints_vec[jj].atom_index_3)) {
			   
			   if ((restraints_vec[jj].atom_index_2 == restraints_vec[ii].atom_index_1) ||
			       (restraints_vec[jj].atom_index_2 == restraints_vec[ii].atom_index_3)) {
			      
			      is_1_4 = true;
			      break;
			   }
			} 
		     }
		  }
	       }
	    }
	 }
      }
      if (is_1_4)
	 break;
   }
   // std::cout << "debug:: check_for_1_4_relation() " << idx_1 << " " << idx_2 << " is " << is_1_4 << std::endl;
   return is_1_4;
}

// check either way round.
//
// this can be static.
bool
coot::restraints_container_t::check_for_O_C_1_5_relation(mmdb::Atom *at_1, mmdb::Atom *at_2) {

   // PDBv3 FIXME.
   
   bool match = false;
   if (at_2->residue != at_1->residue) {

      // std::cout << "debug check_for_O_C_1_5_relation " << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << std::endl;

      // Check first at_1 is O(n) and at_2 is C(n+1)
      // 
      if ((at_1->GetSeqNum() + 1) == at_2->GetSeqNum()) {
	 std::string atom_name_1 = at_1->GetAtomName();
	 std::string atom_name_2 = at_2->GetAtomName();

	 if (atom_name_1 == " O  ") { 
	    if (atom_name_2 == " C  ") { 
	 
	       std::string chain_id_1 = at_1->GetChainID();
	       std::string chain_id_2 = at_2->GetChainID();
	       
	       if (chain_id_2 == chain_id_1) {
		  match = true;
	       } 
	    }
	 }
      }

      if (match) return match;

      // Check now that at_1 is C(n+1) and at_2 is O(n)
      // 
      if ((at_2->GetSeqNum() + 1) == at_1->GetSeqNum()) {
	 std::string atom_name_1 = at_1->GetAtomName();
	 std::string atom_name_2 = at_2->GetAtomName();

	 if (atom_name_1 == " C  ") { 
	    if (atom_name_2 == " O  ") { 
	 
	       std::string chain_id_1 = at_1->GetChainID();
	       std::string chain_id_2 = at_2->GetChainID();
	       
	       if (chain_id_2 == chain_id_1) {
		  match = true;
	       } 
	    }
	 }
      }
   }
   return match;
}


void
coot::restraints_container_t::symmetry_non_bonded_contacts(bool print_table) {

   int n_non_symmetric = 0;
   int n_ele = 0;
   int idx;
   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
      n_ele += filtered_non_bonded_atom_indices[i].size();
      for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) {
	 idx = filtered_non_bonded_atom_indices[i][j];
	 // is i in the idx set?
	 if (std::find(filtered_non_bonded_atom_indices[idx].begin(),
		       filtered_non_bonded_atom_indices[idx].end(),
		       i) == filtered_non_bonded_atom_indices[idx].end()) {
	    // it wasn't - i.e. non-symmetry
	    if (0) {
	       std::cout << "   " << atom_spec_t(atom[idx]) << " was an unreciprocated neighbour of "
			 << atom_spec_t(atom[i]) << std::endl;
	       std::cout << "  to  " << idx << " added " << i << std::endl;
	    }
	    int prev_size = filtered_non_bonded_atom_indices[idx].size();
	    filtered_non_bonded_atom_indices[idx].push_back(i);
	    if (0)
	       std::cout << "  filtered_non_bonded_atom_indices[" << idx << "] was of size "
			 << prev_size << " and now " << filtered_non_bonded_atom_indices[idx].size()
			 << std::endl;
	    n_non_symmetric++;
	 }
      }
   }

   if (print_table) { 
      for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
	 std::cout << "  " << i << " : ";
	 for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++)
	    std::cout << " " << filtered_non_bonded_atom_indices[i][j];
	 std::cout << "\n";
      }
   }
}


// fill the member data filtered_non_bonded_atom_indices
// 
void
coot::restraints_container_t::construct_non_bonded_contact_list(const coot::bonded_pair_container_t &bpc,
								const coot::protein_geometry &geom) {

   if (from_residue_vector)
      construct_non_bonded_contact_list_by_res_vec(bpc, geom);
   else 
      construct_non_bonded_contact_list_conventional();

}

std::pair<bool, double>
coot::simple_restraint::get_nbc_dist(const std::string &atom_1_type,
				     const std::string &atom_2_type, 
				     const protein_geometry &geom) {

   // This function reduces NBC distance for ring systems types and
   // H-B donor&acceptor combinations.  There is no provision of 1-4
   // distances and I dont know about 1-3 distances (which is a bit of
   // a worry).
   // 
   return geom.get_nbc_dist(atom_1_type, atom_2_type);
}

double
coot::simple_restraint::torsion_distortion(double model_theta) const {

   if ((restraint_type != TORSION_RESTRAINT) && (restraint_type != TRANS_PEPTIDE_RESTRAINT)) return 0;

   // this functions needs to mirror distortion_score_torsion()
   double diff = 99999.9; 
   double tdiff; 
   double trial_target; 
   int per = periodicity;
   for(int i=0; i<per; i++) { 
      // trial_target = torsion_restraint.target_value + double(i)*360.0/double(per);  ??
      trial_target = target_value + double(i)*360.0/double(per); 
      if (trial_target >= 360.0) trial_target -= 360.0; 
      tdiff = model_theta - trial_target;
      if (tdiff < -180) tdiff += 360;
      if (tdiff >  180) tdiff -= 360;
      if (fabs(tdiff) < fabs(diff)) { 
	 diff = tdiff;
      }
   }
   if (diff < -180.0) { 
      diff += 360.; 
   } else { 
      if (diff > 180.0) { 
	 diff -= 360.0; 
      }
   }
   if (false)
      std::cout << "in torsion_distortion() per: " << per << " diff: " << diff << " sigma: " << sigma
		<< "   penalty: " << diff*diff/(sigma * sigma) << std::endl;
   return diff*diff/(sigma * sigma);
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
		 ierr =  res_selection_local_inner[jat]->GetUDData(udd_atom_index_handle, 
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
			    it!=bonded_atom_indices[atom_index].end(); it++) {
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
   std::cout << "INFO:: nbc computation " << "elapsed time: " << elapsed_seconds.count() << "s\n";

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

bool
coot::restraints_container_t::is_a_moving_residue_p(mmdb::Residue *r) const {

   bool ret = false;

// Hmm! Was this even the right test?
//
//    for (unsigned int i=0; i<residues_vec.size(); i++) {
//       if (residues_vec[i].second == r) {
// 	 ret = 1;
// 	 break;
//       }
//    }
//    return ret;

   return (residues_vec_moving_set.find(r) != residues_vec_moving_set.end());

}

int
coot::restraints_container_t::get_CA_index(mmdb::Residue *residue_p) const {

   return get_atom_index(std::string(" CA "), residue_p);
}


int
coot::restraints_container_t::get_N_index(mmdb::Residue *residue_p) const {

   return get_atom_index(std::string(" N  "), residue_p);
}

int
coot::restraints_container_t::get_atom_index(const std::string &atom_name_in,
					      mmdb::Residue *residue_p) const {

   int idx = -2; // not here initally 
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      mmdb::Atom *at = residue_atoms[i];
      std::string atom_name(at->GetAtomName());
      if (atom_name == atom_name_in) { // PDBv3 FIXME
	 idx = i;
	 break;
      }
   }

   return idx;
}


// this function should only be called for residues that are standard amino acids.
coot::restraints_container_t::restraint_counts_t
coot::restraints_container_t::add_N_terminal_residue_bonds_and_angles_to_hydrogens(mmdb::Residue *residue_p) {

   restraint_counts_t rc;
   int n_bond_restraints = 0;
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   int N_index = -1; // residue-based. -2 is "checked and not here", -1 is "not check"
   int CA_index = -1;  // ditto

   // we need the map to deal with the alt-confs
   std::map<std::string, int> h1s;
   std::map<std::string, int> h2s;
   std::map<std::string, int> h3s;

   for (int i=0; i<n_residue_atoms; i++) {
      mmdb::Atom *at = residue_atoms[i];
      std::string atom_name(at->GetAtomName());
      if (atom_name == " H1 " || atom_name == " H2 " || atom_name == " H3 ") {  // PDBv3 FIXME
	 if (N_index == -1) // unset
	    N_index  = get_N_index(residue_p);
	 if (CA_index == -1)
	    CA_index = get_CA_index(residue_p);

	 // std::cout << "i " << i << " N_index " << N_index << " CA_index " << CA_index << std::endl;

	 if (N_index >= 0) {
	    if (CA_index >= 0) {
	       int atom_index_1 = -1;
	       int atom_index_2 = -1;
	       int atom_index_3 = -1;
	       int udd_get_data_status_1 = residue_atoms[i       ]->GetUDData(udd_atom_index_handle, atom_index_1);
	       int udd_get_data_status_2 = residue_atoms[N_index ]->GetUDData(udd_atom_index_handle, atom_index_2);
	       int udd_get_data_status_3 = residue_atoms[CA_index]->GetUDData(udd_atom_index_handle, atom_index_3);
	       std::vector<bool> fixed_flags_b = make_fixed_flags(atom_index_1, atom_index_2);
	       std::vector<bool> fixed_flags_a = make_fixed_flags(atom_index_1, atom_index_2, atom_index_3);
	       add(BOND_RESTRAINT,  atom_index_1, atom_index_2, fixed_flags_b, 0.86, 0.02, 1.2);
	       // std::cout << "debug:: add_bond " << atom_index_1 << " " << atom_index_2 << " "
	       // << fixed_flags_b[0] << " " << fixed_flags_b[1] << std::endl;
	       bool is_single_H_atom_angle_restraint = true;
	       add(ANGLE_RESTRAINT, atom_index_1, atom_index_2, atom_index_3, fixed_flags_a, 109.5, 2.0,
		   is_single_H_atom_angle_restraint);
	       n_bond_restraints++;
	       rc.n_bond_restraints++;
	       rc.n_angle_restraints++;
	       // bonded_atom_indices[atom_index_1].push_back(atom_index_2);
	       // bonded_atom_indices[atom_index_2].push_back(atom_index_1);
	       // bonded_atom_indices[atom_index_1].push_back(atom_index_3);
	       // bonded_atom_indices[atom_index_3].push_back(atom_index_1);
	       bonded_atom_indices[atom_index_1].insert(atom_index_2);
	       bonded_atom_indices[atom_index_2].insert(atom_index_1);
	       bonded_atom_indices[atom_index_1].insert(atom_index_3);
	       bonded_atom_indices[atom_index_3].insert(atom_index_1);
	    }
	 }

	 // PDBv3 FIXME

	 // store atoms for inter-hydrogen angle restraints
	 if (atom_name == " H1 ") {
	    int ai;
	    at->GetUDData(udd_atom_index_handle, ai);
	    h1s[at->altLoc] = ai;
	 }
	 if (atom_name == " H2 ") {
	    int ai;
	    at->GetUDData(udd_atom_index_handle, ai);
	    h2s[at->altLoc] = ai;
	 }
	 if (atom_name == " H3 ") {
	    int ai;
	    at->GetUDData(udd_atom_index_handle, ai);
	    h3s[at->altLoc] = ai;
	 }
      }
   }

   // Now do the inter-hydrogen angle restraints

   if (N_index >= 0) {
      std::map<std::string, int>::const_iterator it_1, it_2, it_3;
      for(it_1=h1s.begin(); it_1!=h1s.end(); it_1++) {
	 const std::string &key_alt_conf = it_1->first;
	 it_2 = h2s.find(key_alt_conf);
	 it_3 = h3s.find(key_alt_conf);
	 if (it_2 != h2s.end()) {
	    std::vector<bool> fixed_flags_a12 = make_fixed_flags(it_1->second, N_index, it_2->second);
	    add(ANGLE_RESTRAINT, it_1->second, N_index, it_2->second, fixed_flags_a12, 109.5, 2.0, false);
	 }

	 if (it_3 != h3s.end()) {
	    std::vector<bool> fixed_flags_a13 = make_fixed_flags(it_1->second, N_index, it_3->second);
	    add(ANGLE_RESTRAINT, it_1->second, N_index, it_3->second, fixed_flags_a13, 109.5, 2.0, false);
	 }

	 if (it_2 != h2s.end() && it_3 != h3s.end()) {
	    std::vector<bool> fixed_flags_a23 = make_fixed_flags(it_2->second, N_index, it_3->second);
	    add(ANGLE_RESTRAINT, it_2->second, N_index, it_3->second, fixed_flags_a23, 109.5, 2.0, false);
	 }
      }
   }

   return rc;
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

   for (unsigned int ib=0; ib<geom[idr].second.bond_restraint.size(); ib++) {
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name1(res_selection[iat]->name);

	 if (debug)
	    std::cout << "comparing first (pdb) :" << pdb_atom_name1
		      << ": with (dict) :"
		      << geom[idr].second.bond_restraint[ib].atom_id_1_4c()
		      << ":" << std::endl; 

	 if (pdb_atom_name1 == geom[idr].second.bond_restraint[ib].atom_id_1_4c()) {
	    for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

	       std::string pdb_atom_name2(res_selection[iat2]->name);

	       if (debug)
		  std::cout << "comparing second (pdb) :" << pdb_atom_name2
			    << ": with (dict) :"
			    << geom[idr].second.bond_restraint[ib].atom_id_2_4c()
			    << ":" << std::endl;
	       
	       if (pdb_atom_name2 == geom[idr].second.bond_restraint[ib].atom_id_2_4c()) {

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
			   }
			   if (is_hydrogen(atom[index2])) {
			      mmdb::Atom *H_at = atom[index2];
			      mmdb::Atom *parent_at = atom[index1];
			      std::string atom_name(parent_at->name);
			      std::string te = dict.type_energy(atom_name);
			      hb_t hbt = geom.get_h_bond_type(te);
			      H_atom_parent_energy_type_atom_map[H_at] = hbt;
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

int
coot::restraints_container_t::add_torsions(int idr, mmdb::PPAtom res_selection,
					   int i_no_res_atoms,
					   mmdb::PResidue SelRes,
					   const coot::protein_geometry &geom) {

   int n_torsion_restr = 0; 

   for (unsigned int ib=0; ib<geom[idr].second.torsion_restraint.size(); ib++) {

      // Joel Bard fix: Don't add torsion restraints for torsion that
      // have either s.d. or period 0

      if (geom[idr].second.torsion_restraint[ib].periodicity() > 0) { // we had this test most inner
	 if (geom[idr].second.torsion_restraint[ib].esd() > 0.000001) { // new test
	 
	    // now find the atoms
	    for (int iat=0; iat<i_no_res_atoms; iat++) {
	       std::string pdb_atom_name1(res_selection[iat]->name);

	       if (pdb_atom_name1 == geom[idr].second.torsion_restraint[ib].atom_id_1_4c()) {
		  for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

		     std::string pdb_atom_name2(res_selection[iat2]->name);
		     if (pdb_atom_name2 == geom[idr].second.torsion_restraint[ib].atom_id_2_4c()) {
				    
			// 		  std::cout << "atom match 1 " << pdb_atom_name1;
			// 		  std::cout << " atom match 2 " << pdb_atom_name2
			// 			    << std::endl;

			for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		     
			   std::string pdb_atom_name3(res_selection[iat3]->name);
			   if (pdb_atom_name3 == geom[idr].second.torsion_restraint[ib].atom_id_3_4c()) {
		  
			      for (int iat4=0; iat4<i_no_res_atoms; iat4++) {
		     
				 std::string pdb_atom_name4(res_selection[iat4]->name);
				 if (pdb_atom_name4 == geom[idr].second.torsion_restraint[ib].atom_id_4_4c()) {
		  
				    // now we need the indices of
				    // pdb_atom_name1 and
				    // pdb_atom_name2 in asc.atom_selection:

				    int index1;
				    int index2;
				    int index3;
				    int index4;

				    res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
				    res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
				    res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);
				    res_selection[iat4]->GetUDData(udd_atom_index_handle, index4);

				    double torsion_angle = geom[idr].second.torsion_restraint[ib].angle();
				    if (torsion_angle < 0)
				       torsion_angle += 360;
				    if (torsion_angle > 360)
				       torsion_angle -= 360;

				    std::vector<bool> fixed_flags =
				       make_fixed_flags(index1, index2, index3, index4);
				    add(TORSION_RESTRAINT, index1, index2, index3, index4,
					fixed_flags,
					torsion_angle,
					geom[idr].second.torsion_restraint[ib].esd(),
					1.2,  // junk value
					geom[idr].second.torsion_restraint[ib].periodicity());
				    if (0) // debug
				       std::cout << "Adding monomer torsion restraint: "
						 << index1 << " "
						 << index2 << " "
						 << index3 << " "
						 << index4 << " angle "
						 << geom[idr].second.torsion_restraint[ib].angle() << " esd " 
						 << geom[idr].second.torsion_restraint[ib].esd() << " period " 
						 << geom[idr].second.torsion_restraint[ib].periodicity()
						 << std::endl;
				    n_torsion_restr++;
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
				       restraints_vec.push_back(simple_restraint(CHIRAL_VOLUME_RESTRAINT, indexc,
										 index1, index2, index3,
										 geom[idr].second.chiral_restraint[ic].volume_sign,
										 geom[idr].second.chiral_restraint[ic].target_volume(),
										 geom[idr].second.chiral_restraint[ic].volume_sigma(),
										 fixed_flags, chiral_hydrogen_index));
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
	    std::cout << "INFO:: Adding Ramachandran restraint "
		      << "type " << std::setw(6) << zort << " for " << residue_spec_t(this_res)
		      << " " << this_res->GetResName() << " "
// 		      << coot::atom_spec_t(atom[atom_indices[0]]) << " "
// 		      << coot::atom_spec_t(atom[atom_indices[1]]) << " "
// 		      << coot::atom_spec_t(atom[atom_indices[2]]) << " "
// 		      << coot::atom_spec_t(atom[atom_indices[3]]) << " "
// 		      << coot::atom_spec_t(atom[atom_indices[4]])
		      << "fixed: "
		      << fixed_flag[0] << " " << fixed_flag[1] << " "
		      << fixed_flag[2] << " " << fixed_flag[3] << " "
		      << fixed_flag[4]
		      << std::endl;

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

coot::atom_spec_t
coot::restraints_container_t::get_atom_spec(int atom_index) const {

   if (atom)
      return atom_spec_t(atom[atom_index]);
   else
      return atom_spec_t();
}


int
coot::restraints_container_t::get_asc_index(const coot::atom_spec_t &spec) const {

   return get_asc_index_new(spec.atom_name.c_str(), spec.alt_conf.c_str(), spec.res_no,
			    spec.ins_code.c_str(), spec.chain_id.c_str());
}


int
coot::restraints_container_t::get_asc_index(const char *at_name,
					    const char *alt_loc,
					    int resno,
					    const char *ins_code,
					    const char *chain_id) const {

   return get_asc_index_new(at_name, alt_loc, resno, ins_code, chain_id);
   
}

int
coot::restraints_container_t::get_asc_index(mmdb::Atom *at) {

   int idx = -1;
   at->GetUDData(udd_atom_index_handle, idx);
   return idx;
}

int
coot::restraints_container_t::get_asc_index_new(const char *at_name,
						const char *alt_loc,
						int resno,
						const char *ins_code,
						const char *chain_id) const {

   int index = -1;

   if (mol) { 
      int SelHnd = mol->NewSelection(); // d
      mol->SelectAtoms(SelHnd,
		       0,
		       chain_id,
		       resno, ins_code,
		       resno, ins_code,
		       "*",      // resnames
		       at_name,  // anames
		       "*",      // elements
		       alt_loc  // altLocs 
		       );

      int nSelAtoms;
      mmdb::PPAtom SelAtom = NULL;
      mol->GetSelIndex(SelHnd, SelAtom, nSelAtoms);

      if (nSelAtoms > 0) {
	 if (udd_atom_index_handle >= 0) { 
	    SelAtom[0]->GetUDData(udd_atom_index_handle, index); // sets index
	 } else { 
	    index = get_asc_index_old(at_name, resno, chain_id);
	 } 
      }
      mol->DeleteSelection(SelHnd);
   }
   return index;
}

int
coot::restraints_container_t::get_asc_index_old(const std::string &at_name,
						int resno,
						const char *chain_id) const {

   int index = -1;
   int SelHnd = mol->NewSelection();
   
   mol->SelectAtoms(SelHnd,
			0,
			chain_id,
			resno, "*",
			resno, "*",
			"*", // rnames
			at_name.c_str(), // anames
			"*", // elements
			"*" // altLocs 
			);

   int nSelAtoms;
   mmdb::PPAtom SelAtom;
   mol->GetSelIndex(SelHnd, SelAtom, nSelAtoms);

   if (nSelAtoms > 0) {
      // now check indices.
      // 
      // Sigh. Is mmdb really this shit or am I missing something?
      //
      // It's not shit, you are not using it as it is supposed to be
      // used, I think. Instead of passing around the index to an
      // atom selection, you should simply be passing a pointer to an
      // atom.
      for (int i=0; i<n_atoms; i++) {
	 if (atom[i] == SelAtom[0]) {
	    index = i;
	    break;
	 }
      }
   }
   mol->DeleteSelection(SelHnd);

   if (index == -1 ) { 
      std::cout << "ERROR:: failed to find atom index for "
		<< at_name << " " << resno << " " << chain_id
		<< std::endl;
   } 
   return index; 
}

void
coot::restraints_container_t::post_add_new_restraint() {

   // adjust just one - used for target position (atom pull) restraints

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   unsigned int idx_rest = restraints_vec.size() -1;
   // if idx_rest is not in restraints_indices, add it to the back of restraints_indices.

   // If this is slow, we can go through this list backwards.
   bool found = false;
   for (std::size_t i=0; i<restraints_indices.size(); i++) {
      const std::vector<std::size_t> &v = restraints_indices[i];
      for (std::size_t j=0; j<v.size(); j++) {
	 if (v[j] == idx_rest) {
	    found = true;
	    break;
	 }
      }
      if (found) break;
   }

   if (! found)
      restraints_indices.back().push_back(idx_rest);

#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

}

void
coot::restraints_container_t::post_add_new_restraints() {

   // regenerate them all

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   make_df_restraints_indices();
#endif
}


// now setup the gsl_vector with initial values
//
// We presume that the atoms in mmdb::PPAtom are exactly the same 
// order as they are in the pdb file that refmac/libcheck uses
// to generate the restraints. 
//  
void 
coot::restraints_container_t::setup_gsl_vector_variables() {

   // recall that x is a class variable, 
   // (so are n_atoms and atom, which were set in the constructor)
   //  

   x = gsl_vector_alloc(3*n_atoms);

   // If atom is going out of date, check how atom is handled when the
   // restraints (this object) goes out of scope.  the destructor.
   
//    std::cout << "DEBUG:: using atom array pointer " << atom << std::endl;
//    std::cout << "DEBUG:: Top few atoms" << std::endl;
//    for (int i=0; i<10; i++) {
//       std::cout << "Top atom " << i << " "
// 		<< atom[i] << "   "
// 		<< atom[i]->x << " "
// 		<< atom[i]->y << " "
// 		<< atom[i]->z << " "
// 		<< std::endl;
//    } 

   for (int i=0; i<n_atoms; i++) {
      int idx = 3*i;
      gsl_vector_set(x, idx,   atom[i]->x);
      gsl_vector_set(x, idx+1, atom[i]->y);
      gsl_vector_set(x, idx+2, atom[i]->z);
   }

}



void 
coot::restraints_container_t::update_atoms(gsl_vector *s) { 

   int idx;

   if (false) { 
      std::cout << "update_atom(0): from " << atom[0]->x  << " " << atom[0]->y << " " << atom[0]->z
		<< std::endl;
      double xx = gsl_vector_get(s, 0);
      double yy = gsl_vector_get(s, 1);
      double zz = gsl_vector_get(s, 2);
      std::cout << "                  to " << xx  << " " << yy << " " << zz << std::endl;
   }
 
   if (!s) {
      std::cout << "ERROR:: in update_atoms() s has disappeared! - skip update " << std::endl;
   } else {
      for (int i=0; i<n_atoms; i++) { 
         idx = 3*i;
         atom[i]->x = gsl_vector_get(s,idx);
         atom[i]->y = gsl_vector_get(s,idx+1);
         atom[i]->z = gsl_vector_get(s,idx+2);
      }
   }
}


void
coot::restraints_container_t::position_OXT() { 

   if (oxt_reference_atom_pos.size()== 4) { 
      // std::cout << "DEBUG:: Positioning OXT by dictionary" << std::endl;
      double tors_o = 
	 clipper::Coord_orth::torsion(oxt_reference_atom_pos[0], 
				      oxt_reference_atom_pos[1], 
				      oxt_reference_atom_pos[2], 
				      oxt_reference_atom_pos[3]);
      double angl_o = clipper::Util::d2rad(120.8);
      clipper::Coord_orth oxt_pos(oxt_reference_atom_pos[0], 
				  oxt_reference_atom_pos[1], 
				  oxt_reference_atom_pos[2], 
				  1.231, angl_o, tors_o + M_PI);
      atom[oxt_index]->x = oxt_pos.x();
      atom[oxt_index]->y = oxt_pos.y();
      atom[oxt_index]->z = oxt_pos.z();
   }
} 

int
coot::restraints_container_t::write_new_atoms(std::string pdb_file_name) { 

   //
   int status = -1;
   if (mol != NULL) {
      // return 0 on success, non-zero on failure.
      status = mol->WritePDBASCII(pdb_file_name.c_str());
      if (status == 0)
	 std::cout << "INFO:: output file: " << pdb_file_name
		   << " written." << std::endl;
      else
	 std::cout << "WARNING:: output file: " << pdb_file_name
		   << " not written." << std::endl;
   } else { 
      std::cout << "not constructed from asc, not writing coords" << std::endl;
   }
   return status;
}

void
coot::restraints_container_t::info() const {

   std::cout << "INFO:: There are " << restraints_vec.size() << " restraints" << std::endl;

   for (unsigned int i=0; i< restraints_vec.size(); i++) {
      if (restraints_vec[i].restraint_type == coot::TORSION_RESTRAINT) {
	 std::cout << "INFO:: restraint " << i << " is of type "
		   << restraints_vec[i].restraint_type << std::endl;

	 std::cout << restraints_vec[i].atom_index_1 << " "
		   << restraints_vec[i].atom_index_2 << " "
		   << restraints_vec[i].atom_index_3 << " "
		   << restraints_vec[i].atom_index_4 << " "
		   << restraints_vec[i].target_value << " "
		   << restraints_vec[i].sigma << " " << std::endl
		   << " with "
  		   << restraints_vec[i].plane_atom_index.size() << " vector atoms " << std::endl
		   << " with periodicity "
  		   << restraints_vec[i].periodicity << std::endl;
      }

      std::cout << "restraint number " << i << " is restraint_type " <<
	 restraints_vec[i].restraint_type << std::endl;
   }
} 

void
coot::simple_refine(mmdb::Residue *residue_p,
		    mmdb::Manager *mol,
		    const coot::dictionary_residue_restraints_t &dict_restraints) {

   if (residue_p) {
      if (mol) {

	 int imol = 0; // shouldn't matter
	 protein_geometry geom;
	 geom.replace_monomer_restraints(residue_p->GetResName(), imol, dict_restraints);
   
	 short int have_flanking_residue_at_start = 0;
	 short int have_flanking_residue_at_end = 0;
	 short int have_disulfide_residues = 0;
	 std::string altloc("");
	 std::vector<coot::atom_spec_t> fixed_atom_specs;

	 char *chain_id = residue_p->GetChainID();
	 int istart_res = residue_p->GetSeqNum();
	 int iend_res   = istart_res;
	 clipper::Xmap<float> dummy_xmap;

	 coot::restraints_container_t restraints(istart_res,
						 iend_res,
						 have_flanking_residue_at_start,
						 have_flanking_residue_at_end,
						 have_disulfide_residues,
						 altloc,
						 chain_id,
						 mol,
						 fixed_atom_specs,
						 &dummy_xmap);
   
	 // restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
	 restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
	 pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	 bool do_internal_torsions = true;
	 bool do_trans_peptide_restraints = true;
	 restraints.make_restraints(imol, geom, flags, do_internal_torsions,
				    do_trans_peptide_restraints, 0, 0, true, true, pseudos);
	 restraints.minimize(flags, 3000, 1);
      }
   }
}

void
coot::restraints_container_t::copy_from(int i) {

}

// this was an experiment when trying to work with const ref data xmap.
// It didn't work this way.
void
coot::restraints_container_t::copy_from(const coot::restraints_container_t &rest_in) {

   restraints_vec = rest_in.restraints_vec;

   std::cout << "in copy_from we now have " << restraints_vec.size() << " restraints "<< std::endl;
   atom = rest_in.atom;
   from_residue_vector = rest_in.from_residue_vector;
   SelHnd_atom = rest_in.SelHnd_atom;
		      
   par = rest_in.par;
   n_atoms = rest_in.n_atoms;
   x = rest_in.x;
   are_all_one_atom_residues = rest_in.are_all_one_atom_residues;
   mol = rest_in.mol;
      
   residues_vec = rest_in.residues_vec;
   residues_vec_moving_set = rest_in.residues_vec_moving_set;

   udd_bond_angle = rest_in.udd_bond_angle;
   udd_atom_index_handle = rest_in.udd_atom_index_handle;

   SelResidue_active = rest_in.SelResidue_active;
   nSelResidues_active= rest_in.nSelResidues_active;

   filtered_non_bonded_atom_indices = rest_in.filtered_non_bonded_atom_indices;

   istart_res = rest_in.istart_res;
   iend_res = rest_in.iend_res;
      
   istart_minus_flag = rest_in.istart_minus_flag;
   iend_plus_flag = rest_in.iend_plus_flag;
   chain_id_save = rest_in.chain_id_save;

   previous_residue = rest_in.previous_residue;
   next_residue = rest_in.next_residue;
      
   verbose_geometry_reporting = rest_in.verbose_geometry_reporting;
   
   initial_position_params_vec = rest_in.initial_position_params_vec;

   multimin_func = rest_in.multimin_func;

   include_map_terms_flag = rest_in.include_map_terms_flag;

   lograma = rest_in.lograma;
   zo_rama = rest_in.zo_rama;
   rama_plot_weight = rest_in.rama_plot_weight;
   rama_type = rest_in.rama_type;

   map_weight = rest_in.map_weight;

   non_bonded_neighbour_residues = rest_in.non_bonded_neighbour_residues;

   have_oxt_flag = rest_in.have_oxt_flag;
   oxt_index = rest_in.oxt_index;
   residues_with_OXTs = rest_in.residues_with_OXTs;

   oxt_reference_atom_pos = rest_in.oxt_reference_atom_pos;
   do_numerical_gradients_flag = rest_in.do_numerical_gradients_flag;

   bonded_atom_indices = rest_in.bonded_atom_indices;

   // public:
   fixed_atom_indices = rest_in.fixed_atom_indices;
   restraints_usage_flag = rest_in.restraints_usage_flag;

   use_map_gradient_for_atom = rest_in.use_map_gradient_for_atom;
   atom_z_occ_weight = rest_in.atom_z_occ_weight;
   geman_mcclure_alpha = rest_in.geman_mcclure_alpha;

   cryo_em_mode = rest_in.cryo_em_mode;

#ifdef HAVE_CXX_THREAD
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

   // thread pool!
      //
   thread_pool_p = rest_in.thread_pool_p;
   n_threads = rest_in.n_threads;
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#endif // HAVE_CXX_THREAD
   
}

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
unsigned int 
coot::restraints_container_t::df_by_thread_results_size() const {
   return df_by_thread_results.size();
}
#endif

#endif // HAVE_GSL
