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

#ifdef ANALYSE_REFINEMENT_TIMING
#ifdef _MSC_VER
#include <time.h>
#else
#include <sys/time.h>
#endif
#endif // ANALYSE_REFINEMENT_TIMING


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


zo::rama_table_set coot::restraints_container_t::zo_rama;
std::atomic<bool> coot::restraints_container_t::print_lock(false);

void
coot::restraints_container_t::clear() {
   bool unlocked = false;
   while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
      unlocked = false;
   }
   restraints_vec.clear();
   init(); // resets lock, fwiw
}

void
coot::restraints_container_t::get_restraints_lock() {

   bool unlocked = false;
   while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::nanoseconds(10));
      unlocked = false;
   }
}

void
coot::restraints_container_t::release_restraints_lock() {

   restraints_lock = false;
}

// static
void
coot::restraints_container_t::get_print_lock() {

   bool unlocked = false;
   while (! print_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::microseconds(1));
      unlocked = false;
   }
}

// static
void
coot::restraints_container_t::release_print_lock() {
   print_lock = false;
}




coot::restraints_container_t::~restraints_container_t() {

   unset_fixed_during_refinement_udd();

   // 20240228-PE If the restraints_container_t is closed/finished with before the refinement
   // has terminated (i.e. continue_status == GSL_CONTINUE), then free_delete_reset()
   // doesn't get called. So call it now.
   // Fixes memory leak.
   //
   // 20240229-PE Nope. This causes a crash is test_peptide_omega()
   // free_delete_reset();

   if (from_residue_vector) {
      if (atom) {
	 // this is constructed manually.

	 // Oh we can't do this here because we copy the
	 // restraints in simple_refine_residues() and that
	 // shallow copies the atom pointer - the original
	 // restraints go out of scope and call this destructor.
	 //
	 // We need a new way to get rid of atom - c.f. the
	 // linear/conventional way?
	 //
         // 20200820-PE I don't know what simple_refine_residues() is
         // I am going to ignore the above message and delete the atoms
         // now.
         // atom needs to be a shared_ptr so that I can copy
         // restraints containers.
	 // delete [] atom;
	 // atom = NULL;

         // 20240228-PE another memory leak: atom
         // It seems to me at the moment that atom can be assigned in the constructor
         // in which case, we don't want to delete it.
         // Or it can be assigned in init_from_residue_vec(), in which case, we do
         // want to delete it (atom_array_needs_to_be_deleted_at_end is set there)

         if (atom_array_needs_to_be_deleted_at_end) {
            delete [] atom;
            atom = nullptr;
         }

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


// Used in omega distortion graph
// 
coot::restraints_container_t::restraints_container_t(atom_selection_container_t asc_in,
						     const std::string &chain_id,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {
   init();
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

   for (int ires=0; ires<nSelResidues; ires++) {
      int resno = SelResidues[ires]->GetSeqNum();
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

   init();
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

void
coot::restraints_container_t::link_restraints_counts::report() const {

   // std::cout << "   Made " << n_link_bond_restr    << " " << link_type << " bond restraints\n";
   // std::cout << "   Made " << n_link_angle_restr   << " " << link_type << " angle restraints\n";
   // std::cout << "   Made " << n_link_plane_restr   << " " << link_type << " plane restraints\n";
   // std::cout << "   Made " << n_link_trans_peptide << " " << link_type << " trans-peptide restraints\n";

   logger.log(log_t::INFO, {"created", n_link_bond_restr,    link_type, "bond restraints"});
   logger.log(log_t::INFO, {"created", n_link_angle_restr,   link_type, "angle restraints"});
   logger.log(log_t::INFO, {"created", n_link_plane_restr,   link_type, "plane restraints"});
   logger.log(log_t::INFO, {"created", n_link_trans_peptide, link_type, "trans-peptide restraints"});
}

void
coot::restraints_container_t::restraint_counts_t::report(bool do_residue_internal_torsions) const {

   // std::cout << "created " << n_bond_restraints   << " bond       restraints " << std::endl;
   // std::cout << "created " << n_angle_restraints  << " angle      restraints " << std::endl;
   // std::cout << "created " << n_plane_restraints  << " plane      restraints " << std::endl;
   // std::cout << "created " << n_chiral_restr      << " chiral vol restraints " << std::endl;
   // std::cout << "created " << n_improper_dihedral_restr << " improper dihedral restraints " << std::endl;
   // if (do_residue_internal_torsions)
   //    std::cout << "created " << n_torsion_restr << " torsion restraints " << std::endl;

   logger.log(log_t::INFO, {"created", n_bond_restraints,  "bond       restraints"});
   logger.log(log_t::INFO, {"created", n_angle_restraints, "angle      restraints"});
   logger.log(log_t::INFO, {"created", n_plane_restraints, "plane      restraints"});
   logger.log(log_t::INFO, {"created", n_chiral_restr,     "chiral vol restraints"});
   logger.log(log_t::INFO, {"created", n_improper_dihedral_restr, "improper dihedral restraints"});
   if (do_residue_internal_torsions)
      logger.log(log_t::INFO, {"created", n_torsion_restr, "torsion restraints"});

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
//
// currently links are ignored.
coot::restraints_container_t::restraints_container_t(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						     const std::vector<mmdb::Link> &links_in,
						     const coot::protein_geometry &geom,
						     mmdb::Manager *mol_in,
						     const std::vector<atom_spec_t> &fixed_atom_specs,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {

   istart_minus_flag = false; // used in make_flanking_atoms_rama_restraints
   iend_plus_flag = false;

   init();
   from_residue_vector = 1;
   are_all_one_atom_residues = false;

   // int n_resiudes_in_mol = util::number_of_residues_in_molecule(mol_in);
   // std::cout << "################## n_resiudes_in_mol " << n_resiudes_in_mol << std::endl;

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
	 std::cout << "DEBUG:: restraints_container_t() constructor: passed-residue-vec residue: "
                   << residue_spec_t(residues_local[i].second)
		   << " has index " << residues_local[i].second->index << std::endl;

   residues_vec = residues_local;
   init_from_residue_vec(residues_local, geom, mol_in, fixed_atom_specs);
   fill_links(mol_in);

}

coot::restraints_container_t::restraints_container_t(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						     const coot::protein_geometry &geom,
						     mmdb::Manager *mol_in,
						     const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {

   init();
   from_residue_vector = 1;
   are_all_one_atom_residues = false;

   std::vector<atom_spec_t> dummy_fixed_atom_specs;
   std::vector<std::pair<bool,mmdb::Residue *> > residues_local;
   residues_local.reserve(residues.size());

   for(unsigned int i=0; i<residues.size(); i++)
      if (residues[i].second)
         residues_local.push_back(residues[i]);

   std::sort(residues_local.begin(), residues_local.end(), residue_sorter);
   residues_vec = residues_local;

   if (false) {
      std::cout << "debug:: in restraints_container_t() constructor with local residue size " << residues_local.size() << std::endl;
      for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
         if (residues_vec[ir].second) {
           // std::cout << "INFO:: starting init_from_residue_vec() residue " << residues_vec[ir].second << std::endl;
           logger.log(log_t::INFO, "starting init_from_residue_vec() residue", residues_vec[ir].second);
         } else {
           std::cout << "ERROR:: starting init_from_residue_vec() NUll residue " << ir << std::endl;
         }
      }
   }


   //std::cout << "##################################### constructor D  "

   init_from_residue_vec(residues_local, geom, mol_in, dummy_fixed_atom_specs); // sets mol
   fill_links(mol_in);

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
      for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); ++it)
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
   lennard_jones_epsilon = 1.0; // close enough to 0.997 kJ/mol
   cryo_em_mode = true;
   n_times_called = 0;
   n_small_cycles_accumulator = 0;
   m_s = 0;
   x = 0;
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   n_threads = 0;
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   log_cosh_target_distance_scale_factor = 3000.0;
   convert_plane_restraints_to_improper_dihedral_restraints_flag = false; // as it was in 2019.

   use_proportional_editing = false;
   pull_restraint_neighbour_displacement_max_radius = 10.0; // make this a member of the class

   init_neutron_occupancies();
}

double
coot::restraints_container_t::get_distortion_score() const {

   // Yummy mixing C and C++ APIs...
   return distortion_score(x, const_cast<void *>(static_cast<const void *>(this)));
}

void
coot::restraints_container_t::set_use_proportional_editing(bool state) {
   use_proportional_editing = state;
}


void
coot::restraints_container_t::set_has_hydrogen_atoms_state() {

   // in init model_has_hydrogens is set to true;

   bool found = false;
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      if (is_hydrogen(at)) {
         found = true;
         break;
      }
   }
   if (! found)
      model_has_hydrogen_atoms = false;

}

// pass the formal charge also?
double
coot::restraints_container_t::neutron_occupancy(const std::string &element, int formal_charge) const {

   std::string mod_ele = coot::util::remove_whitespace(element);
   if (mod_ele.length() > 1)
      mod_ele = coot::util::capitalise(mod_ele);
   if (formal_charge != 0)
      mod_ele += coot::util::int_to_string(formal_charge);

   std::map<std::string, double>::const_iterator it = neutron_occupancy_map.find(mod_ele);
   if (it != neutron_occupancy_map.end())
      return it->second;
   else
      return 0.0;

}

void
coot::restraints_container_t::set_z_occ_weights() {

   // z weights:
   //
   atom_z_occ_weight.resize(n_atoms);
   std::vector<std::pair<std::string, int> > atom_list = coot::util::atomic_number_atom_list();
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      if (! at->isTer()) {
         std::string element = at->element;
	 double z = coot::util::atomic_number(at->element, atom_list);
	 double weight = 1.0;
	 double occupancy = atom[i]->occupancy;
	 if (occupancy > 1.0) occupancy = 1.0;
         if (do_neutron_refinement) {
            int formal_charge = 0;
            occupancy = neutron_occupancy(element, formal_charge);
         }
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


}


void
coot::restraints_container_t::init_shared_post(const std::vector<atom_spec_t> &fixed_atom_specs) {

   // std::cout << "##################################### init_shared_post() " << fixed_atom_specs.size() << std::endl;

   bonded_atom_indices.resize(n_atoms);

   set_has_hydrogen_atoms_state();

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
	    if (! is_hydrogen(atom[i]) || do_hydrogen_atom_refinement)
	       use_map_gradient_for_atom[i] = true;
	 } else {
	    // std::cout << "blanking out density for atom " << i << std::endl;
	    use_map_gradient_for_atom[i] = false;
	 }
      }
   }

   set_z_occ_weights();

   // the fixed atoms:   
   // 
   assign_fixed_atom_indices(fixed_atom_specs); // convert from std::vector<atom_spec_t>
   				                // to std::vector<int> fixed_atom_indices;

   // blank out those atoms from seeing electron density map gradients

   std::set<int>::const_iterator it;
   for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); ++it)
      use_map_gradient_for_atom[*it] = false;

   if (verbose_geometry_reporting == VERBOSE)
      for (int i=0; i<n_atoms; i++)
	 std::cout << atom[i]->name << " " << atom[i]->residue->seqNum << " "
		   << use_map_gradient_for_atom[i] << std::endl;

}

// uses fixed_atom_indices
void
coot::restraints_container_t::set_fixed_during_refinement_udd() {

   if (! mol) {
      std::cout << "ERROR:: in set_fixed_during_refinement_udd() mol is null" << std::endl;
      return;
   }
   int uddHnd = mol->RegisterUDInteger(mmdb::UDR_ATOM , "FixedDuringRefinement");
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];

      //std::cout << "  setting fixed udd flag on atom " << atom_spec_t(at) << " residue "
      // << at->residue << " atom " << at << " mol " << mol << std::endl;
      // if (std::find(fixed_atom_indices.begin(), fixed_atom_indices.end(), i) == fixed_atom_indices.end())

      if (fixed_atom_indices.find(i) == fixed_atom_indices.end())
	 at->PutUDData(uddHnd, 0);
      else
	 at->PutUDData(uddHnd, 1);
   }
}

void
coot::restraints_container_t::unset_fixed_during_refinement_udd() {

   if (! mol) {
      // if the mol has been deleted and reset before destruction of a restraints_container_t
      // then it is OK for mol to be null.
      // std::cout << "ERROR:: in unset_fixed_during_refinement_udd() mol is null" << std::endl;
      return;
   }
   int uddHnd = mol->GetUDDHandle(mmdb::UDR_ATOM , "FixedDuringRefinement");
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      at->PutUDData(uddHnd, 0);
   }
}


void
coot::restraints_container_t::init_from_residue_vec(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						    const coot::protein_geometry &geom,
						    mmdb::Manager *mol_in,
						    const std::vector<atom_spec_t> &fixed_atom_specs) {


   // std::cout << "##################################### start of init_from_residue_vec() "
   //              << fixed_atom_specs.size() << std::endl;
   
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
           // std::cout << "INFO:: starting init_from_residue_vec() residue " << residues_vec[ir].second << std::endl;
           logger.log(log_t::INFO, "starting init_from_residue_vec() residue", residues_vec[ir].second);
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

   // dist_crit = 8.0

   // fill neighbour_set from rnr (excluding residues of residues_vec):
   
   std::vector<std::pair<bool,mmdb::Residue *> > active_residues_vec;
   active_residues_vec.reserve(residues_vec.size()/2 + 1);
   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
      if (!residues_vec[ir].first)
         active_residues_vec.push_back(residues_vec[ir]);
   }
   
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnr = residues_near_residues(active_residues_vec, mol, dist_crit);
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > neighbour_set;
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it_map;

   for(it_map=rnr.begin(); it_map!=rnr.end(); ++it_map) {
      mmdb::Residue *r = it_map->first;
      const std::set<mmdb::Residue *> &s = it_map->second;
      std::set<mmdb::Residue *>::const_iterator it_set;
      for (it_set=s.begin(); it_set!=s.end(); ++it_set) {
	 bool found = false;
	 for (std::size_t i=0; i<residues_vec.size(); i++) {
	    if (*it_set == residues_vec[i].second) {
               if (! residues_vec[i].first) { // moving residue
                  found = true;
                  break;
               }
	    }
	 }
	 if (! found) {
	    neighbour_set[r].insert(*it_set);
	 }
      }
   }


   bonded_pair_container_t bpc = bonded_flanking_residues_by_residue_vector(neighbour_set, geom);

   // std::cout << "bonded_pair_container_t bpc for flanking: " << std::endl;
   // std::cout << bpc << std::endl;

   // internal variable non_bonded_neighbour_residues is set by this
   // function:
   set_non_bonded_neighbour_residues_by_residue_vector(neighbour_set, bpc, geom);

   if (false) { // debug

      std::cout << "############## rnr map set: " << std::endl;
      for(it_map=rnr.begin(); it_map!=rnr.end(); ++it_map) {
	 const std::set<mmdb::Residue *> &s = it_map->second;
	 std::cout << "::: rnr " << residue_spec_t(it_map->first) << std::endl;
	 std::set<mmdb::Residue *>::const_iterator it_set;
	 for (it_set=s.begin(); it_set!=s.end(); ++it_set) {
	    std::cout << "     rnr " << residue_spec_t(*it_set) << std::endl;
	 }
      }

      std::cout << "############## neighbour set: " << std::endl;
      for(it_map=neighbour_set.begin(); it_map!=neighbour_set.end(); ++it_map) {
	 const std::set<mmdb::Residue *> &s = it_map->second;
	 std::cout << "::: " << residue_spec_t(it_map->first) << std::endl;
	 std::set<mmdb::Residue *>::const_iterator it_set;
	 for (it_set=s.begin(); it_set!=s.end(); ++it_set) {
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
   for (unsigned int i=0; i<residues.size(); i++) {
      n_atoms_limit_for_nbc += residues[i].second->GetNumberOfAtoms();
      // std::cout << "Here in init_from_residue_vec() with n_atoms_limit_for_nbc " << n_atoms_limit_for_nbc << std::endl;
   }

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

   if (debug) {
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
   atom_array_needs_to_be_deleted_at_end = true;
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
      for(it_map=rnr.begin(); it_map!=rnr.end(); ++it_map) {
	 mmdb::Residue *r = it_map->first;
	 std::cout << "###### debugging rnr: Residue " << coot::residue_spec_t(r) << std::endl;
	 const std::set<mmdb::Residue *> &s = it_map->second;
	 std::set<mmdb::Residue *>::const_iterator it_set;
	 for (it_set=s.begin(); it_set!=s.end(); ++it_set) {
	    mmdb::Residue *residue_neighb = *it_set;
	    std::cout << "###### debugging rnr:    Neighb: " << coot::residue_spec_t(residue_neighb)
		      << " " << residue_neighb << std::endl;
	 }
      }
   }

   for(it_map=rnr.begin(); it_map!=rnr.end(); ++it_map) {
      mmdb::Residue *r = it_map->first;
      const std::set<mmdb::Residue *> &s = it_map->second;
      std::set<mmdb::Residue *>::const_iterator it_set;
      for (it_set=s.begin(); it_set!=s.end(); ++it_set) {
	 mmdb::Residue *residue_neighb = *it_set;
	 // if residue_neigh is fixed and r is not then add residue_neigh as a neighbour of r
	 if (residues_vec_moving_set.find(residue_neighb) == residues_vec_moving_set.end()) {
	    if (residues_vec_moving_set.find(r) != residues_vec_moving_set.end()) {
	       fixed_neighbours_set[r].insert(residue_neighb);
	    }
	 }
      }
   }


   // std::cout << "##################################### init_from_residue_vec() calling init_shared_post() B  "
   // << fixed_atom_specs.size() << std::endl;

   init_shared_post(fixed_atom_specs); // use n_atoms, fills fixed_atom_indices

   if (false) { //debugging
      std::cout << "---- after init_shared_post(): here are the "<< fixed_atom_indices.size()
		<< " fixed atoms " << std::endl;
      std::set<int>::const_iterator it_fixed;
      for (it_fixed=fixed_atom_indices.begin(); it_fixed!=fixed_atom_indices.end(); ++it_fixed)
	 std::cout << "    " << atom_spec_t(atom[*it_fixed]) << std::endl;
   }

   add_fixed_atoms_from_flanking_residues(bpc);

   if (false) { // debug
      std::cout << "in init_from_residue_vec() here are the non_bonded_neighbour_residues " << non_bonded_neighbour_residues.size() << std::endl;
      for (unsigned int i=0; i<non_bonded_neighbour_residues.size(); i++) {
         std::cout << "    " << residue_spec_t(non_bonded_neighbour_residues[i]) << std::endl;
      }
   }

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

   // std::cout << "##################################### assign_fixed_atom_indices() " << fixed_atom_specs.size() << std::endl;

   fixed_atom_indices.clear();

   for (unsigned int i=0; i<fixed_atom_specs.size(); i++) {
      for (int iat=0; iat<n_atoms; iat++) {
         if (fixed_atom_indices.find(iat) == fixed_atom_indices.end()) { // hopefully a speed-up not a slow-down
            if (fixed_atom_specs[i].matches_spec(atom[iat])) {
               fixed_atom_indices.insert(iat); // hello grep: in assign_fixed_atom_indices()
               break;
            }
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
      if (fixed_atom_indices.find(iat) != fixed_atom_indices.end()) is_fixed = true;
      std::cout << std::setw(3) << iat << " " << atom_spec_t(atom[iat]) << "  "
		<< std::right << std::setw(10) << std::fixed << std::setprecision(3) << atom[iat]->x << " "
		<< std::right << std::setw(10) << std::fixed << std::setprecision(3) << atom[iat]->y << " "
		<< std::right << std::setw(10) << std::fixed << std::setprecision(3) << atom[iat]->z
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
		  // std::cout << "INFO:: Pre-Sanitize Minimum found (iteration number " << iter << ") at ";
		  // std::cout << m_s->f << "\n";
		  logger.log(log_t::INFO, "Pre-Sanitize Minimum found (iteration number", iter, ") at", m_s->f);
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
   for (int i=0; i<size(); i++) {
      const simple_restraint &rest = restraints_vec[i];
      if (rest.restraint_type == TARGET_POS_RESTRAINT)
	 n++;
   }

   return n;
}


// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
//
// n_steps_max = 1000 default arg
//
coot::refinement_results_t
coot::restraints_container_t::minimize(restraint_usage_Flags usage_flags, int n_steps_max) {

   short int print_chi_sq_flag = 1;
   refinement_results_t rr = minimize(usage_flags, n_steps_max, print_chi_sq_flag);
   // std::cout << "debug:: minimize() returns " << rr.progress << std::endl;
   return rr;

}


#include <gsl/gsl_blas.h> // for debugging norm of gradient

// this does a reset/setup
//
void
coot::restraints_container_t::setup_minimize() {

   if (m_s) {
      gsl_multimin_fdfminimizer_free(m_s);
      m_s = nullptr;
   }
   if (x) {
      gsl_vector_free(x);
      x = nullptr;
   }

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
   // T = gsl_multimin_fdfminimizer_vector_bfgs2;

   // setup before the first minimize (n_times_called == 1)

   setup_gsl_vector_variables();  //initial positions
   setup_multimin_func(); // provide functions for f, df, fdf

   m_s = gsl_multimin_fdfminimizer_alloc(T, n_variables());

   double step_size_multiplier = 1.0;
   // std::cout << "setting step_size_multiplier with n_atoms " << n_atoms << std::endl;

   // this is a bit "heuristic" - actually I want the number of non-fixed atoms.
   if (n_atoms < 100)
      step_size_multiplier = 0.1;

   m_initial_step_size = step_size_multiplier * gsl_blas_dnrm2(x); // how about just 0.1?

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

   // std::cout << "------------ ::minimize() basic " << size() << std::endl;
   n_times_called++;
   if (n_times_called == 1 || needs_reset)
      setup_minimize();

   refinement_results_t rr = minimize_inner(usage_flags, nsteps_max);

   return rr;
}

// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
//
coot::refinement_results_t
coot::restraints_container_t::minimize(int imol, restraint_usage_Flags usage_flags, int nsteps_max,
				       short int print_initial_chi_sq_flag,
				       const coot::protein_geometry &geom) {


   unsigned int n_steps_per_recalc_nbcs = 300000; // Hmm.
   // n_steps_per_relcalc_nbcs *= 10;

   n_times_called++;
   n_small_cycles_accumulator += nsteps_max;

   // std::cout << "------------ ::minimize() " << size() << " " << n_small_cycles_accumulator << std::endl;

   if (n_times_called == 1 || needs_reset)
      setup_minimize();

#if 1

   // make_non_bonded_contact_restraints_ng() works using the residues of the
   // input molecule - and they are not going to change as the atoms in the atom
   // vector move around.  Ideally, we need a function that takes the atoms in atom and
   // finds neighbours in mol - which is extremely not how refinement works in Coot,
   // because we find the residues to refine before refinement.
   //
   // It will take quite a rewrite to fix this - pass the moving residues and mol
   // and let a function in this class work out what the "residues near residues"
   // are. This is quite doable, but not for now - updating NBCs is not as
   // important as OpenGLv3 (shader) graphics.

   // So, for now we can update non-bonded contacts between that atoms that
   // we have for refinement (including the non-moving atoms).
   // That should do a lot of work for us in domain-refine and coot-cuda-refine
   // but not in wonky-N-terminus refine.

   if (n_small_cycles_accumulator >= n_steps_per_recalc_nbcs) {
      auto tp_0 = std::chrono::high_resolution_clock::now();
      unsigned int n_new = make_non_bonded_contact_restraints_ng(imol, geom);
      auto tp_1 = std::chrono::high_resolution_clock::now();
      // setup_minimize(); // needed? (seems not)
      if (false) {
         auto tp_2 = std::chrono::high_resolution_clock::now();
         auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
         auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
         unsigned int n_restraints = size();
         std::cout << "minimize() nbc updated made " << n_new << " (new) nbc restraints: " << n_restraints
                   << " total - timings: " << d10 << " " << d21 << " milliseconds" << std::endl;
      }
      n_small_cycles_accumulator = 0;
   }
#endif

   refinement_results_t rr = minimize_inner(usage_flags, nsteps_max);

   // std::cout << "debug:: minimize() returns " << rr.progress << std::endl;

   return rr;
}

// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
//
coot::refinement_results_t
coot::restraints_container_t::minimize_inner(restraint_usage_Flags usage_flags,
					     int nsteps_max) {

   // Was just checking that minimize() was not being called multiple times concurrently
   // (it wasn't)
   // std::cout << "DEBUG:: incrementing n_refiners_refining to " << n_refiners_refining+1 << std::endl;
   // n_refiners_refining++;

   // check that we have restraints before we start to minimize:
   if (size() == 0) {
      if (restraints_usage_flag != NO_GEOMETRY_RESTRAINTS) {
	 std::cout << "SPECIFICATION ERROR:  There are no restraints. ";
	 std::cout << "No minimization will happen" << std::endl;
	 return refinement_results_t(0, 0, "No Restraints!");
      }
   }

   if (do_numerical_gradients_flag) {
      std::cout << "debug:: minimize_inner called with usage_flags " << usage_flags << std::endl;
      debug_atoms();
   }

   // debug_atoms(); // for diagnosing refinement (atom index) issues // DEBUG-ATOMS

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

	 bool unlocked = false;
	 // while (! restraints_lock.compare_exchange_weak(unlocked, true) && !unlocked) {
	 while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
	    std::this_thread::sleep_for(std::chrono::microseconds(10));
	    unlocked = false;
	 }

         if (m_s == 0) {
            std::cout << "ERROR:: !! m_s has disappeared! " << std::endl;
            break;
         } else {
	    status = gsl_multimin_fdfminimizer_iterate(m_s);
         }

	 // this might be useful for debugging rama restraints

	 conjugate_pr_state_t *state = static_cast<conjugate_pr_state_t *>(m_s->state);
	 double pnorm  = state->pnorm;
	 double g0norm = state->g0norm;
	 //
	 if (false)
	    std::cout << "iter: " << iter << " f " << m_s->f << " " << gsl_multimin_fdfminimizer_minimum(m_s)
		      << " pnorm " << pnorm << " g0norm " << g0norm << std::endl;

	 if (status != GSL_SUCCESS) {

            if (verbose_geometry_reporting != QUIET)
               std::cout << "Unexpected error from gsl_multimin_fdfminimizer_iterate at iter " << iter << std::endl;
	    if (status == GSL_ENOPROG) {
               if (verbose_geometry_reporting != QUIET)
                  std::cout << "Error:: in gsl_multimin_fdfminimizer_iterate() result was GSL_ENOPROG" << std::endl; 
	       if (verbose_geometry_reporting != QUIET)
		  std::cout << "Error:: iter: " << iter << " f " << m_s->f << " "
			    << gsl_multimin_fdfminimizer_minimum(m_s)
			    << " pnorm " << pnorm << " g0norm " << g0norm << "\n";

               if (verbose_geometry_reporting != QUIET) {
                  // write out gradients here - with numerical gradients for comparison
                  lights_vec = chi_squareds("Final Estimated RMS Z Scores (ENOPROG)", m_s->x);
                  analyze_for_bad_restraints();
               }

	       done_final_chi_squares = true;
	       refinement_lights_info_t::the_worst_t worst_of_all = find_the_worst(lights_vec);
	       if (worst_of_all.is_set) {
		  const simple_restraint &baddie_restraint = restraints_vec[worst_of_all.restraints_index];
                  if (verbose_geometry_reporting != QUIET)
                     // std::cout << "INFO:: Most dissatisfied restraint (refine no-progress): "
                     //          << baddie_restraint.format(atom, worst_of_all.value) << std::endl;
                     logger.log(log_t::INFO, "Most dissatisfied restraint (refine no-progress):", baddie_restraint.format(atom, worst_of_all.value));
	       } else {
                  if (verbose_geometry_reporting != QUIET)
                     // std::cout << "INFO:: somehow the worst restraint was not set (no-progress)" << std::endl;
                     logger.log(log_t::INFO, "somehow the worst restraint was not set (no-progress)");
	       }

               // debugging/analysis
               if (verbose_geometry_reporting != QUIET)
                  std::cout << "----------------------- FAIL, ENOPROG --------------- " << std::endl;

               // follwing is useful - but not for everyone
               // gsl_vector *non_const_v = const_cast<gsl_vector *> (m_s->x); // because there we use gls_vector_set()
               // void *params = static_cast<void *>(this);
               // numerical_gradients(non_const_v, params, m_s->gradient, "failed-gradients.tab");
            }
            restraints_lock = false;
	    break;
	 }

         // std::cout << "Debug:: before gsl_multimin_test_gradient, status is " << status << std::endl;

         if (status == GSL_SUCCESS || status == GSL_CONTINUE) // probably just GSL_SUCCESS is what I want
            status = gsl_multimin_test_gradient(m_s->gradient, m_grad_lim);

         // std::cout << "Debug:: after gsl_multimin_test_gradient, status is " << status << std::endl;

	 if (status == GSL_SUCCESS) {
	    if (verbose_geometry_reporting != QUIET) {
	       // std::cout << "Minimum found (iteration number " << iter << ") at ";
	       // std::cout << m_s->f << "\n";
               logger.log(log_t::INFO, "Minimum found at iteration number", iter, "at", m_s->f);
            }
	    std::string title = "Final Estimated RMS Z Scores:";
	    std::vector<coot::refinement_lights_info_t> results = chi_squareds(title, m_s->x);
	    lights_vec = results;
            if (verbose_geometry_reporting != QUIET) {
               // std::cout << "-------- Results ---------" << std::endl; // should this go into analyze_for_bad_restraints()?
               logger.log(log_t::INFO, "-------- Results ---------");
               update_atoms(m_s->x); // needed for simple_refine() (maybe other times too to catch the last round)
               analyze_for_bad_restraints();
            }
	    done_final_chi_squares = true;
	 }

	 if (verbose_geometry_reporting == VERBOSE)
	    std::cout << "iteration number " << iter << " " << m_s->f << std::endl;

	 restraints_lock = false; // unlock

      }
   while ((status == GSL_CONTINUE) && (iter < nsteps_max));

   if (! done_final_chi_squares) {
      if (status != GSL_CONTINUE) {

         analyze_for_bad_restraints();
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
   // std::cout << "After update atoms " << std::endl;

   // (we don't get here unless restraints were found)
   bool found_restraints_flag = true;
   refinement_results_t rr(found_restraints_flag, status, lights_vec);
   if (refinement_results_add_details) {
      // this may be slowing thing down for big molecules. Needs investigation.
      // std::cout << "mimize_inner() calling add_details_to_refinement_results() " << std::endl;
      add_details_to_refinement_results(&rr);
   }

   // std::cout << "After rr" << std::endl;

   if (status != GSL_CONTINUE) {

      // protection so that clearing of the vectors doesn't coincide with geometric_distortions()
      // evaluation

      bool unlocked = false;
      while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
         // std::cout << "debug:: in minimize_inner() waiting for restraints_lock" << std::endl;
	 std::this_thread::sleep_for(std::chrono::microseconds(10));
	 unlocked = false;
      }

      free_delete_reset();

      // std::cout << "debug:: unlocking restraints in minimize_inner()"  << std::endl;
      restraints_lock = false; // unlock
   }

   // the bottom line from the timing test is the only thing that matters
   // is the time spend in the core minimization iterations

   if (false)
      std::cout << "-------------- Finally returning from minimize_inner() with rr with status "
                << rr.progress << std::endl;

   n_refiners_refining--;
   return rr;
}

coot::refinement_results_for_rama_t::refinement_results_for_rama_t(mmdb::Atom *at_1,
                                                                   mmdb::Atom *at_2,
                                                                   mmdb::Atom *at_3,
                                                                   mmdb::Atom *at_4,
                                                                   mmdb::Atom *at_5,
                                                                   float distortion_in) {
   distortion = distortion_in;
   atom_spec_CA = atom_spec_t(at_3);
   ball_pos_x = 0; ball_pos_y = 0; ball_pos_z = 0;
   if (at_3) {
      ball_pos_x = at_3->x + 0.5;
      ball_pos_y = at_3->y;
      ball_pos_z = at_3->z;
   }
   if (at_1 && at_2 && at_3 && at_4 && at_5) {
      clipper::Coord_orth p2 = co(at_2);
      clipper::Coord_orth p3 = co(at_3);
      clipper::Coord_orth p4 = co(at_4);
      clipper::Coord_orth v1(p3 - p2);
      clipper::Coord_orth v2(p3 - p4);
      clipper::Coord_orth v3(p4 - p2);
      clipper::Coord_orth v1_uv(v1.unit());
      clipper::Coord_orth v2_uv(v2.unit());
      clipper::Coord_orth v3_uv(v3.unit());
      clipper::Coord_orth v4(clipper::Coord_orth::cross(v2_uv, v1_uv));
      clipper::Coord_orth p2p24_mid_point(0.5 * (p4+p2));
      clipper::Coord_orth mid_point_to_CA(p3 - p2p24_mid_point);
      clipper::Coord_orth delta = 0.2 * mid_point_to_CA + 0.4 * v4;
      ball_pos_x = delta.x() + at_3->x;
      ball_pos_y = delta.y() + at_3->y;
      ball_pos_z = delta.z() + at_3->z;
   }
}

void
coot::restraints_container_t::free_delete_reset()  {

   if (false)
      std::cout << "DEBUG:: ---- free/delete/reset m_s and x" << std::endl;
   if (m_s)
      gsl_multimin_fdfminimizer_free(m_s);
   gsl_vector_free(x);
   m_s = 0;
   x = 0;
   needs_reset = true;
}

coot::refinement_results_t
coot::restraints_container_t::get_refinement_results() {

   bool found_restraints_flag = true;
   int status = GSL_SUCCESS;
   setup_minimize();
   std::vector<coot::refinement_lights_info_t> lights_vec =
      chi_squareds("Final Estimated RMS Z Scores", m_s->x);
   refinement_results_t rr(found_restraints_flag, status, lights_vec);
   add_details_to_refinement_results(&rr);
   free_delete_reset();
   return rr;
}


void
coot::restraints_container_t::add_details_to_refinement_results(coot::refinement_results_t *rr) const {

   auto tp_1 = std::chrono::high_resolution_clock::now();
   int n_restraints = size();

   // close in terms of seqnum that is
   auto is_close_main_chain_nbc = [] (const simple_restraint &restraint, mmdb::PPAtom atom) {

                                     std::vector<std::string> main_chain_atom_names = { " CA ", " C  ", " N  ", " H  ", " O  ", " HA" };
                                     mmdb::Atom *at_1 = atom[restraint.atom_index_1];
                                     mmdb::Atom *at_2 = atom[restraint.atom_index_2];
                                     int seq_num_1 = at_1->GetSeqNum();
                                     int seq_num_2 = at_2->GetSeqNum();
                                     if (abs(seq_num_1 - seq_num_2) < 2) {
                                        if (at_1->GetChain() == at_2->GetChain()) {
                                           std::string atom_name_1(at_1->GetAtomName());
                                           std::string atom_name_2(at_2->GetAtomName());
                                           if (std::find(main_chain_atom_names.begin(),
                                                         main_chain_atom_names.end(), atom_name_1) != main_chain_atom_names.end() ||
                                               std::find(main_chain_atom_names.begin(),
                                                         main_chain_atom_names.end(), atom_name_2) != main_chain_atom_names.end())
                                              return true;
                                        }
                                     }
                                     return false;
                                  };

   class nbc_baddie_atom_index_pair_t {
   public:
      int index_1;
      int index_2;
      float distortion;
      nbc_baddie_atom_index_pair_t(const simple_restraint &restraint, float dist_in) : distortion(dist_in) {
         index_1 = restraint.atom_index_1;
         index_2 = restraint.atom_index_2;
      }
   };
   std::vector<nbc_baddie_atom_index_pair_t> nbc_baddie_atom_index_pair_vec;

   std::map<int, float> nbc_baddies; // atom index to badness-score
   std::map<int, float> rama_baddies;
   unsigned int n_non_bonded_restraints = 0;
   unsigned int n_rama_restraints = 0;
   double nbc_distortion_score_sum = 0;
   double rama_distortion_score_sum = 0;

   std::vector<refinement_results_for_chiral_t> chiral_baddies;

   if (! m_s) {
      std::cout << "m_s is null - returning early from add_details_to_refinement_results() " << std::endl;
      return;
   }
   const gsl_vector *v = m_s->x;
   std::vector<refinement_results_for_rama_t> all_ramas;
   all_ramas.reserve(100);

   rr->n_restraints = n_restraints;
   rr->nbc_baddies_atom_index_map.clear();

   unsigned int n_hydrogen_bond_restraints = 0;
   unsigned int n_geman_mcclure_distance_restraints = 0;
   for (int i=0; i<n_restraints; i++) {
      const simple_restraint &restraint = restraints_vec[i];

      if (restraint.restraint_type == coot::TARGET_POS_RESTRAINT) {
          double dist = distortion_score_target_pos(restraint, 1.0, v);
          mmdb::Atom *at = atom[restraint.atom_index_1];
          std::pair<atom_spec_t, float> p(atom_spec_t(at), dist);
          // not sorted yet
          rr->sorted_atom_pulls.push_back(p);
          rr->overall_atom_pull_score += dist;
      }

      if (restraints_usage_flag & coot::NON_BONDED_MASK) {
         if (restraint.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
            n_non_bonded_restraints++;
            double dist = distortion_score_non_bonded_contact(restraint, lennard_jones_epsilon, v);
            // std::cout << "nbc " << dist << std::endl;  Vast majority < -0.05
            const double &dist_crit = 10.0; // 20230521-PE was 5.0, too many
                                            // 20230215-PE was 1.95;
                                            // 20221127-PE was 0.55
                                            // 20220924-PE was 0.25 - but that made too many
            if (dist > dist_crit) { // 20220503-PE was 0.05 - we want to see angry diego only when the
                                    // atom are really too close

               // if this is slow, add the result of this test as a boolean as the restraint
               // is being created, is_close_main_chain_nbc_flag is part of a simple_restraint;
               //
               if (! is_close_main_chain_nbc(restraint, atom)) {
                  nbc_distortion_score_sum += dist;
                  nbc_baddies[restraint.atom_index_1] += 0.5 * dist;
                  nbc_baddies[restraint.atom_index_2] += 0.5 * dist;
                  nbc_baddie_atom_index_pair_t bip(restraint, dist);
                  nbc_baddie_atom_index_pair_vec.push_back(bip);
                  rr->nbc_baddies_atom_index_map[restraint.atom_index_1].push_back(restraint.atom_index_2);
               }
            }
         }
      }

      if (restraints_usage_flag & coot::BONDS_MASK) {
         if (restraint.restraint_type == coot::BOND_RESTRAINT) {
            if (restraint.is_hydrogen_bond) {
               // std::cout << "hyrogen bond restraint " << restraint.atom_index_1 << " " << restraint.atom_index_2 << std::endl;
               n_hydrogen_bond_restraints++;
               rr->hydrogen_bond_atom_index_vec.push_back(std::make_pair(restraint.atom_index_1, restraint.atom_index_2));
            }
         }
      }

      if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) {
         if (restraint.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
            double dist = distortion_score_chiral_volume(restraint, v);
            double chiral_volume_distortion_limit = 6.0; // c.f. limiit in dynamic-valiation.cc make_chiral_volume_buttons()
            if (dist > chiral_volume_distortion_limit) {
               mmdb:: Atom *at = atom[restraint.atom_index_centre];
               clipper::Coord_orth pos(at->x, at->y, at->z);
               refinement_results_for_chiral_t cb(atom_spec_t(at), pos, dist);
               chiral_baddies.push_back(cb);
            }
         }
      }

      if (restraint.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
         if (restraints_usage_flag & coot::GEMAN_MCCLURE_DISTANCE_MASK)
            n_geman_mcclure_distance_restraints++;
      }

      if (restraints_usage_flag & coot::RAMA_PLOT_MASK) {
         if (restraint.restraint_type == coot::RAMACHANDRAN_RESTRAINT) {
            n_rama_restraints++;
            if (rama_type == restraints_container_t::RAMA_TYPE_ZO) {
               // std::cout << "----------------------- ZO type RAMA! " << std::endl;
               double dd_raw = distortion_score_rama(restraint, v, ZO_Rama(), get_rama_plot_weight());
               double dd = dd_raw / rama_plot_weight;
               dd *= 50.0; // scale to non-ZO non weighted
               if (false) // range -13 to 0 with weight 1.4, and 100 times that with weight 140
                  std::cout << "PATH A zo-rama distortion for restraint " << i << " distortion is "
                            << dd_raw << " " << " (post-mod) " << dd << " "
                            << atom_spec_t(atom[restraint.atom_index_3])
                            << std::endl;
               rama_distortion_score_sum += dd;
               if (dd > -200.01) {
                  rama_baddies[restraint.atom_index_3] += dd;
                  // 20200412-PE-merge-complexity: A bug?
                  // mmdb::Atom *at = atom[restraint.atom_index_3];
                  // std::string res_type(at->residue->GetResName());
                  // std::cout << atom_spec_t(at) << " " << res_type << " zo-rama " << dd << std::endl;
                  rama_baddies[restraint.atom_index_1] += 0.5 * dd;
               }
               refinement_results_for_rama_t rp(atom[restraint.atom_index_1],
                                                atom[restraint.atom_index_2],
                                                atom[restraint.atom_index_3],
                                                atom[restraint.atom_index_4],
                                                atom[restraint.atom_index_5], dd);
               all_ramas.push_back(rp);
            } else {
               double w = get_rama_plot_weight();
               double dd = distortion_score_rama(restraint, v, LogRama(), w) / rama_plot_weight;
               rama_distortion_score_sum += dd;
               if (false)
                  std::cout << "PATH B rama for restraint " << i << " distortion " << dd << " "
                            << atom_spec_t(atom[restraint.atom_index_3]) << " using rama plot weight " << rama_plot_weight
                            << std::endl;

               refinement_results_for_rama_t rp(atom[restraint.atom_index_1],
                                                atom[restraint.atom_index_2],
                                                atom[restraint.atom_index_3],
                                                atom[restraint.atom_index_4],
                                                atom[restraint.atom_index_5], dd);
               all_ramas.push_back(rp);
               // this cutoff (or dd) should take account of the rama weight
               // When we are plotting coloured balls and bars, we are not interested
               // in the contribution of this rama to the target function, we want to
               // know how probable this rama value is. So needs unweighting.
               //
               // Actually, only the ZO rama scores are weighted.
               //
               if (dd > -200.0) {
                  // GLY have naturally lower probabilities densities, hence higher -logPr
                  std::string rn(atom[restraint.atom_index_3]->residue->GetResName());
                  if (rn == "GLY") {
                     dd -= 50.0; // utter guess
                  }
                  rama_baddies[restraint.atom_index_3] += dd;
               }
            }
         }
      }
   }

   // --- non-bonded contacts ---

   auto atom_to_coord_orth = [] (mmdb::Atom *at) { return clipper::Coord_orth(at->x, at->y, at->z); };

   std::map<int, float>::const_iterator it;
   unsigned int idx = 0;

   auto sorter = [] (const std::pair<int, float> &v1,
                     const std::pair<int, float> &v2) {
                    return v2.second < v1.second;
                 };

   // now sort the baddie index pairs (which are needed to find the positions for the bad NBC markers)
   //
   auto sorter_2 = [] (const nbc_baddie_atom_index_pair_t &b1, const nbc_baddie_atom_index_pair_t &b2) {
                      return b2.distortion < b1.distortion;
                   };
   std::sort(nbc_baddie_atom_index_pair_vec.begin(), nbc_baddie_atom_index_pair_vec.end(), sorter_2);

   if (false) {
      for (unsigned int i=0; i<nbc_baddie_atom_index_pair_vec.size(); i++) {
         const auto &bip = nbc_baddie_atom_index_pair_vec[i];
         std::cout << "sorted baddie index pair " << std::setw(4) << bip.index_1 << " "
                   << std::setw(4) << bip.index_2 << " " << bip.distortion << std::endl;
      }
   }

   std::vector<refinement_results_nbc_baddie_t> nbc_baddies_with_spec_vec(nbc_baddie_atom_index_pair_vec.size());

   // 20220503-PE now fill nbc_baddies_with_spec_vec using nbc_baddie_atom_index_pair_vec
   unsigned int n_baddies = nbc_baddie_atom_index_pair_vec.size();
   for (unsigned int i=0; i<n_baddies; i++) {
      const auto &bip = nbc_baddie_atom_index_pair_vec[i];
      mmdb::Atom *at_1 = atom[bip.index_1];
      mmdb::Atom *at_2 = atom[bip.index_2];
      clipper::Coord_orth p_1 = atom_to_coord_orth(at_1);
      clipper::Coord_orth p_2 = atom_to_coord_orth(at_2);
      clipper::Coord_orth p_m = 0.5 * (p_1 + p_2);
      nbc_baddies_with_spec_vec[i].atom_spec_1 = atom_spec_t(at_1);
      nbc_baddies_with_spec_vec[i].atom_spec_2 = atom_spec_t(at_2);
      nbc_baddies_with_spec_vec[i].mid_point = p_m;
      nbc_baddies_with_spec_vec[i].atom_1_pos = p_1;
      nbc_baddies_with_spec_vec[i].atom_2_pos = p_2;
      nbc_baddies_with_spec_vec[i].score = bip.distortion;
   }
   
   // std::cout << "done filling nbc baddies" << std::endl;

   rr->overall_nbc_score = nbc_distortion_score_sum;
   rr->sorted_nbc_baddies = nbc_baddies_with_spec_vec;
   rr->refinement_results_contain_overall_nbc_score = true;

   // --- rama ---

   if (n_rama_restraints > 0) {
      std::vector<std::pair<int, float> > rama_baddies_vec(rama_baddies.size());
      idx = 0;
      for (it=rama_baddies.begin(); it!=rama_baddies.end(); ++it)
         rama_baddies_vec[idx++] = std::pair<int, float>(it->first, it->second);
      std::sort(rama_baddies_vec.begin(), rama_baddies_vec.end(), sorter);
      if (rama_baddies_vec.size() > 20)
         rama_baddies_vec.resize(20);
      std::vector<std::pair<atom_spec_t, float> > rama_baddies_with_spec_vec(rama_baddies_vec.size());
      for (unsigned int i=0; i<rama_baddies_vec.size(); i++) {
         int atom_index = rama_baddies_vec[i].first;
         rama_baddies_with_spec_vec[i].first  = atom_spec_t(atom[atom_index]);
         rama_baddies_with_spec_vec[i].second = rama_baddies_vec[i].second;
         // set user data meaning "is_in_a_moving_atoms_residue"
         if (fixed_atom_indices.find(atom_index) != fixed_atom_indices.end())
            rama_baddies_with_spec_vec[i].first.int_user_data = 0;
         else
            rama_baddies_with_spec_vec[i].first.int_user_data = 1;
         if (false) // debug
            std::cout << "debug:: rama_baddies_with_spec_vec " << i << " "
                      << rama_baddies_with_spec_vec[i].first << " "
                      << rama_baddies_with_spec_vec[i].second << std::endl;
      }
      rr->refinement_results_contain_overall_rama_plot_score = true;
      rr->all_ramas = all_ramas;
      rr->sorted_rama_baddies = rama_baddies_with_spec_vec;
      rr->overall_rama_plot_score = rama_distortion_score_sum;

   }

   // --- chiral ---

   auto chiral_sorter = [] (const refinement_results_for_chiral_t &r1,
                            const refinement_results_for_chiral_t &r2) {
      return (r2.distortion < r1.distortion);
   };
   std::sort(chiral_baddies.begin(), chiral_baddies.end(), chiral_sorter);
   rr->sorted_chiral_volume_baddies = chiral_baddies;

   if (false) {
      for (unsigned int ii=0; ii<rr->sorted_chiral_volume_baddies.size(); ii++)
         std::cout << "chiral-vol-baddie " << ii
                   << " " << rr->sorted_chiral_volume_baddies[ii].atom_spec
                   << " " << rr->sorted_chiral_volume_baddies[ii].distortion << std::endl;
   }

   // --- atom pulls ---

   {
      auto sorter_ap = [] (const std::pair<atom_spec_t, float> &v1,
                           const std::pair<atom_spec_t, float> &v2) {
                          return v2.second < v1.second;
                    };
      std::sort(rr->sorted_atom_pulls.begin(), rr->sorted_atom_pulls.end(), sorter_ap);
   }
      
   if (false) {
      // ~1ms for 100 residues
      auto tp_2 = std::chrono::high_resolution_clock::now();
      auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
      std::cout << "info:: add_details_to_refinement_results(): " << d21 << " microseconds\n";
   }

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

void
coot::refinement_results_t::show() const {

   std::string progress_str = "other";
   if (progress == GSL_CONTINUE) progress_str = "GSL_CONTINUE";
   if (progress == GSL_SUCCESS)  progress_str = "GSL_SUCCESS";
   if (progress == GSL_ENOPROG)  progress_str = "GSL_NO_PROGRESS";
   if (progress == GSL_FAILURE)  progress_str = "GSL_FAILURE";
   std::cout << "Refinement Ressults: " << info_text
             << " n_restraints " << n_restraints
             << " found_restraints_flag: " << found_restraints_flag
             << " progress_status " << progress_str << std::endl;
   for (const auto &l : lights) {
      std::cout << " " << l.name << " " << l.label << " " << l.value << std::endl;
   }
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
   double dist = 0; 

   // if (v->size != restraints->init_positions_size() ) {

   for (int i=0; i<restraints->init_positions_size(); i++) { 
      double d = restraints->initial_position(i) - gsl_vector_get(v, i);
      dist += 0.01*d*d;
   }
   std::cout << "starting_structure_diff_score: " << dist << std::endl; 
   return dist;
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
	    resultant.add(delta, clipper::Coord_orth(b_uv_abs * delta));
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

   if (restraints_p->thread_pool_p) {

      double results[1024]; // we will always have less than 1024 threads

      unsigned int n_ranges = ranges.size(); // clang scan-build fix.
      for(unsigned int i=0; i<n_ranges; i++) {
         results[i] = 0.0;
	 restraints_p->thread_pool_p->push(electron_density_score_from_restraints_using_atom_index_range,
					   v, std::cref(ranges[i]), restraints_p, &(results[i]),
					   std::ref(done_count_for_restraints_sets));
      }
      while (done_count_for_restraints_sets < ranges.size())
	 std::this_thread::sleep_for(std::chrono::nanoseconds(300));

      // consolidate
      for(unsigned int i=0; i<n_ranges; i++)
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

   // auto tp_1 = std::chrono::high_resolution_clock::now();

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

               //               std::cout << "ao:" << ao.format() << std::endl; // prograam terminated before
                                                                                // sphere-refine had finished.

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

   // auto tp_2 = std::chrono::high_resolution_clock::now();
   // auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
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



#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
void coot::restraints_container_t::make_distortion_electron_density_ranges() {

   // std::vector<std::pair<unsigned int, unsigned int> > ranges =
   // atom_index_ranges(n_atoms, restraints_p->n_threads);

   unsigned int nt = n_threads;
   if (nt == 0) nt = 1;
   m_atom_index_ranges = atom_index_ranges(n_atoms, nt);

}
#endif



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
										  const coot::bonded_pair_container_t &bonded_flanking_pairs,
                                                                                  const coot::protein_geometry &geom) {

   // non_bonded_neighbour_residues becomes this:
   //
   std::vector<mmdb::Residue *> nbr; // non-bonded residues

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it_map;

   // don't iterate like this:
   // for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
   // std::vector<mmdb::Residue *> neighbours =
   // coot::residues_near_residue(residues_vec[ir].second, mol, dist_crit);

   for(it_map=neighbour_set.begin(); it_map!=neighbour_set.end(); ++it_map) {

      const std::set<mmdb::Residue *> &neighbours = it_map->second;
      std::set<mmdb::Residue *>::const_iterator it_set;

      for (it_set=neighbours.begin(); it_set!=neighbours.end(); ++it_set) {
	 mmdb::Residue *test_res = *it_set;
	 if (std::find(nbr.begin(), nbr.end(), test_res) == nbr.end()) {
	    // not already there...
	    bool found = false;

	    if (false) // debug
	       std::cout << ".... about to compare " << residue_spec_t(test_res) << " to "
			 << residues_vec.size() << " refining residues " << std::endl;
	    for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	       if (test_res == residues_vec[ires].second) {
                  if (! residues_vec[ires].first) {
                     found = true;
                     break;
                  }
	       }
	    }

	    if (! found) {
	       // OK, so this neighbour was not in the passed set of
	       // moving residues (and not already in nbr)... it can
	       // be a flanking residue then...

	       // check that it is not a bonded flanking residue...
	       for (unsigned int iflank=0; iflank<bonded_flanking_pairs.size(); iflank++) { 
		  if (bonded_flanking_pairs[iflank].res_1 == test_res) {
		     found = true;
		     break;
		  } 
		  if (bonded_flanking_pairs[iflank].res_2 == test_res) {
		     found = true;
		     break;
		  }
	       }

	       if (! found) {
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

// delete this
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

      mmdb::Atom *at_1 = atom[i];

      for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) {

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
               // Hmm... Is that sensible?
	    
	       std::string res_name_1 = at_1->GetResName();
	       std::string res_name_2 = at_2->GetResName();

	       if (false)
		  std::cout << "DEBUG:: here with " << atom_spec_t(at_1) << " " << atom_spec_t(at_2)
			    << " res_names " << res_name_1 << " " << res_name_2 << " "
			    << at_1->GetAtomName() << " " << at_2->GetAtomName()
			    << std::endl;

               if (false) {
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
	       bool is_1_4_related = ai.is_1_4(i, filtered_non_bonded_atom_indices[i][j], fixed_atom_flags);

	       if (false)
		  std::cout << "here C with at_1 " << atom_spec_t(at_1) << " at_2 " << atom_spec_t(at_2)
			    << " is_1_4_related " << is_1_4_related << std::endl;

	       if (is_1_4_related) {
                  if (in_same_ring_flag)
                     dist_min = 2.64; // was 2.7 but c.f. guanine ring distances
                  else
                     dist_min = 3.8;
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
				  is_H_non_bonded_contact,
				  fixed_atom_flags, dist_min);

	       if (false)
		  std::cout << "Adding NBC " << i << " " << filtered_non_bonded_atom_indices[i][j]
			    << " " << type_1 << " " << type_2 << " " 
			    << is_H_non_bonded_contact << " "
			    << fixed_atom_flags[0] << " " << fixed_atom_flags[1] << " "
			    << dist_min <<  "\n";

	       r.n_atoms_from_all_restraints = n_atoms; // for debugging crash in non-bonded contact
	                                                // restraints
	       r.restraints_index = size(); // likewise
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
      std::pair<bool, dictionary_residue_restraints_t> rest = geom.get_monomer_restraints(res_name, imol);
      const auto &dict = rest.second;
      if (rest.first) {
	 std::vector<std::vector<std::string> > ri = dict.get_ligand_ring_list();
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
   std::vector<bool> fixed_atom_flags = {false, false};
   is_1_4 = ai.is_1_4(idx_1, idx_2, fixed_atom_flags);
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

      if (r[ii].restraint_type == BOND_RESTRAINT) {
	 bonds[r[ii].atom_index_1].insert(r[ii].atom_index_2);
	 bonds[r[ii].atom_index_2].insert(r[ii].atom_index_1);
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
coot::restraints_container_t::reduced_angle_info_container_t::is_1_4(int indx_1, int indx_2,
								     const std::vector<bool> &fixed_atom_flags) const {

   // this function can be const because we don't use [] operator on the angles map.

   // This doesn't find 1-4 related main-chain when one (or more) of the atoms is in the residue
   // are fixed atoms (because fixed atoms don't have angle restraints)

   // We could catch some of those by looking to see if one atom (C(n-1)) is fixed and the
   // others are not (C(n), CA(n), N(n)) and that there is a bond between N(n) and C(n+1).

   bool f = false;

   bool has_a_fixed_atom = false;

   if (fixed_atom_flags.size() != 2) {
      std::cout << "ERROR:: in reduced_angle_info_container_t is_1_4(): " << fixed_atom_flags.size()
		<< std::endl;
      return false;
   }

   if (fixed_atom_flags[0]) has_a_fixed_atom = true;
   if (fixed_atom_flags[1]) has_a_fixed_atom = true;

   if (! has_a_fixed_atom) {
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
   } else {

      // *Does* have a fixed atom index

      bool fixed_1 = false;
      bool fixed_2 = false;
      if (fixed_atom_flags[0]) fixed_1 = true;
      if (fixed_atom_flags[1]) fixed_2 = true;

      if (fixed_2 && ! fixed_1) {
	 // key: atom index, data: a vector of the other 2 atoms in the angle restraints
	 std::map<int, std::vector<std::pair<int, int> > >::const_iterator it = angles.find(indx_1);
	 if (it != angles.end()) {
	    const std::vector<std::pair<int, int> > &v = it->second;
	    for (unsigned int ii=0; ii<v.size(); ii++) {
	       int idx_other_end = v[ii].second;
	       // is there a bond between idx_other_end and idx_2?
	       std::map<int, std::set<int> >::const_iterator it_bonds = bonds.find(idx_other_end);
	       if (it_bonds != bonds.end()) {
		  const std::set<int> &s = it_bonds->second;
		  if (s.find(indx_2) != s.end()) {
		     f = true;
		     break;
		  }
	       }
	    }
	 }
      }

      // same as above, but reversed indices
      if (fixed_1 && ! fixed_2) {
	 std::map<int, std::vector<std::pair<int, int> > >::const_iterator it = angles.find(indx_2);
	 if (it != angles.end()) {
	    const std::vector<std::pair<int, int> > &v = it->second;
	    for (unsigned int ii=0; ii<v.size(); ii++) {
	       int idx_other_end = v[ii].second;
	       // is there a bond between idx_other_end and idx_1?
	       std::map<int, std::set<int> >::const_iterator it_bonds = bonds.find(idx_other_end);
	       if (it_bonds != bonds.end()) {
		  const std::set<int> &s = it_bonds->second;
		  if (s.find(indx_1) != s.end()) {
		     f = true;
		     break;
		  }
	       }
	    }
	 }
      }
   }
   return f;
}

bool
coot::restraints_container_t::check_for_1_4_relation(int idx_1, int idx_2) const {

   bool is_1_4 = false;

   for (unsigned int ii=0; ii<restraints_vec.size(); ii++) {

      { // double indexing angle_1
	 const simple_restraint &restraint_i = restraints_vec[ii];
	 if (restraint_i.restraint_type == coot::ANGLE_RESTRAINT) {

	    if (idx_1 == restraint_i.atom_index_1 ||
		idx_1 == restraint_i.atom_index_3 ||
		idx_2 == restraint_i.atom_index_1 ||
		idx_2 == restraint_i.atom_index_3) { 

	       for (unsigned int jj=ii; jj<restraints_vec.size(); jj++) {
		  const simple_restraint &restraint_j = restraints_vec[jj];

		  if (jj != ii) { // check both indices
		     if (restraint_j.restraint_type == coot::ANGLE_RESTRAINT) {

			if (idx_2 == restraint_j.atom_index_1 ||
			    idx_2 == restraint_j.atom_index_3 ||
			    idx_1 == restraint_j.atom_index_1 ||
			    idx_1 == restraint_j.atom_index_3) {

			   if (false)
			      std::cout << "check_for_1_4_relation() indices "
					<< idx_1 << " " << idx_2
					<< " examining angle restraint pair "
					<< restraint_i.atom_index_1 << " "
					<< restraint_i.atom_index_2 << " "
					<< restraint_i.atom_index_3 << " and "
					<< restraint_j.atom_index_1 << " "
					<< restraint_j.atom_index_2 << " "
					<< restraint_j.atom_index_3 << std::endl;

			   if ((restraint_i.atom_index_2 == restraint_j.atom_index_1) ||
			       (restraint_i.atom_index_2 == restraint_j.atom_index_3)) {

			      if ((restraint_j.atom_index_2 == restraint_i.atom_index_1) ||
				  (restraint_j.atom_index_2 == restraint_i.atom_index_3)) {

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

std::pair<double, double>
coot::simple_restraint::distortion(mmdb::PAtom *atoms, const double &lj_epsilon) const {

   // distortion and bond length delta - bone-headed recalculation - because interface
   // is from atoms, not gsl vector x

   std::pair<double, double> distortion_pair(-1.0, -1.0);

   if (restraint_type == CHIRAL_VOLUME_RESTRAINT) {
      mmdb::Atom *at_c = atoms[atom_index_centre];
      mmdb::Atom *at_1 = atoms[atom_index_1];
      mmdb::Atom *at_2 = atoms[atom_index_2];
      mmdb::Atom *at_3 = atoms[atom_index_3];

      clipper::Coord_orth centre = co(at_c);
      clipper::Coord_orth a1 = co(at_1);
      clipper::Coord_orth a2 = co(at_2);
      clipper::Coord_orth a3 = co(at_3);

      clipper::Coord_orth a = a1 - centre;
      clipper::Coord_orth b = a2 - centre;
      clipper::Coord_orth c = a3 - centre;

      double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
      double delta = cv - target_chiral_volume;
      distortion_pair.first = delta*delta/(sigma * sigma);
      distortion_pair.second = delta;
   }

   if (restraint_type == BOND_RESTRAINT) {
      mmdb::Atom *at_1 = atoms[atom_index_1];
      mmdb::Atom *at_2 = atoms[atom_index_2];
      if (at_1 && at_2) {
         clipper::Coord_orth p1 = co(at_1);
         clipper::Coord_orth p2 = co(at_2);
         double d = sqrt((p2-p1).lengthsq());
	 if (false)
	    std::cout << atom_spec_t(at_1) << " " << atom_spec_t(at_2)
		      << " d " << d << " target_value " << target_value << std::endl;
         double delta = d - target_value;
         double z = delta/sigma;
         double pen_score = z * z;
         distortion_pair.first = pen_score;
         distortion_pair.second = delta;
      }
   }

   if (restraint_type == GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
      mmdb::Atom *at_1 = atoms[atom_index_1];
      mmdb::Atom *at_2 = atoms[atom_index_2];
      if (at_1 && at_2) {
         clipper::Coord_orth p1 = co(at_1);
         clipper::Coord_orth p2 = co(at_2);
         double d = sqrt((p2-p1).lengthsq());
	 if (false)
	    std::cout << atom_spec_t(at_1) << " " << atom_spec_t(at_2)
		      << " d " << d << " target_value " << target_value << std::endl;
         double alpha = 0.01; // pass this
         double delta = d - target_value;
         double z = delta/sigma;
         double distortion = z*z/(1+alpha*z*z);
         distortion_pair.first = distortion;
         distortion_pair.second = delta;
      }
   }

   if (restraint_type == ANGLE_RESTRAINT) {
      mmdb::Atom *at_1 = atoms[atom_index_1];
      mmdb::Atom *at_2 = atoms[atom_index_2];
      mmdb::Atom *at_3 = atoms[atom_index_3];
      if (at_1 && at_2 && at_3) {
         double angle_deg = coot::angle(at_1, at_2, at_3);
         double delta = angle_deg - target_value;
         double z = delta/sigma;
         double distortion = z*z;
         distortion_pair.first = distortion;
         distortion_pair.second = delta;
      }
   }

   if (restraint_type == NON_BONDED_CONTACT_RESTRAINT) {
      mmdb::Atom *at_1 = atoms[atom_index_1];
      mmdb::Atom *at_2 = atoms[atom_index_2];
      if (at_1 && at_2) {
         clipper::Coord_orth p1 = co(at_1);
         clipper::Coord_orth p2 = co(at_2);

         double dist_sq = clipper::Coord_orth(p2-p1).lengthsq();
         double dist = sqrt(dist_sq);
         double dist_delta = dist - target_value;

         // float bl = sqrt(bl_sq);
         // float delta =  bl - target_value;

         float V_lj = 0;

         // the value lj_sigma is r when is the potential is 0.
         // for lj_sigma + delta the potential is negative
         // for lj_sigma - delta the potential is positive
         // the potential is at a minimum at lj_r_min
         // 
         // so if target_value is say 3.4A, lj_sigma is 3.4 sigma
         // and lj_r_min ~ 3.4 * 1.122 = 3.82
         double lj_sigma = target_value;
         // double lj_r_min = std::pow(2.0, 1.0/6.0) * lj_sigma;
         double lj_r_min = 1.122462048309373 * lj_sigma;

         double max_dist = 2.5 * lj_sigma; // r_max

         if (dist_sq < 0.81) dist_sq = 0.81; // 0.9^2
         double alpha_sqrd = lj_r_min*lj_r_min/dist_sq;
         double alpha_up_6  = alpha_sqrd * alpha_sqrd * alpha_sqrd;
         double alpha_up_12 = alpha_up_6 * alpha_up_6;
         V_lj = lj_epsilon * (alpha_up_12 - 2.0 * alpha_up_6);

         // offset the Vlj so that it is zero at r_max (beyond which we no longer
         // consider contributions to the distortion)

         double Vlj_at_rmax = -0.016316891136 * lj_epsilon; // see Lennard-Jones truncated and shifted for

         V_lj += Vlj_at_rmax;

         distortion_pair.first = V_lj;
         distortion_pair.second = dist_delta;
      }
   }
   return distortion_pair;
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

   // 20191005-PE we add the bonded_atom_indices for angles between H so that NBCs
   // are not generted between these Hs.

   if (N_index >= 0) {
      std::map<std::string, int>::const_iterator it_1, it_2, it_3;
      for(it_1=h1s.begin(); it_1!=h1s.end(); it_1++) {
	 const std::string &key_alt_conf = it_1->first;
	 it_2 = h2s.find(key_alt_conf);
	 it_3 = h3s.find(key_alt_conf);
	 if (it_2 != h2s.end()) {
	    std::vector<bool> fixed_flags_a12 = make_fixed_flags(it_1->second, N_index, it_2->second);
	    add(ANGLE_RESTRAINT, it_1->second, N_index, it_2->second, fixed_flags_a12, 109.5, 2.0, false);
            bonded_atom_indices[it_1->second].insert(it_2->second);
            bonded_atom_indices[it_2->second].insert(it_1->second);
	 }

	 if (it_3 != h3s.end()) {
	    std::vector<bool> fixed_flags_a13 = make_fixed_flags(it_1->second, N_index, it_3->second);
	    add(ANGLE_RESTRAINT, it_1->second, N_index, it_3->second, fixed_flags_a13, 109.5, 2.0, false);
            bonded_atom_indices[it_1->second].insert(it_3->second);
            bonded_atom_indices[it_3->second].insert(it_1->second);
	 }

	 if (it_2 != h2s.end() && it_3 != h3s.end()) {
	    std::vector<bool> fixed_flags_a23 = make_fixed_flags(it_2->second, N_index, it_3->second);
	    add(ANGLE_RESTRAINT, it_2->second, N_index, it_3->second, fixed_flags_a23, 109.5, 2.0, false);
            bonded_atom_indices[it_2->second].insert(it_3->second);
            bonded_atom_indices[it_3->second].insert(it_2->second);
	 }
      }
   }

   return rc;
}


std::vector<unsigned int>
coot::restraints_container_t::make_torsion_restraint_indices_vector() const {

   std::vector<unsigned int> v;
   v.reserve(20);
   unsigned int n = restraints_vec.size(); 
   for (unsigned int ir=0; ir<n; ir++) {
      if (restraints_vec[ir].restraint_type == TORSION_RESTRAINT) {
         v.push_back(ir);
      }
   }
   return v;
}


bool
coot::restraints_container_t::add_or_replace_torsion_restraints_with_closest_rotamer_restraints(const std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > &rotamer_torsions) {

   bool status = false;
   bool debug = false;

   // we have a book-keeping problem: the torsions in the monomer library (ie. the definintion of chi1 (for example)
   // is not the same in the monomer libary and the atoms from rotamer_atoms() (and presumably molprobity dictionary).
   // So, if there is an active monomer library torsion - set from the dictionary, then we need to delete
   // the current torsion restraints for that residue.

   std::vector<unsigned int> tri = make_torsion_restraint_indices_vector();

   for (unsigned int ir=0; ir<rotamer_torsions.size(); ir++) {
      mmdb::Residue *rotamer_residue = rotamer_torsions[ir].first;
      for (unsigned int i=0; i<residues_vec.size(); i++) {
         if (! residues_vec[i].first) {
            mmdb::Residue *residue_p = residues_vec[i].second;
            if (residue_p == rotamer_residue) {
               for (unsigned int j=0; j<rotamer_torsions[ir].second.size(); j++) {
                  const dict_torsion_restraint_t &new_torsion_restraint = rotamer_torsions[ir].second[j];
                  if (debug)
                     std::cout << "debug:: in add_or_replace_torsion_restraints_with_closest_rotamer_restraints() consider "
                               << new_torsion_restraint << std::endl;
                  mmdb::Atom **residue_atoms = 0;
                  int n_residue_atoms;
                  residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                  bool replaced = replace_torsion_restraint(new_torsion_restraint, residue_atoms, n_residue_atoms, tri);
                  if (debug)
                     std::cout << "debug:: in add_or_replace_torsion_restraints_with_closest_rotamer_restraints() replaced flag " << replaced << " for "
                               << new_torsion_restraint << std::endl;
                  if (! replaced)
                     status = add_torsion_internal(new_torsion_restraint, residue_atoms, n_residue_atoms, torsion_restraints_weight);
               }
            }
         }
      }
   }

   if (debug)
      std::cout << "add_or_replace_torsion_restraints_with_closest_rotamer_restraints() returning " << status << std::endl;

   return status;
}

int
coot::restraints_container_t::get_atom_index_for_restraint_using_alt_conf(const std::string &atom_name,
                                                                          const std::string &alt_conf,
                                                                          mmdb::PPAtom res_selection, int num_res_atoms) const {

   int idx = -1;
   for (int i=0; i<num_res_atoms; i++) {
      mmdb::Atom *at = res_selection[i];
      std::string n(at->GetAtomName());
      if (n == atom_name) {
         std::string a(at->altLoc);
         if (a.empty() || a == alt_conf) {
            at->GetUDData(udd_atom_index_handle, idx);
         }
      }
   }
   return idx;
}

bool
coot::restraints_container_t::replace_torsion_restraint(const coot::dict_torsion_restraint_t &new_torsion_restraint,
                                                        mmdb::PPAtom res_selection, int num_res_atoms,
                                                        const std::vector<unsigned int> &torsion_restraint_indices) {

   // I don't want to add a rotamer torsion restraint if it's already there

   bool replaced = false;

   unsigned int n = torsion_restraint_indices.size(); 

   std::string alt_conf; // run through all alt confs in the residue ideally.
   if (true) {
      int idx_1 = get_atom_index_for_restraint_using_alt_conf(new_torsion_restraint.atom_id_1_4c(), alt_conf, res_selection, num_res_atoms);
      if (idx_1 >= 0) {
         int idx_2 = get_atom_index_for_restraint_using_alt_conf(new_torsion_restraint.atom_id_2_4c(), alt_conf, res_selection, num_res_atoms);
         if (idx_2 >= 0) {
            int idx_3 = get_atom_index_for_restraint_using_alt_conf(new_torsion_restraint.atom_id_3_4c(), alt_conf, res_selection, num_res_atoms);
            if (idx_3 >= 0) {
               int idx_4 = get_atom_index_for_restraint_using_alt_conf(new_torsion_restraint.atom_id_4_4c(), alt_conf, res_selection, num_res_atoms);
               if (idx_4 >= 0) {
                  // OK, so we have real atoms for a restraints, does a torsion for that set of atoms exist already?
                  for (unsigned int it=0; it<n; it++) {
                     simple_restraint &rest = restraints_vec[torsion_restraint_indices[it]];
                     if (rest.restraint_type == TORSION_RESTRAINT) {
                        // std::cout << " comparing atom indices " << idx_1 << " " << idx_2 << " " << idx_3 << " " << idx_4 << "  vs " << rest.atom_index_1 << " " << rest.atom_index_2 << " " << rest.atom_index_3 << " " << rest.atom_index_4 << std::endl;
                        if (idx_1 == rest.atom_index_1) {
                           if (idx_2 == rest.atom_index_2) {
                              if (idx_3 == rest.atom_index_3) {

                                 if (idx_4 != rest.atom_index_4)
                                    rest.atom_index_4 = idx_4;
                                 rest.target_value = new_torsion_restraint.angle(); // angle mean "torsion" - hmm!
                                 replaced = true;
                                 if (true)
                                    std::cout << "debug:: in replace_torsion_restraint() replacing restraints with " << new_torsion_restraint << std::endl;
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
   }
   return replaced;
}


void
coot::restraints_container_t::set_torsion_restraints_weight(double w) {
   // we make the changes as the refinement is running possibly. Might be dangerous.
   torsion_restraints_weight = w;
   for (auto &r : restraints_vec)
      if (r.restraint_type == TORSION_RESTRAINT)
         r.torsion_restraint_weight = w;
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
   unsigned int idx_rest = size() -1;
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

   if (! found) {
      restraints_indices.back().push_back(idx_rest);
   }

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

void
coot::restraints_container_t::set_use_harmonic_approximations_for_nbcs(bool flag) {

   bool needs_reset_flag = false;
   simple_restraint::nbc_function_t nbc_func_type = simple_restraint::LENNARD_JONES;
   if (flag)
      nbc_func_type = simple_restraint::HARMONIC;
   for (unsigned int i=0; i< restraints_vec.size(); i++) {
      simple_restraint &restraint = restraints_vec[i];
      if (restraint.nbc_function != nbc_func_type) {
         restraint.nbc_function = nbc_func_type;
         needs_reset_flag = true;
      }
   }
   if (needs_reset_flag)
      set_needs_reset();
}

int
coot::restraints_container_t::write_new_atoms(std::string pdb_file_name) { 

   //
   int status = -1;
   if (mol != NULL) {
      // return 0 on success, non-zero on failure.
      status = mol->WritePDBASCII(pdb_file_name.c_str());
      if (status == 0)
	 // std::cout << "INFO:: output file: " << pdb_file_name
	 //          << " written." << std::endl;
	 logger.log(log_t::INFO, "output file:", pdb_file_name, "written.");
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

   // std::cout << "INFO:: There are " << n_atoms << " atoms" << std::endl;
   logger.log(log_t::INFO, "There are", n_atoms, "atoms");
   // std::cout << "INFO:: There are " << size() << " restraints" << std::endl;
   logger.log(log_t::INFO, "There are", size(), "restraints");

   for (unsigned int i=0; i< restraints_vec.size(); i++) {
      const simple_restraint &restraint = restraints_vec[i];
      if (restraint.restraint_type == coot::TORSION_RESTRAINT) {
	 // std::cout << "INFO:: restraint " << i << " is of type "
	 //          << restraint.restraint_type << std::endl;
	 logger.log(log_t::INFO, "restraint", i, "is of type", restraint.restraint_type);

	 std::cout << restraint.atom_index_1 << " "
		   << restraint.atom_index_2 << " "
		   << restraint.atom_index_3 << " "
		   << restraint.atom_index_4 << " "
		   << restraint.target_value << " "
		   << restraint.sigma << " " << std::endl
		   << " with "
  		   << restraint.plane_atom_index.size() << " vector atoms " << std::endl
		   << " with periodicity "
  		   << restraint.periodicity << std::endl;
      }
      // std::cout << "INFO:: restraint number " << i << " is restraint_type " <<
      //    restraint.restraint_type << std::endl;
      logger.log(log_t::INFO, "restraint number", i, "is restraint_type", restraint.restraint_type);
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
   
	 clipper::Xmap<float> dummy_xmap;

         std::vector<std::pair<bool,mmdb::Residue *> > residues;
         residues.push_back(std::pair<bool,mmdb::Residue *>(false, residue_p));

         coot::restraints_container_t restraints(residues, geom, mol, &dummy_xmap);
   
	 // restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
	 restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
	 pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	 bool do_internal_torsions = true;
	 bool do_trans_peptide_restraints = true;
         int n_threads = coot::get_max_number_of_threads();
         ctpl::thread_pool thread_pool(n_threads);
         restraints.thread_pool(&thread_pool, n_threads);
	 restraints.make_restraints(imol, geom, flags, do_internal_torsions,
				    do_trans_peptide_restraints, 0, 0, true, true, false, pseudos);

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


// maybe this should have its own file

// Because this is used for "celebration", I think that the Rama score
// should be analysed (used as a filter) even if it is not used for
// refinement.  Otherwise we could get celebration for a model with
// Rama outliers - and that would be misleading for most people (not
// me).
//
// Along similar lines, the fit to density of the residues in the atom selection should
// be better than ... something.. average?
//
std::pair<bool, std::string>
coot::refinement_results_t::hooray() const {

   bool status = true;
   std::string message;
   for (unsigned int i=0; i<lights.size(); i++) {
      const refinement_lights_info_t &light = lights[i];
      if (false)
         // std::cout << "INFO:: for lights index " << i << " " << light.name << " " << light.value << std::endl;
         logger.log(log_t::INFO, "for lights index", i, light.name, light.value);
      float crit_value = 1.0; // 20210906-PE was 1.4
      if (light.name == "Trans_peptide")
         crit_value = 2.0; // 20210906-PE  was 6.0. Should Trans-peptide even be tested in hooray()?
      if (light.value > crit_value) {
         // std::cout << "Boo for lights index " << i << " " << light.name << " " << light.value << std::endl;
         status = false;
      }
   }

   // to prevent hooray() returning true "too often" only make it return true if there
   // have been pull restraints
   //

   int n_pull_restraints = sorted_atom_pulls.size();
   if (n_pull_restraints == 0)
      status = false;

   return std::pair<bool, std::string> (status, message);
}
