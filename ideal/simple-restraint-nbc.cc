/* ideal/simple-restraint-nbc.cc
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


// Fixed-atom flags and non-bonded contact restraints for restraints_container_t

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

