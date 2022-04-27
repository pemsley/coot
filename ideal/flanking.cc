/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2016 by Medical Research Council
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


// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL

#include "simple-restraint.hh"

#include "coot-utils/coot-coord-extras.hh"  // is_nucleotide_by_dict

// also uses class variable non_bonded_neighbour_residues
void
coot::restraints_container_t::add_fixed_atoms_from_flanking_residues(const coot::bonded_pair_container_t &bpc) {

   // std::cout << "debug:: add_fixed_atoms_from_flanking_residues() called " << bpc.size() << std::endl;

   std::vector<mmdb::Residue *> residues_for_fixed_atoms;

   if (false) { // debug
      std::cout << "------ in add_fixed_atoms_from_flanking_residues() here are the elements of the bpc" << std::endl;
      for (unsigned int i=0; i<bpc.size(); i++) {
         std::cout << " bpc " << i << " " << bpc[i] << std::endl;
      }
   }

   for (unsigned int i=0; i<bpc.size(); i++) {
      if (bpc[i].is_fixed_first)
	 residues_for_fixed_atoms.push_back(bpc[i].res_1);
      if (bpc[i].is_fixed_second)
	 residues_for_fixed_atoms.push_back(bpc[i].res_2);
   }

   if (false) { // debug
      std::cout << "----in add_fixed_atoms_from_flanking_residues() here are the residues_for_fixed_atoms"
                << std::endl;
      for (unsigned int i=0; i<residues_for_fixed_atoms.size(); i++) {
         std::cout << "             " << residue_spec_t(residues_for_fixed_atoms[i]) << std::endl;
      }
   }

   for (unsigned int i=0; i<residues_for_fixed_atoms.size(); i++) {
      int idx;
      mmdb::PPAtom residue_atoms = 0;
      int n_res_atoms;
      residues_for_fixed_atoms[i]->GetAtomTable(residue_atoms, n_res_atoms);
      for (int iat=0; iat<n_res_atoms; iat++) { 
	 mmdb::Atom *at = residue_atoms[iat];
	 if (! (at->GetUDData(udd_atom_index_handle, idx) == mmdb::UDDATA_Ok)) {
	    std::cout << "ERROR:: bad UDD for atom " << iat << std::endl;
	 } else {
	    fixed_atom_indices.insert(idx); // hello grep: in add_fixed_atoms_from_flanking_residues()
	 }
      }
   }
}


// Not a "flanking function" - maybe it shouldn't be in this file
//
// use non_bonded_neighbour_residues to extend the set of fixed atoms.
// (these are not flanking residues)
void
coot::restraints_container_t::add_fixed_atoms_from_non_bonded_neighbours() {

   // std::cout << "####################### in add_fixed_atoms_from_non_bonded_neighbours() "
   // << non_bonded_neighbour_residues.size() << std::endl;

   for (std::size_t jj=0; jj<non_bonded_neighbour_residues.size(); jj++) {
      // std::cout << "debug:: in add_fixed_atoms_from_non_bonded_neighbours "
      // << residue_spec_t(non_bonded_neighbour_residues[jj]) << std::endl;
      mmdb::Residue *residue_p = non_bonded_neighbour_residues[jj];
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 int idx = -1;
	 mmdb::Atom *at = residue_atoms[iat];
	 if (! (at->GetUDData(udd_atom_index_handle, idx) == mmdb::UDDATA_Ok)) {
	    std::cout << "ERROR:: in add_fixed_atoms_from_non_bonded_neighbours() "
		      << " bad UDD for atom " << iat << std::endl;
	 } else {
	    if (std::find(fixed_atom_indices.begin(),
			  fixed_atom_indices.end(), idx) == fixed_atom_indices.end())
	       fixed_atom_indices.insert(idx); // hello grep: in add_fixed_atoms_from_non_bonded_neighbours()
	 }
      }
   }
}


void
coot::restraints_container_t::add_fixed_atoms_from_flanking_residues(bool have_flanking_residue_at_start,
								     bool have_flanking_residue_at_end,
								     int iselection_start_res,
								     int iselection_end_res) {

   if (false)
      std::cout << "debug:: in add_fixed_atoms_from_flanking_residues() "
		<< have_flanking_residue_at_start << " "
		<< have_flanking_residue_at_end << " "
		<< std::endl;

   if (have_flanking_residue_at_start || have_flanking_residue_at_end) {
      for (int iat=0; iat<n_atoms; iat++) {
	 mmdb::Atom *at = atom[iat];
	 if (have_flanking_residue_at_start) {
	    if (at->residue->GetSeqNum() == iselection_start_res) {
	       // perhaps this should be a set - yes.
	       fixed_atom_indices.insert(iat);
	    }
	 }
	 if (have_flanking_residue_at_end) {
	    if (at->residue->GetSeqNum() == iselection_end_res) {
	       fixed_atom_indices.insert(iat);
	    }
	 }
      }
   }
}


coot::bonded_pair_container_t
coot::restraints_container_t::make_flanking_atoms_restraints(const coot::protein_geometry &geom,
							     bool do_rama_plot_restraints,
							     bool do_trans_peptide_restraints) {

   coot::bonded_pair_container_t bonded_residue_pairs = bonded_flanking_residues(geom);
   int iv = make_link_restraints_by_pairs(geom, bonded_residue_pairs, do_trans_peptide_restraints, "Flanking residue");

   int n_rama_restraints = -1; // unset, don't output an info line if
			       // do_rama_plot_restraints is not set.
   if (do_rama_plot_restraints) {
      // e.g 1 free 2 free 3 flanking (fixed).
      n_rama_restraints = make_flanking_atoms_rama_restraints(geom);  // returns 0 or something.
   }
   return bonded_residue_pairs;
}


int coot::restraints_container_t::make_flanking_atoms_rama_restraints(const protein_geometry &geom) {
   int n_rama_restraints = 0;

   if (istart_minus_flag && iend_plus_flag) {  // have flanking residues

      std::vector<coot::ramachandran_restraint_flanking_residues_helper_t> vrrfr;
      coot::ramachandran_restraint_flanking_residues_helper_t rrfr_1;
      rrfr_1.resno_first = istart_res-1;
      rrfr_1.resno_third = istart_res+1;
      rrfr_1.is_fixed[0] = 1;
      if (istart_res == iend_res) // i.e. just one moving residue
	 rrfr_1.is_fixed[2] = 1;
      vrrfr.push_back(rrfr_1);

      // we don't want to add 2 sets of flanking ramas for when
      // refining just one residue (with 2 flanking residues)
      if (istart_res != iend_res) { 
	 coot::ramachandran_restraint_flanking_residues_helper_t rrfr_2;
	 rrfr_2.resno_first = iend_res-1;
	 rrfr_2.resno_third = iend_res+1;
	 rrfr_2.is_fixed[2] = 1;
	 vrrfr.push_back(rrfr_2);
      }

      for (unsigned int iround=0; iround<vrrfr.size(); iround++) { 
      
	 int selHnd = mol->NewSelection();
	 mmdb::PPResidue SelResidue = NULL;
	 int nSelResidues;
	 mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		      chain_id_save.c_str(), // Chain(s)
		      vrrfr[iround].resno_first,   "*",  // starting res
		      vrrfr[iround].resno_third,   "*",  // ending res
		      "*",  // residue name
		      "*",  // Residue must contain this atom name?
		      "*",  // Residue must contain this Element?
		      "*",  // altLocs
		      mmdb::SKEY_NEW); // selection key 
	 mol->GetSelIndex ( selHnd, SelResidue,nSelResidues );
	 // std::cout << "DEBUG:: GetSelIndex (make_flanking_atoms_rama_restraints) returned " 
	 // << nSelResidues << " residues (for flanking rama restraints)" << std::endl;
      
	 if (nSelResidues == 3) {
	    // super careful would mean that we check the link type of
	    // both pairs before calling this function:

	    if (0) { // debugging fixed atoms
	       for (int i=0; i<3; i++)
		  std::cout << "   make_flanking_atoms_rama_restraints() calling add_rama() with index "
			    << i << " resno " 
			    << coot::residue_spec_t(SelResidue[i]) << " Fixed: "
			    << vrrfr[iround].is_fixed[i] << std::endl;
	    }

	    add_rama("TRANS",
		     SelResidue[0], SelResidue[1], SelResidue[2],
		     vrrfr[iround].is_fixed[0],
		     vrrfr[iround].is_fixed[1],
		     vrrfr[iround].is_fixed[2], geom);
	 }
      
	 mol->DeleteSelection(selHnd);
      }
   }
      
   return n_rama_restraints;
}

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_flanking_residues(const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t bpc;

   // residue n is at the end of the active selection.  What is
   // residue n+1 from the mol? We will make a bonded pair of residues
   // n and n+1, and make the is_fixed_residue true for the n+1
   // (flanking) residue (and False for residue n, obviously).
   //
   // We need to ignore (don't add) a [n,n+1] pair if residue n+1 is
   // in the vector of residues that we pass.
   // 
   // We need to do this for [n,n+1] pairs and [n,n-1] pairs.
   //
   // First, what residue n?  I suppose that that would be a member of
   // the residue vector passed to the constructor.  But what if we
   // are not using the residues vector constructor?
   // 

   if (from_residue_vector)
      bpc = bonded_flanking_residues_by_residue_vector(geom);
   else
      bpc = bonded_flanking_residues_by_linear(geom);

   return bpc;
}



coot::bonded_pair_container_t
coot::restraints_container_t::bonded_flanking_residues_by_linear(const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t bpc;
   std::string link_type = "TRANS";
   mmdb::PPResidue SelResidue = NULL;
   int nSelResidues;
   int selHnd = mol->NewSelection();
   mol->Select (selHnd,mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res-1,   "*",  // starting res
		istart_res,     "*",  // ending res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",  // altLocs
		mmdb::SKEY_NEW); // selection key 
   mol->GetSelIndex (selHnd, SelResidue, nSelResidues);
   std::cout << "INFO:: GetSelIndex (make_flanking_atoms_restraints) returned " 
	     << nSelResidues << " residues (flanking restraints)" << std::endl;
   if (nSelResidues > 1) {
      link_type = find_link_type(SelResidue[0], SelResidue[1], geom);
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage
   
      coot::bonded_pair_t bp(SelResidue[0], SelResidue[1], 1, 0, link_type);
      bpc.try_add(bp);
   }
   mol->DeleteSelection(selHnd);
   
   // And now again for the C-terminal flanking residue:
   // 
   selHnd = mol->NewSelection();
   mol->Select (selHnd,mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		iend_res,   "*",  // starting res
		iend_res+1, "*",  // ending res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",  // altLocs
		mmdb::SKEY_NEW); // selection key 
   mol->GetSelIndex (selHnd, SelResidue, nSelResidues);
   std::cout << "INFO:: GetSelIndex (make_flanking_atoms_restraints) returned " 
	     << nSelResidues << " residues (flanking restraints)" << std::endl;
   if (nSelResidues > 1) {
      link_type = find_link_type(SelResidue[0], SelResidue[1], geom);
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage
      coot::bonded_pair_t bp(SelResidue[0], SelResidue[1], 0, 1, link_type);
      bpc.try_add(bp);
   }
   mol->DeleteSelection(selHnd);

   // std::cout << "DEBUG:: bonded_flanking_residues_by_linear() reutrns " << bpc;
   return bpc;
}



coot::bonded_pair_container_t
coot::restraints_container_t::bonded_flanking_residues_by_residue_vector(const std::map<mmdb::Residue *, std::set<mmdb::Residue *> > &neighbour_set,
									 const coot::protein_geometry &geom) const {

   // Don't make flanking residue restraints if both residues are fixed!

   coot::bonded_pair_container_t bpc;

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;

   // 20180104 2.0 is a terrible distance.  Sometimes it will find a disulfide bond
   //              and the next it will not.  Use 2.3
   float dist_crit = 2.3; // 20170924-PE was 3.0 but this made a horrible link in a tight turn
                          // (which I suspect is not uncommon) crazy-neighbour-refine-519.pdb
                          // for EMDB 6224.
                          // 520 was bonded to 522 in a neighb (3-residue) refine on 519.
                          // This function is called by init (and (I think) make_restraints)
                          // init doesn't set bonded_pairs_container (make_restraints does that).

   // Don't forget to consider the case were we refine one residue,
   // and that is the ASN for a glycosylation.  So the neighbouring
   // residues and the GLC/NAG (say) will be flanking residues.
   //
   // But that is a really hard thing - isn't it? To find a GLC that
   // is connected to this residue?

   if (false)
      std::cout << "DEBUG:: here in bonded_flanking_residues_by_residue_vector() bonded_pairs_container has size "
		<< bonded_pairs_container.size() << std::endl;

   if (false) { // debug
      for (it=neighbour_set.begin(); it != neighbour_set.end(); ++it) {
	 std::cout << "Residue " << residue_spec_t(it->first) << " has neighbours ";
	 const std::set<mmdb::Residue *> &neighbours = it->second;
	 std::set<mmdb::Residue *>::const_iterator it_set;
	 for (it_set=neighbours.begin(); it_set!=neighbours.end(); ++it_set) {
	    std::cout << " " << residue_spec_t(*it_set);
	 }
	 std::cout << std::endl;
      }
   }

   // 20180220 no longer iterate on the residue_vec. Elinor Breiner
   //
   // Use this instead:
   for (it=neighbour_set.begin(); it != neighbour_set.end(); ++it) {

      const std::set<mmdb::Residue *> &neighbours = it->second;
      std::set<mmdb::Residue *>::const_iterator it_set;
      for (it_set=neighbours.begin(); it_set!=neighbours.end(); ++it_set) {

	 if (false) {
	    std::cout << "base residue " << it->first << " " << residue_spec_t(it->first) << std::endl;
	    std::cout << "checking for set member " << *it_set << " " << residue_spec_t(*it_set) << " in "
		      << residues_vec.size() << " residue_vec residues " << std::endl;
	 }

	 bool found = false;
	 for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	    // pointer comparison
	    if (*it_set == residues_vec[ires].second) {
               found = true;
               break;
	    }
	 }

	 if (! found) {
	    // OK, so this neighbour was not in the passed set of
	    // moving residues, it can be a flanking residue then...
	    std::pair<bool, float> d = closest_approach(*it_set, it->first);

	    if (d.first) {
	       if (d.second < dist_crit) {

                  unsigned int n_fixed_residues = 0;
                  for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
                     if (*it_set == residues_vec[ires].second)
                        if (residues_vec[ires].first) {
                           n_fixed_residues++;
                        }
                     if (it->first == residues_vec[ires].second)
                        if (residues_vec[ires].first) {
                           n_fixed_residues++;
                        }
                  }

                  if (false)
                     std::cout << "DEBUG:: here in bonded_flanking_residues_by_residue_vector() "
                               << " considering residues " << residue_spec_t(*it_set) << " "
                               << residue_spec_t(it->first) << " with n_fixed_residues "
                               << n_fixed_residues << std::endl;

                  if (n_fixed_residues != 2) {
                     std::pair<std::string, bool> l = find_link_type_complicado(*it_set, it->first, geom);
                     const std::string &link_type = l.first;
                     if (! link_type.empty()) {
                        const bool &order_switch_flag = l.second;
                        if (! order_switch_flag) {
                           coot::bonded_pair_t bp(*it_set, it->first, true, false, link_type);
                           bpc.try_add(bp);
                        } else {
                           coot::bonded_pair_t bp(it->first, *it_set, false, true, link_type);
                           bpc.try_add(bp);
                        }
                     }
                  }
	       }
	    }
	 }
      }
   }

   // does your linking problem lie in here?
   //
   bpc.filter();

   return bpc;
}


// old (pre-Weizmann)
coot::bonded_pair_container_t
coot::restraints_container_t::bonded_flanking_residues_by_residue_vector(const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t bpc;

   // 20180104 2.0 is a terrible distance.  Sometimes it will find a disulfide bond
   //              and the next it will not.  Use 2.3
   float dist_crit = 2.3; // 20170924-PE was 3.0 but this made a horrible link in a tight turn
                          // (which I suspect is not uncommon) crazy-neighbour-refine-519.pdb
                          // for EMDB 6224.
                          // 520 was bonded to 522 in a neighb (3-residue) refine on 519.
                          // This function is called by init (and (I think) make_restraints)
                          // init doesn't set bonded_pairs_container (make_restraints does that).

   // Don't forget to consider the case were we refine one residue,
   // and that is the ASN for a glycosylation.  So the neighbouring
   // residues and the GLC/NAG (say) will be flanking residues.
   //
   // But that is a really hard thing - isn't it? To find a GLC that
   // is connected to this residue?

   if (false)
      std::cout << "DEBUG:: here in bonded_flanking_residues_by_residue_vector() bonded_pairs_container has size "
		<< bonded_pairs_container.size() << std::endl;

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > neighbour_set = residues_near_residues(residues_vec, mol, dist_crit);
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;

   if (false) { // debug
      for (it=neighbour_set.begin(); it != neighbour_set.end(); it++) {
	 std::cout << "Residue " << residue_spec_t(it->first) << " has neighbours ";
	 const std::set<mmdb::Residue *> &neighbours = it->second;
	 std::set<mmdb::Residue *>::const_iterator it_set;
	 for (it_set=neighbours.begin(); it_set!=neighbours.end(); it_set++) {
	    std::cout << " " << residue_spec_t(*it_set);
	 }
	 std::cout << std::endl;
      }
   }

   // 20180220 no longer iterate on the residue_vec. Elinor Breiner
   //
   // Use this instead:
   for (it=neighbour_set.begin(); it != neighbour_set.end(); it++) {

      const std::set<mmdb::Residue *> &neighbours = it->second;
      std::set<mmdb::Residue *>::const_iterator it_set;
      for (it_set=neighbours.begin(); it_set!=neighbours.end(); it_set++) {

	 // std::cout << "base residue " << it->first << " " << residue_spec_t(it->first) << std::endl;
	 // std::cout << "checking for set member " << *it_set << " " << residue_spec_t(*it_set) << " in "
	 // << residues_vec.size() << " residue_vec residues " << std::endl;

	 bool found = false;
	 for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	    // pointer comparison
	    if (*it_set == residues_vec[ires].second) {
	       found = true;
	       break;
	    }
	 }

	 if (! found) {
	    // OK, so this neighbour was not in the passed set of
	    // moving residues, it can be a flanking residue then...
	    std::pair<bool, float> d = closest_approach(*it_set, it->first);

	    if (d.first) {
	       if (d.second < dist_crit) {

		  std::pair<std::string, bool> l = find_link_type_complicado(*it_set, it->first, geom);
		  const std::string &link_type = l.first;
		  if (! link_type.empty()) {
		     const bool &order_switch_flag = l.second;
		     if (! order_switch_flag) {
			coot::bonded_pair_t bp(*it_set, it->first, 1, 0, link_type);
			bpc.try_add(bp);
		     } else {
			coot::bonded_pair_t bp(it->first, *it_set, 0, 1, link_type);
			bpc.try_add(bp);
		     }
		  }
	       }
	    }
	 }
      } 
   }
   return bpc;
}


#endif // HAVE_GSL
