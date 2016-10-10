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

#include <string.h> // for strcmp


// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL


// #include <fstream>
// #include <algorithm> // for sort
// #include <stdexcept>

#include "simple-restraint.hh"

#include "coot-utils/coot-coord-extras.hh"  // is_nucleotide_by_dict

// also uses class variable non_bonded_neighbour_residues
void
coot::restraints_container_t::add_fixed_atoms_from_flanking_residues(const coot::bonded_pair_container_t &bpc) {

   std::vector<mmdb::Residue *> residues_for_fixed_atoms;
   
   for (unsigned int i=0; i<bpc.size(); i++) { 
      if (bpc[i].is_fixed_first)
	 residues_for_fixed_atoms.push_back(bpc[i].res_1);
      if (bpc[i].is_fixed_second)
	 residues_for_fixed_atoms.push_back(bpc[i].res_2);
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
	    if (std::find(fixed_atom_indices.begin(),
			  fixed_atom_indices.end(), idx) == fixed_atom_indices.end())
	       fixed_atom_indices.push_back(idx);
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
coot::restraints_container_t::bonded_flanking_residues_by_residue_vector(const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t bpc;
   float dist_crit = 3.0;

   // Don't forget to consider the case were we refine one residue,
   // and that is the ASN for a glycosylation.  So the neighbouring
   // residues and the GLC/NAG (say) will be flanking residues.
   //
   // But that is a really hard thing - isn't it? To find a GLC that
   // is connected to this residue?

   
   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {

      std::vector<mmdb::Residue *> neighbours = coot::residues_near_residue(residues_vec[ir].second,
								       mol, dist_crit);
//       std::cout << " DEBUG:: bonded_flanking_residues_by_residue_vector using reference res "
//  		<< ir << " of " << residues_vec.size() << " " 
//  		<< residues_vec[ir].second << " with " << neighbours.size() << " neighbours"
// 		<< std::endl;
      
      for (unsigned int ineighb=0; ineighb<neighbours.size(); ineighb++) {

	 bool found = 0;
	 for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	    // pointer comparison
	    if (neighbours[ineighb] == residues_vec[ires].second) {
	       found = 1;
	       break;
	    } 
	 }

	 if (! found) {
	    // OK, so this neighbour was not in the passed set of
	    // moving residues, it can be a flanking residue then...
	    std::pair<bool, float> d = closest_approach(neighbours[ineighb], residues_vec[ir].second);
// 	    std::cout << " not in residue vec... " << d.first << " " << d.second
// 		      << " for " << coot::residue_spec_t(neighbours[ineighb]) << " to "
// 		      << coot::residue_spec_t(residues_vec[ir].second) << std::endl;
	    if (d.first) {
	       if (d.second < dist_crit) {

		  std::pair<std::string, bool> l =
		     find_link_type_complicado(neighbours[ineighb], residues_vec[ir].second, geom);
		  std::string link_type = l.first;
		  if (link_type != "") {
		     bool order_switch_flag = l.second;
		     if (! order_switch_flag) {
			coot::bonded_pair_t bp(neighbours[ineighb], residues_vec[ir].second, 1, 0, link_type);
			bpc.try_add(bp);
		     } else {
			coot::bonded_pair_t bp(residues_vec[ir].second, neighbours[ineighb], 0, 1, link_type);

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
