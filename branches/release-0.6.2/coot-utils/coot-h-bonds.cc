/* coot-utils/coot-h-bonds.cc
 * 
 * Copyright 2010 by The University of Oxford
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

#include "coot-h-bonds.hh"


// typically the atom selection selHnd_1 is for the ligand and
// selHnd_2 is for everything (else).
// 
std::vector<coot::h_bond>
coot::h_bonds::get(int selHnd_1, int selHnd_2, CMMDBManager *mol, const coot::protein_geometry &geom) {

   bool debug = 0;
   
   realtype min_dist = 2.4; // H-bonds are longer than this
   realtype max_dist = 3.9; // H-bonds are shorter than this
   
   std::vector<coot::h_bond> v;

   int hb_type_udd_handle = mark_donors_and_acceptors(selHnd_1, selHnd_2, mol, geom); // using UDD data

   // What is the nearest neighbour of the atoms in mol?
   // 
   std::map<CAtom *, std::vector<std::pair<CAtom *, float> > > neighbour_map =
      make_neighbour_map(selHnd_1, selHnd_2, mol);

   PPCAtom sel_1_atoms = 0;
   PPCAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mat44 my_matt;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   
   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);

   PSContact pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;

   mol->SeekContacts(sel_1_atoms, n_sel_1_atoms,
		     sel_2_atoms, n_sel_2_atoms,
		     min_dist, max_dist,
		     0, // seqDist 0 -> also in same res.
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   std::cout << "h_bonds:: get(): found n_contacts between atom selections: " << n_contacts << std::endl;

   if (n_contacts > 0) {
      if (pscontact) {
	 for (int i_contact=0; i_contact<n_contacts; i_contact++) {
	    CAtom *at_1 = sel_1_atoms[pscontact[i_contact].id1];
	    CAtom *at_2 = sel_2_atoms[pscontact[i_contact].id2];

	    // are they donor and acceptor?
	    //
	    int hb_type_1 = coot::energy_lib_atom::HB_UNASSIGNED;
	    int hb_type_2 = coot::energy_lib_atom::HB_UNASSIGNED;

	    // hb_type_1 = 999; // debugging
	    // hb_type_2 = 999;

	    at_1->GetUDData(hb_type_udd_handle, hb_type_1);
	    at_2->GetUDData(hb_type_udd_handle, hb_type_2);

	    if (0) { 
	       std::cout << "atom 1st selection: " << coot::atom_spec_t(at_1) << " has hb_type "
			 << hb_type_1 << std::endl;
	       std::cout << "atom 2nd selection: " << coot::atom_spec_t(at_2) << " has hb_type "
			 << hb_type_2 << std::endl;
	    }

	    bool match = 0;
	    bool swap = 0;
	    bool ligand_atom_is_donor = 1;

	    if ((hb_type_1 == coot::energy_lib_atom::HB_DONOR ||
		 hb_type_1 == coot::energy_lib_atom::HB_BOTH) &&
		(hb_type_2 == coot::energy_lib_atom::HB_ACCEPTOR ||
		 hb_type_2 == coot::energy_lib_atom::HB_BOTH)) {
	       match = 1;
	    }

	    if ((hb_type_1 == coot::energy_lib_atom::HB_ACCEPTOR ||
		 hb_type_1 == coot::energy_lib_atom::HB_BOTH) &&
		(hb_type_2 == coot::energy_lib_atom::HB_DONOR ||
		 hb_type_2 == coot::energy_lib_atom::HB_BOTH)) {
	       match = 1;
	       swap = 1;
	       ligand_atom_is_donor = 0;
	    }

	    if (match) {

	       std::cout << "MATCH (pre-swap) H-bond donor acceptor match "
			 << coot::atom_spec_t(at_1) << " " << at_1->GetResName() 
			 << " to " << coot::atom_spec_t(at_2) << " " << at_2->GetResName() << "\n";
	       
	       if (swap) { 
		  std::swap(at_1, at_2);
	       }

	       // donor first: xxx_1
	       // 
	       std::vector<std::pair<CAtom *, float> > nm_1 = neighbour_map[at_1];
	       std::vector<std::pair<CAtom *, float> > nm_2 = neighbour_map[at_2];
	       
	       std::string res_type_1 = at_1->GetResName();
	       std::string res_type_2 = at_2->GetResName();

	       double angle_1_0 = -1; // assigned later
	       double angle_2_0 = -1; // assigned later
	       
	       std::cout << "   Here 1 #donor-neighbs: " << nm_1.size()
			 << " #acceptor-neighbs: " << nm_2.size() << std::endl;
	       
	       if (((nm_1.size() > 0) || (res_type_1 == "HOH")) &&
		   ((nm_2.size() > 0) || (res_type_2 == "HOH"))) {

		  // only check the angles if there is not already a
		  // bond in the list between these atoms.
		  // 
		  coot::h_bond hb(at_1, at_2);
		  std::vector<coot::h_bond>::const_iterator it = std::find(v.begin(), v.end(), hb);
		  if (it == v.end()) {
		     
		     bool neighbour_angles_are_good = 1;  // initially
		     double dist = coot::distance(at_1, at_2);
		     for (unsigned int ii=0; ii<nm_1.size(); ii++) {

			double angle_1 = coot::angle(nm_1[ii].first, at_1, at_2);

			if (angle_1_0 < 0) {  // as yet unset
			   coot::residue_spec_t donor_res_spec(at_1);
			   coot::residue_spec_t neigh_res_spec(nm_1[ii].first);
			   if (donor_res_spec == neigh_res_spec)
			      angle_1_0 = angle_1; 
			}
		     
			if (angle_1 < 90) {
			   neighbour_angles_are_good = 0;
			   if (debug)
			      std::cout << "Rejecting (donor) "
					<< coot::atom_spec_t(at_1) << "..."
					<< coot::atom_spec_t(at_2) << "  because angle "
					<< coot::atom_spec_t(nm_1[ii].first) << "-"
					<< coot::atom_spec_t(at_1) << "-"
					<< coot::atom_spec_t(at_2) << " is " << angle_1
					<< "\n";

			   break;
			}
		     }

		     // We don't loop the second angle.  Is that a bad thing?

		     double angle_2 = -1;
		     if (nm_2.size()) // not a HOH
			angle_2 = coot::angle(at_1, at_2, nm_2[0].first);
		     
		     if (angle_2_0 < 0) {  // as yet unset
			angle_2_0 = angle_2;
		     }

		     // Reject if bad angles (but don't reject if the
		     // bonded residue is a HOH (in that case there is
		     // no angle_2)).
		     // 
		     if ((angle_2 < 90) && ((res_type_1 != "HOH") && (res_type_2 != "HOH"))) {
			neighbour_angles_are_good = 0;
			if (debug) 
			   std::cout << "Rejecting (acceptor) "
				     << coot::atom_spec_t(at_1) << "..."
				     << coot::atom_spec_t(at_2) << "  because angle "
				     << coot::atom_spec_t(at_1) << "-"
				     << coot::atom_spec_t(at_2)
				     << coot::atom_spec_t(nm_2[0].first)
				     << " is " << angle_2
				     << "\n";
		     }
		  
		     std::cout << "   Here 4 angles_good: " << neighbour_angles_are_good << std::endl;
		     
		     if (neighbour_angles_are_good) {

			if (nm_1.size())
			   hb.donor_neigh = nm_1[0].first;
			if (nm_2.size())
			   hb.acceptor_neigh = nm_2[0].first;
			hb.dist = dist;
			hb.angle_1 = angle_1_0;
			hb.angle_2 = angle_2_0;
			hb.ligand_atom_is_donor = ligand_atom_is_donor; // selHnd_1 is presumed "ligand"
			if (1) 
			   std::cout << "Adding bond " << coot::atom_spec_t(hb.donor) << " to "
				     << coot::atom_spec_t(hb.acceptor) << std::endl;
			v.push_back(hb);
		     }
		  }
	       }
	    } else {
	       // std::cout << "No H-bond donor acceptor match " << coot::atom_spec_t(at_1)
	       // << " to " << coot::atom_spec_t(at_2) << "\n";
	    } 
	 }
      }
   }

   std::sort(v.begin(), v.end());
   return v;
}


// using UDD data
// 
// return the UDD handle
int
coot::h_bonds::mark_donors_and_acceptors(int selHnd_1, int selHnd_2, CMMDBManager *mol,
					 const coot::protein_geometry &geom) {

   bool debug = 0;
   PPCAtom sel_1_atoms = 0;
   PPCAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);
   int udd_h_bond_type_handle = mol->RegisterUDInteger(UDR_ATOM, "hb_type");

   for (unsigned int i=0; i<n_sel_1_atoms; i++) { 
      std::string name = sel_1_atoms[i]->name;
      std::string res_name = sel_1_atoms[i]->GetResName();
      int h_bond_type = geom.get_h_bond_type(name, res_name);
      sel_1_atoms[i]->PutUDData(udd_h_bond_type_handle, h_bond_type);
      if (debug)
	 std::cout << sel_1_atoms[i]->GetChainID() << " "
		   << sel_1_atoms[i]->GetSeqNum() << " "
		   << sel_1_atoms[i]->GetResName() << " "
		   << "name: " << name << " marked as " << h_bond_type << "\n";
   }

   if (selHnd_1 != selHnd_2) {
      for (unsigned int i=0; i<n_sel_2_atoms; i++) { 
	 std::string name = sel_2_atoms[i]->name;
	 std::string res_name = sel_2_atoms[i]->GetResName();
	 int h_bond_type = geom.get_h_bond_type(name, res_name);
	 sel_2_atoms[i]->PutUDData(udd_h_bond_type_handle, h_bond_type);
      if (debug)
	 std::cout << sel_2_atoms[i]->GetChainID() << " "
		   << sel_2_atoms[i]->GetSeqNum() << " "
		   << sel_2_atoms[i]->GetResName() << " "
		   << "name: " << name << " marked as " << h_bond_type << "\n";
      }
   }

   return udd_h_bond_type_handle;
}


// What is the nearest neighbour of the atoms in mol?
// 
std::map<CAtom *, std::vector<std::pair<CAtom *, float> > >
coot::h_bonds::make_neighbour_map(int selHnd_1, int selHnd_2, CMMDBManager *mol) {

   std::map<CAtom *, std::vector<std::pair<CAtom *, float> > > atom_map;
   PPCAtom sel_1_atoms = 0;
   PPCAtom sel_2_atoms = 0;
   int n_sel_1_atoms;
   int n_sel_2_atoms;
   mol->GetSelIndex   (selHnd_1, sel_1_atoms, n_sel_1_atoms);
   mol->GetSelIndex   (selHnd_2, sel_2_atoms, n_sel_2_atoms);

   mat44 my_matt;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   
   PSContact pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;

   mol->SeekContacts(sel_1_atoms, n_sel_1_atoms,
		     sel_1_atoms, n_sel_1_atoms,
		     0.1, 1.7,
		     0, // seqDist 0 -> also in same res.
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   if (n_contacts) {
      if (pscontact) {
	 for (int i_contact=0; i_contact<n_contacts; i_contact++) {
	    clipper::Coord_orth pt_1(sel_1_atoms[pscontact[i_contact].id1]->x,
				     sel_1_atoms[pscontact[i_contact].id1]->y,
				     sel_1_atoms[pscontact[i_contact].id1]->z);
	    clipper::Coord_orth pt_2(sel_1_atoms[pscontact[i_contact].id2]->x,
				     sel_1_atoms[pscontact[i_contact].id2]->y,
				     sel_1_atoms[pscontact[i_contact].id2]->z);
	    float d = clipper::Coord_orth::length(pt_1, pt_2);
	    coot::residue_spec_t res_1(sel_1_atoms[pscontact[i_contact].id1]->GetResidue());
	    coot::residue_spec_t res_2(sel_1_atoms[pscontact[i_contact].id2]->GetResidue());
	    
	    std::pair<CAtom *, float> p(sel_1_atoms[pscontact[i_contact].id2], d);
	    atom_map[sel_1_atoms[pscontact[i_contact].id1]].push_back(p);
	 }
      }
   }
   delete [] pscontact;
   pscontact = NULL;


   mol->SeekContacts(sel_2_atoms, n_sel_2_atoms,
		     sel_2_atoms, n_sel_2_atoms,
		     0.1, 1.7,
		     0, // seqDist 0 -> also in same res.
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   if (n_contacts) {
      if (pscontact) {
	 for (int i_contact=0; i_contact<n_contacts; i_contact++) {
	    clipper::Coord_orth pt_1(sel_2_atoms[pscontact[i_contact].id1]->x,
				     sel_2_atoms[pscontact[i_contact].id1]->y,
				     sel_2_atoms[pscontact[i_contact].id1]->z);
	    clipper::Coord_orth pt_2(sel_2_atoms[pscontact[i_contact].id2]->x,
				     sel_2_atoms[pscontact[i_contact].id2]->y,
				     sel_2_atoms[pscontact[i_contact].id2]->z);
	    float d = clipper::Coord_orth::length(pt_1, pt_2);
	    coot::residue_spec_t res_1(sel_2_atoms[pscontact[i_contact].id1]->GetResidue());
	    coot::residue_spec_t res_2(sel_2_atoms[pscontact[i_contact].id2]->GetResidue());
	    
	    std::pair<CAtom *, float> p(sel_2_atoms[pscontact[i_contact].id2], d);

	    // only add p if is not already in the atom map vector for this atom:
	    // (this relies on the doubles matching :) but it seems to work...
	    // 
	    std::vector<std::pair<CAtom *, float> >::const_iterator it = 
	       std::find(atom_map[sel_2_atoms[pscontact[i_contact].id1]].begin(),
			 atom_map[sel_2_atoms[pscontact[i_contact].id1]].end(), p);
	    if (it == atom_map[sel_2_atoms[pscontact[i_contact].id1]].end())
	       atom_map[sel_2_atoms[pscontact[i_contact].id1]].push_back(p);
	 }
      }
   }
   delete [] pscontact;
   pscontact = NULL;
   

   // sort the neighbour by distance 
   //
   std::map<CAtom *, std::vector<std::pair<CAtom *, float> > >::iterator it;
   for (it=atom_map.begin(); it != atom_map.end(); it++) {

      //  You can't use a const_iterator to sort this vector's innards.
      // 
      coot::h_bonds::atom_sorter as(it->first);
      std::sort(it->second.begin(), it->second.end(), as);
   }

   bool debug = 0;

   // were they sorted correctly?  Debug
   //
   if (debug) {
      for (it=atom_map.begin(); it != atom_map.end(); it++) {
	 std::cout << coot::atom_spec_t(it->first) << "\n";
	 for (unsigned int i=0; i<it->second.size(); i++) { 
	    std::cout << "    DEUBG:: " << coot::atom_spec_t(it->second[i].first) << " "
		      << it->second[i].second << "\n";
	 }
      }
   }

   return atom_map;
}

