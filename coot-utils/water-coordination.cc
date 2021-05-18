/* coot-utils/coot-coord-utils.hh
 * 
 * Copyright 2010 by The University of Oxford
 * Copyright 2015 by Medical Research Council
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

#include <algorithm>
#include "coot-coord-utils.hh"

coot::util::water_coordination_t::water_coordination_t(mmdb::Manager *mol, mmdb::realtype radius) {

   init_internal(mol, radius, false);
}


coot::util::water_coordination_t::water_coordination_t(mmdb::Manager *mol,
                                                       mmdb::realtype radius,
                                                       bool do_metals_only_flag) {

   init_internal(mol, radius, do_metals_only_flag);

}

void
coot::util::water_coordination_t::init_internal(mmdb::Manager *mol,
                                                mmdb::realtype radius,
                                                bool do_metals_only_flag) {

   if (! mol)
      return; 

   mmdb::realtype min_dist = 0.5;
   mmdb::realtype max_dist = radius;
   
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
         my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   mmdb::PAtom *water_selection = NULL;
   int n_water_atoms;
   mmdb::PAtom *atom_selection = NULL;
   int n_selected_atoms;
   int SelHnd =        mol->NewSelection(); // d
   int SelHnd_waters = mol->NewSelection(); // d

   if (! do_metals_only_flag) 
      mol->SelectAtoms (SelHnd_waters, 0, "*",
                        mmdb::ANY_RES, // starting resno, an int
                        "*", // any insertion code
                        mmdb::ANY_RES, // ending resno
                        "*", // ending insertion code
                        "HOH", // residue name
                        "*",   // atname
                        "*",   // elements
                        "*"    // alt loc.
                        );

   mol->SelectAtoms (SelHnd_waters, 0, "*",
                     mmdb::ANY_RES, // starting resno, an int
                     "*", // any insertion code
                     mmdb::ANY_RES, // ending resno
                     "*", // ending insertion code
                     "*", // residue name
                     "*",   // atname
                     "MG,CA,K,NA,LI,RB,BE,BA,FR,CS,SR,RA,SC.TI,V,CR,MN,FE,CO,NI,CU,ZN,ZR,NB,MO,RU,RH,Ag,Cd,W,OS,IR,PT,AU,HG",   // metals
                     "*" ,   // alt loc.
                     mmdb::SKEY_OR);
   

   // c.f. addSymmetry_whole_chain() in Bond-lines.cc 
   
   mol->SelectAtoms (SelHnd, 0, "*",
                     mmdb::ANY_RES, // starting resno, an int
                     "*", // any insertion code
                     mmdb::ANY_RES, // ending resno
                     "*", // ending insertion code
                     "*", // any residue name
                     "*",
                     "*", // elements
                     "*"  // alt loc.
                     );

   mol->GetSelIndex(SelHnd_waters, water_selection, n_water_atoms);
   mol->GetSelIndex(SelHnd,         atom_selection, n_selected_atoms);

   mmdb::mat44 test_mat;
   int i_symm_err = mol->GetTMatrix(test_mat, 0, 0, 0, 0);

   if (i_symm_err) { // no symmetry fallback
      add_contacts(mol, water_selection, n_water_atoms, atom_selection, n_selected_atoms,
                   min_dist, max_dist, my_matt);
   } else {
      for (int ix = -1; ix < 2; ix++) { 
         for (int iy = -1; iy < 2; iy++) { 
            for (int iz = -1; iz < 2; iz++) {
               for (int isym = 0; isym < mol->GetNumberOfSymOps(); isym++) { 
                  mol->GetTMatrix(my_matt, isym, ix, iy, iz);
                  add_contacts(mol, water_selection, n_water_atoms, atom_selection, n_selected_atoms,
                               min_dist, max_dist, my_matt);
               }
            }
         }
      }
   }
   mol->DeleteSelection(SelHnd);
   mol->DeleteSelection(SelHnd_waters);
   

} 



coot::util::contact_atoms_info_t::contact_atom_t::contact_atom_t(mmdb::Atom *contactor, mmdb::Atom *central_atom) {

   clipper::Coord_orth co_1(   contactor->x,    contactor->y,    contactor->z);
   clipper::Coord_orth co_2(central_atom->x, central_atom->y, central_atom->z);
   dist = clipper::Coord_orth::length(co_1, co_2);
   at = contactor;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++)
         mat[i][j] = 0;
   for (int i=0; i<4; i++) 
         mat[i][i] = 1;
}


coot::util::contact_atoms_info_t::contact_atom_t::contact_atom_t(mmdb::Atom *contactor,
                                                                 mmdb::Atom *central_atom,
                                                                 const mmdb::mat44 &m) {

   clipper::Coord_orth co_1(   contactor->x,    contactor->y,    contactor->z);
   clipper::Coord_orth co_2(central_atom->x, central_atom->y, central_atom->z);
   dist = clipper::Coord_orth::length(co_1, co_2);
   at = contactor;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++)
         mat[i][j] = m[i][j];
}


// But don't add them if they are not matching alt confs
// 
void
coot::util::water_coordination_t::add_contact(mmdb::Atom *atom_central,
                                              mmdb::Atom *atom_contactor,
                                              const mmdb::mat44 &mat) {

   std::string alt_conf_1 = atom_contactor->altLoc;
   std::string alt_conf_2 = atom_central->altLoc;

   if ((alt_conf_1 == alt_conf_2) || (alt_conf_1 == "") || (alt_conf_2 == "")) { 

      // filter out H water contacts.
      std::string ele(atom_contactor->element);
      if (ele != " H") { 
         coot::util::contact_atoms_info_t::contact_atom_t con_at(atom_contactor, atom_central, mat);
         //
         // Now where does con_at go?
         // is atom_central already in atom_contacts?
         //
         bool found = 0;
         for (unsigned int i=0; i<atom_contacts.size(); i++) {
            if (atom_contacts[i].matches_atom(atom_central)) {
               atom_contacts[i].add(con_at);
               found = 1; 
               break;
            } 
         }

         if (! found) {
            coot::util::contact_atoms_info_t cai(atom_central, con_at);
            atom_contacts.push_back(cai);
         }
      }
   }
}



std::vector<coot::util::contact_atoms_info_t>
coot::util::water_coordination_t::get_highly_coordinated_waters(int n_contacts,  // at least n_contacts
                                                                double dist_max) const { // within dist_max

   std::vector<coot::util::contact_atoms_info_t> v;
   //
   for (unsigned int i=0; i<atom_contacts.size(); i++) {
      if (atom_contacts[i].has_contacts(n_contacts, dist_max)) {
         v.push_back(atom_contacts[i]);
      }
   }
   sort_contacts(&v); // fiddle with v;
   return v;
}

std::vector<std::pair<coot::util::contact_atoms_info_t, coot::util::contact_atoms_info_t::ele_index_t> >
coot::util::water_coordination_t::metals() const {

   std::vector<std::pair<coot::util::contact_atoms_info_t, coot::util::contact_atoms_info_t::ele_index_t> > v;
   std::vector<coot::util::contact_atoms_info_t::ele_index_t> eles;

   // search for these in order, ie. if we find/hit a MG2, don't test same atom for NA.
   // 
   eles.push_back(coot::util::contact_atoms_info_t::ELE_CA2);
   eles.push_back(coot::util::contact_atoms_info_t::ELE_MG2);
   eles.push_back(coot::util::contact_atoms_info_t::ELE_NA);
   eles.push_back(coot::util::contact_atoms_info_t::ELE_K);
   eles.push_back(coot::util::contact_atoms_info_t::ELE_LI);
   
   for (unsigned int i=0; i<atom_contacts.size(); i++) {
      for (unsigned int j =0; j<eles.size(); j++) {
         if (atom_contacts[i].test_for_ele(eles[j])) {
            std::pair<coot::util::contact_atoms_info_t, coot::util::contact_atoms_info_t::ele_index_t>
               p(atom_contacts[i], eles[j]);
            v.push_back(p);
            break; // ele loop
         }
      }
   }


   // sort_contacts(&v); // fiddle with v;
   return v;
}

void
coot::util::water_coordination_t::transform_atom(int i, int j) {

   mmdb::mat44 &m = atom_contacts[i][j].mat;
   // atom_contacts[i][j].at->Transform(atom_contacts[i][j].mat);
   atom_contacts[i][j].at->Transform(m);

}




void
coot::util::water_coordination_t::sort_contacts(std::vector<coot::util::contact_atoms_info_t> *v) const {
   std::sort(v->begin(), v->end(), sort_contacts_func);
}

// static 
bool
coot::util::water_coordination_t::sort_contacts_func(const coot::util::contact_atoms_info_t &first,
                                                     const coot::util::contact_atoms_info_t &second) {
   try { 
      double smallest_dist_1 = first.smallest_contact_dist();
      double smallest_dist_2 = second.smallest_contact_dist();
      return (smallest_dist_1 < smallest_dist_2);
   }
   catch (const std::runtime_error &rte) {
      return 0;
   }
}

void
coot::util::water_coordination_t::add_contacts(mmdb::Manager *mol, 
                                               mmdb::PAtom *water_selection, int n_water_atoms, 
                                               mmdb::PAtom *atom_selection, int n_selected_atoms,
                                               mmdb::realtype min_dist, mmdb::realtype max_dist,
                                               const mmdb::mat44 &my_mat) {

   

   mmdb::Contact *pscontact = NULL;
   int n_contacts = 0;
   long i_contact_group = 1;
   mmdb::mat44 other_mat;

   // It's a bit grim that this is needed.  I can't see how to use
   // my_mat in SeekContacts() directly.
   // 
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++)
         other_mat[i][j] = my_mat[i][j];

   mol->SeekContacts(water_selection, n_water_atoms, 
                     atom_selection, n_selected_atoms,
                     min_dist, max_dist, // min, max distances
                     0,        // seqDist 0 -> in same res also
                     pscontact, n_contacts,
                     0, &other_mat, i_contact_group);

   if (n_contacts > 0) {
      for (int i=0; i< n_contacts; i++) {
         
         // e.g. the NE2 in atom_selection contacts O in an HOH in the water_selection
         // (the test for matching altconfs id done in add_contact().
         //
         add_contact(water_selection[pscontact[i].id1], atom_selection[pscontact[i].id2], other_mat);
      }
   } 
}

// Only test for Na like this if the water_coordination_t max_dist was
// 4.0A (and min_dist 0.1 or so).
// 
bool
coot::util::contact_atoms_info_t::test_for_na() const {

   return test_for_ele(coot::util::contact_atoms_info_t::ELE_NA);
   
}

// Only test for Na like this if the water_coordination_t max_dist was
// 4.0A (and min_dist 0.1 or so).
// 
bool
coot::util::contact_atoms_info_t::test_for_ele(coot::util::contact_atoms_info_t::ele_index_t ele_index) const {

   bool r = 0; // return this result (boolean hit)
   
   // For the calculation of V_{Na^+}
   double R0 = 1.622; // Brown and Wu, 1976; Nayal & Di Cera, 1996
   double  N = 4.29;
   double v_ele_target = 0.95; // 95% of +1 charge compensation (a bit of wiggle room)

   // More Brown and Wu values (Mg2+ same as Na)
   
   if (ele_index == coot::util::contact_atoms_info_t::ELE_K) {
      R0 = 2.276;
      N = 9.1;
   }

   if (ele_index == coot::util::contact_atoms_info_t::ELE_LI) {
      R0 = 1.378;
      N = 4.065;
   }

   if (ele_index == coot::util::contact_atoms_info_t::ELE_MG2) {
      v_ele_target *= 2.0; // double charge
   }

   if (ele_index == coot::util::contact_atoms_info_t::ELE_CA2) {
      R0 = 1.909;
      N = 5.4;
      v_ele_target *= 2.0; // double charge
   }


   double sum_v_j = 0.0;
   int n_protein_contacts = 0;
   int n_contacts = 0;
   if (contact_atoms.size() > 2) { 
      for (unsigned int j=0; j<contact_atoms.size(); j++) {
         std::string ele(contact_atoms[j].at->element);
         double Rj = contact_atoms[j].dist;
         if (Rj < R0) { // it's too close to a neighbor to be the required metal type.
            sum_v_j = 0; 
            break;
         }
         if (ele == " O") {
            double v_j = pow(Rj/R0, -N);
            double occ = contact_atoms[j].at->occupancy;  // symmetry overlaps handled.
            sum_v_j += v_j * occ;
            // was it a protein atom contact < 3.5A?  We need at least 1 of them.
            std::string resname(contact_atoms[j].at->GetResName());
            if (Rj < 3.5) { 
               if (resname != "HOH") {
                  n_protein_contacts++;
               }
               n_contacts++;
            }
         }
         if (ele == " N") {
            if (Rj < 3.4) { 
               sum_v_j = 0;
               break; // a metal does not have close nitrogens.
            }
         }
         if (ele == " C") {
            if (Rj < 2.9) {  // it was a bump.  Possibly/probably not even a good water position.
               sum_v_j = 0;
               break; // a metal does not have close carbons
            }
         }
         // and similarly, no metals next to metals
         if (ele == "NA"  || ele == " K" || ele == "CA" || ele == "MG") {
            if (Rj < 3.0) {
               sum_v_j = 0;
               break; // a metal does not have close metals.
            }
         }
      }
      if (n_protein_contacts > 1) { 
         if (n_contacts > 3) { 
            r = (sum_v_j > v_ele_target);
            if (r) {
               std::cout << "Hit for " << coot::atom_spec_t(central_atom()) << " "
                         << ele_index << std::endl;
               sum_v_j = 0.0;
               for (unsigned int j=0; j<contact_atoms.size(); j++) {
                  std::string ele(contact_atoms[j].at->element);
                  double Rj = contact_atoms[j].dist;
                  if (ele == " O") {
                     double v_j = pow(Rj/R0, -N);
                     sum_v_j += v_j;
                     std::cout << "   contribution " << v_j << " makes "<< sum_v_j << std::endl; 
                  }
                  if (ele == " N") {
                     if (Rj < 3.4) { 
                        sum_v_j = 0;
                        break; // a metal does not have close nitrogens.
                     }
                  }
               }
               std::cout << "   DEBUG:: got sum_v_j " << sum_v_j << " comparing to "
                         << v_ele_target << std::endl;
               std::cout << std::endl;
            }
         }
      }
   }
   return r;
}

