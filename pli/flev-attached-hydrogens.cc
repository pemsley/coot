/* pli/flev-attached-hydrogens.cc
 * 
 * Copyright 2013 by Medical Research Council
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

#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "flev-annotations.hh" // break out flev_attached_hydrogens_t from here?

// examine the dictionary and find the atoms to which the hydrogens
// are attached.  Is the hydrogen riding or rotatable?
// 
//            ... riding hydrogens do not have torsions for them in
// the dictionary.  Rotable hydrogens do.  Note though, that PRODRG
// hydrogen names (in the cif file) are problematic.  You can't rely
// on name matching of hydrogens between the PDB and the dictionary.
//
// So here is what we have to do:
//
// i) From the dictionary, look through the list of non-Hydrogen atoms
// to see if there are hydrogens bonded to the atom.
// 
//    ii) If so, does this non-Hydrogen atom appear as atom 2 in a 
//    torsion restraint with a hydrogen as atom 1 or appear as atom 3
//    in a torsion restraint with a hydrogen as atom 4?
//
//        iii) If so, then this is a torsionable hydrogen, we need to
//        generate vectors by random sampling from the probability
//        distribution of this torsion.
//
//        However, if it is not, then this is a riding hydrogen.  Add
//        the vector from the ligand atom to the hydrogen as a
//        cannonball vector associated with that ligand atom.
//
pli::flev_attached_hydrogens_t::flev_attached_hydrogens_t(const coot::dictionary_residue_restraints_t &restraints) {

   for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
      std::string atom_name_1 = restraints.bond_restraint[ibond].atom_id_1_4c();
      std::string atom_name_2 = restraints.bond_restraint[ibond].atom_id_2_4c();
      if ((restraints.is_hydrogen(atom_name_1)) && (! restraints.is_hydrogen(atom_name_2))) {
	 std::swap(atom_name_1, atom_name_2);
      }
      
      if ((! restraints.is_hydrogen(atom_name_1)) && (restraints.is_hydrogen(atom_name_2))) {
	 // a heavy atom connected to a hydrogen.
	 // Does it exist in the torsions?

	 bool found = 0;
	 std::pair<std::string, std::string> p(atom_name_1, atom_name_2);
	 for (unsigned int itor=0; itor<restraints.torsion_restraint.size(); itor++) { 
	    if (((restraints.torsion_restraint[itor].atom_id_1_4c() == atom_name_2) &&
		 (restraints.torsion_restraint[itor].atom_id_2_4c() == atom_name_1)) ||
		((restraints.torsion_restraint[itor].atom_id_4_4c() == atom_name_2) &&
		 (restraints.torsion_restraint[itor].atom_id_3_4c() == atom_name_1))) { 
	       if (restraints.torsion_restraint[itor].is_const()) {
		  atoms_with_riding_hydrogens.push_back(p);
	       } else {
		  atoms_with_rotating_hydrogens.push_back(p);
	       }
	       found = 1;
	       break;
	    }
	 }
	 if (! found) {
	    atoms_with_riding_hydrogens.push_back(p);
	 }
      }
   }
}

void
pli::flev_attached_hydrogens_t::cannonballs(mmdb::Residue *ligand_residue_3d,
					     const std::string &prodrg_3d_ligand_file_name,
					     const coot::dictionary_residue_restraints_t &restraints) {

   atom_selection_container_t asc = get_atom_selection(prodrg_3d_ligand_file_name, false, true, false);
   if (asc.read_success) {
      cannonballs(ligand_residue_3d, asc.mol, restraints);
   }
}

void
pli::flev_attached_hydrogens_t::cannonballs(mmdb::Residue *ligand_residue_3d,
					     mmdb::Manager *mol, 
					     const coot::dictionary_residue_restraints_t &restraints) {

   if (! mol)
      return;
   
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;


   int SelHnd_H = mol->NewSelection();
   int SelHnd_non_H = mol->NewSelection();

   mmdb::PPAtom hydrogen_selection = 0;
   mmdb::PPAtom non_hydrogen_selection = 0;
   int n_hydrogen_atoms;
   int n_non_hydrogen_atoms;
      
      
   mol->SelectAtoms(SelHnd_H,     0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", " H", "*");
   mol->SelectAtoms(SelHnd_non_H, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "!H", "*");
      
   mol->GetSelIndex(SelHnd_H, hydrogen_selection, n_hydrogen_atoms);
   mol->GetSelIndex(SelHnd_non_H, non_hydrogen_selection, n_non_hydrogen_atoms);
      
   std::cout << "Found " << n_hydrogen_atoms << " Hydrogens " << std::endl;
   std::cout << "Found " << n_non_hydrogen_atoms << " non Hydrogens " << std::endl;

   if (n_hydrogen_atoms == 0) {
      std::cout << "WARNING:: Oops found no hydrogens for cannonballs" << std::endl;
      return;
   }
   if (n_non_hydrogen_atoms == 0) {
      std::cout << "WARNING:: Oops found no non-hydrogens for cannonballs" << std::endl;
      return;
   }

   mol->SeekContacts(hydrogen_selection, n_hydrogen_atoms,
		     non_hydrogen_selection, n_non_hydrogen_atoms,
		     0.1, 1.5,
		     0, // in same res also
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   std::cout << "Found " << n_contacts << " contacts to Hydrogens " << std::endl;

   // We need to find the torsion of the hydrogen,
   // A torsion (that can be mapped to the reference ligand) is:
   //
   // At_name_base At_name_2 At_name_bond_to_H bond_length bond_angle torsion_angle
   //
   // that is, we work from "inside" ligand atoms out to the hydrogen
   // 
   // 
   if (n_contacts > 0) {
      for (int i=0; i< n_contacts; i++) {
	 mmdb::Atom *at = non_hydrogen_selection[pscontact[i].id2];
	 std::string atom_name_bonded_to_H(at->name);

	 bool found_torsion_for_this_H = 0;

	 // riding hydrogens:
	 // 
	 for (unsigned int iat=0; iat<atoms_with_riding_hydrogens.size(); iat++) { 
	    if (atom_name_bonded_to_H == atoms_with_riding_hydrogens[iat].first) {
	       mmdb::Atom *h_at = hydrogen_selection[pscontact[i].id1];
	       found_torsion_for_this_H = add_named_torsion(h_at, at, restraints, mol, H_IS_RIDING);
	    }
	    if (found_torsion_for_this_H)
	       break;
	 }

	 // rotating hydrogens:
	 // 
	 for (unsigned int iat=0; iat<atoms_with_rotating_hydrogens.size(); iat++) { 
	    if (atom_name_bonded_to_H == atoms_with_rotating_hydrogens[iat].first) {
	       mmdb::Atom *h_at = hydrogen_selection[pscontact[i].id1];
	       found_torsion_for_this_H = add_named_torsion(h_at, at, restraints, mol, H_IS_ROTATABLE);
	    }
	    if (found_torsion_for_this_H)
	       break;
	 }
      }
   }

   mol->DeleteSelection(SelHnd_H);
   mol->DeleteSelection(SelHnd_non_H);

   named_hydrogens_to_reference_ligand(ligand_residue_3d, restraints);

}


// RDKit version (well, the version that is used when the hydrogens
// are correctly named (according to the dictionary) and placed on the
// ligand of interest.
//
// This fills the atom_bashes vector.
//
void
pli::flev_attached_hydrogens_t::distances_to_protein_using_correct_Hs(mmdb::Residue *ligand_residue,
                                                                      mmdb::Manager *mol,
                                                                      const coot::protein_geometry &geom) {

   // the constructor (called just before this) should fill
   // atoms_with_rotating_hydrogens and atoms_with_riding_hydrogens
   // vectors (using the restraints).

   float radius = 6.0;
   std::vector<mmdb::Residue *> env_residues =
      coot::residues_near_residue(ligand_residue, mol, radius);

   // -------------------------------------------------------------------
   //                    riding hydrogens
   // -------------------------------------------------------------------
   //
   mmdb::PPAtom residue_atoms = 0;
   int n_ligand_atoms;
   ligand_residue->GetAtomTable(residue_atoms, n_ligand_atoms);
   for (unsigned int irh=0; irh<atoms_with_riding_hydrogens.size(); irh++) {
      mmdb::Atom *lig_at = NULL;
      mmdb::Atom *H_at = NULL;
      for (int iat=0; iat<n_ligand_atoms; iat++) {
         std::string atom_name(residue_atoms[iat]->name);
         if (atom_name == atoms_with_riding_hydrogens[irh].first)
            lig_at = residue_atoms[iat];
         if (atom_name == atoms_with_riding_hydrogens[irh].second)
            H_at = residue_atoms[iat];
         if (lig_at && H_at)
            break;
      }
      if (lig_at && H_at) {
         clipper::Coord_orth H_pt(H_at->x, H_at->y, H_at->z);
         clipper::Coord_orth lig_atom_pt(lig_at->x, lig_at->y, lig_at->z);

         std::vector<mmdb::Atom *> atoms = close_atoms(H_pt, env_residues);
         coot::bash_distance_t bash = find_bash_distance(lig_atom_pt, H_pt, atoms);
         atom_bashes[atoms_with_riding_hydrogens[irh].first].push_back(bash);
         if (true)
            std::cout << " adding bash distance " << bash << " to atom "
                      << atoms_with_riding_hydrogens[irh].first << std::endl;
      }
   }

   // -------------------------------------------------------------------
   //                 rotatable hydrogens (more complex)
   // -------------------------------------------------------------------
   //

   for (unsigned int irh=0; irh<atoms_with_rotating_hydrogens.size(); irh++) {
      mmdb::Atom *lig_at = NULL;
      mmdb::Atom *H_at = NULL;
      for (int iat=0; iat<n_ligand_atoms; iat++) {
         std::string atom_name(residue_atoms[iat]->name);
         if (atom_name == atoms_with_rotating_hydrogens[irh].first)
            lig_at = residue_atoms[iat];
         if (atom_name == atoms_with_rotating_hydrogens[irh].second)
            H_at = residue_atoms[iat];
         if (lig_at && H_at)
            break;
      }
      if (lig_at && H_at) {
         clipper::Coord_orth H_pt(H_at->x, H_at->y, H_at->z);
         clipper::Coord_orth lig_atom_pt(lig_at->x, lig_at->y, lig_at->z);

         std::vector<mmdb::Atom *> atoms = close_atoms(H_pt, env_residues);

         try {
            clipper::Coord_orth vector_pt = get_atom_pos_bonded_to_atom(lig_at, H_at, // not H_at
                                                                        ligand_residue, geom);
            clipper::Coord_orth base_ref_pt(0,0,0);
            double tors = clipper::Coord_orth::torsion(base_ref_pt, vector_pt, lig_atom_pt, H_pt);
            double dist = sqrt((lig_atom_pt - H_pt).lengthsq());
            double angle = clipper::Coord_orth::angle(vector_pt, lig_atom_pt, H_pt);

            int n_torsion_samples = 8;
            for (int itor=0; itor<n_torsion_samples; itor++) {

               double tmp_tor_d =  double(itor) * 360.0/double(n_torsion_samples);
               double tmp_tor = clipper::Util::d2rad(tmp_tor_d);
               tmp_tor += tors;
               clipper::Coord_orth new_pt =
                  clipper::Coord_orth(base_ref_pt, vector_pt, lig_atom_pt, dist, angle, tmp_tor);
               coot::bash_distance_t bash = find_bash_distance(lig_atom_pt, new_pt, atoms);
               atom_bashes[atoms_with_rotating_hydrogens[irh].first].push_back(bash);
            }
         }
         catch (const std::runtime_error &rte) {
            std::cout << rte.what() << std::endl;
         }
      }
   }
}

// apply those cannonball directions onto the real reference ligand:
void
pli::flev_attached_hydrogens_t::distances_to_protein(mmdb::Residue *residue_reference, 
                                                      mmdb::Manager *mol_reference) {

   float radius = 6.0;
   std::vector<mmdb::Residue *> env_residues =
      coot::residues_near_residue(residue_reference, mol_reference, radius);
   bool debug = false;

   if (debug) { 
      std::cout << "named torsion: " << named_torsions.size() << std::endl;
      for (unsigned int i=0; i<named_torsions.size(); i++) {
         std::cout << "  " << i
                   << named_torsions[i].base_atom_name << "  "
                   << named_torsions[i].atom_name_2 << "  "
                   << named_torsions[i].atom_name_bonded_to_H << "  "
                   << named_torsions[i].dist << "  "
                   << named_torsions[i].angle << "  "
                   << named_torsions[i].torsion << "  "
                   << std::endl;
      }
   }

   for (unsigned int i=0; i<named_torsions.size(); i++) {
      try {
         if (named_torsions[i].hydrogen_type == H_IS_RIDING) {

            if (debug)
               std::cout << "hydrogen on " << named_torsions[i].atom_name_bonded_to_H
                         << " is riding " << std::endl;

            std::pair<clipper::Coord_orth, clipper::Coord_orth> pt_base_and_H =
               hydrogen_pos(named_torsions[i], residue_reference);

            if (debug)
               std::cout << "pt_base_and_H: ligand atom at: " << pt_base_and_H.first.format()
                         << " H atom at: " << pt_base_and_H.second.format() << std::endl;

            std::vector<mmdb::Atom *> atoms = close_atoms(pt_base_and_H.second, env_residues);

            // bash is one of a (potentially) number of bash distances
            // for a given ligand (non-hydrogen) atom.
            // (passing the H coord, the lig-atom coord, a vector of mmdb::Atoms *s.)
            //
            coot::bash_distance_t bash = find_bash_distance(pt_base_and_H.first,
                                                            pt_base_and_H.second,
                                                            atoms);
            if (debug)
               std::cout << "   found bash: " << bash << std::endl;
            atom_bashes[named_torsions[i].atom_name_bonded_to_H].push_back(bash);
         }

         if (named_torsions[i].hydrogen_type == H_IS_ROTATABLE) {

            if (debug)
               std::cout << "hydrogen on " << named_torsions[i].atom_name_bonded_to_H
                         << " is rotatable " << std::endl;

            int n_torsion_samples = 8;
            std::pair<clipper::Coord_orth, clipper::Coord_orth> pt_base_and_H =
               hydrogen_pos(named_torsions[i], residue_reference);
            std::vector<mmdb::Atom *> atoms = close_atoms(pt_base_and_H.second, env_residues);
            for (int itor=0; itor<n_torsion_samples; itor++) {

               pli::named_torsion_t tmp_tor = named_torsions[i];
               tmp_tor.torsion += double(itor) * 360.0/double(n_torsion_samples);
               if (tmp_tor.torsion > 360)
                  tmp_tor.torsion -= 360;
               pt_base_and_H = hydrogen_pos(tmp_tor, residue_reference);

               coot::bash_distance_t bash = find_bash_distance(pt_base_and_H.first,
                                                               pt_base_and_H.second,
                                                               atoms);
               if (debug)
                  std::cout << "Adding bash distance " << bash << " to atom "
                            << named_torsions[i].atom_name_bonded_to_H
                            << std::endl;
               atom_bashes[named_torsions[i].atom_name_bonded_to_H].push_back(bash);
            }
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << rte.what() << std::endl;
      }
   }
}




coot::bash_distance_t
pli::flev_attached_hydrogens_t::find_bash_distance(const clipper::Coord_orth &ligand_atom_pos,
                                                    const clipper::Coord_orth &hydrogen_pos,
                                                    const std::vector<mmdb::Atom *> &close_residue_atoms) const {

   // find the residue from the close residue atoms and cache the
   // dictionaries here so that we can then call
   // dictionary_map[residue_p].type_energy(atom_name) and use that
   // type energy (checking for not "") to find if the atom is a
   // hydrogen bond accetor (or both) to use
   // energy_lib_t::some_new_accessor_function_hb_type(type_energy).
   // If hb_type is acceptor, then decrease bash distance,
   // atom_radius_plus_cbr by 0.8A or so.
   //
   std::map<mmdb::Residue *, coot::dictionary_residue_restraints_t> dictionary_map;

   double cannonball_radius = 0.8; // radius of the cannonball, c.f. at least a hydrogen.

   double max_dist = 4.05; // if we have travelled 4A without bashing
                           // into a protein atom then this has
                           // practically unlimited substitution distance.

   double max_dist_squared = max_dist * max_dist;
   clipper::Coord_orth h_vector((hydrogen_pos - ligand_atom_pos).unit());

   if (0)
      std::cout << "h_vector: " << h_vector.format() << " from hydrogen pos: "
                << hydrogen_pos.format() << "and ligand atom pos: " << ligand_atom_pos.format()
                << std::endl;

   // set the atomic radii:
   //
   std::vector<double> radius(close_residue_atoms.size());
   for (unsigned int iat=0; iat<close_residue_atoms.size(); iat++) {
      std::string ele(close_residue_atoms[iat]->element);
      radius[iat] = get_radius(ele);
   }

   coot::bash_distance_t dd;

   std::vector<clipper::Coord_orth> atom_positions(close_residue_atoms.size());
   // likewise set the atom positions so that we don't have to keep doing it.
   for (unsigned int i=0; i<close_residue_atoms.size(); i++)
      atom_positions[i] = clipper::Coord_orth(close_residue_atoms[i]->x,
                                              close_residue_atoms[i]->y,
                                              close_residue_atoms[i]->z);

   for (double slide=0; slide<=max_dist; slide+=0.04) {
      clipper::Coord_orth test_pt = ligand_atom_pos + slide * h_vector;
      if (true)
         std::cout << "   bash distance for ligand atom at " << ligand_atom_pos.format() << " "
                   << "determined from " << atom_positions.size() << " atom positions"
                   << std::endl;
      for (unsigned int iat=0; iat<atom_positions.size(); iat++) {
         double atom_radius_plus_cbr = radius[iat] + cannonball_radius;
         double d_squared = (test_pt - atom_positions[iat]).lengthsq();
         if (true)
            std::cout << "   atom " << iat << " "
                      << close_residue_atoms[iat]->GetChainID() << " "
                      << close_residue_atoms[iat]->GetSeqNum() << " "
                      << close_residue_atoms[iat]->GetAtomName() << " "
                      << " slide: " << slide
                      << " comparing " << sqrt(d_squared) << "^2  and "
                    << atom_radius_plus_cbr << "^2" << std::endl;
         if (d_squared < atom_radius_plus_cbr*atom_radius_plus_cbr) {
            dd = coot::bash_distance_t(slide);
            break;
         }
      }
      if (dd.limited) {
         break;
      }
   }
   return dd;
}



// What are the atoms that are close (distance < 6A) to pt?
//
// waters are not counted as close atoms.
//
std::vector<mmdb::Atom *>
pli::flev_attached_hydrogens_t::close_atoms(const clipper::Coord_orth &pt,
                                            const std::vector<mmdb::Residue *> &env_residues) const {

   std::vector<mmdb::Atom *> v;
   double dist_crit = 6.0;
   double dist_crit_squared = dist_crit * dist_crit;

   for (unsigned int i=0; i<env_residues.size(); i++) {
      mmdb::Residue *residue_p = env_residues[i];
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         clipper::Coord_orth atom_pos(residue_atoms[iat]->x, residue_atoms[iat]->y, residue_atoms[iat]->z);
         double d_squared = (pt - atom_pos).lengthsq();
         if (d_squared < dist_crit_squared) {
            std::string rn(residue_atoms[iat]->GetResName());
            if (rn != "HOH")
               v.push_back(residue_atoms[iat]);
         }
      }
   }
   return v;
}


double
pli::flev_attached_hydrogens_t::get_radius(const std::string &ele) const {

   double radius = 1.70;
   if (ele == " H")
      radius = 1.20;
   if (ele == " N")
      radius = 1.55;
   if (ele == " O")
      radius = 1.52;
   if (ele == " S")
      radius = 1.8;
   return radius;
}

// find an atom (the atom, perhaps) bonded to lig_at that is not H_at.
// Return its position.
//
// Can throw a std::runtime_error if not found.
//
clipper::Coord_orth
pli::flev_attached_hydrogens_t::get_atom_pos_bonded_to_atom(mmdb::Atom *lig_at, mmdb::Atom *H_at, // not H_at
                                                            mmdb::Residue *ligand_residue,
                                                            const coot::protein_geometry &geom) const {
   int imol = 0; // FIXME needs checking
   std::string res_name(lig_at->residue->GetResName());
   std::pair<bool, coot::dictionary_residue_restraints_t> p =
      geom.get_monomer_restraints_at_least_minimal(res_name, imol);

   if (! p.first) {
      std::string m = "No monomer type ";
      m += res_name;
      m += " found in dictionary";
      throw std::runtime_error(m);
   } else {
      mmdb::Atom *bonded_atom = NULL;
      std::string bonded_atom_name;
      std::string lig_at_name = lig_at->name;
      std::string H_at_name = H_at->name;
      for (unsigned int ibond=0; ibond<p.second.bond_restraint.size(); ibond++) {
         std::string atom_name_1 = p.second.bond_restraint[ibond].atom_id_1_4c();
         std::string atom_name_2 = p.second.bond_restraint[ibond].atom_id_2_4c();
         if (atom_name_1 == lig_at_name) {
            if (atom_name_2 != H_at_name) {
               bonded_atom_name = atom_name_2;
               break;
            }
         }
         if (atom_name_2 == lig_at_name) {
            if (atom_name_1 != H_at_name) {
               bonded_atom_name = atom_name_1;
               break;
            }
         }
      }
      if (bonded_atom_name != "") {
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         ligand_residue->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms ; iat++) {
            std::string atom_name = residue_atoms[iat]->name;
            if (atom_name == bonded_atom_name) {
               bonded_atom = residue_atoms[iat];
               break;
            }
         }
      }

      if (! bonded_atom) {
         std::string m = "No atom bonded to ";
         m += lig_at_name;
         m += " found in dictionary for ";
         m += res_name;
         throw std::runtime_error(m);
      } else {
         // good
         return clipper::Coord_orth(bonded_atom->x, bonded_atom->y, bonded_atom->z);
      }
   }

}

// For (each?) of the atoms in our real reference residue
// ligand_residue_3d (that should have hydrogens attached) give us a
// unit vector from the bonding atom in the direction of (each of, if
// there are more than one) the hydrogen(s).
// 
std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > >
pli::flev_attached_hydrogens_t::named_hydrogens_to_reference_ligand(mmdb::Residue *ligand_residue_3d,
								     const coot::dictionary_residue_restraints_t &restraints) const {

   std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > > v;

   for (unsigned int i=0; i<named_torsions.size(); i++) { 
      if (named_torsions[i].hydrogen_type == H_IS_RIDING) {
	 mmdb::Atom *atom_base = NULL;
	 mmdb::Atom *atom_2 = NULL;
	 mmdb::Atom *atom_bonded_to_H = NULL;
	 
	 mmdb::PPAtom residue_atoms = 0;
	 int n_residue_atoms;
	 ligand_residue_3d->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) { 
	    std::string atom_name(residue_atoms[iat]->name);
	    if (atom_name == named_torsions[i].base_atom_name) {
	       atom_base = residue_atoms[iat];
	    }
	    if (atom_name == named_torsions[i].atom_name_2) {
	       atom_2 = residue_atoms[iat];
	    }
	    if (atom_name == named_torsions[i].atom_name_bonded_to_H) {
	       atom_bonded_to_H = residue_atoms[iat];
	    }
	 }

	 if (atom_base && atom_2 && atom_bonded_to_H) {
	    clipper::Coord_orth pos_atom_base(atom_base->x, atom_base->y, atom_base->z);
	    clipper::Coord_orth pos_atom_2(atom_2->x, atom_2->y, atom_2->z);
	    clipper::Coord_orth pos_atom_bonded_to_H(atom_bonded_to_H->x, atom_bonded_to_H->y, atom_bonded_to_H->z);

	    clipper::Coord_orth new_pt(pos_atom_base, pos_atom_2, pos_atom_bonded_to_H,
				       1.0, // unit vector
				       clipper::Util::d2rad(named_torsions[i].angle),
				       clipper::Util::d2rad(named_torsions[i].torsion));
	    clipper::Coord_orth vect = new_pt = pos_atom_bonded_to_H;


	    // add that to the pile
	    bool found_atom = 0;
	    for (unsigned int irv=0; irv<v.size(); irv++) { 
	       if (v[irv].first == atom_bonded_to_H) {
		  v[irv].second.push_back(vect);
		  found_atom = 1;
	       }
	    }
	    if (! found_atom) {
	       std::vector<clipper::Coord_orth> cov;
	       cov.push_back(vect);
	       std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > p(atom_bonded_to_H, cov);
	       v.push_back(p);
	    } 
	 }
      }
   }

   return v;
}

// Can throw an exception
// 
// Return the position of the H-ligand atom (the atom to which the H
// is attached) and the hydrogen position - in that order.
// 
std::pair<clipper::Coord_orth, clipper::Coord_orth>
pli::flev_attached_hydrogens_t::hydrogen_pos(const pli::named_torsion_t &named_tor,
                                             mmdb::Residue *residue_p) const {
   clipper:: Coord_orth pt(0,0,0);
   clipper:: Coord_orth pt_ligand_atom(0,0,0);
   mmdb::Atom *at_1 = NULL;
   mmdb::Atom *at_2 = NULL;
   mmdb::Atom *at_3 = NULL;

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == named_tor.base_atom_name)
	 at_1 = residue_atoms[i];
      if (atom_name == named_tor.atom_name_2)
	 at_2 = residue_atoms[i];
      if (atom_name == named_tor.atom_name_bonded_to_H)
	 at_3 = residue_atoms[i];
   }

   if (! (at_1 && at_2 && at_3)) {
      throw(std::runtime_error("missing atoms in residue"));
   } else {
      clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
      clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
      clipper::Coord_orth pt_3(at_3->x, at_3->y, at_3->z);
      clipper::Coord_orth p4_h(pt_1, pt_2, pt_3,
			       named_tor.dist,
			       clipper::Util::d2rad(named_tor.angle),
			       clipper::Util::d2rad(named_tor.torsion));
      pt = p4_h;
      pt_ligand_atom = pt_3;
      
//       std::cout << "in hydrogen_pos() constructed H pos " << pt.format()
// 		<< " from atom at_1: " << at_1 << " "
// 		<< "atom at_2: " << at_2 << " "
// 		<< "atom at_3: " << at_3 << " "
// 		<< " dist: " << named_tor.dist << " " 
// 		<< " angle: " << named_tor.angle << " " 
// 		<< " torsion: " << named_tor.torsion << std::endl;
   } 
   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (pt_ligand_atom, pt);
}

// hydrogen_type is either H_IS_RIDING or H_IS_ROTATABLE
//
bool
pli::flev_attached_hydrogens_t::add_named_torsion(mmdb::Atom *h_at, mmdb::Atom *at,
                                                  const coot::dictionary_residue_restraints_t &restraints,
                                                  mmdb::Manager *mol, // 3d prodrg ligand mol
                                                  int hydrogen_type)  {

   bool found_torsion_for_this_H = 0;
   std::string atom_name_bonded_to_H(at->name);
   clipper::Coord_orth p_h(h_at->x, h_at->y, h_at->z);
   clipper::Coord_orth p_1(at->x, at->y, at->z);

   // now we work back through the restraints, finding
   // an atom that bonds to at/atom_name_bonded_to_H,
   // and one that bonds to that (by name).
   for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
      std::string atom_name_1 = restraints.bond_restraint[ibond].atom_id_1_4c();
      std::string atom_name_2 = restraints.bond_restraint[ibond].atom_id_2_4c();

      if (atom_name_bonded_to_H == atom_name_2)
	 std::swap(atom_name_1, atom_name_2);

      if (atom_name_bonded_to_H == atom_name_1) {

	 if (! restraints.is_hydrogen(atom_name_2)) { 
	    std::string At_name_2 = atom_name_2;

	    // now for the base...
	    //
	    for (unsigned int jbond=0; jbond<restraints.bond_restraint.size(); jbond++) {
	       std::string atom_name_b_1 = restraints.bond_restraint[jbond].atom_id_1_4c();
	       std::string atom_name_b_2 = restraints.bond_restraint[jbond].atom_id_2_4c();

	       if (At_name_2 == atom_name_b_2)
		  std::swap(atom_name_b_1, atom_name_b_2); // same trick

	       if (At_name_2 == atom_name_b_1) {

		  // atom_name_b_2 is the base then (maybe)

		  if (atom_name_b_2 != atom_name_bonded_to_H) {
		     if (! restraints.is_hydrogen(atom_name_b_1)) {
			std::string base_atom_name = atom_name_b_2;

			// now, where are those atoms (At_name_2 and base_atom_name)?
			mmdb::Atom *At_2 = NULL;
			mmdb::Atom *base_atom = NULL;

			int imod = 1;
			mmdb::Model *model_p = mol->GetModel(imod);
			mmdb::Chain *chain_p;
			int nchains = model_p->GetNumberOfChains();
			for (int ichain=0; ichain<nchains; ichain++) {
			   chain_p = model_p->GetChain(ichain);
			   int nres = chain_p->GetNumberOfResidues();
			   mmdb::Residue *residue_p;
			   mmdb::Atom *residue_at;
			   for (int ires=0; ires<nres; ires++) { 
			      residue_p = chain_p->GetResidue(ires);
			      int n_atoms = residue_p->GetNumberOfAtoms();
			      for (int iat=0; iat<n_atoms; iat++) {
				 residue_at = residue_p->GetAtom(iat);
				 std::string res_atom_name(residue_at->name);
				 if (res_atom_name == At_name_2)
				    At_2 = residue_at;
				 if (res_atom_name == base_atom_name)
				    base_atom = residue_at;
			      }
			   }
			}

			if (!base_atom || !At_2) {

			   if (!base_atom)
			      std::cout << "Failed to find base in 3d prodrg residue "
					<< base_atom_name << std::endl;
			   if (!At_2)
			      std::cout << "Failed to find base or At_2 in 3d prodrg residue "
					<< At_name_2 << std::endl;
			} else {
			   try { 
			      clipper::Coord_orth p_2(At_2->x, At_2->y, At_2->z);
			      clipper::Coord_orth p_base(base_atom->x, base_atom->y, base_atom->z);
					     
			      double tors_r = clipper::Coord_orth::torsion(p_base, p_2, p_1, p_h);
			      double tors = clipper::Util::rad2d(tors_r);
			      double angle = coot::angle(h_at, at, At_2);
			      double dist = clipper::Coord_orth::length(p_h, p_1);

			      pli::named_torsion_t torsion(base_atom_name,
							    At_name_2,
							    atom_name_bonded_to_H,
							    dist, angle, tors, hydrogen_type);

			      // std::cout << "  Yeah!!! adding named torsion " << std::endl;
			      named_torsions.push_back(torsion);
			      found_torsion_for_this_H = 1;
			   }
			   catch (const std::runtime_error &rte) {
			      std::cout << "WARNING:: " << rte.what() << std::endl;
			   } 
			} 
		     }
		  }
	       }
	       if (found_torsion_for_this_H)
		  break;
	    }
	 }
      }
      if (found_torsion_for_this_H)
	 break;
   }

   return found_torsion_for_this_H;
}


