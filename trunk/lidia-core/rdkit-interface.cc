/* lidia-core/rdkit-interface.cc
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
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

#ifdef MAKE_ENTERPRISE_TOOLS

#include "coot-utils.hh"
#include "rdkit-interface.hh"
#include "coot-coord-utils.hh" // after rdkit-interface.hh to avoid ::strchr problems


// This can throw an runtime_error exception (residue not in
// dictionary).
//
// Return an RDKit::RWMol with one conformer (perhaps we should do a
// conformer for each alt conf).
//
// alt_conf is an optional arg.
// 
RDKit::RWMol
coot::rdkit_mol(CResidue *residue_p, const coot::protein_geometry &geom) {

   std::string res_name = residue_p->GetResName();
   if (0)
      std::cout << "====================  here in rdkit_mol() with geometry with res_name \""
		<< res_name << "\"" << std::endl;
   
   std::pair<bool, coot::dictionary_residue_restraints_t> p = 
      geom.get_monomer_restraints_at_least_minimal(res_name);
   if (! p.first) {

      std::string m = "rdkit_mol(): residue type ";
      m += res_name;
      m += " not in dictionary";
      throw(std::runtime_error(m));

   } else {
      if (0)
	 std::cout << "......... calling rdkit_mol() with restraints that have "
		   << p.second.bond_restraint.size() << " bond restraints"
		   << std::endl;
      return rdkit_mol(residue_p, p.second);
   } 
}

RDKit::RWMol
coot::rdkit_mol(CResidue *residue_p,
		const coot::dictionary_residue_restraints_t &restraints,
		const std::string &alt_conf) {

   bool debug = false;
   
   if (debug)
      std::cout << "==================== here in rdkit_mol() with restraints that have "
		<< restraints.bond_restraint.size() << " bond restraints" << std::endl;
   
   RDKit::RWMol m;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   // this is so that we don't add multiple copies of an atom with
   // the same name (that is, only add the first atom of a given
   // atom name in a residue with alt confs).
   // 
   // We also need added_atoms so that that is what we iterate over
   // when handling bond restraints (we don't want to consider B
   // conformer for restraints if they are ignored when we are adding
   // atoms).
   // 
   std::vector<std::string> added_atom_names;
   std::vector<CAtom *>     added_atoms; // gets added to as added_atom_names gets added to.
   std::map<std::string, int> atom_index;
   int current_atom_id = 0;
   std::vector<int> bonded_atoms; // vector of the atoms that we will
				  // add to the rdkit molecule.  We
				  // don't want to add atoms that are
				  // not bonded to anything
				  // (e.g. hydrogens with mismatching names)
   
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (unsigned int iat_1=0; iat_1<n_residue_atoms; iat_1++) {
      CAtom *at_1 = residue_atoms[iat_1];
      if (! at_1->Ter) {
	 std::string atom_name_1(at_1->name);
	 std::string atom_alt_conf(at_1->altLoc);
	 if (atom_alt_conf == alt_conf) { 
	    bool found_a_bonded_atom = false;
	    for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
	       if (restraints.bond_restraint[ib].atom_id_1_4c() == atom_name_1) {
		  // is the atom to wich atom_name_1 is bonded in the
		  // atoms of the residue?
		  for (unsigned int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
		     CAtom *at_2 = residue_atoms[iat_2];
		     std::string atom_name_2 = at_2->name;
		     if (atom_name_2 == restraints.bond_restraint[ib].atom_id_2_4c()) {
			found_a_bonded_atom = true;
			break;
		     } 
		  }
	       } 
	       if (restraints.bond_restraint[ib].atom_id_2_4c() == atom_name_1) {
		  // is the atom to wich atom_name_1 is bonded in the
		  // atoms of the residue?
		  for (unsigned int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
		     CAtom *at_2 = residue_atoms[iat_2];
		     std::string atom_name_2 = at_2->name;
		     if (atom_name_2 == restraints.bond_restraint[ib].atom_id_1_4c()) {
			found_a_bonded_atom = true;
			break;
		     } 
		  }
	       }
	    }

	    if (found_a_bonded_atom)
	       bonded_atoms.push_back(iat_1);
	 }
      }
   }

   if (debug)
      std::cout << "DEBUG:: number of bonded atoms with alt conf \"" << alt_conf << "\" found: "
		<< bonded_atoms.size() << std::endl;


   for (unsigned int iat=0; iat<bonded_atoms.size(); iat++) {

      CAtom *at = residue_atoms[bonded_atoms[iat]];
      std::string atom_name(at->name);
      if (debug)
	 std::cout << "   handling atom " << iat << " of " << n_residue_atoms << " " 
		   << atom_name << std::endl;
	 
      // only add the atom if the atom_name is not in the list of
      // already-added atom names.
      if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name) != added_atom_names.end()) {
	 std::cout << "!!!! Problem? atom name \"" << atom_name
		   << "\" was already added" << std::endl;
	    
      } else { 
	 RDKit::Atom *rdkit_at = new RDKit::Atom;
	 try {
	    std::string ele_capped =
	       coot::util::capitalise(coot::util::remove_leading_spaces(at->element));
	    int atomic_number = tbl->getAtomicNumber(ele_capped);
	    rdkit_at->setAtomicNum(atomic_number);
	    rdkit_at->setMass(tbl->getAtomicWeight(atomic_number));
	    rdkit_at->setProp("name", atom_name);


	    // set the valence from they type energy.  Abstract?
	    //
	    std::string type_energy = restraints.type_energy(at->name);
	    if (type_energy != "") {
	       if (type_energy == "NT") {
		  rdkit_at->setFormalCharge(1);
	       }
		  
	       // other NT*s will drop hydrogens in RDKit, no need to
	       // fix up formal charge (unless there is a hydrogen! Hmm).
	       
	       // now that phosphates are undelocalized, we don't
	       // need to make this hack:
	       // 
	       // if (type_energy == "P") {
	       //    rdkit_at->setFormalCharge(1);
	       // }
		  
	    }

	    // set the chirality
	    // (if this atom is chiral)
	    //
	    for (unsigned int ichi=0; ichi<restraints.chiral_restraint.size(); ichi++) { 
	       if (restraints.chiral_restraint[ichi].atom_id_c_4c() == atom_name) {
		  if (!restraints.chiral_restraint[ichi].has_unassigned_chiral_volume()) {
		     if (!restraints.chiral_restraint[ichi].is_a_both_restraint()) {
			// e.g. RDKit::Atom::CHI_TETRAHEDRAL_CCW;
			RDKit::Atom::ChiralType chiral_tag =
			   get_chiral_tag(residue_p, restraints, residue_atoms[iat]);
			rdkit_at->setChiralTag(chiral_tag);
		     } 
		  } 
	       } 
	    }
	    
	    m.addAtom(rdkit_at);
	       
	    if (debug) 
	       std::cout << "      about to add atom with name \"" << atom_name
			 << "\" to added_atom_names which is currently of size "
			 << added_atom_names.size() << std::endl;
	       
	    added_atom_names.push_back(atom_name);
	    added_atoms.push_back(residue_atoms[iat]);
	    atom_index[atom_name] = current_atom_id;
	    current_atom_id++; // for next round
	 }
	 catch (std::exception rte) {
	    std::cout << rte.what() << std::endl;
	 }
      }
   }

   if (debug) {
      std::cout << "DEBUG::number of atoms in rdkit mol: " << m.getNumAtoms() << std::endl;
      std::cout << "DEBG:: number of bond restraints: " << restraints.bond_restraint.size()
		<< std::endl;
   } 
   
   for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
      if (debug)
	 std::cout << "   handling bond " << ib << " of " << restraints.bond_restraint.size()
		   << " :" << restraints.bond_restraint[ib].atom_id_1_4c() << ": " 
		   << " :" << restraints.bond_restraint[ib].atom_id_2_4c() << ": " 
		   << std::endl;
      RDKit::Bond::BondType type = convert_bond_type(restraints.bond_restraint[ib].type());
      RDKit::Bond *bond = new RDKit::Bond(type);
	    
      std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
      std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
      std::string ele_1 = restraints.element(atom_name_1);
      std::string ele_2 = restraints.element(atom_name_2);
      int idx_1 = -1; // unset initially
      int idx_2 = -1; // unset

      // this block sets idx_1 and idx_2
      //
      for (unsigned int iat=0; iat<added_atom_names.size(); iat++) {
	 if (added_atom_names[iat] == atom_name_1) { 
	    idx_1 = iat;
	    break;
	 }
      }
      for (unsigned int iat=0; iat<added_atom_names.size(); iat++) { 
	 if (added_atom_names[iat] == atom_name_2) { 
	    idx_2 = iat;
	    break;
	 }
      }

      if (idx_1 != -1) { 
	 if (idx_2 != -1) {	 
	    bond->setBeginAtomIdx(idx_1);
	    bond->setEndAtomIdx(  idx_2);
	    if (type == RDKit::Bond::AROMATIC) { 
	       bond->setIsAromatic(true);
	       m[idx_1]->setIsAromatic(true);
	       m[idx_2]->setIsAromatic(true);
	    }
	    m.addBond(bond); // worry about ownership or memory leak.
	    
	 } else {
	    if (ele_2 != " H") { 
	       std::cout << "WARNING:: oops, bonding in rdkit_mol() "
			 << "failed to get atom index idx_2 for atom name: "
			 << atom_name_2 << " ele :" << ele_2 << ":" << std::endl;
	       // give up trying to construct this thing then.
	       std::string message = "Failed to get atom index for atom name \"";
	       message += atom_name_2;
	       message += "\"";
	       throw std::runtime_error(message);
	    }
	 }
      } else {
	 if (ele_1 != " H") { 
	    std::cout << "WARNING:: oops, bonding in rdkit_mol() "
		      << "failed to get atom index idx_1 for atom name: \""
		      << atom_name_1 << "\" ele :" << ele_1 << ":" << std::endl;
	    // give up trying to construct this thing then.
	    std::string message = "Failed to get atom index for atom name \"";
	    message += atom_name_1;
	    message += "\"";
	    throw std::runtime_error(message);
	 }
      }
   }

   // Now all the bonds are in place.  We can now try to add an extra H to the ring N that
   // need them.  We do this after, because all the molecule bonds have to be in place for
   // us to know if the H is there already (tested in add_H_to_ring_N_as_needed()).
   
   std::vector<int> Hs_added_list; // a list of atoms to which extra Hs have been added
				      // (this is needed, because we only want to add a
				      // one H to an aromatic N and as we run though the
				      // bond list, we typically find N-C and C-N bonds,
				      // which would result in the H being added twice -
				      // not good.

   for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) { 

      RDKit::Bond::BondType type = convert_bond_type(restraints.bond_restraint[ib].type());
      if (type == RDKit::Bond::AROMATIC) { 
	 std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
	 std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
	 std::string ele_1 = restraints.element(atom_name_1);
	 std::string ele_2 = restraints.element(atom_name_2);
	 int idx_1 = -1; // unset
	 int idx_2 = -1; // unset

	 // we can't run through n_residue_atoms because the atom in the
	 // CResidue may not have been added to the atom in the rdkit
	 // molecule (as is the case for an alt conf).

	 for (unsigned int iat=0; iat<m.getNumAtoms(); iat++) {
	    try {
	       std::string name;
	       RDKit::ATOM_SPTR at_p = m[iat];
	       at_p->getProp("name", name);
	       if (name == atom_name_1)
		  idx_1 = iat;
	       if (name == atom_name_2)
		  idx_2 = iat;
	    }
	    catch (KeyErrorException &err) {
	       // this happens for alt conf
	       // std::cout << "caught no-name exception in rdkit_mol H-block" << std::endl;
	    }
	 }

	 if (debug) { 
	    std::cout << "idx_1 " << idx_1 << " " << n_residue_atoms << std::endl;
	    std::cout << "idx_2 " << idx_2 << " " << n_residue_atoms << std::endl;
	 }
      
	 if (idx_1 != -1) { 
	    if (idx_2 != -1) {	 
   
	       // special edge case for aromatic ring N that may need an H attached for
	       // kekulization (depending on energy type).

	       if (m[idx_1]->getAtomicNum() == 7) {
		  if (std::find(Hs_added_list.begin(), Hs_added_list.end(), idx_1) == Hs_added_list.end()) { 
		     std::string n = add_H_to_ring_N_as_needed(&m, idx_1, atom_name_1, restraints);
		     if (0)
			std::cout << "testing 1 idx_1 " << idx_1 << " idx_2 " << idx_2
				  << " n was :" << n << ":" << std::endl;
		     if (n != "")
			added_atom_names.push_back(n);
		     Hs_added_list.push_back(idx_1);
		  }
	       }
	       if (m[idx_2]->getAtomicNum() == 7) {
		  if (std::find(Hs_added_list.begin(), Hs_added_list.end(), idx_2) == Hs_added_list.end()) { 
		     std::string n = add_H_to_ring_N_as_needed(&m, idx_2, atom_name_2, restraints);
		     // std::cout << "testing 2 idx_1 " << idx_1 << " idx_2 " << idx_2
		     // << " n was :" << n << ":" << std::endl;
		     if (n != "")
			added_atom_names.push_back(n);
		     Hs_added_list.push_back(idx_2);
		  }
	       }
	    }
	 }
      }
   }


    if (debug)
       std::cout << "=============== calling undelocalise() " << &m << std::endl;
    coot::undelocalise(&m);

    if (debug)
       std::cout << "---------------------- calling assign_formal_charges() -----------" << std::endl;
    coot::assign_formal_charges(&m);
   
    if (debug)
       std::cout << "---------------------- calling cleanUp() -----------" << std::endl;
    RDKit::MolOps::cleanUp(m);

   
   // OK, so cleanUp() doesn't fix the N charge problem our prodrg molecule
   // 
   if (debug) { // debug, formal charges
      std::cout << "::::::::::::::::::::::::::: after cleanup :::::::::::::::::"
		<< std::endl;
      int n_mol_atoms = m.getNumAtoms();
      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 RDKit::ATOM_SPTR at_p = m[iat];
	 std::string name = "";
	 try {
	    at_p->getProp("name", name);
	 }
	 catch (KeyErrorException kee) {
	    std::cout << "caught no-name for atom exception in make_molfile_molecule(): "
		      <<  kee.what() << std::endl;
	 }
	 int formal_charge = at_p->getFormalCharge();
	 std::cout << " " << iat << " " << name << " formal charge: "
		   << formal_charge << std::endl;
      }
   }

   if (debug)
      std::cout << "DEBUG:: sanitizeMol() " << std::endl;
   RDKit::MolOps::sanitizeMol(m);

   if (debug)
      std::cout << "in constructing rdk molecule now adding a conf" << std::endl;
   RDKit::Conformer *conf = new RDKit::Conformer(added_atom_names.size());
   conf->set3D(true);
      
   // Add positions to the conformer (only the first instance of an
   // atom with a particular atom name).
   // 
   for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::string atom_alt_conf(residue_atoms[iat]->altLoc);
      if (atom_alt_conf == alt_conf) { 
	 std::map<std::string, int>::const_iterator it = atom_index.find(atom_name);
	 if (it != atom_index.end()) {
	    RDGeom::Point3D pos(residue_atoms[iat]->x,
				residue_atoms[iat]->y,
				residue_atoms[iat]->z);
	    conf->setAtomPos(it->second, pos);
	    if (debug) 
	       std::cout << "in construction of rdkit mol: making a conformer atom "
			 << iat << " " << it->second << " " << atom_name << " at pos "
			 << pos << std::endl;
	 }
      }
   }
   m.addConformer(conf);
   if (debug) 
      std::cout << "ending construction of rdkit mol: n_atoms " << m.getNumAtoms()
		<< std::endl;
   return m;
}

// can throw a std::runtime_error or std::exception.
//
// should kekulize flag be an argmuent?
// 
RDKit::RWMol
coot::rdkit_mol_sanitized(CResidue *residue_p, const protein_geometry &geom) {

   RDKit::RWMol mol = coot::rdkit_mol(residue_p, geom);
   rdkit_mol_sanitize(mol);
   return mol;
}

void
coot::rdkit_mol_sanitize(RDKit::RWMol &mol) {
   
   // clear out any cached properties
   mol.clearComputedProps();
   // clean up things like nitro groups
   RDKit::MolOps::cleanUp(mol);
   // update computed properties on atoms and bonds:
   mol.updatePropertyCache();
   RDKit::MolOps::Kekulize(mol);
   RDKit::MolOps::assignRadicals(mol);
	    
   // then do aromaticity perception
   RDKit::MolOps::setAromaticity(mol);
    
   // set conjugation
   RDKit::MolOps::setConjugation(mol);
	       
   // set hybridization
   RDKit::MolOps::setHybridization(mol); // non-linear ester bonds.

   // remove bogus chirality specs:
   RDKit::MolOps::cleanupChirality(mol);

}



RDKit::Bond::BondType
coot::convert_bond_type(const std::string &t) {
   RDKit::Bond::BondType bt = RDKit::Bond::UNSPECIFIED;
   if (t == "single")
      bt = RDKit::Bond::SINGLE;
   if (t == "double")
      bt = RDKit::Bond::DOUBLE;
   if (t == "triple")
      bt = RDKit::Bond::TRIPLE;
   if (t == "coval")
      bt = RDKit::Bond::SINGLE;
   if (t == "deloc")
      bt = RDKit::Bond::ONEANDAHALF;
   if (t == "aromatic")
      bt = RDKit::Bond::AROMATIC;
   
   if (0) { // debug
      std::cout << "created RDKit bond type " << bt;
      if (bt == RDKit::Bond::AROMATIC)
	 std::cout << " (aromatic)";
      std::cout << std::endl;
   }
   
   return bt;
}

// used in the rdkit_mol() "constructor".
// 
RDKit::Atom::ChiralType
coot::get_chiral_tag(CResidue *residue_p,
		     const dictionary_residue_restraints_t &restraints,
		     CAtom *atom_p) {

   RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_UNSPECIFIED; // as yet
   
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   std::string atom_name = atom_p->name;
   
   bool atom_orders_match = 0; 
   // does the order of the restraints match the order of the atoms?
   //
   for (unsigned int ichi=0; ichi<restraints.chiral_restraint.size(); ichi++) { 
      if (restraints.chiral_restraint[ichi].atom_id_c_4c() == atom_name) {
	 const coot::dict_chiral_restraint_t &chiral_restraint = restraints.chiral_restraint[ichi];

	 int n_neigbours_found = 0;
	 unsigned int i_next = 0;
	 for (unsigned int iat=0; iat<n_residue_atoms; iat++) { 
	    std::string atom_name_local = residue_atoms[iat]->name;
	    if (atom_name_local == chiral_restraint.atom_id_1_4c()) {
	       n_neigbours_found++;
	       i_next = iat+1;
	       break;
	    }
	 }
	 for (unsigned int iat=i_next; iat<n_residue_atoms; iat++) {
	    std::string atom_name_local = residue_atoms[iat]->name;
	    if (atom_name_local == chiral_restraint.atom_id_2_4c()) {
	       n_neigbours_found++;
	       i_next = iat+1;
	       break;
	    }
	 }
	 for (unsigned int iat=i_next; iat<n_residue_atoms; iat++) {
	    std::string atom_name_local = residue_atoms[iat]->name;
	    if (atom_name_local == chiral_restraint.atom_id_3_4c()) {
	       n_neigbours_found++;
	       break;
	    }
	 }

	 if (n_neigbours_found == 3) {
	    // yes, they match
	    atom_orders_match = 1;
	 }

	 // This bit needs checking
	 // 
	 if (atom_orders_match) {
	    if (chiral_restraint.volume_sign == 1)
	       chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
	    else 
	       chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
	 } else {
	    if (chiral_restraint.volume_sign == 1)
	       chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
	    else 
	       chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
	 } 
	 break;
      } 
   }

   return chiral_tag;
} 



// tweaking function used by rdkit mol construction function.
// (change mol maybe).
//
// Don't try to add an H if there is already an H on this atom (typically, this ligand had
// hydrogens (everywhere) already).
//
// Try to find the name of the Hydrogen from the bond restraints.
// If not found, add an atom called "-".
// 
// @return the added hydrogen name - or "" if nothing was added.
// 
std::string
coot::add_H_to_ring_N_as_needed(RDKit::RWMol *mol,
				int idx, const std::string &atom_name,
				const coot::dictionary_residue_restraints_t &restraints) {


   std::string r = "";
   std::string energy_type = restraints.type_energy(atom_name);
   
   if ((energy_type == "NR15") || (energy_type == "NR16")) {

      // first, is there an H bonded to this atom already?
      bool already_there = 0;
      int n_bonds = mol->getNumBonds();
      for (unsigned int ib=0; ib<n_bonds; ib++) {
	 const RDKit::Bond *bond_p = mol->getBondWithIdx(ib);
	 int idx_1 = bond_p->getBeginAtomIdx();
	 int idx_2 = bond_p->getEndAtomIdx();
	 if (idx_1 == idx)
	    if ((*mol)[idx_2]->getAtomicNum() == 1) {
	       already_there = 1;
	       break;
	    }
	 if (idx_2 == idx)
	    if ((*mol)[idx_1]->getAtomicNum() == 1) {
	       already_there = 1;
	       break;
	    }
      }

      if (! already_there) { 
      
	 // -------------  add an H atom --------------------------
	 // 
	 RDKit::Atom *at = new RDKit::Atom;
	 at->setAtomicNum(1);
	 int idx_for_H = mol->addAtom(at);

	 // std::string name_H = "-";
	 // what is the Name of this H?

	 std::string name_H = "-";
	 for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) { 
	    if (restraints.bond_restraint[ib].atom_id_1_4c() == atom_name) {
	       if (restraints.element(restraints.bond_restraint[ib].atom_id_2_4c()) == " H") { 
		  name_H = restraints.bond_restraint[ib].atom_id_2_4c();
		  break;
	       }
	    } 
	    if (restraints.bond_restraint[ib].atom_id_2_4c() == atom_name) {
	       if (restraints.element(restraints.bond_restraint[ib].atom_id_1_4c()) == " H") { 
		  name_H = restraints.bond_restraint[ib].atom_id_1_4c();
		  break;
	       }
	    } 
	 }
      
      
	 if (name_H != "-") {
	    at->setProp("name", name_H);
	 }
	 r = name_H;

	 // -------------  now add a bond --------------------------

	 RDKit::Bond *bond = new RDKit::Bond(RDKit::Bond::SINGLE);
	 bond->setBeginAtomIdx(idx);
	 bond->setEndAtomIdx(idx_for_H);
	 mol->addBond(bond);
      }
   }
   return r;
}


lig_build::molfile_molecule_t
coot::make_molfile_molecule(const RDKit::ROMol &rdkm, int iconf) {

   lig_build::molfile_molecule_t mol;
   int n_conf  = rdkm.getNumConformers();

   if (n_conf) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
      
      RDKit::Conformer conf = rdkm.getConformer(iconf);
      int n_mol_atoms = rdkm.getNumAtoms();

      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 RDKit::ATOM_SPTR at_p = rdkm[iat];
	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
	 std::string name = "";
	 try {
	    at_p->getProp("name", name);
	 }
	 catch (KeyErrorException kee) {
	    std::cout << "caught no-name for atom exception in make_molfile_molecule(): "
		      <<  kee.what() << std::endl;
	 } 
	 clipper::Coord_orth pos(r_pos.x, r_pos.y, r_pos.z);
	 int n = at_p->getAtomicNum();
	 std::string element = tbl->getElementSymbol(n);
	 int charge = at_p->getFormalCharge();
	 lig_build::molfile_atom_t mol_atom(pos, element, name);
	 mol.add_atom(mol_atom);
      }

      int n_bonds = rdkm.getNumBonds();
      for (unsigned int ib=0; ib<n_bonds; ib++) {
	 const RDKit::Bond *bond_p = rdkm.getBondWithIdx(ib);
	 int idx_1 = bond_p->getBeginAtomIdx();
	 int idx_2 = bond_p->getEndAtomIdx();
	 lig_build::bond_t::bond_type_t bt = convert_bond_type(bond_p->getBondType());
	 if (0) 
	    std::cout << "   make_molfile_molecule() " << idx_1 << " "
		      << idx_2 << "      "
		      << bond_p->getBondType() << " to "
		      << bt 
		      << std::endl;
	 lig_build::molfile_bond_t mol_bond(idx_1, idx_2, bt);
	 RDKit::Bond::BondDir bond_dir = bond_p->getBondDir();
	 if (bond_dir != RDKit::Bond::NONE) {
	    if (bond_dir == RDKit::Bond::BEGINWEDGE)
	       mol_bond.bond_type = lig_build::bond_t::OUT_BOND;
	    if (bond_dir == RDKit::Bond::BEGINDASH)
	       mol_bond.bond_type = lig_build::bond_t::IN_BOND;
	 }
	 mol.add_bond(mol_bond);
      }
   }
   
   return mol;
}

// returns NULL on fail. Caller deletes.
//
CResidue *
coot::make_residue(const RDKit::ROMol &rdkm, int iconf, const std::string &res_name) {

   CResidue *residue_p = NULL;
   lig_build::molfile_molecule_t mol = coot::make_molfile_molecule(rdkm, iconf);

   // now convert mol to a CResidue *
   // 
   if (mol.atoms.size()) {
      residue_p = new CResidue;
      residue_p->seqNum = 1;
      residue_p->SetResName(res_name.c_str());
      CChain *chain_p = new CChain;
      chain_p->SetChainID("");
      chain_p->AddResidue(residue_p);
      for (unsigned int iat=0; iat<mol.atoms.size(); iat++) { 
	 CAtom *at = new CAtom;
	 std::string atom_name = mol.atoms[iat].name; // overridden hopefully
	 at->SetAtomName(atom_name.c_str());
	 at->SetElementName(mol.atoms[iat].element.c_str());
	 at->SetCoordinates(mol.atoms[iat].atom_position.x(),
			    mol.atoms[iat].atom_position.y(),
			    mol.atoms[iat].atom_position.z(),
			    1.0, 30.0);
	 at->Het = 1;
	 residue_p->AddAtom(at);
      }
   }

   return residue_p;
}


lig_build::bond_t::bond_type_t
coot::convert_bond_type(const RDKit::Bond::BondType &type) {

   lig_build::bond_t::bond_type_t t = lig_build::bond_t::SINGLE_BOND; // Hmmm..
   
   if (type == RDKit::Bond::SINGLE)
      t = lig_build::bond_t::SINGLE_BOND;
   if (type == RDKit::Bond::DOUBLE)
      t = lig_build::bond_t::DOUBLE_BOND;
   if (type == RDKit::Bond::TRIPLE)
      t = lig_build::bond_t::TRIPLE_BOND;
   if (type == RDKit::Bond::AROMATIC)
      t = lig_build::bond_t::AROMATIC_BOND;
    
    // PUTS_IT_BACK?
    // enabling this makes phosphate bonds disconnected in the
    // thumbnail and Residue->2D
    // 
    //     if (type == RDKit::Bond::ONEANDAHALF)
    //        t = lig_build::bond_t::DELOC_ONE_AND_A_HALF;

   return t;
}

// fiddle with rdkm.  This can throw a std::exception
// 
void
coot::remove_non_polar_Hs(RDKit::RWMol *rdkm) {

   int n_bonds = rdkm->getNumBonds();
   // std::vector<RDKit::ATOM_SPTR> atoms_to_be_deleted;
   std::vector<RDKit::Atom *> atoms_to_be_deleted;
   for (unsigned int ib=0; ib<n_bonds; ib++) {
      RDKit::Bond *bond_p = rdkm->getBondWithIdx(ib);
      int idx_1 = bond_p->getBeginAtomIdx();
      int idx_2 = bond_p->getEndAtomIdx();
      RDKit::ATOM_SPTR at_p_1 = (*rdkm)[idx_1];
      RDKit::ATOM_SPTR at_p_2 = (*rdkm)[idx_2];
      // If this was a bond for a hydrogen attached to a carbon, delete it.
      if ((at_p_1->getAtomicNum() == 1) && (at_p_2->getAtomicNum() == 6)) {
	 // rdkm->removeBond(idx_1, idx_2);
	 atoms_to_be_deleted.push_back(at_p_1.get());
      }
      if ((at_p_2->getAtomicNum() == 1) && (at_p_1->getAtomicNum() == 6)) { 
	 // rdkm->removeBond(idx_1, idx_2);
	 atoms_to_be_deleted.push_back(at_p_2.get());
      }
      if ((at_p_1->getAtomicNum() == 1) && (at_p_2->getAtomicNum() == 5)) {
	 // rdkm->removeBond(idx_1, idx_2);
	 atoms_to_be_deleted.push_back(at_p_1.get());
      }
      if ((at_p_2->getAtomicNum() == 1) && (at_p_1->getAtomicNum() == 5)) { 
	 // rdkm->removeBond(idx_1, idx_2);
	 atoms_to_be_deleted.push_back(at_p_2.get());
      }
   }
   for (unsigned int i=0; i<atoms_to_be_deleted.size(); i++) { 
      rdkm->removeAtom(atoms_to_be_deleted[i]);
   }

   std::cout << "DEBUG:: remove_non_polar_Hs() clearComputedProps() " << std::endl;
   rdkm->clearComputedProps();
   // clean up things like nitro groups
   std::cout << "DEBUG:: remove_non_polar_Hs() cleanUp() " << std::endl;
   RDKit::MolOps::cleanUp(*rdkm);
   // update computed properties on atoms and bonds:

   coot::assign_formal_charges(rdkm);

   std::cout << "DEBUG:: remove_non_polar_Hs() updatePropertyCache() " << std::endl;
   rdkm->updatePropertyCache();
   
   std::cout << "DEBUG:: remove_non_polar_Hs() kekulize() " << std::endl;
   RDKit::MolOps::Kekulize(*rdkm);
   std::cout << "DEBUG:: remove_non_polar_Hs() assignRadicals() " << std::endl;
   RDKit::MolOps::assignRadicals(*rdkm);
   // set conjugation
   std::cout << "DEBUG:: remove_non_polar_Hs() setConjugation() " << std::endl;
   RDKit::MolOps::setConjugation(*rdkm);
   // set hybridization
   std::cout << "DEBUG:: remove_non_polar_Hs() setHybridization() " << std::endl;
   RDKit::MolOps::setHybridization(*rdkm); 
   // remove bogus chirality specs:
   std::cout << "DEBUG:: remove_non_polar_Hs() cleanupChirality() " << std::endl;
   RDKit::MolOps::cleanupChirality(*rdkm);
   
}


// Delete a hydrogen (if possible) from a N with valence 4.
void
coot::delete_excessive_hydrogens(RDKit::RWMol *rdkm) {

   int n_mol_atoms = rdkm->getNumAtoms();   
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
      if (at_p->getAtomicNum() == 7) {
	 
	 int e_valence = at_p->getExplicitValence();

	 // std::cout << " atom N has explicit valence: " << e_valence << std::endl;

	 if (e_valence == 4) { 

	    RDKit::ROMol::OEDGE_ITER current, end;
	    boost::tie(current, end) = rdkm->getAtomBonds(at_p.get());
	    RDKit::Atom *last_hydrogen_p = NULL;
	    while (current != end) {
	       RDKit::BOND_SPTR bond=(*rdkm)[*current];
	       // is this a bond to a hydrogen?
	       int idx = bond->getOtherAtomIdx(iat);
	       RDKit::ATOM_SPTR at_other_p = (*rdkm)[iat];
	       if (at_other_p->getAtomicNum() == 1)
		  last_hydrogen_p = at_other_p.get();
	       current++;
	    }

	    if (last_hydrogen_p) {
	       // delete it then
	       rdkm->removeAtom(last_hydrogen_p);
	    }
	 }
      }
   }
}

// Put a +1 charge on Ns with 4 bonds, 
void
coot::assign_formal_charges(RDKit::RWMol *rdkm) {

   
   int n_mol_atoms = rdkm->getNumAtoms();
   if (0)
      std::cout << "---------------------- in assign_formal_charges() with " << n_mol_atoms
		<< " atoms -----------" << std::endl;

   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
      // debug
      if (0)
	 std::cout << "in assign_formal_charges() calcExplicitValence on atom "
		   << iat << "/" << n_mol_atoms
		   << "  " << at_p->getAtomicNum() << std::endl;
      at_p->calcExplicitValence(false);
   }
   
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
      if (0) 
	 std::cout << "atom " << iat << "/" << n_mol_atoms << "  " << at_p->getAtomicNum()
		   << " with valence " << at_p->getExplicitValence()
		   << std::endl;
      if (at_p->getAtomicNum() == 7) { // N
	 int e_valence = at_p->getExplicitValence();
	 if (0) 
	    std::cout << " atom N has explicit valence: " << e_valence << std::endl;
	 if (e_valence == 4) {
	    if (0)
	       std::cout << ".......... assign_formal_charges: found one! "
			 << at_p << std::endl;
	    at_p->setFormalCharge(1);
	 }
      }
      if (at_p->getAtomicNum() == 12) { // Mg
	 at_p->setFormalCharge(2);
      }
   }
   if (0) 
      std::cout << "----------- normal completion of assign_formal_charges()" << std::endl;
}

// a wrapper for the above, matching hydrogens names to the
// dictionary.  Add atoms to residue_p, return success status.
// 
std::pair<bool, std::string>
coot::add_hydrogens_with_rdkit(CResidue *residue_p,
			      const coot::dictionary_residue_restraints_t &restraints) {

   bool r = 0;
   std::string error_message;
   if (residue_p) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
      try {

	 // first save the existing hydrogen names. We don't want to
	 // add a hydrogen with the same name as an atom we already
	 // have.
	 //
	 std::vector<std::string> existing_H_names;
	 PPCAtom residue_atoms = 0;
	 int n_residue_atoms;
	 residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	    if (! residue_atoms[iat]->isTer()) {
	       std::string ele = residue_atoms[iat]->element;
	       if (ele == " H")
		  existing_H_names.push_back(residue_atoms[iat]->name);
	    }
	 }

	 // get out now if thre are existing hydrogens,
	 // MolOps::addHs() throws an exception (not clear to me why)
	 // if there are hydrogens on the molecule already.
	 if (existing_H_names.size()) {
	    return std::pair<bool, std::string> (0, "Ligand contains (some) hydrogens already");
	 }


	 std::vector<std::string> residue_alt_confs = util::get_residue_alt_confs(residue_p);

	 for (unsigned int i_alt_conf=0; i_alt_conf<residue_alt_confs.size(); i_alt_conf++) { 

	    std::string alt_conf = residue_alt_confs[i_alt_conf];
	    RDKit::RWMol m_no_Hs = rdkit_mol(residue_p, restraints, residue_alt_confs[i_alt_conf]);
	    int n_mol_atoms = m_no_Hs.getNumAtoms();
	    // std::cout << ".......  pre H addition mol has " << n_mol_atoms
	    // << " atoms" << std::endl;

	    coot::undelocalise(&m_no_Hs);
	 
	    for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	       RDKit::ATOM_SPTR at_p = m_no_Hs[iat];
	       at_p->calcImplicitValence(true);
	    }

	    // std::cout << "MolOps:: adding hydrogens " << std::endl;
	    RDKit::ROMol *m_pre = RDKit::MolOps::addHs(m_no_Hs, 0, 1);
	    RDKit::RWMol m(*m_pre);
	    // std::cout << "....... post H addition mol has " << m_pre->getNumAtoms()
	    // << std::endl;
	    int n_atoms_new = m_pre->getNumAtoms();
	    int n_conf = m_pre->getNumConformers();

	    double vdwThresh=10.0;
	    int confId = 0;
	    bool ignoreInterfragInteractions=true;
	    int maxIters = 500;

	    // std::cout << "==== constructForceField()......" << std::endl;
	 
	    ForceFields::ForceField *ff =
	       RDKit::UFF::constructForceField(m, vdwThresh, confId,
					       ignoreInterfragInteractions);

	    for (unsigned int iat=0; iat<n_mol_atoms; iat++)
	       ff->fixedPoints().push_back(iat);

	    // std::cout << "==== ff->initialize() ......" << std::endl;
	    ff->initialize();
	    // std::cout << "==== ff->minimize() ......" << std::endl;
	    int res=ff->minimize(maxIters);
	    // std::cout << "rdkit minimize() returns " << res << std::endl;
	    delete ff;
	 

	    if (! n_conf) {
	       std::cout << "ERROR:: mol with Hs: no conformers" << std::endl;
	    } else { 
	       RDKit::Conformer conf = m_pre->getConformer(0);
	       std::vector<std::string> H_names_already_added;
	 
	       for (unsigned int iat=0; iat<n_atoms_new; iat++) {
		  RDKit::ATOM_SPTR at_p = (*m_pre)[iat];
		  RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
		  std::string name = "";
		  try {
		     at_p->getProp("name", name);
		     CAtom *res_atom = residue_p->GetAtom(name.c_str());
		     if (res_atom) {
			std::cout << "setting heavy atom " << name << " to "
				  << r_pos << std::endl;
			res_atom->x = r_pos.x;
			res_atom->y = r_pos.y;
			res_atom->z = r_pos.z;
		     }
		     
		  }
		  catch (KeyErrorException kee) {

		     // OK...
		     //
		     // typically when we get here, that's because the
		     // atom is a new one, generated by RDKit.
		  
		     std::string name = coot::infer_H_name(iat, at_p, &m, restraints,
							   H_names_already_added);
		     if (name != "") {

			// add atom if the name is not already there:
			// 
			if (std::find(existing_H_names.begin(),
				      existing_H_names.end(),
				      name) == existing_H_names.end()) {
			
			   H_names_already_added.push_back(name);
			
			   int n = at_p->getAtomicNum();
			   std::string element = tbl->getElementSymbol(n);
			
			   CAtom *at = new CAtom;
			   at->SetAtomName(name.c_str());
			   // at->SetElementName(element.c_str()); // FIXME?
			   at->SetElementName(" H");  // PDBv3 FIXME
			   at->SetCoordinates(r_pos.x, r_pos.y, r_pos.z, 1.0, 30.0);
			   at->Het = 1;
			   if (alt_conf != "") {
			      strncpy(at->altLoc, alt_conf.c_str(), alt_conf.length()+1);
			   }
			   residue_p->AddAtom(at);
			   r = 1;
			}
		     }
		  }
	       }
	    }
	 }
	 
	 // delete m;
      }
      catch (std::runtime_error e) {
	 std::cout << e.what() << std::endl;
      }
      catch (std::exception rdkit_error) {
	 std::cout << rdkit_error.what() << std::endl;
      }
   }
   return std::pair<bool, std::string> (r, error_message);
}

// atom_p is a hydrogen atom we presume, of degree 1.  This is tested
// before the restraints are checked.
// 
std::string
coot::infer_H_name(int iat,
		   RDKit::ATOM_SPTR atom_p,
		   const RDKit::ROMol *mol,
		   const dictionary_residue_restraints_t &restraints,
		   const std::vector<std::string> &H_names_already_added) {

   std::string r = "";

   unsigned int deg = mol->getAtomDegree(atom_p.get());
   if (deg == 1) {
      RDKit::ROMol::OEDGE_ITER current, end;
      boost::tie(current, end) = mol->getAtomBonds(atom_p.get());
      while (current != end) {
	 RDKit::BOND_SPTR bond=(*mol)[*current];
	 int idx = bond->getOtherAtomIdx(iat);
	 RDKit::ATOM_SPTR other_atom_p = (*mol)[idx];
	 std::string bonding_atom_name;
	 try {
	    other_atom_p->getProp("name", bonding_atom_name);
	    // in the restraints, what is the name of the hydrogen
	    // bonded to atom with name bonding_atom_name? (if any).
	    std::vector<std::string> nv = 
	       restraints.get_attached_H_names(bonding_atom_name);
	    for (unsigned int i=0; i<nv.size(); i++) { 
	       if (std::find(H_names_already_added.begin(),
			     H_names_already_added.end(), nv[i]) ==
		   H_names_already_added.end()) {
		  r = nv[i];
		  break;
	       }
	    }
	 }
	 catch (KeyErrorException kee) {
	    // this should not happen, there should be no way we get
	    // here where we have a hydrogen (check in calling
	    // function) with no name attached to another atom with no
	    // name.
	    std::cout << "ERROR:: in infer_H_name() bonding atom with no name "
		      << std::endl;
	 }
	 current++;
      }
   }
   // std::cout << "returning infered H name :" << r << ":" << std::endl;
   return r;
}


// tweak rdkmol
int
coot::add_2d_conformer(RDKit::ROMol *rdk_mol, double weight_for_3d_distances) {

   bool debug = false;

   int icurrent_conf = 0; // the conformer number from which the
                          // distance matrix is generated.  Should this
			  // be passed?

   int n_conf  = rdk_mol->getNumConformers();
   if (n_conf == 0) {
      std::cout << "WARNING:: no conformers in add_2d_conformer() - aborting"
		<< std::endl;
   }

   int n_mol_atoms = rdk_mol->getNumAtoms();   

   if (debug) 
      std::cout << "::::: add_2d_conformer before compute2DCoords n_atoms: "
		<< rdk_mol->getConformer(0).getNumAtoms()
		<< " n_bonds " << rdk_mol->getNumBonds() << std::endl;

   // We must call calcImplicitValence() before getNumImplictHs()
   // [that is to say that compute2DCoords() calls getNumImplictHs()
   // without calling calcImplicitValence() and that results in problems].
   //
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR at_p = (*rdk_mol)[iat];
      at_p->calcImplicitValence(true);
   }

   int n_items = n_mol_atoms * (n_mol_atoms - 1)/2;

   double *cData = new double[n_items]; // handled by smart pointer, I think.
   for (unsigned int i=0; i<n_items; i++)
      cData[i] = -1;

   // fill cData with distances (don't include hydrogen distance
   // metrics (makes layout nicer?)).
   // 
   RDKit::Conformer conf = rdk_mol->getConformer(icurrent_conf);
   int ic_index = 0;
   for (unsigned int iat=1; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR iat_p = (*rdk_mol)[iat];
      if (iat_p->getAtomicNum() != 1) { 
	 RDGeom::Point3D &pos_1 = conf.getAtomPos(iat);
	 // std::cout << "   in 3d conformer: pos " << iat << " is " << pos_1 << std::endl;
	 for (unsigned int jat=0; jat<iat; jat++) {
	    RDKit::ATOM_SPTR jat_p = (*rdk_mol)[jat];
	    if (jat_p->getAtomicNum() != 1) { 
	       RDGeom::Point3D &pos_2 = conf.getAtomPos(jat);
	       RDGeom::Point3D diff = pos_1 - pos_2;

	       // (thanks to JED for useful discussions)
	       ic_index = iat*(iat - 1)/2 + jat;

 	       if (iat < jat)
 		  ic_index = jat*(jat -1)/2 + iat;

	       if (ic_index >= n_items)
		  std::cout << "indexing problem! " << ic_index << " but limit "
			    << n_items << std::endl;
	       if (debug) 
		  std::cout << "mimic: atoms " << iat << " " << jat
			    << " ic_index " << ic_index << " for max " << n_items
			    << " dist " << diff.length() << std::endl;
	       cData[ic_index] = diff.length();
	    }
	 }
      }
   }

   RDDepict::DOUBLE_SMART_PTR dmat(cData);
   
   // AFAICS, the number of conformers before and after calling
   // compute2DCoords() is the same (1).
   // 
   // Does it clear first? Yes, optional parameter clearConfs=true.
   // all atoms, long along x-axis, flip 2 rotatable bonds, try it for
   // 20 samples
   // 
   // int iconf = RDDepict::compute2DCoords(*rdk_mol, NULL, 1, 0, 2, 20);
   //

   int nRB = RDKit::Descriptors::calcNumRotatableBonds(*rdk_mol);

   // The can screw up the coordinates of rdk_mol.  I don't know
   // why. It is not due to distance matrix (dmat) I am pretty sure of
   // that.
    int iconf =
       RDDepict::compute2DCoordsMimicDistMat(*rdk_mol, &dmat, true, true,
					     weight_for_3d_distances, nRB, 200);

   conf = rdk_mol->getConformer(iconf);
   RDKit::WedgeMolBonds(*rdk_mol, &conf);

   if (debug) { // .................... debug ...................
      std::cout << "::::: add_2d_conformer after  compute2DCoords n_atoms: "
		<< rdk_mol->getConformer(0).getNumAtoms()
		<< " n_bonds " << rdk_mol->getNumBonds() << std::endl;
   
      std::cout << ":::::: in add_2d_conformer here are the coords: "
		<< std::endl;
      conf = rdk_mol->getConformer(iconf);
      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 RDKit::ATOM_SPTR at_p = (*rdk_mol)[iat];
	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
	 std::string name = "";
	 try {
	    at_p->getProp("name", name);
	 }
	 catch (const KeyErrorException &kee) {
	 }
	 std::cout << "   " << iat << " " << name << "  " << r_pos << std::endl;
      }
   }

   return iconf;
}

// undelocalise if you can, i.e. there is are 2 deloc bonds to a atom
// (a carbon), make one of these a double (the second, non-N bound)
// and the first (N-C) a single bond.
// 
void
coot::undelocalise(RDKit::RWMol *rdkm) {

   undelocalise_aminos(rdkm);
   undelocalise_methyl_carboxylates(rdkm);
   undelocalise_phosphates(rdkm); 
}

void
coot::undelocalise_aminos(RDKit::RWMol *rdkm) {


   // Plan:
   //
   // 1) Identify a bond with ONEANDAHALF bond order (bond_1)
   // 2) find a nitrogen in the bond
   // 3)   if yes, then is the other atom a carbon
   // 4)      if yes, then find another bond to this carbon that
   //            is deloc (bond_2)
   // 5)         if found, then make bond_1 single, bond_2 double

   bool debug = false;

   int n_bonds = rdkm->getNumBonds();
   RDKit::ROMol::BondIterator bondIt;
   RDKit::ROMol::BondIterator bondIt_inner;
   for(bondIt=rdkm->beginBonds(); bondIt!=rdkm->endBonds(); ++bondIt) {
      if ((*bondIt)->getBondType() == RDKit::Bond::ONEANDAHALF) {
	 // was one of these atoms a Nitrogen?
	 RDKit::Atom *atom_1 = (*bondIt)->getBeginAtom();
	 RDKit::Atom *atom_2 = (*bondIt)->getEndAtom();
	 bool do_it = 0;
	 if (atom_1->getAtomicNum() == 7) {
	    if (atom_2->getAtomicNum() == 6) {
	       do_it = 1;
	    }
	 }
	 if (atom_2->getAtomicNum() == 7) {
	    if (atom_1->getAtomicNum() == 6) {
	       do_it = 1;
	       std::swap(atom_1, atom_2);
	    }
	 }

	 if (do_it) {

	    // atom_1 is a Nitrogen, atom_2 is a Carbon.  Does the
	    // carbon have a bond (not this one) that is delocalised?
	    // 
	    for(bondIt_inner=rdkm->beginBonds(); bondIt_inner!=rdkm->endBonds(); ++bondIt_inner) {
	       RDKit::Atom *atom_1_in = (*bondIt_inner)->getBeginAtom();
	       RDKit::Atom *atom_2_in = (*bondIt_inner)->getEndAtom();
	       if (atom_1_in == atom_2) {
		  if (atom_2_in != atom_1) {
		     if ((*bondIt_inner)->getBondType() == RDKit::Bond::ONEANDAHALF) {

			// this can throw exception:
			// getExplicitValence() called without call to calcExplicitValence()
			// int e_valence_pre = atom_1->getExplicitValence();
			
			// swap them then
			(*bondIt)->setBondType(RDKit::Bond::SINGLE);
			(*bondIt_inner)->setBondType(RDKit::Bond::DOUBLE);
			if (debug)
			   std::cout << "^^^^^^^^^^^^^^^^^ fixed up a bond! 1" << std::endl;
		     }
		  }
	       }
	       if (atom_2_in == atom_2) {
		  if (atom_1_in != atom_1) {
		     if ((*bondIt_inner)->getBondType() == RDKit::Bond::ONEANDAHALF) {
			// swap them then
			(*bondIt)->setBondType(RDKit::Bond::SINGLE);
			(*bondIt_inner)->setBondType(RDKit::Bond::DOUBLE);
			if (debug)
			   std::cout << "^^^^^^^^^^^^^^^^^ fixed up a alternative path bond!"
				     << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
}

void
coot::undelocalise_methyl_carboxylates(RDKit::RWMol *rdkm) {

   // The valence of 2 1/2 on a methyl on one of the oxygens of a carboxylate
   // 
   RDKit::ROMol::BondIterator bondIt;
   RDKit::ROMol::BondIterator bondIt_inner;
   for(bondIt=rdkm->beginBonds(); bondIt!=rdkm->endBonds(); ++bondIt) {
      if ((*bondIt)->getBondType() == RDKit::Bond::ONEANDAHALF) {
	 RDKit::Atom *atom_1 = (*bondIt)->getBeginAtom();
	 RDKit::Atom *atom_2 = (*bondIt)->getEndAtom();

	 if (atom_1->getAtomicNum() == 6) {
	    if (atom_2->getAtomicNum() == 8) {

	       // rename for clarity
	       RDKit::Atom *central_C = atom_1;
	       RDKit::Atom *O1 = atom_2;
	       
	       for(bondIt_inner=rdkm->beginBonds(); bondIt_inner!=rdkm->endBonds(); ++bondIt_inner) {
		  if ((*bondIt_inner)->getBondType() == RDKit::Bond::ONEANDAHALF) {
		     RDKit::Atom *atom_1_in = (*bondIt_inner)->getBeginAtom();
		     RDKit::Atom *atom_2_in = (*bondIt_inner)->getEndAtom();
		     if (atom_1_in == central_C) {
			if (atom_2_in != O1) {
			   if (atom_2_in->getAtomicNum() == 8) {

			      // OK, we have a carbon (atom_1) bonded to two Os via delocs -
			      // the oxygens are atom_2 and atom_2_in
			      // 
			      // rename for clarity
			      //
			      RDKit::Atom *O2 = atom_2_in;
			      
			      // bondIt and bondIt_inner are the bonds that we will ultimately modify
			      // 
			      deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);

			   }
			}
		     }

		     // The central carbon was the other atom?
		     if (atom_2_in == central_C) {
			if (atom_1_in != O1) {
			   if (atom_1_in->getAtomicNum() == 8) {

			      // OK, we have a carbon (atom_1) bonded to two Os via delocs -
			      // the oxygens are atom_2 and atom_2_in
			      // 
			      // rename for clarity
			      //
			      RDKit::Atom *O2 = atom_1_in;
			      deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);
			   } 
			} 
		     }
		  }
	       }
	    }
	 }
	    
	 if (atom_1->getAtomicNum() == 8) {
	    if (atom_2->getAtomicNum() == 6) {
	       // rename for clarity
	       RDKit::Atom *central_C = atom_2;
	       RDKit::Atom *O1 = atom_1;
	       
	       for(bondIt_inner=rdkm->beginBonds(); bondIt_inner!=rdkm->endBonds(); ++bondIt_inner) {
		  if ((*bondIt_inner)->getBondType() == RDKit::Bond::ONEANDAHALF) {
		     RDKit::Atom *atom_1_in = (*bondIt_inner)->getBeginAtom();
		     RDKit::Atom *atom_2_in = (*bondIt_inner)->getEndAtom();
		     if (atom_1_in == central_C) {
			if (atom_2_in != O1) {
			   if (atom_2_in->getAtomicNum() == 8) {

			      // again, we have detected a central carbon bonded to two
			      // oxygens with deloc bonds.
			      RDKit::Atom *O2 = atom_2_in;
			      deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);
			   }
			}
		     }
		     // The central carbon was the other atom?
		     if (atom_2_in == central_C) {
			if (atom_1_in != O1) {
			   if (atom_1_in->getAtomicNum() == 8) {
			      RDKit::Atom *O2 = atom_1_in;
			      deloc_O_check_inner(rdkm, central_C, O1, O2, *bondIt, *bondIt_inner);
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



// fiddle with the bonds in rdkm as needed.
void
coot::deloc_O_check_inner(RDKit::RWMol *rdkm, RDKit::Atom *central_C,
			  RDKit::Atom *O1, RDKit::Atom *O2,
			  RDKit::Bond *b1, RDKit::Bond *b2) {

   RDKit::ROMol::BondIterator bondIt_in_in;
   // OK, so was there something attached to either of the Oxygens?
   // 
   for(bondIt_in_in=rdkm->beginBonds(); bondIt_in_in!=rdkm->endBonds(); ++bondIt_in_in) {
      if ((*bondIt_in_in)->getBondType() == RDKit::Bond::SINGLE) {
	 RDKit::Atom *atom_1_in_in = (*bondIt_in_in)->getBeginAtom();
	 RDKit::Atom *atom_2_in_in = (*bondIt_in_in)->getEndAtom();
				    
	 // check atom_1_in_in vs the first oxygen
	 if (atom_1_in_in == O1) {
	    if (atom_2_in_in != central_C) {

	       // OK, so O1 was bonded to something else
	       // 
	       b1->setBondType(RDKit::Bond::SINGLE);
	       b2->setBondType(RDKit::Bond::DOUBLE);
	    }
	 }
				    
	 // check vs the second oxygen
	 if (atom_1_in_in == O2) {
	    if (atom_2_in_in != central_C) {
	       // OK, so O2 was bonded to something else
	       b1->setBondType(RDKit::Bond::DOUBLE);
	       b2->setBondType(RDKit::Bond::SINGLE);
	    }
	 }

	 

	 // check atom_2_in_in vs the first oxygen
	 if (atom_2_in_in == O1) {
	    if (atom_1_in_in != central_C) {
	       // OK, so O1 was bonded to something else
	       // 
	       b1->setBondType(RDKit::Bond::SINGLE);
	       b2->setBondType(RDKit::Bond::DOUBLE);
	    }
	 }
				    
	 // check atom_2_in_in vs the second oxygen
	 if (atom_2_in_in == O2) { 
	    if (atom_1_in_in != central_C) {
	       // OK, so O2 was bonded to something else
	       b1->setBondType(RDKit::Bond::DOUBLE);
	       b2->setBondType(RDKit::Bond::SINGLE);
	    }
	 }
      }
   }
}

void
coot::undelocalise_phosphates(RDKit::ROMol *rdkm) {

   RDKit::ROMol::AtomIterator ai;
   for(ai=rdkm->beginAtoms(); ai!=rdkm->endAtoms(); ai++) {

      if ((*ai)->getAtomicNum() == 15) {
	 RDKit::Atom *P_at = *ai;
	 int idx_1 = P_at->getIdx();
	 std::vector<RDKit::Bond *> deloc_O_bonds;
	 
	 RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
	 boost::tie(nbrIdx, endNbrs) = rdkm->getAtomNeighbors(P_at);
	 while(nbrIdx != endNbrs) {
	    const RDKit::ATOM_SPTR at = (*rdkm)[*nbrIdx];
	    RDKit::Bond *bond = rdkm->getBondBetweenAtoms(idx_1, *nbrIdx);
	    if (bond) {
	       if (bond->getBondType() == RDKit::Bond::ONEANDAHALF)
		  deloc_O_bonds.push_back(bond);
	    } 
	    ++nbrIdx;
	 }
	 
	 if (deloc_O_bonds.size() == 3) {
	    // make 2 single and one double.  Handle formal charge too?
	    deloc_O_bonds[0]->setBondType(RDKit::Bond::SINGLE);
	    deloc_O_bonds[1]->setBondType(RDKit::Bond::SINGLE);
	    deloc_O_bonds[2]->setBondType(RDKit::Bond::DOUBLE);
	 }

	 if (deloc_O_bonds.size() == 2) {
	    // make 1 single, one double. Handle formal charge too?
	    deloc_O_bonds[0]->setBondType(RDKit::Bond::SINGLE);
	    deloc_O_bonds[1]->setBondType(RDKit::Bond::DOUBLE);
	 } 
      }
   }
}




// now update the atom positions of the conformer iconf in
// rdkit_molecule using the atom positions in residue_p (perhaps this
// should be in rdkit-interface.hh/cc?)
// 
// ignore alt confs.
void coot::update_coords(RDKit::RWMol *mol_p, int iconf, CResidue *residue_p) {

   int n_atoms;
   PPCAtom residue_atoms = NULL;
   residue_p->GetAtomTable(residue_atoms, n_atoms);
   RDKit::Conformer &conf = mol_p->getConformer(iconf);
   for (unsigned int iat=0; iat<n_atoms; iat++) {
      std::string residue_atom_name(residue_atoms[iat]->name);
      CAtom *r_at = residue_atoms[iat];
      for (unsigned int jat=0; jat<n_atoms; jat++) { 
	 RDKit::ATOM_SPTR at_p = (*mol_p)[jat];
	 try {
	    std::string rdkit_atom_name;
	    at_p->getProp("name", rdkit_atom_name);
	    if (rdkit_atom_name == residue_atom_name) {
	       RDGeom::Point3D r_pos(r_at->x, r_at->y, r_at->z);
	       conf.setAtomPos(jat, r_pos);
	    }
	 }
	 catch (KeyErrorException kee) {
// 	    std::cout << "caught no-name for atom exception in update_coords(): "
// 		      <<  kee.what() << std::endl;
	 }
      }
   }
}


// return something interesting in due course
//
// bond_index is the index of the bond to be replace.  atom_index
// is the index of the atom (i.e. one of the atoms in bond of
// bond_index), the R-group of which the atoms are marked for
// deletion and the atom_index atom is changed to "*" for later
// modification with various R-groups.
// 
RDKit::ROMol *
coot::split_molecule(const RDKit::ROMol &mol, int bond_index, int atom_index) {

   RDKit::ROMol *ret_mol = NULL;
   RDKit::RWMol *working_mol = new RDKit::RWMol(mol);

   if (bond_index < 0 || bond_index >= working_mol->getNumBonds()) {
      std::cout << "split_molecule() bad bond index " << bond_index << " vs " << working_mol->getNumBonds()
		<< "  in split_molecule()" << std::endl;
   } else {
      // happy path
      const RDKit::Bond *bond_p = working_mol->getBondWithIdx(bond_index);
      if (bond_p) {
	 int idx_1 = bond_p->getBeginAtomIdx();
	 int idx_2 = bond_p->getEndAtomIdx();
	 if (atom_index == idx_2)
	    std::swap(idx_1, idx_2);
	 if (atom_index != idx_1) {
	    std::cout << "bond index and atom index inconsistent - fail" << std::endl;
	 } else {
	    // OK, so far so good
	    std::queue<int> q;
	    std::vector<int> considered;
	    std::vector<int> R_group_atoms;
	    considered.push_back(idx_1); // not the picked atom;
	    considered.push_back(idx_2); // not the first neighbour of idx_1.
	    // what are the neighbours of idx_1 that is not idx_2?
	    RDKit::ATOM_SPTR at_p = (*working_mol)[idx_1];
	    RDKit::ROMol::OEDGE_ITER current, end;
	    boost::tie(current, end) = working_mol->getAtomBonds(at_p.get());

	    // add some atoms to the queue
	    while (current != end) {
	       RDKit::BOND_SPTR bond= (*working_mol)[*current];
	       int idx = bond->getOtherAtomIdx(idx_1);
	       if (idx != idx_2) {
		  q.push(idx);
		  considered.push_back(idx);
		  R_group_atoms.push_back(idx);
	       }
	       current++;
	    }

	    while (q.size()) {
	       int current_atom_idx = q.front();
	       q.pop();
	       RDKit::ATOM_SPTR at_p = (*working_mol)[current_atom_idx];
	       boost::tie(current, end) = working_mol->getAtomBonds(at_p.get());
	       // std::cout << "current and end: " << current << " " << end << std::endl;
	       while (current != end) {
		  RDKit::BOND_SPTR bond=(*working_mol)[*current];
		  int idx = bond->getOtherAtomIdx(current_atom_idx);
		  if (std::find(considered.begin(),
				considered.end(),
				idx) == considered.end()) {
		     q.push(idx);
		     considered.push_back(idx);
		     R_group_atoms.push_back(idx);
		  }
		  current++;
	       }
	    }

	    std::cout << "R-group atoms has " << R_group_atoms.size() << " atoms"
		      << std::endl;
	    if (1)
	       for (unsigned int iat=0; iat<R_group_atoms.size(); iat++)
		  std::cout << "    " << R_group_atoms[iat] << std::endl;

	    if (R_group_atoms.size()) { 

	       // make a list of atoms to deleted, then delete them
	       // and set the returned molecule pointer to non-null
	       // 
	       std::vector<RDKit::Atom *> atoms_to_be_deleted;
	       for (unsigned int iat=0; iat<R_group_atoms.size(); iat++) {
		  std::cout << "... deleting atom " << R_group_atoms[iat] << std::endl;
		  RDKit::ATOM_SPTR at_p = (*working_mol)[R_group_atoms[iat]];
		  atoms_to_be_deleted.push_back(at_p.get());
	       }
	       for (unsigned int iat=0; iat<R_group_atoms.size(); iat++)
		  working_mol->removeAtom(atoms_to_be_deleted[iat]);

	       ret_mol = working_mol;
	    }
	 }
      }
   }
   
   return ret_mol;
}


// atom_index is the atom in mol that is the "any" atom that is part of the fragment.
//
// This atoms is/becomes the '*' atom of the trial fragment (having atomic number of 0)
//
// Caller deletes the molecules of the vector.
// 
std::vector<RDKit::ROMol *>
coot::join_molecules(const RDKit::ROMol &mol, int atom_index, const RDKit::ROMol &trial_fragment) {

   std::vector<RDKit::ROMol *> v;

   bool i_joining_atom_found = false;
   for (unsigned int iat=0; iat<trial_fragment.getNumAtoms(); iat++) {
      RDKit::ATOM_SPTR at_p = trial_fragment[iat];

      // Was it the '*' atom of the trial_fragment?
      //
      if (at_p->getAtomicNum() == 0) {
	 i_joining_atom_found = true;
	 std::cout << "in join_molecules() found '*' atom " << iat << std::endl;
	 RDKit::RWMol working_mol = mol;
	 // make a coordMap of the working mol atoms for later use in 2d coords generation
	 //
	 RDGeom::INT_POINT2D_MAP coordMap;
	 try {
	    int iconf = 0;
	    int n_conf = working_mol.getNumConformers();
	    std::cout << "------------ working mol has " << n_conf << " conformers" << std::endl;
	    RDKit::Conformer conf = working_mol.getConformer(0);
	    for (unsigned int iwk_at=0; iwk_at<working_mol.getNumAtoms(); iwk_at++) {
	       RDGeom::Point3D p = conf.getAtomPos(iwk_at);
	       coordMap[iwk_at] = RDGeom::Point2D(p.x, p.y);
	    }
	 }
	 catch (...) {
	    std::cout << "caught something on working mol conformer details" << std::endl;
	 } 
	 
	 std::map<unsigned int, unsigned int> atom_idx_map;
	 std::cout << "working_mol has " << working_mol.getNumAtoms() << " atoms." << std::endl;

	 //make copies of the atoms of trial_fragment and add them to working_mol
	 for (unsigned int iat_inner=0; iat_inner<trial_fragment.getNumAtoms(); iat_inner++) {
	    // we don't want to add the '*' atom of the trial_fragment
	    if (iat_inner != iat) { 
	       at_p = trial_fragment[iat_inner];
	       unsigned int new_n_atoms = working_mol.addAtom(at_p); // copies
	       unsigned int working_idx = new_n_atoms; // I think
	       atom_idx_map[iat_inner] = working_idx;
	    }
	 }

	 if (0) { 
	    std::cout << "in join_molecules() iterating over " << trial_fragment.getNumBonds()
		      << " trial_fragment bonds " << std::endl;
	    
	    std::cout << "------ atom idx map ---------" << std::endl;
	    std::map<unsigned int, unsigned int>::const_iterator it;
	    for (it=atom_idx_map.begin(); it!=atom_idx_map.end(); it++) {
	       std::cout << "    atom_name " << it->first << " -> " << it->second << std::endl;
	    } 
	    std::cout << "------ end of atom idx map ---------" << std::endl;
	 }

	 
	 for (unsigned int ibond=0; ibond<trial_fragment.getNumBonds(); ibond++) {
	    const RDKit::Bond *bond = trial_fragment.getBondWithIdx(ibond);
	    unsigned int itb = bond->getBeginAtomIdx();
	    unsigned int ite = bond->getEndAtomIdx();
	    std::map<unsigned int, unsigned int>::const_iterator it_b = atom_idx_map.find(itb);
	    std::map<unsigned int, unsigned int>::const_iterator it_e = atom_idx_map.find(ite);
	    if (it_b != atom_idx_map.end()) {
	       if (it_e != atom_idx_map.end()) {
		  unsigned int iwb = it_b->second;
		  unsigned int iwe = it_e->second;
		  working_mol.addBond(iwb, iwe, bond->getBondType());
	       } else {
		  // OK, so ite could have been the '*' atom
		  if (ite == iat) {
		     unsigned int iwb = it_b->second;
		     working_mol.addBond(iwb, atom_index, bond->getBondType());
		  }
	       }
	    } else { 
	       // OK, so itb could have been the '*' atom
	       if (itb == iat) {
		  if (it_e != atom_idx_map.end()) { // just for safety sake
		     unsigned int iwe = it_e->second;
		     working_mol.addBond(atom_index, iwe, bond->getBondType());
		  }
	       }
	    }
	 }

	 rdkit_mol_sanitize(working_mol);
	    
	 RDDepict::compute2DCoords(working_mol, &coordMap);
	 RDKit::ROMol *romol = new RDKit::ROMol(working_mol);
	 v.push_back(romol);
      }

      
   }
   return v;
}


#endif // MAKE_ENTERPRISE_TOOLS   
