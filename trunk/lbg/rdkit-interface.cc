
#ifdef HAVE_GOOCANVAS

#include "coot-utils.hh"
#include "lbg.hh"

#ifdef MAKE_ENTERPRISE_TOOLS

#include "rdkit-interface.hh"


// Caller deletes
// 
RDKit::RWMol
lbg_info_t::rdkit_mol(const widgeted_molecule_t &mol) const {

   RDKit::RWMol m;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      if (! mol.atoms[iat].is_closed()) { 
	 RDKit::Atom *at = new RDKit::Atom;
	 at->setAtomicNum(tbl->getAtomicNumber(mol.atoms[iat].element));

	 // add the name to at too here (if you can).
	 // 
	 try {
	    at->setProp("name", std::string(mol.atoms[iat].get_atom_name()));
	 }
	 catch (std::exception rte) {
	    std::cout << rte.what() << std::endl;
	 }
	 m.addAtom(at);
      }
   }

   for (unsigned int ib=0; ib<mol.bonds.size(); ib++) {
      if (! mol.bonds[ib].is_closed()) {
	 RDKit::Bond::BondType type = convert_bond_type(mol.bonds[ib].get_bond_type());
	 RDKit::Bond *bond = new RDKit::Bond(type);
	 int idx_1 = mol.bonds[ib].get_atom_1_index();
	 int idx_2 = mol.bonds[ib].get_atom_2_index();
	 if (!mol.atoms[idx_1].is_closed() && !mol.atoms[idx_2].is_closed()) { 
	    bond->setBeginAtomIdx(idx_1);
	    bond->setEndAtomIdx(  idx_2);
	    if (type == RDKit::Bond::AROMATIC) { 
	       bond->setIsAromatic(true);
	       m[idx_1]->setIsAromatic(true);
	       m[idx_2]->setIsAromatic(true);
	    } 
	    m.addBond(bond);
	 }
      }
   }
   return m;
}

// This can throw an runtime_error exception (residue not in
// dictionary).
//
// Return an RDKit::RWMol with one conformer (perhaps we should do a
// conformer for each alt conf).
// 
RDKit::RWMol
coot::rdkit_mol(CResidue *residue_p, const coot::protein_geometry &geom) {

   std::string res_name = residue_p->GetResName();
   
   std::pair<bool, coot::dictionary_residue_restraints_t> p = 
      geom.get_monomer_restraints_at_least_minimal(res_name);
   if (! p.first) {

      std::string m = "rdkit_mol(): residue type ";
      m += res_name;
      m += " not in dictionary";
      throw(std::runtime_error(m));

   } else { 
      return rdkit_mol(residue_p, p.second);
   } 
}

RDKit::RWMol
coot::rdkit_mol(CResidue *residue_p,
		const coot::dictionary_residue_restraints_t &restraints) {

   RDKit::RWMol m;
   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   // this is so that we don't add multiple copies of an atom with
   // the same name (that is, only add the first atom of a given
   // atom name in a residue with alt confs).
   std::vector<std::string> added_atom_names;
   std::map<std::string, int> atom_index;
   int current_atom_id = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      // only add the atom if the atom_name is not in the list of
      // already-added atom names.
      if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name) == added_atom_names.end()) {
	 RDKit::Atom *at = new RDKit::Atom;
	 try {
	    at->setAtomicNum(tbl->getAtomicNumber(coot::util::capitalise(coot::util::remove_leading_spaces(residue_atoms[iat]->element))));
	    // for debugging, delete if you feel like it.
	    // if (std::string(residue_atoms[iat]->element) == " H")
	    //    at->setAtomicNum(tbl->getAtomicNumber(std::string("XXX")));
	       
	    at->setProp("name", std::string(residue_atoms[iat]->name));
	    m.addAtom(at);
	    if (0) 
	       std::cout << "adding atom with name " << atom_name << " to added_atom_names"
			 << " which is now size " << added_atom_names.size() << std::endl;
	    added_atom_names.push_back(atom_name);
	    atom_index[atom_name] = current_atom_id;
	    current_atom_id++; // for next round
	 }
	 catch (std::exception rte) {
	    std::cout << rte.what() << std::endl;
	 } 
      }
   }

   for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) { 
      RDKit::Bond::BondType type = convert_bond_type(restraints.bond_restraint[ib].type());
      RDKit::Bond *bond = new RDKit::Bond(type);
	    
      std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
      std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
      std::string ele_1 = restraints.element(atom_name_1);
      std::string ele_2 = restraints.element(atom_name_2);
      int idx_1 = -1; // unset
      int idx_2 = -1; // unset
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 if (! residue_atoms[iat]->Ter) { 
	    std::string atom_name(residue_atoms[iat]->name);
	    if (atom_name == atom_name_1)
	       if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name)
		   != added_atom_names.end())
		  idx_1 = iat;
	    if (atom_name == atom_name_2)
	       if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name)
		   != added_atom_names.end())
		  idx_2 = iat;
	 }
      }
      if (idx_1 != -1) { 
	 if (idx_2 != -1) {	 
	    if (type == RDKit::Bond::AROMATIC) { 
	       bond->setIsAromatic(true);
	       m[idx_1]->setIsAromatic(true);
	       m[idx_2]->setIsAromatic(true);

	    }
	    bond->setBeginAtomIdx(idx_1);
	    bond->setEndAtomIdx(  idx_2);
	    m.addBond(bond);
	 } else {
	    if (ele_2 != " H") 
	       std::cout << "WARNING:: oops, bonding in making rdkit mol "
			 << "failed to get atom index idx_2 for atom name: "
			 << atom_name_2 << " ele :" << ele_2 << ":" << std::endl;
	 }
      } else {
	 if (ele_1 != " H") 
	    std::cout << "WARNING:: oops, bonding in making rdkit mol "
		      << "failed to get atom index idx_1 for atom name: "
		      << atom_name_1 << " ele :" << ele_1 << ":" << std::endl;
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
      RDKit::Bond *bond = new RDKit::Bond(type);
	    
      std::string atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
      std::string atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();
      std::string ele_1 = restraints.element(atom_name_1);
      std::string ele_2 = restraints.element(atom_name_2);
      int idx_1 = -1; // unset
      int idx_2 = -1; // unset
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 if (! residue_atoms[iat]->Ter) { 
	    std::string atom_name(residue_atoms[iat]->name);
	    if (atom_name == atom_name_1)
	       if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name)
		   != added_atom_names.end())
		  idx_1 = iat;
	    if (atom_name == atom_name_2)
	       if (std::find(added_atom_names.begin(), added_atom_names.end(), atom_name)
		   != added_atom_names.end())
		  idx_2 = iat;
	 }
      }
      if (idx_1 != -1) { 
	 if (idx_2 != -1) {	 
	    if (type == RDKit::Bond::AROMATIC) { 
   
	       // special edge case for aromatic ring N that may need an H attached for
	       // kekulization (depending on energy type).

	       if (m[idx_1]->getAtomicNum() == 7) {
		  if (std::find(Hs_added_list.begin(), Hs_added_list.end(), idx_1) == Hs_added_list.end()) { 
		     std::string n = add_H_to_ring_N_as_needed(&m, idx_1, atom_name_1, restraints);
		     std::cout << "testing 1 idx_1 " << idx_1 << " idx_2 " << idx_2 << " n was :" << n << ":"
			       << std::endl;
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

   RDKit::MolOps::cleanUp(m);

   // OK, so cleanUp() doesn't fix the N charge problem our prodrg molecule
   if (0) { // debug, formal charges
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
      std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
		<< std::endl;
   }

   RDKit::MolOps::sanitizeMol(m);

   std::cout << "in constructing rdk molecule now adding a conf" << std::endl;
   RDKit::Conformer *conf = new RDKit::Conformer(added_atom_names.size());
   conf->set3D(true);
      
   // Add positions to the conformer (only the first instance of an
   // atom with a particular atom name).
   // 
   for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::map<std::string, int>::const_iterator it = atom_index.find(atom_name);
      if (it != atom_index.end()) {
	 RDGeom::Point3D pos(residue_atoms[iat]->x,
			     residue_atoms[iat]->y,
			     residue_atoms[iat]->z);
	 conf->setAtomPos(it->second, pos);
	 if (0) 
	    std::cout << "in construction of rdkit mol: making a conformer "
		      << iat << " " << it->second << " " << atom_name << " "
		      << pos << std::endl;
      } 
   }
   m.addConformer(conf);
   std::cout << "ending construction of rdkit mol: n_atoms " << m.getNumAtoms() << std::endl;
   return m;
}


RDKit::Bond::BondType
lbg_info_t::convert_bond_type(const lig_build::bond_t::bond_type_t &t) const {

   // There are lots more in RDKit::Bond::BondType!
   // 
   RDKit::Bond::BondType bt = RDKit::Bond::UNSPECIFIED;
   if (t == lig_build::bond_t::SINGLE_BOND)
      bt = RDKit::Bond::SINGLE;
   if (t == lig_build::bond_t::DOUBLE_BOND)
      bt = RDKit::Bond::DOUBLE;
   if (t == lig_build::bond_t::TRIPLE_BOND)
      bt = RDKit::Bond::TRIPLE;
   return bt;
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


// this can throw a std::exception
// 
std::string
lbg_info_t::get_smiles_string_from_mol_rdkit() const {

   RDKit::RWMol rdkm = rdkit_mol(mol);
   RDKit::ROMol *rdk_mol_with_no_Hs = RDKit::MolOps::removeHs(rdkm);
   std::string s = RDKit::MolToSmiles(*rdk_mol_with_no_Hs);
   delete rdk_mol_with_no_Hs;

   return s;
}


lig_build::molfile_molecule_t
coot::make_molfile_molecule(const RDKit::ROMol &rdkm, int iconf) {

   lig_build::molfile_molecule_t mol;
   int n_conf  = rdkm.getNumConformers();

   if (n_conf) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
      
      RDKit::Conformer conf = rdkm.getConformer(iconf);
      int n_conf_atoms = conf.getNumAtoms();
      int n_mol_atoms = rdkm.getNumAtoms();

      std::cout << "in make_molfile_molecule() conformer has " << n_conf_atoms
		<< " atoms" << std::endl;
      
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
	 lig_build::molfile_bond_t mol_bond(idx_1, idx_2, bt);
	 mol.add_bond(mol_bond);
      }
   }
   
   return mol;
}

// returns NULL on fail. Caller deletes.
//
CResidue *
coot::make_residue(const RDKit::ROMol &rdkm, int iconf) {

   CResidue *residue_p = NULL;
   lig_build::molfile_molecule_t mol = coot::make_molfile_molecule(rdkm, iconf);

   // now convert mol to a CResidue *
   // 
   if (mol.atoms.size()) {
      residue_p = new CResidue;
      residue_p->seqNum = 1;
      CChain *chain_p = new CChain;
      chain_p->SetChainID("");
      chain_p->AddResidue(residue_p);
      for (unsigned int iat=0; iat<mol.atoms.size(); iat++) { 
	 CAtom *at = new CAtom;
	 at->SetAtomName(mol.atoms[iat].name.c_str());
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
// BL says:: dont have it now (old rdkit?)
#ifndef WINDOWS_MINGW
   RDKit::MolOps::assignRadicals(*rdkm);
#endif // WINDOWS
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
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR at_p = (*rdkm)[iat];
      if (at_p->getAtomicNum() == 7) {
	 
	 int e_valence = at_p->getExplicitValence();

	 std::cout << " atom N has explicit valence: " << e_valence << std::endl;

	 if (e_valence == 4) { 
	    std::cout << ".......... assign_formal_charges: found one! "
		      << at_p << std::endl;
	    at_p->setFormalCharge(1);
	 }
      }
   }
}

// a wrapper for the above, matching hydrogens names to the
// dictionary.  Add atoms to residue_p, return success status.
// 
bool
coot::add_hydrogens_with_rdkit(CResidue *residue_p,
			      const coot::dictionary_residue_restraints_t &restraints) {

   bool r = 0;
   if (residue_p) {
      const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
      try {
	 RDKit::RWMol m_no_Hs = rdkit_mol(residue_p, restraints);
	 int n_mol_atoms = m_no_Hs.getNumAtoms();
	 std::cout << ".......  pre H addition mol has " << n_mol_atoms
		   << " atoms" << std::endl;
	 
	 for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	    RDKit::ATOM_SPTR at_p = m_no_Hs[iat];
	    at_p->calcImplicitValence(true);
	 }
	 
	 RDKit::ROMol *m_pre = RDKit::MolOps::addHs(m_no_Hs, 0, 1);
	 RDKit::RWMol m(*m_pre);
	 std::cout << "....... post H addition mol has " << m_pre->getNumAtoms()
		   << std::endl;
	 int n_atoms_new = m_pre->getNumAtoms();
	 int n_conf = m_pre->getNumConformers();

	 double vdwThresh=10.0;
	 int confId = 0;
	 bool ignoreInterfragInteractions=true;
	 int maxIters = 500;
	 
 	 ForceFields::ForceField *ff =
 	    RDKit::UFF::constructForceField(m, vdwThresh, confId,
 					    ignoreInterfragInteractions);

 	 for (unsigned int iat=0; iat<n_mol_atoms; iat++)
 	    ff->fixedPoints().push_back(iat);

 	 ff->initialize();
 	 int res=ff->minimize(maxIters);
	 std::cout << "rdkit minimize() returns " << res << std::endl;
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
		     H_names_already_added.push_back(name);
		  
		     int n = at_p->getAtomicNum();
		     std::string element = tbl->getElementSymbol(n);
		  
		     CAtom *at = new CAtom;
		     at->SetAtomName(name.c_str());
		     // at->SetElementName(element.c_str()); // FIXME?
		     at->SetElementName(" H");
		     at->SetCoordinates(r_pos.x, r_pos.y, r_pos.z, 1.0, 30.0);
		     at->Het = 1;
		     residue_p->AddAtom(at);
		     r = 1;
		  }
	       }
	    }
	 }
	 
	 // delete m;
      }
      catch (std::runtime_error e) {
	 std::cout << e.what() << std::endl;
      }
   }
   return r;
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

   int icurrent_conf = 0; // the conformer number from which the
                          // distance matrix is generated.  Should this
			  // be passed?

   int n_conf  = rdk_mol->getNumConformers();
   if (n_conf == 0) {
      std::cout << "WARNING:: no conformers in add_2d_conformer() - aborting"
		<< std::endl;
   }

   int n_mol_atoms = rdk_mol->getNumAtoms();   

   if (0) 
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

   double *cData = new double[n_items];
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
	       if (0) 
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
   int iflip = 3;
   int iconf =
      RDDepict::compute2DCoordsMimicDistMat(*rdk_mol, &dmat, 1, 1, weight_for_3d_distances, iflip, 200);

   if (0) { // .................... debug ...................
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
	 catch (KeyErrorException kee) {
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

   // Plan:
   //
   // 1) Identify a bond with ONEANDAHALF bond order (bond_1)
   // 2) find a nitrogen in the bond
   // 3)   if yes, then is the other atom a carbon
   // 4)      if yes, then find another bond to this carbon that
   //            is deloc (bond_2)
   // 5)         if found, then make bond_1 single, bond_2 double
    

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
			// swap them then
			(*bondIt)->setBondType(RDKit::Bond::SINGLE);
			(*bondIt_inner)->setBondType(RDKit::Bond::DOUBLE);
			// std::cout << "^^^^^^^^^^^^^^^^^ fixed up a bond! 1" << std::endl;
		     }
		  }
	       }
	       if (atom_2_in == atom_2) {
		  if (atom_1_in != atom_1) {
		     if ((*bondIt_inner)->getBondType() == RDKit::Bond::ONEANDAHALF) {
			// swap them then
			(*bondIt)->setBondType(RDKit::Bond::SINGLE);
			(*bondIt_inner)->setBondType(RDKit::Bond::DOUBLE);
			// std::cout << "^^^^^^^^^^^^^^^^^ fixed up a bond! 2" << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
}


#endif // MAKE_ENTERPRISE_TOOLS   
#endif // HAVE_GOOCANVAS
