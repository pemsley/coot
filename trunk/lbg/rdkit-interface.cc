
#include "coot-utils.hh"
#include "lbg.hh"
#include "enterprise.hh"

#ifdef MAKE_ENTERPRISE_TOOLS

void
lbg_info_t::update_statusbar_smiles_string() const {

   std::string status_string;
   try {
      std::string s = get_smiles_string_from_mol();
      std::cout << "SMILES string: " << s << std::endl;
      status_string = " SMILES:  ";
      status_string += s;

   }
   catch (std::exception rte) {
      std::cout << rte.what() << std::endl;
   }
   if (lbg_statusbar) { 
      guint statusbar_context_id =
	 gtk_statusbar_get_context_id(GTK_STATUSBAR(lbg_statusbar), status_string.c_str());
      gtk_statusbar_push(GTK_STATUSBAR(lbg_statusbar),
			 statusbar_context_id,
			 status_string.c_str());
   }
}

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
coot::rdkit_mol(CResidue *residue_p, const protein_geometry &geom) {


   std::string res_name = residue_p->GetResName();
   
   std::pair<bool, dictionary_residue_restraints_t> p = 
      geom.get_monomer_restraints(res_name);
   if (! p.first) {

      throw(std::runtime_error("residue type not in dictionary"));

   } else { 
   
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
	       at->setAtomicNum(tbl->getAtomicNumber(coot::util::remove_leading_spaces(residue_atoms[iat]->element)));
	       m.addAtom(at);
	       added_atom_names.push_back(atom_name);
	       atom_index[atom_name] = current_atom_id;
	    }
	    catch (std::exception rte) {
	       std::cout << rte.what() << std::endl;
	    } 
	    current_atom_id++; // for next round
	 }
      }
      for (unsigned int ib=0; ib<p.second.bond_restraint.size(); ib++) { 
	 RDKit::Bond::BondType type = convert_bond_type(p.second.bond_restraint[ib].type());
	 RDKit::Bond *bond = new RDKit::Bond(type);
	    
	 std::string atom_name_1 = p.second.bond_restraint[ib].atom_id_1_4c();
	 std::string atom_name_2 = p.second.bond_restraint[ib].atom_id_2_4c();
	 int idx_1 = -1;
	 int idx_2 = -1;
	 for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
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
	 if (idx_1 != -1) { 
	    if (idx_2 != -1) {	 
	       if (type == RDKit::Bond::AROMATIC) { 
		  bond->setIsAromatic(true);
		  m[idx_1]->setIsAromatic(true);
		  m[idx_2]->setIsAromatic(true);
	       }
	       bond->setBeginAtomIdx(idx_1);
	       bond->setEndAtomIdx(  idx_2);
	       std::cout << "Adding rdkit bond with atom indices " << idx_1 << " and " << idx_2 << std::endl;
	       m.addBond(bond);
	    }
	 }
      }

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
	 } 
      }
      m.addConformer(conf);
      return m;
   } 
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
   std::cout << "created RDKit bond type " << bt;
   if (bt == RDKit::Bond::AROMATIC)
      std::cout << " (aromatic)";
   std::cout << std::endl;
   return bt;
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
      for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
	 RDKit::ATOM_SPTR at_p = rdkm[iat];
	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
	 clipper::Coord_orth pos(r_pos.x, r_pos.y, r_pos.z);
	 int n = at_p->getAtomicNum();
	 std::string element = tbl->getElementSymbol(n);
	 int charge = at_p->getFormalCharge();
	 lig_build::molfile_atom_t mol_atom(pos, element, "");
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

// tweak rdkmol
int
coot::add_2d_conformer(RDKit::ROMol *rdk_mol) {

   int iconf = RDDepict::compute2DCoords(*rdk_mol);
   return iconf;
} 

#endif // MAKE_ENTERPRISE_TOOLS   


