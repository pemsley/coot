

#include "lbg.hh"

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
	    m.addBond(bond);
	 }
      }
   }
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
