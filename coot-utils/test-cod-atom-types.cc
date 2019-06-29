/* lidia-core/test-cod-atom-types.cc
 * 
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

#include <limits>
#ifndef MAKE_ENHANCED_LIGAND_TOOLS
int main(int argc, char **argv) {return 0;}
#else 
#include "utils/coot-utils.hh"

#include "geometry/residue-and-atom-specs.hh"
#include "lidia-core/cod-atom-types.hh"
#include "lidia-core/rdkit-interface.hh"
#include "lidia-core/bond-record-container-t.hh"
#include "lidia-core/primes.hh"

#include "atom-selection-container.hh"

mmdb::Residue *
get_residue(coot::residue_spec_t residue_spec, mmdb::Manager *mol) {

   // Fill this FIXME
   return nullptr;

}

mmdb::Residue *
get_first_residue(mmdb::Manager *mol) {

   // Fill this FIXME
   return nullptr;

}


// rdkit_mol is not const because there is no const beginAtoms() operator.
// 
void write_types(RDKit::RWMol &rdkm) {

   if (false) {
      for (unsigned int iat=0; iat<rdkm.getNumAtoms(); iat++) {
	 try {
	    std::string name;
	    RDKit::ATOM_SPTR at_p = rdkm[iat];
	    at_p->getProp("name", name);
	    std::cout << iat << "   \"" << name << "\"\n";
	 }
	 catch (const KeyErrorException &err) {
	    std::cout << "caught no-name exception in rdkit_mol H-block" << std::endl;
	 }
      }
   }

   cod::atom_types_t t;
   std::vector<cod::atom_type_t> v = t.get_cod_atom_types(rdkm);

   const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
   for (unsigned int iat=0; iat<v.size(); iat++) {
      std::string name;
      try {
	 RDKit::ATOM_SPTR at_p = rdkm[iat];
	 at_p->getProp("name", name);

	 int n = at_p->getAtomicNum();
	 std::string atom_ele = tbl->getElementSymbol(n);
      
	 // Acedrg writes out: index element name COD-type, we should match that.
	 // 
	 //
	 if (false)
	    std::cout << " " << std::right << std::setw(3) << iat << "    "
		      << atom_ele << "    " << name << "     " << v[iat].level_4
		      << "\n";

	 // Now I want to check the tier-2 types too
	 // 
	 std::cout << " " << std::right << std::setw(3) << iat << "    "
		   << atom_ele << "    " << name << "     "
		   << v[iat].level_2.string() << "  "
		   << v[iat].level_4 << "  "
		   << "\n";
      }
      catch (const KeyErrorException &err) { }
   }
}

void molecule_from_ccd_pdbx(const std::string &comp_id,
			    const std::string &file_name) {

   coot::protein_geometry geom;
   int read_number = 0;
   geom.init_refmac_mon_lib(file_name, read_number++);

   bool idealised_flag = true;
   int imol = 0; // dummy
   mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, imol, idealised_flag);

   if (! mol) {
      std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
   } else {
      mmdb::Residue *residue_p = get_first_residue(mol);

      if (! residue_p) {
	 // pretty strange
	 std::cout << "Null residue from mol from mol_from_dictionary() for "
		   << comp_id << std::endl;
      } else {
	 try {
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, geom);
	    write_types(rdkm);
	 }

	 catch(const std::runtime_error &rte) {
	    std::cout << "error::" << rte.what() << std::endl;
	 }
      }
   }
}

void molecule_from_comp_id(const std::string &comp_id) {
   try {
      coot::protein_geometry geom;
      // geom.init_standard();

      // std::string three_letter_code = "001"; // nice test
      // std::string three_letter_code = "0PY"; // simple
      // std::string three_letter_code = "06C"; // simplest

      int imol = 0; // dummy
      bool idealised_flag = true;
      mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, imol, idealised_flag);

      if (! mol) {
	 std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
      } else {
	 
	 mmdb::Residue *residue_p = get_first_residue(mol);

	 if (! residue_p) {
	    // pretty strange
	    std::cout << "Null residue from mol from mol_from_dictionary() for "
		      << comp_id << std::endl;
	 } else { 

	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol, geom);
	    coot::debug_rdkit_molecule(&rdkm);

	    write_types(rdkm);
	       
	 }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "test-cod-atom-types caught exception: " << rte.what() << std::endl;
   }
}

void molecule_from_SMILES(const std::string &smiles_string) {

   std::string s = smiles_string;
   // std::string s = "COc1ccc(cc1O[C@H]1C[C@@H]2CC[C@H]1C2)C1CNC(=O)NC1";
   // s = "CC12CC[C@@]3(C1)CC(=C)CC3C(C2C(C)(C)C)C(C)(C)C"; // bird
   // s = "O[V]1(O)(O)ONC(=[O]1)c1ccccc1";
   // s = "C#O";
   // s = "[H]/N=C(\\N)/NCCC[C@@H](C(=O)N[C@@H](CCc1ccccc1)/C=C/S(=O)(=O)c2ccccc2)NC(=O)N3CC[NH+](CC3)C";

   try {
   
      RDKit::RWMol *rdkmp = RDKit::SmilesToMol(s);
      if (rdkmp) {
	 coot::mogulify_mol(*rdkmp);

	 { 
	    bool includeStereo = true;
	    bool kekulize = false;
	    std::string sdf_file_name = "smiles-string.mdl";
	    int confId = -1;
	    RDDepict::compute2DCoords(*rdkmp, NULL, true);
	    std::cout << "confId: " << confId << std::endl;
	    RDKit::MolToMolFile(*rdkmp, sdf_file_name, includeStereo, confId, kekulize);
	 }
   
	 coot::debug_rdkit_molecule(rdkmp);
	 RDKit::RWMol rdkm(*rdkmp);

	 cod::atom_types_t t;
	 std::vector<cod::atom_type_t> v = t.get_cod_atom_types(rdkm);

	 std::cout << "PE-TYPES:: -------- got " << v.size() << " atoms " << std::endl;
	 for (unsigned int i=0; i<v.size(); i++)
	    std::cout << "   " << i << " " << v[i].level_4 << "" << std::endl;
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
   }

}

void
write_bonds_by_type(RDKit::RWMol *rdkm, const std::vector<cod::atom_type_t> &v) {

   if (rdkm) {
      unsigned int n_atoms = rdkm->getNumAtoms();
      if (n_atoms == v.size()) {
	 unsigned int n_bonds = rdkm->getNumBonds();
	 for (unsigned int ib=0; ib<n_bonds; ib++) {
	    const RDKit::Bond *bond_p = rdkm->getBondWithIdx(ib);
	    int idx_1 = bond_p->getBeginAtomIdx();
	    int idx_2 = bond_p->getEndAtomIdx();
	    std::string t1 = v[idx_1].level_4;
	    std::string t2 = v[idx_2].level_4;
	    if (t1 > t2) {
	       // std::swap(t1, t2);
	    }
	    std::cout << "BOND level-4 " << t1 << "     " << t2 << std::endl;
	    // std::cout << "BOND level-3 " << v[idx_1].level_3 << "     " << v[idx_2].level_3 << std::endl;
	    // std::cout << "BOND level-2 " << v[idx_1].level_2.string() << "     "
	    // << v[idx_2].level_2.string() << std::endl;
	 }
      }
   }
}

void molecule_from_mdl_mol(const std::string &file_name) {

   try {
      RDKit::RWMol *rdkm = RDKit::MolFileToMol(file_name);
      if (rdkm) {
	 cod::atom_types_t t;
	 std::vector<cod::atom_type_t> v = t.get_cod_atom_types(*rdkm);
	 std::cout << "PE-TYPES:: -------- " << v.size() << " atoms " << " from "
		   << file_name << std::endl;
	 for (unsigned int i=0; i<v.size(); i++)
	    std::cout << "   " << i << " " << v[i].level_4 << "" << std::endl;
	 write_bonds_by_type(rdkm, v);
      } else {
	 std::cout << "WARNING:: file " << file_name << " null molecule" << std::endl;
      }
   }
   catch (const RDKit::MolSanitizeException &e) {
      std::cout << "WARNING:: file " << file_name << " "  << e.what() << std::endl;
   }
}

// tables_dir_name is xxx/tables/AlllOrgBondTables
void read_tables(const std::string &tables_dir_name, bool make_db = false) {

   cod::bond_record_container_t brc(tables_dir_name);

   if (false) {  // debug
      unsigned int n_records = 10;
      for (unsigned int i=0; i<n_records; i++)
	 std::cout << "   " << brc.bonds[i] << std::endl;
   }

   if (make_db)
      brc.make_db("at.db");
   else
      brc.write("atom-indices.tab", "bonds.tab"); // writes atom-indices.tab too atm

   if (false) {
      std::cout << "----------- reading " << std::endl;
      cod::bond_record_container_t brc_read;
      brc_read.read("atom-indices.tab", "bonds.tab");
      std::cout << "----------- done reading " << std::endl;
      brc_read.check();
   }
}

void
validate(const std::string &comp_id,
	 const std::string &chain_id,
	 int res_no,
	 const std::string &pdb_file_name,
	 const std::string &cif_file_name,
	 const std::string &acedrg_install_dir) {

   cod::bond_record_container_t brc; // perhaps we need a better constructor here

   // this uses bonds.tab, only the level-4 types, but this is not
   // good enough to find all the hits with generalization, we need
   // information in the acedrg tables/*.table files
   // brc.read("atom-indices.tab", "bonds.tab");

   std::string dir = acedrg_install_dir + "/share/acedrg/tables/allOrgBondTables";
   brc.read_acedrg_table_dir(dir);
   
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);

   coot::residue_spec_t residue_spec(chain_id, res_no, "");

   mmdb::Residue *res = get_residue(residue_spec, asc.mol);

   if (res) { 

      coot::protein_geometry geom;
      int read_number = 0;
      geom.init_refmac_mon_lib(cif_file_name, read_number++);

      int imol = 0; // dummy
      std::pair<bool, coot::dictionary_residue_restraints_t> p = 
	 geom.get_monomer_restraints(comp_id, imol);

      if (p.first) {

	 coot::dictionary_residue_restraints_t rest = p.second;
	 brc.validate(res, rest);
      }
   }
}

void
test_primes() {

   cod::primes primes(200);
   std::vector<unsigned int> pr = primes.get_primes();
   for (unsigned int i=0; i<pr.size(); i++) { 
      std::cout << "   " << pr[i];
   }
   std::cout << "" << std::endl;
   std::cout << "max unsigned int: " << std::numeric_limits<unsigned int>::max()
	     << std::endl;
}

void proc_mols(const std::string &mol_dir) {

   std::vector<std::string> files = coot::util::glob_files(mol_dir, "*.mol");
   for (unsigned int i=0; i<files.size(); i++) {
      molecule_from_mdl_mol(files[i]);
   }
}


int main(int argc, char **argv) {

   int status = 0;

   if (argc > 1) {
      if (argc == 2) {
	 std::string s = argv[1];

	 if (s == "primes") {
	    test_primes();
	 } else {
	    if (s.length() == 3)
	       molecule_from_comp_id(s);
	    else
	       molecule_from_SMILES(s);
	 }
      }
      
      if (argc == 3) {

	 std::string comp_id   = argv[1];
	 std::string file_name = argv[2];

	 if (comp_id == "tables") {
	    bool make_db = true;
	    read_tables(file_name, make_db); // dir-name in this case
	 } else {
	    if (comp_id == "mol-dir") {
	       std::string mol_dir = file_name;
	       proc_mols(mol_dir);
	    } else {
	       molecule_from_ccd_pdbx(comp_id, file_name);
	    }
	 }
      }

      if (argc == 8) {

	 std::string residue = argv[1];
	 if (residue == "residue") {
	    std::string comp_id       = argv[2];
	    std::string chain_id      = argv[3];
	    std::string res_no_str    = argv[4];
	    std::string pdb_file_name = argv[5];
	    std::string cif_file_name = argv[6];
	    std::string acedrg_dir    = argv[7];

	    try {
	       int res_no = coot::util::string_to_int(res_no_str);
	       validate(comp_id, chain_id, res_no, pdb_file_name,
			cif_file_name, acedrg_dir);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "Failure:: " << rte.what() << std::endl;
	    }
	 }
      }
   }
   return status;
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS

