
#ifndef MAKE_ENHANCED_LIGAND_TOOLS
int main(int argc, char **argv) {return 0;}
#else 
#include "cod-types.hh"
#include "rdkit-interface.hh"
#include "coot-utils/coot-coord-utils.hh"

// rdkit_mol is not const because there is no const beginAtoms() operator.
// 
void write_types(RDKit::RWMol &rdkm) {

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

   std::vector<std::string> v = cod::get_cod_atom_types(rdkm);

   for (unsigned int iat=0; iat<v.size(); iat++) {
      std::string name;
      try {
	 RDKit::ATOM_SPTR at_p = rdkm[iat];
	 at_p->getProp("name", name);
      }
      catch (const KeyErrorException &err) { }
      std::cout << " " << std::right << std::setw(3) << iat << " " << name << " " << v[iat]
		<< "\n";
   }
}

void molecule_from_ccd_pdbx(const std::string &comp_id,
			    const std::string &file_name) {

   coot::protein_geometry geom;
   int read_number = 0;
   geom.init_refmac_mon_lib(file_name, read_number++);

   bool idealised_flag = true;
   mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, idealised_flag);

   if (! mol) {
      std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
   } else {
      mmdb::Residue *residue_p = coot::util::get_first_residue(mol);

      if (! residue_p) {
	 // pretty strange
	 std::cout << "Null residue from mol from mol_from_dictionary() for "
		   << comp_id << std::endl;
      } else {
	 try {
	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, geom);
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
	 
      bool idealised_flag = true;
      mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, idealised_flag);

      if (! mol) {
	 std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
      } else {
	 
	 mmdb::Residue *residue_p = coot::util::get_first_residue(mol);

	 if (! residue_p) {
	    // pretty strange
	    std::cout << "Null residue from mol from mol_from_dictionary() for "
		      << comp_id << std::endl;
	 } else { 

	    RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, geom);
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
   s = "[H]/N=C(\\N)/NCCC[C@@H](C(=O)N[C@@H](CCc1ccccc1)/C=C/S(=O)(=O)c2ccccc2)NC(=O)N3CC[NH+](CC3)C";
   
   RDKit::RWMol *rdkmp = RDKit::SmilesToMol(s);
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
	       
   std::vector<std::string> v = cod::get_cod_atom_types(rdkm);

   std::cout << "PE-TYPES:: -------- got " << v.size() << " atoms " << std::endl;
   for (unsigned int i=0; i<v.size(); i++)
      std::cout << "   " << i << " " << v[i] << "" << std::endl;

}

int main(int argc, char **argv) {

   int status = 0;

   if (argc > 1) {
      if (argc == 2) {
	 std::string s = argv[1];
	 if (s.length() == 3)
	    molecule_from_comp_id(s);
	 else
	    molecule_from_SMILES(s);
      }
      
      if (argc == 3) {
	 std::string comp_id = argv[1];
	 std::string file_name = argv[2];
	 molecule_from_ccd_pdbx(comp_id, file_name);
      }
   }
   return status;
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS

