
#ifndef MAKE_ENTERPRISE_TOOLS
int main(int argc, char **argv) {return 0;}
#else 
#include "cod-types.hh"
#include "rdkit-interface.hh"
#include "coot-coord-utils.hh"

int main(int argc, char **argv) {

   int status = 0;

   if (argc > 1) {
      std::string comp_id = argv[1];

      try {
	 CResidue *residue_p = 0;
	 coot::protein_geometry geom;
	 // geom.init_standard();

	 // std::string three_letter_code = "001"; // nice test
	 // std::string three_letter_code = "0PY"; // simple
	 // std::string three_letter_code = "06C"; // simplest
	 
	 bool idealised_flag = true;
	 CMMDBManager *mol = geom.mol_from_dictionary(comp_id, idealised_flag);

	 if (! mol) {
	    std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
	 } else {

	    CResidue *residue_p = coot::util::get_first_residue(mol);

	    if (! residue_p) {
	       // pretty strange
	       std::cout << "Null residue from mol from mol_from_dictionary() for "
			 << comp_id << std::endl;
	    } else { 

	       std::string s = "COc1ccc(cc1O[C@H]1C[C@@H]2CC[C@H]1C2)C1CNC(=O)NC1";
	       // s = "CC12CC[C@@]3(C1)CC(=C)CC3C(C2C(C)(C)C)C(C)(C)C"; // bird
	       // RDKit::RWMol *rdkmp = RDKit::SmilesToMol(s);
	       // RDKit::RWMol rdkm(*rdkmp);
	       
	       RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, geom);
	       
	       for (unsigned int iat=0; iat<rdkm.getNumAtoms(); iat++) {
		  try {
		     std::string name;
		     RDKit::ATOM_SPTR at_p = rdkm[iat];
		     at_p->getProp("name", name);
		     std::cout << iat << "   \"" << name << "\"" << std::endl;
		  }
		  catch (KeyErrorException &err) {
		     std::cout << "caught no-name exception in rdkit_mol H-block" << std::endl;
		  }
	       }
	       
	       std::vector<std::string> v = cod::get_cod_atom_types(rdkm);
	       std::cout << "PE-TYPES:: -------- got " << v.size() << " atoms " << std::endl;
	       for (unsigned int i=0; i<v.size(); i++) { 
		  std::cout << "   " << i << " " << v[i] << "" << std::endl;
	       }
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "test-cod-atom-types caught exception: " << rte.what() << std::endl;
      }
   } 
   return status;
}

#endif // MAKE_ENTERPRISE_TOOLS

