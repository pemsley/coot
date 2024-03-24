
#include "lucrezia-tests.hh"

int test_multiligands_lig_bonding(molecules_container_t &mc) {
   starting_test(__FUNCTION__);
   int status = 0;
   int imol_paul_prot = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_lucr_prot = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-2.pdb"));
   mc.import_cif_dictionary(reference_data("00Z_for_prot1.cif"), imol_paul_prot);
   mc.import_cif_dictionary(reference_data("00Z_for_prot2.cif"), imol_lucr_prot);
   coot::protein_geometry &geom= mc.get_geometry();
   
   std::pair<bool, coot::dictionary_residue_restraints_t> lig1_pair = geom.get_monomer_restraints("00Z", imol_paul_prot);
   std::pair<bool, coot::dictionary_residue_restraints_t> lig2_pair = geom.get_monomer_restraints("00Z", imol_lucr_prot);
   if (lig1_pair.first) {
      if (lig2_pair.first) {
         std::cout << "found both dictionaries\n";
         // std::cout << "     " << lig1_pair.second << std::endl;
         // std::cout << "     " << lig2_pair.second << std::endl;
         int imol_lig1= mc.get_monomer_from_dictionary("00Z", imol_paul_prot, true);
         int imol_lig2= mc.get_monomer_from_dictionary("00Z", imol_lucr_prot, true);
         mc.write_coordinates(imol_lig1, "00Z_from_dict_for_paul_prot.pdb");
         mc.write_coordinates(imol_lig2, "00Z_from_dict_for_lucr_prot.pdb");
         mc.merge_molecules(imol_paul_prot, std::to_string(imol_lig1));
         mc.merge_molecules(imol_lucr_prot, std::to_string(imol_lig2));
         mc.write_coordinates(imol_paul_prot, "p_test.pdb");
         mc.write_coordinates(imol_lucr_prot, "l_test.pdb");
         bool check_hydrogens_too_flag = false;
         bool apply_bond_distance_check = false;
         coot::protein_geometry & geom = mc.get_geom();
         std::pair<bool, coot::dictionary_residue_restraints_t> pair = geom.get_monomer_restraints("00Z", imol_lucr_prot);
         if (pair.first) {
            coot::dictionary_residue_restraints_t restraints = pair.second; 
            // std::cout << "xxxxxxx" << restraints << std::endl;
            coot:: residue_spec_t res_spec("B", 1, "");
            mmdb::Residue *res = mc.get_residue(imol_lucr_prot, res_spec);
            if (res) {
                mmdb::Atom **atoms = nullptr;
                int natoms = 0;
                res->GetAtomTable(atoms, natoms);
                for(int iat = 0; iat < natoms; iat++) {
                    mmdb::Atom *at = atoms[iat];
                    // std::cout << "  " << iat << "  " << at->GetAtomName() << std::endl;
                } 
                std::pair<bool, std::vector<std::string> >
                matchers = geom.atoms_match_dictionary(res, check_hydrogens_too_flag, apply_bond_distance_check, restraints); 
                if (matchers.first) {
                   // std::cout << "found " << matchers.second.size() << " matchers" << std::endl;
                   if (matchers.second.size() == 15) status = 1;
                } else {
                   std::cout << "invalid match" << std::endl;
                }
            } else {
               std::cout << "failed to find " << res_spec << "in molecule" << imol_lucr_prot << std::endl;
            }
         }
      }
   }
   return status;
}
