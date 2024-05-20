
#include "lucrezia-tests.hh"

int test_import_ligands_with_same_name(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   // Fragment185/dict.cif is not in the test data, so this test
   // can't/won't properly pass

   int coord_mol_no1 = mc.read_pdb(reference_data("5a3h.pdb"));
   int coord_mol_no2 = mc.read_pdb(reference_data("5fjj.pdb"));
   int map_mol_no1 = mc.read_mtz(reference_data("5a3h_sigmaa.mtz"),"FWT","PHWT","", false, false);
   mc.import_cif_dictionary(reference_data("Fragment185/dict.cif"), coord_mol_no1);
   std::string tlc = "LIG";
   int lig_mol_no1 = mc.get_monomer_and_position_at(tlc, coord_mol_no1, 0, 0, 0);
   auto merge_info_1 = mc.merge_molecules(coord_mol_no1, std::to_string(lig_mol_no1));

   int copy_mol_no1 = mc.copy_fragment_for_refinement_using_cid(coord_mol_no1, "/1/C/1/*");
   mc.copy_dictionary("LIG", coord_mol_no1, copy_mol_no1);
   mc.init_refinement_of_molecule_as_fragment_based_on_reference(copy_mol_no1, coord_mol_no1, map_mol_no1);
   auto refine_results_1 = mc.refine(copy_mol_no1, 30);
   // mc.write_coordinates(copy_mol_no1, "copy_mol_no1-refine-30.pdb");
   auto &instanced_mesh_1 = refine_results_1.second;
   auto &geom_vec_1 = instanced_mesh_1.geom;
   int geom_vec_1_size = geom_vec_1.size();

   std::vector <glm::vec3> results_1;
   for (int i = 0; i < geom_vec_1_size; i++) {
      auto &geom = geom_vec_1[i];
      auto &inst_data_B_vec = geom.instancing_data_B;
      int inst_data_B_vec_size =  inst_data_B_vec.size();
      for (int j = 0; j < inst_data_B_vec_size; j++) {
         auto &inst_data_B = inst_data_B_vec[j];
         results_1.push_back(inst_data_B.size);
         // check that all the bonds are less than 2.2A
      }
   }
   mc.clear_refinement(coord_mol_no1);

   std::cout << "results_1 size: " << results_1.size() << std::endl;

   for (const auto &bond_size : results_1)
      std::cout << "results_1 " << bond_size.x <<  " " << bond_size.y << " " << bond_size.z << std::endl;

   status = 1; 
   return status;
}

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
