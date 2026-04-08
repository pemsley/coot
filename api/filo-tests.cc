/*
* @author Filomeno Sanchez
* 
* Copyright 2023, University of York
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation; either version 3 of the License, or (at
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

#include "filo-tests.hh"

int test_something_filo(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::cout << "BEFORE " << mc.map_sampling_rate << "\n";

   mc.set_map_sampling_rate(1.8);

   std::cout << "AFTER " << mc.map_sampling_rate << "\n";

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         coot::Cartesian atom_pos = atom_to_cartesian(at_1);
         double dd = coot::Cartesian::lengthsq(atom_pos, atom_pos);
         double d = std::sqrt(dd);
         std::cout << "test_ d " << d << std::endl;

         coot::validation_information_t dca = mc.density_correlation_analysis(imol, imol_map);
         for (const auto &chain : dca.cviv) {
            for (const auto &res : chain.rviv) {
               if (res.residue_spec.res_no == 62) {
                  std::cout << "function value " << res.function_value << std::endl;
                  if (res.function_value > 0.6) {
                     status = 1;
                  }
               }
            }
         }
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   return status;
}


int test_get_diff_map_peaks(molecules_container_t &mc) {

   // this test needs a different mtz file. Come back later.

   starting_test(__FUNCTION__);
   int status = 0;

   int coordMolNo   = mc.read_pdb(reference_data("5a3h.pdb"));
   int mapMolNo     = mc.read_mtz(reference_data("5a3h_sigmaa.mtz"), "FWT",    "PHWT",    "FOM", false, false);
   int diffMapMolNo = mc.read_mtz(reference_data("5a3h_sigmaa.mtz"), "DELFWT", "PHDELWT", "FOM", false, true);

   mc.associate_data_mtz_file_with_map(mapMolNo, reference_data("5a3h_sigmaa.mtz"), "FP", "SIGFP", "FREE");

   mc.connect_updating_maps(coordMolNo, mapMolNo, mapMolNo, diffMapMolNo);
   // if sfcalc_genmaps_using_bulk_solvent() fails stats.r_factor is -1

   coot::util::sfcalc_genmap_stats_t stats = mc.sfcalc_genmaps_using_bulk_solvent(coordMolNo, mapMolNo, diffMapMolNo, mapMolNo);
   std::cout << "   stats 0 : r_factor " << stats.r_factor << std::endl;

   molecules_container_t::r_factor_stats rfs_1 = mc.get_r_factor_stats();
   std::cout << "   stats 1 : r_factor " << rfs_1.r_factor << std::endl;
   mc.get_map_contours_mesh(mapMolNo,  77.501,  45.049,  22.663,  13,  0.48);

   mc.delete_using_cid(coordMolNo, "/1/A/300/*", "LITERAL");
   mc.get_map_contours_mesh(mapMolNo,  77.501,  45.049,  22.663,  13,  0.48);
   molecules_container_t::r_factor_stats rfs_2 = mc.get_r_factor_stats();
   std::cout << "   stats 2 : r_factor " << rfs_2.r_factor << std::endl;

   auto diff_diff_map_peaks = mc.get_diff_diff_map_peaks(diffMapMolNo,  77.501,  45.049,  22.663);

   if (diff_diff_map_peaks.size() >  0) {
      status = 1;
   }

   mc.close_molecule(coordMolNo);
   mc.close_molecule(mapMolNo);
   mc.close_molecule(diffMapMolNo);

   return status;

}

int test_non_drawn_bond_multi_cid_2(molecules_container_t &mc) {

   // test.only("Test non-drawn bonds and multi CID selection mesh --second", () => {
   starting_test(__FUNCTION__);
   int status = 0;

   mc.set_use_gemmi(false);
   int coordMolNo = mc.read_pdb(reference_data("5a3h.pdb"));
   if (mc.is_valid_model_molecule(coordMolNo)) {
      auto instanceMesh_1 = mc.get_bonds_mesh_for_selection_instanced(
         coordMolNo, "//A/10-20||//A/25-30", "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, false, false, false, true, 1);

      mc.add_to_non_drawn_bonds(coordMolNo, "//A/26-30");
      auto instanceMesh_2 = mc.get_bonds_mesh_for_selection_instanced(
         coordMolNo, "//A/10-20||//A/25-30", "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, false, false, false, true, 1);

      // mc.print_non_drawn_bonds(coordMolNo);

      // std::cout << "debug:: mesh sizes " << instanceMesh_1.geom[1].instancing_data_B.size() << " "
      //                                    << instanceMesh_2.geom[1].instancing_data_B.size()<< std::endl;

      std::cout << "debug:: A sizes "
                << instanceMesh_1.geom[0].instancing_data_A.size() << " "
                << instanceMesh_2.geom[0].instancing_data_A.size()<< std::endl;
      std::cout << "debug:: B sizes "
                << instanceMesh_1.geom[1].instancing_data_B.size() << " "
                << instanceMesh_2.geom[1].instancing_data_B.size()<< std::endl;
      if (instanceMesh_2.geom[1].instancing_data_B.size() < instanceMesh_1.geom[1].instancing_data_B.size()) {
         if (instanceMesh_2.geom[0].instancing_data_A.size() < instanceMesh_1.geom[0].instancing_data_A.size()) {
            status = 1;
         } else {
            std::cout << "bad A mesh size match" << std::endl;
         }
      } else {
         std::cout << "bad B mesh size match" << std::endl;
      }
   } else {
      std::cout << "Bad read for 5a3h.pdb" << std::endl;
   }

   return status;

}

int test_change_chain_id_1(molecules_container_t &molecules_container) {

   // test.only("Test change chain ID --first", () => {

   starting_test(__FUNCTION__);
   int status = 0;

   molecules_container.set_use_gemmi(false);
   int coordMolNo = molecules_container.read_pdb(reference_data("5a3h.pdb"));
   molecules_container.delete_colour_rules(coordMolNo);

   // let colourMap = new cootModule.MapIntFloat3();
   // let indexedResiduesVec = new cootModule.VectorStringUInt_pair();

   std::map<unsigned int, std::array<float, 4> > colourMap;
   std::vector<std::pair<std::string, unsigned int> > indexedResiduesVec;

   std::vector<std::pair<std::string, std::array<float, 4> > > colours = {
      std::make_pair("//A", std::array<float,4>({1, 0, 0, 1})),
      std::make_pair("//B", std::array<float,4>({0, 0, 1, 1}))};

   // colours.forEach((colour, index) => {
   //    colourMap.set(index + 51, colour.rgb)
   //    const i = { first: colour.cid, second: index + 51 }
   //    indexedResiduesVec.push_back(i)
   // });

   for(unsigned int i=0; i<colours.size(); i++) {
         colourMap[i+51] = colours[i].second;
         indexedResiduesVec.push_back(std::make_pair(colours[i].first, i+51));
   }

   molecules_container.set_user_defined_bond_colours(coordMolNo, colourMap);
   molecules_container.set_user_defined_atom_colour_by_selection(coordMolNo, indexedResiduesVec, false);
   molecules_container.add_colour_rule(coordMolNo, "//A", "#ff0000");
   molecules_container.add_colour_rule(coordMolNo, "//B", "#0000ff");

   auto instanceMesh_1 = molecules_container.get_bonds_mesh_instanced(
      coordMolNo, "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, false, false, false, true, 1);

   std::vector<std::string> original_chains;
   auto original_chains_vec = molecules_container.get_chains_in_model(coordMolNo);
   unsigned int original_chains_vec_size = original_chains_vec.size();
   for (unsigned int i = 0; i < original_chains_vec_size; i++) {
      const auto &chain_name = original_chains_vec.at(i);
      original_chains.push_back(chain_name);
   }

   molecules_container.change_chain_id(coordMolNo, "A", "X", false, 0, 0);

   std::vector<std::string> new_chains;
   auto new_chains_vec = molecules_container.get_chains_in_model(coordMolNo);
   int new_chains_vec_size = new_chains_vec.size();
   for (int i = 0; i < new_chains_vec_size; i++) {
      const auto &chain_name = new_chains_vec.at(i);
      new_chains.push_back(chain_name);
   }

   // expect(new_chains).not.toEqual(original_chains)
   // expect(original_chains.includes('X')).toBeFalsy()
   // expect(new_chains.includes('A')).toBeFalsy()
   // expect(new_chains.includes('X')).toBeTruthy()

   // ----> IF HERE I DELETE THE COLOUR RULES AND I DEFINE THEM AGAIN WITHOUT SPECIFYING A COLOUR FOR CHAIN X I GET BUS ERROR
   // BUT IF I DON'T DELETE THEM THEN EVERYHTING IS FINE ?????
   molecules_container.delete_colour_rules(coordMolNo);

   // let colourMap_2 = new cootModule.MapIntFloat3()
   // let indexedResiduesVec_2 = new cootModule.VectorStringUInt_pair()

   std::map<unsigned int, std::array<float, 4> > colourMap_2;
   std::vector<std::pair<std::string, unsigned int> > indexedResiduesVec_2;

   // const colours_2 = [
   //    { cid: "//A", rgb: [1, 0, 0] },
   //    { cid: "//B", rgb: [0, 0, 1] },
   //    //{ cid: "//X", rgb: [0, 1, 0] }   ----> IF I COMMENT THIS OUT I GET BUS ERROR
   // ]

   std::vector<std::pair<std::string, std::array<float, 4> > > colours_2 = {
      std::make_pair("//A", std::array<float,4>({1, 0, 0, 1})),
      std::make_pair("//B", std::array<float,4>({0, 0, 1, 1}))
   };

   // colours_2.forEach((colour, index) => {
   //    colourMap_2.set(index + 51, colour.rgb)
   //    const i = { first: colour.cid, second: index + 51 }
   //    indexedResiduesVec_2.push_back(i)
   // });

   for(unsigned int i=0; i<colours_2.size(); i++) {
         colourMap[i+51] = colours_2[i].second;
         indexedResiduesVec.push_back(std::make_pair(colours_2[i].first, i+51));
   }

   molecules_container.set_user_defined_bond_colours(coordMolNo, colourMap_2);
   molecules_container.set_user_defined_atom_colour_by_selection(coordMolNo, indexedResiduesVec_2, false);
   molecules_container.add_colour_rule(coordMolNo, "//A", "#ff0000");
   molecules_container.add_colour_rule(coordMolNo, "//B", "#0000ff");
   //molecules_container.add_colour_rule(coordMolNo, '//X', '#0000ff').  ----> IF I COMMENT THIS OUT I GET BUS ERROR

   auto instanceMesh_2 = molecules_container.get_bonds_mesh_instanced(
      coordMolNo, "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, false, false, false, true, 1);

   // expect(
   //    instanceMesh_2.geom.get(1).instancing_data_B.size()
   // ).toBe(
   //    instanceMesh_1.geom.get(1).instancing_data_B.size()
   // )

   if (instanceMesh_2.geom.size() > 1)
      if (instanceMesh_2.geom.at(1).instancing_data_B.size() == instanceMesh_1.geom.at(1).instancing_data_B.size())
         status = 1;

   return status;
}


int test_change_rotamer(molecules_container_t &molecules_container) {


   starting_test(__FUNCTION__);

   int status = 0;
   int imol_molecule = molecules_container.read_pdb(reference_data("5a3h.pdb"));

   // Create a fragment and change rotamer
   int imol_fragment = molecules_container.copy_fragment_using_cid(imol_molecule, "//A/179");

   auto change_info = molecules_container.change_to_next_rotamer(imol_fragment, "//A/179", "");

   std::cout << "change_info: rank " << change_info.rank << std::endl;
   std::cout << "change_info: name " << change_info.name << std::endl;
   std::cout << "change_info: richardson_probability " << change_info.richardson_probability << std::endl;
   std::cout << "change_info: status " << change_info.status << std::endl;

   // Get the OG atom for that new rotamer (still in the fragment)
   auto resSpec = coot::residue_spec_t("A", 179, "");
   mmdb::Residue *res_fragment = molecules_container.get_residue(imol_fragment, resSpec);
   mmdb::Atom *atom_fragment = res_fragment->GetAtom(5);

   // Replace fragment back into the molecule and get new OG atom
   int status_replace = molecules_container.replace_fragment(imol_molecule, imol_fragment, "//A/179");
   std::cout << "replace_fragment() status " << status_replace << std::endl;
   mmdb::Residue *res_new = molecules_container.get_residue(imol_molecule, resSpec);
   mmdb::Atom *atom_new = res_new->GetAtom(5);

   std::cout << "atom_fragment pos: " << atom_fragment->x << " " << atom_fragment->y << " " << atom_fragment->z << std::endl;
   std::cout << "atom_new pos:      " << atom_new->x      << " " << atom_new->y      << " " << atom_new->z      << std::endl;

   // This fails...
   // expect(atom_new.x).toBe(atom_fragment.x)
   // expect(atom_new.y).toBe(atom_fragment.y)
   // expect(atom_new.z).toBe(atom_fragment.z)

   clipper::Coord_orth co_1 = coot::co(atom_fragment);
   clipper::Coord_orth co_2 = coot::co(atom_new);
   double d2 = (co_2-co_1).lengthsq();
   double d = std::sqrt(d2);
   if (d < 0.0001) status = 1;

   return status;
 }


// test.skip("Test import ligands with same name and animated refinement", () => {
int test_import_ligands_with_same_name_and_animated_refinement(molecules_container_t &molecules_container) {

   int status = 0;
   int coordMolNo_1 = molecules_container.read_pdb(reference_data("./5a3h.pdb"));
   int coordMolNo_2 = molecules_container.read_pdb(reference_data("./5fjj.pdb"));
   int mapMolNo = molecules_container.read_mtz(reference_data("./5a3h_sigmaa.mtz"), "FWT", "PHWT", "", false, false);

   molecules_container.import_cif_dictionary(     reference_data("./benzene.cif"), coordMolNo_1);
   molecules_container.import_cif_dictionary(reference_data("./nitrobenzene.cif"), coordMolNo_2);

   const std::string tlc = "LIG";

   int ligandMolNo_1 = molecules_container.get_monomer_and_position_at(tlc, coordMolNo_1, 0,0,0);
   auto merge_info_1 = molecules_container.merge_molecules(coordMolNo_1, std::to_string(ligandMolNo_1));
   std::cout << "merge_info_1.second.size() " << merge_info_1.second.size() << std::endl;
   if (merge_info_1.second.size() == 1) {

      int ligandMolNo_2 = molecules_container.get_monomer_and_position_at(tlc, coordMolNo_2, 0,0,0);
      auto merge_info_2 = molecules_container.merge_molecules(coordMolNo_2, std::to_string(ligandMolNo_2));
      // the first indicates that the merge actually happened.
      std::cout << "merge_info_2.first " << merge_info_2.first << std::endl;
      std::cout << "merge_info_2.second.size() " << merge_info_2.second.size() << std::endl;
      if (merge_info_2.second.size() == 1) {

         int copyMolNo_1 = molecules_container.copy_fragment_for_refinement_using_cid(coordMolNo_1, "/1/C/1/*");
         molecules_container.init_refinement_of_molecule_as_fragment_based_on_reference(copyMolNo_1, coordMolNo_1, mapMolNo);
         bool copy_status_1 = molecules_container.copy_dictionary("LIG", coordMolNo_1, copyMolNo_1);
         std::cout << "debug:: copy_dictionary() copy_status_1: " << copy_status_1 << std::endl;
         // let result_1 = [];
         molecules_container.display_molecule_names_table();
         std::cout << "debug:: copyMolNo_1 is " << copyMolNo_1 << std::endl;
         std::vector<glm::vec3> result_1;
         auto refine_result_1 = molecules_container.refine(copyMolNo_1, 5000); // returns a mesh in a pair
         const coot::instanced_mesh_t &instanced_mesh_1 = refine_result_1.second;
         const auto &geom_vec_1 = instanced_mesh_1.geom;
         unsigned int geom_vec_1_size = geom_vec_1.size();
         for (unsigned int i = 0; i < geom_vec_1_size; i++) {
            const auto &geom = geom_vec_1.at(i);
            const auto &inst_data_B_vec = geom.instancing_data_B;
            unsigned int inst_data_B_vec_size = inst_data_B_vec.size();
            for (unsigned int j = 0; j < inst_data_B_vec_size; j++) {
               const auto &inst_data_B = inst_data_B_vec.at(j);
               result_1.push_back(inst_data_B.size);
            }
         }
         molecules_container.clear_refinement(coordMolNo_1);

         int copyMolNo_2 = molecules_container.copy_fragment_for_refinement_using_cid(coordMolNo_2, "/1/j/1/*");
         molecules_container.init_refinement_of_molecule_as_fragment_based_on_reference(copyMolNo_2, coordMolNo_2, mapMolNo);
         bool copy_status_2 = molecules_container.copy_dictionary("LIG", coordMolNo_2, copyMolNo_2);
         std::cout << "debug:: copy_dictionary() copy_status_2: " << copy_status_2 << std::endl;
         // let result_2 = []
         std::vector<glm::vec3> result_2;
         auto refine_result_2 = molecules_container.refine(copyMolNo_2, 5000);
         const auto &instanced_mesh_2 = refine_result_2.second;
         const auto &geom_vec_2 = instanced_mesh_2.geom;
         unsigned int geom_vec_2_size = geom_vec_2.size();
         for (unsigned int i = 0; i < geom_vec_2_size; i++) {
            const auto &geom = geom_vec_2.at(i);
            const auto &inst_data_B_vec = geom.instancing_data_B;
            unsigned inst_data_B_vec_size = inst_data_B_vec.size();
            for (unsigned int j = 0; j < inst_data_B_vec_size; j++) {
               const auto &inst_data_B = inst_data_B_vec.at(j);
               result_2.push_back(inst_data_B.size);
            }
         }
         molecules_container.clear_refinement(coordMolNo_2);

         if (result_1.size() == 15) {
            if (result_2.size() == 22) {

               status = 1;

               // expect(result_1.every(size => size <= 2.25)).toBeTruthy()
               // expect(result_2.every(size => size <= 2.25)).toBeTruthy()

               for (const auto &r1 : result_1) {
                  std::cout << "r1.z " << r1.z << std::endl;
                  if (r1.z > 2.25) { status = 0; break; }
               }
               for (const auto &r2 : result_2) {
                  std::cout << "r2.z " << r2.z << std::endl;
                  if (r2.z > 2.25) { status = 0; break; }
               }
            }
         }
      }
   }
   return status;
}
   
