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
         coordMolNo, "//A/10-20||//A/25-30", "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, 1);

      mc.add_to_non_drawn_bonds(coordMolNo, "//A/26-30");
      auto instanceMesh_2 = mc.get_bonds_mesh_for_selection_instanced(
         coordMolNo, "//A/10-20||//A/25-30", "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, 1);

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

   std::map<unsigned int, std::array<float, 3> > colourMap;
   std::vector<std::pair<std::string, unsigned int> > indexedResiduesVec;

   std::vector<std::pair<std::string, std::array<float, 3> > > colours = {
      std::make_pair("//A", std::array<float,3>({1, 0, 0})),
      std::make_pair("//B", std::array<float,3>({0, 0, 1}))};

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
      coordMolNo, "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, 1);

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

   std::map<unsigned int, std::array<float, 3> > colourMap_2;
   std::vector<std::pair<std::string, unsigned int> > indexedResiduesVec_2;

   // const colours_2 = [
   //    { cid: "//A", rgb: [1, 0, 0] },
   //    { cid: "//B", rgb: [0, 0, 1] },
   //    //{ cid: "//X", rgb: [0, 1, 0] }   ----> IF I COMMENT THIS OUT I GET BUS ERROR
   // ]

   std::vector<std::pair<std::string, std::array<float, 3> > > colours_2 = {
      std::make_pair("//A", std::array<float,3>({1, 0, 0})),
      std::make_pair("//B", std::array<float,3>({0, 0, 1}))
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
      coordMolNo, "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1, 1, 1);

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
