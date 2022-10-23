
#include <iostream>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include "molecules_container.hh"

void starting_test(const char *func) {
   std::cout << "Starting test " << func << "()" << std::endl;
}

// wrap with the directory for the test data.
std::string
reference_data(const std::string &file) {
   return file; // 20221016-PE for now.
}

int test_auto_fit_rotamer(molecules_container_t &mc_in) {

   starting_test(__FUNCTION__);
   int status = 0; // initially fail status

   molecules_container_t mc;
   mc.geometry_init_standard();
   int imol = mc.read_pdb(reference_data("tm-A.pdb"));
   int imol_map = mc.read_mtz(reference_data("rnasa-1.8-all_refmac1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {
         coot::residue_spec_t res_spec("A", 89, "");
         mmdb::Residue *r = coot::util::get_residue(res_spec, mc[imol].atom_sel.mol);
         if (r) {
            mmdb::Atom *cz = r->GetAtom(" CZ ");
            coot::Cartesian pt_1(cz->x, cz->y, cz->z);
            status = mc.auto_fit_rotamer(imol, "A", 89, "", "", imol_map);
            coot::Cartesian pt_2(cz->x, cz->y, cz->z);
            double dd = coot::Cartesian::lengthsq(pt_1, pt_2);
            double d = std::sqrt(dd);
            std::cout << "d " << d << std::endl;
            if (d > 7.0)
               status = 1; // yay.
         } else {
            std::cout << "residue not found" << res_spec << std::endl;
         }
      } else {
         std::cout << "Non-valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "Non-valid model molecule " << imol << std::endl;
   }
   return status;
}

int test_pepflips(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   // A 14, 18, 20, 49
   //
   std::vector<coot::residue_spec_t> residues_for_flipping;
   std::vector<int> res_nos = {14, 18, 20, 49};
   for (const auto &rn : res_nos) {
      coot::residue_spec_t rs("A", rn, "");
      residues_for_flipping.push_back(rs);
   }

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));

   unsigned int n_flipped = 0;
   for (const auto &res_spec : residues_for_flipping) {
      coot::atom_spec_t atom_spec(res_spec.chain_id, res_spec.res_no, res_spec.ins_code, " O  ","");
      mmdb::Atom *at = mc.get_atom(imol, atom_spec);
      if (at) {
         coot::Cartesian pt_1(at->x, at->y, at->z);
         mc.flip_peptide(imol, atom_spec, "");
         coot::Cartesian pt_2(at->x, at->y, at->z);
         double dd = coot::Cartesian::lengthsq(pt_1, pt_2);
         double d = std::sqrt(dd);
         std::cout << "debug:: in test_pepflips() for " << atom_spec << " d is " << d << std::endl;
         if (d > 3.0)
            n_flipped++;
      } else {
         std::cout << "ERROR:: in test_pepflips() failed to find atom " << atom_spec << std::endl;
      }
   }

   // This flip A100 is OK to begin with so flipping it sould make the GruPoints
   // worse - and it does.

   if (n_flipped == res_nos.size()) {

      // test Atom cid and that flipping the N atom flips the previous residue
      std::string atom_cid = "//A/100/N";
      auto atom_spec = mc.atom_cid_to_atom_spec(imol, atom_cid);
      if (! atom_cid.empty()) {
         coot::atom_spec_t atom_spec_of_moving_O("A", 99, "", " O  ", "");
         mmdb::Atom *at = mc.get_atom(imol, atom_spec_of_moving_O);
         if (at) {
            coot::Cartesian pt_1(at->x, at->y, at->z);
            mc.flip_peptide_using_cid(imol, atom_cid, "");
            coot::Cartesian pt_2(at->x, at->y, at->z);
            double dd = coot::Cartesian::lengthsq(pt_1, pt_2);
            double d = std::sqrt(dd);
            if (d > 3.0) {
               status = 1;
            }
         }
      }
   }
   return status;
}

int test_updating_maps(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol =        mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   int imol_map      = mc.read_mtz(reference_data("gideondoesntapprove.mtz"), "FWT",    "PHWT",    "W", false, false);
   int imol_diff_map = mc.read_mtz(reference_data("gideondoesntapprove.mtz"), "DELFWT", "PHDELWT", "W", false, true);
   mc.associate_data_mtz_file_with_map(imol_map, reference_data("gideondoesntapprove.mtz"), "F", "SIGF", "FREER");

   mc.display_molecule_names_table();

   // set to the clipper map, overwriting the refmac map.
   //
   mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map);
   // After you have changed maps the firs time, add a starting point for the gru score:
   mc.calculate_new_gru_points(imol_diff_map);

   // modify the model by flipping a peptide.
   //
   std::string atom_cid = "//A/14/CA";
   auto as = mc.atom_cid_to_atom_spec(imol, atom_cid);
   if (! as.empty()) {
      coot::atom_spec_t atom_spec(as.chain_id, as.res_no, as.ins_code, " O  ","");
      mmdb::Atom *at = mc.get_atom(imol, atom_spec);
      if (at) {
         mc.flip_peptide_using_cid(imol, atom_cid, "");
      }
   }

   // now update the maps
   mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map);
   float new_gru_points = mc.calculate_new_gru_points(imol_diff_map);
   float gpt = mc.gru_points_total();
   std::cout << "###### GruPoints gained: " << new_gru_points << " gru points total " << gpt << std::endl;

   if (new_gru_points > 4.0)
      status = 1;

   return status;

}

int test_undo_and_redo(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   std::string atom_cid = "//A/14/CA";
   coot::atom_spec_t atom_spec("A", 14, "", " O  ", "");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   coot::Cartesian pt_1(at_1->x, at_1->y, at_1->z);
   mc.flip_peptide_using_cid(imol, atom_cid, "");
   coot::Cartesian pt_2(at_1->x, at_1->y, at_1->z);
   mc.undo(imol); // deletes atoms so now at_1 is out of date
   mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec);
   coot::Cartesian pt_3(at_2->x, at_2->y, at_2->z);

   double dd_1 = coot::Cartesian::lengthsq(pt_1, pt_2);
   double dd_2 = coot::Cartesian::lengthsq(pt_1, pt_3);
   double d_1 = std::sqrt(dd_1);
   double d_2 = std::sqrt(dd_2);

   std::cout << "test_undo(): debug distances " << d_1 << " " << d_2 << std::endl;
   if (d_1 > 3.0) {
      if (d_2 < 0.01) {

         // now let's test redo
         mc.redo(imol); // deletes atoms so now at_1 is out of date
         mmdb::Atom *at_3 = mc.get_atom(imol, atom_spec);
         coot::Cartesian pt_4(at_3->x, at_3->y, at_3->z);
         // modified and redone should be the same:
         double dd_3 = coot::Cartesian::lengthsq(pt_2, pt_4);
         double d_3 = std::sqrt(dd_3);
         std::cout << "test_undo(): debug distance d3 " << d_3 << std::endl;
         if (d_3 < 0.01)
            status = 1;
      }
   }
   return status;
}


int test_rama_balls_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string coords_fn = "tm-A.pdb";
   int imol = mc.read_pdb(coords_fn);

   coot::simple_mesh_t rvmm = mc.ramachandran_validation_markup_mesh(imol);
   std::cout << "debug:: rama mesh: " << rvmm.vertices.size() << " vertices and " << rvmm.triangles.size()
             << " triangles" << std::endl;

   // Let's look at the colours of the balls.
   if (false) // let's not.
      for (unsigned int i=0; i<rvmm.vertices.size(); i+=100)
         std::cout << i << " " << glm::to_string(rvmm.vertices[i].color) << std::endl;

   if (rvmm.vertices.size() > 2000) status = 1;
   return status;

}

int test_rota_dodecs_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   coot::simple_mesh_t rota_mesh = mc.get_rotamer_dodecs(imol);
   std::cout << "rota mesh: " << rota_mesh.vertices.size() << " vertices and " << rota_mesh.triangles.size()
             << " triangles" << std::endl;

   if (rota_mesh.vertices.size() > 2000)
      status = 1;

   return status;

}

int test_density_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   // this could be any mtz file I suppose
   int imol_map = mc.read_mtz("rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", false, false);

   clipper::Coord_orth p(55, 10, 10);
   float radius = 12;
   float contour_level = 0.13;
   coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
   std::cout << "DEBUG:: test_density_mesh(): " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
             << " triangles" << std::endl;

   if (map_mesh.vertices.size() > 30000)
      status = 1;

   return status;
}

int test_delete_atom(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   std::string atom_cid = "//A/14/O";
   mmdb::Atom *at_1 = mc.get_atom_using_cid(imol, atom_cid);
   if (at_1) {
      mc.delete_atom(imol, "A", 14, "", " O  ", "");
      mmdb::Atom *at_2 = mc.get_atom_using_cid(imol, atom_cid);
      if (at_2) {
         // bad, it was not deleted
      } else {
         mc.delete_using_cid(imol, "/1/B/200/CA", "ATOM");
         status = 1;
      }
   }
   return status;
}

int test_delete_residue(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   std::string residue_cid = "//A/14";
   mmdb::Residue *r_1 = mc.get_residue_using_cid(imol, residue_cid);
   if (r_1) {
      mc.delete_residue(imol, "A", 14, "");
      mmdb::Residue *r_2 = mc.get_residue_using_cid(imol, residue_cid);
      if (r_2) {
         // bad, it was not deleted
      } else {
         status = 1;
      }
   }
   return status;
}

coot::Cartesian atom_to_cartesian(mmdb::Atom *at) {

   return coot::Cartesian(at->x, at->y, at->z);

}

int test_rsr(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   int imol_map = mc.read_mtz("gideondoesntapprove.mtz", "FWT", "PHWT", "W", false, false);

   // int refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
   // const std::string &alt_conf, coot::molecule_t::refine_residues_mode mode);

   mc.set_imol_refinement_map(imol_map);

   std::string chain_id = "A";
   int res_no = 14; // this residue is problematic in gideondoesntapprove.pdb
   std::string ins_code;

   coot::atom_spec_t atom_spec_N_1(chain_id, res_no, ins_code, " N  ","");
   coot::atom_spec_t atom_spec_N_2(chain_id, res_no, ins_code, " CG2",""); // not a nitrogen atom
   mmdb::Atom *at_n_1 = mc.get_atom(imol, atom_spec_N_1);
   mmdb::Atom *at_n_2 = mc.get_atom(imol, atom_spec_N_2);
   coot::Cartesian pt_n_1_pre = atom_to_cartesian(at_n_1);
   coot::Cartesian pt_n_2_pre = atom_to_cartesian(at_n_2);
   
   std::string mode = "SPHERE";
   mc.refine_residues(imol, "A", 14, "", "", mode);
   coot::Cartesian pt_n_1_post = atom_to_cartesian(at_n_1);
   coot::Cartesian pt_n_2_post = atom_to_cartesian(at_n_2);
   mc.write_coordinates(imol, "refined-with-big-sphere.pdb");

   double dd_n_1 = coot::Cartesian::lengthsq(pt_n_1_pre, pt_n_1_post);
   double dd_n_2 = coot::Cartesian::lengthsq(pt_n_2_pre, pt_n_2_post);
   double d_1 = std::sqrt(dd_n_1);
   double d_2 = std::sqrt(dd_n_2);

   std::cout << "debug:: rsr distances " << d_1 << " " << d_2 << std::endl;

   if (d_1 > 0.1)
      if (d_2 > 0.1)
         status = true;

   return status;
}

int test_rsr_using_atom_cid(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   int imol_map = mc.read_mtz("gideondoesntapprove.mtz", "FWT", "PHWT", "W", false, false);

   mc.set_imol_refinement_map(imol_map);

   std::string chain_id = "A";
   int res_no = 14; // this residue is problematic in gideondoesntapprove.pdb
   std::string ins_code;

   std::string cid = "//A/14/CA";

   coot::atom_spec_t atom_spec_N_1(chain_id, res_no, ins_code, " N  ","");
   coot::atom_spec_t atom_spec_N_2(chain_id, res_no, ins_code, " CG2",""); // not a nitrogen atom
   mmdb::Atom *at_n_1 = mc.get_atom(imol, atom_spec_N_1);
   mmdb::Atom *at_n_2 = mc.get_atom(imol, atom_spec_N_2);
   coot::Cartesian pt_n_1_pre = atom_to_cartesian(at_n_1);
   coot::Cartesian pt_n_2_pre = atom_to_cartesian(at_n_2);
   
   std::string mode = "SPHERE";
   mc.refine_residues_using_atom_cid(imol, cid, mode);
   coot::Cartesian pt_n_1_post = atom_to_cartesian(at_n_1);
   coot::Cartesian pt_n_2_post = atom_to_cartesian(at_n_2);

   double dd_n_1 = coot::Cartesian::lengthsq(pt_n_1_pre, pt_n_1_post);
   double dd_n_2 = coot::Cartesian::lengthsq(pt_n_2_pre, pt_n_2_post);
   double d_1 = std::sqrt(dd_n_1);
   double d_2 = std::sqrt(dd_n_2);

   std::cout << "debug:: rsr distances " << d_1 << " " << d_2 << std::endl;

   if (d_1 > 0.1)
      if (d_2 > 0.1)
         status = true;

   return status;
}

int test_add_terminal_residue(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   int imol_map = mc.read_mtz("gideondoesntapprove.mtz", "FWT", "PHWT", "W", false, false);
   mc.set_imol_refinement_map(imol_map);

   // test adding to the N-terminus
   bool part_one_done = false;
   coot::atom_spec_t atom_spec_in_new_residue("A", 284, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec_in_new_residue);
   if (! at_1) { // it's not there to begin with
      mc.add_terminal_residue_directly(imol, "A", 285, "");
      mc.write_coordinates(imol, "with-added-terminal-residue.pdb");
      mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec_in_new_residue);
      if (at_2) {
         // now test that it is there
         coot::Cartesian reference_pos(76.3, 59, 23);
         coot::Cartesian atom_pos = atom_to_cartesian(at_2);
         double dd = coot::Cartesian::lengthsq(reference_pos, atom_pos);
         double d = std::sqrt(dd);
         std::cout << "debug d1 in test_add_terminal_residue " << d << std::endl;
         if (d < 1.0)
            part_one_done = true;
      } else {
         std::cout << "ERROR:: failed to find atom " << atom_spec_in_new_residue << std::endl;
      }
   } else {
      std::cout << "ERROR:: atom already exists " << atom_spec_in_new_residue << std::endl;
   }


   // test adding to the C-terminus
   bool part_two_done = false;
   atom_spec_in_new_residue = coot::atom_spec_t("A", 279, "", " O  ","");
   at_1 = mc.get_atom(imol, atom_spec_in_new_residue);
   if (! at_1) { // it's not there to begin with
      mc.add_terminal_residue_directly(imol, "A", 278, "");
      mc.write_coordinates(imol, "with-added-terminal-residue.pdb");
      mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec_in_new_residue);
      if (at_2) {
         // now test that it is there
         coot::Cartesian reference_pos(67.3, 60, 23);
         coot::Cartesian atom_pos = atom_to_cartesian(at_2);
         double dd = coot::Cartesian::lengthsq(reference_pos, atom_pos);
         double d = std::sqrt(dd);
         std::cout << "debug d2 in test_add_terminal_residue " << d << std::endl;
         if (d < 2.0)
            part_two_done = true;
      } else {
         std::cout << "ERROR:: failed to find atom " << atom_spec_in_new_residue << std::endl;
      }
   } else {
      std::cout << "ERROR:: atom already exists " << atom_spec_in_new_residue << std::endl;
   }

   if (part_one_done && part_two_done)
      status = 1;

   return status;
}

int test_delete_chain(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));

   coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      mc.delete_using_cid(imol, "/1/A/20/CA", "CHAIN");
      mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec);
      if (at_2) {
         // fail - it should not be there
      } else {
         status = 1;
      }
   } else {
      std::cout << "ERROR:: existing atom in imol not found " << atom_spec << std::endl;
   }
   return status;
}

int test_delete_molecule(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      mc.delete_using_cid(imol, "/1/A/20/CA", "MOLECULE");
      mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec);
      if (at_2) {
         // fail - it should not be there
      } else {
         status = 1;
      }
   } else {
      std::cout << "ERROR:: existing atom in imol not found " << atom_spec << std::endl;
   }
   return status;
}

int test_mutate(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   coot::atom_spec_t atom_spec_ser("A", 270, "", " OG ",""); // current SER
   coot::atom_spec_t atom_spec_tyr("A", 270, "", " OH ",""); // TYR
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec_ser);
   if (at_1) {
      // std::cout << "debug:: found atom at_1 " << at_1 << std::endl;
      std::string cid = "//A/270/CA";
      int mutate_status = mc.mutate(imol, cid, "TYR");
      mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec_tyr);
      if (at_2) {
         if (mutate_status == 1)
            status = 1;
      } else {
         std::cout << "failed to find at_2" << std::endl;
      }
   } else {
      std::cout << "failed to find at_1" << std::endl;
   }
   return status;
}


int test_weird_delete(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   coot::atom_spec_t atom_spec("A", 151, "", " N  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      mc.delete_using_cid(imol, "//A/151", "RESIDUE");
      coot::residue_spec_t res_spec("A", 151, "");
      mmdb::Residue *r = mc.get_residue(imol, res_spec);
      if (!r)
         status = 1;
   } else {
      std::cout << "ERROR:: failed to find test atom " << atom_spec << std::endl;
   }
   return status;
}

int test_side_chain_180(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   coot::atom_spec_t atom_spec("A", 268, "", " OD1","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
      std::string atom_cid = "//A/268/CA";
      mc.side_chain_180(imol, atom_cid);
      coot::Cartesian atom_pos_2 = atom_to_cartesian(at_1);
      double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
      double d = std::sqrt(dd);
      // std::cout << "debug test_side_chain_180 d " << d << std::endl;
      if (d > 2.0) // its 2.09
         status = 1.0;
   }
   return status;
}

int test_bonds_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));

   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   coot::simple_mesh_t mesh = mc.get_bonds_mesh(imol, mode);

   if (mesh.vertices.size() > 1000)
      if (mesh.triangles.size() > 1000)
         status = 1;

   std::cout << "INFO:: bonds mesh: " << mesh.vertices.size() << " vertices and " << mesh.triangles.size()
             << " triangles" << std::endl;

   return status;
}

int test_template(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      coot::Cartesian atom_pos = atom_to_cartesian(at_1);
   }
   return status;
}

int n_tests = 0;

int
run_test(int (*test_func) (molecules_container_t &mc), const std::string &test_name, molecules_container_t &mc) {

   n_tests++;
   int status = test_func(mc);
   std::string status_string = "FAIL: ";
   if (status == 1)
      status_string = "PASS: ";
   std::cout << status_string << std::setw(40) << std::left << test_name << " status " << status << std::endl;

   return status;
}

int main(int argc, char **argv) {

   int status = 0;

   bool last_test_only = false;
   if (argc > 1) {
      std::string arg(argv[1]);
      if (arg == "last-test-only")
         last_test_only = true;
   }

   molecules_container_t mc;

   mc.fill_rotamer_probability_tables();

   if (! last_test_only) {

      status += run_test(test_rama_balls_mesh,    "rama balls mesh",          mc);
      status += run_test(test_density_mesh,       "density mesh",             mc);
      status += run_test(test_pepflips,           "pepflips",                 mc);
      status += run_test(test_auto_fit_rotamer,   "auto-fit rotamer",         mc);
      status += run_test(test_updating_maps,      "updating maps",            mc);
      status += run_test(test_undo_and_redo,      "undo and redo",            mc);
      status += run_test(test_delete_residue,     "delete residue",           mc);
      status += run_test(test_delete_chain,       "delete chain",             mc);
      status += run_test(test_rota_dodecs_mesh,   "rotamer dodecahedra mesh", mc);
      status += run_test(test_rsr,                "rsr",                      mc);
      status += run_test(test_rsr_using_atom_cid, "rsr using atom cid",       mc);
      status += run_test(test_delete_molecule,    "delete_moelcule",          mc);
      status += run_test(test_add_terminal_residue, "add terminal residue",   mc);
      status += run_test(test_mutate,              "mutate",                  mc);
      status += run_test(test_delete_atom,        "delete atom",              mc);
      status += run_test(test_weird_delete,        "delete II",               mc);
      status += run_test(test_rsr,                "rsr",                      mc);
      status += run_test(test_side_chain_180,     "side-chain 180",           mc);
   }

   status += run_test(test_bonds_mesh,        "bonds mesh",      mc);


   int all_tests_status = 1; // fail!
   if (status == n_tests) all_tests_status = 0;

   return all_tests_status;

}
