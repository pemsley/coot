
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
   char *env = getenv("MOORHEN_TEST_DATA_DIR");
   if (env) {
      std::string joined = coot::util::append_dir_file(env, file);
      return joined;
   } else {
      return file;
   }
}

int test_auto_fit_rotamer(molecules_container_t &mc_in) {

   starting_test(__FUNCTION__);
   int status = 0; // initially fail status

   molecules_container_t mc;
   mc.geometry_init_standard();
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {
         coot::residue_spec_t res_spec("A", 61, "");
         mmdb::Residue *r = coot::util::get_residue(res_spec, mc[imol].atom_sel.mol);
         if (r) {
            mmdb::Atom *cz = r->GetAtom(" CZ ");
            coot::Cartesian pt_1(cz->x, cz->y, cz->z);
            status = mc.auto_fit_rotamer(imol, "A", 61, "", "", imol_map);
            coot::Cartesian pt_2(cz->x, cz->y, cz->z);
            double dd = coot::Cartesian::lengthsq(pt_1, pt_2);
            double d = std::sqrt(dd);
            std::cout << "d " << d << std::endl;
            if (d > 6.0)
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

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

   int imol =        mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map      = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT",    "PHWT",    "W", false, false);
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "DELFWT", "PHDELWT", "W", false, true);
   mc.associate_data_mtz_file_with_map(imol_map, reference_data("moorhen-tutorial-map-number-1.mtz"), "F", "SIGF", "FREER");

   mc.display_molecule_names_table();

   // set to the clipper map, overwriting the refmac map.
   //
   mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map, imol_map);
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
   mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map, imol_map);
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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
         if (d_3 < 0.01) {

            mc.flip_peptide_using_cid(imol, "//A/20/CA", "");
            mc.flip_peptide_using_cid(imol, "//A/22/CA", "");

            coot::atom_spec_t atom_spec_b("A", 24, "", " O  ", "");
            mmdb::Atom *at_5 = mc.get_atom(imol, atom_spec_b);
            coot::Cartesian pt_5(at_5->x, at_5->y, at_5->z);
            mc.flip_peptide_using_cid(imol, "//A/24/CA", "");
            mmdb::Atom *at_6 = mc.get_atom(imol, atom_spec_b);
            coot::Cartesian pt_6(at_6->x, at_6->y, at_6->z);
            mc.undo(imol);
            mmdb::Atom *at_7 = mc.get_atom(imol, atom_spec_b);
            coot::Cartesian pt_7(at_7->x, at_7->y, at_7->z);
            mc.redo(imol);
            mmdb::Atom *at_8 = mc.get_atom(imol, atom_spec_b);
            coot::Cartesian pt_8(at_8->x, at_8->y, at_8->z);

            if (true) { // debugging
               std::cout << "pt_4 " << pt_4 << std::endl;
               std::cout << "pt_5 " << pt_5 << std::endl;
               std::cout << "pt_6 " << pt_6 << std::endl;
               std::cout << "pt_7 " << pt_7 << std::endl;
               std::cout << "pt_8 " << pt_8 << std::endl;
            }

            double dd_4 = coot::Cartesian::lengthsq(pt_6, pt_7);
            double dd_5 = coot::Cartesian::lengthsq(pt_6, pt_8);
            double d_4 = std::sqrt(dd_4);
            double d_5 = std::sqrt(dd_5);

            std::cout << "debug d_4 " << d_4 << " dd_5 " << d_5 << std::endl;
            if (dd_4 > 3.0)
               if (dd_5 < 0.01)
                  status = 1;
         }
      }
   }
   return status;
}

int test_rama_validation(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string coords_fn = reference_data("moorhen-tutorial-structure-number-1.pdb");
   int imol = mc.read_pdb(coords_fn);

   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > rv = mc.ramachandran_validation(imol);

   bool r_285 = false;
   bool r_286 = false;
   for (const auto &r : rv) {
      // std::cout << " " << r.first << " " << r.second << std::endl;
      if (r.second.residue_number == 285) r_285 = true;
      if (r.second.residue_number == 286) r_286 = true;
   }

   if (r_286 && ! r_285)
      status = 1;

   return status;
}


int test_rama_balls_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string coords_fn = reference_data("moorhen-tutorial-structure-number-1.pdb");
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   // int refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
   // const std::string &alt_conf, coot::molecule_t::refine_residues_mode mode);

   mc.set_imol_refinement_map(imol_map);

   std::string chain_id = "A";
   int res_no = 14; // this residue is problematic in moorhen-tutorial-structure-number-1.pdb
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   mc.set_imol_refinement_map(imol_map);

   std::string chain_id = "A";
   int res_no = 14; // this residue is problematic in moorhen-tutorial-structure-number-1.pdb
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

int test_rsr_using_residue_range(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   mc.set_imol_refinement_map(imol_map);

   coot::atom_spec_t atom_spec_N_1("A", 130, "", " N  ","");
   coot::atom_spec_t atom_spec_N_2("A", 133, "", " N  ","");
   coot::atom_spec_t atom_spec_N_3("A", 137, "", " N  ","");
   mmdb::Atom *at_N_1 = mc.get_atom(imol, atom_spec_N_1);
   mmdb::Atom *at_N_2 = mc.get_atom(imol, atom_spec_N_2);
   mmdb::Atom *at_N_3 = mc.get_atom(imol, atom_spec_N_3);
   coot::Cartesian atom_pos_N_1_1 = atom_to_cartesian(at_N_1);
   coot::Cartesian atom_pos_N_2_1 = atom_to_cartesian(at_N_2);
   coot::Cartesian atom_pos_N_3_1 = atom_to_cartesian(at_N_3);
   float w = mc.get_map_weight();
   mc.set_map_weight(w * 100.0);
   mc.refine_residue_range(imol, "A", 131, 136);
   mc.set_map_weight(w); // restore sanity.
   coot::Cartesian atom_pos_N_1_2 = atom_to_cartesian(at_N_1);
   coot::Cartesian atom_pos_N_2_2 = atom_to_cartesian(at_N_2);
   coot::Cartesian atom_pos_N_3_2 = atom_to_cartesian(at_N_3);
   double dd_N_1 = coot::Cartesian::lengthsq(atom_pos_N_1_1, atom_pos_N_1_2);
   double dd_N_2 = coot::Cartesian::lengthsq(atom_pos_N_2_1, atom_pos_N_2_2);
   double dd_N_3 = coot::Cartesian::lengthsq(atom_pos_N_3_1, atom_pos_N_3_2);
   double d_N_1 =  std::sqrt(dd_N_1);
   double d_N_2 =  std::sqrt(dd_N_2);
   double d_N_3 =  std::sqrt(dd_N_3);
   // std::cout << "ds: " << d_N_1 << " " << d_N_2 << " " << d_N_3 << std::endl;
   if (d_N_1 < 0.0001)  // no move
      if (d_N_3 < 0.0001) // no move
         if (d_N_2 > 0.08) // move a bit
            status = 1;
   return status;
}

int test_add_terminal_residue(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
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
      // mc.write_coordinates(imol, "with-added-terminal-residue.pdb");
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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   coot::simple_mesh_t mesh = mc.get_bonds_mesh(imol, mode);

   if (mesh.vertices.size() > 1000)
      if (mesh.triangles.size() > 1000)
         status = 1;

   std::cout << "INFO:: bonds mesh: " << mesh.vertices.size() << " vertices and " << mesh.triangles.size()
             << " triangles" << std::endl;

   return status;
}

int test_copy_fragment_using_residue_range(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      int res_no_start = 131;
      int res_no_end   = 140;
      int imol_new = mc.copy_fragment_using_residue_range(imol, "A", res_no_start, res_no_end);
      if (mc.is_valid_model_molecule(imol_new)) {
         std::vector<mmdb::Residue *> residues;
         mmdb::Manager *mol = mc.get_mol(imol_new);
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     residues.push_back(residue_p);
                  }
               }
            }
         }
         if (residues.size() == 10) {
            coot::atom_spec_t atom_spec_2("A", 136, "", " O  ","");
            mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec_2);
            if (at_2)
               status = 1;
         }
      }
   }
   return status;
}

int test_copy_fragment_using_cid(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      std::string cid = "//A/131-140";
      int imol_new = mc.copy_fragment_using_cid(imol, cid);
      if (mc.is_valid_model_molecule(imol_new)) {
         std::vector<mmdb::Residue *> residues;
         mmdb::Manager *mol = mc.get_mol(imol_new);
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     residues.push_back(residue_p);
                  }
               }
            }
         }
         if (residues.size() == 10) {
            coot::atom_spec_t atom_spec_2("A", 136, "", " O  ","");
            mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec_2);
            if (at_2)
               status = 1;
         }
      }
   }
   return status;
}

int test_move_molecule_here(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::Cartesian pos(11,22,33);
   mc.move_molecule_to_new_centre(imol, pos.x(), pos.y(), pos.z());
   coot::Cartesian molecule_centre = mc.get_molecule_centre(imol);
   double dd = coot::Cartesian::lengthsq(pos, molecule_centre);
   double d = std::sqrt(dd);
   std::cout << "test_move_molecule_here d " << d << std::endl;
   if (d < 0.001)
      status = 1;
   return status;
}

int test_difference_map_contours(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "DELFWT", "PHDELWT", "", false, true);

   clipper::Coord_orth p(65, 50, 30);
   float radius = 10;
   float contour_level = 0.03;
   coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
   std::cout << "DEBUG:: test_density_mesh(): " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
             << " triangles" << std::endl;

   for (const auto &vertex : map_mesh.vertices) {
      // std::cout << "vertex " <<  glm::to_string(vertex.pos) << " " << glm::to_string(vertex.color) << std::endl;
      if (vertex.color[0] > vertex.color[1])
         if (vertex.color[0] > vertex.color[2])
            status = 1;
   }
   
   return status;
}

int test_jed_flip(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 225, "", " CZ ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
      std::string atom_cid = "//A/225/CB";
      bool invert_selection = true;
      mc.jed_flip(imol, atom_cid, invert_selection);
      coot::Cartesian atom_pos_2 = atom_to_cartesian(at_1);
      double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
      double d = std::sqrt(dd);
      std::cout << "test_jed_flip d " << d << std::endl;
      mc.write_coordinates(imol, "jed-flip.pdb");
      if (d > 0.9)
         status = true;
   }
   return status;
}


int test_sequence_generator(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::vector<std::string> chain_ids = mc.get_chains_in_model(imol);
   for (const auto &chain_id : chain_ids) {
      std::vector<std::pair<coot::residue_spec_t, std::string> > seq = mc.get_single_letter_codes_for_chain(imol, chain_id);
      std::string sp; // accumlate the sequence of the protein chain.
      unsigned int n_wat = 0;
      for (unsigned int i=0; i<seq.size(); i++) {
         const std::string &l = seq[i].second;
         if (l != "~")
            sp += l;
         else
            n_wat++;
      }
      std::cout << sp << std::endl;
      if (sp.length() > 50)
         if (n_wat == 150)
            status = 1;
   }
   return status;
}

int test_eigen_flip(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.get_monomer("ARG");
   coot::atom_spec_t atom_spec("A", 1, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      coot::Cartesian atom_pos_0 = atom_to_cartesian(at_1);
      mc.eigen_flip_ligand(imol, "A", 1, "");
      coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
      mc.eigen_flip_ligand(imol, "A", 1, "");
      coot::Cartesian atom_pos_2 = atom_to_cartesian(at_1);
      mc.eigen_flip_ligand(imol, "A", 1, "");
      coot::Cartesian atom_pos_3 = atom_to_cartesian(at_1);
      mc.eigen_flip_ligand(imol, "A", 1, "");
      coot::Cartesian atom_pos_4 = atom_to_cartesian(at_1);
      double dd_1 = coot::Cartesian::lengthsq(atom_pos_0, atom_pos_1);
      double dd_2 = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
      double dd_3 = coot::Cartesian::lengthsq(atom_pos_2, atom_pos_3);
      double dd_4 = coot::Cartesian::lengthsq(atom_pos_3, atom_pos_4);
      double dd_5 = coot::Cartesian::lengthsq(atom_pos_4, atom_pos_0);

      double d1 = std::sqrt(dd_1);
      double d2 = std::sqrt(dd_2);
      double d3 = std::sqrt(dd_3);
      double d4 = std::sqrt(dd_4);
      double d5 = std::sqrt(dd_5);

      std::cout << "test_eigen_flip ds " << d1 << " " << d2 << " " << d3 << " " << d4 << " " << d5 << std::endl;
      if (d1 > 3.0)
         if (d2 > 3.0)
            if (d3 > 3.0)
               if (d4 > 3.0)
                  if (d5 < 0.01)
                     status = true;
      mc.close_molecule(imol);
   }
   return status;
}

int test_non_standard_types(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::vector<std::string> nst = mc.get_residue_names_with_no_dictionary(imol);

   // weak test
   if (nst.empty())
      status = 1;

   return status;
}

int test_import_cif_dictionary(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   if (coot::file_exists("ATP.cif")) {

      mc.import_cif_dictionary("ATP.cif", coot::protein_geometry::IMOL_ENC_ANY);
      int imol = mc.get_monomer("ATP");
      if (mc.is_valid_model_molecule(imol)) {
         mc.close_molecule(imol);

         imol = mc.get_monomer_and_position_at("ATP", -999999, 5, 6, 7);
         std::cout << "debug:: in test_import_cif_dictionary() imol in centre test: " << imol << std::endl;
         coot::Cartesian ligand_centre = mc.get_molecule_centre(imol);
         coot::Cartesian expected(5,6,7);

         double dd = coot::Cartesian::lengthsq(ligand_centre, expected);
         double d = std::sqrt(dd);
         std::cout << "debug:: in test_import_cif_dictionary() ligand_centre is " << ligand_centre
                   << " d is " << d << std::endl;
         if (d < 0.001)
            status = 1;
      }

   } else {

      std::cout << "SKIP_TEST: test_import_cif_dictionary no ATP.cif in directory" << std::endl;
      status = 1; // can't test

   }

   return status;
}

int test_difference_map_peaks(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "DELFWT", "PHDELWT", "W", false, true);
   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_diff_map)) {
         float n_rmsd = 6.0; // lots of high peaks in the tutorial data. 6.0 is not a normal limit
         auto sites = mc.difference_map_peaks(imol_diff_map, imol, n_rmsd);
         if (sites.size() > 20)
            status = 1;
         if (false)
            for (const auto &site : sites)
               std::cout << "site " << site.feature_type << " " << site.button_label << " "
                         << site.x << " " << site.y << " " << site.z << " residue " << site.residue_spec
                         << " height " << site.feature_value << " badness " << site.badness << std::endl;
      }
   }
   return status;
}

int test_pepflips_using_difference_map(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "DELFWT", "PHDELWT", "W", false, true);
   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_diff_map)) {
         float n_rmsd = 4.0;
         std::vector<coot::molecule_t::interesting_place_t> flips = mc.pepflips_using_difference_map(imol, imol_diff_map, n_rmsd);
         if (flips.size() > 2)
            status = 1;
         if (true)
            for (const auto &flip : flips)
               std::cout << "flip: " << flip.feature_type << " " << flip.button_label << " " << flip.x << " " << flip.y << " " << flip.z
                         << " badness " << flip.badness << std::endl;
         mc.close_molecule(imol);
         mc.close_molecule(imol_diff_map);
      }
   }
   return status;
}

int test_template(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      coot::Cartesian atom_pos = atom_to_cartesian(at_1);
      double dd = coot::Cartesian::lengthsq(atom_pos, atom_pos);
      double d = std::sqrt(dd);
      std::cout << "test_ d " << d << std::endl;
   }
   mc.close_molecule(imol);
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
      status += run_test(test_difference_map_contours, "difference map density mesh", mc);
      status += run_test(test_pepflips,           "pepflips",                 mc);
      status += run_test(test_auto_fit_rotamer,   "auto-fit rotamer",         mc);
      status += run_test(test_updating_maps,      "updating maps",            mc);
      status += run_test(test_delete_residue,     "delete residue",           mc);
      status += run_test(test_delete_chain,       "delete chain",             mc);
      status += run_test(test_rota_dodecs_mesh,   "rotamer dodecahedra mesh", mc);
      status += run_test(test_rsr,                "rsr",                      mc);
      status += run_test(test_rsr_using_atom_cid, "rsr using atom cid",       mc);
      status += run_test(test_rsr_using_residue_range, "rsr using residue range", mc);
      status += run_test(test_delete_molecule,    "delete_moelcule",          mc);
      status += run_test(test_mutate,              "mutate",                  mc);
      status += run_test(test_delete_atom,        "delete atom",              mc);
      status += run_test(test_weird_delete,       "delete II",                mc);
      status += run_test(test_side_chain_180,     "side-chain 180",           mc);
      status += run_test(test_bonds_mesh,         "bonds mesh",               mc);
      status += run_test(test_undo_and_redo,      "undo and redo",            mc);
      status += run_test(test_add_terminal_residue, "add terminal residue",   mc);
      status += run_test(test_rama_validation,    "rama validation",          mc);
      status += run_test(test_copy_fragment_using_residue_range, "copy-fragment using residue range", mc);
      status += run_test(test_copy_fragment_using_cid, "copy-fragment using cid", mc);
      status += run_test(test_move_molecule_here,    "move_molecule_here",    mc);
      status += run_test(test_jed_flip,             "JED Flip",               mc);
      status += run_test(test_sequence_generator, "Make a sequence string",   mc);
      status += run_test(test_non_standard_types, "non-standard residue types in molecule",   mc);
      status += run_test(test_pepflips_using_difference_map,         "Pepflips from Difference Map",               mc);
      status += run_test(test_difference_map_peaks, "Difference Map Peaks",   mc);
      status += run_test(test_eigen_flip,         "Eigen Flip",               mc);
   }



      status += run_test(test_import_cif_dictionary, "import cif dictionary", mc);

      // change the autofit_rotamer test so that it tests the change of positions of the atoms of the neighboring residues.

   int all_tests_status = 1; // fail!
   if (status == n_tests) all_tests_status = 0;

   return all_tests_status;

}
