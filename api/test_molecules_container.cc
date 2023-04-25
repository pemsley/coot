
#include <iostream>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include "molecules_container.hh"

void starting_test(const char *func) {
   std::cout << "Starting " << func << "()" << std::endl;
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

int test_utils(molecules_container_t &mc_in) {

   int status = 0; // initially fail status
   starting_test(__FUNCTION__);

   std::string test_string = "aaaaa||bbbb||c";
   std::vector<std::string> parts = coot::util::split_string(test_string, "||");
   if (parts.size() == 3) {
      if (parts[0] == "aaaaa") {
         if (parts[1] == "bbbb") {
            if (parts[2] == "c") {
               status = 1;
            }
         }
      }
   }

   return status;
}

int test_auto_fit_rotamer_1(molecules_container_t &mc_in) {

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
            if (d > 6.0) {
               status = 1;
            } else {
               std::cout << "bad d " << d << std::endl;
            }
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

int test_auto_fit_rotamer_2(molecules_container_t &mc_in) {

   starting_test(__FUNCTION__);
   int status = 0; // initially fail status

   molecules_container_t mc;
   mc.geometry_init_standard();
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {

         mc.mutate(imol, "//A/62/CA", "ARG");
         mc.write_coordinates(imol, "test-mc-post-mutate.pdb");
         mc.auto_fit_rotamer(imol, "A", 62, "", "", imol_map);
         mc.write_coordinates(imol, "test-mc-post-auto-fit-rotamer.pdb");
         // did it fit?
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

   int imol          = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map      = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT",    "PHWT",    "W", false, false);
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "DELFWT", "PHDELWT", "W", false, true);
   mc.associate_data_mtz_file_with_map(imol_map, reference_data("moorhen-tutorial-map-number-1.mtz"), "F", "SIGF", "FREER");

   // debugging
   mc.display_molecule_names_table();

   // set to the clipper map, overwriting the refmac map.
   //
   mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map, imol_map);
   mc.imol_difference_map = imol_diff_map; // happens for you in connect_updating_maps() (but we are not using that here).
   // After you have changed maps the firs time, add a starting point for the gru score:
   mc.calculate_new_rail_points();

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
   float new_rail_points = mc.calculate_new_rail_points();
   float gpt = mc.rail_points_total();
   std::cout << "###### RailPoints gained: " << new_rail_points << " rail points total " << gpt << std::endl;

   if (new_rail_points > 4.0)
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

int test_undo_and_redo_2(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   coot::atom_spec_t atom_spec("A", 61, "", " CZ ", "");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   coot::Cartesian pt_1(at_1->x, at_1->y, at_1->z);

   int status_af = mc.auto_fit_rotamer(imol, "A", 61, "", "", imol_map);
   if (status_af == 1) {
      coot::Cartesian pt_2(at_1->x, at_1->y, at_1->z);
      double dd = coot::Cartesian::lengthsq(pt_1, pt_2);
      double d = std::sqrt(dd);
      if (d > 6.0) {
         // OK, it moved (fitted)
         mc.undo(imol);
         mmdb::Atom *at_3 = mc.get_atom(imol, atom_spec);
         coot::Cartesian pt_3(at_3->x, at_3->y, at_3->z);
         dd = coot::Cartesian::lengthsq(pt_1, pt_3);
         d = std::sqrt(dd);
         std::cout << "debug:: in test_undo_and_redo_2() d " << d << std::endl;
         if (d < 0.001)
            status = 1;
      }
   }
   return status;
}


int test_ramachandran_analysis(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      unsigned int n_res = 0;
      coot::validation_information_t ra = mc.ramachandran_analysis(imol);
      for (const auto &chain : ra.cviv) {
         for (const auto &res : chain.rviv) {
            if (res.function_value > 0.5)
               n_res++;
         }
      }
      std::cout << "debug:: in test_ramachandran_validation n_res: " << n_res << std::endl;
      if (n_res > 100)
         status = 1;
   }
   mc.close_molecule(imol);
   return status;
}


int test_rama_validation(molecules_container_t &mc) {

   // test  that 285 is not there

   starting_test(__FUNCTION__);
   int status = 0;

   std::string coords_fn = reference_data("moorhen-tutorial-structure-number-1.pdb");
   int imol = mc.read_pdb(coords_fn);

   std::vector<coot::phi_psi_prob_t> rv = mc.ramachandran_validation(imol);

   bool r_285 = false;
   bool r_286 = false;
   for (const auto &r : rv) {
      // std::cout << " " << r.first << " " << r.second << std::endl;
      if (r.phi_psi.residue_number == 285) r_285 = true;
      if (r.phi_psi.residue_number == 286) r_286 = true;
   }

   if (r_286 && ! r_285)
      status = 1;

   return status;
}


int test_rama_balls_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string coords_fn = reference_data("moorhen-tutorial-structure-number-1.pdb");
   // std::string coords_fn = reference_data("frag.pdb");
   int imol = mc.read_pdb(coords_fn);

   coot::simple_mesh_t rvmm = mc.get_ramachandran_validation_markup_mesh(imol);
   std::cout << "debug:: rama mesh: " << rvmm.vertices.size() << " vertices and " << rvmm.triangles.size()
             << " triangles" << std::endl;

   // Let's look at the colours of the balls.
   if (false) // let's not.
      for (unsigned int i=0; i<rvmm.vertices.size(); i+=1) {
         // std::cout << i << " " << glm::to_string(rvmm.vertices[i].color) << std::endl;
         // const auto &v = rvmm.vertices[i];
         // std::cout << i << " " << glm::to_string(v.pos) << " " << glm::to_string(v.normal) << " " << glm::to_string(v.color)  << std::endl;
      }

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

   auto mmdb_to_cartesian = [] (mmdb::Atom *at) {
      return coot::Cartesian(at->x, at->y, at->z);
   };

   auto glm_to_cartesian = [] (const glm::vec3 &gp) {
      return coot::Cartesian(gp[0], gp[1], gp[2]);
   };

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
      mc.write_coordinates(imol, "test-add-terminal-residue-with-added-terminal-residue.pdb");
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

   // part three - is the mesh updated correctly to give a bond between 278 and 279?
   bool part_three_done = false;

   // what is the mid-point along the new peptide bond?
   mmdb::Atom *at_3 = mc.get_atom(imol, coot::atom_spec_t("A", 278, "", " C  ", ""));
   mmdb::Atom *at_4 = mc.get_atom(imol, coot::atom_spec_t("A", 279, "", " N  ", ""));
   if (at_3 && at_4) {

      std::cout << "got peptide atoms " << std::endl;
      coot::Cartesian pt_3 = mmdb_to_cartesian(at_3);
      coot::Cartesian pt_4 = mmdb_to_cartesian(at_4);

      double ddt = coot::Cartesian::lengthsq(pt_3, pt_4);
      double dt = std::sqrt(ddt);
      std::cout << " distance between peptide atoms " << dt << std::endl;
      coot::Cartesian mid_point = pt_3.mid_point(pt_4);
      std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
      auto mesh = mc.get_bonds_mesh(imol, mode, true, 0.1, 1.0, 1);
      double d_crit = 0.15;
      unsigned int n_vertex = 0;
      for (const auto &vert : mesh.vertices) {
         const auto &pos = vert.pos;
         coot::Cartesian pos_cart = glm_to_cartesian(pos);
         double dd = coot::Cartesian::lengthsq(mid_point, pos_cart);
         double d = std::sqrt(dd);
         if (d < d_crit) {
            n_vertex++;
         }
      }
      if (n_vertex > 2)
         part_three_done = true;
   }

   // part four - add to a residue that has an OXT atom
   bool part_four_done = false;
   mmdb::Atom *oxt_1 = mc.get_atom(imol, coot::atom_spec_t("A", 303, "", " OXT", ""));
   if (oxt_1) {
      mc.add_terminal_residue_directly(imol, "A", 303, "");
      mmdb::Atom *oxt_2 = mc.get_atom(imol, coot::atom_spec_t("A", 303, "", " OXT", ""));
      if (oxt_2 == nullptr) // not there
         part_four_done = true;
   }

   if (part_one_done && part_two_done && part_three_done && part_four_done)
      status = 1;

   mc.close_molecule(imol);
   mc.close_molecule(imol_map);

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
   coot::simple_mesh_t mesh = mc.get_bonds_mesh(imol, mode, true, 0.1, 1.0, 1);
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

int test_no_dictionary_residues(molecules_container_t &mc) {

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
         float rmsd = mc.get_map_rmsd_approx(imol_diff_map);
         float n_rmsd = 4.5; // lots of high peaks in the tutorial data. 6.0 is not a normal limit
         auto sites = mc.difference_map_peaks(imol_diff_map, imol, n_rmsd);
         if (sites.size() > 8)
            status = 1;
         if (false)
            for (const auto &site : sites) {
               float peak_n_rmsd = site.feature_value / rmsd;
               std::cout << "site " << site.feature_type << " " << site.button_label << " "
                         << site.x << " " << site.y << " " << site.z << " residue " << site.residue_spec
                         << " height " << site.feature_value
                         << " n-rmsd " << peak_n_rmsd
                         << " badness " << site.badness
                         << std::endl;
            }
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

int test_dictionary_bonds(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("pdb2sar-part.ent"));
   mc.import_cif_dictionary("ATP.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("3GP.cif", imol_1);
   int imol_2 = mc.get_monomer("ATP");
   int imol_3 = mc.read_pdb(reference_data("pdb2sar-part.ent"));

   std::cout << ":::: test_dictionary_bonds() imol_2: " << imol_2 << std::endl;
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");

   glm::vec3 atom_ligand_C4_position(53.4, 9.7, 20.3);

   coot::simple_mesh_t mesh = mc.get_bonds_mesh(imol_3, mode, true, 0.1, 1.0, 1);

   // there is no dictionary, but we should see vertices for the atoms
   //
   unsigned int n_ligand_vertices = 0;
   for (const auto &vert : mesh.vertices) {
      double d = glm::distance(vert.pos, atom_ligand_C4_position);
      if (d < 1.0)
         n_ligand_vertices++;
   }
   std::cout << "debug:: test_dictionary_bonds n_ligand_vertices: " << n_ligand_vertices << std::endl;
   if (n_ligand_vertices > 0)
      status = 1;

   return status;
}

int test_transformation_for_atom_selection(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 20, "", " CA ", "");
   
   //! get the atom position (function for testing)
   std::pair<bool, coot::Cartesian> pos = mc.get_atom_position(imol, atom_spec);
   if (pos.first) {
      coot::Cartesian pos_pre = pos.second;
      int n_atoms = 7; // VAL
      int n_atoms_moved  =  mc.apply_transformation_to_atom_selection(imol, "//A/20", n_atoms,
                                                1,0,0,
                                                0,1,0,
                                                0,0,1,
                                                0,0,0, // centre of rotation
                                                2,3,4);
      pos = mc.get_atom_position(imol, atom_spec);
      if (n_atoms_moved == 7) {
         if (pos.first) {
            coot::Cartesian pos_post = pos.second;
            double dd = coot::Cartesian::lengthsq(pos_pre, pos_post);
            double d = std::sqrt(dd);
            std::cout << "debug:: in test_transformation_for_atom_selection() d " << d << std::endl;
            if (d > 4.0)
               if (d < 10)
                  status = 1;
         }
      }
   }
   return status;
}

int test_new_position_for_atoms(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 20, "", " CA ", "");
   
   //! get the atom position (function for testing)
   std::pair<bool, coot::Cartesian> pos = mc.get_atom_position(imol, atom_spec);
   if (pos.first) {
      coot::Cartesian pos_pre = pos.second;

      std::vector<coot::molecule_t::moved_atom_t> moved_atom_positions;
      std::vector<std::string> atom_names = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2"};
      for (const auto &name : atom_names) {
         coot::atom_spec_t as("A", 20, "", name, "");
         moved_atom_positions.push_back(coot::molecule_t::moved_atom_t(name, "", 2.1, 3.2, 4.3));
      }

      mc.new_positions_for_residue_atoms(imol, "//A/20", moved_atom_positions);
      pos = mc.get_atom_position(imol, atom_spec);
      if (pos.first) {
         coot::Cartesian pos_post = pos.second;
         double dd = coot::Cartesian::lengthsq(pos_pre, pos_post);
         double d = std::sqrt(dd);
         std::cout << "debug:: in test_new_position_for_atoms() d " << d << std::endl;
         if (d > 70.0)
            status = 1;
      }
   }
   return status;
}

int test_new_position_for_atoms_in_residues(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 20, "", " CA ", "");
   
   //! get the atom position (function for testing)
   std::pair<bool, coot::Cartesian> pos = mc.get_atom_position(imol, atom_spec);
   if (pos.first) {
      coot::Cartesian pos_pre = pos.second;

      std::vector<coot::molecule_t::moved_atom_t> moved_atom_positions;
      std::vector<std::string> atom_names = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2"};
      std::vector<coot::molecule_t::moved_residue_t> moved_atoms_residue_vec;
      coot::molecule_t::moved_residue_t mr("A", 20, "");
      for (const auto &name : atom_names) {
         coot::atom_spec_t as("A", 20, "", name, "");
         mr.add_atom(coot::molecule_t::moved_atom_t(name, "", 2.1, 3.2, 4.3));
      }
      moved_atoms_residue_vec.push_back(mr);

      mc.new_positions_for_atoms_in_residues(imol, moved_atoms_residue_vec);
      pos = mc.get_atom_position(imol, atom_spec);
      if (pos.first) {
         coot::Cartesian pos_post = pos.second;
         double dd = coot::Cartesian::lengthsq(pos_pre, pos_post);
         double d = std::sqrt(dd);
         std::cout << "debug:: in test_new_position_for_atoms() d " << d << std::endl;
         if (d > 70.0)
            status = 1;
      }
   }
   return status;
}



int test_merge_molecules(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   mc.import_cif_dictionary("ATP.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("3GP.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("NUT.cif", coot::protein_geometry::IMOL_ENC_ANY);
  
   int imol_2 = mc.get_monomer_and_position_at("ATP", coot::protein_geometry::IMOL_ENC_ANY, 60, 50, 30);
   int imol_3 = mc.get_monomer_and_position_at("3GP", coot::protein_geometry::IMOL_ENC_ANY, 80, 55, 20);
   int imol_4 = mc.get_monomer_and_position_at("NUT", coot::protein_geometry::IMOL_ENC_ANY, 10, 15, 10);

   // imol_2 and imol_3 should be given chain_id "A" because they are close to "A"
   // and imol_4 should have chain-id "B" because it is far from "A"

   std::vector<std::string> chains_ids_pre = mc.get_chains_in_model(imol_1);
   std::string ls = std::to_string(imol_2) + std::string(":") + std::to_string(imol_3) + std::string(":") + std::to_string(imol_4);

   std::pair<int, std::vector<merge_molecule_results_info_t> > merge_results =
      mc.merge_molecules(imol_1, ls);
   std::vector<std::string> chain_ids_post = mc.get_chains_in_model(imol_1);

   // std::cout << "in test_merge_molecules() chain_ids_post size " << chain_ids_post.size() << std::endl;
   mc.write_coordinates(imol_1, "post-merge.pdb");
   std::set<std::string> chain_ids_set;
   for (const auto &ch : chain_ids_post)
      chain_ids_set.insert(ch);

   if (chain_ids_post.size() == 2) {
         status = 1;
   } else {
      std::cout << "in test_merge_molecules() failed to make unique chain ids " << chain_ids_post.size() << std::endl;
   }

   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   coot::simple_mesh_t mesh = mc.get_bonds_mesh(imol_1, mode, true, 0.1, 1.0, 1);

   return status;
}

int test_density_correlation_validation(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, true);
   bool bad_correls = false;
   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {
         unsigned int n_res = 0;
         coot::validation_information_t dca = mc.density_correlation_analysis(imol, imol_map);
         for (const auto &chain : dca.cviv) {
            for (const auto &res : chain.rviv) {
               if (res.function_value > 0.5)
                  n_res++;
               if (res.function_value < -1.0) bad_correls = true;
               if (res.function_value >  1.0) bad_correls = true;

               std::cout << "correl " << res.residue_spec.res_no << " " << res.function_value << std::endl;
            }
         }
         std::cout << "debug:: in test_density_correlation_validation n_res: " << n_res << std::endl;
         if (bad_correls == false) {
            if (n_res > 400)
               status = 1;
         } else {
            std::cout << "debug:: in test_density_correlation_validation() bad correls! " << std::endl;
         }
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   return status;
}

int test_rotamer_validation(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      unsigned int n_res = 0;
      coot::validation_information_t dca = mc.rotamer_analysis(imol);
      for (const auto &chain : dca.cviv) {
         for (const auto &res : chain.rviv) {
            if (res.function_value > 0.5)
               n_res++;
         }
      }
      std::cout << "debug:: in test_rotamer_validation n_res: " << n_res << std::endl;
      if (n_res > 150)
         status = 1;
   }
   mc.close_molecule(imol);
   return status;
}

int test_add_water(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {
         unsigned int n_waters = mc.add_waters(imol, imol_map);
         if (n_waters > 10)
            status = 1;
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   return status;
}

int test_read_a_map(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   bool is_diff_map = false;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_ccp4_map(reference_data("test.map"), is_diff_map);
   std::cout << "Here in test_read_a_map() with imol_map " << imol_map << std::endl;
   if (mc.is_valid_map_molecule(imol_map)) {

      float radius = 20;
      float contour_level = 0.013;
      coot::Cartesian p(88.25823211669922, 69.19033813476562, 89.1391372680664);
      coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
      std::cout << "DEBUG:: test_density_mesh(): " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
                << " triangles" << std::endl;

      if (map_mesh.vertices.size() > 30000)
         status = 1;
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   return status;

}

int test_ligand_fitting_here(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   int imol_ligand = mc.get_monomer("TRS");

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_model_molecule(imol_ligand)) {
         float n_rmsd = 1.0;
         coot::Cartesian pos(62.4, 40.7, 26.6);
         std::vector<int> fits = mc.fit_ligand_right_here(imol, imol_map, imol_ligand, pos.x(), pos.y(), pos.z(), n_rmsd, true, 10);
         std::cout << "debug:: in test_ligand_fitting_here() found " << fits.size() << " fits" << std::endl;
         if (fits.size() == 1) {
            int imol_fit = fits[0];
            mc.display_molecule_names_table();
            std::cout << "debug:: in test_fit_ligand() imol_ligand: " << imol_ligand << " imol_fit: " << imol_fit << std::endl;
            coot::validation_information_t dca_1 = mc.density_correlation_analysis(imol_ligand, imol_map);
            coot::validation_information_t dca_2 = mc.density_correlation_analysis(imol_fit,    imol_map);
            std::cout << "debug:: in test_ligand_fitting_here() dca_2 has " << dca_2.cviv.size() << " chains" << std::endl;
            if (! dca_2.empty()) {
               double cc_1 = dca_1.cviv[0].rviv[0].function_value;
               double cc_2 = dca_2.cviv[0].rviv[0].function_value;
               std::cout << "debug:: in test_ligand_fitting_here() cc: " << cc_1 << " " << cc_2 << std::endl;
               // this is not a good test, because this model doesn't have a ligand. Here we basically just test
               // that the function added some atoms.
               if (cc_2 > 0.2)
                  status = 1;
            }
         }
      } else {
         std::cout << "debug:: test_fit_ligand() failed to get model for ligand " << imol_ligand << std::endl;
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   mc.close_molecule(imol_ligand);
   return status;
}

int test_jiggle_fit(molecules_container_t &mc) {

   // 20221119-PE this needs a better test. I need to construct a problem
   // where there is a good solution with a real ligand.

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);


   if (mc.is_valid_model_molecule(imol)) {
      mc.imol_refinement_map = imol_map;
      coot::atom_spec_t atom_spec("A", 61, "", " CZ ","");
      coot::residue_spec_t residue_spec("A", 61, "");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
         mc.fit_to_map_by_random_jiggle(imol, residue_spec, 9110, -1);
         // you can test that the density for the CZ has improved when that function is available.
         // or maybe "density_fit_for_residue()" ?
         // or maybe "density_correlation_for_residue()" ?
         // or maybe "density_correlation_for_residues()" with a list of a single residue
         mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec);
         coot::Cartesian atom_pos_2 = atom_to_cartesian(at_2);
         double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
         double d = std::sqrt(dd);
         std::cout << "test_jiggle_fit d " << d << std::endl;
         if (d > 0.4)
            status = true;
      } else {
         std::cout << "ERROR: test_jiggle_fit() missing atom " << atom_spec << std::endl;
      }
   } else {
      std::cout << "ERROR: test_jiggle_fit() invalid model molecule " << imol << std::endl;
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);


   if (status == 1) {

      status = 0;
      std::cout << "Second jiggle-fit test: using atom selection ------------------------------" << std::endl;
      imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
      int imol_start = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
      int imol_other = mc.read_pdb("weird-orientation.pdb");
      if (mc.is_valid_model_molecule(imol_other)) {
         // now test that we stared with bad fit to density
         coot::validation_information_t vi_0 = mc.density_correlation_analysis(imol_start, imol_map);
         coot::validation_information_t vi_1 = mc.density_correlation_analysis(imol_other, imol_map);
         // fit!
         int imol_blur = mc.sharpen_blur_map(imol_map, 100, false);
         mc.imol_refinement_map = imol_blur;
         mc.fit_to_map_by_random_jiggle_using_cid(imol_other, "//A", 10100, 1);
         // now test that we have good fit to density.
         coot::validation_information_t vi_2 = mc.density_correlation_analysis(imol_other, imol_map);

         coot::stats::single s_0 = vi_0.get_stats();
         coot::stats::single s_1 = vi_1.get_stats();
         coot::stats::single s_2 = vi_2.get_stats();

         // 20230402-PE These results are disappointing - they are not as good as doing it interactively.
         // I wonder what the difference is.

         std::cout << "orig:     " << std::fixed << s_0.mean() << " " << std::fixed << std::sqrt(s_0.variance()) << std::endl;
         std::cout << "pre-fit:  " << std::fixed << s_1.mean() << " " << std::fixed << std::sqrt(s_1.variance()) << std::endl;
         std::cout << "post-fit: " << std::fixed << s_2.mean() << " " << std::fixed << std::sqrt(s_2.variance()) << std::endl;

         float d1 = s_0.mean() - s_2.mean();
         float d2 = s_2.mean() - s_1.mean();

         if (d2 > d1) // we got most of the way there
            status = 1;

         mc.write_coordinates(imol_other, "jiggled.pdb");
      }
   }
   return status;
}

int test_peptide_omega(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      auto vi = mc.peptide_omega_analysis(imol);
      for (const auto &chain : vi.cviv) {
         for (const auto &res : chain.rviv) {
            // std::cout << "in test_peptide_omega() " << res.residue_spec << " " << res.function_value << std::endl;
            if (res.residue_spec.res_no == 32) {
               if (res.function_value > 32.0) { // not planned
                  status = 1;
               }
            }
         }
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_delete_literal(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {

      int status_1 = mc.delete_using_cid(imol, "//A/10-20", "LITERAL").first;

      if (status_1) {

         unsigned int n_found = 0;
         for (unsigned int ires=9; ires<=21; ires++) {
            coot::atom_spec_t atom_spec("A", ires, "", " O  ","");
            mmdb::Atom *at = mc.get_atom(imol, atom_spec);
            if (at)
               n_found++;
         }

         if (n_found == 2)
            status = 1;
      }
   }

   // add another test here that it doesn't delete anything if the atom selection doesn't match
   
   mc.close_molecule(imol);
   return status;
}

int test_cis_trans(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 263, "", " N  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
         int status_l = mc.cis_trans_convert(imol, "//A/262/O");
         if (status_l) {
            mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec);
            coot::Cartesian atom_pos_2 = atom_to_cartesian(at_2);
            double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
            double d = std::sqrt(dd);
            std::cout << "test_cis_trans d " << d << std::endl;
            if (d > 1.0)
               status = true;
            // mc.write_coordinates(imol, "cis-trans-converted.pdb");
         }
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_add_compound(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         coot::Cartesian atom_pos = atom_to_cartesian(at_1);
         double dd = coot::Cartesian::lengthsq(atom_pos, atom_pos);
         double d = std::sqrt(dd);

         coot::Cartesian lig_pos(1,23,3);
         mc.add_compound(imol, "GOL", coot::protein_geometry::IMOL_ENC_ANY, imol_map, lig_pos.x(), lig_pos.y(), lig_pos.z());

         coot::validation_information_t dca = mc.density_correlation_analysis(imol, imol_map);
         for (const auto &chain : dca.cviv) {
            for (const auto &res : chain.rviv) {
               if (res.residue_spec.res_no == 62) {
                  std::cout << "function value " << res.function_value << std::endl;
                  if (res.function_value > 0.6) {
                     status = 1;
                     mc.write_coordinates(imol, "post-add-compound.pdb");
                  }
               }
            }
         }
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_non_standard_residues(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string fn = "moorhen-tutorial-structure-number-1.pdb";
   int imol     = mc.read_pdb(reference_data(fn));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      coot::Cartesian lig_pos_1(1,23,3);
      coot::Cartesian lig_pos_2(63,40,27);
      mc.add_compound(imol, "GOL", coot::protein_geometry::IMOL_ENC_ANY, imol_map, lig_pos_1.x(), lig_pos_1.y(), lig_pos_1.z());
      // should the MPD be added to the A chain?
      int status_l = mc.add_compound(imol, "MPD", coot::protein_geometry::IMOL_ENC_ANY, imol_map, lig_pos_2.x(), lig_pos_2.y(), lig_pos_2.z());
      std::cout << "status_l " << status_l << std::endl;
      if (status_l) {
         //mc.write_coordinates(imol, "test_non_standard_residues.pdb");
         auto specs = mc.get_non_standard_residues_in_molecule(imol);
         std::cout << "DEBUG:: there were " << specs.size() << " non-standard residues in " << fn << std::endl;
         for (const auto &spec : specs) {
            std::cout << "    " << spec << " " << spec.string_user_data << std::endl;
         }
         if (specs.size() == 2)
            status = 1;
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_molecular_representation(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      std::string selection = "//";
      std::string colour = "colorRampChainsScheme";
      std::string style = "Ribbon";

      colour = "Chain";
      colour = "RampChains";
      // style  = "MolecularSurface";

      coot::simple_mesh_t mesh = mc.get_molecular_representation_mesh(imol, selection, colour, style);
      if (mesh.vertices.size() > 10) {

         std::cout << "test_molecular_representation() Ribbons OK" << std::endl;

         style = "MolecularSurface";
         coot::simple_mesh_t surface_mesh = mc.get_molecular_representation_mesh(imol, selection, colour, style);

         status = true;
      }
   }
   mc.close_molecule(imol);
   return status;
}


int test_replace_fragment(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 15, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {

         coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
         int imol_frag = mc.copy_fragment_using_cid(imol, "//A/10-20");
         int n_atoms = 77;
         int n_atoms_moved = mc.apply_transformation_to_atom_selection(imol_frag, "//", n_atoms,
                                                                       1, 0, 0,
                                                                       0, 1, 0,
                                                                       0, 0, 1,
                                                                       0, 0, 0,
                                                                       1, 2, 3);
         if (n_atoms_moved > 20) {
            mmdb::Atom *at_2 = mc.get_atom(imol_frag, atom_spec);
            coot::Cartesian atom_pos_2 = atom_to_cartesian(at_2);
            double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
            double d = std::sqrt(dd);
            std::cout << "test_replace_fragment d " << d << std::endl;
            int status_2 = mc.replace_fragment(imol, imol_frag, "//A/12-17");
            if (status_2 == 1) {
               mmdb::Atom *at_3 = mc.get_atom(imol, atom_spec);
               coot::Cartesian atom_pos_3 = atom_to_cartesian(at_2);
               dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_3);
               d = std::sqrt(dd);
               if (d > 3.0)
                  status = 1;
            } else {
               std::cout << "test_replace_fragment() bad status_2 " << std::endl;
            }
         } else {
            std::cout << "test_replace_fragment() n_atoms_moved " << n_atoms_moved << std::endl;
         }
         
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_instanced_rota_markup(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::simple_mesh_t mz = mc.get_rotamer_dodecs(imol);
      coot::instanced_mesh_t m = mc.get_rotamer_dodecs_instanced(imol);
      if (! m.geom.empty()) {
         const coot::instanced_geometry_t &ig = m.geom[0];
         if (ig.vertices.size() > 30) {
            if (ig.triangles.size() > 30) {
               if (ig.instancing_data_A.size() > 30) {
                  status = 1;

                  if (false) {
                     for (unsigned int i=0; i<ig.vertices.size(); i++) {
                        const auto &vert = ig.vertices[i];
                        std::cout << "   " << i << " " << glm::to_string(vert.pos) << std::endl;
                     }
                  }

                  for (unsigned int i=0; i<ig.instancing_data_A.size(); i++) {
                     if (i > 20) continue;
                     const auto &item = ig.instancing_data_A[i];
                     std::cout << i << " "
                               << glm::to_string(item.position) << " "
                               << glm::to_string(item.colour)   << " "
                               << glm::to_string(item.size)     << std::endl;
                  }
               } else {
                  std::cout << "error:: in test_instanced_rota_markup() instancing_data_A size " << ig.instancing_data_A.size() << std::endl;
               }
            } else {
               std::cout << "error:: in test_instanced_rota_markup() triangles size " << ig.triangles.size() << std::endl;
            }
         } else {
            std::cout << "error:: in test_instanced_rota_markup() vertices size " << ig.vertices.size() << std::endl;
         }
      } else {
         std::cout << "error:: in test_instanced_rota_markup() geom is empty " << std::endl;
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_gaussian_surface(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      float sigma = 4.4;
      float contour_level = 4.0;
      float box_radius = 5.0;
      float grid_scale = 0.7;
      coot::simple_mesh_t mesh = mc.get_gaussian_surface(imol, sigma, contour_level, box_radius, grid_scale);
      std::cout << "in test_gaussian_surface() " << mesh.vertices.size() << " " << mesh.triangles.size() << std::endl;
      if (mesh.vertices.size() > 0)
         if (mesh.triangles.size() > 0)
            status = 1;
   }
   mc.close_molecule(imol);
   return status;
}

int test_instanced_bonds_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   if (mc.is_valid_model_molecule(imol)) {
      coot::instanced_mesh_t im = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, 1);
      std::cout << "instanced mesh has " << im.geom.size()  << " geoms" << std::endl;
      if (im.geom.size() > 3) {
         if (false) {
            for (unsigned int i=0; i<im.geom[0].instancing_data_A.size(); i++) {
               std::cout << "   " << i << " " << glm::to_string(im.geom[0].instancing_data_A[i].position) << std::endl;
            }
         }
         if (im.geom[0].instancing_data_A.size() > 1000)
            status = 1;
      }

      for (unsigned int i=0; i<im.geom.size(); i++) {
         const auto &geom = im.geom[i];
         std::cout << "geom i = " << i << " has " << geom.vertices.size() << " vertices" << std::endl;

         if (i == 2) {
            for (unsigned int j=0; j<geom.instancing_data_B.size(); j++) {
               if (j >= 10) continue;
               std::cout << j << " " << glm::to_string(geom.instancing_data_B[j].orientation) << std::endl;
            }
         }
      }
   }

   std::string cid("/*/A/270");
   coot::instanced_mesh_t im_lig = mc.get_bonds_mesh_for_selection_instanced(imol, cid, mode, true, 0.1, 1.0, 1);
   unsigned int n_geoms = im_lig.geom.size();
   for (unsigned int i=0; i<n_geoms; i++) {
      std::cout << "test_instanced_bonds_mesh()) im_lig " << im_lig.geom[i].name << " " << i << " has A " << im_lig.geom[i].instancing_data_A.size() << std::endl;
      std::cout << "test_instanced_bonds_mesh()) im_lig " << im_lig.geom[i].name << " " << i << " has B " << im_lig.geom[i].instancing_data_B.size() << std::endl;
      for (unsigned int j=0; j<im_lig.geom[i].instancing_data_B.size(); j++) {
         const auto &item = im_lig.geom[i].instancing_data_B[j];
         std::cout << "   instanced bond: " << j << " has colour " << glm::to_string(item.colour) << std::endl;
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_add_alt_conf(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      mc.add_alternative_conformation(imol, "//A/72");
      mc.write_coordinates(imol, "with-alt-conf.pdb");

      coot::residue_spec_t res_spec("A", 72, "");
      mmdb::Residue *r = coot::util::get_residue(res_spec, mc[imol].atom_sel.mol);

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      r->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::cout << iat << " " << coot::atom_spec_t(at) << " " << at->x << " " << at->y << " " << at->z << std::endl;
         }
      }
      if (n_residue_atoms > 22)
         status = 1;
   }

   return status;
}

int test_fill_partial(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int part_1 = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   mc.set_imol_refinement_map(imol_map);

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " OG ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         mc.delete_atom_using_cid(imol, "//A/270/OG");
         at_1 = mc.get_atom(imol, atom_spec);
         if (at_1) {
            std::cout << "failed to delete " << std::endl;
         } else {
            mc.fill_partial_residue_using_cid(imol, "//A/270");
            at_1 = mc.get_atom(imol, atom_spec);
            if (at_1) {
               std::cout << "Found the OG" << std::endl;
               part_1 = 1;
            }
         }
      }
   }

   // now test the all-molecule function

   if (part_1 == 1) {
      // there are more things to be filled than just these 2 residues
      mc.delete_atom_using_cid(imol, "//A/43/CG");
      mc.delete_atom_using_cid(imol, "//A/44/CG1");
      mmdb::Atom *at_1 = mc.get_atom_using_cid(imol, "//A/43/CG");
      mmdb::Atom *at_2 = mc.get_atom_using_cid(imol, "//A/44/CG1");
      if (at_1) {
         std::cout << "fail to delete 1" << std::endl;
      } else {
         if (at_2) {
            std::cout << "fail to delete 2 " << std::endl;
         } else {
            mc.fill_partial_residues(imol);
            at_1 = mc.get_atom_using_cid(imol, "//A/43/CG");
            at_2 = mc.get_atom_using_cid(imol, "//A/44/CG1");
            if (at_1) {
               if (at_2) {
                  status = 1;
               }
            }
         }
      }
   }
   return status;
}

int test_missing_atoms_info(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      mc.delete_atom_using_cid(imol, "//A/43/CG");
      mc.delete_atom_using_cid(imol, "//A/44/CG1");
      coot::util::missing_atom_info mai = mc.missing_atoms_info_raw(imol);
      std::cout << "missing_atom_info: residues with no dictionary: size " << mai.residues_with_no_dictionary.size() << std::endl;
      for (unsigned int i=0; i<mai.residues_with_no_dictionary.size(); i++) {
         std::cout << "   " << mai.residues_with_no_dictionary[i] << std::endl;
      }
      std::cout << "missing_atom_info: residues with missing atoms: size " << mai.residues_with_missing_atoms.size() << std::endl;
      for (unsigned int i=0; i<mai.residues_with_missing_atoms.size(); i++) {
         mmdb::Residue *r = mai.residues_with_missing_atoms[i];
         std::cout << "   with missing atoms: " << coot::residue_spec_t(r)<< " " << r->GetResName() << std::endl;
      }
      if (mai.residues_with_missing_atoms.size() > 1)
         status = 1;
   }
   return status;
}

int test_editing_session_tutorial_1(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol          = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map      = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "DELFWT", "PHDELWT", "W", false, true);

   if (mc.is_valid_model_molecule(imol)) {
      coot::simple_mesh_t map_mesh; // just throw it away for testing

      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);

      mc.set_imol_refinement_map(imol_map);
      mc.associate_data_mtz_file_with_map(imol_map, reference_data("moorhen-tutorial-map-number-1.mtz"), "FP", "SIGFP", "FREE");
      mc.connect_updating_maps(imol, imol_map, imol_map, imol_diff_map);
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      molecules_container_t::r_factor_stats stats_1 = mc.get_r_factor_stats();
      std::cout << "stats_1: " << mc.r_factor_stats_as_string(stats_1) << std::endl;

      // debugging
      mc.display_molecule_names_table();
      mc.add_waters(imol, imol_map);

      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      molecules_container_t::r_factor_stats stats_2 = mc.get_r_factor_stats();
      std::cout << "stats_2: " << mc.r_factor_stats_as_string(stats_2) << std::endl;
      if (stats_2.rail_points_total > 500)
         status = 1;
   }
   return status;
}

int test_editing_session_tutorial_4(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol          = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map      = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "DELFWT", "PHDELWT", "W", false, true);

   if (mc.is_valid_model_molecule(imol)) {
      mc.set_imol_refinement_map(imol_map);
      mc.associate_data_mtz_file_with_map(imol_map, reference_data("moorhen-tutorial-map-number-4.mtz"), "F", "SIGF", "FREER");
      mc.connect_updating_maps(imol, imol_map, imol_map, imol_diff_map);

      // debugging
      mc.display_molecule_names_table();

      coot::simple_mesh_t map_mesh; // just throw it away for testing

      mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map, imol_map);
      int rpn_1 = mc.calculate_new_rail_points();
      int rpt_1 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points A latest_move: " << rpn_1 << " total: " << rpt_1 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.flip_peptide_using_cid(imol, "//A/20/C", "");
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_2 = mc.calculate_new_rail_points();
      int rpt_2 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points B: latest_move: " << rpn_2 << " total: " << rpt_2 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.flip_peptide_using_cid(imol, "//A/262/C", "");
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_3 = mc.calculate_new_rail_points();
      int rpt_3 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points C: latest_move: " << rpn_3 << " total: " << rpt_3 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.flip_peptide_using_cid(imol, "//A/33/C", "");
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_4 = mc.calculate_new_rail_points();
      int rpt_4 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points D: latest_move: " << rpn_4 << " total: " << rpt_4 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.fill_partial_residue_using_cid(imol, "//A/205");
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_5 = mc.calculate_new_rail_points();
      int rpt_5 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points E: latest_move: " << rpn_5 << " total: " << rpt_5 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.cis_trans_convert(imol, "//A/262/CA");
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_6 = mc.calculate_new_rail_points();
      int rpt_6 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points F: latest_move: " << rpn_6 << " total: " << rpt_6 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.refine_residues_using_atom_cid(imol, "//A/262/CA", "QUINTUPLE");
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_7 = mc.calculate_new_rail_points();
      int rpt_7 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points G: latest_move: " << rpn_7 << " total: " << rpt_7 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.auto_fit_rotamer(imol, "A", 58, "", "", imol_map);
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_8 = mc.calculate_new_rail_points();
      int rpt_8 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points G: latest_move: " << rpn_8 << " total: " << rpt_8 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      mc.auto_fit_rotamer(imol, "A", 61, "", "", imol_map);
      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      int rpn_9 = mc.calculate_new_rail_points();
      int rpt_9 = mc.rail_points_total();
      std::cout << "::::::::::::::::::::::::::::::::::: Rail points G: latest_move: " << rpn_9 << " total: " << rpt_9 << std::endl;
      std::cout << "::::::::::::::::::::::::::::::::::: R-factor " << mc.get_latest_sfcalc_stats().r_factor << std::endl;

      if (rpt_9 > 500)
         status = 1;
   }

   return status;
}

int test_ligand_contact_dots(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      unsigned int num_subdivisions = 1;
      coot::instanced_mesh_t im = mc.contact_dots_for_ligand(imol, "262", num_subdivisions);
      if (im.geom.size() > 1) {
         if (im.geom[0].instancing_data_A.size() > 10) status = 1;
         if (im.geom[1].instancing_data_A.size() > 10) status = 1;
      }
      for (unsigned int i=0; i<im.geom.size(); i++) {
         std::cout << "geom " << i
                   << " A: " << im.geom[i].instancing_data_A.size()
                   << " B: " << im.geom[i].instancing_data_B.size()
                   << std::endl;
      }
   }
   return status;
}

int test_broken_function(molecules_container_t &mc) {

   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   auto mesh = mc.get_bonds_mesh(imol, mode, true, 0.1, 1.0, 1);
   return status;

}

int test_delete_side_chain(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec_1("A", 270, "", " O  ","");
      coot::atom_spec_t atom_spec_2("A", 270, "", " OG ","");
      coot::atom_spec_t atom_spec_3("A",  51, "", " CG ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec_1);
      mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec_2);
      if (at_1 && at_2) {
         mc.delete_side_chain(imol, "A", 270, "");
         mc.delete_side_chain_using_cid(imol, "//A/51");
         at_1 = mc.get_atom(imol, atom_spec_1);
         at_2 = mc.get_atom(imol, atom_spec_2);
         mmdb::Atom *at_3 = mc.get_atom(imol, atom_spec_3);
         if (at_1)
            if (! at_2)
               if (! at_3)
                  status = 1;
      }
   }
   return status;
}

int test_colour_rules(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string pdb_fn = "pdb7vvl.ent";
   unsigned int n_colour_rules = 6; // number expected
   int imol_0 = mc.read_pdb(reference_data(pdb_fn));
   if (imol_0 == -1) {
      pdb_fn = "moorhen-tutorial-structure-number-1.pdb";
      imol_0 = mc.read_pdb(reference_data(pdb_fn));
      n_colour_rules = 1;
   }
   int imol_1 = mc.read_pdb(reference_data(pdb_fn));
   int imol_2 = mc.read_pdb(reference_data(pdb_fn));

   auto v = mc.get_colour_rules(imol_0);

   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   auto mesh = mc.get_bonds_mesh_instanced(imol_0, mode, true, 0.1, 1.0, 1);

   if (true) {
      std::cout << "colour rules: " << std::endl;
      std::cout << "-------------" << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
         std::cout << i << " " << v[i].first << " " << v[i].second << std::endl;
      }
      std::cout << "-------------" << std::endl;
   }

   if (v.size() == n_colour_rules) status = 1;

   return status;
}


int test_multi_colour_rules(molecules_container_t &mc) {

   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::string crs = "//A/1^#cc0000|//A/2^#cb0002|//A/3^#c00007";
   mc.add_colour_rules_multi(imol, crs);
   auto v = mc.get_colour_rules(imol);
   if (v.size() == 4)
      status = 1;
   std::cout << "n colour rules " << v.size() << std::endl;
   return status;
}



int test_add_hydrogen_atoms(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int n_atoms_0 = mc.get_number_of_atoms(imol);
   auto mesh_0 = mc.get_bonds_mesh(imol, mode, true, 0.1, 1.0, 1);
   mc.add_hydrogen_atoms(imol);
   auto mesh_1 = mc.get_bonds_mesh(imol, mode, true, 0.1, 1.0, 1);
   int n_atoms_1 = mc.get_number_of_atoms(imol);
   mc.delete_hydrogen_atoms(imol);
   int n_atoms_2 = mc.get_number_of_atoms(imol);
   auto mesh_2 = mc.get_bonds_mesh(imol, mode, true, 0.1, 1.0, 1);

   if (n_atoms_1 > n_atoms_0) {
      if (n_atoms_2 == n_atoms_0) {
         int imol_lig = mc.get_monomer("GLC");
         auto mesh_lig_1 = mc.get_bonds_mesh(imol_lig, mode, true, 0.1, 1.0, 1);
         int n_atom_pre = mc.get_number_of_atoms(imol_lig);
         mc.delete_hydrogen_atoms(imol_lig);
         int n_atoms_post = mc.get_number_of_atoms(imol_lig);
         std::cout << "test_add_hydrogen_atoms(): pre: " << n_atom_pre << " n_atom_post " << n_atoms_post << std::endl;
         auto mesh_lig_2 = mc.get_bonds_mesh(imol_lig, mode, true, 0.1, 1.0, 1);
         status = 1;
      }
   }

   return status;
}

#include <fstream>

int test_mmrrcc(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      std::string chain_id = "A";
      std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
                std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> > results =
         mc.mmrrcc(imol, chain_id, imol_map);
      auto mc = results.first;
      auto sc = results.second;

      if (mc.size() > 90)
         status = 1;

      std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
      for (it=mc.begin(); it!=mc.end(); ++it)
         std::cout << "   " << it->first << " " << it->second.correlation() << std::endl;
   }
   mc.close_molecule(imol);


   if (false) {
      // Filo's example 11729 and 7adk
      imol = mc.read_pdb("pdb7adk.ent");
      imol_map = mc.read_ccp4_map("emd_11729.map", 0);
      auto results = mc.mmrrcc(imol, "B", imol_map);
      auto mc = results.first;
      auto sc = results.second;
      std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
      std::ofstream f("7adk-b-chain-all-atom.table");
      for (it=mc.begin(); it!=mc.end(); ++it)
         f << "   " << it->first << " " << it->second.correlation() << std::endl;
      std::ofstream fs("7adk-b-chain-side-chain.table");
      for (it=sc.begin(); it!=sc.end(); ++it)
         fs << "   " << it->first << " " << it->second.correlation() << std::endl;
   }
   return status;
}

int test_auto_read_mtz(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   std::vector<int> imol_maps = mc.auto_read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"));

   if (imol_maps.size() == 2) {
      float rmsd_0 = mc.get_map_rmsd_approx(imol_maps[0]);
      float rmsd_1 = mc.get_map_rmsd_approx(imol_maps[1]);
      std::cout << "rmsds " << rmsd_0 << " " << rmsd_1 << std::endl;
      if (rmsd_0 > 0.4) // test that the FWT map is the first of the pair
         if (rmsd_1 > 0.2)
      status = 1;
   }

   return status;
}

int test_svg(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_2 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));

   mc.import_cif_dictionary("ATP.cif", imol_1);
   mc.import_cif_dictionary("ATP.cif", imol_2);
   bool dark_bg = false;
   std::string s = mc.get_svg_for_residue_type(imol_1, "ATP", dark_bg);

   if (s.length() > 0) {

      if (true) {
         std::ofstream f("ATP.svg");
         f << s;
         f.close();
      }

      status = 1;
   }
   return status;

}

int test_superpose(molecules_container_t &mc) {

   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_2 = mc.read_pdb(reference_data("3pzt.pdb"));

   unsigned int n_pre = mc.get_number_of_molecules();

   coot::atom_spec_t atom_spec_1("A", 227, "", " CA ","");
   coot::atom_spec_t atom_spec_2("B", 256, "", " CA ","");
   mmdb::Atom *at_1 = mc.get_atom(imol_1, atom_spec_1);
   mmdb::Atom *at_2 = mc.get_atom(imol_2, atom_spec_2);

   coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
   coot::Cartesian atom_pos_2 = atom_to_cartesian(at_2);

   double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
   double d1 = std::sqrt(dd);
   std::cout << "test d1 " << d1 << std::endl;

   // std::pair<std::string, std::string> ss_result_pair = mc.SSM_superpose(imol_1, "A", imol_2, "B");
   superpose_results_t ss_results = mc.SSM_superpose(imol_1, "A", imol_2, "B");

   std::cout << "ss_result:\n" << ss_results.suppose_info << std::endl;
   std::cout << "ss_result:\n" << ss_results.alignment.first  << std::endl;
   std::cout << "ss_result:\n" << ss_results.alignment.second << std::endl;

   if (true) {
      for (const auto &chain : ss_results.alignment_info.cviv) {
         for (const auto &res : chain.rviv) {
            std::cout << res.residue_spec << " " << res.function_value << std::endl;
         }
      }
   }

   mc.write_coordinates(imol_2, "superposed.pdb");

   mmdb::Atom *at_3 = mc.get_atom(imol_2, atom_spec_2);
   mmdb::Atom *at_4 = mc.get_atom(imol_1, atom_spec_1);
   coot::Cartesian atom_pos_3 = atom_to_cartesian(at_3);
   coot::Cartesian atom_pos_4 = atom_to_cartesian(at_4);

   dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_3);
   double d2 = std::sqrt(dd);
   std::cout << "test d2 " << d2 << std::endl;

   unsigned int n_post = mc.get_number_of_molecules();
   if (n_pre == n_post)
      if (d1 > 50.0)
         if (d2 < 1.0)
            status = 1;

   std::cout << "debug:: n_mol_pre " << n_pre << " n_mol_post " << n_post << std::endl;
   std::cout << "debug:: atom_pos_1 " << atom_pos_1 << " atom_pos_2 " << atom_pos_2
             << " atom_pos_3 " << atom_pos_3 << " atom_pos_4 " << atom_pos_4 << std::endl;

   return status;
}

int test_non_drawn_atoms(molecules_container_t &mc) {

   auto atom_in_mesh = [] (const coot::instanced_mesh_t &mesh, const glm::vec3 &ca_pos) {
      bool status = false;
      for (unsigned int igeom=0; igeom<mesh.geom.size(); igeom++) {
         const coot::instanced_geometry_t &ig = mesh.geom[0];
         for (unsigned int j=0; j<ig.instancing_data_A.size(); j++) {
            if (glm::distance(ig.instancing_data_A[j].position, ca_pos) < 0.01) return true;
         }
      }
      return status;
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   glm::vec3 ca_pos(70.328, 51.541, 38.560);
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   auto mesh_1 = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, 1);
   mc.add_to_non_drawn_bonds(imol, "//A/270");
   auto mesh_2 = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, 1);

   // the first one should have the atom, the second should not.
   bool f1 = atom_in_mesh(mesh_1, ca_pos);
   bool f2 = atom_in_mesh(mesh_2, ca_pos);
   std::cout << "the f1 and f2 " << f1 << " " << f2 << std::endl;
   if (f1 == true)
      if (f2 == false)
         status = 1;

   return status;
}

int test_rigid_body_fit(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   //! Rigid-body fitting
   //!
   //! `multi_cids" is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||/B/56-66"
   std::string multi_cids = "//A/20-22||//A/60-66";
   status = mc.rigid_body_fit(imol, multi_cids, imol_map);

   if (status == 1) {
      int imol_lig = mc.read_pdb("misplaced-moorhen-tutorial-1-ligand.pdb");
      coot::validation_information_t vi_0 = mc.density_correlation_analysis(imol_lig, imol_map);
      status = mc.rigid_body_fit(imol_lig, "//A/301", imol_map);
      coot::validation_information_t vi_1 = mc.density_correlation_analysis(imol_lig, imol_map);
      coot::stats::single s_0 = vi_0.get_stats();
      coot::stats::single s_1 = vi_1.get_stats();
      // std::cout << "orig:      " << std::fixed << s_0.mean() << " " << std::fixed << std::sqrt(s_0.variance()) << std::endl;
      // std::cout << "post-fit:  " << std::fixed << s_1.mean() << " " << std::fixed << std::sqrt(s_1.variance()) << std::endl;
      if (s_0.mean() < 0.6)
         if (s_1.mean() > 0.8)
            status = 1;
   }

   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   return status;
}

int test_symmetry(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::Cartesian pos(1,0,10);
   coot::symmetry_info_t si = mc.get_symmetry(imol, 10.0, pos.x(), pos.y(), pos.z());
   std::vector<std::pair<symm_trans_t, Cell_Translation> > v = si.symm_trans;

   std::cout << "Got " << v.size() << " symmetry-related molecules" << std::endl;

   if (! v.empty()) // not a good test
      status = 1;

   std::cout << "cell: " << si.cell.a << " " << si.cell.b << " " << si.cell.c << " "
             << si.cell.alpha << " " << si.cell.beta << " " << si.cell.gamma
             << std::endl;

   for (unsigned int i=0; i<v.size(); i++) {
      mmdb::mat44 my_matt;
      const symm_trans_t &st     = v[i].first;
      const Cell_Translation &ct = v[i].second;
      int err = mc.get_mol(imol)->GetTMatrix(my_matt, st.isym(), st.x(), st.y(), st.z());
      std::cout << "  " << i << " " << st << " " << ct << std::endl;
      std::cout
         << "   " << st.mat[0][0] << " " << st.mat[0][1] << " " << st.mat[0][2] << " " << st.mat[0][3] << "\n"
         << "   " << st.mat[1][0] << " " << st.mat[1][1] << " " << st.mat[1][2] << " " << st.mat[1][3] << "\n"
         << "   " << st.mat[2][0] << " " << st.mat[2][1] << " " << st.mat[2][2] << " " << st.mat[2][3] << "\n"
         << "   " << st.mat[3][0] << " " << st.mat[3][1] << " " << st.mat[3][2] << " " << st.mat[3][3] << "\n"
         << std::endl;
   }

   return status;
}

int test_read_file(molecules_container_t &mc) {

   int status = 0;
   std::string s = mc.file_name_to_string("1x8b.pdb");
   if (s.size() == 225099)
      status = 1;

   return status;

}

int test_set_rotamer(molecules_container_t &mc) {

   // Actually this is a test for change to next rotamer

   starting_test(__FUNCTION__);
   int status = 0;
   int status_1 = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec_CB("A", 270, "", " CB ","");
   coot::atom_spec_t atom_spec_CG("A", 270, "", " CG ","");
   mmdb::Atom *at_start_CB = mc.get_atom(imol, atom_spec_CB);
   mmdb::Atom *at_start_CG = mc.get_atom(imol, atom_spec_CG);
   if (at_start_CG) {
      // std::cout << "debug:: " << __FUNCTION__ << " got at_start_CG " << std::endl;
      if (at_start_CB) {
         // std::cout << "debug:: " << __FUNCTION__ << " got at_start_CB " << std::endl;
         coot::Cartesian atom_pos_start_CB = atom_to_cartesian(at_start_CB);
         coot::Cartesian atom_pos_start_CG = atom_to_cartesian(at_start_CG);
         mc.change_to_next_rotamer(imol, "//A/270", "");
         coot::Cartesian atom_pos_done_CB = atom_to_cartesian(at_start_CB);
         coot::Cartesian atom_pos_done_CG = atom_to_cartesian(at_start_CG);
         double dd_CB = coot::Cartesian::lengthsq(atom_pos_start_CB, atom_pos_done_CB);
         double dd_OG = coot::Cartesian::lengthsq(atom_pos_start_CG, atom_pos_done_CG);
         double d_CB = std::sqrt(dd_CB);
         double d_CG = std::sqrt(dd_OG);
         std::cout << "d_CB " << d_CB << " d_CG " << d_CG << std::endl;
         if (d_CG > 0.3)
            if (d_CB < 0.1)
               status_1 = 1;

         if (status_1) {

            int status_2 = 0;
            mc.change_to_previous_rotamer(imol, "//A/270", "");
            mc.change_to_previous_rotamer(imol, "//A/270", "");
            coot::Cartesian atom_pos_done_CB = atom_to_cartesian(at_start_CB);
            coot::Cartesian atom_pos_done_CG = atom_to_cartesian(at_start_CG);
            double dd_CB = coot::Cartesian::lengthsq(atom_pos_start_CB, atom_pos_done_CB);
            double dd_OG = coot::Cartesian::lengthsq(atom_pos_start_CG, atom_pos_done_CG);
            double d_CB = std::sqrt(dd_CB);
            double d_CG = std::sqrt(dd_OG);
            std::cout << "d_CB " << d_CB << " d_CG " << d_CG << std::endl;
            if (d_CG > 0.3)
               if (d_CB < 0.1)
                  status_2 = 1;

            if (status_2 == 1) status = 1; // all gooe

         }
      }
   }

   return status;
}

int test_replace_model_from_file(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         std::string rn_1 = at_1->GetResidue()->GetResName();
         mc.replace_molecule_by_model_from_file(imol, "3pzt.pdb");
         mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec);
         if (at_2) {
            std::string rn_2 = at_2->GetResidue()->GetResName();

            if (rn_1 == "SER")
               if (rn_2 == "PHE")
                  status = 1;
         }
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_user_defined_bond_colours(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 1, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         std::cout << "...................... here A " << std::endl;
         std::map<unsigned int, std::array<float, 3> > colour_map;
         colour_map[0] = std::array<float, 3> {0.42222222, 0.7, 0.4};
         colour_map[2] = std::array<float, 3> {0.42222222, 0.4, 0.7};
         colour_map[1] = std::array<float, 3> {0.7, 0.4, 0.42222222};
         std::vector<std::pair<std::string, unsigned int> > indexed_residues_cids;
         indexed_residues_cids.push_back(std::make_pair("//A",2));
         indexed_residues_cids.push_back(std::make_pair("//A/100-200",1));
         indexed_residues_cids.push_back(std::make_pair("//A/130-150",0));
         std::string mode("USER-DEFINED-COLOURS");
         mc.set_user_defined_bond_colours(imol, colour_map);
         mc.set_user_defined_atom_colour_by_residue(imol, indexed_residues_cids);
         coot::instanced_mesh_t im = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, 1);
         std::cout << "...................... here B " << im.geom.size() << std::endl;
         if (im.geom.size() > 3) {
            if (im.geom[0].instancing_data_A.size() > 1000)
               status = 1;
            const auto& g = im.geom[1].instancing_data_B;
            std::cout << "........... g.size() " << g.size() << std::endl;
            for (unsigned int j=0; j<g.size(); j++) {
               if (j > 300) continue;
               const auto &d = g[j];
               std::cout << j << " " << glm::to_string(d.colour) << std::endl;
            }
         }
      }
   }
   return status;
}

int test_replace_map(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   mc.replace_map_by_mtz_from_file(imol_map, reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false);
   status = 1;
   return status;
}

int test_residue_name_group(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   mc.get_monomer("BGC");

   std::string g1 = mc.get_group_for_monomer("PHE");
   std::string g2 = mc.get_group_for_monomer("BGC");

   std::cout << "g1: " << g1 << std::endl;
   std::cout << "g2: " << g2 << std::endl;

   if (g1 == "peptide")
      if (g2 == "pyranose")
         status = 1;

   return status;
}

int test_alt_conf_and_rotamer(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 131, "", " O  ", "A");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         mc.imol_refinement_map = imol_map;
         mc.refine_residues_using_atom_cid(imol, "//A/131", "TRIPLE");
         mc.write_coordinates(imol, "alt-conf-and-rotamer-and-refine.pdb");
         status = 1;
      }
   }
   return status;
}

int test_alt_conf_and_rotamer_v2(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("mol-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   std::cout << "........... test_alt_conf_and_rotamer_v2() imol_map " << imol_map << std::endl;

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec_A("A", 131, "", " NE2","A");
      coot::atom_spec_t atom_spec_B("A", 131, "", " NE2","B");
      mmdb::Atom *at_A = mc.get_atom(imol, atom_spec_A);
      mmdb::Atom *at_B = mc.get_atom(imol, atom_spec_B);
      if (at_A) {
         coot::Cartesian atom_pos_start_A = atom_to_cartesian(at_A);
         coot::Cartesian atom_pos_start_B = atom_to_cartesian(at_B);
         mc.imol_refinement_map = imol_map;
         mc.refine_residues_using_atom_cid(imol, "//A/131", "TRIPLE");
         mc.change_to_next_rotamer(imol, "//A/131", "A");
         mc.change_to_next_rotamer(imol, "//A/131", "A");
         mc.change_to_next_rotamer(imol, "//A/131", "A");
         coot::Cartesian atom_pos_done_A = atom_to_cartesian(at_A);
         coot::Cartesian atom_pos_done_B = atom_to_cartesian(at_B);
         mc.write_coordinates(imol, "alt-conf-and-rotamer.pdb");
         double dd_1 = coot::Cartesian::lengthsq(atom_pos_done_A, atom_pos_done_B);
         std::cout << "dd_1 " << dd_1 << std::endl;
         if (dd_1 > 9.0)
            status = 1;
      } else {
         std::cout << "In test_alt_conf_and_rotamer_v2() No atom found " << atom_spec_A << std::endl;
      }

      
   }

   std::cout << "done test_alt_conf_and_rotamer_v2()" << std::endl;
   return status;
}


int test_moorhen_h_bonds(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   mc.add_hydrogen_atoms(imol); // no hydrogen bonds found without hydrogens in the model
   const std::string &cid_str = "//A/270";
   std::vector<moorhen::h_bond> h_bonds = mc.get_h_bonds(imol, cid_str);

   std::cout << "INFO:: in test_moorhen_h_bonds() we got " << h_bonds.size() << " H-bonds" << std::endl;
   for (unsigned int i=0; i<h_bonds.size(); i++) {
      std::cout << " "
                << h_bonds[i].donor.chain    << " " << h_bonds[i].donor.res_no    << " " << h_bonds[i].donor.name
                << h_bonds[i].acceptor.chain << " " << h_bonds[i].acceptor.res_no << " " << h_bonds[i].acceptor.name
                << std::endl;
   }

   // test one.
   if (h_bonds.size() == 3) {
      if (h_bonds[2].acceptor.name == " OD1")
         if (h_bonds[2].donor.name == " N  ")
            status = 1;
   }
   return status;
}

int test_template(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

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

   molecules_container_t mc(false); // quiet

   mc.fill_rotamer_probability_tables();

   if (! last_test_only) {

      status += run_test(test_new_position_for_atoms_in_residues, "new positions for atoms in residues", mc);
      status += run_test(test_transformation_for_atom_selection, "transformation for atoms",             mc);
      status += run_test(test_copy_fragment_using_residue_range, "copy-fragment using residue range",    mc);
      status += run_test(test_density_correlation_validation, "density correlation validation",          mc);
      status += run_test(test_pepflips_using_difference_map, "Pepflips from Difference Map",             mc);
      status += run_test(test_difference_map_contours, "difference map density mesh", mc);
      status += run_test(test_rsr_using_residue_range, "rsr using residue range", mc);
      status += run_test(test_copy_fragment_using_cid, "copy-fragment using cid", mc);
      status += run_test(test_no_dictionary_residues,  "no-dictionary residues", mc);
      status += run_test(test_rsr_using_atom_cid,    "rsr using atom cid",       mc);
      status += run_test(test_auto_fit_rotamer_1,    "auto-fit rotamer",         mc);
      status += run_test(test_auto_fit_rotamer_2,    "auto-fit rotamer t2",      mc);
      status += run_test(test_rota_dodecs_mesh,      "rotamer dodecahedra mesh", mc);
      status += run_test(test_delete_molecule,       "delete_moelcule",          mc);
      status += run_test(test_rama_balls_mesh,       "rama balls mesh",          mc);
      status += run_test(test_density_mesh,          "density mesh",             mc);
      status += run_test(test_updating_maps,         "updating maps",            mc);
      status += run_test(test_delete_residue,        "delete residue",           mc);
      status += run_test(test_delete_chain,          "delete chain",             mc);
      status += run_test(test_delete_atom,           "delete atom",              mc);
      status += run_test(test_jiggle_fit,            "Jiggle-fit",               mc);
      status += run_test(test_pepflips,              "pepflips",                 mc);
      status += run_test(test_mutate,                "mutate",                   mc);
      status += run_test(test_rsr,                   "rsr",                      mc);
      status += run_test(test_jed_flip,              "JED Flip",                 mc);
      status += run_test(test_add_water,             "add waters",               mc);
      status += run_test(test_cis_trans,             "cis_trans conversion",     mc);
      status += run_test(test_bonds_mesh,            "bonds mesh",               mc);
      status += run_test(test_eigen_flip,            "Eigen Flip",               mc);
      status += run_test(test_read_a_map,            "read a map",               mc);
      status += run_test(test_add_compound,          "add compound",             mc);
      status += run_test(test_weird_delete,          "delete II",                mc);
      status += run_test(test_fill_partial,    "fill partially-filled residues", mc);
      status += run_test(test_add_alt_conf,          "add alt conf",             mc);
      status += run_test(test_delete_literal,        "delete literal",           mc);
      status += run_test(test_side_chain_180,        "side-chain 180",           mc);
      status += run_test(test_peptide_omega,         "peptide omega",            mc);
      status += run_test(test_undo_and_redo,         "undo and redo",            mc);
      status += run_test(test_undo_and_redo_2,       "undo/redo 2",              mc);
      status += run_test(test_merge_molecules,       "merge molecules",          mc);
      status += run_test(test_dictionary_bonds,      "dictionary bonds",         mc);
      status += run_test(test_replace_fragment,      "replace fragment",         mc);
      status += run_test(test_gaussian_surface,      "Gaussian surface",         mc);
      status += run_test(test_missing_atoms_info,    "missing atom info",        mc);
      status += run_test(test_move_molecule_here,    "move_molecule_here",       mc);
      status += run_test(test_sequence_generator,    "Make a sequence string",   mc);
      status += run_test(test_rotamer_validation,    "rotamer validation",       mc);
      status += run_test(test_ligand_fitting_here,   "Ligand fitting here",      mc);
      status += run_test(test_ligand_contact_dots,   "lgiand contact dots",      mc);
      status += run_test(test_rama_validation,       "rama validation 2",        mc); // for the plot, not the graph
      status += run_test(test_ramachandran_analysis, "ramachandran analysis",    mc); // for the graph, not the plot
      status += run_test(test_difference_map_peaks,  "Difference Map Peaks",     mc);
      status += run_test(test_non_standard_residues, "non-standard residues",    mc);
      status += run_test(test_import_cif_dictionary, "import cif dictionary",    mc);
      status += run_test(test_add_terminal_residue,  "add terminal residue",     mc);
      status += run_test(test_instanced_rota_markup, "instanced rotamer mesh",   mc);
      status += run_test(test_new_position_for_atoms,"new positions for atoms",  mc);
      status += run_test(test_molecular_representation, "molecular representation mesh", mc);
   }

   // status += run_test(test_undo_and_redo, "undo and redo", mc);

   // status += run_test(test_alt_conf_and_rotamer,            "Alt Conf then rotamer", mc);

   // status += run_test(test_rigid_body_fit, "rigid-body fit", mc);

   // status += run_test(test_jiggle_fit,            "Jiggle-fit",               mc);

   // status += run_test(test_editing_session_tutorial_1, "an Tutorial 1 editing session",         mc);

   // status += run_test(test_broken_function, "Something was broken",         mc);

   // status += run_test(test_molecular_representation, "molecular representation mesh", mc);

   // status += run_test(test_delete_side_chain, "delete side chain", mc);

   // status += run_test(test_colour_rules, "colour rules", mc);

   // status += run_test(test_mmrrcc, "MMRRCC", mc);

   //status += run_test(test_auto_read_mtz, "auto-read MTZ", mc);

   // status += run_test(test_instanced_bonds_mesh, "insta bonds mesh", mc);

   // status = run_test(test_utils, "utils", mc);

   // status = run_test(test_instanced_bonds_mesh, "instanced_bonds", mc);

   // status = run_test(test_svg, "svg string", mc);

   // status = run_test(test_superpose, "SSM superpose ", mc);

   // status = run_test(test_multi_colour_rules, "multi colour rules ", mc);

   // status = run_test(test_non_drawn_atoms, "non-drawn atoms", mc);

   // status = run_test(test_add_terminal_residue, "add terminal residue", mc);

   // status = run_test(test_symmetry, "symmetry", mc);

   // status += run_test(test_add_hydrogen_atoms, "add hydrogen atoms", mc);

   // status = run_test(test_set_rotamer, "set rotamer ", mc);

   // status = run_test(test_alt_conf_and_rotamer_v2, "alt-conf and rotamer v2 ", mc);

   status = run_test(test_moorhen_h_bonds, "moorhen H-bonds ", mc);

   // Note to self:
   //
   // change the autofit_rotamer test so that it tests the change of positions of the atoms of the neighboring residues.

   // status = run_test(test_replace_model_from_file, "replace model from file", mc);

   // status = run_test(test_user_defined_bond_colours, "user-defined bond colours", mc);

   // status = run_test(test_replace_map, "replace map from mtz", mc);

   // status = run_test(test_residue_name_group, "residue name group", mc);


   int all_tests_status = 1; // fail!
   if (status == n_tests) all_tests_status = 0;

   return all_tests_status;

}
