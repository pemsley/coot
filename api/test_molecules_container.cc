
#include <iostream>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include "molecules_container.hh"
#include "filo-tests.hh"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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

void colour_analysis(const coot::simple_mesh_t &mesh) {

   auto is_near_colour = [] (const glm::vec4 &col_1, const glm::vec4 &col_2) {
      float cf = 0.04;
      if (std::fabs(col_2.r - col_1.r) < cf)
         if (std::fabs(col_2.g - col_1.g) < cf)
            if (std::fabs(col_2.b - col_1.b) < cf)
               if (std::fabs(col_2.a - col_1.a) < cf)
                  return true;
      return false;
   };

   auto sorter = [] (const std::pair<glm::vec4, unsigned int> &p1,
                     const std::pair<glm::vec4, unsigned int> &p2) {
      if (p1.first[0] == p2.first[0]) {
         return (p1.first[1] > p2.first[1]);
      } else {
         return (p1.first[0] > p2.first[0]);
      }
   };

   std::vector<std::pair<glm::vec4, unsigned int> > colour_count;
   for (unsigned int i=0; i<mesh.vertices.size(); i++) {
      const auto &vertex = mesh.vertices[i];
      const glm::vec4 &col = vertex.color;
      bool found_col = false;
      for (unsigned int j=0; j<colour_count.size(); j++) {
         if (is_near_colour(col, colour_count[j].first)) {
            colour_count[j].second ++;
            found_col = true;
            break;
         }
      }
      if (! found_col) {
         colour_count.push_back(std::make_pair(col, 1));
      }
   }

   std::sort(colour_count.begin(), colour_count.end(), sorter);

   std::cout << "INFO:: " << colour_count.size() << " colours" << std::endl;
   for (unsigned int i=0; i<colour_count.size(); i++)
      std::cout << "    " << glm::to_string(colour_count[i].first) << " "
                << std::setw(7) << std::right << colour_count[i].second << std::endl;

}

void colour_analysis(const coot::instanced_mesh_t &mesh) {

   auto is_near_colour = [] (const glm::vec4 &col_1, const glm::vec4 &col_2) {
      float cf = 0.04;
      if (std::fabs(col_2.r - col_1.r) < cf)
         if (std::fabs(col_2.g - col_1.g) < cf)
            if (std::fabs(col_2.b - col_1.b) < cf)
               if (std::fabs(col_2.a - col_1.a) < cf)
                  return true;
      return false;
   };

   auto sorter = [] (const std::pair<glm::vec4, unsigned int> &p1,
                     const std::pair<glm::vec4, unsigned int> &p2) {
      if (p1.first[0] == p2.first[0]) {
         return (p1.first[1] > p2.first[1]);
      } else {
         return (p1.first[0] > p2.first[0]);
      }
   };

   std::vector<std::pair<glm::vec4, unsigned int> > colour_count;

   for (unsigned int i=0; i<mesh.geom.size(); i++) {
      const coot::instanced_geometry_t &ig = mesh.geom[i];
      for (unsigned int jj=0; jj<ig.instancing_data_A.size(); jj++) {
         const auto &col =  ig.instancing_data_A[jj].colour;
         bool found_col = false;
         for (unsigned int j=0; j<colour_count.size(); j++) {
            if (is_near_colour(col, colour_count[j].first)) {
               colour_count[j].second ++;
               found_col = true;
               break;
            }
         }
         if (! found_col) {
            colour_count.push_back(std::make_pair(col, 1));
         }
      }

      for (unsigned int jj=0; jj<ig.instancing_data_B.size(); jj++) {
         const auto &col =  ig.instancing_data_B[jj].colour;
         bool found_col = false;
         for (unsigned int j=0; j<colour_count.size(); j++) {
            if (is_near_colour(col, colour_count[j].first)) {
               colour_count[j].second ++;
               found_col = true;
               break;
            }
         }
         if (! found_col) {
            colour_count.push_back(std::make_pair(col, 1));
         }
      }

   }


   for (unsigned int i=0; i<mesh.markup.vertices.size(); i++) {
      const auto &vertex = mesh.markup.vertices[i];
      const glm::vec4 &col = vertex.color;
      bool found_col = false;
      for (unsigned int j=0; j<colour_count.size(); j++) {
         if (is_near_colour(col, colour_count[j].first)) {
            colour_count[j].second ++;
            found_col = true;
            break;
         }
      }
      if (! found_col) {
         colour_count.push_back(std::make_pair(col, 1));
      }
   }

   std::sort(colour_count.begin(), colour_count.end(), sorter);

   std::cout << "INFO:: " << colour_count.size() << " colours" << std::endl;
   for (unsigned int i=0; i<colour_count.size(); i++)
      std::cout << "    " << glm::to_string(colour_count[i].first) << " "
                << std::setw(7) << std::right << colour_count[i].second << std::endl;


}

int test_auto_fit_rotamer_1(molecules_container_t &mc_in) {

   starting_test(__FUNCTION__);
   int status = 0; // initially fail status

   molecules_container_t mc;
   mc.geometry_init_standard();
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {

         coot::residue_spec_t res_spec("A", 61, "");
         mmdb::Residue *r = coot::util::get_residue(res_spec, mc[imol].atom_sel.mol);
         if (r) {
            mmdb::Atom *cz = r->GetAtom(" CZ ");
            if (cz) {
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
               std::cout << "in test_auto_fit_rotamer_1() CZ atom not found " << std::endl;
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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);

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

   int imol          = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map      = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT",    "PHWT",    "W", false, false);
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "DELFWT", "PHDELWT", "W", false, true);
   mc.associate_data_mtz_file_with_map(imol_map, reference_data("moorhen-tutorial-map-number-4.mtz"), "F", "SIGF", "FREER");

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
   } else {
      std::cout << "in test_updating_maps() atom spec was null" << std::endl;
   }

   // now update the maps
   mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map, imol_map);
   float new_rail_points = mc.calculate_new_rail_points();
   float gpt = mc.rail_points_total();
   std::cout << "###### RailPoints gained: " << new_rail_points << " rail points total " << gpt << std::endl;

   std::vector<std::pair<clipper::Coord_orth, float> > ddmp = mc.get_diff_diff_map_peaks(imol_diff_map, 70, 50, 30);

   std::cout << "test_updating_maps(): We got " << ddmp.size() << " difference map peaks" << std::endl;

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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);

   coot::atom_spec_t atom_spec("A", 61, "", " CZ ", "");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
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
   } else {
      std::cout << "in test_undo_and_redo_2() failed to find atom " << std::endl;
   }
   return status;
}


int test_ramachandran_analysis(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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

   std::string coords_fn = reference_data("moorhen-tutorial-structure-number-4.pdb");
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
   float radius = 22;
   float contour_level = 0.13;
   mc.set_map_is_contoured_with_thread_pool(true);
   coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);

   // std::cout << "DEBUG:: test_density_mesh(): " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
   // << " triangles" << std::endl;

   unsigned int size_1 = map_mesh.vertices.size();
   if (map_mesh.vertices.size() > 30000)
      status = 1;

   mc.set_map_is_contoured_with_thread_pool(false);
   coot::simple_mesh_t map_mesh_2 = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
   unsigned int size_2 = map_mesh_2.vertices.size();

   std::cout << "compare sizes " << size_1 << " " << size_2 << std::endl;

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

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {
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

         int n_cycles = 1000;
         std::string mode = "SPHERE";
         mc.refine_residues(imol, "A", 14, "", "", mode, n_cycles);
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
      }
   }

   return status;
}

int test_rsr_using_atom_cid(molecules_container_t &mc_in) {

   starting_test(__FUNCTION__);
   int status = 0;

   molecules_container_t mc;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   std::cout << "imol: " << imol << std::endl;
   std::cout << "imol_map: " << imol_map << std::endl;

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {

         std::cout << "calling set_imol_refinement_map with imol_map " << imol_map << std::endl;
         mc.set_imol_refinement_map(imol_map);
         std::cout << "debug:: A imol_refinement_map is now " << mc.imol_refinement_map << std::endl;

         std::string chain_id = "A";
         int res_no = 14; // this residue is problematic in moorhen-tutorial-structure-number-1.pdb
         std::string ins_code;

         // std::string cid = "//A/14/CA";
         std::string cid = "//A/187/CA";

         coot::atom_spec_t atom_spec_N_1(chain_id, res_no, ins_code, " N  ","");
         coot::atom_spec_t atom_spec_N_2(chain_id, res_no, ins_code, " CG2",""); // not a nitrogen atom
         mmdb::Atom *at_n_1 = mc.get_atom(imol, atom_spec_N_1);
         mmdb::Atom *at_n_2 = mc.get_atom(imol, atom_spec_N_2);
         coot::Cartesian pt_n_1_pre = atom_to_cartesian(at_n_1);
         coot::Cartesian pt_n_2_pre = atom_to_cartesian(at_n_2);

         int n_cycles = 1000;
         std::string mode = "SPHERE";
         std::cout << "debug:: B imol_refinement_map is now " << mc.imol_refinement_map << std::endl;
         float f = mc.get_map_weight();
         std::cout << "debug:: map weight " << f << std::endl;

         mc.add_to_non_drawn_bonds(imol, cid);
         mc.refine_residues_using_atom_cid(imol, cid, mode, n_cycles);
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

         mc.clear_non_drawn_bonds(imol);
      }
   }
   return status;
}

int test_rsr_using_residue_range(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
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
   // mc.set_map_weight(w * 100.0);
   int n_cycles = 500;
   mc.refine_residue_range(imol, "A", 131, 136, n_cycles);
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
   mc.write_coordinates(imol, "post-refine-using-residue-range.pdb");
   return status;
}

int test_rsr_using_multi_atom_cid(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   // this can be a coot::utils function
   //
   auto get_n_diffs = [] (mmdb::Manager *mol, mmdb::Manager *mol_orig) {

      int n_checked = 0;
      int n_diffs = 0;
      int imod = 1;
      mmdb::Model *model_p_1 = mol->GetModel(imod);
      if (model_p_1) {
         mmdb::Model *model_p_2 = mol_orig->GetModel(imod);
         if (model_p_2) {
            int n_chains_1 = model_p_1->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains_1; ichain++) {
               mmdb::Chain *chain_p_1 = model_p_1->GetChain(ichain);
               int n_chains_2 = model_p_2->GetNumberOfChains();
               for (int jchain=0; jchain<n_chains_2; jchain++) {
                  mmdb::Chain *chain_p_2 = model_p_2->GetChain(jchain);
                  int n_res_1 = chain_p_1->GetNumberOfResidues();
                  int n_res_2 = chain_p_2->GetNumberOfResidues();

                  for (int ires=0; ires<n_res_1; ires++) {
                     mmdb::Residue *residue_p_1 = chain_p_1->GetResidue(ires);
                     if (residue_p_1) {
                        int seqnum_1 = residue_p_1->GetSeqNum();
                        std::string ins_code_1 = residue_p_1->GetInsCode();

                        for (int jres=0; jres<n_res_2; jres++) {
                           mmdb::Residue *residue_p_2 = chain_p_2->GetResidue(jres);
                           if (residue_p_2) {
                              int seqnum_2 = residue_p_2->GetSeqNum();
                              std::string ins_code_2 = residue_p_2->GetInsCode();

                              if (seqnum_1 ==  seqnum_2) {
                                 if (ins_code_1 == ins_code_2) {
                                    std::string rn_1 = residue_p_1->GetResName();
                                    std::string rn_2 = residue_p_2->GetResName();
                                    if (rn_1 == rn_2) {

                                       int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
                                       for (int iat=0; iat<n_atoms_1; iat++) {
                                          mmdb::Atom *at_1 = residue_p_1->GetAtom(iat);
                                          if (! at_1->isTer()) {
                                             std::string atom_name_1(at_1->GetAtomName());
                                             std::string alt_conf_1(at_1->altLoc);

                                             int n_atoms_2 = residue_p_2->GetNumberOfAtoms();
                                             for (int jat=0; jat<n_atoms_2; jat++) {
                                                mmdb::Atom *at_2 = residue_p_2->GetAtom(jat);
                                                if (! at_2->isTer()) {
                                                   std::string atom_name_2(at_2->GetAtomName());
                                                   std::string alt_conf_2(at_2->altLoc);
                                                   if (atom_name_1 == atom_name_2) {
                                                      if (alt_conf_1 == alt_conf_2) {
                                                         n_checked++;
                                                         float dx = at_1->x - at_2->x;
                                                         float dy = at_1->y - at_2->y;
                                                         float dz = at_1->z - at_2->z;
                                                         if ((fabsf(dx) + fabsf(dy) + fabsf(dz)) > 0.01)
                                                            n_diffs++;
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      // std::cout << "n_diffs: " << n_diffs << " n_checked: " << n_checked << std::endl;
      return n_diffs;
   };

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
   mc.set_imol_refinement_map(imol_map);
   int n_cycles = 1000;
   std::string mode = "SINGLE";
   std::string multi_cid = "//A/10-12";
   mmdb::Manager *mol = mc.get_mol(imol);
   mmdb::Manager *mol_orig = coot::util::copy_molecule(mol);
   mc.refine_residues_using_atom_cid(imol, multi_cid, mode, n_cycles);
   int n_diffs_1 = get_n_diffs(mol, mol_orig);
   multi_cid = "//A/10-12||//A/20-22";
   mc.refine_residues_using_atom_cid(imol, multi_cid, mode, n_cycles);
   int n_diffs_2 = get_n_diffs(mol, mol_orig);
   if (n_diffs_2 > (n_diffs_1 + 10)) status = 1;
   std::cout << "n_diffs_1 " << n_diffs_1 << " n_diffs_2 " << n_diffs_2 << std::endl;

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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
   mc.set_imol_refinement_map(imol_map);

   // test adding to the N-terminus
   bool part_one_done = false;
   coot::atom_spec_t atom_spec_in_new_residue("A", 284, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec_in_new_residue);
   std::cout << "Here 1 A" << std::endl;
   if (! at_1) { // it's not there to begin with
      std::cout << "Here 1 B" << std::endl;
      mc.add_terminal_residue_directly(imol, "A", 285, "");
      mc.write_coordinates(imol, "test-add-terminal-residue-with-added-terminal-residue.pdb");
      mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec_in_new_residue);
      if (at_2) {
         std::cout << "Here 1 C" << std::endl;
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

   std::cout << "part_one_done "   << part_one_done   << std::endl;
   std::cout << "part_two_done "   << part_two_done   << std::endl;
   std::cout << "part_three_done " << part_three_done << std::endl;
   std::cout << "part_four_done "  << part_four_done  << std::endl;

   if (part_one_done && part_two_done && part_three_done && part_four_done)
      status = 1;

   mc.close_molecule(imol);
   mc.close_molecule(imol_map);

   return status;
}

// this was converted from Filo's test
int test_add_terminal_residue_v2(molecules_container_t &molecules_container) {

   int status = 0;
   int coordMolNo = molecules_container.read_pdb("./5a3h.mmcif");
   int mapMolNo = molecules_container.read_mtz("./5a3h_sigmaa.mtz", "FWT", "PHWT", "", false, false);
   molecules_container.set_imol_refinement_map(mapMolNo);

   int atom_count_1 = molecules_container.get_number_of_atoms(coordMolNo);
   molecules_container.delete_using_cid(coordMolNo, "A/4-100/*", "LITERAL");
   molecules_container.delete_using_cid(coordMolNo, "A/105-200/*", "LITERAL");
   molecules_container.delete_using_cid(coordMolNo, "A/205-303/*", "LITERAL");
   std::string mmcifString_1 = molecules_container.molecule_to_mmCIF_string(coordMolNo);
   int atom_count_2 = molecules_container.get_number_of_atoms(coordMolNo);
   molecules_container.write_coordinates(coordMolNo, "pre-addition.cif");
   molecules_container.write_coordinates(coordMolNo, "pre-addition-2.cif");

   int result = molecules_container.add_terminal_residue_directly_using_cid(coordMolNo, "A/104");
   if (result == -1) {
       std::cout << "test_add_terminal_residue_v2() fail: Result is -1" << std::endl;
   }
   int atom_count_3 = molecules_container.get_number_of_atoms(coordMolNo);
   std::string mmcifString_2 = molecules_container.molecule_to_mmCIF_string(coordMolNo);
   if (mmcifString_1 != mmcifString_2) {
       std::cout << "Error: MMCIF Strings are not equal" << std::endl;
   }

   molecules_container.write_coordinates(coordMolNo, "post-addition.cif");
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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   coot::atom_spec_t atom_spec("A", 151, "", " N  ", "");
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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

int test_copy_fragment_for_refinement_using_cid(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   if (at_1) {
      std::string cid = "//A/131-140";
      int imol_new = mc.copy_fragment_for_refinement_using_cid(imol, cid);
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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
         if (d < 0.001) {
            std::string fn = mc.get_cif_file_name("ATP", coot::protein_geometry::IMOL_ENC_ANY);
            if (fn == "ATP.cif")
               status = 1;
         }
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "DELFWT", "PHDELWT", "W", false, true);
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

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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

               // std::cout << "correl " << res.residue_spec.res_no << " " << res.function_value << std::endl;
            }
         }
         std::cout << "debug:: in test_density_correlation_validation n_res: " << n_res << std::endl;
         if (bad_correls == false) {
            if (n_res > 250)
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

int test_read_a_missing_map(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map_1 = mc.read_ccp4_map("a-map-that-just-isnt-there.map", false);
   if (mc.is_valid_map_molecule(imol_map_1))
      status = 0;
   else
      status = 1;
   return status;
}



int test_ligand_fitting_here(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
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
         std::cout << "debug:: test_ligand_fitting_here() failed to get model for ligand " << imol_ligand << std::endl;
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   mc.close_molecule(imol_ligand);
   return status;
}

int test_ligand_fitting_in_map(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   auto largest_eigenvalue = [] (const std::vector<double> &evs) {
      double l = 0.0;
      for (unsigned int i=0; i<evs.size(); i++)
         if (evs[i] > l) l = evs[i];
      return l;
   };

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
   int imol_ligand = mc.get_monomer("GLC");
   mc.write_coordinates(imol_ligand, "ligand.pdb");

   mc.write_map(imol_map, "number-4.map");

   mc.write_map(imol_map, "moorhen-test.map");

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_model_molecule(imol_ligand)) {
         if (mc.is_valid_map_molecule(imol_map)) {
            float n_rmsd = 1.0;
            bool make_conformers = true;
            unsigned int n_conformers = 8;
            std::vector<molecules_container_t::fit_ligand_info_t> solutions =
               mc.fit_ligand(imol, imol_map, imol_ligand, n_rmsd, make_conformers, n_conformers);
            std::cout << "found " << solutions.size() << " ligand fitting solutions" << std::endl;

            // check that these solutions have different eigen values (because they
            // are different conformers)
            std::vector<double> ligands_largest_eigenvector;
            std::vector<molecules_container_t::fit_ligand_info_t>::const_iterator it;
            for (it=solutions.begin(); it!=solutions.end(); ++it) {
               const auto &fli(*it);
               if (false) { // let's write out those solutions
                  std::string fn("Ligand-sol-" + coot::util::int_to_string(fli.imol) + ".pdb");
                  mc.write_coordinates(fli.imol, fn);
               }
               auto eigenvalues = mc.get_eigenvalues(fli.imol, "A", 1, "");
               double f = largest_eigenvalue(eigenvalues);
               ligands_largest_eigenvector.push_back(f);
            }

            if (false)
               for (unsigned int ii=0; ii<ligands_largest_eigenvector.size(); ii++)
                  std::cout << "Largest-ev: " << ii << " " << ligands_largest_eigenvector[ii] << std::endl;

            coot::stats::single ss(ligands_largest_eigenvector);
            double sd = std::sqrt(ss.variance());
            std::cout << "EV sd " << sd << std::endl;
            if (sd > 0.001)
                if (solutions.size()< 5)
                   status = 1;

         } else {
            std::cout << "Not a valid map molecule for moorhen-tutorial-map-number-4.mtz" << std::endl;
         }

      } else {
         std::cout << "Not a valid model molecule for GLC get_monomer() " << std::endl;
      }
   } else {
      std::cout << "Not a valid molecule molecule for moorhen-tutorial-structure-number-4.pdb" << std::endl;
   }

   return status;

}


int test_write_map_is_sane(molecules_container_t &mc) {


   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
   mc.write_map(imol_map, "test_for_map.map");
   struct stat buf_2;
         int istat_2 = stat("test_for_map.map", &buf_2);
         if (istat_2 == 0)
            if (buf_2.st_size > 1000000)
               status = 1;
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
         mc.fit_to_map_by_random_jiggle(imol, residue_spec, 110, -1);
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

   std::cout << "in test_jiggle_fit() A with status " << status   << std::endl;

   if (status == 1) {

      status = 0;
      std::cout << "Second jiggle-fit test: using atom selection ------------------------------" << std::endl;
      imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
      int imol_start = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
      int imol_other = mc.read_pdb(reference_data("weird-orientation-tut-4.pdb"));
      if (mc.is_valid_model_molecule(imol_other)) {
         // now test that we started with bad fit to density
         // fit!
         int imol_blur = mc.sharpen_blur_map(imol_map, 200, false);
         mc.write_map(imol_blur, "blurred.map");
         mc.imol_refinement_map = imol_blur;
         coot::validation_information_t vi_0 = mc.density_correlation_analysis(imol_start, imol_blur);
         coot::validation_information_t vi_1 = mc.density_correlation_analysis(imol_other, imol_blur);
         mc.fit_to_map_by_random_jiggle_using_cid(imol_other, "//A", 5000, 2);
         // now test that we have good fit to density.
         coot::validation_information_t vi_2 = mc.density_correlation_analysis(imol_other, imol_blur);

         coot::stats::single s_0 = vi_0.get_stats();
         coot::stats::single s_1 = vi_1.get_stats();
         coot::stats::single s_2 = vi_2.get_stats();

         // 20230402-PE These results are disappointing - they are not as good as doing it interactively.
         // I wonder what the difference is.

         std::cout << "orig:     mean " << std::fixed << std::right << std::setw(10) << s_0.mean() << " sd " << std::fixed << std::sqrt(s_0.variance()) << std::endl;
         std::cout << "pre-fit:  mean " << std::fixed << std::right << std::setw(10) << s_1.mean() << " sd " << std::fixed << std::sqrt(s_1.variance()) << std::endl;
         std::cout << "post-fit: mean " << std::fixed << std::right << std::setw(10) << s_2.mean() << " sd " << std::fixed << std::sqrt(s_2.variance()) << std::endl;

         float d1 = s_0.mean() - s_2.mean();
         float d2 = s_2.mean() - s_1.mean();

         if (d2 > d1) // we got most of the way there
            status = 1;

         mc.write_coordinates(imol_other, "jiggled.pdb");
      }
   }
   return status;
}

int test_jiggle_fit_with_blur(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4-unfit.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
   if (mc.is_valid_model_molecule(imol)) {
      mc.fit_to_map_by_random_jiggle_with_blur_using_cid(imol, imol_map, "//", 200, 2000, 3.0);
      mc.write_coordinates(imol, "jiggled-with-blur.pdb");
   }
   return status;
}

int test_jiggle_fit_params(molecules_container_t &mc) {

   auto get_score_vs_true_solution = [&mc] (int imol_testing, int imol_ref) {

      // ref is the true solution

      std::map<int, clipper::Coord_orth> ca_map;

      std::cout << "imol_testing " << imol_testing << " imol_ref " << imol_ref << std::endl;

      mmdb::Manager *mol_ref = mc.get_mol(imol_ref);
      int imod = 1;
      mmdb::Model *model_p = mol_ref->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id(chain_p->GetChainID());
            if (chain_id == "N") {
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at->isTer()) {
                           std::string atom_name(at->GetAtomName());
                           if (atom_name == " CA ") {
                              int res_no = residue_p->GetSeqNum();
                              clipper::Coord_orth co = coot::co(at);
                              ca_map[res_no] = co;
                           }
                        }
                     }
                  }
               }
            }
         }
      }


      // now find the CAs of the "hypothesis" chain
      mmdb::Manager *mol_testing = mc.get_mol(imol_testing);
      model_p = mol_testing->GetModel(imod);
      int n_found = 0;
      double sum_dist = 0.0;
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id(chain_p->GetChainID());
            if (chain_id == "E") {
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at->isTer()) {
                           std::string atom_name(at->GetAtomName());
                           if (atom_name == " CA ") {
                              int res_no = residue_p->GetSeqNum();
                              clipper::Coord_orth co = coot::co(at);

                              std::map<int, clipper::Coord_orth>::const_iterator it;
                              it = ca_map.find(res_no);
                              if (it != ca_map.end()) {
                                 double dd = (co - it->second).lengthsq();
                                 double d = std::sqrt(dd);
                                 sum_dist += d;
                                 n_found++;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      std::cout << "n_found " << n_found << " sum-dist: " << sum_dist << std::endl;
      return sum_dist;
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_1   = mc.read_pdb(reference_data("7vvl.pdb"));
   int imol_2   = mc.read_pdb(reference_data("6gdg.cif"));
   int imol_map = mc.read_ccp4_map(reference_data("emd_32143.map"), false);

   int imol_3 = mc.copy_fragment_using_cid(imol_2, "//E");
   int n_atoms = mc.get_number_of_atoms(imol_3);
   if (n_atoms > 100) {
      clipper::Coord_orth ref_centre(121, 71, 102);
      std::vector<float> b_variants = {0, 11, 24, 50, 100, 200, 400, 800};
      std::vector<int> n_trials_variants = {10, 24, 50, 100, 200, 400, 800, 2000, 4000};
      std::vector<float> tr_variants = {0, 1, 2, 4, 8, 16};
      int r = 8;
      std::vector<float>::const_iterator it_1;
      std::vector<int>::const_iterator it_2;
      std::vector<float>::const_iterator it_3;

      for (it_1=b_variants.begin(); it_1!=b_variants.end(); ++it_1) {
         for (it_2=n_trials_variants.begin(); it_2!=n_trials_variants.end(); ++it_2) {
            for (it_3=tr_variants.begin(); it_3!=tr_variants.end(); ++it_3) {
               float b = *it_1;
               float tr = *it_3;
               int n_tr = *it_2;

               for (float xo=static_cast<float>(-r); xo<=static_cast<float>(r); xo += 4.0) {
                  for (float yo=static_cast<float>(-r); yo<=static_cast<float>(r); yo += 4.0) {
                     for (float zo=static_cast<float>(-r); zo<=static_cast<float>(r); zo += 4.0) {
                        clipper::Coord_orth offset(xo, yo, zo);
                        clipper::Coord_orth this_centre = ref_centre + offset;
                        int imol_E_copy = mc.copy_fragment_using_cid(imol_3, "/");
                        mc.move_molecule_to_new_centre(imol_3, this_centre.x(), this_centre.y(), this_centre.z());
                        float score = mc.fit_to_map_by_random_jiggle_with_blur_using_cid(imol_E_copy, imol_map, "//", b, n_tr, tr);
                        float score_vs_true = get_score_vs_true_solution(imol_E_copy, imol_1);
                        std::cout << "--- score:: " << xo << " " << yo << " " << zo << " "
                                  << b << " " << n_tr << " " << tr << " score: " << score
                                  << " score_vs_true: " << score_vs_true << std::endl;
                        mc.write_coordinates(imol_3, "fitted.pdb");
                     }
                  }
               }
            }
         }
      }

   }
   return status;
}

int test_peptide_omega(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));

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

int test_omega_5tig_cif(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("5tig.cif"));
   coot::validation_information_t oa = mc.peptide_omega_analysis(imol);
   status = 1; // it didn't crash

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

int test_electro_molecular_representation(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      std::string selection = "//";
      std::string colour = "ByOwnPotential";
      std::string style = "MolecularSurface";

      coot::simple_mesh_t mesh = mc.get_molecular_representation_mesh(imol, selection, colour, style);
      if (mesh.vertices.size() > 10) {
         status = true;
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_replace_fragment(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));

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
                                                                       0, 0, 0, // centre of rotation
                                                                       1, 2, 3); // translation
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

int test_replace_large_fragment(molecules_container_t &mc) {

   // speed test the replacement
   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb("pdb8oie.ent");
   int imol_map = mc.read_ccp4_map("emd_16890.map", false);
   mc.import_cif_dictionary("CLF.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("HCA.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("S5Q.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("ADP.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary( "MG.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("AF3.cif", coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary("SF4.cif", coot::protein_geometry::IMOL_ENC_ANY);

   int imol_new =  mc.copy_fragment_for_refinement_using_cid(imol, "/");
   mc.init_refinement_of_molecule_as_fragment_based_on_reference(imol_new, imol, imol_map);

   std::pair<int, coot::instanced_mesh_t> im_2 = mc.refine(imol_new, 20);
   int refine_status = im_2.first;
   if (refine_status == GSL_CONTINUE) {
      std::pair<int, coot::instanced_mesh_t> im_3 = mc.refine(imol_new, 20);
      refine_status = im_3.first;
   }
   if (refine_status == GSL_CONTINUE) {
      std::pair<int, coot::instanced_mesh_t> im_4 = mc.refine(imol_new, 20);
      refine_status = im_4.first;
   }
   mc.clear_refinement(imol);
   mc.replace_fragment(imol, imol_new, "//");
   status = 1;

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
      float b_factor = 20.0;
      coot::simple_mesh_t mesh = mc.get_gaussian_surface(imol, sigma, contour_level, box_radius, grid_scale, b_factor);
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

int test_instanced_bonds_mesh_v2(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("pdb8ox7.ent"));
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   std::string selection_cid = "//A/1301"; // "//A/1301||//A/456";

   int imol_frag = mc.copy_fragment_using_cid(imol, selection_cid);
   coot::instanced_mesh_t im_frag = mc.get_bonds_mesh_instanced(imol, mode, true, 0.2, 1.0, 2);
   colour_analysis(im_frag);

   coot::instanced_mesh_t im = mc.get_bonds_mesh_for_selection_instanced(imol_frag, selection_cid, mode, true, 0.2, 1.0, 2);
   colour_analysis(im);

   unsigned int n_geoms = im.geom.size();
   for (unsigned int i=0; i<n_geoms; i++) {
      for (unsigned int j=0; j<im.geom[i].instancing_data_A.size(); j++) {
         const auto &item = im.geom[i].instancing_data_A[j];
         const auto &col = item.colour;
         // test that the colour is more (than a bit) more green than it is blue
         // (because grey is the wrong colour, this is a useful test)
         std::cout << "    atom selection (Mg) colour " << glm::to_string(col)  << std::endl;
         if (col[1] > (col[2] + 0.2))
            status = true;
      }
   }
   return status;
}


int test_add_alt_conf(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));

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

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);

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

      mc.refine_residues_using_atom_cid(imol, "//A/262/CA", "QUINTUPLE", 400);
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
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
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
      imol_map = mc.read_ccp4_map(reference_data("emd_11729.map"), 0);
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

int test_map_histogram(molecules_container_t &mc) {

   auto print_hist = [&mc] (int imol_map, float hist_scale_factor) {

      if (mc.is_valid_map_molecule(imol_map)) {
         unsigned int n_bins = 200;
         float zoom_factor = 18.0; // 10 is fine
         coot::molecule_t::histogram_info_t hist = mc.get_map_histogram(imol_map, n_bins, zoom_factor);
         std::cout << "STATS:: mean: " << hist.mean << " sd: " << std::sqrt(hist.variance) << std::endl;
         for (unsigned int i=0; i<hist.counts.size(); i++) {
            float range_start = hist.base + static_cast<float>(i)   * hist.bin_width;
            float range_end   = hist.base + static_cast<float>(i+1) * hist.bin_width;
            std::cout << "    "
                      << std::setw(10) << std::right << range_start << " - "
                      << std::setw(10) << std::right << range_end   << "  "
                      << std::setw(10) << std::right << hist.counts[i] << " ";
            unsigned int n_stars = static_cast<int>(static_cast<float>(hist.counts[i]) * hist_scale_factor);
            for (unsigned int jj=0; jj<n_stars; jj++)
               std::cout << "*";
            std::cout << std::endl;
         }
         return static_cast<int>(hist.counts.size());
      }
      return -1;
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map_1 = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   int imol_map_2 = mc.read_ccp4_map(reference_data("emd_16890.map"), false);
   std::cout << "map_1:" << std::endl;
   int counts_1 = print_hist(imol_map_1, 0.03);
   std::cout << "map_2:" << std::endl;
   int counts_2 = print_hist(imol_map_2, 0.00004);

   if (counts_1 > 10)
      if (counts_2 > 10) status = 1;

   return status;
}

int test_auto_read_mtz(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::vector<molecules_container_t::auto_read_mtz_info_t> imol_maps_5a3h
      = mc.auto_read_mtz(reference_data("5a3h_sigmaa.mtz"));

   for (const auto &item : imol_maps_5a3h) {
      std::cout << "auto-read: map-idx: " << item.idx << " Fobs: " << item.F_obs << " sigFobs: " << item.sigF_obs << " "
                << "Rfree: " << item.Rfree << std::endl;
   }

   std::vector<molecules_container_t::auto_read_mtz_info_t> imol_maps
      = mc.auto_read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"));

   // one of these (the last one) should be observed data without an imol
   if (imol_maps.size() == 3) {
      float rmsd_0 = mc.get_map_rmsd_approx(imol_maps[0].idx);
      float rmsd_1 = mc.get_map_rmsd_approx(imol_maps[1].idx);
      std::cout << "test_auto_read_mtz() rmsds " << rmsd_0 << " " << rmsd_1 << std::endl;
      if (rmsd_0 > 0.3) { // test that the FWT map is the first of the pair
         if (rmsd_1 > 0.1) {
            // what observed data did we find?
            unsigned int n_fobs_found = 0;
            for (unsigned int i=0; i<imol_maps.size(); i++) {
               const auto &mtz_info = imol_maps[i];
               if (! mtz_info.F_obs.empty()) {
                  n_fobs_found++;
                  if (! mtz_info.sigF_obs.empty()) {
                     if (mtz_info.F_obs == "/2vtq/1/FP") {
                        if (mtz_info.sigF_obs == "/2vtq/1/SIGFP") {
                           status = 1;
                        }
                     }
                  }
               }
               if (mtz_info.Rfree == "/HKL_base/HKL_base/FREE") {
                  // we are good.
               } else {
                  status = 0;
               }
            }
            if (n_fobs_found != 1) {
               std::cout << "Too many: " << n_fobs_found << std::endl;
               status = 0;
            }
         }
      }
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
   bool use_rdkit_svg = false;
   std::string s = mc.get_svg_for_residue_type(imol_1, "ATP", use_rdkit_svg, dark_bg);

   if (s.length() > 0) {

      std::ofstream f("ATP.svg");
      f << s;
      f.close();

      {
         mc.import_cif_dictionary("G37.cif", coot::protein_geometry::IMOL_ENC_ANY);
         s = mc.get_svg_for_residue_type(imol_1, "G37", use_rdkit_svg, dark_bg);
         std::ofstream f2("G37.svg");
         f2 << s;
         f2.close();
      }

      {
         mc.import_cif_dictionary("GLC.cif", coot::protein_geometry::IMOL_ENC_ANY);
         s = mc.get_svg_for_residue_type(imol_1, "GLC", use_rdkit_svg, dark_bg);
         std::ofstream f2("GLC.svg");
         f2 << s;
         f2.close();
      }

      status = 1;
   }
   return status;

}

int test_superpose(molecules_container_t &mc) {

   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1-with-gap.pdb"));
   int imol_2 = mc.read_pdb(reference_data("1phk-with-gap.pdb"));

   unsigned int n_pre = mc.get_number_of_molecules();

   coot::atom_spec_t atom_spec_1("A", 227, "", " CA ","");
   coot::atom_spec_t atom_spec_2("A", 256, "", " CA ","");
   mmdb::Atom *at_1 = mc.get_atom(imol_1, atom_spec_1);
   mmdb::Atom *at_2 = mc.get_atom(imol_2, atom_spec_2);

   coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
   coot::Cartesian atom_pos_2 = atom_to_cartesian(at_2);

   double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
   double d1 = std::sqrt(dd);
   std::cout << "test d1 " << d1 << std::endl;

   // std::pair<std::string, std::string> ss_result_pair = mc.SSM_superpose(imol_1, "A", imol_2, "B");
   superpose_results_t ss_results = mc.SSM_superpose(imol_1, "A", imol_2, "A");

   std::cout << "ss_result: info:\n" << ss_results.superpose_info << std::endl;
   std::cout << "ss_result: alnR\n" << ss_results.alignment.first  << std::endl;
   std::cout << "ss_result: alnM\n" << ss_results.alignment.second << std::endl;

   if (false) {
      for (unsigned int i=0; i<ss_results.alignment_info_vec.size(); i++) {
         for (const auto &chain : ss_results.alignment_info_vec[i].cviv) {
            for (const auto &res : chain.rviv) {
               std::cout << res.residue_spec << " " << res.function_value << std::endl;
            }
         }
      }
   }

   if (true) {
      const auto &pairs = ss_results.aligned_pairs;
      for (unsigned int i=0; i<pairs.size(); i++) {
         const auto &r1 = pairs[i].first;
         const auto &r2 = pairs[i].second;
         std::cout << "   " << r1.residue_spec << " " << r2.residue_spec << std::endl;
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
         bool colour_applies_to_non_carbon_atoms_also = true;
         mc.set_user_defined_atom_colour_by_selection(imol, indexed_residues_cids, colour_applies_to_non_carbon_atoms_also);
         coot::instanced_mesh_t im = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, 1);
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
         mc.refine_residues_using_atom_cid(imol, "//A/131", "TRIPLE", 400);
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
         mc.refine_residues_using_atom_cid(imol, "//A/131", "TRIPLE", 400);
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
   // if mcdonald_and_thornton is true, then we need to add hydrogen atoms.
   // if mcdonald_and_thornton is false then we don't need hydrogen atoms - if there are hydrogen atoms
   // then no H-bonds will be made.
   bool mcdonald_and_thornton = true;
   mc.add_hydrogen_atoms(imol); // no hydrogen bonds found without hydrogens in the model
   const std::string &cid_str = "//A/270";
   std::vector<moorhen::h_bond> h_bonds = mc.get_h_bonds(imol, cid_str, mcdonald_and_thornton);

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

int test_bespoke_carbon_colour(molecules_container_t &mc) {

   auto close_float = [] (float a, float b) {
      return fabsf(a - b) < 0.001;
   };

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.get_monomer("LZA");

   if (mc.is_valid_model_molecule(imol)) {
      coot::colour_t col(0.0999, 0.0888, 0.0777);
      mc.set_use_bespoke_carbon_atom_colour(imol, true);
      mc.set_bespoke_carbon_atom_colour(imol, col);
      std::string mode("VDW-BALLS");
      coot::instanced_mesh_t im = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, 1);

      std::cout << "There are " << im.geom.size() << " geoms " << std::endl;
      for (unsigned int ig=0; ig<im.geom.size(); ig++) {
         const auto &g = im.geom[ig];
         std::cout << "geom[" << ig << "] info: "
                   << " n-vertices: "   << g.vertices.size()
                   << " n-triangles: "  << g.triangles.size()
                   << " instancing-A: " << g.instancing_data_A.size()
                   << " instancing-B: " << g.instancing_data_B.size()
                   << std::endl;
      }

      const auto &g = im.geom[0];
      for (unsigned int j=0; j<g.instancing_data_A.size(); j++) {
         const auto &idA = g.instancing_data_A[j];
         std::cout << "object-indx: " << j << " col: " << glm::to_string(idA.colour) << std::endl;
         if (close_float(static_cast<float>(idA.colour.r), 0.0999f))
            if (close_float(static_cast<float>(idA.colour.g), 0.0888f))
               if (close_float(static_cast<float>(idA.colour.b), 0.0777f)) {
                  status = 1;
                  break;
               }
      }
   }
   return status;
}


int test_dark_mode_colours(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.get_monomer("LZA");
   if (mc.is_valid_model_molecule(imol)) {
      std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";
      auto mesh_light = mc.get_bonds_mesh_instanced(imol, mode, false, 0.2, 1.0, 1);
      auto mesh_dark  = mc.get_bonds_mesh_instanced(imol, mode, true,  0.2, 1.0, 1);
      colour_analysis(mesh_light);
      colour_analysis(mesh_dark);
   }

   return status;

}


int test_number_of_hydrogen_atoms(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
     mc.add_hydrogen_atoms(imol);
     int n_hydrogen_atoms = mc.get_number_of_hydrogen_atoms(imol);
     if (n_hydrogen_atoms > 100)
       status = 1;
   }
   return status;
}

int test_cell(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   auto print_cell = [] (const api::cell_t &c) {
      std::cout << c.a << " " << c.b << " " << c.c << " " << c.alpha << " " << c.beta << " " << c.gamma
                << std::endl;
   };

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   api::cell_t c1 = mc.get_cell(imol);
   api::cell_t c2 = mc.get_cell(imol_map);

   print_cell(c1);
   print_cell(c2);

   if (c1.a > 10) status = 1;

   return status;
}

int test_map_centre(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   std::cout << "in test_map_centre() imol_map is " << imol_map << std::endl;
   coot::util::map_molecule_centre_info_t mci = mc.get_map_molecule_centre(imol_map);
   std::cout << "new centre: " << mci.updated_centre.format() << std::endl;

   if (mci.success == true) {
      std::cout << "map centre success " << std::endl;
      if (mci.updated_centre.x() > 10.0)
         status = 1;
   }

   return status;
}

int test_dragged_atom_refinement(molecules_container_t &mc_in) {

   starting_test(__FUNCTION__);
   int status = 0;

   molecules_container_t mc;
   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         coot::Cartesian atom_pos = atom_to_cartesian(at_1);
         coot::Cartesian atom_pos_new = atom_pos + coot::Cartesian(2,2,2); // say
         int imol_new = mc.copy_fragment_for_refinement_using_cid(imol, "//A/265-275"); // make molten molecule
         mc.init_refinement_of_molecule_as_fragment_based_on_reference(imol_new, imol, imol_map);

         mc.write_coordinates(imol_new, "rtest-0.pdb");

         int refine_status = 0;
         std::pair<int, coot::instanced_mesh_t> im_0 = mc.refine(imol_new, 10);
         refine_status = im_0.first;
         mc.write_coordinates(imol_new, "rtest-1.pdb");

         // keep calling this, and displaying the subsequent mesh again as the mouse moves...
         coot::instanced_mesh_t im_1 =
            mc.add_target_position_restraint_and_refine(imol_new, "//A/270/O",
                                                        atom_pos_new.x(), atom_pos_new.y(), atom_pos_new.z(),
                                                        1000);

         mc.write_coordinates(imol_new, "rtest-2.pdb");

         // mouse button release
         mc.clear_target_position_restraints(imol_new);

         std::pair<int, coot::instanced_mesh_t> im_2 = mc.refine(imol_new, 10);
         refine_status = im_2.first;
         if (refine_status == GSL_CONTINUE) {
            std::pair<int, coot::instanced_mesh_t> im_3 = mc.refine(imol_new, 10);
            refine_status = im_3.first;
         }
         if (refine_status == GSL_CONTINUE) {
            std::pair<int, coot::instanced_mesh_t> im_4 = mc.refine(imol_new, 10);
            refine_status = im_4.first;
         }
         if (refine_status == GSL_CONTINUE) {
            std::pair<int, coot::instanced_mesh_t> im_5 = mc.refine(imol_new, 10);
            refine_status = im_5.first;
         }
         if (refine_status == GSL_CONTINUE) {
            std::pair<int, coot::instanced_mesh_t> im_6 = mc.refine(imol_new, 100);
            refine_status = im_6.first;
         }

         // finished mousing:
         mc.clear_refinement(imol);

         // Let's say they liked it:
         mc.replace_fragment(imol, imol_new, "//");

         if (refine_status == GSL_CONTINUE) status = 1; // the atoms still are moving (a bit) (that's what I want
                                                        // for success of this test).

         mc.pop_back(); // remove imol_new
      }
   }
   mc.close_molecule(imol);
   return status;
}

int test_bucca_ml_growing(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   // int imol     = mc.read_pdb(reference_data("1gwd-large-C-terminal-fragment-missing.pdb"));
   std::string fn = "1gwd-large-C-terminal-fragment-missing.pdb";
   fn = "1gwd-118-chop.pdb";
   fn = "1gwd-2-chop.pdb";
   int imol     = mc.read_pdb(reference_data(fn));
   int imol_map = mc.read_mtz(reference_data("1gwd_map.mtz"), "FWT", "PHWT", "FOM", false, false);

   mc.set_refinement_is_verbose();
   mc.geometry_init_standard();
   mc.get_monomer("CL");
   mc.get_monomer("IOD");
   mc.get_monomer("CMO");

   std::string chain_id = "A";
   int res_no = 2; // was 118
   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {
         mc.set_imol_refinement_map(imol_map);
         ::api::cell_t c = mc.get_cell(imol_map);
         std::cout << "debug cell: "
                   << c.a << " " << c.b << " " << c.c << " " << c.alpha << " " << c.beta << " " << c.gamma << std::endl;
         int add_status = 1; // get started.
         while (add_status == 1) {
            std::string atom_cid = "//A/" + std::to_string(res_no);
            std::cout << "-------- building based on " << atom_cid << std::endl;
            add_status = mc.add_terminal_residue_directly_using_bucca_ml_growing_using_cid(imol, atom_cid);

            if (add_status == 1) {

               std::string fnp = "test-pre-ref-" + std::to_string(res_no) + std::string(".pdb");
               mc.write_coordinates(imol, fnp);
               bool do_refine = true;
               if (do_refine) {
                  mc.refine_residue_range(imol, chain_id, res_no, res_no+1, 500);
               }

               std::string fn = "test-" + std::to_string(res_no) + std::string(".pdb");
               mc.write_coordinates(imol, fn);
            }
            res_no += 1; // for next residue
            if (res_no > 129)
               add_status = 0; // force stop.
         }
      } else {
         std::cout << "Not a valid map " << imol_map << std::endl;
      }
   } else {
      std::cout << "Not a valid model " << imol << std::endl;
   }
   return status;
}

int test_user_defined_bond_colours_v2(molecules_container_t &mc) {

   auto close_float = [] (float a, float b) {
      return fabsf(a - b) < 0.001;
   };

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   auto v = mc.get_colour_rules(imol);
   if (false) {
      std::cout << "colour rules: " << std::endl;
      std::cout << "-------------" << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
         std::cout << i << " " << v[i].first << " " << v[i].second << std::endl;
      }
      std::cout << "-------------" << std::endl;
   }

   std::map<unsigned int, std::array<float, 3> > colour_index_map;
   colour_index_map[12] = {1, 0, 1};
   colour_index_map[13] = {0, 1, 1};
   colour_index_map[14] = {0, 0, 1};
   colour_index_map[15] = {0.7, 0.7, 0};

   mc.set_user_defined_bond_colours(imol, colour_index_map);
   std::vector<std::pair<std::string, unsigned int> > indexed_cids;
   indexed_cids.push_back(std::make_pair("//A/1-20",   12));
   indexed_cids.push_back(std::make_pair("//A/22-36",  13));
   indexed_cids.push_back(std::make_pair("//A/46-80",  14));
   indexed_cids.push_back(std::make_pair("//A/90-180", 15));
   bool colour_applies_to_non_carbon_atoms_also = false;
   mc.set_user_defined_atom_colour_by_selection(imol, indexed_cids, colour_applies_to_non_carbon_atoms_also);

   std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";

   {
     auto colour_table = mc.get_colour_table(imol, false);
     for (unsigned int i = 0; i < colour_table.size(); i++) {
       std::cout << "   pre " << i << " " << glm::to_string(colour_table[i])
                 << std::endl;
     }
   }

   auto bonds = mc.get_bonds_mesh_instanced(imol, mode, false, 0.2, 1.0, 1);

   auto colour_table = mc.get_colour_table(imol, false);
   for (unsigned int i=0; i<colour_table.size(); i++) {
      // std::cout << "   post " << i << " " << glm::to_string(colour_table[i]) << std::endl;
   }

   // test a couple of these
   if (close_float(colour_table[12][0], 1.0))
      if (close_float(colour_table[12][1], 0.0))
         if (close_float(colour_table[12][2], 1.0))
            if (close_float(colour_table[13][0], 0.0))
               if (close_float(colour_table[13][1], 1.0))
                  if (close_float(colour_table[13][2], 1.0))
                     status = 1;

   coot::instanced_mesh_t im = mc.get_bonds_mesh_for_selection_instanced(imol, "//A/1-3", "VDW-BALLS", false, 0.1, 1.0, 1);
   if (! im.geom.empty()) {
      const coot::instanced_geometry_t &ig = im.geom[0]; // 0 is spheres
      std::cout << "debug:: in im type A data size: " << ig.instancing_data_A.size() << std::endl;
      if (ig.instancing_data_A.size() == 25) {
         // as it should be!

         // The colour of the [11]th atom ("CA") should be 1,0,1

         for (unsigned int i=0; i<25; i++) {
            const auto &sphere = ig.instancing_data_A[i];
            if (false)
               std::cout << "sphere " << i << " pos " << glm::to_string(sphere.position)
                         << " colour " << glm::to_string(sphere.colour) << std::endl;
            if (i == 11) { // the is "CA" the CA in the first residue (strangely)
               status = 0;
               if (close_float(sphere.colour[0], 1.0))
                  if (close_float(sphere.colour[1], 0.0))
                     if (close_float(sphere.colour[2], 1.0))
                        status = 1;
            }
         }
         
      } else {
         status = 0;
      }
   }
   return status;
}

int test_user_defined_bond_colours_v3(molecules_container_t &mc) {

      // from Filo:

      // const imol = molecules_container.read_pdb('./4ri2.pdb')
      // let colourMap = new cootModule.MapIntFloat3()
      // let indexedResiduesVec = new cootModule.VectorStringUInt_pair()
      // colourMap.set(51, [0.627, 0.529, 0.400])
      // indexedResiduesVec.push_back( { first: '//A', second: 51 })
      // colourMap.set(52, [0.424, 0.627, 0.400])
      // indexedResiduesVec.push_back( { first: '//B', second: 52 })
      // colourMap.set(53, [0.957, 0.263, 0.212])
      // indexedResiduesVec.push_back( { first: '//', second: 53 })
      // molecules_container.set_user_defined_bond_colours(imol, colourMap)
      // molecules_container.set_user_defined_atom_colour_by_selection(imol, indexedResiduesVec, applyColourToNonCarbonAtoms)
      // const bonds = molecules_container.get_bonds_mesh_for_selection_instanced(imol, '//', 'COLOUR-BY-CHAIN-AND-DICTIONARY')

   auto close_float = [] (float a, float b) {
      return fabsf(a - b) < 0.001;
   };

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("pdb4ri2.ent"));

   if (mc.is_valid_model_molecule(imol)) {
      std::map<unsigned int, std::array<float, 3> > colour_map;
      colour_map[51] = {0.627, 0.529, 0.400};
      colour_map[52] = {0.424, 0.627, 0.400};
      colour_map[53] = {0.957, 0.263, 0.212};
      bool C_only = false;
      std::vector<std::pair<std::string, unsigned int> > indexed_residues_cids;
      indexed_residues_cids.push_back(std::make_pair("//A", 51));
      indexed_residues_cids.push_back(std::make_pair("//B", 52));
      indexed_residues_cids.push_back(std::make_pair("//",  53));
      mc.set_user_defined_bond_colours(imol, colour_map);
      mc.set_user_defined_atom_colour_by_selection(imol, indexed_residues_cids, C_only);
      std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";

      // now test the colours:
      auto bonds = mc.get_bonds_mesh_for_selection_instanced(imol, "/", mode, false, 0.2, 1.0, 1);
      auto &geom = bonds.geom;
      auto &vb   = geom[1].instancing_data_B; // bonds
      colour_analysis(bonds);
   }

   return status;
}



int test_is_em_map(molecules_container_t &mc) {

   auto close_float = [] (float a, float b) {
      return fabsf(a - b) < 0.001;
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map = mc.read_ccp4_map(reference_data("emd_25074.map"), 0);
   bool is_EM_map = mc.is_EM_map(imol_map);
   float cl = mc.get_suggested_initial_contour_level(imol_map);
   float rmsd = mc.get_map_rmsd_approx(imol_map);

   std::cout << "test_is_em_map(): EM: " << is_EM_map << " cl " << cl << " with rmsd: " << rmsd
             << " ratio " << cl/rmsd << std::endl;

   if (is_EM_map) {
      if (close_float(cl, 4.0 * rmsd)) {
         coot::util::map_molecule_centre_info_t mci = mc.get_map_molecule_centre(imol_map);
         std::cout << "mci.suggested_radius " << mci.suggested_radius << std::endl;
         if (mci.suggested_radius > 50.0)
            status = 1;
      }
   }

   return status;
}

int test_other_user_define_colours_other(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         std::map<unsigned int, std::array<float, 3> > colour_index_map;
         colour_index_map[21] = {1.11111111, 1.111111, 0};
         std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
         mc.set_user_defined_bond_colours(imol, colour_index_map);
         std::vector<std::pair<std::string, unsigned int> > indexed_cids;
         indexed_cids.push_back(std::make_pair("//A/1-5", 21));
         bool non_carbon_atoms_also_flag = false;
         mc.set_user_defined_atom_colour_by_selection(imol, indexed_cids, non_carbon_atoms_also_flag);
         auto bonds_1 = mc.get_bonds_mesh_for_selection_instanced(imol, "/", mode, false, 0.2, 1.0, 1);
         auto bonds_2 = mc.get_bonds_mesh_instanced(imol, mode, false, 0.2, 1.0, 1);
         auto bonds_3 = mc.get_bonds_mesh_for_selection_instanced(imol, "/", mode, false, 0.2, 1.0, 1);
         auto &geom_1 = bonds_1.geom;
         auto &geom_3 = bonds_3.geom;
 
         auto &vb_1 = geom_1[1].instancing_data_B;
         auto &vb_3 = geom_3[1].instancing_data_B;

         for (unsigned int i=0; i<vb_1.size(); i++) {
             if (i > 10) continue;
             auto col = vb_1[i].colour;
             std::cout << "instancing colour_1: " << i << " " << glm::to_string(col) << "\n";
         }

         for (unsigned int i=0; i<vb_3.size(); i++) {
            // if (i > 10) continue;
            auto col = vb_3[i].colour;
            std::cout << "instancing colour_3: " << i << " " << glm::to_string(col) << "\n";
         }
      }
   }
   return status;
}


int test_self_restraints(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   mc.generate_self_restraints(imol, 5.0);
   coot::instanced_mesh_t im = mc.get_extra_restraints_mesh(imol, 0);
   unsigned int size = im.geom[0].instancing_data_B.size();
   if (size > 10) {
      int imol_2 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
      mc.generate_local_self_restraints(imol_2, 5, "//A/10");
      im = mc.get_extra_restraints_mesh(imol_2, 0);
      unsigned int size_1 = im.geom[0].instancing_data_B.size();
      mc.generate_local_self_restraints(imol_2, 5, "//A/11-22");
      im = mc.get_extra_restraints_mesh(imol_2, 0);
      unsigned int size_2 = im.geom[0].instancing_data_B.size();
      std::cout << "size-1 " << size_1 << " size-2: " << size_2 << std::endl;
      if (size_1 > 10) {
         if (size_2 > 2 * size_1) {
            status = 1;
         }
      }
   }
   return status;
}

int test_read_extra_restraints(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_2 = mc.read_pdb(reference_data("3pzt.pdb"));
   mc.read_extra_restraints(imol_1, reference_data("moorhen-tutorial-structure-number-1-prosmart.txt"));
   coot::instanced_mesh_t im = mc.get_extra_restraints_mesh(imol_1, 0);
   if (! im.geom.empty()) {
      std::cout << "instancing_data_B size " << im.geom[0].instancing_data_B.size() << std::endl;
      if (im.geom[0].instancing_data_B.size() > 10) status = 1;
   } else {
      std::cout << "ERROR:: im geom is empty" << std::endl;
   }
   return status;
}




int test_colour_map_by_other_map(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_map_1 = mc.read_ccp4_map(reference_data("emd_16890.map"), false);
   int imol_map_2 = mc.read_ccp4_map("scale_res_emd_16890.mrc", false);
   if (mc.is_valid_map_molecule(imol_map_1)) {
      if (mc.is_valid_map_molecule(imol_map_2)) {
         coot::simple_mesh_t mesh = mc.get_map_contours_mesh_using_other_map_for_colours(imol_map_1, imol_map_2,
                                                                                         160, 160, 160,
                                                                                         100, 0.16,
                                                                                         0.3, 0.9, false);
         std::cout << "test: mesh v and t: " << mesh.vandt() << std::endl;
         colour_analysis(mesh);
         if (mesh.vertices.size() > 1000) status = true;
      }
   }
   return status;

}

int test_residues_near_residues(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      std::string residue_cid = "//A/60";
      std::vector<coot::residue_spec_t> neighbour_specs = mc.get_residues_near_residue(imol, residue_cid, 5);
      std::cout << "debug:: found " << neighbour_specs.size() << " neighbouring residues" << std::endl;
      if (neighbour_specs.size() > 3) status = 1;
   }
   return status;
}

int test_ncs_chains(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("7bqx-assembly1.cif"));
   if (mc.is_valid_model_molecule(imol)) {
      int n_chains = 0;
      auto vvc = mc.get_ncs_related_chains(imol);
      std::cout << "found " << vvc.size() << " NCS-chain groups" << std::endl;
      for (const auto &vc : vvc) {
         for (const auto &c : vc) {
            // std::cout << " " << c;
            n_chains++;
         }
         // std::cout << std::endl;
      }
      std::cout << "Found " << n_chains << " chains in total" << std::endl;
      if (n_chains == 95) status = 1;
   }
   return status;
}

int test_pdbe_dictionary_depiction(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   // this test doesn't have a good/correct success criterion.
   // Just that the file is written. It is up to us to look at the image.

   mc.import_cif_dictionary(reference_data("MOI.restraints.cif"), coot::protein_geometry::IMOL_ENC_ANY); // from Oliver Smart

   mc.write_png("MOI", coot::protein_geometry::IMOL_ENC_ANY, "MOI-depiction.png");

   // if (coot::file_exists("MOI-depiction.png")) status = 1; // not a good test.

   bool use_rdkit_rendering = true;
   bool dark_background = false;
   std::string svg = mc.get_svg_for_residue_type(coot::protein_geometry::IMOL_ENC_ANY, "MOI", use_rdkit_rendering, dark_background);
   std::ofstream f("MOI.svg");
   f << svg;
   f.close();
   return status;
}


int test_cif_writer(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   mc.import_cif_dictionary(reference_data("HEM.restraints.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   std::string s1 = mc.get_cif_restraints_as_string("xHEMx", coot::protein_geometry::IMOL_ENC_ANY);
   std::string s2 = mc.get_cif_restraints_as_string("HEM",   coot::protein_geometry::IMOL_ENC_ANY);
   if (s1.length() == 0)
      if (s2.length() > 10)
         status = 1;
   if (false) {
      std::cout << "debug s2 length " << s2.length() << std::endl;
      std::ofstream f("s2.out");
      f << s2;
      f.close();
   }
   return status;
}

int test_cif_gphl_chem_comp_info(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   mc.import_cif_dictionary(reference_data("HEM.restraints.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   const auto &info = mc.get_gphl_chem_comp_info("HEM", coot::protein_geometry::IMOL_ENC_ANY);
   if (info.size() > 6) status = 1;
   for (unsigned int i=0; i<info.size(); i++) {
      std::cout << "   " << i << " " << info[i].first << " " << info[i].second << std::endl;
   }
   return status;
}

int test_pdb_as_string(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   std::string s1 = mc.molecule_to_PDB_string(imol);
   mc.auto_fit_rotamer(imol, "A", 61, "", "", imol_map);
   std::string s2 = mc.molecule_to_PDB_string(imol);

   std::cout << "test_pdb_as_string(): lengths " << s1.length() << " " << s2.length() << std::endl;

   if (s1.length() == s2.length()) status = 1;

   if (false) {
      std::ofstream f("test_pdb_as_string.pdb");
      f << s2;
      f.close();
   }

   return status;

}

int test_mmcif_as_string(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("2vtq.cif"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   std::string s1 = mc.molecule_to_mmCIF_string(imol);
   mc.auto_fit_rotamer(imol, "A", 61, "", "", imol_map);
   std::string s2 = mc.molecule_to_mmCIF_string(imol);

   std::cout << "test_mmcif_as_string(): lengths " << s1.length() << " " << s2.length() << std::endl;

   if (s1.length() == s2.length()) status = 1;

   if (true) {
      std::ofstream f1("test_mmcif_as_string_1.mmcif");
      f1 << s1;
      f1.close();
      std::ofstream f2("test_mmcif_as_string_2.mmcif");
      f2 << s2;
      f2.close();
   }

   return status;

}

int test_mmcif_atom_selection(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string fn = "1ej6-assembly1.cif";
   std::cout << "reading " << fn << std::endl;
   int imol = mc.read_pdb(reference_data(fn));
   mmdb::Manager *mol = mc.get_mol(imol);
   if (mol) {
      int n_selected_atoms_1 = 0;
      int n_selected_atoms_2 = 0;
      int n_selected_atoms_3 = 0;
      mmdb::Atom **selected_atoms_1 = 0;
      mmdb::Atom **selected_atoms_2 = 0;
      mmdb::Atom **selected_atoms_3 = 0;
      int selHnd_1 = mol->NewSelection();
      int selHnd_2 = mol->NewSelection();
      int selHnd_3 = mol->NewSelection();
      std::cout << "selecting //A" << std::endl;
      mol->Select(selHnd_1, mmdb::STYPE_ATOM, "//A",     mmdb::SKEY_NEW);
      std::cout << "selecting //A-1" << std::endl;
      mol->Select(selHnd_2, mmdb::STYPE_ATOM, "//A-1",   mmdb::SKEY_NEW);
      std::cout << "selecting //A-1,A" << std::endl;
      mol->Select(selHnd_3, mmdb::STYPE_ATOM, "//A-1,A", mmdb::SKEY_NEW);
      mol->GetSelIndex(selHnd_1, selected_atoms_1, n_selected_atoms_1);
      mol->GetSelIndex(selHnd_2, selected_atoms_2, n_selected_atoms_2);
      mol->GetSelIndex(selHnd_3, selected_atoms_3, n_selected_atoms_3);
      std::cout << "n-selected " << n_selected_atoms_1 << " " << n_selected_atoms_2 << " " << n_selected_atoms_3
                << std::endl;
      // there should be nothing in A-1 selection that is in //A
      unsigned int n_matcher = 0;
      for (int i=0; i<n_selected_atoms_1; i++) {
         if (i >= 100) break;
         mmdb:: Atom *at_1 = selected_atoms_1[i];
         for (int j=0; i<n_selected_atoms_2; j++) {
            mmdb:: Atom *at_2 = selected_atoms_2[j];
            if (at_1 == at_2) {
               n_matcher++;
               break;
            }
         }
      }
      std::cout << "Looked for 100 atoms and found " << n_matcher << " matchers" << std::endl;
      if (n_matcher == 0) status = 1;
   }
   return status;
}

int test_contouring_timing(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      float contour_level = 0.12;

      clipper::Coord_orth p(55, 10, 10);
      for (unsigned int i=0; i<80; i++) {
         float radius = i;
         coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
         double t = mc.get_contouring_time();
         std::cout << "contouring time: " << i << " " << t << std::endl;
         if (t > 10) status = true;
      }
   }

   return status;
}

int test_test_the_threading(molecules_container_t &mc) {

   int status = 1; // no faiiure

   for (unsigned int i=0; i<50; i++) {
      double r = mc.test_the_threading(i);
      std::cout << " test_threading: " << i << " " << r << std::endl;
   }

   return status;
}

int test_thread_launching(molecules_container_t &mc) {

   int status = 1; // no failure
   for (unsigned int i=1; i<50; i++) {
      for (unsigned int j=1; j<50; j++) {
         double t = mc.test_launching_threads(j, i);
         std::cout << " launching " << j << " " << i << " : " << t << " micro-seconds " << std::endl;
      }
   }
   return status;
}

int test_thread_pool(molecules_container_t &mc) {

   int status = 1; // no failure
   for (unsigned int j=1; j<50; j++) {
      double t = mc.test_thread_pool_threads(j);
      std::cout << " launching " << j << " " << t << " micro-seconds " << std::endl;
   }
   return status;
}

int test_long_name_ligand_cif_merge(molecules_container_t &mc) {

   int status = 0;
   int imol = mc.read_pdb(reference_data("8a2q.cif"));
   mc.import_cif_dictionary(reference_data("7ZTVU.cif"), coot::protein_geometry::IMOL_ENC_ANY);

   return status;
}

#ifdef USE_GEMMI
#include "gemmi/mmread.hpp"
#include "gemmi/mmdb.hpp"

int test_disappearing_ligand(molecules_container_t &mc) {

   int status = 0;
   // int imol = mc.read_pdb(reference_data("6ttq.cif")); // needs gemmi
   int imol = mc.read_pdb(reference_data("8a2q.cif"));
   mmdb::Manager *mol = mc.get_mol(imol);
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            std::cout << "    " << chain_p->GetChainID() << " " << n_res << " residues " << std::endl;
         }
      }
   }
   mc.import_cif_dictionary(reference_data("MOI.restraints.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   int imol_lig = mc.get_monomer("MOI");
   std::string sl = std::to_string(imol_lig);
   std::pair<int, std::vector<merge_molecule_results_info_t> > ss = mc.merge_molecules(imol, sl);
   mc.write_coordinates(imol, "merged.cif");
   gemmi::Structure st = gemmi::read_structure_file("merged.cif");

   return status;
}
#endif

// 20240205-PE this is not a good name for this test. The failure is not in merging, the failure
// is in writing out the cif file.
//
int test_ligand_merge(molecules_container_t &mc) {

   auto test_mmdb = [] () {
      // int imol_2 = mc.read_pdb(reference_data("2vtq.cif"));
      // mc.write_coordinates(imol_2, "2vtq-just-input-output.cif");
      mmdb::Manager *mol = new mmdb::Manager;
      mol->ReadCoorFile("2vtq.cif");
      mol->WriteCIFASCII("2vtq-input-output-pure-mmdb.cif");
      delete mol;
   };

   int status = 0;

   test_mmdb();
   int imol = mc.read_pdb(reference_data("2vtq-sans-ligand.cif"));
   mc.write_coordinates(imol, "2vtq-sans-ligand-just-input-output.cif");
   mmdb::Manager *mol = mc.get_mol(imol);
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            std::cout << "    " << chain_p->GetChainID() << " " << n_res << " residues " << std::endl;
         }
      }
   }
   mc.import_cif_dictionary(reference_data("MOI.restraints.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   int imol_lig = mc.get_monomer("MOI");
   std::string sl = std::to_string(imol_lig);
   std::pair<int, std::vector<merge_molecule_results_info_t> > ss = mc.merge_molecules(imol, sl);
   mc.write_coordinates(imol, "2vtq-sans-ligand-with-merged-MOI.cif");

   // we test that the output is sane by looking for an atom that has unset/default/"." for atom and element
   // columns 3 and 4 (starting from 1).
   std::string cif_for_testing = "2vtq-sans-ligand-just-input-output.cif";
   cif_for_testing = "2vtq-sans-ligand-just-input-output.cif"; // it's not the merge, it's the writing.
   if (coot::file_exists(cif_for_testing)) {
      bool found_bogus_atom = false;
      std::ifstream f(cif_for_testing.c_str());
      if (f) {
         std::string line;
         while (std::getline(f, line)) {
            std::vector<std::string> parts = coot::util::split_string_no_blanks(line);
            if (parts.size() > 10) {
               if (parts[0] == "ATOM") {
                  if (parts[2] == ".") {
                     if (parts[3] == ".") {
                        std::cout << "found bogus null atom in cif output " << cif_for_testing << std::endl;
                        found_bogus_atom = true;
                     }
                  }
               }
            }
         }
      }
      if (! found_bogus_atom) status = 1; // OK then I suppose
   }
   return status;
}


int test_gltf_export(molecules_container_t &mc) {

   auto make_multi_cid = [] (const std::vector<coot::residue_spec_t> &neighbs) {

      std::string multi_cid;
      if (neighbs.size() == 1) {
         multi_cid = "//" + neighbs[0].chain_id + "/" + std::to_string(neighbs[0].res_no);
      }
      if (neighbs.size() > 1) {
         unsigned int m = neighbs.size() - 1;
         for (unsigned int i=0; i<m; i++) {
            const auto &n = neighbs[i];
            std::string rs = "//" + n.chain_id + "/" + std::to_string(n.res_no);
            multi_cid += rs;
            multi_cid += "||";
         }
         multi_cid += "//" + neighbs.back().chain_id + "/" + std::to_string(neighbs.back().res_no);
      }
      return multi_cid;
   };

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("2vtq.cif"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   clipper::Coord_orth p(25, 4, 62);
   float radius = 10;
   float contour_level = 0.4;
   std::cout << "-------------------------------------------------- map mesh " << std::endl;
   coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
   map_mesh.export_to_gltf("map-around-ligand.glb", true);

   std::cout << "-------------------------------------------------- ligand mesh " << std::endl;

   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   int imol_lig = mc.get_monomer("LZA");
   int imol_frag = mc.copy_fragment_using_cid(imol, "//A/1299");
   std::cout << "test_gltf_export() imol_frag " << imol_frag << std::endl;
   coot::instanced_mesh_t im    = mc.get_bonds_mesh_instanced(imol_frag, mode, true, 0.1, 1.0, 1);
   coot::simple_mesh_t sm_lig = coot::instanced_mesh_to_simple_mesh(im);
   sm_lig.export_to_gltf("lig.glb", true);

   std::cout << "-------------------------------------------------- neighbour mesh " << std::endl;
   std::vector<coot::residue_spec_t> neighbs = mc.get_residues_near_residue(imol, "//A/1299", 4.2);
   std::string multi_cid = make_multi_cid(neighbs);
   mc.set_draw_missing_residue_loops(false);
   coot::instanced_mesh_t im_neighbs = mc.get_bonds_mesh_for_selection_instanced(imol, multi_cid, mode, true, 0.15, 1.0, 1);
   coot::simple_mesh_t sm_neighbs = coot::instanced_mesh_to_simple_mesh(im_neighbs);
   sm_neighbs.export_to_gltf("neighbs.glb", true);

   return status;
}


int test_gltf_export_via_api(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("2vtq.cif"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   clipper::Coord_orth p(25, 4, 62);
   float radius = 10;
   float contour_level = 0.4;
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   int imol_lig = mc.get_monomer("LZA");
   mc.export_map_molecule_as_gltf(imol_map, p.x(), p.y(), p.z(), radius, contour_level, "map-around-ligand.glb");
   mc.export_model_molecule_as_gltf(imol, "//A/1299", mode, true, 0.2, 1.4, 2, true, true, "fat-ligand.glb");

   struct stat buf_1;
   int istat_1 = stat("map-around-ligand.glb", &buf_1);
   if (istat_1 == 0) {
      if (buf_1.st_size > 1000000) {

         struct stat buf_2;
         int istat_2 = stat("fat-ligand.glb", &buf_2);
         if (istat_2 == 0) {
            if (buf_2.st_size > 100000) {
               status = 1;
            }
         }
      }
   }
   return status;
}


int test_5char_ligand_merge(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   mc.import_cif_dictionary(reference_data("acedrg-7z-new.cif"), imol_enc);
   int imol_lig = mc.get_monomer("7ZTVU");
   if (mc.is_valid_model_molecule(imol)) {
      mc.merge_molecules(imol, std::to_string(imol_lig));
      mc.write_coordinates(imol, "5-char-ligand-merged.cif");
      status = 1;
   }
   return status;
}

int test_mask_atom_selection(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      int imol_masked = mc.mask_map_by_atom_selection(imol, imol_map, "//A/1-20||//A/50-70||//A/100-120", 4.0, true);
      mc.write_map(imol_masked, "multi-cid-masked.map");
      // if the masking worked there will be zero density at the CA of residue 30.
      mmdb::Manager *mol = mc.get_mol(imol);
      int sel_hnd = mol->NewSelection();
      mol->SelectAtoms(sel_hnd, 1, "A", 30, "", 30, "", "VAL", " CA ", "*", "");
      mmdb::Atom **atom_selection = 0; // member data - cleared on destruction
      int n_selected_atoms = 0;
      mol->GetSelIndex(sel_hnd, atom_selection, n_selected_atoms);
      std::cout << "------------- got n_selected_atoms " << n_selected_atoms << std::endl;
      if (n_selected_atoms > 0) {
         for (int i=0; i<n_selected_atoms; i++) {
            mmdb::Atom *at = atom_selection[i];
            clipper::Coord_orth pos(at->x, at->y, at->z);
            std::cout << "in test_mask_atom_selection() found atom "
                      << at->GetResName() << " " << at->GetSeqNum() << " "
                      << ":" << at->GetAtomName() << ": " << pos.format() << std::endl;
            float f = mc.get_density_at_position(imol_masked, at->x, at->y, at-> z);
            if (f < 0.00001) {
               status = 1;
            }
         }
      }
      mol->DeleteSelection(sel_hnd);
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
      status += run_test(test_rota_dodecs_mesh,      "rotamer dodecahedra mesh", mc);
      status += run_test(test_rsr_using_residue_range, "rsr using residue range", mc);
      status += run_test(test_copy_fragment_using_cid, "copy-fragment using cid", mc);
      status += run_test(test_no_dictionary_residues,  "no-dictionary residues", mc);
      status += run_test(test_cis_trans,             "cis_trans conversion",     mc);
      status += run_test(test_rsr_using_atom_cid,    "rsr using atom cid",       mc);
      status += run_test(test_auto_fit_rotamer_1,    "auto-fit rotamer",         mc);
      status += run_test(test_auto_fit_rotamer_2,    "auto-fit rotamer t2",      mc);
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
      status += run_test(test_bonds_mesh,            "bonds mesh",               mc);
      status += run_test(test_eigen_flip,            "Eigen Flip",               mc);
      status += run_test(test_read_a_map,            "read a map",               mc);
      status += run_test(test_add_compound,          "add compound",             mc);
      status += run_test(test_weird_delete,          "delete II",                mc);
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
      status += run_test(test_rotamer_validation,    "rotamer validation",       mc);
      status += run_test(test_ligand_fitting_here,   "Ligand fitting here",      mc);
      status += run_test(test_ligand_contact_dots,   "ligand contact dots",      mc);
      status += run_test(test_difference_map_peaks,  "Difference Map Peaks",     mc);
      status += run_test(test_rama_validation,       "rama validation 2",        mc); // for the plot, not the graph
      status += run_test(test_ramachandran_analysis, "ramachandran analysis",    mc); // for the graph, not the plot
      status += run_test(test_non_standard_residues, "non-standard residues",    mc);
      status += run_test(test_import_cif_dictionary, "import cif dictionary",    mc);
      status += run_test(test_add_terminal_residue,  "add terminal residue",     mc);
      status += run_test(test_sequence_generator,    "Make a sequence string",   mc);
      status += run_test(test_instanced_rota_markup, "instanced rotamer mesh",   mc);
      status += run_test(test_new_position_for_atoms,"new positions for atoms",  mc);
      status += run_test(test_molecular_representation, "molecular representation mesh", mc);
      status += run_test(test_fill_partial,          "fill partially-filled residues", mc);
   }

   // status += run_test(test_multiligands_lig_bonding, "some multiligands bonding", mc);

   // status += run_test(test_gltf_export_via_api, "gltf via api", mc);

   // status += run_test(test_multi_ligand_ligands, "multi-ligand ligands", mc);

   // status += run_test(test_updating_maps, "updating maps", mc);

   // status += run_test(test_disappearing_ligand, "disappearning ligand", mc);

   // status += run_test(test_long_name_ligand_cif_merge, "Long-name ligand cif merge", mc);

   // status += run_test(test_pdbe_dictionary_depiction, "PDBe dictionary depiction", mc);

   // status += run_test(test_user_defined_bond_colours_v3, "user-defined colours v3", mc);

   // status += run_test(test_gltf_export, "glTF export", mc);

   // status += run_test(test_5char_ligand_merge, "5-char ligand merge", mc);

   // status += run_test(test_density_mesh,          "density mesh",             mc);

   // status += run_test(test_thread_pool, "thread pool",    mc);

   // status += run_test(test_thread_launching, "thread launching",    mc);

   // status += run_test(test_cif_gphl_chem_comp_info, "extracting gphl info",    mc);

   // status += run_test(test_test_the_threading, "threading speed test",    mc);

   // status += run_test(test_contouring_timing, "contouring timing",    mc);

   // status += run_test(test_mmcif_atom_selection, "mmCIF atom selection",    mc);

   // status += run_test(test_mmcif_as_string, "mmCIF as string",    mc);

   // status += run_test(test_pdb_as_string, "PDB as string",    mc);

   // status += run_test(test_cif_writer, "mmCIF dictionary writer",    mc);

   // status += run_test(test_pdbe_dictionary_depiction, "PDBe dictionary depiction",    mc);

   // status += run_test(test_rsr_using_multi_atom_cid, "multi-atom-cid RSR",    mc);

   //status += run_test(test_rsr_using_atom_cid, "atom-cid RSR",    mc);

   // status += run_test(test_residues_near_residues, "residues near residues",    mc);

   // status += run_test(test_import_cif_dictionary, "import cif dictionary",    mc);

   // status += run_test(test_electro_molecular_representation, "electro molecular representation mesh", mc);

   // status += run_test(test_replace_fragment,      "replace fragment",         mc);

   // status += run_test(test_ncs_chains,      "NCS chains",         mc);

   // status += run_test(test_omega_5tig_cif,      "Omega for 5tig cif",         mc);

   // status += run_test(test_jiggle_fit_params, "actually testing for goodness pr params", mc);

   // status += run_test(test_dark_mode_colours, "light vs dark mode colours", mc);

   // status += run_test(test_read_extra_restraints, "read extra restraints", mc);

   // status += run_test(test_electro_molecular_representation, "electro molecular representation mesh", mc);

   // status += run_test(test_map_histogram, "map histogram", mc);

   // status += run_test(test_auto_read_mtz, "auto-read-mtz", mc);

   // status += run_test(test_read_a_missing_map, "read a missing map file ", mc);

   // status += run_test(test_colour_map_by_other_map, "colour-map-by-other-map", mc);

   // status += run_test(test_jiggle_fit_with_blur, "Jiggle-fit-with-blur", mc);

   // status += run_test(test_something_filo, "Self something filo", mc);

   // status += run_test(test_self_restraints, "Self restraints mesh", mc);

   // status += run_test(test_other_user_define_colours_other, "New colour test", mc);

   // status += run_test(test_is_em_map, "test if EM map flag is correctly set", mc);

   // status += run_test(test_user_defined_bond_colours_v2, "user-defined bond colours v2", mc);

   // status += run_test(test_ramachandran_analysis, "--- current_test ---", mc);

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

   // status += run_test(test_auto_read_mtz, "auto-read MTZ", mc);

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

   // status = run_test(test_moorhen_h_bonds, "moorhen H-bonds ", mc);

   // status = run_test(test_number_of_hydrogen_atoms, "number of hydrogen atoms ", mc);

   // status += run_test(test_molecular_representation, "molecular representation mesh", mc);

   // status += run_test(test_cell, "cell", mc);

   // status += run_test(test_map_centre, "map centre", mc);

   // status += run_test(test_dragged_atom_refinement, "dragged atom refinement", mc);

   // status = run_test(test_bespoke_carbon_colour, "bespoke carbon colours ", mc);

   // Note to self:
   //
   // change the autofit_rotamer test so that it tests the change of positions of the atoms of the neighboring residues.

   // status = run_test(test_replace_model_from_file, "replace model from file", mc);

   // status = run_test(test_user_defined_bond_colours, "user-defined bond colours", mc);

   // status = run_test(test_replace_map, "replace map from mtz", mc);

   // status = run_test(test_residue_name_group, "residue name group", mc);

   // status = run_test(test_superpose, "SSM superpose ", mc);

   // status = run_test(test_rsr_using_residue_range, "test_rsr using residue range", mc);

   // status = run_test(test_bucca_ml_growing, "Bucca ML growing", mc);

   // status += run_test(test_mask_atom_selection, "mask atom selection", mc);

   // status += run_test(test_instanced_bonds_mesh_v2, "test instanced bond selection v2", mc);

   // status += run_test(test_ligand_merge, "test ligand merge", mc);

   // status += run_test(test_add_terminal_residue_v2, "test add terminal residue v2", mc);

   // status += run_test(test_auto_read_mtz, "test auto_read_mtz", mc);

   // status += run_test(test_replace_fragment, "replace fragment",         mc);

   // status += run_test(test_dragged_atom_refinement, "dragged atom refinement", mc);

   // status += run_test(test_rsr_using_atom_cid,    "rsr using atom cid",       mc);

   // status += run_test(test_replace_large_fragment,      "refine and replace large fragment",         mc);

   // status += run_test(test_auto_read_mtz, "test ------ ", mc);

   // status += run_test(test_ligand_fitting_in_map, "ligand fitting in map",    mc);

   // status += run_test(test_write_map_is_sane, "write map is sane",    mc);

   status += run_test(test_ligand_fitting_in_map, "ligand fitting in map",    mc);

   int all_tests_status = 1; // fail!
   if (status == n_tests) all_tests_status = 0;

   return all_tests_status;

}
