
#include <iostream>
#include <iomanip>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#ifdef USE_GEMMI
#include <gemmi/mmread.hpp>
#include <gemmi/mmdb.hpp>
#include <gemmi/enumstr.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/polyheur.hpp>
#endif
#include "MoleculesToTriangles/CXXClasses/MyMolecule.h"
#include "molecules-container.hh"
#include "coot-utils/acedrg-types-for-residue.hh"
#include "filo-tests.hh"
#include "lucrezia-tests.hh"

void starting_test(const char *func) {
   std::cout << "\nStarting " << func << "()" << std::endl;
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

std::vector<std::pair<glm::vec4, unsigned int> >
colour_analysis(const coot::simple_mesh_t &mesh) {

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

   return colour_count;
}

std::vector<std::pair<glm::vec4, unsigned int> >
colour_analysis(const coot::instanced_mesh_t &mesh) {

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

   return colour_count;
}

class colour_analysis_row {
public:
   colour_analysis_row(const glm::vec4 &v4, unsigned int cc) : col(v4), count(cc) {};
   glm::vec4 col;
   unsigned int count;
};

std::vector<colour_analysis_row> get_colour_analysis(const coot::instanced_mesh_t &mesh) {

   auto is_near_colour = [] (const glm::vec4 &col_1, const glm::vec4 &col_2) {
      float cf = 0.04;
      if (std::fabs(col_2.r - col_1.r) < cf)
         if (std::fabs(col_2.g - col_1.g) < cf)
            if (std::fabs(col_2.b - col_1.b) < cf)
               if (std::fabs(col_2.a - col_1.a) < cf)
                  return true;
      return false;
   };

   auto sorter = [] (const colour_analysis_row &p1,
                     const colour_analysis_row &p2) {
      if (p1.col[0] == p2.col[0]) {
         return (p1.col[1] > p2.col[1]);
      } else {
         return (p1.col[0] > p2.col[0]);
      }
   };

   std::vector<colour_analysis_row> colour_count;

   for (unsigned int i=0; i<mesh.geom.size(); i++) {
      const coot::instanced_geometry_t &ig = mesh.geom[i];
      for (unsigned int jj=0; jj<ig.instancing_data_A.size(); jj++) {
         const auto &col =  ig.instancing_data_A[jj].colour;
         bool found_col = false;
         for (unsigned int j=0; j<colour_count.size(); j++) {
            if (is_near_colour(col, colour_count[j].col)) {
               colour_count[j].count ++;
               found_col = true;
               break;
            }
         }
         if (! found_col) {
            colour_count.push_back(colour_analysis_row(col, 1));
         }
      }

      for (unsigned int jj=0; jj<ig.instancing_data_B.size(); jj++) {
         const auto &col =  ig.instancing_data_B[jj].colour;
         bool found_col = false;
         for (unsigned int j=0; j<colour_count.size(); j++) {
            if (is_near_colour(col, colour_count[j].col)) {
               colour_count[j].count++;
               found_col = true;
               break;
            }
         }
         if (! found_col) {
            colour_count.push_back(colour_analysis_row(col, 1));
         }
      }

   }


   for (unsigned int i=0; i<mesh.markup.vertices.size(); i++) {
      const auto &vertex = mesh.markup.vertices[i];
      const glm::vec4 &col = vertex.color;
      bool found_col = false;
      for (unsigned int j=0; j<colour_count.size(); j++) {
         if (is_near_colour(col, colour_count[j].col)) {
            colour_count[j].count ++;
            found_col = true;
            break;
         }
      }
      if (! found_col) {
         colour_count.push_back(colour_analysis_row(col, 1));
      }
   }

   std::sort(colour_count.begin(), colour_count.end(), sorter);

   std::cout << "INFO:: get_colour_analysis(): " << colour_count.size() << " colours" << std::endl;
   for (unsigned int i=0; i<colour_count.size(); i++)
      std::cout << "    " << glm::to_string(colour_count[i].col) << " "
                << std::setw(7) << std::right << colour_count[i].count << std::endl;

   return colour_count;
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
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);
   return status;
}

int test_auto_fit_rotamer_2(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0; // initially fail status

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
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);
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

   mc.close_molecule(imol_diff_map);
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);
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
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);
   return status;
}

// Test for set_residue_to_rotamer_number()
//
int test_set_residue_to_rotamer_number(molecules_container_t &mc) {

    starting_test(__FUNCTION__);
    int status = 0;

    // Load test structure model
    int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
    if (!mc.is_valid_model_molecule(imol)) {
        std::cout << "Failed to load model molecule" << std::endl;
        return status;
    }

    // Pick a known residue (e.g. chain A, residue 270)
    int rotamer_number = 2; // Try setting to rotamer #2

    // Get a key atom to track the change (e.g. CG for ASP)
    std::string chain_id = "A";
    int res_no = 269;
    std::string residue_cid = "//A/269";
    std::string alt_conf;
    coot::atom_spec_t atom_spec(chain_id, res_no, "", " CG ", alt_conf);
    mmdb::Atom* at_start = mc.get_atom(imol, atom_spec);
    if (!at_start) {
        std::cout << "Failed to find atom CG in starting residue" << std::endl;
        mc.close_molecule(imol);
        return status;
    }
    coot::Cartesian pos_start = coot::Cartesian(at_start->x, at_start->y, at_start->z);

    int result = mc.set_residue_to_rotamer_number(imol, residue_cid, alt_conf, rotamer_number);

    // Find the atom after the function call
    mmdb::Atom* at_end = mc.get_atom(imol, atom_spec);
    if (!at_end) {
        std::cout << "Failed to find atom CG after rotamer set" << std::endl;
        mc.close_molecule(imol);
        return status;
    }
    coot::Cartesian pos_end = coot::Cartesian(at_end->x, at_end->y, at_end->z);

    double dist_moved = std::sqrt(coot::Cartesian::lengthsq(pos_start, pos_end));
    std::cout << "CG atom moved " << dist_moved << " Å by set_residue_to_rotamer_number()" << std::endl;

    // Success: function returned 1 and atom moved
    if (result == 1 && dist_moved > 0.3) {
        status = 1;
    }

    mc.close_molecule(imol);
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
   mc.close_molecule(imol_map);

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
   mc.close_molecule(imol_map);

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

         std::string cid = "//A/14/CA";
         // std::string cid = "//A/187/CA";

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

         mc.write_coordinates(imol, "pre-refine.pdb");
         mc.add_to_non_drawn_bonds(imol, cid);
         int refine_status = mc.refine_residues_using_atom_cid(imol, cid, mode, n_cycles);
         std::cout << "debug:: refine_status " << refine_status << std::endl;
         coot::Cartesian pt_n_1_post = atom_to_cartesian(at_n_1);
         coot::Cartesian pt_n_2_post = atom_to_cartesian(at_n_2);
         mc.write_coordinates(imol, "post-refine.pdb");

         double dd_n_1 = coot::Cartesian::lengthsq(pt_n_1_pre, pt_n_1_post);
         double dd_n_2 = coot::Cartesian::lengthsq(pt_n_2_pre, pt_n_2_post);
         double d_1 = std::sqrt(dd_n_1);
         double d_2 = std::sqrt(dd_n_2);

         std::cout << "debug:: rsr distances " << d_1 << " " << d_2 << std::endl;

         if (d_1 > 0.1)
            if (d_2 > 0.1)
               status = true;

         std::cout << "debug:: rsr set status " << status << std::endl;
         mc.clear_non_drawn_bonds(imol);
      }
   }
   std::cout << "debug:: rsr pre close_molecule(imol) " << std::endl;
   mc.close_molecule(imol);
   std::cout << "debug:: rsr pre close_molecule(imol_map) " << std::endl;
   mc.close_molecule(imol_map);
   return status;
}

int test_rsr_using_residue_range(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-4.mtz"), "FWT", "PHWT", "W", false, false);
   mc.set_imol_refinement_map(imol_map);

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {

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
         mc.set_map_weight(w * 10.0);
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
         std::cout << "DEBUG:: ds: " << d_N_1 << " " << d_N_2 << " " << d_N_3 << std::endl;
         if (d_N_1 < 0.0001)  // no move
            if (d_N_3 < 0.0001) // no move
               if (d_N_2 > 0.08) // move a bit
                  status = 1;
         mc.write_coordinates(imol, "post-refine-using-residue-range.pdb");
      }
   }
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);
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
   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_map_molecule(imol_map)) {
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
         mc.close_molecule(imol_map);
      }
   }
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
      std::cout << "DEBUG:: distance between peptide atoms " << dt << std::endl;
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

   if (false) {
      std::cout << "part_one_done "   << part_one_done   << std::endl;
      std::cout << "part_two_done "   << part_two_done   << std::endl;
      std::cout << "part_three_done " << part_three_done << std::endl;
      std::cout << "part_four_done "  << part_four_done  << std::endl;
   }

   if (part_one_done && part_two_done && part_three_done && part_four_done)
      status = 1;

   mc.close_molecule(imol);
   mc.close_molecule(imol_map);

   return status;
}

// this was converted from Filo's test
int test_add_terminal_residue_v2(molecules_container_t &molecules_container) {

   // this is far from working. Needs new cif writer.

   starting_test(__FUNCTION__);
   int status = 0;
   int coordMolNo = molecules_container.read_pdb(reference_data("./5a3h.mmcif"));
   int mapMolNo = molecules_container.read_mtz(reference_data("./5a3h_sigmaa.mtz"), "FWT", "PHWT", "", false, false);
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
       std::ofstream f1("mmcifString_1");
       std::ofstream f2("mmcifString_2");
       f1 << mmcifString_1 << std::endl;
       f2 << mmcifString_2 << std::endl;
       f1.close();
       f2.close();
   }

   molecules_container.write_coordinates(coordMolNo, "post-addition.cif");
   molecules_container.close_molecule(mapMolNo);
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
   std::cout << "DEBUG:: test_difference_map_contoursh(): " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
             << " triangles" << std::endl;

   for (const auto &vertex : map_mesh.vertices) {
      // std::cout << "vertex " <<  glm::to_string(vertex.pos) << " " << glm::to_string(vertex.color) << std::endl;
      if (vertex.color[0] > vertex.color[1])
         if (vertex.color[0] > vertex.color[2])
            status = 1;
   }
   mc.close_molecule(imol_map);

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
      if (d > 0.9) {

	 mc.import_cif_dictionary(reference_data("NUT.cif"), coot::protein_geometry::IMOL_ENC_ANY);
         // now test an altconf ligand
         int imol_lig = mc.get_monomer("NUT");
         mc.delete_hydrogen_atoms(imol_lig);
         mc.add_alternative_conformation(imol_lig, "//A/1");
         mc.write_coordinates(imol_lig, "NUT-with-alt-conf.pdb");
         coot::atom_spec_t spun_atom_spec("A", 1, "", " C7 ", "A");
         mmdb:: Atom *at_2 = mc.get_atom(imol_lig, spun_atom_spec);
         if (at_2) {
            coot::Cartesian atom_pos_3 = atom_to_cartesian(at_2);
            mc.jed_flip(imol_lig, "//A/1/O1:A", false);
            mc.write_coordinates(imol_lig, "NUT-with-alt-conf-and-jed-flip.pdb");
            coot::Cartesian atom_pos_4 = atom_to_cartesian(at_2);
            dd = coot::Cartesian::lengthsq(atom_pos_3, atom_pos_4);
            d = std::sqrt(dd);
            if (d > 2.0)
               status = 1;
         } else {
            std::cout << "failed to select atom " << spun_atom_spec << std::endl;
         }
      }
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
   if (mc.is_valid_model_molecule(imol)) {
      std::vector<std::string> nst = mc.get_residue_names_with_no_dictionary(imol);
      // weak test
      if (nst.empty())
	 status = 1;
   }
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
   mc.close_molecule(imol_diff_map);
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
   mc.close_molecule(imol_diff_map);
   return status;
}

int test_dictionary_bonds(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("pdb2sar-part.ent"));
   mc.import_cif_dictionary(reference_data("ATP.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("3GP.cif"), imol_1);
   int imol_2 = mc.get_monomer("ATP");
   int imol_3 = mc.read_pdb(reference_data("pdb2sar-part.ent"));

   if (mc.is_valid_model_molecule(imol_2)) {
      if (mc.is_valid_model_molecule(imol_3)) {

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
      }
   }

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
   mc.close_molecule(imol);
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

      std::vector<coot::api::moved_atom_t> moved_atom_positions;
      std::vector<std::string> atom_names = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2"};
      for (const auto &name : atom_names) {
         coot::atom_spec_t as("A", 20, "", name, "");
         moved_atom_positions.push_back(coot::api::moved_atom_t(name, "", 2.1, 3.2, 4.3));
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
   mc.close_molecule(imol);
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

      std::vector<std::string> atom_names = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2"};
      std::vector<coot::api::moved_residue_t> moved_atoms_residue_vec;
      coot::api::moved_residue_t mr("A", 20, "");
      for (const auto &name : atom_names) {
         coot::atom_spec_t as("A", 20, "", name, "");
         mr.add_atom(coot::api::moved_atom_t(name, "", 2.1, 3.2, 4.3));
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
   mc.close_molecule(imol);
   return status;
}



int test_merge_molecules(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   mc.import_cif_dictionary(reference_data("ATP.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("3GP.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("NUT.cif"), coot::protein_geometry::IMOL_ENC_ANY);

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
      std::cout << "DEBUG:: test_read_a_map(): " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
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

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_model_molecule(imol_ligand)) {
         if (mc.is_valid_map_molecule(imol_map)) {
            float n_rmsd = 1.0;
            bool make_conformers = true;
            unsigned int n_conformers = 80;
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
            std::cout << "Eigenvector size std. dev.: " << sd << std::endl;
            if (sd > 0.001) {
               if (solutions.size() < 5) {

                  // tell me about the solutions:
                  for (unsigned int i=0; i<solutions.size(); i++) {
                     const auto &sol(solutions[i]);
                     std::cout << "    Solution " << i << " : "
                               << " volume " << sol.get_cluster_volume() << " "
                               << sol.imol << " "
                               << sol.cluster_idx << " "
                               << sol.ligand_idx << " "
                               << " correl " << sol.get_fitting_score() << " "
                               << std::endl;
                  }
                  status = 1;
               }
            }

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


int test_ligand_fitting_in_map_LZA(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   int imol_ligand = mc.get_monomer("LZA");

   if (mc.is_valid_model_molecule(imol)) {
      if (mc.is_valid_model_molecule(imol_ligand)) {
         if (mc.is_valid_map_molecule(imol_map)) {
            float n_rmsd = 1.0;
            bool make_conformers = false;
            unsigned int n_conformers = 8;
            std::vector<molecules_container_t::fit_ligand_info_t> solutions =
               mc.fit_ligand(imol, imol_map, imol_ligand, n_rmsd, make_conformers, n_conformers);
            std::cout << "DEBUG:: in test_ligand_fitting_in_map_LZA(): found "
                      << solutions.size() << " ligand fitting solutions" << std::endl;

            // tell me about the solutions:
            for (unsigned int i=0; i<solutions.size(); i++) {
               const auto &sol(solutions[i]);
               std::cout << "    LZA Solution " << i << " : "
                         << " volume: " << sol.get_cluster_volume()
                         << " imol: " << sol.imol
                         << " cluster-idx: " << sol.cluster_idx
                         << " ligand-idx: " << sol.ligand_idx
                         << " correl: " << sol.get_fitting_score() << " " << std::endl;
               std::string fn("Ligand-sol-" + coot::util::int_to_string(sol.imol) + ".pdb");
               mc.write_coordinates(sol.imol, fn);
            }
         }
      }
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
   mc.close_molecule(imol_map);
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
         mc.close_molecule(imol_blur);
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
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
   mc.close_molecule(imol_1);
   mc.close_molecule(imol_2);
   mc.close_molecule(imol_map);
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
   mc.close_molecule(imol_map);
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

      coot::simple_mesh_t mesh = mc.get_molecular_representation_mesh(imol, selection, colour, style, CALC_SECONDARY_STRUCTURE);
      if (mesh.vertices.size() > 10) {

         std::cout << "test_molecular_representation() Ribbons OK" << std::endl;

         style = "MolecularSurface";
         coot::simple_mesh_t surface_mesh = mc.get_molecular_representation_mesh(imol, selection, colour, style, CALC_SECONDARY_STRUCTURE);

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

      coot::simple_mesh_t mesh = mc.get_molecular_representation_mesh(imol, selection, colour, style, CALC_SECONDARY_STRUCTURE);
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
   int imol = mc.read_pdb(reference_data("pdb8oie.ent"));
   int imol_map = mc.read_ccp4_map(reference_data("emd_16890.map"), false);
   mc.set_refinement_is_verbose(false);
   mc.import_cif_dictionary(reference_data("CLF.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("HCA.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("S5Q.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("ADP.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data( "MG.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("AF3.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   mc.import_cif_dictionary(reference_data("SF4.cif"), coot::protein_geometry::IMOL_ENC_ANY);

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
   mc.close_molecule(imol_map);
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
                     if (i > 5) continue;
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

int test_instanced_goodsell_style_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   mc.set_use_gemmi(false);
   int imol = mc.read_pdb(reference_data("pdb8ox7.ent"));
   float cwr = 97.0;
   coot::instanced_mesh_t im = mc.get_goodsell_style_mesh_instanced(imol, cwr, 0.8, 0.6);
   std::vector<std::pair<glm::vec4, unsigned int> > ca = colour_analysis(im);

   // check that we have colour variation
   unsigned int n_above = 0;
   unsigned int n_below = 0;
   for (unsigned int i=0; i<im.geom.size(); i++) {
      const coot::instanced_geometry_t &ig = im.geom[i];
      for (unsigned int jj=0; jj<ig.instancing_data_A.size(); jj++) {
         const auto &col =  ig.instancing_data_A[jj].colour;
         if (col.r > 0.8) n_above++;
         if (col.r < 0.5) n_below++;
      }
   }
   if (n_above > 1)
      if (n_below > 1)
         status = 1;
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
      mc.add_colour_rule(imol, "//A", "#66666666");
      coot::simple_mesh_t mesh = mc.get_gaussian_surface(imol, sigma, contour_level, box_radius, grid_scale, b_factor);
      auto colours = colour_analysis(mesh);
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
      coot::instanced_mesh_t im = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, false, false, true, 1);
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
   coot::instanced_mesh_t im_lig = mc.get_bonds_mesh_for_selection_instanced(imol, cid, mode, true, 0.1, 1.0, false, false, true, 1);
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
   coot::instanced_mesh_t im_frag = mc.get_bonds_mesh_instanced(imol, mode, true, 0.2, 1.0, false, false, true, 2);
   colour_analysis(im_frag);

   coot::instanced_mesh_t im = mc.get_bonds_mesh_for_selection_instanced(imol_frag, selection_cid, mode, true, 0.2, 1.0, false, false, true, 2);
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
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
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
      // mc.display_molecule_names_table();
      mc.add_waters(imol, imol_map);

      map_mesh = mc.get_map_contours_mesh(imol_map, 40,40,40, 6, 0.8);
      molecules_container_t::r_factor_stats stats_2 = mc.get_r_factor_stats();
      std::cout << "stats_2: " << mc.r_factor_stats_as_string(stats_2) << std::endl;
      if (stats_2.rail_points_total > 500)
         status = 1;
   }

   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   mc.close_molecule(imol_diff_map);

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

   // what does this test?

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   auto mesh = mc.get_bonds_mesh(imol, mode, true, 0.1, 1.0, 1);
   status = 1;
   return status;

}

int test_delete_side_chain(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec_1("A", 270, "", " O  ","");
      coot::atom_spec_t atom_spec_2("A", 270, "", " OD1","");
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
   auto mesh = mc.get_bonds_mesh_instanced(imol_0, mode, true, 0.1, 1.0, false, false, true, 1);

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

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      std::string crs = "//A/1^#cc0000|//A/2^#cb0002|//A/3^#c00007";
      mc.add_colour_rules_multi(imol, crs);
      auto v = mc.get_colour_rules(imol);
      for (const auto &cr : v) {
         std::cout << " colour rule " << cr.first << " " << cr.second << std::endl;
      }
      if (v.size() == 5)
         status = 1;
      std::cout << "n colour rules: " << v.size() << std::endl;
   }
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
         std::cout << "DEBUG:: test_add_hydrogen_atoms(): pre: " << n_atom_pre << " n_atom_post " << n_atoms_post << std::endl;
         auto mesh_lig_2 = mc.get_bonds_mesh(imol_lig, mode, true, 0.1, 1.0, 1);
         if (n_atoms_post < n_atom_pre)
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
   unsigned int n_residue_per_residue_range = 11;

   if (mc.is_valid_model_molecule(imol)) {
      std::string chain_id = "A";
      std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
                std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> > results =
         mc.mmrrcc(imol, chain_id, imol_map);
      auto mcc = results.first;
      if (mcc.size() > 90)
         status = 1;

      std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
      if (false) // remove noise
         for (it=mcc.begin(); it!=mcc.end(); ++it)
            std::cout << "   " << it->first << " " << it->second.correlation() << std::endl;
   }
   mc.close_molecule(imol);


   if (false) {
      // Filo's example 11729 and 7adk
      imol = mc.read_pdb(reference_data("pdb7adk.ent"));
      imol_map = mc.read_ccp4_map(reference_data("emd_11729.map"), 0);
      auto results = mc.mmrrcc(imol, "B", imol_map);
      auto mcc = results.first;
      auto scc = results.second;
      std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
      std::ofstream f("7adk-b-chain-all-atom.table");
      for (it=mcc.begin(); it!=mcc.end(); ++it)
         f << "   " << it->first << " " << it->second.correlation() << std::endl;
      std::ofstream fs("7adk-b-chain-side-chain.table");
      for (it=scc.begin(); it!=scc.end(); ++it)
         fs << "   " << it->first << " " << it->second.correlation() << std::endl;
      mc.close_molecule(imol);
   }

   mc.close_molecule(imol_map);
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

   mc.close_molecule(imol_map_1);
   mc.close_molecule(imol_map_2);
   return status;
}

int test_auto_read_mtz(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::vector<molecules_container_t::auto_read_mtz_info_t> imol_maps_5a3h
      = mc.auto_read_mtz(reference_data("5a3h_sigmaa.mtz"));

   for (const auto &item : imol_maps_5a3h) {
      std::cout << "    auto-read: map-idx: " << std::setw(2) << item.idx << " F: \"" << item.F << "\" phi: \"" << item.phi
                << "\" Fobs: \"" << item.F_obs << "\" sigFobs: \"" << item.sigF_obs
                << "\" Rfree: \"" << item.Rfree << "\"" << std::endl;
   }

   std::vector<molecules_container_t::auto_read_mtz_info_t> imol_maps
      = mc.auto_read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"));

   // one of these (the last one) should be observed data without an imol
   if (imol_maps.size() == 3) {
      int imol_idx_1 = imol_maps[1].idx;
      int imol_idx_2 = imol_maps[2].idx;
      float rmsd_1 = mc.get_map_rmsd_approx(imol_idx_1);
      float rmsd_2 = mc.get_map_rmsd_approx(imol_idx_2);
      if (mc.is_valid_map_molecule(imol_idx_1)) {
         if (mc.is_valid_map_molecule(imol_idx_1)) {
            std::cout << "test_auto_read_mtz() rmsds " << rmsd_1 << " " << rmsd_2 << std::endl;
            if (rmsd_1 > 0.3) { // test that the FWT map is the first of the pair
               if (rmsd_2 > 0.1) {
                  // what observed data did we find?
                  for (unsigned int i=0; i<imol_maps.size(); i++) {
                     const auto &mtz_info = imol_maps[i];
                     if (! mtz_info.F_obs.empty()) {
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
               }
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

   mc.import_cif_dictionary(reference_data("ATP.cif"), imol_1);
   mc.import_cif_dictionary(reference_data("ATP.cif"), imol_2);
   bool use_rdkit_svg = false;
   std::string bg = "dark-bonds/opaque-bg";
   std::string s = mc.get_svg_for_residue_type(imol_1, "ATP", use_rdkit_svg, bg);

   if (s.length() > 0) {

      std::ofstream f("ATP.svg");
      f << s;
      f.close();
      {
         mc.import_cif_dictionary("G37.cif", coot::protein_geometry::IMOL_ENC_ANY);
         s = mc.get_svg_for_residue_type(imol_1, "G37", use_rdkit_svg, bg);
         std::ofstream f2("G37.svg");
         f2 << s;
         f2.close();
      }

      {
         mc.import_cif_dictionary("GLC.cif", coot::protein_geometry::IMOL_ENC_ANY);
         s = mc.get_svg_for_residue_type(imol_1, "GLC", use_rdkit_svg, bg);
         std::ofstream f2("GLC.svg");
         f2 << s;
         f2.close();
      }

      status = 1;
   }
   return status;

}

int test_superpose(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1-with-gap.pdb"));
   int imol_2 = mc.read_pdb(reference_data("1phk-with-gap.pdb"));

   unsigned int n_pre = mc.get_number_of_molecules();

   if (mc.is_valid_model_molecule(imol_1)) {
      if (mc.is_valid_model_molecule(imol_2)) {

         coot::atom_spec_t atom_spec_1("A", 284, "", " CA ","");
         coot::atom_spec_t atom_spec_2("A", 285, "", " CA ","");
         mmdb::Atom *at_1 = mc.get_atom(imol_1, atom_spec_1);
         mmdb::Atom *at_2 = mc.get_atom(imol_2, atom_spec_2);

         coot::Cartesian atom_pos_1 = atom_to_cartesian(at_1);
         coot::Cartesian atom_pos_2 = atom_to_cartesian(at_2);

         double dd = coot::Cartesian::lengthsq(atom_pos_1, atom_pos_2);
         double d1 = std::sqrt(dd);
         std::cout << "test d1 " << d1 << std::endl;

         // std::pair<std::string, std::string> ss_result_pair = mc.SSM_superpose(imol_1, "A", imol_2, "B");
         // ref ref mov mov
         superpose_results_t ss_results = mc.SSM_superpose(imol_1, "A", imol_2, "A");

         std::cout << "ss_result: info:\n" << ss_results.superpose_info   << std::endl;
         std::cout << "ss_result: alnR\n"  << ss_results.alignment.first  << std::endl;
         std::cout << "ss_result: alnM\n"  << ss_results.alignment.second << std::endl;

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
      }
   }
   std::cout << "done superpose test" << std::endl;
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
   glm::vec3 ca_pos(26.83, 3.43, 31.43);
   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   auto mesh_1 = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, false, false, true, 1);
   mc.add_to_non_drawn_bonds(imol, "//A/270");
   auto mesh_2 = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, false, false, true, 1);

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
      int imol_lig = mc.read_pdb(reference_data("misplaced-moorhen-tutorial-1-ligand.pdb"));
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
         mc.replace_molecule_by_model_from_file(imol, reference_data("3pzt.pdb"));
         mmdb::Atom *at_2 = mc.get_atom(imol, atom_spec);
         if (at_2) {
            std::string rn_2 = at_2->GetResidue()->GetResName();

            if (rn_1 == "ASP")
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
         std::map<unsigned int, std::array<float, 4> > colour_map;
         colour_map[0] = std::array<float, 4> {0.42222222, 0.7, 0.4, 1.0};
         colour_map[2] = std::array<float, 4> {0.42222222, 0.4, 0.7, 1.0};
         colour_map[1] = std::array<float, 4> {0.7, 0.4, 0.42222222, 1.0};
         std::vector<std::pair<std::string, unsigned int> > indexed_residues_cids;
         indexed_residues_cids.push_back(std::make_pair("//A",2));
         indexed_residues_cids.push_back(std::make_pair("//A/100-200",1));
         indexed_residues_cids.push_back(std::make_pair("//A/130-150",0));
         std::string mode("USER-DEFINED-COLOURS");
         mc.set_user_defined_bond_colours(imol, colour_map);
         bool colour_applies_to_non_carbon_atoms_also = true;
         mc.set_user_defined_atom_colour_by_selection(imol, indexed_residues_cids, colour_applies_to_non_carbon_atoms_also);
         coot::instanced_mesh_t im = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, false, false, true, 1);
         if (im.geom.size() > 3) {
            if (im.geom[0].instancing_data_A.size() > 1000)
               status = 1;
            const auto& g = im.geom[1].instancing_data_B;
            std::cout << "........... g.size() " << g.size() << std::endl;
            for (unsigned int j=0; j<g.size(); j++) {
               if (j > 3) continue;
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
   mc.close_molecule(imol_map);

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
      std::cout << "  "                             << h_bonds[i].donor.res_no    << " " << h_bonds[i].donor.name
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

   int import_status = mc.import_cif_dictionary(reference_data("LZA.cif"), -999999);
   if (import_status == 0)
      return status;
   int imol = mc.get_monomer("LZA");

   if (mc.is_valid_model_molecule(imol)) {
      coot::colour_t col(0.0999, 0.0888, 0.0777);
      mc.set_use_bespoke_carbon_atom_colour(imol, true);
      mc.set_bespoke_carbon_atom_colour(imol, col);
      std::string mode("VDW-BALLS");
      coot::instanced_mesh_t im = mc.get_bonds_mesh_instanced(imol, mode, true, 0.1, 1.0, false, false, true, 1);

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

   int import_status = mc.import_cif_dictionary(reference_data("LZA.cif"), -999999);
   if (import_status == 0)
      return status;
   int imol = mc.get_monomer("LZA");
   if (mc.is_valid_model_molecule(imol)) {
      std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";
      auto mesh_light = mc.get_bonds_mesh_instanced(imol, mode, false, 0.2, 1.0, false, false, true, 1);
      auto mesh_dark  = mc.get_bonds_mesh_instanced(imol, mode, true,  0.2, 1.0, false, false, true, 1);
      std::cout << "starting colour analysis for mesh_light" << std::endl;
      colour_analysis(mesh_light);
      std::cout << "starting colour analysis for mesh_dark" << std::endl;
      colour_analysis(mesh_dark);
      std::cout << "done analyses" << std::endl;
      status = 1; // this should be a better test
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

   mc.close_molecule(imol_map);
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
      if (mci.updated_centre.z() > 10.0)
         status = 1;
   }
   mc.close_molecule(imol_map);

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
   mc.close_molecule(imol_map);
   return status;
}

int test_bucca_ml_growing(molecules_container_t &mc) {

   // ignore this - for now at least

   starting_test(__FUNCTION__);
   int status = 0;

   // int imol     = mc.read_pdb(reference_data("1gwd-large-C-terminal-fragment-missing.pdb"));
   std::string fn = "1gwd-large-C-terminal-fragment-missing.pdb";
   fn = "1gwd-118-chop.pdb";
   fn = "1gwd-2-chop.pdb";
   int imol     = mc.read_pdb(reference_data(fn));
   int imol_map = mc.read_mtz(reference_data("1gwd_map.mtz"), "FWT", "PHWT", "FOM", false, false);

   mc.set_refinement_is_verbose(true);
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
   mc.set_refinement_is_verbose(false);
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
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

   std::map<unsigned int, std::array<float, 4> > colour_index_map;
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

   auto bonds = mc.get_bonds_mesh_instanced(imol, mode, false, 0.2, 1.0, false, false, true, 1);

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

   coot::instanced_mesh_t im = mc.get_bonds_mesh_for_selection_instanced(imol, "//A/1-3", "VDW-BALLS", false, 0.1, 1.0, false, false, true, 1);
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

   auto is_near_colour = [] (const glm::vec4 &col_1, const std::array<float, 4> &col_2) {
      float cf = 0.04;
      if (std::fabs(col_2[0] - col_1.r) < cf)
         if (std::fabs(col_2[1] - col_1.g) < cf)
            if (std::fabs(col_2[2] - col_1.b) < cf)
               return true;
      return false;
   };


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
      std::map<unsigned int, std::array<float, 4> > colour_map;
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
      auto bonds = mc.get_bonds_mesh_for_selection_instanced(imol, "/", mode, false, 0.2, 1.0, false, false, true, 1);
      auto &geom = bonds.geom;
      auto ca = get_colour_analysis(bonds);

      // colour 53 supercedes/replaces the others
      //
      bool col_51 = false;
      bool col_52 = false;
      bool col_53 = false;
      for (unsigned int i=0; i<ca.size(); i++) {
         const auto &ca_row = ca[i];
         if (is_near_colour(ca_row.col, colour_map[51])) col_51 = true;
         if (is_near_colour(ca_row.col, colour_map[52])) col_52 = true;
         if (is_near_colour(ca_row.col, colour_map[53])) col_53 = true;
      }
      if (col_51 == false)
         if (col_52 == false)
            if (col_53 == true)
               status = true;
   }

   return status;
}



int test_is_em_map(molecules_container_t &mc) {

   auto close_float = [] (float a, float b) {
      return fabsf(a - b) < 0.001;
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map = mc.read_ccp4_map(reference_data("emd_32143.map"), 0);
   bool is_EM_map = mc.is_EM_map(imol_map);
   float cl = mc.get_suggested_initial_contour_level(imol_map);
   float rmsd = mc.get_map_rmsd_approx(imol_map);

   std::cout << "test_is_em_map(): EM: " << is_EM_map << " cl " << cl << " with rmsd: " << rmsd
             << " ratio " << cl/rmsd << std::endl;

   if (is_EM_map) {
      if (close_float(cl, 4.0 * rmsd)) {
         coot::util::map_molecule_centre_info_t mci = mc.get_map_molecule_centre(imol_map);
         std::cout << "mci.suggested_radius " << mci.suggested_radius << std::endl;
         if (mci.suggested_radius > 30.0)
            status = 1;
      }
   }
   mc.close_molecule(imol_map);
   return status;
}

int test_other_user_defined_colours_other(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         std::map<unsigned int, std::array<float, 4> > colour_index_map;
         colour_index_map[21] = {1.11111111, 1.111111, 0};
         std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
         mc.set_user_defined_bond_colours(imol, colour_index_map);
         std::vector<std::pair<std::string, unsigned int> > indexed_cids;
         indexed_cids.push_back(std::make_pair("//A/1-5", 21));
         bool non_carbon_atoms_also_flag = false;
         auto bonds_1 = mc.get_bonds_mesh_for_selection_instanced(imol, "/", mode, false, 0.2, 1.0, false, false, true, 1);
         auto bonds_2 = mc.get_bonds_mesh_instanced(imol, mode, false, 0.2, 1.0, false, false, true, 1);
         mc.set_user_defined_atom_colour_by_selection(imol, indexed_cids, non_carbon_atoms_also_flag);
         auto bonds_3 = mc.get_bonds_mesh_for_selection_instanced(imol, "/", mode, false, 0.2, 1.0, false, false, true, 1);
         auto &geom_1 = bonds_1.geom;
         auto &geom_3 = bonds_3.geom;

         auto &vb_1 = geom_1[1].instancing_data_B;
         auto &vb_3 = geom_3[1].instancing_data_B;

         for (unsigned int i=0; i<vb_1.size(); i++) {
             if (i > 5) continue;
             auto col = vb_1[i].colour;
             std::cout << "instancing colour_1: " << i << " " << glm::to_string(col) << "\n";
         }

         for (unsigned int i=0; i<vb_3.size(); i++) {
            if (i > 5) continue;
            auto col = vb_3[i].colour;
            std::cout << "instancing colour_3: " << i << " " << glm::to_string(col) << "\n";
         }
         std::vector<colour_analysis_row> ca_1 = get_colour_analysis(bonds_1);
         std::vector<colour_analysis_row> ca_3 = get_colour_analysis(bonds_3);
         // different vec indices because UD colour becomes the first colour.
         // This is weird. When run in "single test only" the indices need to
         // be 4 and 3. It might have something to do with reading the ATP at the start.
         //
         if (ca_3[3].count == (ca_1[2].count - 85))
            status = 1;

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

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_2 = mc.read_pdb(reference_data("3pzt.pdb"));
   if (mc.is_valid_model_molecule(imol_1)) {
      if (mc.is_valid_model_molecule(imol_2)) {
         // it's actually for moleecule 4 (1 was renamed to 4)
         int n_extra =
	    mc.read_extra_restraints(imol_1,
				     reference_data("moorhen-tutorial-structure-number-1-prosmart.txt"));
	 std::cout << "test_read_extra_restraints made " << n_extra << " extra restraints" << std::endl;
	 if (n_extra > 0) {
	    coot::instanced_mesh_t im = mc.get_extra_restraints_mesh(imol_1, 0);
	    if (! im.geom.empty()) {
	       std::cout << "instancing_data_B size " << im.geom[0].instancing_data_B.size() << std::endl;
	       if (im.geom[0].instancing_data_B.size() > 10)
		  status = 1;
	    } else {
	       std::cout << "ERROR:: im geom is empty" << std::endl;
	    }
	 }
      }
   }
   return status;
}




int test_colour_map_by_other_map(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_map_1 = mc.read_ccp4_map(reference_data("emd_16890.map"), false);
   int imol_map_2 = mc.read_ccp4_map(reference_data("scale_res_emd_16890.mrc"), false);
   if (mc.is_valid_map_molecule(imol_map_1)) {
      if (mc.is_valid_map_molecule(imol_map_2)) {
         coot::simple_mesh_t mesh = mc.get_map_contours_mesh_using_other_map_for_colours(imol_map_1, imol_map_2,
                                                                                         160, 160, 160,
                                                                                         100, 0.16,
                                                                                         0.3, 0.9, false);
         std::cout << "test: mesh v and t: " << mesh.vandt() << std::endl;
         // colour_analysis(mesh);
         if (mesh.vertices.size() > 1000) status = true;
      }
   }
   mc.close_molecule(imol_map_1);
   mc.close_molecule(imol_map_2);
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
   std::string bg = "dark-bonds/opaque-bg";
   std::string svg = mc.get_svg_for_residue_type(coot::protein_geometry::IMOL_ENC_ANY, "MOI", use_rdkit_rendering, bg);
   std::ofstream f("MOI.svg");
   f << svg;
   f.close();
   if (svg.length() > 100)
      status = 1;
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
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);

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
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);

   return status;

}

int test_mmcif_atom_selection(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::string fn = reference_data("1ej6-assembly1.cif");
   std::cout << "reading " << fn << std::endl;
   int imol = mc.read_pdb(reference_data(fn));
   mmdb::Manager *mol = mc.get_mol(imol); // testing get_mol()
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
      mol->DeleteSelection(selHnd_3);
      mol->DeleteSelection(selHnd_2);
      mol->DeleteSelection(selHnd_1);
      std::cout << "Looked for 100 atoms and found " << n_matcher << " matchers" << std::endl;
      if (n_matcher == 0) status = 1;
   }
   mc.close_molecule(imol);
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

   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   return status;
}

int test_test_the_threading(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 1; // no faiiure

   for (unsigned int i=0; i<20; i++) { // was 50, but that caused the process to be Killed
      double r = mc.test_the_threading(i);
      std::cout << " test_threading: " << i << " " << r << std::endl;
   }

   return status;
}

int test_thread_launching(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
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

   starting_test(__FUNCTION__);
   int status = 1; // no failure
   for (unsigned int j=1; j<20; j++) {
      double t = mc.test_thread_pool_threads(j);
      std::cout << " launching " << j << " " << t << " micro-seconds " << std::endl;
   }
   return status;
}

int test_long_name_ligand_cif_merge(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   mc.set_use_gemmi(true);
   int imol = mc.read_pdb(reference_data("8a2q.cif"));
   mc.import_cif_dictionary(reference_data("7ZTVU.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   int imol_lig = mc.get_monomer("7ZTVU");
   std::string sl = std::to_string(imol_lig);
   std::pair<int, std::vector<merge_molecule_results_info_t> > ss = mc.merge_molecules(imol, sl);
   mc.write_coordinates(imol, "8a2q-with-new-ligand.cif");
   std::ifstream f("8a2q-with-new-ligand.cif");
   std::string line;
   std::string ss1("7ZTVU");
   std::string ss2("7ZTVU");
   std::string ss3("C10");
   bool found_it = false;
   while (std::getline(f, line)) {
      if (line.find(ss1) != std::string::npos) {
         if (line.find(ss2) != std::string::npos) {
            if (line.find(ss3) != std::string::npos) {
               found_it = true;
            }
         }
      }
   }
   if (found_it) status = 1;
   return status;
}
#undef USE_GEMMI
#ifdef USE_GEMMI
#include "gemmi/mmread.hpp"
#include "gemmi/mmdb.hpp"

int test_disappearing_ligand(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   // int imol = mc.read_pdb(reference_data("6ttq.cif")); // needs gemmi
   // mc.set_use_gemmi(true);
   int imol = mc.read_coordinates(reference_data("8a2q.cif"));
   if (mc.is_valid_model_molecule(imol)) {
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
      gemmi::Structure structure = gemmi::read_structure_file("merged.cif");
      std::string target_residue_name = "MOI";
      for (auto& model : structure.models) {
         for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
               if (residue.name == target_residue_name) {
                    std::cout << "Found residue " << target_residue_name << " in chain " << chain.name << std::endl;
                    status = 1;
                    break;
               }
            }
         }
      }
   } else {
      std::cout << "Failed to correctly read 8a2q.cif" << std::endl;
   }

   return status;
}
#endif

// 20240205-PE this is not a good name for this test. The failure is not in merging, the failure
// is in writing out the cif file.
//
// The fix for this is to use gemmi for the output
//
// This fails - waiting for gemmi-based fix.
//
int test_ligand_merge(molecules_container_t &mc) {

   auto test_mmdb = [] () {
      // int imol_2 = mc.read_pdb(reference_data("2vtq.cif"));
      // mc.write_coordinates(imol_2, "2vtq-just-input-output.cif");
      mmdb::Manager *mol = new mmdb::Manager;
      mol->ReadCoorFile(reference_data("2vtq.cif").c_str());
      mol->WriteCIFASCII("2vtq-input-output-pure-mmdb.cif");
      delete mol;
   };

   starting_test(__FUNCTION__);
   int status = 0;

   test_mmdb();
   int imol = mc.read_pdb(reference_data("2vtq-sans-ligand.cif"));
   mc.write_coordinates(imol, "2vtq-sans-ligand-just-input-output.cif");
   if (mc.is_valid_model_molecule(imol)) {
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
   map_mesh.export_to_gltf("map-around-ligand.glb", 0.5, 0.5, true);

   std::cout << "-------------------------------------------------- ligand mesh " << std::endl;

   std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
   int import_status = mc.import_cif_dictionary(reference_data("LZA.cif"), -999999);
   if (import_status == 0)
      return status;
   int imol_lig = mc.get_monomer("LZA");
   int imol_frag = mc.copy_fragment_using_cid(imol, "//A/1299");
   std::cout << "test_gltf_export() imol_frag " << imol_frag << std::endl;
   coot::instanced_mesh_t im    = mc.get_bonds_mesh_instanced(imol_frag, mode, true, 0.1, 1.0, false, false, true, 1);
   coot::simple_mesh_t sm_lig = coot::instanced_mesh_to_simple_mesh(im);
   sm_lig.export_to_gltf("lig.glb", 0.5, 0.5, true);

   std::cout << "-------------------------------------------------- neighbour mesh " << std::endl;
   std::vector<coot::residue_spec_t> neighbs = mc.get_residues_near_residue(imol, "//A/1299", 4.2);
   std::string multi_cid = make_multi_cid(neighbs);
   mc.set_draw_missing_residue_loops(false);
   coot::instanced_mesh_t im_neighbs = mc.get_bonds_mesh_for_selection_instanced(imol, multi_cid, mode, true, 0.15, 1.0, false, false, true, 1);
   coot::simple_mesh_t sm_neighbs = coot::instanced_mesh_to_simple_mesh(im_neighbs);
   sm_neighbs.export_to_gltf("neighbs.glb", 0.5f, 0.5f, true);

   struct stat buf_1;
   int istat_1 = stat("lig.glb", &buf_1);
   if (istat_1 == 0) {
      if (buf_1.st_size > 100000) {
         struct stat buf_2;
         int istat_2 = stat("neighbs.glb", &buf_2);
         if (istat_2 == 0) {
            if (buf_2.st_size > 100000) {
               status = 1;
            }
         }
      }
   }

   mc.close_molecule(imol_map);
   mc.close_molecule(imol);

   return status;
}


int test_gltf_export_via_api(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   mc.set_use_gemmi(true); // 20240727-PE there seems to be a memory problem when using gemmi atm
                           // so for now, let's not use gemmi for the tests.

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
   std::cout << "stat for buf_1: " << istat_1 << std::endl;
   if (istat_1 == 0) {
      std::cout << "buf_1 size: " << buf_1.st_size << std::endl;
      // In the github action this is 856136. On my PC it's more than 1,000,000
      // I don't know if the github action version is broken. Maybe
      // related to map sampling?
      if (buf_1.st_size > 800000) {

         struct stat buf_2;
         int istat_2 = stat("fat-ligand.glb", &buf_2);
         std::cout << "stat for buf_2: " << istat_2 << std::endl;
         if (istat_2 == 0) {
            std::cout << "buf_2 size: " << buf_2.st_size << std::endl;
            if (buf_2.st_size > 100000) {
               status = 1;
            }
         }
      }
   }
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);
   return status;
}


int test_5char_ligand_merge(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
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
   mc.close_molecule(imol_map);
   mc.close_molecule(imol);
   return status;
}

int test_molecule_diameter(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_2 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-4.pdb"));
   int imol_3 = mc.get_monomer("GOL");
   float d1 = mc.get_molecule_diameter(imol_1);
   float d2 = mc.get_molecule_diameter(imol_2);
   float d3 = mc.get_molecule_diameter(imol_3);
   if (d1 > 30.0)
      if (d1 < 300.0)
         if (d2 > 30.0)
            if (d2 < 300.0)
               if (d3 > 3.0)
                  if (d3 < 30.0)
               status = 1;

   return status;
}

int test_B_factor_multiply(molecules_container_t &mc) {

   auto close_float = [] (float a, float b) {
      return fabsf(a - b) < 0.001;
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.get_monomer("GOL");

   std::cout << "test_B_factor_multiply imol " << imol << std::endl;

   coot::residue_spec_t rs("A", 1, "");
   mmdb::Residue *residue_p = mc.get_residue(imol, rs);
   if (residue_p) {
      std::vector<float> B_pre;
      std::vector<float> B_post;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            B_pre.push_back(at->tempFactor);
         }
      }

      mc.multiply_residue_temperature_factors(imol, "//", 2.0);

      residue_p = mc.get_residue(imol, rs);
      if (residue_p) {
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               B_post.push_back(at->tempFactor);
            }
         }
      }

      if (B_pre.size() > 12) {
         if (B_pre.size() == B_post.size()) {
            bool clean = true;
            for (unsigned int iat=0; iat<B_pre.size(); iat++) {
               if (close_float(B_pre[iat] * 2.0, B_post[iat])) {
               } else {
                  std::cout << "fail for iat " << iat << " " << B_pre[iat] << " " << B_post[iat] << std::endl;
                  clean = false;
               }
            }
            if (clean) status = 1;
         }
      }
   }

   return status;
}

int test_change_chain_id(molecules_container_t &mc) {

   auto get_min_max_in_chain = [] (mmdb::Manager *mol, const std::string &chain_id_in) {

      int min_res_no_J =  9999;
      int max_res_no_J = -9999;

      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               std::string chain_id(chain_p->GetChainID());
               if (chain_id == chain_id_in) {
                  int n_res = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<n_res; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        int res_no = residue_p->GetSeqNum();
                        if (res_no < min_res_no_J) min_res_no_J = res_no;
                        if (res_no > max_res_no_J) max_res_no_J = res_no;
                     }
                  }
               }
            }
         }
      }
      return std::pair<int, int>(min_res_no_J, max_res_no_J);
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (!mc.is_valid_model_molecule(imol)) return 0;

   std::pair<int, std::string> r_1 = mc.change_chain_id(imol, "A", "J", true, 38, 42); // non-existant
   std::cout << "change_chain_id result-1: " << r_1.first << " \"" << r_1.second << "\"" << std::endl;
   std::pair<int, std::string> r_2 = mc.change_chain_id(imol, "A", "J", true,  2, 22); // exists
   std::cout << "change_chain_id result-2: " << r_2.first << " \"" << r_2.second << "\"" << std::endl;

   int min_res_no_J = 9999;
   int max_res_no_J = 0;

   mmdb::Manager *mol = mc.get_mol(imol);
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id(chain_p->GetChainID());
            if (chain_id == "J") {
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     if (res_no < min_res_no_J) min_res_no_J = res_no;
                     if (res_no > max_res_no_J) max_res_no_J = res_no;
                  }
               }
            }
         }
      }
   }
   std::cout << "min_res_no_J " << min_res_no_J << std::endl;
   std::cout << "max_res_no_J " << max_res_no_J << std::endl;

   // now it in again and try to change A to B. It should fail
   imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::pair<int, std::string> r_3 = mc.change_chain_id(imol, "A", "B", false, -1, -1); // crash
   std::cout << "change_chain_id result-3: " << r_3.first << " \"" << r_3.second << "\"" << std::endl;
   // moving to the C chain is fine though
   std::pair<int, std::string> r_4 = mc.change_chain_id(imol, "A", "C", false, -1, -1); // OK
   std::cout << "change_chain_id result-4: " << r_4.first << " \"" << r_4.second << "\"" << std::endl;

   mol = mc.get_mol(imol);
   std::pair<int,int> C_min_max = get_min_max_in_chain(mol, "C");
   // std::cout << "C_min_max " << C_min_max.first << " " << C_min_max.second << std::endl;

   if (r_1.first == 0) {
      if (r_2.first == 1) {
         if (min_res_no_J == 2) {
            if (max_res_no_J == 22) {
               if (r_3.first == 0) {
                  if (r_4.first == 1) {
                     if (C_min_max.first == 1) {
                        if (C_min_max.second == 298) {
                           status = 1;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return status;
}

int test_non_drawn_CA_bonds(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   if (mc.is_valid_model_molecule(imol)) {
      int imol_frag = mc.copy_fragment_using_cid(imol, "//A/101-111");
      mc.add_to_non_drawn_bonds(imol_frag, "//A/103-111");
      std::string mode = "CA+LIGANDS";
      auto bonds = mc.get_bonds_mesh_for_selection_instanced(imol_frag, "//A", mode, false, 0.2, 1.0, false, false, true, 1);
      auto &geom = bonds.geom;
      // should be size 2 of course, if we don't add the range to the non-drawn bond
      // not 4
      std::cout << ":::::::::::::::::::::::: bonds geom was of size " << geom.size() << std::endl;

      if (geom.empty()) {
         std::cout << "geom empty" << std::endl;
      } else {
         if (geom.size() == 1) {
            const std::vector<coot::instancing_data_type_B_t> &idB = geom[0].instancing_data_B;
            std::cout << "idB size " << idB.size() << std::endl;
            // print the instancing data here. You should see a duplicate/reverse
            status = 1;
         }
      }
   }

   return  status;
}

int test_17257(molecules_container_t &mc) {
   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map = mc.read_ccp4_map(reference_data("emd_17257.map.gz"), false);
   std::cout << "imol_map is " << imol_map << std::endl;
   if (mc.is_valid_map_molecule(imol_map)) {
      status = 1;
   }
   // however, until there is a gzip map reader for CCP4, we expect the
   // return value to be -3:
   if (imol_map == -3) {
      status = 1;
   }
   mc.close_molecule(imol_map);
   return status;
}

int test_shiftfield_b_factor_refinement(molecules_container_t &mc) {

   auto get_average_b_factor = [] (mmdb::Residue *residue_p) {
      float sum = 0.0;
      if (!residue_p) {
         std::cout << "Null residue " << std::endl;
         return 0.0f;
      }
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      int count = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            sum += at->tempFactor;
            count++;
         }
      }
      return sum/static_cast<float>(count);
   };

   starting_test(__FUNCTION__);
   int status = 0;

   mc.set_use_gemmi(true); // 20240211-PE crash if set_use_gemmi(true) (the default).
   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   int imol_diff_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "DELFWT", "PHDELWT", "W", false, true);
   mc.associate_data_mtz_file_with_map(imol_map, reference_data("moorhen-tutorial-map-number-1.mtz"), "FP", "SIGFP", "FREE");

   if (mc.is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec_1("A", 10, "");
      coot::residue_spec_t res_spec_2("A", 66, "");
      mmdb::Residue *r_1 = mc.get_residue(imol, res_spec_1);
      mmdb::Residue *r_2 = mc.get_residue(imol, res_spec_2);
      // set up some weird B-factors on some atoms
      float b_orig_1 = get_average_b_factor(r_1);
      float b_orig_2 = get_average_b_factor(r_2);
      mc.multiply_residue_temperature_factors(imol, "//A/3-13",  2.0);
      mc.multiply_residue_temperature_factors(imol, "//A/63-69", 0.2);
      float b_pre_1 = get_average_b_factor(r_1);
      float b_pre_2 = get_average_b_factor(r_2);
      bool shiftfield_status = mc.shift_field_b_factor_refinement(imol, imol_map);
      if (shiftfield_status) {
         auto stats = mc.sfcalc_genmaps_using_bulk_solvent(imol, imol_map, imol_diff_map, imol_map);
         std::cout << "DEBUG:: in test_shiftfield_b_factor_refinement() with r-factor " << stats.r_factor << std::endl;
         float b_post_1 = get_average_b_factor(r_1);
         float b_post_2 = get_average_b_factor(r_2);
         std::cout << "B-factors: orig " << b_orig_1 << " " << b_orig_2
                   << " " << b_pre_1 << " " << b_pre_2 << " post " << b_post_1 << " " << b_post_2 << std::endl;
         if (b_post_1 < 66.0)
            if (b_post_2 > 12.0)
               status = 1;
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   mc.close_molecule(imol_diff_map);
   return status;
}

int test_split_model(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("203d.pdb"));
   std::vector<int> new_mol_indices = mc.split_multi_model_molecule(imol);
   std::cout << "new_mol_indices was of size " << new_mol_indices.size() << std::endl;
   if (new_mol_indices.size() == 40) {
      unsigned int n_models = 0;
      for (int i : new_mol_indices) {
         mmdb::Manager *mol = mc.get_mol(i);
         if (mol) {
            mmdb::Model *model_p = mol->GetModel(1);
            // std::cout << "MODEL 1 for molecule " << i << " " << model_p << std::endl;
            if (model_p) n_models++;
         }
      }
      std::cout << "DEBUG:: in test_split_model() n_models is " << n_models << std::endl;
      if (n_models == 40) status = 1;
   } else {
      std::cout << "DEBUG:: in test_split_model() new_models size was " << new_mol_indices.size() << std::endl;
   }

   for (const auto &idx : new_mol_indices)
      mc.close_molecule(idx);

   return status;
}

int test_copy_molecule_memory_leak(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   mc.set_use_gemmi(true);
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   const unsigned int n_new_mols = 200;

   if (mc.is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = mc.get_mol(imol);
      std::vector<mmdb::Manager *> mol_copies;
      for (unsigned int i=0; i<n_new_mols; i++) {
         mmdb::Manager *mol_copy = coot::util::copy_molecule(mol);
         mol_copies.push_back(mol_copy);
      }
      for (unsigned int i=0; i<n_new_mols; i++) {
         delete mol_copies[i];
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(500));

      // now create coot::molecule_t from that:
      std::vector<int> new_molecule_vec;
      for (unsigned int i=0; i<n_new_mols; i++) {
         int imol_new = mc.copy_fragment_using_cid(imol, "/");
         new_molecule_vec.push_back(imol_new);
      }
      for (unsigned int i=0; i<n_new_mols; i++) {
         mc.close_molecule(new_molecule_vec[i]);
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
      status = 1;
   }
   mc.close_molecule(imol);
   return status;
}

int test_make_ensemble(molecules_container_t &mc) {

   auto make_molecule_string_list = [] (const std::vector<int> &mols) {
      std::string s;
      for (unsigned int i=0; i<mols.size(); i++) {
         s += std::to_string(mols[i]);
         if (i<(mols.size()-1)) s += ":";
      }
      return s;
   };

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("203d.pdb"));
   std::vector<int> new_mol_indices = mc.split_multi_model_molecule(imol);
   std::cout << "new_mol_indices was of size " << new_mol_indices.size() << std::endl;
   if (new_mol_indices.size() == 40) {
      std::string s = make_molecule_string_list(new_mol_indices);
      std::cout << "s: " << s << std::endl;
      int imol_new = mc.make_ensemble(s);
      if (mc.is_valid_model_molecule(imol_new)) {
         status = 1;
         mc.write_coordinates(imol_new, "ensemble.pdb");
      }
      mc.close_molecule(imol_new);
   }
   mc.close_molecule(imol);
   for (const auto &idx : new_mol_indices)
      mc.close_molecule(idx);
   mc.end_delete_closed_molecules();
   return status;
}

int test_ligand_torsions(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_lig = mc.read_pdb(reference_data("LZA-wiggled.pdb"));
   int imol_ref = mc.get_monomer("LZA");

   if (mc.is_valid_model_molecule(imol_lig)) {
      if (mc.is_valid_model_molecule(imol_ref)) {

         mc.match_ligand_torsions_and_position_using_cid(imol_lig, imol_ref, "//A/1");
         mc.write_coordinates(imol_ref, "LZA-ref.pdb");
         mc.write_coordinates(imol_lig, "LZA-unwiggled.pdb");

         // now a different molecule:
         mc.import_cif_dictionary(reference_data("acedrg-LZB.cif"), coot::protein_geometry::IMOL_ENC_ANY);
         int imol_lzb = mc.get_monomer("LZB");
         mc.match_ligand_torsions_and_position_using_cid(imol_lzb, imol_ref, "//A/1");
         mc.write_coordinates(imol_lzb, "LZB-fit.pdb");

         // test atom position difference here
         mmdb::Atom *at_ref_F1 = mc.get_atom(imol_ref, coot::atom_spec_t("A", 1, "", " F19", ""));
         mmdb::Atom *at_ref_F2 = mc.get_atom(imol_ref, coot::atom_spec_t("A", 1, "", " F24", ""));
         mmdb::Atom *at_ref_N1 = mc.get_atom(imol_ref, coot::atom_spec_t("A", 1, "", " N11", ""));
         mmdb::Atom *at_ref_N2 = mc.get_atom(imol_ref, coot::atom_spec_t("A", 1, "", " N6 ", ""));

         mmdb::Atom *at_lza_F1 = mc.get_atom(imol_lig, coot::atom_spec_t("A", 1, "", " F19", ""));
         mmdb::Atom *at_lza_F2 = mc.get_atom(imol_lig, coot::atom_spec_t("A", 1, "", " F24", ""));
         mmdb::Atom *at_lza_N1 = mc.get_atom(imol_lig, coot::atom_spec_t("A", 1, "", " N11", ""));
         mmdb::Atom *at_lza_N2 = mc.get_atom(imol_lig, coot::atom_spec_t("A", 1, "", " N6 ", ""));

         mmdb::Atom *at_lzb_F1 = mc.get_atom(imol_lzb, coot::atom_spec_t("A", 1, "", " F1 ", ""));
         mmdb::Atom *at_lzb_F2 = mc.get_atom(imol_lzb, coot::atom_spec_t("A", 1, "", " F2 ", ""));
         mmdb::Atom *at_lzb_N1 = mc.get_atom(imol_lzb, coot::atom_spec_t("A", 1, "", " N4 ", ""));
         mmdb::Atom *at_lzb_N2 = mc.get_atom(imol_lzb, coot::atom_spec_t("A", 1, "", " N1 ", ""));

         std::cout << "at_ref_F1 " << at_ref_F1 << std::endl;
         std::cout << "at_ref_F2 " << at_ref_F2 << std::endl;
         std::cout << "at_ref_N1 " << at_ref_N1 << std::endl;
         std::cout << "at_ref_N2 " << at_ref_N2 << std::endl;

         std::cout << "at_lza_F1 " << at_lza_F1 << std::endl;
         std::cout << "at_lza_F2 " << at_lza_F2 << std::endl;
         std::cout << "at_lza_N1 " << at_lza_N1 << std::endl;
         std::cout << "at_lza_N2 " << at_lza_N2 << std::endl;

         std::cout << "at_lzb_F1 " << at_lzb_F1 << std::endl;
         std::cout << "at_lzb_F2 " << at_lzb_F2 << std::endl;
         std::cout << "at_lzb_N1 " << at_lzb_N1 << std::endl;
         std::cout << "at_lzb_N2 " << at_lzb_N2 << std::endl;

         clipper::Coord_orth pos_at_ref_F1 = coot::co(at_ref_F1);
         clipper::Coord_orth pos_at_ref_F2 = coot::co(at_ref_F2);
         clipper::Coord_orth pos_at_ref_N1 = coot::co(at_ref_N1);
         clipper::Coord_orth pos_at_ref_N2 = coot::co(at_ref_N2);

         clipper::Coord_orth pos_at_lza_F1 = coot::co(at_lza_F1);
         clipper::Coord_orth pos_at_lza_F2 = coot::co(at_lza_F2);
         clipper::Coord_orth pos_at_lza_N1 = coot::co(at_lza_N1);
         clipper::Coord_orth pos_at_lza_N2 = coot::co(at_lza_N2);

         clipper::Coord_orth pos_at_lzb_F1 = coot::co(at_lzb_F1);
         clipper::Coord_orth pos_at_lzb_F2 = coot::co(at_lzb_F2);
         clipper::Coord_orth pos_at_lzb_N1 = coot::co(at_lzb_N1);
         clipper::Coord_orth pos_at_lzb_N2 = coot::co(at_lzb_N2);

         double d_lza_F1 = std::sqrt((pos_at_lza_F1-pos_at_ref_F1).lengthsq());
         double d_lza_F2 = std::sqrt((pos_at_lza_F2-pos_at_ref_F2).lengthsq());
         double d_lza_N1 = std::sqrt((pos_at_lza_N1-pos_at_ref_N1).lengthsq());
         double d_lza_N2 = std::sqrt((pos_at_lza_N2-pos_at_ref_N2).lengthsq());

         double d_lzb_F1 = std::sqrt((pos_at_lzb_F1-pos_at_ref_F1).lengthsq());
         double d_lzb_F2 = std::sqrt((pos_at_lzb_F2-pos_at_ref_F2).lengthsq());
         double d_lzb_N1 = std::sqrt((pos_at_lzb_N1-pos_at_ref_N1).lengthsq());
         double d_lzb_N2 = std::sqrt((pos_at_lzb_N2-pos_at_ref_N2).lengthsq());

         // the piperidine rings don't superpose correctly sadly. so these number
         // are bigger than I'd like them to be
         //
         double dc = 0.7;

         std::cout << "d_lza_F1 " << d_lza_F1 << std::endl;
         std::cout << "d_lza_F2 " << d_lza_F2 << std::endl;
         std::cout << "d_lza_N1 " << d_lza_N1 << std::endl;
         std::cout << "d_lza_N2 " << d_lza_N2 << std::endl;

         std::cout << "d_lzb_F1 " << d_lzb_F1 << std::endl;
         std::cout << "d_lzb_F2 " << d_lzb_F2 << std::endl;
         std::cout << "d_lzb_N1 " << d_lzb_N1 << std::endl;
         std::cout << "d_lzb_N2 " << d_lzb_N2 << std::endl;

         if ((d_lza_F1 < dc) && (d_lza_F2 < dc))
            if ((d_lza_N1 < dc) && (d_lza_N2 < dc))
               if ((d_lzb_F1 < dc) && (d_lzb_F2 < dc))
                  if ((d_lzb_N1 < dc) && (d_lzb_N2 < dc))
                     status = 1;
      }
   }
   return status;

}

int test_end_delete_closed_molecules(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_2 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_3 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_4 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_5 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_6 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_7 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_8 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));

   int n_molecules_1 = mc.get_number_of_molecules();

   mc.close_molecule(imol_8);
   mc.close_molecule(imol_2);
   mc.close_molecule(imol_6);
   mc.close_molecule(imol_7);

   mc.display_molecule_names_table();
   mc.end_delete_closed_molecules();
   mc.display_molecule_names_table();

   int n_molecules_2 = mc.get_number_of_molecules();

   if ((n_molecules_1 - n_molecules_2) == 3)
      status = 1;

   mc.close_molecule(imol_1);
   mc.close_molecule(imol_3);
   mc.close_molecule(imol_4);
   mc.close_molecule(imol_5);

   return status;
}

int test_texture_as_floats(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (!mc.is_valid_map_molecule(imol_map)) std::cout << "Failed to read moorhen-tutorial-map-number-1.mtz" << std::endl;

   texture_as_floats_t tf = mc.get_map_section_texture(imol_map, 6, 0, -0.1, 0.2);

   std::cout << " image data size " << tf.image_data.size() << std::endl;

   if (tf.x_size > 10.0) {
      if (tf.y_size > 10.0) {
         if (tf.image_data.size() > 10000) {
            for (unsigned int i=0; i<15000; i+=1000) {
               // std::cout << "   " << i << " " << tf.image_data[i] << std::endl;
               if (tf.image_data[i] > 0.5)
                  status = 1;
            }
            for (unsigned int i=0; i<15000; i+=1000) {
               // std::cout << "   " << i << " " << tf.image_data[i] << std::endl;
               if (tf.image_data[i] > 1.0)
                  status = 0;
            }
         }
      }
   }

   return status;
}

int test_n_map_sections(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map = mc.read_ccp4_map(reference_data("emd_16890.map"), false);

   int n = mc.get_number_of_map_sections(imol_map, 2);

   if (n == 380) status = 1;

   return status;
}

#if 0
int test_rdkit_mol(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   mc.import_cif_dictionary(reference_data("ATP.cif"), coot::protein_geometry::IMOL_ENC_ANY);
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   RDKit::RWMol m = mc.get_rdkit_mol("ATP", coot::protein_geometry::IMOL_ENC_ANY);
   std::string smiles = RDKit::MolToSmiles(m);
   std::cout << smiles << std::endl;
   if (smiles == "[H]O[C@@]1([H])[C@@]([H])(O[H])[C@]([H])(n2c([H])nc3c(N([H])[H])nc([H])nc32)O[C@]1([H])C([H])([H])OP(~O)(~O)OP(~O)(~O)O[P+](~O)(~O)~O") status = 1;
   if (smiles == "[H]O[C@]1([H])[C@@]([H])(O[H])[C@@]([H])(C([H])([H])OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O[C@]1([H])n1c([H])nc2c(N([H])[H])nc([H])nc21") status = 1;
#endif
   return status;
}
#endif



int test_lsq_superpose(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_2 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   mc.clear_lsq_matches();
   mc.add_lsq_superpose_match("A", 185, 195, "A", 105, 115, 1);
   bool summary_to_screen = false;
   auto tm = mc.get_lsq_matrix(imol_1, imol_2, summary_to_screen);
   mc.lsq_superpose(imol_1, imol_2);
   for(unsigned int i=0; i<tm.rotation_matrix.size(); i++)
      std::cout << " " << tm.rotation_matrix[i];
   std::cout << std::endl;
   for(unsigned int i=0; i<tm.translation.size(); i++)
      std::cout << " " << tm.translation[i];
   std::cout << std::endl;
   mc.write_coordinates(imol_2, "superposed.pdb");
   mmdb::Atom *at_1 = mc.get_atom(imol_1, coot::atom_spec_t("A", 190, "", " CA ", ""));
   mmdb::Atom *at_2 = mc.get_atom(imol_2, coot::atom_spec_t("A", 110, "", " CA ", ""));
   clipper::Coord_orth co_1 = coot::co(at_1);
   clipper::Coord_orth co_2 = coot::co(at_2);
   double d2 = (co_2-co_1).lengthsq();
   double d = std::sqrt(d2);
   std::cout << "d " << d << std::endl;
   if (d < 0.32)
      status = 1;
   return status;
}

int test_alpha_in_colour_holder(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   std::string col = "#aabbccdd";
   coot::colour_holder ch(col);
   std::cout << "alpha: " << ch.alpha << std::endl;
   int status = 0;
   if (ch.alpha > 0.7 && ch.alpha < 0.9) {

      std::vector<std::pair<std::string, unsigned int> > indexed_residues_cids;
      std::map<unsigned int, std::array<float, 4> > colour_map;
      colour_map[0] = std::array<float, 4> {0.42222222, 0.7, 0.4, 0.5};
      colour_map[2] = std::array<float, 4> {0.42222222, 0.4, 0.7, 0.5};
      colour_map[1] = std::array<float, 4> {0.7, 0.4, 0.42222222, 0.5};
      indexed_residues_cids.push_back(std::make_pair("//A",2));
      indexed_residues_cids.push_back(std::make_pair("//A/100-200",1));
      indexed_residues_cids.push_back(std::make_pair("//A/130-150",0));
      int imol_1 = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
      mc.set_user_defined_bond_colours(imol_1, colour_map);
      mc.set_user_defined_atom_colour_by_selection(imol_1, indexed_residues_cids, true);
      std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";
      auto mesh = mc.get_bonds_mesh_instanced(imol_1, mode, true,  0.2, 1.0, false, false, true, 1);
      std::vector<std::pair<glm::vec4, unsigned int> > colour_count = colour_analysis(mesh);
      unsigned int n_transparent = 0;
      for(const auto &cc : colour_count) {
         if (cc.first[3] < 0.6)
            n_transparent += cc.second;
      }
      std::cout << "................. n_transparent " << n_transparent << std::endl;
      if (n_transparent > 5000) status = 1;
   }
   return status;
}

int test_Q_Score(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      mc.set_use_gemmi(true);
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         coot::validation_information_t vi = mc.get_q_score(imol, imol_map);
         coot::stats::single s = vi.get_stats();
         double q_mean = s.mean();
         std::cout << "q_mean " << q_mean << std::endl;
         if (q_mean > 0.8)
            if (q_mean < 1.0)
               status = 1;
      }
   }
   mc.close_molecule(imol);
   mc.close_molecule(imol_map);
   return status;
}

int test_assign_sequence(molecules_container_t &mc) {

   std::cout << "------------------ test_assign_sequence() " << std::endl;

   starting_test(__FUNCTION__);
   int status = 0;

   mc.set_use_gemmi(true);
   int imol = mc.read_pdb(reference_data("pdb7vvl.ent"));
   if (mc.is_valid_model_molecule(imol)) {
      int imol_map = mc.read_ccp4_map(reference_data("emd_32143.map"), false);
      if (mc.is_valid_map_molecule(imol_map)) {
         // strip the side-chains fro the hormone
         int imol_hormone = mc.copy_fragment_using_cid(imol, "//P");
         mmdb::Manager *mol = mc.get_mol(imol_hormone);
         if (mol) {
            int imod = 1;
            mmdb::Model *model_p = mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int n_res = chain_p->GetNumberOfResidues();
                  std::string chain_id = chain_p->GetChainID();
                  if (chain_id == "P") {
                     int nres = 0;
                     mmdb::PResidue *residue_table = 0;
                     chain_p->GetResidueTable(residue_table, nres);
                     for (int ires=0; ires<nres; ires++) {
                        int rn = residue_table[ires]->GetSeqNum();
                        mc.delete_side_chain(imol_hormone, chain_id, rn, "");
                     }
                     // now add the sidechain
                     std::string seq = "SVSEIQLMHNLGKHLNSMERVEWLRKKLQDVHNF";
                     mc.associate_sequence(imol_hormone, "P", seq);
                     mc.assign_sequence(imol_hormone, imol_map);

                     // residue 4 should be a GLU. Let's test for the
                     // presence of an OE1
                     coot::residue_spec_t rs("P", 4, "");
                     mmdb::Residue *r_test = mc.get_residue(imol_hormone, rs);
                     if (r_test) {
                        int n_res_atoms = 0;
                        mmdb::Atom **res_atoms;
                        r_test->GetAtomTable(res_atoms, n_res_atoms);
                        for (int i=0; i<n_res_atoms; i++) {
                           mmdb::Atom *at = res_atoms[i];
                           std::string atom_name(at->GetAtomName());
                           if (atom_name == " OE1") {
                              status = 1;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   mc.write_coordinates(imol, "para-sans-side-chains.pdb");

   return status;
}

int test_dictionary_conformers(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::vector<int> new_mols = mc.get_dictionary_conformers("TYR", coot::protein_geometry::IMOL_ENC_ANY, true);

   if (new_mols.size() == 108) status = 1;

   return status;
}

int test_ligand_distortion(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.get_monomer("ATP");

   mc.get_ligand_validation_vs_dictionary(imol, "//A/1", true);

   // needs more work

   return status;
}

int test_for_long_bonds(molecules_container_t &mc, int imol) {

   int state = -1; // unset
   if (mc.is_valid_model_molecule(imol)) {
      auto instanced_mesh = mc.get_bonds_mesh_instanced(imol, "COLOUR-BY-CHAIN-AND-DICTIONARY", false, 0.1f, 1.0f, false, false, true, 1);
      const auto &geom_vec = instanced_mesh.geom;
      unsigned int geom_vec_size = geom_vec.size();
      for (unsigned int i = 0; i < geom_vec_size; i++) {
         const auto &geom = geom_vec.at(i);
         const auto &inst_data_B_vec = geom.instancing_data_B;
         unsigned inst_data_B_vec_size = inst_data_B_vec.size();
         unsigned int n_long = 0;
         for (unsigned int j = 0; j < inst_data_B_vec_size; j++) {
            float z = inst_data_B_vec[j].size.z;
            if (z > 1.8f) {
               std::cout << j << " " << z << std::endl;
               n_long++;
            }
         }
         if (n_long == 0) state = 0;
         if (n_long > 0)  state = 1;
      }
   }
   return state;
}

int test_import_LIG_dictionary(molecules_container_t &mc) {

   int status = 0;
   starting_test(__FUNCTION__);
   mc.import_cif_dictionary(reference_data("LIG.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   int imol_pdb = mc.read_pdb("7vvl.pdb");
   status = mc.import_cif_dictionary(reference_data("LIG.cif"), imol_pdb);
   return status;
}

int test_tricky_ligand_problem(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("tutorial-1-with-nitrobenzene.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   mc.set_imol_refinement_map(imol_map);
   mc.import_cif_dictionary(reference_data("YXG-as-LIG.cif"), imol);
   mc.import_cif_dictionary(reference_data("nitrobenzene.cif"), imol);
   std::cout << "------------- nitrobenzene.cif had been read --------------" << std::endl;
   mc.refine_residues_using_atom_cid(imol, "//A/301", "SPHERE", 1000);
   mc.write_coordinates(imol, "nitrobenzene-refined.pdb");

   int state = test_for_long_bonds(mc, imol);
   if (state == 0) status = 1;

   return status;
}

int test_dictionary_acedrg_atom_types(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   mc.import_cif_dictionary(reference_data("YXG-as-LIG.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   std::vector<std::pair<std::string, std::string> > v = mc.get_acedrg_atom_types("LIG", coot::protein_geometry::IMOL_ENC_ANY);

   if (v.size() > 10) {
      std::cout << "types vector size: " << v.size() << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
         if (v[i].first == "C21")
            if (v[i].second == "C[6a](C[6a]C[6a]N[5a])(C[6a]C[6a]N[6a])(C[6a]C[6a]H){1|Cl<1>,1|N<2>,2|H<1>,4|C<3>}")
               status = 1;
         // std::cout << v[i].first << " " << v[i].second << std::endl;
      }
   }
   return status;
}

int test_dictionary_acedrg_atom_types_for_ligand(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   mc.import_cif_dictionary(reference_data("YXG-as-LIG.cif"), coot::protein_geometry::IMOL_ENC_ANY);
   int imol = mc.get_monomer_from_dictionary("LIG", coot::protein_geometry::IMOL_ENC_ANY, false);
   coot::acedrg_types_for_residue_t types = mc.get_acedrg_atom_types_for_ligand(imol, "//A/1");

   std::cout << "------------- types (bond) ----------- " << std::endl;
   for (unsigned int i=0; i<types.bond_types.size(); i++) {
      const auto &bond_types = types.bond_types[i];
      std::cout << "bond_type "
                << bond_types.atom_id_1 << " " << bond_types.atom_type_1 << " "
                << bond_types.atom_id_2 << " " << bond_types.atom_type_2
                << " bond-length: " << bond_types.bond_length << std::endl;
   }
   if (types.bond_types.size() > 10) status = 1;
   return status;
}


int test_links_in_model_read_via_gemmi(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   mc.set_use_gemmi(true);
   int imol = mc.read_pdb(reference_data("6dgd.cif"));
   if (mc.is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = mc.get_mol(imol);
      for (int imod = 1; imod <= mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_links = model_p->GetNumberOfLinks();
            std::cout << "Found n_links: " << n_links << std::endl;
            for (int i_link = 0; i_link < n_links; i_link++) {
               mmdb::Link *link_p = model_p->GetLink(i_link);
               // std::cout << "Link " << i_link << " " << link_p << std::endl;
            }
            if (n_links > 4) status = 1;
         }
      }
   }

   return status;
}

int test_delete_two_add_one_using_gemmi(molecules_container_t &mc) {


   starting_test(__FUNCTION__);
   int status = 0;
   mc.set_use_gemmi(true);

   int coordMolNo = mc.read_pdb(reference_data("./5a3h.pdb"));
   int mapMolNo = mc.read_mtz(reference_data("./5a3h_sigmaa.mtz"), "FWT", "PHWT", "", false, false);
   mc.set_imol_refinement_map(mapMolNo);

   mc.delete_using_cid(coordMolNo, "A/100-101/*", "LITERAL");
   mc.write_coordinates(coordMolNo, "add_terminal_test_tmp_1.pdb");

   int result = mc.add_terminal_residue_directly_using_cid(coordMolNo, "A/99");

   mc.write_coordinates(coordMolNo, "add_terminal_test_tmp_2.pdb");

   std::vector<mmdb::Residue *> A100s;
   mmdb::Manager *mol = mc.get_mol(coordMolNo);
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            std::cout << "    " << chain_p->GetChainID() << " " << n_res << " residues " << std::endl;
            std::string chain_id = chain_p->GetChainID();
            if (chain_id == "A") {
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  int seq_num = residue_p->GetSeqNum();
                  if (seq_num == 100)
                     A100s.push_back(residue_p);
               }
            }
         }
      }
   }
   std::cout << " here with " << A100s.size() << " residues A 100s" << std::endl;
   if (A100s.size() == 1)
      status = 1;
   return status;
}

int test_merge_ligand_and_gemmi_parse_mmcif(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

#ifdef USE_GEMMI_REALLY

  auto read_structure_from_string = [] (const std::string &data, const std::string& path){
    char *c_data = (char *)data.c_str();
    size_t size = data.length();
    // return gemmi::read_structure_from_char_array(c_data,size,path);
    return nullptr;
  };

   auto coordMolNo_1 = mc.read_pdb(reference_data("5a3h.mmcif"));
   // expect(coordMolNo_1).toBe(0)
   if (! mc.is_valid_model_molecule(coordMolNo_1))
     return status;

   int import_status = mc.import_cif_dictionary(reference_data("LZA.cif"), -999999);
   // expect(result_import_dict).toBe(1)
   if (import_status == 0)
     return status;

   auto ligandMolNo = mc.get_monomer_and_position_at("LZA", -999999, 0, 0, 0);
   // expect(ligandMolNo).toBe(1);

   auto merge_info = mc.merge_molecules(coordMolNo_1, std::to_string(ligandMolNo));
   // expect(merge_info.second.size()).toBe(1)

   mc.write_coordinates(coordMolNo_1, "my-mol-with-ligand.cif");

   bool merge_chain_parts = true;
   gemmi::Structure st = gemmi::read_structure_file("my-mol-with-ligand.cif");
   if (merge_chain_parts)
      st.merge_chain_parts();

   // std::string mmcifString = mc.get_molecule_atoms(coordMolNo_1, "mmcif");
   // auto st = read_structure_from_string(mmcifString, "test-molecule");
   gemmi::setup_entities(st);
   gemmi::add_entity_types(st, true);

   auto model = st.first_model();
   auto chains = model.chains;
   auto chain = chains[2];
   gemmi::ResidueSpan ligands = chain.get_ligands();
   // expect(ligands.length()).toBe(1);

   std::cout << "we found " << model.chains.size() << " model chains" << std::endl;

   std::cout << "debug:: test_merge_ligand_and_gemmi_parse_mmcif() merge_info second length "
	     << merge_info.second.size() << std::endl;

   std::cout << "debug:: test_merge_ligand_and_gemmi_parse_mmcif() ligands length() "
	     << ligands.size() << std::endl;

   std::cout << "chain things " << chain.whole().size() << std::endl;
   std::cout << "chain 0 things " << chains[0].name << " " << chains[0].whole().size() << std::endl;
   std::cout << "chain 1 things " << chains[1].name << " " << chains[1].whole().size() << std::endl;
   std::cout << "chain 2 things " << chains[2].name << " " << chains[2].whole().size() << std::endl;
   // std::cout << "chain 3 things " << chains[3].name << " " << chains[3].whole().size() << std::endl;

   std::cout << "get_ligands(): size    " << chains[2].get_ligands().size()   << std::endl;
   std::cout << "get_ligands(): length  " << chains[2].get_ligands().length() << std::endl;
   std::cout << "get polymer " << chains[2].get_polymer().length() << std::endl;
   std::cout << "get waters  " << chains[2].get_waters().length()  << std::endl;

   gemmi::EntityType et = chains[2].whole().at(0).entity_type;
   std::cout << "get_ligands(): whole " << gemmi::entity_type_to_string(et) << std::endl;

   gemmi::ResidueSpan whole = chains[2].whole();
   std::cout << "whole span length: " << whole.length() << std::endl;

   if (model.chains.size() == 3)
      if (chains[2].get_ligands().size() == 1)
         status = 1;

#endif

   return status;

}

int test_dictionary_atom_name_match(molecules_container_t &mc) {

   int status = 0;
   starting_test(__FUNCTION__);

   std::string comp_id_1 = "TYR";
   std::string comp_id_2 = "PTR";
   int imol_1 = coot::protein_geometry::IMOL_ENC_ANY;
   int imol_2 = coot::protein_geometry::IMOL_ENC_ANY;

   // 20241031-PE these are relatively modern dictionaries (from CCP4-9)
   mc.import_cif_dictionary(reference_data("TYR.cif"), imol_1);
   mc.import_cif_dictionary(reference_data("PTR.cif"), imol_2);

   int imol_pty = mc.get_monomer(comp_id_2);  // I want the dictionary, not the molecule

   if (mc.is_valid_model_molecule(imol_pty)) {
      std::map<std::string, std::string> m =
         mc.dictionary_atom_name_map(comp_id_1, imol_1, comp_id_2, imol_2);

      std::cout << "got a name map of size " << m.size() << std::endl;
      for (const auto &name : m)
         std::cout << "    " << name.first << " " << name.second << std::endl;
      if (m.size() > 8) status = 1;
   }

   return status;
}

int test_average_position_functions(molecules_container_t &mc) {

   auto close_position = [] (const std::vector<double> &p,
                             const clipper::Coord_orth &r) {
      double t = 0.9;
      if (abs(p[0] - r.x()) < t) {
         if (abs(p[1] - r.y()) < t) {
            if (abs(p[2] - r.z()) < t) {
               return true;
            }
         }
      }
      return false;
   };

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::vector<double> ap_1 = mc.get_residue_average_position(imol, "//A/15");
   std::vector<double> ap_2 = mc.get_residue_sidechain_average_position(imol, "//A/15");
   std::vector<double> ap_3 = mc.get_residue_CA_position(imol, "//A/15");

   int n_pass = 0;
   if (close_position(ap_1, clipper::Coord_orth(16.4, 10.1, 57.5))) n_pass += 1;
   if (close_position(ap_2, clipper::Coord_orth(16.2,  8.9, 56.0))) n_pass += 2;
   if (close_position(ap_3, clipper::Coord_orth(16.4, 11.5, 58.7))) n_pass += 4;

   std::cout << "here with n_pass " << n_pass << std::endl;
   if (n_pass == 7) status = 1;
   return status;
}

int test_set_occupancy(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int n_zero = 0;
   mc.set_occupancy(imol, "//A/18-19", 0.0);

   mmdb::Manager *mol = mc.get_mol(imol);
   for (int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        if (at->occupancy < 0.00001) n_zero++;
                     }
                  }
               }
            }
         }
      }
   }
   if (n_zero == 19) status = 1;
   return status;
}

int test_missing_residues(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      mc.delete_residue(imol, "A", 10, "");
      mc.delete_residue(imol, "A", 11, "");
      mc.delete_residue(imol, "A", 12, "");
      mc.delete_residue(imol, "A", 15, "");

      auto vec = mc.get_missing_residue_ranges(imol);
      for (const auto &rr : vec) {
         std::cout << "   " << rr.chain_id << " " << rr.res_no_start << " " << rr.res_no_end
	  	   << std::endl;
       }
       if (vec.size() == 5) {
          if (vec[0].res_no_start ==  10 && vec[0].res_no_end ==  12)
          if (vec[1].res_no_start ==  15 && vec[1].res_no_end ==  15)
          if (vec[2].res_no_start ==  37 && vec[2].res_no_end ==  45)
          if (vec[3].res_no_start ==  73 && vec[3].res_no_end ==  75)
          if (vec[4].res_no_start == 147 && vec[4].res_no_end == 165)
 	     status = 1;
       }
   }
   return status;
}

int test_mutation_info(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      mc.delete_residue(imol, "A", 10, "");
      mc.delete_residue(imol, "A", 11, "");
      mc.delete_residue(imol, "A", 12, "");
      mc.delete_residue(imol, "A", 15, "");

      std::string t("MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL");

      mc.associate_sequence(imol, "A", t);
      auto mi = mc.get_mutation_info(imol);
      std::cout << "mutation-info: "
		<< mi.mutations.size()  << " mutations "
		<< mi.insertions.size() << " insertions "
		<< mi.deletions.size()  << " deletions "
		<< std::endl;

      for (const auto &m : mi.mutations)
	 std::cout << "   mutation " << m.first << " " << m.second << std::endl;
      for (const auto &m : mi.insertions) {
	 std::cout << "   insertions " << m.start_resno << std::endl << "  ";
	 for (const auto &t : m.types)
	    std::cout << " " << t;
	 std::cout << std::endl;
      }
      for (const auto &m : mi.deletions)
	 std::cout << "   deletions " << m << std::endl;

      if (mi.mutations.size() > 2)
	 if (mi.insertions.size() > 2)
	    status = 1;
   }
   return status;
}

int test_scale_map(molecules_container_t &mc) {

   auto close_float = [] (float a, float b) {
      return fabsf(a - b) < 0.001;
   };

   starting_test(__FUNCTION__);
   int status = 0;

   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"),
                              "FWT", "PHWT", "W", false, false);
   if (mc.is_valid_map_molecule(imol_map)) {
      float sf = 2.4;
      float r_1 = mc.get_map_rmsd_approx(imol_map);
      mc.scale_map(imol_map, sf);
      float r_2 = mc.get_map_rmsd_approx(imol_map);
      float f = r_2 / r_1;
      if (close_float(f, sf))
         status = 1;
   } else {
      std::cout << "ERROR:: failed to read moorhen-tutorial-map-number-1.mtz" << std::endl;
   }
   return status;
}

int test_add_RNA_residue(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_coordinates(reference_data("5bjo-needs-E6.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      mc.add_terminal_residue_directly_using_cid(imol, "//E/5");
      mmdb::Residue *residue_p = mc.get_residue_using_cid(imol, "//E/6");
      if (residue_p) {
         status = true;
      }
   } else {
      std::cout << "failed to read 5bjo-needs-E6.pdb" << std::endl;
   }
   return status;
}

int test_HOLE(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_coordinates(reference_data("pdb2y5m.ent"));
   // get restraints for 15P, DVA, FVA, DLE and, of course ETA.
   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   std::vector<std::string> new_types = {"15P", "DVA", "FVA", "DLE", "ETA"};
   for (const auto &type : new_types) {
      int imol_type = mc.get_monomer_from_dictionary(type, imol_enc, true);
      if (! mc.is_valid_model_molecule(imol_type))
         std::cout << "Failed to get_monomer_from_dictionary() for " << type << std::endl;
   }
   if (mc.is_valid_model_molecule(imol)) {
      clipper::Coord_orth s(1.6, 14.9, -14.2);
      clipper::Coord_orth e(8.4, -16.8, -17.0);
      coot::instanced_mesh_t im = mc.get_HOLE(imol, s.x(), s.y(), s.z(), e.x(), e.y(), e.z());
      if (! im.geom.empty()) {
         const coot::instanced_geometry_t &ig = im.geom[0];
         unsigned int n_spots = ig.instancing_data_A.size();
         std::cout << "n_spots " << n_spots << std::endl;
         if (n_spots > 1000) {
            status = 1;
         }
      }
   } else {
      std::cout << "Failed to read pdb2y5m.ent" << std::endl;
   }
   return status;
}

int test_is_nucleic_acid(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   bool fail = false;
   int imol_1 = mc.read_coordinates(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_2 = mc.read_coordinates(reference_data("5bjo-needs-E6.pdb"));
   for (unsigned int res_no=0; res_no<=298; res_no++) {
      std::string cid = "//A/" + std::to_string(res_no);
      bool r = mc.residue_is_nucleic_acid(imol_1, cid);
      std::string rn = "--unset--";
      mmdb::Residue *residue_p = mc.get_residue_using_cid(imol_1, cid);
      if (residue_p)
         rn = residue_p->GetResName();
      if (r) {
         std::cout << "fail for imol_1 " << cid << " type " << rn << std::endl;
         fail = true;
         break;
      }
   }

   for (unsigned int res_no=1; res_no<=16; res_no++) {
      if (res_no == 6) continue;
      if (res_no == 7) continue;
      std::string cid = "//E/" + std::to_string(res_no);
      bool r = mc.residue_is_nucleic_acid(imol_2, cid);
      mmdb::Residue *residue_p = mc.get_residue_using_cid(imol_2, cid);
      std::string rn = "--unset--";
      if (residue_p)
         rn = residue_p->GetResName();
      if (! r) {
         std::cout << "fail for imol_2 " << cid << " type " << rn << std::endl;
         fail = true;
      } else {
         std::cout << "pass for imol_2 " << cid << " type " << rn << std::endl;
      }
   }

   if (fail)
      status = 0;
   else
      status = 1;
   return status;
}

int test_delete_all_carbohydrate(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("pdb8ox7.ent"));
   if (mc.is_valid_model_molecule(imol)) {
      mc.delete_all_carbohydrate(imol);
      mmdb::Manager *mol = mc.get_mol(imol);
      if (mol) {
         // check for remaining NAGs
         int n_nags = 0;
         for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int n_res = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<n_res; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        std::string rn = residue_p->GetResName();
                        if (rn == "NAG") n_nags++;
                     }
                  }
               }
            }
         }
         if (n_nags == 0) status = 1;
      }
   }
   return status;
}

int test_map_vertices_histogram(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   auto make_n_stars = [] (int counts) {
      std::string n_stars = "";
      for (int i=0; i<counts; i+=2000)
	 n_stars += "*";
      return n_stars;
   };

   int imol_map_1 = mc.read_ccp4_map(reference_data("emd_16890.map"), false);
   int imol_map_2 = mc.read_ccp4_map(reference_data("scale_res_emd_16890.mrc"), false);
   if (mc.is_valid_map_molecule(imol_map_1)) {
      if (mc.is_valid_map_molecule(imol_map_2)) {
	 unsigned int n_bins = 40;
	 coot::molecule_t::histogram_info_t histo =
	    mc.get_map_vertices_histogram(imol_map_1, imol_map_2,
					  160, 160, 160,
					  100, 0.16, n_bins);
	 unsigned int n_bins_hist = histo.counts.size();
	 std::cout << "n_bins_hist " << n_bins_hist << std::endl;
	 for (unsigned int i=0; i<n_bins_hist; i++) {
	    int counts = histo.counts[i];
	    float f_min = histo.base + static_cast<float>(i)   * histo.bin_width;
	    float f_max = histo.base + static_cast<float>(i+1) * histo.bin_width;
	    std::string n_stars = make_n_stars(counts);
	    std::cout << "   " << std::setw(9) << f_min << " " << std::setw(9) << f_max
		      << " " << std::setw(6) << counts << " " << n_stars << std::endl;
	    if (n_bins_hist >= 40)
	       if (counts > 100000)
		  status = 1;
	 }
      }
   }
   return status;
}

int test_non_XYZ_EM_map_status(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   // the test here is that we shouldd be seeing zero vertices outside the unit cell box
   // when the map is an EM map.

   int imol = mc.read_ccp4_map("initial_map_clement.ccp4", 0);
   if (mc.is_valid_map_molecule(imol)) {
      bool em_status = mc.is_EM_map(imol);
      if (em_status) {
         float radius = 20.0;
         float contour_level = 0.1;
         coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol, -111, 111, 111, radius, contour_level);
         std::cout << "n-vertices: " << map_mesh.vertices.size() << std::endl;
         std::cout << "n-triangles: " << map_mesh.triangles.size() << std::endl;
         if (map_mesh.triangles.size() == 0)
            status = 1;
      } else {
         std::cout << "ERROR:: in " << __FUNCTION__ << " the EM map is marked as a non-EM map" << std::endl;
      }
   } else {
      std::cout << "Failed to read initial_map_clement.ccp4 " << std::endl;
   }

   return status;
}

int test_radius_of_gyration(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {
      double rg = mc.get_radius_of_gyration(imol);
      std::cout << "Radius of gyration: " << rg << std::endl;
      // Expect a positive value for a valid molecule
      if (rg > 0.0 && rg < 100.0)
         status = 1;
      else
         std::cout << "Unexpected radius of gyration value: " << rg << std::endl;
   } else {
      std::cout << "Invalid model molecule for radius of gyration test." << std::endl;
   }
   mc.close_molecule(imol);
   return status;
}

int test_temperature_factor_of_atom(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   float b1 = mc.get_temperature_factor_of_atom(imol, "//A/8/CB");
   float b2 = mc.get_temperature_factor_of_atom(imol, "//x/8/CB");

   std::cout << "debug b1 " << b1 << " b2 " << b2 << std::endl;

   if (b2 < 0.0)
      if (b1 < 40.0)
         if (b1 > 39.0)
            status = 1;

   return status;
}

int test_water_spherical_variance(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   if (mc.is_valid_model_molecule(imol)) {
      unsigned int n_waters = mc.add_waters(imol, imol_map);
      std::cout << "DEBUG:: test_water_spherical_variance(): n_waters: " << n_waters << std::endl;
      mc.write_coordinates(imol, "with-waters.pdb");
      bool all_pass = true;
      for (unsigned int rn=1; rn<=50; rn++) {
         std::string atom_cid = "//B/" + std::to_string(rn) + "/O";
         std::pair<float,float> mv = mc.get_mean_and_variance_of_density_for_non_water_atoms(imol, imol_map);
         if (mv.first > 0) { // this test that the function worked as expected/hoped
            float sv = mc.get_spherical_variance(imol_map, imol, atom_cid, mv.first);
            // std::cout << "spherical variance: " << atom_cid << " " << sv << std::endl;
            if (sv <= 0.0) all_pass = false;
         }
      }
      if (all_pass) status = 1;
   }
   return status;
}


int test_dedust(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol_map = mc.read_ccp4_map(reference_data("emd_16890.map"), false);
   int imol_new = mc.dedust_map(imol_map);
   mc.write_map(imol_new, "dedust-16890.map");
   return status;
}

int test_atom_overlaps(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   std::vector<coot::plain_atom_overlap_t> aov = mc.get_atom_overlaps(imol);
   for (unsigned int i=0; i<aov.size(); i++) {
      if (i > 10) continue;
      const auto &ao = aov[i];
      std::cout << "Overlapping atom " << ao.atom_spec_1 << " " << ao.atom_spec_2
                << " with overlap volume " << ao.overlap_volume << std::endl;
      if (ao.overlap_volume > 2.0) status = 1;
   }

   return status;
}

#include "coot-utils/json.hpp"
using json = nlohmann::json;

int test_pucker_info(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("2pwt.cif"));
   std::string pucker_info_json = mc.get_pucker_analysis_info(imol);
   // parsign an empty json string causes a crash
   if (pucker_info_json.size() > 10) {
      json j = json::parse(pucker_info_json);
      unsigned int count = 0;
      for (json::iterator it=j.begin(); it!=j.end(); ++it) {
         count += 1;
         if (count > 5) continue;
         json &item = *it;
         std::string s = item.dump(4);
         std::cout << s << std::endl;
      }
      if (count > 10) status = 1;
   }
   return status;
}

int test_inner_bond_kekulization(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   int imol = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   if (mc.is_valid_model_molecule(imol)) {

      // first delete the water
      // status,atom_count
      std::pair<int, unsigned int> dw = mc.delete_residue(imol, "B", 1, "");
      if (dw.first == 1) {

         std::string mode("COLOUR-BY-CHAIN-AND-DICTIONARY");
         std::string selection_cid_1 = "//*";
         std::string selection_cid_2 = "//*/(!HOH)";
         coot::instanced_mesh_t m_1 = mc.get_bonds_mesh_for_selection_instanced(imol, selection_cid_1, mode,
                                                                                false, 0.1, 1.0, false, false, false, 2);
         coot::instanced_mesh_t m_2 = mc.get_bonds_mesh_for_selection_instanced(imol, selection_cid_2, mode,
                                                                                false, 0.1, 1.0, false, false, false, 2);
         std::cout << "--------------------------------- mesh 1 ---------------------------" << std::endl;
         std::vector<std::pair<glm::vec4, unsigned int> > r_1 = colour_analysis(m_1);
         std::cout << "--------------------------------- mesh 2 ---------------------------" << std::endl;
         std::vector<std::pair<glm::vec4, unsigned int> > r_2 = colour_analysis(m_2);
         status = 1; // failure on mismatch
         for (unsigned int i=0; i<r_1.size(); i++) {
            const auto &cp_1 = r_1[i];
            const auto &cp_2 = r_2[i];
            if (cp_1.second != cp_2.second) status = 0;
         }
      }
   }
   return status;
}

int test_gaussian_surface_to_map_molecule(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   //!
   //! @param imol is the model molecule index
   //! @param cid is the atom selection CID
   //! @param sigma default 4.4
   //! @param contour_level default 4.0
   //! @param box_radius default 5.0
   //! @param grid_scale default 0.7
   //! @param b_factor default 100.0 (use 0.0 for no FFT-B-factor smoothing)
   //!
   //! @return a new molecule index for the map or -1 on failur
   float sigma = 4.0;
   float box_radius = 5.0;
   float grid_scale = 1.0;
   float fft_b_factor = 30.0;
   std::string cid  = "//A";
   int imol_new = mc.gaussian_surface_to_map_molecule(imol, cid, sigma, box_radius, grid_scale, fft_b_factor);
   if (mc.is_valid_map_molecule(imol_new)) {
      status = 1;
      if (true) {
         mc.write_map(imol_new, "gaussian-surface-as-map.map");
         int imol_mask = mc.make_mask(imol_map, imol, "//A", 13.0f);
         mc.write_map(imol_mask, "mask-map.map");
      }
   }
   return status;
}

int test_template(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      mc.set_use_gemmi(true);
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

int n_tests = 0;
static std::vector<std::pair<std::string, int> > test_results;

void
write_test_name(const std::string &test_name) {

   std::ofstream f(".current-test");
   f << "\"" << test_name << "\"" << "\n";
   f.close();
}

int
run_test(int (*test_func) (molecules_container_t &mc), const std::string &test_name, molecules_container_t &mc) {

   n_tests++;
   write_test_name(test_name);
   int status = test_func(mc);
   std::string status_string = "FAIL: ";
   std::string uncol = "[m";
   std::string col = "[31m";
   if (status == 1) {
      status_string = "PASS: ";
      col = "[32m";
   }
   std::cout << status_string << std::setw(40) << std::left << test_name << col << " ⬤ " << uncol << std::endl;
   test_results.push_back(std::make_pair(test_name, status));

   return status;
}

void
print_results_summary() {

   std::cout << "n_tests: " << n_tests << std::endl;
   unsigned int n_failed = 0;
   // for (const auto &result[function_name, status] : test_results) { // structured binding
   std::cout << "LIGHTS: ";
   unsigned int count = 0;
   for (const auto &result : test_results) {
      count++;
      const auto &status = result.second;
      if (status == 0) {
         n_failed++;
         std::cout << "[31m⬤ ";
      } else {
         std::cout << "[32m⬤ ";
      }
      if (count%40 == 0) {
         std::cout << "\n,       ";
      }
   }
   std::cout << "[m  failures: " << n_failed << "/" << n_tests << std::endl;

   if (n_failed > 0) {
      std::cout << "Test summary: " << n_failed << " failed tests of " << n_tests << std::endl;
      for (const auto &result : test_results) {
         const auto &name   = result.first;
         const auto &status = result.second;
         if (status == 0) {
            std::cout << "FAIL:   " << name << std::endl;
         }
      }
   } else {
      std::cout << "Test summary: all " << n_tests << " tests passed "  << std::endl;
   }

}

int main(int argc, char **argv) {

   int status = 0;
   write_test_name("---");

   bool last_test_only = false;
   if (argc > 1) {
      std::string arg(argv[1]);
      if (arg == "last-test-only")
         last_test_only = true;
   }

   int all_tests_status = 1; // fail!

   {
      // we use this scope so that molecule_container_t mc goes out of scope before
      // the end of the program. Hopefully this can allow us to distinguish
      // between memory consumed by the mc, and memory that we can no longer
      // access (memory allocated but not deleted - e.g. Atom Selections)
      //
      molecules_container_t mc(true); // quiet
      mc.set_use_gemmi(true);

      // now check that the monomer library was read OK
      int imol = mc.get_monomer("ATP");
      if (! mc.is_valid_model_molecule(imol)) {
         std::cout << "Failed to read the monomer library - exit now" << std::endl;
         exit(1);
      } else {
         mc.close_molecule(imol);
      }

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
         status += run_test(test_rsr_using_multi_atom_cid, "multi-atom-cid RSR",    mc);
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
         status += run_test(test_fill_partial,          "Fill partially-filled residues", mc);
         status += run_test(test_delete_atom,           "delete atom",              mc);
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
         status += run_test(test_ligand_contact_dots,   "ligand contact dots",      mc);
         status += run_test(test_difference_map_peaks,  "Difference Map Peaks",     mc);
         status += run_test(test_rama_validation,       "rama validation 2",        mc); // for the plot, not the graph
         status += run_test(test_ramachandran_analysis, "ramachandran analysis",    mc); // for the graph, not the plot
         status += run_test(test_non_standard_residues, "non-standard residues",    mc);
         status += run_test(test_import_cif_dictionary, "import cif dictionary",    mc);
         status += run_test(test_add_terminal_residue,  "add terminal residue",     mc);
         status += run_test(test_sequence_generator,    "Make a sequence string",   mc);
         status += run_test(test_instanced_rota_markup, "Instanced rotamer mesh",   mc);
         status += run_test(test_new_position_for_atoms,"New positions for atoms",  mc);
         status += run_test(test_molecular_representation, "Molecular representation mesh", mc);
         // remove these for now - I know why they don't work and they are slow.
         // status += run_test(test_rigid_body_fit,        "Rigid-body fit", mc);
         // status += run_test(test_ligand_fitting_here,   "Ligand fitting here",      mc);
         // status += run_test(test_jiggle_fit,            "Jiggle-fit",               mc);
         // status += run_test(test_jiggle_fit_with_blur,  "Jiggle-fit-with-blur",     mc);
         // status += run_test(test_ligand_fitting_in_map, "ligand fitting in map",    mc);
         status += run_test(test_multiligands_lig_bonding, "Some multiligands bonding", mc);
         status += run_test(test_gltf_export_via_api,   "glTF via api", mc);
         status += run_test(test_long_name_ligand_cif_merge, "Long-name ligand cif merge", mc);
         status += run_test(test_user_defined_bond_colours_v3, "user-defined colours v3", mc);
         status += run_test(test_gltf_export,           "glTF export", mc);
         status += run_test(test_5char_ligand_merge,    "5-char ligand merge", mc);
         status += run_test(test_thread_pool,           "thread pool",    mc);
         // status += run_test(test_thread_launching,      "thread launching",    mc); // this is not a helpful test
         status += run_test(test_cif_gphl_chem_comp_info, "extracting gphl info",    mc);
         // status += run_test(test_test_the_threading,    "threading speed test",    mc); // not helpful
         // status += run_test(test_contouring_timing,     "contouring timing",    mc); // not helpful

         //reinstate this test when mmdb chain selection works
         // status += run_test(test_mmcif_atom_selection,  "mmCIF atom selection",    mc);

         //reinstate this test when gemmi is used for writing cif files
         // status += run_test(test_mmcif_as_string,       "mmCIF as string",    mc);
         status += run_test(test_pdb_as_string,         "PDB as string",    mc);
         status += run_test(test_cif_writer,            "mmCIF dictionary writer",    mc);
         status += run_test(test_residues_near_residues, "residues near residues",    mc);
         status += run_test(test_electro_molecular_representation, "electro molecular representation mesh", mc);
         status += run_test(test_replace_fragment,      "replace fragment",         mc);
         status += run_test(test_ncs_chains,            "NCS chains",         mc);
         status += run_test(test_omega_5tig_cif,        "Omega for 5tig cif",         mc);
         // status += run_test(test_jiggle_fit_params,     "actually testing for goodness pr params", mc); // not useful
         status += run_test(test_dark_mode_colours,     "light vs dark mode colours", mc);
         status += run_test(test_read_extra_restraints, "read extra restraints", mc);
         status += run_test(test_map_histogram,         "map histogram", mc);
         status += run_test(test_auto_read_mtz,         "auto-read-mtz", mc);
         status += run_test(test_read_a_missing_map,    "read a missing map file ", mc);
         status += run_test(test_colour_map_by_other_map, "colour-map-by-other-map", mc);
         status += run_test(test_something_filo,        "Self something filo", mc);
         status += run_test(test_self_restraints,       "Self restraints mesh", mc);
         status += run_test(test_other_user_defined_colours_other, "New colour test", mc);
         status += run_test(test_is_em_map,             "EM map flag is correctly set?", mc);
         status += run_test(test_user_defined_bond_colours_v2, "user-defined bond colours v2", mc);
         // reinstate when add alt conf has been added
         // status += run_test(test_alt_conf_and_rotamer,            "Alt Conf then rotamer", mc);
         status += run_test(test_editing_session_tutorial_1, "an Tutorial 1 editing session",         mc);
         status += run_test(test_broken_function, "Something was broken",         mc);
         status += run_test(test_delete_side_chain, "delete side chain", mc);
         status += run_test(test_colour_rules, "colour rules", mc);
         status += run_test(test_mmrrcc, "MMRRCC", mc);
         status += run_test(test_instanced_bonds_mesh, "insta bonds mesh", mc);
         status += run_test(test_instanced_bonds_mesh_v2, "test instanced bond selection v2", mc);
         status += run_test(test_utils, "utils", mc);
         status += run_test(test_svg, "svg string", mc);
         status += run_test(test_superpose, "SSM superpose ", mc);
         status += run_test(test_multi_colour_rules, "multi colour rules ", mc);
         status += run_test(test_non_drawn_atoms, "non-drawn atoms", mc);
         status += run_test(test_symmetry, "symmetry", mc);
         status += run_test(test_add_hydrogen_atoms, "add hydrogen atoms", mc);
         status += run_test(test_set_rotamer, "set rotamer ", mc);
         status += run_test(test_alt_conf_and_rotamer_v2, "alt-conf and rotamer v2 ", mc);
         status += run_test(test_moorhen_h_bonds, "moorhen H-bonds ", mc);
         status += run_test(test_number_of_hydrogen_atoms, "number of hydrogen atoms ", mc);
         status += run_test(test_cell, "cell", mc);
         status += run_test(test_map_centre, "map centre", mc);
         status += run_test(test_dragged_atom_refinement, "dragged atom refinement", mc);
         status += run_test(test_bespoke_carbon_colour, "bespoke carbon colours ", mc);
         status += run_test(test_replace_model_from_file, "replace model from file", mc);
         status += run_test(test_user_defined_bond_colours, "user-defined bond colours", mc);
         status += run_test(test_replace_map, "replace map from mtz", mc);
         status += run_test(test_residue_name_group, "residue name group", mc);
         status += run_test(test_mask_atom_selection, "mask atom selection", mc);
         status += run_test(test_write_map_is_sane, "write map is sane",    mc);
         status += run_test(test_replace_large_fragment,      "refine and replace large fragment",         mc);
         status += run_test(test_molecule_diameter, "molecule diameter",    mc);
         status += run_test(test_B_factor_multiply, "B-factor multiply",    mc);
         status += run_test(test_change_chain_id, "change chain id",    mc);
         status += run_test(test_17257, "read emd_17257.map.gz",    mc);
         status += run_test(test_get_diff_map_peaks, "get diff map peaks",    mc);
         status += run_test(test_shiftfield_b_factor_refinement, "Shiftfield B",    mc);
         status += run_test(test_non_drawn_CA_bonds,       "non-drawn bonds in CA+LIGANDS", mc);
         status += run_test(test_change_chain_id_1,        "change chain-id filo-1", mc);
         status += run_test(test_split_model,              "Split model", mc);
         status += run_test(test_make_ensemble,            "Make Ensemble", mc);
         status += run_test(test_end_delete_closed_molecules, "end delete close molecules", mc);
         status += run_test(test_moorhen_h_bonds, "moorhen H-bonds ", mc);
         status += run_test(test_texture_as_floats, "Texture as Floats ", mc);
         status += run_test(test_n_map_sections, "N map sections ", mc);
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
         status += run_test(test_pdbe_dictionary_depiction, "pdbe dictionary depiction", mc);
         // status += run_test(test_rdkit_mol, "RDKit mol", mc);
#endif

#ifdef USE_GEMMI
         status += run_test(test_disappearing_ligand,   "Disappearing ligand", mc);
#endif
         // Note to self:
         // change the autofit_rotamer test so that it tests the change of positions of the atoms of the neighboring residues.
      }

      {
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#endif
         // status += run_test(test_lsq_superpose, "LSQ superpose", mc);
         // status += run_test(test_change_rotamer, "Change Rotamer (Filo)", mc);
         // status += run_test(test_alpha_in_colour_holder, "Alpha value in colour holder", mc);
         // status += run_test(test_gaussian_surface, "Gaussian surface", mc);
         // status += run_test(test_Q_Score, "Q Score", mc);
         // status += run_test(test_assign_sequence, "Assign Sequence", mc);
         // status += run_test(test_undo_and_redo_2, "Undo and redo 2", mc);
         // status += run_test(test_gltf_export_via_api,   "glTF via api", mc);
         // status += run_test(test_import_ligands_with_same_name_and_animated_refinement, "Test import ligands with same name and animated refinement", mc);
         // status += run_test(test_dictionary_conformers,   "Dictionary Conformers", mc);
         // status += run_test(test_ligand_distortion,   "Ligand Distortion", mc);
         // status += run_test(test_import_LIG_dictionary,   "Import LIG.cif", mc);
         // status += run_test(test_tricky_ligand_problem,   "Tricky Ligand import/refine", mc);
         // status += run_test(test_dictionary_acedrg_atom_types, "Acedrg atom types", mc);
         // status += run_test(test_dictionary_acedrg_atom_types_for_ligand, "Acedrg atom types for ligand", mc);
         // status += run_test(test_long_name_ligand_cif_merge, "test long name ligand cif merge", mc);
         // status += run_test(test_merge_ligand_and_gemmi_parse_mmcif, "test_merge_ligand_and_gemmi_parse_mmcif", mc);
         // status += run_test(test_delete_two_add_one_using_gemmi, "test_delete_two_add_one_using_gemmi", mc);
         // status += run_test(test_dictionary_atom_name_match, "dictionary atom names match", mc);
         // status += run_test(test_average_position_functions, "average position functions", mc);

         // status += run_test(test_get_torsion, "get_torsion", mc);
         // status += run_test(test_set_occupancy, "set occupancy", mc);
         // status += run_test(test_missing_residues, "missing residues", mc);
         // status += run_test(test_mutation_info, "mutation info", mc);
         // status += run_test(test_scale_map, "scale_map", mc);
         // status += run_test(test_add_RNA_residue, "add RNA residue", mc);
         // status += run_test(test_HOLE, "HOLE", mc);
         // status += run_test(test_is_nucleic_acid, "is nucleic acid?", mc);
         // status += run_test(test_delete_all_carbohydrate, "delete all carbohydrate", mc);
         // status += run_test(test_instanced_goodsell_style_mesh, "instanced goodsell style mesh", mc);
         // status += run_test(test_map_vertices_histogram, "map vertices histogram", mc);
         // status += run_test(test_non_XYZ_EM_map_status, "non-XYZ map status", mc);

         // put these up
         // status += run_test(test_radius_of_gyration, "radius of gyration", mc);
         // status += run_test(test_temperature_factor_of_atom, "temperature factor of atom", mc);
         // status += run_test(test_water_spherical_variance, "water spherical variance", mc);
         // status += run_test(test_dedust, "dedust", mc);  .... maybe not this one
         // status += run_test(test_atom_overlaps, "atom overlaps", mc);
         // status += run_test(test_pucker_info, "pucker info", mc);
         // status += run_test(test_set_residue_to_rotamer_number, "set residue", mc);
         // status += run_test(test_inner_bond_kekulization, "inner-bond kekulization", mc);
         status += run_test(test_gaussian_surface_to_map_molecule, "gaussian-surface to map", mc);
         if (status == n_tests) all_tests_status = 0;

         print_results_summary();
      }
   }

   std::this_thread::sleep_for(std::chrono::milliseconds(500));
   return all_tests_status;

}
