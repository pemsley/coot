
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
         mc.flip_peptide(imol, res_spec, "");
         coot::Cartesian pt_2(at->x, at->y, at->z);
         double dd = coot::Cartesian::lengthsq(pt_1, pt_2);
         double d = std::sqrt(dd);
         if (d > 3.0) {
            n_flipped++;
         }
      }
   }

   // This flip A100 is OK to begin with so flipping it sould make the GruPoints
   // worse - and it does.

   if (n_flipped == res_nos.size()) {

      std::string residue_cid = "//A/100";
      auto rs = mc.residue_cid_to_residue_spec(imol, residue_cid);
      if (! rs.empty()) {
         const auto &res_spec = rs;
         coot::atom_spec_t atom_spec(res_spec.chain_id, res_spec.res_no, res_spec.ins_code, " O  ","");
         mmdb::Atom *at = mc.get_atom(imol, atom_spec);
         if (at) {
            coot::Cartesian pt_1(at->x, at->y, at->z);
            mc.flip_peptide_using_cid(imol, residue_cid, "");
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
   std::string residue_cid = "//A/14";
   auto rs = mc.residue_cid_to_residue_spec(imol, residue_cid);
   if (! rs.empty()) {
      const auto &res_spec = rs;
      coot::atom_spec_t atom_spec(res_spec.chain_id, res_spec.res_no, res_spec.ins_code, " O  ","");
      mmdb::Atom *at = mc.get_atom(imol, atom_spec);
      if (at) {
         mc.flip_peptide_using_cid(imol, residue_cid, "");
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
   std::string residue_cid = "//A/14";
   coot::atom_spec_t atom_spec("A", 14, "", " O  ", "");
   mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
   coot::Cartesian pt_1(at_1->x, at_1->y, at_1->z);
   mc.flip_peptide_using_cid(imol, residue_cid, "");
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
   std::cout << "rama mesh: " << rvmm.vertices.size() << " vertices and " << rvmm.triangles.size()
             << " triangles" << std::endl;

   // Let's look at the colours of the balls.
   if (false)
      for (unsigned int i=0; i<rvmm.vertices.size(); i+=100)
         std::cout << i << " " << glm::to_string(rvmm.vertices[i].color) << std::endl;

   // Rama dodecs
   coot::simple_mesh_t rota_mesh = mc.get_rotamer_dodecs(imol);
   std::cout << "rota mesh: " << rota_mesh.vertices.size() << " vertices and " << rota_mesh.triangles.size()
             << " triangles" << std::endl;

   if (rota_mesh.vertices.size() > 2000) status = 1;

   return status;

}

int test_density_mesh(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;
   // this could be any mtz file I suppose
   int imol_map = mc.read_mtz("rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", false, false);

   std::cout << "................ test_density_mesh imol_map is " << imol_map << std::endl;

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

   molecules_container_t mc;

#if 0
   mc.fill_rotamer_probability_tables();

   // --- rama mesh
   status += run_test(test_rama_balls_mesh, "rama balls mesh", mc);

   // --- density mesh
   status += run_test(test_density_mesh, "density mesh", mc);
 
   // --- auto-fit rotamer
   status += run_test(test_auto_fit_rotamer, "auto-fit rotamer", mc);

   // --- pepflips
   status += run_test(test_pepflips, "pepflips", mc);

   // --- updating maps
   status += run_test(test_updating_maps, "updating maps", mc);

   // --- undo
   status += run_test(test_undo_and_redo, "undo and redo", mc);

#endif

   status += run_test(test_delete_atom, "delete atom", mc);

   status += run_test(test_delete_residue, "delete residue", mc);

   // add a test for:
   // delete_atoms
   // rsr triple

   int all_tests_status = 1;
   if (status == n_tests) all_tests_status = 0;

   return all_tests_status;

}
