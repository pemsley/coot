
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

int test_auto_fit_rotamer() {

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
   std::vector<int> res_nos = {13, 18, 20, 49};
   for (const auto &rn : res_nos) {
      coot::residue_spec_t rs("A", rn, "");
      residues_for_flipping.push_back(rs);
   }

   int imol = mc.read_pdb(reference_data("gideondoesntapprove.pdb"));
   int imol_map = mc.read_mtz(reference_data("gideondoesntapprove.mtz"), "FWT", "PHWT", "W", false, false);

   float rmsd_diff_map_1 = mc.get_map_rmsd_approx(imol_map);

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
   if (n_flipped == res_nos.size())
      status = 1;

   // 20221016-PE now update the maps!

   return status;

}

int test_density_mesh(molecules_container_t &mc) {

   int status = 0;
   int imol_map = mc.read_mtz("rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", false, false);

   clipper::Coord_orth p(55, 10, 10);
   float radius = 12;
   float contour_level = 0.13;
   coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
   std::cout << "density mesh: " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
             << " triangles" << std::endl;

   if (map_mesh.vertices.size() > 30000)
      status = 1;

   return status;
}

int main(int argc, char **argv) {

   int status = 0;
   std::string status_string = "FAIL";

   molecules_container_t mc;
   mc.geometry_init_standard();
   mc.fill_rotamer_probability_tables();

   std::string coords_fn = "tm-A.pdb";
   int imol = mc.read_pdb(coords_fn);

   // --- rama mesh
   
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


   // --- density mesh

   status = test_density_mesh(mc);
   status_string = "FAIL:";
   if (status == 1)
      status_string = "PASS:";
   std::cout << status_string << " test_density_mesh() test status " << status << std::endl;

 
   // --- auto-fit rotamer

   status = test_auto_fit_rotamer();
   status_string = "FAIL:";
   if (status == 1)
      status_string = "PASS:";
   std::cout << status_string << " test_auto_fit_rotamer() test status " << status << std::endl;

   // --- pepflips

   status = test_pepflips(mc);
   status_string = "FAIL:";
   if (status == 1)
      status_string = "PASS:";
   std::cout << status_string << " test_pepflips() test status " << status << std::endl;


   // add a test for:
   // delete_atom
   // delete_atoms
   // delete_residue

   return status;

}
