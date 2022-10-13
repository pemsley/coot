
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include "molecules_container.hh"

void starting_test(const char *func) {

   std::cout << "Starting test " << func << "()" << std::endl;

}

int test_auto_fit_rotamer() {

   starting_test(__FUNCTION__);

   int status = 0; // initially fail status
   molecules_container_t mc;
   mc.geometry_init_standard();
   int imol = mc.read_pdb("tm-A.pdb");
   int imol_map = mc.read_mtz("rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", false, false);

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

   int imol_map = mc.read_mtz("rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", false, false);

   clipper::Coord_orth p(55, 10, 10);
   float radius = 12;
   float contour_level = 0.13;
   coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
   std::cout << "density mesh: " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
             << " triangles" << std::endl;


   // --- auto-fit rotamer

   status = test_auto_fit_rotamer();
   status_string = "FAIL:";
   if (status == 1)
      status_string = "PASS:";
   std::cout << status_string << " auto_fit_rotamer() test status " << status << std::endl;


   // add a test for:
   // delete_atom
   // delete_atoms
   // delete_residue

   return status;

}
