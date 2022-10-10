
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include "molecules_container.hh"

int main(int argc, char **argv) {

   int status = 0;

   molecules_container_t mc;
   int imol = mc.read_pdb("pdb5a3h.ent");

   // --- rama mesh
   
   coot::simple_mesh_t rvmm = mc.ramachandran_validation_markup_mesh(imol);
   std::cout << "rama mesh: " << rvmm.vertices.size() << " vertices and " << rvmm.triangles.size()
             << " triangles" << std::endl;

   // Let's look at the colours of the balls.
   if (false)
      for (unsigned int i=0; i<rvmm.vertices.size(); i+=100)
         std::cout << i << " " << glm::to_string(rvmm.vertices[i].color) << std::endl;


   // --- density mesh

   int imol_map = mc.read_mtz("rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "W", false, false);

   clipper::Coord_orth p(55, 10, 10);
   float radius = 12;
   float contour_level = 0.13;
   coot::simple_mesh_t map_mesh = mc.get_map_contours_mesh(imol_map, p.x(), p.y(), p.z(), radius, contour_level);
   std::cout << "density mesh: " << map_mesh.vertices.size() << " vertices and " << map_mesh.triangles.size()
             << " triangles" << std::endl;

   return status;

}
