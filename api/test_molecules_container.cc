
#include "molecules_container.hh"

int main(int argc, char **argv) {

   int status = 0;

   molecules_container_t mc;
   int imol = mc.read_pdb("test.pdb");
   coot::simple_mesh_t rvmm = mc.ramachandran_validation_markup_mesh(imol);

   std::cout << "mesh: " << rvmm.vertices.size() << " vertices and " << rvmm.triangles.size()
             << " traingles" << std::endl;

   return status;

}
