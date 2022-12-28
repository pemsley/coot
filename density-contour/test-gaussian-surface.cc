
#include <iostream>
#include <string>
#include "gaussian-surface.hh"

int main(int argc, char **argv) {

   if (argc > 1) {
      std::string file_name = argv[1];

      mmdb::Manager *mol = new mmdb::Manager;
      mol->ReadCoorFile(file_name.c_str());

      coot::gaussian_surface_t gauss_surf(mol, "A");
      coot::simple_mesh_t mesh = gauss_surf.get_surface();

      std::cout << "test-gaussian-surface got " << mesh.vertices.size() << " vertices and " << mesh.triangles.size() << std::endl;
   }
   return 0;
}
