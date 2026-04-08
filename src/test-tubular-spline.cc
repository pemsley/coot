
#include "coot-utils/atom-selection-container.hh"
#include "MoleculesToTriangles/CXXClasses/tubes.hh"

int main(int argc, char **argv) {

   int status = 0;
   std::string file_name = "pdb9eih.ent";
   bool use_gemmi = false;
   atom_selection_container_t asc = get_atom_selection(file_name, use_gemmi);
   if (asc.read_success) {
      std::string atom_selection_str = "//";
      std::string colour_scheme = "Chain";
      int secondaryStructureUsageFlag = 0;
      float radius_for_coil = 0.4;
      int Cn_for_coil = 3;
      int accuracy_for_coil = 12;
      unsigned int n_slices_for_coil = 12;
      /* -- removed until tubes can be compiled without coot libs/includes
      coot::simple_mesh_t mesh = make_tubes_representation(asc.mol,
                                                           atom_selection_str,
                                                           colour_scheme,
                                                           radius_for_coil,
                                                           Cn_for_coil, accuracy_for_coil,
                                                           n_slices_for_coil,
                                                           secondaryStructureUsageFlag);
                                                         */
      coot::simple_mesh_t mesh;
      float roughness = 0.2;
      float smoothnesss =- 0.8;
      bool use_binary_format = true;
      mesh.export_to_gltf("pdb9eih-spline.gltf", roughness, smoothnesss, use_binary_format);

   } else {
      status = 1;
   }
   return status;
}
