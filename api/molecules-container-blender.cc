

#include "molecules_container.hh"
#include "blender-mesh.hh"

std::vector<float>
molecules_container_t::get_colour_table_for_blender(int imol) {

   if (imol >=0)
      if (imol < static_cast<int>(molecules.size()))
         return molecules[imol].get_colour_table_for_blender();

   std::vector<float> v;
   return v;
}

std::vector<float>
molecules_container_t::get_vertices_for_blender(int imol) {

   if (imol >=0)
      if (imol < static_cast<int>(molecules.size()))
         return molecules[imol].get_vertices_for_blender();
   std::vector<float> v;
   return v;
}

std::vector<int>
molecules_container_t::get_triangles_for_blender(int imol) {

   if (imol >=0)
      if (imol < static_cast<int>(molecules.size()))
         return molecules[imol].get_triangles_for_blender();
   std::vector<int> v;
   return v;
}

void
molecules_container_t::make_mesh_for_bonds_for_blender(int imol, const std::string &mode) {

   if (is_valid_model_molecule(imol)) {
      // pass other params later
      molecules[imol].make_mesh_for_bonds_for_blender(mode, &geom); // bonds
   }

}

void
molecules_container_t::make_mesh_for_molecular_representation_for_blender(int imol,
                                                                          const std::string &cid,
                                                                          const std::string &colour_scheme,
                                                                          const std::string &style) {
   if (is_valid_model_molecule(imol)) {
      molecules[imol].make_mesh_for_molecular_representation_for_blender(cid, colour_scheme, style); // ribbons, etc
   }
}

void
molecules_container_t::make_mesh_for_map_contours_for_blender(int imol, float x, float y, float z, float level, float radius) {

   if (is_valid_map_molecule(imol)) {
      coot::Cartesian pos(x,y,z);
      molecules[imol].make_mesh_for_map_contours_for_blender(pos, level, radius);
   }

}

