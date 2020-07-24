
#include "molecule-class-info.h"

#include "oct.hh"

// We can think about a more efficient interface when this one works
//
std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >
molecule_class_info_t::make_generic_vertices_for_atoms(const std::vector<glm::vec4> &index_to_colour) const {

   float sphere_radius = 0.085; // how big should atoms be?
   float radius_scale = 0.2 * bond_width; // arbs
   if (is_intermediate_atoms_molecule) radius_scale *= 1.8f;
   // atom_scale = 1.45;

   std::vector<vertex_with_rotation_translation> v1;
   std::vector<g_triangle> v2;

   bool against_a_dark_background = true;
   glm::vec3 origin(0,0,0);
   unsigned int num_subdivisions = 2; // 2 should be the default?
   if (is_intermediate_atoms_molecule)
      num_subdivisions = 1;
   float radius = 1;
   glm::vec4 col(0.5, 0.5, 0.5, 0.5);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > octaball =
      make_octasphere(num_subdivisions, origin, radius, col);

   // first count the number of atoms so that we can resize the vertices v1:
   unsigned int n_atoms = 0;
   for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++)
      n_atoms += bonds_box.consolidated_atom_centres[icol].num_points;
   v1.reserve(n_atoms * octaball.first.size());
   v2.reserve(n_atoms * octaball.second.size());
   glm::mat4 unit_matrix(1.0f);

   unsigned int atom_idx = 0; // running index so that the indices for the triangles can be calculated
   for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
      glm::vec4 atom_col = index_to_colour[icol];
      for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &ai = bonds_box.consolidated_atom_centres[icol].points[i];
         float sphere_scale = radius_scale * ai.radius_scale * 1.18;
         if (ai.is_hydrogen_atom) // this should be set already (in the generator). Is it?
            sphere_scale = 0.5;

         glm::vec3 atom_position = cartesian_to_glm(ai.position);
         unsigned int idx_base = v1.size();

         for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
            vertex_with_rotation_translation vertex(octaball.first[ibv], sphere_radius * sphere_scale);
            vertex.colour = atom_col;
            vertex.model_rotation_matrix = unit_matrix;  // for now
            vertex.model_translation = atom_position;
            v1.push_back(vertex);
         }
         std::vector<g_triangle> octaball_triangles = octaball.second;
         for (unsigned int ii=0; ii<octaball_triangles.size(); ii++)
            octaball_triangles[ii].rebase(idx_base);
         v2.insert(v2.end(), octaball_triangles.begin(), octaball_triangles.end());

         atom_idx++;
      }
   }

   // about 200,000 vertices for tutorial modern with subdivision 2
   // std::cout << "returning from make_generic_vertices_for_atoms() with sizes "
   // << v1.size() << " " << v2.size() << std::endl;

   return std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >(v1, v2);
}
