
#include <glm/gtc/matrix_transform.hpp>

#include "molecule-class-info.h"
#include "oct.hh"
#include "eyes.hh"

// We can think about a more efficient interface when this one works
//
std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >
molecule_class_info_t::make_generic_vertices_for_atoms(const std::vector<glm::vec4> &index_to_colour) const {

   bool fat_atoms_mode = false;

   float sphere_radius = 0.085; // how big should atoms be?
   float radius_scale = 0.2 * bond_width; // arbs
   if (is_intermediate_atoms_molecule) radius_scale *= 1.8f;
   if (fat_atoms_mode)
      radius_scale *= 4.0;

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
      const glm::vec4 &atom_col = index_to_colour[icol];
      // std::cout << "debug index " << icol  << " colour " << glm::to_string(atom_col)  << std::endl;
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
            vertex.model_rotation_matrix = unit_matrix; // for now
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

   if (false) { // a bit of fun

      for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
         for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {
            const graphical_bonds_atom_info_t &ai = bonds_box.consolidated_atom_centres[icol].points[i];
            float sphere_scale = radius_scale * ai.radius_scale * 1.18;
            if (ai.is_hydrogen_atom) continue;
            if (ai.is_hydrogen_atom)
               sphere_scale = 0.5;

            glm::vec3 atom_position = cartesian_to_glm(ai.position);
            float solid_theta = 0.1;
            float atom_radius = sphere_radius * radius_scale;

            float badness = fabs(ai.position.x());
            if (badness > 2.0) badness = 2.0;
            if (badness < 0.2) badness = 0.2;

            float h_scale = 1.0;
            float v_scale = badness * 0.1;

            unsigned int n_slices = 16;
            std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
               spherical_surface_circular_patch_mesh =
               make_spherical_surface_circular_patch(atom_radius * 1.189, solid_theta, h_scale, v_scale, n_slices);
            std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
               spherical_surface_circular_patch_mesh_pupil =
               make_spherical_surface_circular_patch(atom_radius * 1.191, 0.4 * solid_theta, h_scale, h_scale, n_slices);

            unsigned int idx_base = v1.size();
            unsigned int idx_base_tri = v2.size();
            float m_w = 0.6;
            float m_h = 0.1;
            m_h = badness;

            std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
               mouth_patch = make_curved_bar_patch(atom_radius * 1.29, -m_w, m_w, m_h, 16, 0.8);
            for (unsigned int i=0; i<mouth_patch.first.size(); i++) {
               const position_normal_vertex &v = mouth_patch.first[i];
               glm::vec4 c(0.2, 0.2, 0.2, 1.0);

               vertex_with_rotation_translation vrt(v.pos, v.normal, c);
               glm::mat4 rm(1.0f);
               glm::mat3 mm = glm::mat3(rm);
               vrt.model_rotation_matrix = mm;
               vrt.model_translation = atom_position;
               v1.push_back(vrt);
            }
            v2.insert(v2.end(), mouth_patch.second.begin(), mouth_patch.second.end());
            for (unsigned int i=idx_base_tri; i<v2.size(); i++)
               v2[i].rebase(idx_base);

            for (int ii=0; ii<2; ii++) {

               unsigned int idx_base = v1.size();
               unsigned int idx_base_tri = v2.size();

               for (unsigned int i=0; i<spherical_surface_circular_patch_mesh.first.size(); i++) {
                  const position_normal_vertex &v = spherical_surface_circular_patch_mesh.first[i];
                  float whiteness = 0.82;
                  glm::vec4 c(whiteness, whiteness, whiteness, 1.0);
                  vertex_with_rotation_translation vrt(v.pos, v.normal, c);
                  glm::mat4 rm(1.0f);
                  float angle = 0.5 * static_cast<float>(2 * ii - 1);
                  rm = glm::rotate(rm, angle, glm::vec3(0,1,0));
                  glm::mat3 mm = glm::mat3(rm);
                  vrt.model_rotation_matrix = mm;
                  vrt.model_translation     = atom_position;
                  v1.push_back(vrt);
               }
               v2.insert(v2.end(),
                         spherical_surface_circular_patch_mesh.second.begin(),
                         spherical_surface_circular_patch_mesh.second.end());
               for (unsigned int i=idx_base_tri; i<v2.size(); i++)
                  v2[i].rebase(idx_base);

               idx_base = v1.size();
               idx_base_tri = v2.size();
               for (unsigned int i=0; i<spherical_surface_circular_patch_mesh_pupil.first.size(); i++) {
                  const position_normal_vertex &v = spherical_surface_circular_patch_mesh_pupil.first[i];
                  glm::vec4 c(0.1, 0.1, 0.1, 1.0);
                  vertex_with_rotation_translation vrt(v.pos, v.normal, c);
                  glm::mat4 rm(1.0f);
                  float angle = 0.5 * static_cast<float>(2 * ii - 1) + 0.1;
                  rm = glm::rotate(rm, angle, glm::vec3(0,1,0));
                  glm::mat3 mm = glm::mat3(rm);
                  vrt.model_rotation_matrix = mm;
                  vrt.model_translation = atom_position + 0.0f * glm::vec3(ii+1,ii+1,ii+1);
                  v1.push_back(vrt);
               }
               v2.insert(v2.end(),
                         spherical_surface_circular_patch_mesh_pupil.second.begin(),
                         spherical_surface_circular_patch_mesh_pupil.second.end());
               for (unsigned int i=idx_base_tri; i<v2.size(); i++)
                  v2[i].rebase(idx_base);
            }
         }
      }
   }

   return std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >(v1, v2);
}
