
#include <glm/gtc/matrix_transform.hpp>

#include "molecule-class-info.h"
#include "oct.hh"
#include "eyes.hh"

// We can think about a more efficient interface when this one works
//
std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >
molecule_class_info_t::make_generic_vertices_for_atoms(const std::vector<glm::vec4> &index_to_colour,
                                                       float atom_radius_scale_factor) const {


   std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> > v;
   std::vector<vertex_with_rotation_translation> &v1 = v.first;
   std::vector<g_triangle> &v2 = v.second;

   // this is not consistent with the bonds - the hydrogen atoms are too small
   float sphere_radius = 0.085; // how big should atoms be?
   float radius_scale = 0.2 * bond_width; // arbs
   if (is_intermediate_atoms_molecule) radius_scale *= 1.8f;

   radius_scale *= atom_radius_scale_factor;

   bool against_a_dark_background = true;
   glm::vec3 origin(0,0,0);
   unsigned int num_subdivisions = 2; // 2 should be the default?
   if (is_intermediate_atoms_molecule)
      num_subdivisions = 1;
   float radius = 1;
   glm::vec4 col(0.5, 0.5, 0.5, 1.0);
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

         // std::cout << "sphere_scale " << sphere_scale << std::endl;
         if (sphere_scale > 3.5)
            sphere_scale = 3.5; // hacketty-hack for now

         if (ai.is_hydrogen_atom) // this should be set already (in the generator). Is it?
            sphere_scale *= 0.6;

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

   if (false) {
      if (is_intermediate_atoms_molecule) {
         unsigned int idx_base = v1.size();
         unsigned int idx_base_tri = v2.size();
         std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> > v_fun =
            fun(4.0 * radius_scale);
         std::cout << "fun triangles " << v1.size() << " " << v2.size() << std::endl;
         v1.insert(v1.end(), v_fun.first.begin(), v_fun.first.end());
         v2.insert(v2.end(), v_fun.second.begin(), v_fun.second.end());
         for (unsigned int ii=idx_base_tri; ii<v2.size(); ii++)
            v2[ii].rebase(idx_base);
      }
   }

   return v;
}


std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >
molecule_class_info_t::fun(float radius_scale) const {

   // about 200,000 vertices for tutorial modern with subdivision 2
   // std::cout << "returning from make_generic_vertices_for_atoms() with sizes "
   // << v1.size() << " " << v2.size() << std::endl;

   std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> > v;
   std::vector<vertex_with_rotation_translation> &v1 = v.first;
   std::vector<g_triangle> &v2 = v.second;

   for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
      for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &ai = bonds_box.consolidated_atom_centres[icol].points[i];
         float sphere_scale = radius_scale * ai.radius_scale * 1.18;
         if (ai.is_hydrogen_atom) continue;
         if (ai.is_hydrogen_atom)
            sphere_scale = 0.5;

         glm::vec3 atom_position = cartesian_to_glm(ai.position);
         float solid_theta = 0.1;
         float sphere_radius = 0.085; // how big should atoms be?
         float atom_radius = sphere_radius * radius_scale;

         float badness = fabs(ai.position.x());
         if (badness > 2.0) badness = 2.0;
         if (badness < 0.2) badness = 0.2;

         float h_scale = 1.0;
         float v_scale = 1.0;

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
            mouth_patch = make_curved_bar_patch(atom_radius * 0.29, -m_w, m_w, m_h, 16, 0.8);
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
   return v;
}


std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >
molecule_class_info_t::make_generic_vertices_for_bad_CA_CA_distances() const {

   std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> > vp;

   std::vector<vertex_with_rotation_translation> &v1 = vp.first;
   std::vector<g_triangle> &v2 = vp.second;
   glm::vec3 origin(0,0,0);
   unsigned int num_subdivisions = 2;
   if (is_intermediate_atoms_molecule)
      return vp;

   float radius = 1;
   float sphere_radius = 0.13;
   float sphere_scale = 1.0;
   glm::vec4 col(0.99, 0.55, 0.1, 1.0);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > octaball =
      make_octasphere(num_subdivisions, origin, radius, col);
   glm::mat4 unit_matrix(1.0f);

   for (int i=0; i<bonds_box.n_bad_CA_CA_dist_spots; i++) {
      unsigned int idx_base = v1.size();
      glm::vec3 position(bonds_box.bad_CA_CA_dist_spots_ptr[i].x(),
                          bonds_box.bad_CA_CA_dist_spots_ptr[i].y(),
                          bonds_box.bad_CA_CA_dist_spots_ptr[i].z());
      for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
         vertex_with_rotation_translation vertex(octaball.first[ibv], sphere_radius * sphere_scale);
         vertex.colour = col;
         vertex.model_rotation_matrix = unit_matrix; // for now
         vertex.model_translation = position;
         v1.push_back(vertex);
      }
      std::vector<g_triangle> octaball_triangles = octaball.second;
      for (unsigned int ii=0; ii<octaball_triangles.size(); ii++)
         octaball_triangles[ii].rebase(idx_base);
      v2.insert(v2.end(), octaball_triangles.begin(), octaball_triangles.end());
   }

   return std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >(v1, v2);
}


// rama balls.
std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >
molecule_class_info_t::make_generic_vertices_for_rama_balls(float ball_scale_factor,
                                                            const glm::vec3 &screen_up_dir) const {

   // treat them like atoms for now.

   std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> > v;
   std::vector<vertex_with_rotation_translation> &v1 = v.first;
   std::vector<g_triangle> &v2 = v.second;
   float rama_ball_pos_offset_scale = 0.6;

   unsigned int num_subdivisions = 2;
   glm::vec3 origin(0,0,0);
   glm::mat4 unit_matrix(1.0f);
   float radius = 0.62f * ball_scale_factor;
   glm::vec4 col(0.5, 0.5, 0.5, 1.0);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > octaball =
      make_octasphere(num_subdivisions, origin, radius, col);

   for (int i=0; i<bonds_box.n_ramachandran_goodness_spots; i++) {
      coot::Cartesian position = bonds_box.ramachandran_goodness_spots_ptr[i].first;
      const float &prob_raw    = bonds_box.ramachandran_goodness_spots_ptr[i].second;
      double prob(prob_raw);
      if (prob > 0.5) prob = 0.5; // 0.4 and 2.5 f(for q) might be better (not tested)
      double q = (1 - 2.0 * prob);
      q = pow(q, 20);
      coot::colour_holder col = coot::colour_holder(q, 0.0, 1.0, false, std::string(""));
      glm::vec3 atom_position = cartesian_to_glm(position) + rama_ball_pos_offset_scale * screen_up_dir;
      unsigned int idx_base = v1.size();
      for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
         vertex_with_rotation_translation vertex(octaball.first[ibv], radius);
         vertex.colour = glm::vec4(col.red, col.green, col.blue, 1.0f);
         vertex.model_rotation_matrix = unit_matrix;
         vertex.model_translation = atom_position;
         v1.push_back(vertex);
      }
      std::vector<g_triangle> octaball_triangles = octaball.second;
      for (unsigned int ii=0; ii<octaball_triangles.size(); ii++)
         octaball_triangles[ii].rebase(idx_base);
      v2.insert(v2.end(), octaball_triangles.begin(), octaball_triangles.end());
   }

   return v;
}
