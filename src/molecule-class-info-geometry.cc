
#include "molecule-class-info.h"
#include "utils/dodec.hh"

// We can think about a more efficient interface when this one works
//
std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >
molecule_class_info_t::make_generic_vertices_for_atoms(const std::vector<glm::vec4> &index_to_colour) const {

   float sphere_radius = 0.105; // how big should atoms be?
   float radius_scale = 0.2 * bond_width; // arbs
   if (is_intermediate_atoms_molecule) radius_scale *= 1.8f;
   // atom_scale = 1.45;

   std::vector<vertex_with_rotation_translation> v1;
   std::vector<g_triangle> v2;

   bool against_a_dark_background = true;
   pentakis_dodec d(1.0);
   unsigned int atom_idx = 0; // running index so that the indices for the triangles can be calculated
   for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
      // set_bond_colour_by_mol_no(icol, against_a_dark_background); // this function is const ATM.
      for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {

         const graphical_bonds_atom_info_t &ai = bonds_box.consolidated_atom_centres[icol].points[i];
         float sphere_scale = radius_scale * ai.radius_scale * 1.18;

         if (ai.is_hydrogen_atom) // this should be set already (in the generator). Is it?
            sphere_scale = 0.5;

         if (true) {
            const coot::Cartesian &at_pos = ai.position;

            for (unsigned int i=0; i<20; i++) { // dodec vertices
               const clipper::Coord_orth &pt = d.d.get_point(i);
               vertex_with_rotation_translation gv;
               gv.model_rotation_matrix = glm::mat3(1.0f); // identity
               gv.model_translation = glm::vec3(sphere_radius * sphere_scale * pt.x(),
                                                sphere_radius * sphere_scale * pt.y(),
                                                sphere_radius * sphere_scale * pt.z());
               gv.pos = glm::vec3(at_pos.x(), at_pos.y(), at_pos.z());
               gv.normal = glm::normalize(glm::vec3(gv.model_translation));
               gv.colour = index_to_colour[icol];
               v1.push_back(gv);
            }

            for (unsigned int i=0; i<12; i++) { // dodec faces
               const clipper::Coord_orth &pv = d.pyrimid_vertices[i];
               vertex_with_rotation_translation gv;
               gv.model_rotation_matrix = glm::mat3(1.0f); // identity
               gv.model_translation = glm::vec3(sphere_radius * sphere_scale * pv.x(),
                                                sphere_radius * sphere_scale * pv.y(),
                                                sphere_radius * sphere_scale * pv.z());
               gv.pos = glm::vec3(at_pos.x(), at_pos.y(), at_pos.z());
               gv.normal = glm::normalize(glm::vec3(gv.model_translation));
               gv.colour = index_to_colour[icol];
               v1.push_back(gv);
            }

            // make the pyrimid triangles, with the tips being the pyramid vertices
            for (unsigned int i=0; i<12; i++) {

               // pyramid vertices are at the centre of every face
               std::vector<unsigned int> f_indices = d.d.face(i);

               for(unsigned int j=0; j<5; j++) {
                  unsigned int j_next = j+1;
                  if (j == 4) j_next = 0;
                  unsigned int t_idx_1 = 32 * atom_idx + f_indices[j];
                  unsigned int t_idx_2 = 32 * atom_idx + f_indices[j_next];
                  unsigned int t_idx_3 = 32 * atom_idx + 20 + i;
                  g_triangle t(t_idx_1, t_idx_2, t_idx_3);
                  v2.push_back(t);
               }
            }
            atom_idx++;
         }
      }
   }
   return std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> >(v1, v2);
}
