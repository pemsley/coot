
#include <iostream>

#include "simple-mesh.hh"
#include "oct.hh"

std::string
coot::simple_mesh_t::vandt() const {

   std::string s;
   s += " n-vertices: ";
   s += std::to_string(vertices.size());
   s += " n-triangles: ";
   s += std::to_string(triangles.size());
   return s;
}

void
coot::simple_mesh_t::translate(const glm::vec3 &t) {

   for (auto &item : vertices)
      item.pos += t;
}

void
coot::simple_mesh_t::add_submesh(const simple_mesh_t &submesh) {

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangles.size();
   vertices.insert(vertices.end(), submesh.vertices.begin(), submesh.vertices.end());
   triangles.insert(triangles.end(), submesh.triangles.begin(), submesh.triangles.end());
   for (unsigned int i=idx_base_tri; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);
   
}

//! if the colour map is empty then go through the vector of vertices finding colours and putting them
//! into a colour table. This is for Blender - where the colour are assigned to a Material, and a Material
//! is assigned to a face.
void
coot::simple_mesh_t::fill_colour_map() {

#if 0 // 20230110-PE having a colour index in g_triangle messes the setup of the GL indexing buffers
      // So let's keep it out for now - and use the class that is derived from it when we want to use
      // index colours for blender

   // fill this: std::map<int, glm::vec4> colour_index_to_colour_map;

   if (! colour_index_to_colour_map.empty()) return;

   std::map<float, unsigned int> reverse_colour_map;
   unsigned int colour_index = 0;

   for (const auto &vert : vertices) {
      float colour_hash = 1000.0f * vert.color[0] + 100.0f * vert.color[1] + 10.0f * vert.color[2] + vert.color[3];
      if (reverse_colour_map.find(colour_hash) == reverse_colour_map.end()) {
         reverse_colour_map[colour_hash] = colour_index;
         colour_index_to_colour_map[colour_index] = vert.color;
         colour_index++;
      }
   }

   std::cout << "in fill_colour_map() found " << colour_index << " colours" << std::endl;

   unsigned int found_correctly = 0;
   unsigned int fail_type_1 = 0;
   unsigned int fail_type_2 = 0;

   for (auto &tri : triangles) {
      const auto &v_0 = vertices[tri[0]];
      const auto &v_1 = vertices[tri[1]];
      const auto &v_2 = vertices[tri[2]];
      float colour_hash_0 = 1000.0f * v_0.color[0] + 100.0f * v_0.color[1] + 10.0f * v_0.color[2] + v_0.color[3];
      float colour_hash_1 = 1000.0f * v_1.color[0] + 100.0f * v_1.color[1] + 10.0f * v_1.color[2] + v_1.color[3];
      float colour_hash_2 = 1000.0f * v_2.color[0] + 100.0f * v_2.color[1] + 10.0f * v_2.color[2] + v_2.color[3];
      if (colour_hash_0 == colour_hash_1) {
         if (colour_hash_0 == colour_hash_2) {
            tri.colour_index = reverse_colour_map[colour_hash_0];
            found_correctly++;
         } else {
            // std::cout << "tri with different colours 0 2" << std::endl;
            glm::vec4 colour_sum = v_0.color + v_1.color + v_2.color;
            glm::vec4 col_av = colour_sum * 0.33333333f;
            float colour_hash_av = 1000.0f * col_av[0] + 100.0f * col_av[1] + 10.0f * col_av[2] + col_av[3];
            if (reverse_colour_map.find(colour_hash_av) == reverse_colour_map.end()) {
               reverse_colour_map[colour_hash_av] = colour_index;
               tri.colour_index = colour_index;
               colour_index_to_colour_map[colour_index] = col_av;
               colour_index++;
            }
            fail_type_1++;
         }
      } else {
         // std::cout << "tri with different colours 0 1" << std::endl;
         glm::vec4 colour_sum = v_0.color + v_1.color + v_2.color;
         glm::vec4 col_av = colour_sum * 0.33333333f;
         float colour_hash_av = 1000.0f * col_av[0] + 100.0f * col_av[1] + 10.0f * col_av[2] + col_av[3];
         if (reverse_colour_map.find(colour_hash_av) == reverse_colour_map.end()) {
            reverse_colour_map[colour_hash_av] = colour_index;
            tri.colour_index = colour_index;
            colour_index_to_colour_map[colour_index] = col_av;
            colour_index++;
         }
         fail_type_2++;
      }
   }

   std::cout << "post-colour-fill: correct " << found_correctly << " fail-type-1 " << fail_type_1 << " fail_type_2 " << fail_type_2
             << std::endl;
#endif

}

// static
coot::simple_mesh_t
coot::simple_mesh_t::make_sphere() {

   simple_mesh_t m;
   unsigned int num_subdivisions = 2;
   glm::vec3 position(0,0,0);
   float radius = 1.0f;
   glm::vec4 colour(0.5, 0.5, 0.5, 1.0);
   std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > octaball =
      make_octasphere(num_subdivisions, position, radius, colour);

   m.vertices  = octaball.first;
   m.triangles = octaball.second;
   return m;
}

void
coot::simple_mesh_t::change_colour(const glm::vec4 &c) {

   for (auto &v : vertices)
      v.color = c;
}

// scale before translation of course
void
coot::simple_mesh_t::scale(float scale_factor) {

   for (auto &v : vertices)
      v.pos *= scale_factor;

}

