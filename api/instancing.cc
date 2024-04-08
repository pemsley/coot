
#include<iostream>

#include "instancing.hh"


// "unwrap" the instancing
coot::simple_mesh_t
coot::instanced_mesh_to_simple_mesh(const coot::instanced_mesh_t &im) {

   coot::simple_mesh_t sm = im.markup;
   // now add to sm

   std::cout << "debug:: instanced_mesh_to_simple_mesh() geom.size() " << im.geom.size() << std::endl;

   for (unsigned int ig=0; ig<im.geom.size(); ig++) {
      const auto &g = im.geom[ig];
      std::cout << "g: " << ig << " " << g.vertices.size() << " " << g.triangles.size()
                << " A " << g.instancing_data_A.size() << " "
                << " B " << g.instancing_data_B.size() << std::endl;
      for (unsigned int iA=0; iA<g.instancing_data_A.size(); iA++) {
         std::vector<api::vnc_vertex> vertices;
         std::vector<g_triangle> triangles = g.triangles;
         unsigned int idx_base = sm.vertices.size();
         unsigned int idx_tri_base = sm.triangles.size();
         for (unsigned int iv=0; iv<g.vertices.size(); iv++) {
            auto &v_ref = g.vertices[iv];
            auto &col = g.instancing_data_A[iA].colour;
            api::vnc_vertex v(v_ref.pos * g.instancing_data_A[iA].size, v_ref.normal, col);
            v.pos += g.instancing_data_A[iA].position;
            vertices.push_back(v);
         }
         sm.vertices.insert(sm.vertices.end(), vertices.begin(), vertices.end());
         sm.triangles.insert(sm.triangles.end(), triangles.begin(), triangles.end());
         for (unsigned int k=idx_tri_base; k<sm.triangles.size(); k++)
            sm.triangles[k].rebase(idx_base);
      }

      for (unsigned int iB=0; iB<g.instancing_data_B.size(); iB++) {
         std::vector<api::vnc_vertex> vertices;
         std::vector<g_triangle> triangles = g.triangles;
         unsigned int idx_base = sm.vertices.size();
         unsigned int idx_tri_base = sm.triangles.size();
         for (unsigned int iv=0; iv<g.vertices.size(); iv++) {
            auto &v_ref   = g.vertices[iv];
            auto &col     = g.instancing_data_B[iB].colour;
            auto &size    = g.instancing_data_B[iB].size;
            auto &ori     = g.instancing_data_B[iB].orientation;
            auto &instpos = g.instancing_data_B[iB].position;
            api::vnc_vertex v(v_ref.pos * size, v_ref.normal, col);
            glm::vec4 p(v.pos, 1.0);
            glm::vec4 pr = ori * p;
            v.pos = glm::vec3(pr);
            v.pos += instpos;
            vertices.push_back(v);
         }
         sm.vertices.insert(sm.vertices.end(), vertices.begin(), vertices.end());
         sm.triangles.insert(sm.triangles.end(), triangles.begin(), triangles.end());
         for (unsigned int k=idx_tri_base; k<sm.triangles.size(); k++)
            sm.triangles[k].rebase(idx_base);
      }
   }

   std::cout << "debug:: instanced_mesh_to_simple_mesh() vertices " << sm.vertices.size()
             << " triangles " << sm.triangles.size() << std::endl;

   return sm;
}
