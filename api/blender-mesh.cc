
#include <iostream>
#include <functional> // for std::hash

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "blender-mesh.hh"

// static
size_t
coot::blender_mesh_t:: make_colour_hash(const glm::vec4 &col) {

   // Convert the components to integers and combine them
   int r = static_cast<int>(col.r *  255);
   int g = static_cast<int>(col.g *  255);
   int b = static_cast<int>(col.b *  255);
   int a = static_cast<int>(col.a *  255);

   // Combine the components into a single integer
   size_t combined = (r <<  24) | (g <<  16) | (b <<  8) | a;

   // Use std::hash to create a hash from the combined integer
   return std::hash<size_t>{}(combined);
}


coot::blender_mesh_t::blender_mesh_t(const coot::simple_mesh_t &mesh) {

   auto colour_from_index_triple = [&mesh] (const g_triangle &tri) {
      unsigned int idx = tri.point_id[0];
      return mesh.vertices[idx].color;
   };

   std::map<size_t, int > colour_hash_map; // hash code to colour index

   // colour table
   for (unsigned int i=0; i<mesh.triangles.size(); i++) {
      glm::vec4 col = colour_from_index_triple(mesh.triangles[i]);
      size_t hash = make_colour_hash(col);
      std::map<size_t, int >::const_iterator it = colour_hash_map.find(hash);
      if (it == colour_hash_map.end()) {
         int new_idx = colour_hash_map.size();
         colour_hash_map[hash] = new_idx;
         colour_table[new_idx] = col;
      }
   }

   vertices.resize(mesh.vertices.size());
   for (unsigned int i=0; i<mesh.vertices.size(); i++) {
      vertices[i] = mesh.vertices[i].pos;
   }

   normals.resize(mesh.vertices.size());
   for (unsigned int i=0; i<mesh.vertices.size(); i++) {
      normals[i] = mesh.vertices[i].normal;
   }

   triangles.resize(mesh.triangles.size());
   for (unsigned int i=0; i<mesh.triangles.size(); i++) {
      glm::vec4 col = colour_from_index_triple(mesh.triangles[i]);
      size_t hash = make_colour_hash(col);
      int colour_index = colour_hash_map[hash]; // should always find
      blender_triangle_t bt(mesh.triangles[i], colour_index);
      triangles[i] = bt;
   }

}

coot::blender_mesh_t::blender_mesh_t(const coot::instanced_mesh_t &im) {



   std::map<size_t, int > colour_hash_map; // hash code to colour index

   auto get_colour_index = [] (const glm::vec4 &col,
                               std::map<size_t, int > *colour_hash_map_p,
                               std::map<int, glm::vec4> *colour_map_p) {
      int idx = -1;
      size_t hash = make_colour_hash(col);
      std::map<size_t, int >::const_iterator it = colour_hash_map_p->find(hash);
      if (it == colour_hash_map_p->end()) {
         int new_idx = colour_hash_map_p->size();
         // std::cout << "assigning hash colour index " << new_idx << " hash " << hash << " col " << glm::to_string(col)  << std::endl;
         (*colour_hash_map_p)[hash] = new_idx;
         (*colour_map_p)[new_idx] = col;
         idx = new_idx;
      } else {
         idx = it->second;
      }
      return idx;
   };

   // ------------------------------ instanced geometry -----------------------------------

   for (unsigned int ig=0; ig<im.geom.size(); ig++) {

      const auto &instanced_geom(im.geom[ig]);

      // it's one or the other
      const std::vector<instancing_data_type_A_t> &dA = instanced_geom.instancing_data_A;
      const std::vector<instancing_data_type_B_t> &dB = instanced_geom.instancing_data_B;

      for (unsigned int iB=0; iB<dB.size(); iB++) {

         std::vector<glm::vec3> l_vertices;
         std::vector<glm::vec3> l_normals;
         const auto &col = dB[iB].colour;
         int colour_index = get_colour_index(col, &colour_hash_map, &colour_table);
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         const auto &size    = dB[iB].size;
         const auto &ori     = dB[iB].orientation;
         const auto &instpos = dB[iB].position;
         for (unsigned int iv=0; iv<instanced_geom.vertices.size(); iv++) {
            const api::vn_vertex &v_ref = instanced_geom.vertices[iv];
            glm::vec4 p(v_ref.pos * size, 1.0);
            glm::vec4 pr = ori * p;
            glm::vec3 pos = glm::vec3(pr);
            pos += instpos;
            l_vertices.push_back(pos);
            glm::vec4 n4(v_ref.normal, 1.0);
            glm::vec4 n4r = ori * n4;
            glm::vec3 nr(n4r);
            l_normals.push_back(nr);
         }

         vertices.insert(vertices.end(), l_vertices.begin(), l_vertices.end());

         std::vector<blender_triangle_t> l_triangles; // = instanced_geom.triangles;
         for (unsigned int itr=0; itr<instanced_geom.triangles.size(); itr++) {
            blender_triangle_t t(instanced_geom.triangles[itr], colour_index);
            l_triangles.push_back(t);
         }
         triangles.insert(triangles.end(), l_triangles.begin(), l_triangles.end());
         for (unsigned int k=idx_tri_base; k<triangles.size(); ++k)
            triangles[k].rebase(idx_base);
      }

      for (unsigned int iA=0; iA<dA.size(); iA++) {
         std::vector<glm::vec3> l_vertices;
         std::vector<glm::vec3> l_normals;
         const auto &col = dA[iA].colour;
         int colour_index = get_colour_index(col, &colour_hash_map, &colour_table);
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         for (unsigned int iv=0; iv<instanced_geom.vertices.size(); iv++) {
            glm::vec3 vert(instanced_geom.vertices[iv].pos);
            vert *= dA[iA].size;
            vert += dA[iA].position;
            l_vertices.push_back(vert);
            l_normals.push_back(instanced_geom.vertices[iv].normal);
         }

         std::vector<blender_triangle_t> l_triangles; // = instanced_geom.triangles;
         for (unsigned int it=0; it<instanced_geom.triangles.size(); it++) {
            blender_triangle_t t(instanced_geom.triangles[it], colour_index);
            l_triangles.push_back(t);
         }

         vertices.insert(vertices.end(), l_vertices.begin(), l_vertices.end());

         triangles.insert(triangles.end(), l_triangles.begin(), l_triangles.end());
         for (unsigned int k=idx_tri_base; k<triangles.size(); k++)
            triangles[k].rebase(idx_base);
      }
   }

   // ------------------------------ simple mesh geometry -------------------------------


   // ... fill later

}
