#ifndef COOT_API_BLENDER_MESH_HH
#define COOT_API_BLENDER_MESH_HH

#include "instancing.hh"

namespace coot {

   class blender_triangle_t {
   public:
      blender_triangle_t() : colour_index(-1) {};
      blender_triangle_t(const g_triangle &t, int col_idx) : triangle(t), colour_index(col_idx) {}
      g_triangle triangle;
      int colour_index;

      void rebase(const unsigned int &idx_base) {
         triangle.rebase(idx_base);
      }
   };

   class blender_mesh_t {
      static size_t make_colour_hash(const glm::vec4 &col);
   public:
      blender_mesh_t() {}
      explicit blender_mesh_t(const instanced_mesh_t &im);
      explicit blender_mesh_t(const simple_mesh_t &sm);
      std::map<int, glm::vec4> colour_table; // include alpha
      std::vector<glm::vec3> vertices;
      std::vector<glm::vec3> normals;
      std::vector<blender_triangle_t> triangles;
   };

   // the array/lists that we need to send to Blender Python are:
   // PyList_SetItem(r_py, 0, vertices_py);                  // get_vertices_for_blender()
   // PyList_SetItem(r_py, 1, tris_py);                      // get_triangles_for_blender()
   // PyList_SetItem(r_py, 2, face_colours_dict_py);         // get_colour_table_for_blender()

   // where
   // vertices_py:          a list of (x,y,z)
   // tris_py:              a 4-member (int) list of vertex indices for a triangle and a colour_index (int)
   // face_colours_dict_py: a list of index and 3 members (idx, (r,g,b))
}


#endif // BLENDER_MESH_HH
