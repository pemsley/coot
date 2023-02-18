#ifndef SIMPLE_MESH_HH
#define SIMPLE_MESH_HH

#include <vector>
#include <map>
#include "vertex.hh"
#include "g_triangle.hh"

namespace coot {

   //! The basic mesh for transfering mesh geometry and colours
   class simple_mesh_t {
   public:
      std::vector<api::vnc_vertex> vertices;
      std::vector<g_triangle> triangles;
      std::string name;
      //! constructor (for vectors)
      simple_mesh_t() {}
      //! constructor with name
      explicit simple_mesh_t(const std::string &name_in) : name(name_in) {}
      simple_mesh_t(const std::vector<api::vnc_vertex> &vertices_in,
                    const std::vector<g_triangle> &triangles_in) : vertices(vertices_in), triangles(triangles_in) {}
      void translate(const glm::vec3 &t);
      // 20221101-PE blender uses colours/materials for faces. So let's store those too.
      // Now each face (each g_triangle) can have a colour_index (default is -1 (unset)).
      // Maybe a std::vector would be a better/faster container.
      std::map<int, glm::vec4> colour_index_to_colour_map;

      //! utilty function
      void add_submesh(const simple_mesh_t &submesh);

      //! if the colour map is empty then go through the vector of vertices findling colours and putting them
      //! into a colour table. This is for Blender - where the colour are assigned to a Material, and a Material
      //! is assigned to a face
      void fill_colour_map();

      //! export to glTF
      void export_to_gltf(const std::string &file_name, bool use_binary_format) const;

      void set_name(const std::string &n) { name = n; }

      static simple_mesh_t make_sphere();
      void scale(float scale_factor);
      void change_colour(const glm::vec4 &c);
      // for debugging - the string for the number of vertices and triangles
      std::string vandt() const;

   };
}

#endif // SIMPLE_MESH_HH
