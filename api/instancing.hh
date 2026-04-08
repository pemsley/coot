#ifndef INSTANCING_HH
#define INSTANCING_HH

#include <vector>
#include "coot-utils/vertex.hh"
#include "coot-utils/g_triangle.hh"
#include "coot-utils/simple-mesh.hh"

namespace coot {

   //! class for A type instancing data - this does not contain an orientation matrix
   class instancing_data_type_A_t {
   public:
      glm::vec3 position;
      glm::vec4 colour;
      glm::vec3 size;
      instancing_data_type_A_t(const glm::vec3 &position_in, const glm::vec4 &colour_in, const glm::vec3 &size_in) :
         position(position_in), colour(colour_in), size(size_in) {}
      instancing_data_type_A_t() {}
   };

   //! class for B type instancing data - this _does_ contain an orientation matrix
   class instancing_data_type_B_t {
   public:
      glm::vec3 position; // origin offset
      glm::vec4 colour;
      glm::vec3 size;
      //! the orientation matrix rotates the vector away from "z is up"
      glm::mat4 orientation; // 3 sets of vec3 in the shader. 20230114-PE not glm::mat4
      instancing_data_type_B_t(const glm::vec3 &position_in, const glm::vec4 &colour_in, const glm::vec3 &size_in, const glm::mat4 &ori) :
         position(position_in), colour(colour_in), size(size_in), orientation(ori) {}
      instancing_data_type_B_t() {}
   };

   //! instancing container for vertices, triangles and instancing data
   class instanced_geometry_t {
   public:
      //! vertices (containing positions and normals)
      std::vector<api::vn_vertex> vertices;
      //! triangle indices
      std::vector<g_triangle> triangles;
      std::string name;

      instanced_geometry_t() {}
      explicit instanced_geometry_t(const std::string &n) : name(n) {};
      instanced_geometry_t(const std::vector<api::vn_vertex> &v, const std::vector<g_triangle> &t) :
         vertices(v), triangles(t) {}
      bool empty() { return instancing_data_A.empty() && instancing_data_B.empty(); }

      //! a vector of type A instancing
      std::vector<instancing_data_type_A_t> instancing_data_A;
      //! a vector of type B instancing
      std::vector<instancing_data_type_B_t> instancing_data_B;
   };

   //! a simple container for instancing containers - for multiple geometries that are instanced.
   //! (e.g. balls and sticks).
   class instanced_mesh_t {
   public:
      instanced_mesh_t() {};
      std::vector<instanced_geometry_t> geom;
      void add(const instanced_geometry_t &ig) { geom.push_back(ig); }

      //! cis-peptide markup can't be drawn instanced
      simple_mesh_t markup;

      //! message
      std::string message; // if this is non-empty then the mesh generation
                           // failed and the user should see this messaage.

      //! clear
      void clear() { geom.clear(); markup.clear(); }

      void export_to_glTF(const std::string &file_name, bool use_binary_format) const;
   };

   // convert for export. A better exporter would preserve the instancing (but each ball colour
   // would need it's own reference I think) - the instancing is the orientation matrix (only).
   simple_mesh_t instanced_mesh_to_simple_mesh(const instanced_mesh_t &im);

}

#endif // INSTANCING_HH
