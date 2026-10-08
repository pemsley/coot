/*
 * coot-utils/simple-mesh.hh
 *
 * Copyright 2022 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef SIMPLE_MESH_HH
#define SIMPLE_MESH_HH

#include <vector>
#include <map>
#include <string>
#include "vertex.hh"
#include "g_triangle.hh"
#include "cylinder.hh"

namespace coot {

   //! The basic mesh for transfering mesh geometry and colours
   class simple_mesh_t {
   public:
      //! check this status before using a `simple_mesh_t`. 1 is good, 0 is bad (0 is set when when we get a bad_alloc)
      int status;
      //! vertices
      std::vector<api::vnc_vertex> vertices;
      //! vertex index triples
      std::vector<g_triangle> triangles;
      //! mesh name
      std::string name;

      //! Which part of the model each vertex belongs to, for meshes that know.
      //!
      //! Several of the MoleculesToTriangles primitives already record, per vertex, the atom
      //! that generated it - a surface patch comes from one atom's sphere or from a torus
      //! between two, a ribbon vertex from the residue it passes through. Until now that was
      //! dropped when the primitive was flattened into a mesh, because there was nowhere for
      //! it to go. This is that place, one entry per vertex, indexing `owners` below.
      //!
      //! It lets a caller answer "which residue is under the pointer" from a triangle
      //! intersection alone, and light exactly that residue's vertices.
      //!
      //! Either empty, meaning the mesh does not know, or exactly vertices.size() long. A
      //! vertex whose owner could not be determined carries no_owner.
      std::vector<unsigned int> vertex_owner;
      //! A second part the same vertex also belongs to, or no_owner.
      //!
      //! Some vertices genuinely belong to two: the saddle a probe sweeps between a pair of
      //! atoms is shared by both of them, and it is in exactly those grooves that one
      //! residue's territory gives way to the next. Naming only one owner there puts the
      //! boundary on the far rim of the groove instead of down the middle of it.
      //!
      //! Either empty or the same length as vertex_owner.
      std::vector<unsigned int> vertex_owner_other;
      //! How much of the vertex belongs to vertex_owner, the rest belonging to
      //! vertex_owner_other. 1 where there is only one owner.
      //!
      //! Either empty or the same length as vertex_owner.
      std::vector<float> vertex_owner_weight;
      //! The parts themselves, as CID selection strings, indexed by vertex_owner.
      //!
      //! Strings rather than residue_spec_t so that this header, which is included very
      //! widely and currently depends on nothing above coot-utils, does not acquire a
      //! dependency on geometry/ for the sake of a label.
      std::vector<std::string> owners;
      //! The value in vertex_owner for a vertex with no known owner.
      //! constexpr, not const: add_submesh passes it to vector::resize, which takes a reference,
      //! and a plain static const member would then need an out-of-line definition to link.
      static constexpr unsigned int no_owner = 0xFFFFFFFFu;

      //! Does this mesh know where its vertices came from? True only if the arrays are the
      //! right size to be trusted, so a caller can use this in place of its own checks.
      bool has_owners() const { return !owners.empty() && vertex_owner.size() == vertices.size(); }
      //! constructor (for vectors)
      simple_mesh_t() : status(1) {}
      //! constructor with name
      explicit simple_mesh_t(const std::string &name_in) : status(1), name(name_in) {}
      simple_mesh_t(const std::vector<api::vnc_vertex> &vertices_in,
                    const std::vector<g_triangle> &triangles_in) : status(1), vertices(vertices_in), triangles(triangles_in) {}
      explicit simple_mesh_t(const cylinder &cyl) : vertices(cyl.vertices), triangles(cyl.triangles) {
         status = 1;
      }
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
      void export_to_gltf(const std::string &file_name,
                          float roughness_factor, float smoothness_factor,
                          bool use_binary_format) const;

      void set_name(const std::string &n) { name = n; }

      static simple_mesh_t make_sphere();
      void scale(float scale_factor);
      void change_colour(const glm::vec4 &c);
      void clear();
      // for debugging - the string for the number of vertices and triangles
      std::string vandt() const;

   };
}

#endif // SIMPLE_MESH_HH
