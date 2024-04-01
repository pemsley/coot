/*
 * coot-utils/vertex.hh
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
#ifndef VERTEX_HH
#define VERTEX_HH

#include <glm/glm.hpp>

// for standard objects at the origin - typically used in instancing
namespace coot {

   namespace api {

      //! a vertex with (only) postion and normal. Useful for instancing, perhaps
      class vn_vertex {
      public:
         glm::vec3 pos;
         glm::vec3 normal; // normalized on input
         vn_vertex(const glm::vec3 &pos_in,
                   const glm::vec3 &norm_in) :
            pos(pos_in), normal(norm_in) {}
         vn_vertex() {}
      };

      //! a vertex with postion, normal and colour
      class vnc_vertex {
      public:
         glm::vec3 pos;
         glm::vec3 normal; // normalized on input
         glm::vec4 color;  // or colour?
         vnc_vertex(const glm::vec3 &pos_in,
                    const glm::vec3 &norm_in,
                    const glm::vec4 &col_in) : pos(pos_in), normal(norm_in), color(col_in) {}
         explicit vnc_vertex(const vn_vertex &vn) :
            pos(vn.pos), normal(vn.normal), color(glm::vec4(0.5, 0.5, 0.5, 1.0)) {}
         vnc_vertex(const vn_vertex &vn, const glm::vec4 &c) : pos(vn.pos), normal(vn.normal), color(c) {}
         vnc_vertex() {}
      };

      //! 20230108-PE copied from generic-vertex.hh (then edited).
      //!
      //! vertex with rotation and translation (e.g. for oriented bonds)
      class vertex_with_rotation_translation {
      public:
         glm::mat3 model_rotation_matrix; // orientation
         glm::vec3 model_translation; // the coordinates of the first atom of the bond
         glm::vec3 pos;
         glm::vec3 normal; // normalized when set
         glm::vec4 colour;
         //! constructor
         vertex_with_rotation_translation(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &c) : pos(p), normal(n), colour(c) {}
         //! constructor
         vertex_with_rotation_translation(const vnc_vertex &v, const glm::vec3 &atom_position, float scale) :
            model_rotation_matrix(glm::mat3(1.0f)), model_translation(atom_position),
            pos(v.pos * scale), normal(v.normal), colour(v.color) {}
         //! constructor
         vertex_with_rotation_translation() {}
      };
   }
}


#endif // VERTEX_HH
