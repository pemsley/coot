/*
 * src/Material.hh
 *
 * Copyright 2020 by Medical Research Council
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

#ifndef COOT_MATERIAL_HH
#define COOT_MATERIAL_HH

#include <glm/glm.hpp>

class Material {
public:
   glm::vec4 ambient;
   glm::vec4 diffuse;
   glm::vec4 specular;
   bool do_specularity;
   float shininess;
   float specular_strength;
   Material(const glm::vec4 &ambient,
            const glm::vec4 &diffuse,
            const glm::vec4 &specular,
            const float &s) : ambient(ambient), diffuse(diffuse), specular(specular), shininess(s) {
      specular_strength = 1.0;
      do_specularity = false;
   }
   Material() :
      ambient( glm::vec4(0.2, 0.2, 0.2, 1.0)),
      diffuse( glm::vec4(0.5, 0.5, 0.5, 1.0)),
      specular(glm::vec4(0.5, 0.5, 0.5, 1.0)) {
      do_specularity = false;
      shininess = 64.0;
      specular_strength = 0.4;
   }
   void turn_specularity_on(bool state) {
      do_specularity = state;
   }
};


#endif // COOT_MATERIAL_HH
