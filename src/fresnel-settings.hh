/*
 * src/fresnel-settings.hh
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

#ifndef FRESNEL_SETTINGS_T_HH
#define FRESNEL_SETTINGS_T_HH

#include <glm/glm.hpp>

class fresnel_settings_t {
public:
   bool state;
   float bias;
   float scale;
   float power;
   glm::vec4 colour;
   fresnel_settings_t(bool state_in, const float &f1, const float &f2, const float &f3) :
      state(state_in), bias(f1), scale(f2), power(f3), colour(glm::vec4(1,1,1,1)) {}
   fresnel_settings_t() : colour(glm::vec4(1,1,1,1)) {
      state = false;
      bias  = 0.0;
      scale = 0.3;
      power = 4.0;
   }
   void update_settings(bool state_in, const float &f1, const float &f2, const float &f3) {
      state = state_in;
      bias  = f1;
      scale = f2;
      power = f3;
   }
};


#endif // FRESNEL_SETTINGS_T_HH
