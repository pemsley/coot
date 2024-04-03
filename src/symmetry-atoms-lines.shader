/*
 * src/symmetry-atoms-lines.shader
 *
 * Copyright 2021 by Medical Research Council
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


#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec4 colour;

uniform mat4 mvp;

out vec4 line_colour;

void main() {
   gl_Position = mvp * vec4(position, 1.0);
   // gl_Position = vec4(0.5 * position, 1.0);
   line_colour = colour;
}

#shader fragment

#version 330 core

uniform bool do_depth_fog;
uniform vec4 background_colour;

in vec4 line_colour;
out vec4 outputColour;

float get_fog_amount(float depth_in) {

   bool is_perspective_projection = false; // this needs to be passed to draw_symmetry() and set as a uniform

   if (! is_perspective_projection) {
      return depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d * d;
   }

}

void main() {

   if (true) {
      float fog_amount = get_fog_amount(gl_FragCoord.z);
      outputColour = mix(line_colour, background_colour, fog_amount);
   } else {
      outputColour = vec4(0,0.5,0.5,1);
   }
}
