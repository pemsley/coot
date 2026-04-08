/*
 * src/atom-label.shader
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

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;  // quad position (around the origin)
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec3 bitangent;
layout(location = 4) in vec4 colour;
layout(location = 5) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec3 label_position;
uniform float stereo_x_scale;
uniform float stereo_x_offset;

out vec4 colour_transfer;
out vec2 texCoord_transfer;

void main() {

   // float scale = 0.00016;
   // float scale = 0.0001;
   float scale = 0.00008;
   scale = 0.0001;
   mat4 t = transpose(view_rotation);
   vec4 pos_down = scale * vec4(position, 1.0);
   vec4 p = pos_down + vec4(label_position, 1.0) * t;
   // gl_Position = mvp * t * vec4(p);
   gl_Position = mvp * t * vec4(p) * vec4(stereo_x_scale, 1.0f, 1.0f, 1.0f) + vec4(stereo_x_offset, 0.0f, 0.0f, 0.0f);
   colour_transfer = colour;
   texCoord_transfer = texCoord;

}


#shader fragment

#version 330 core

uniform sampler2D text;
uniform bool do_depth_fog;
uniform bool is_perspective_projection;
uniform vec4 background_colour;

in vec4 colour_transfer;
in vec2 texCoord_transfer;

out vec4 outputColor;

float get_fog_amount(float depth_in) {

   // text is a bit wispy already.
   // Do don't fade it out with depth too heavily

   if (! is_perspective_projection) {
      return 0.6 * depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d4;
   }

}

void main() {

   // This is for text in an image
   // vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);

   vec4 sampled_raw = texture(text, texCoord_transfer);
   vec4 sampled = vec4(1.0, 1.0, 1.0, sampled_raw.r);

   float fog_amount = 0.0;
   if (do_depth_fog)
      fog_amount = get_fog_amount(gl_FragCoord.z);
   vec4 col = colour_transfer * sampled;

   // outputColor = (1.0 - fog_amount) * col;
   // outputColor = vec4(fog_amount,0,1,1);
   float r = 1.0 - fog_amount;
   outputColor = r * col;

}
