/*
 * src/camera-facing-quad-shader-for-testing.shader
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
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;
layout(location = 3) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;
out vec4 colour_transfer;
out vec2 texCoord_transfer;

void main() {

   float scale = 0.0001; // for "real" atom label
   scale = 1.0;
   mat4 t = transpose(view_rotation);
   gl_Position = mvp * t * vec4(scale * position, 1.0);
   colour_transfer = colour;
   texCoord_transfer = texCoord;

}


#shader fragment

#version 330 core

uniform sampler2D text;

in vec4 colour_transfer;
in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   // This is for text in an image
   // vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);

   vec4 sampled = texture(text, texCoord_transfer);
   sampled = vec4(0.6, 0.6, 0.1, 1-sampled.r);

   // sampled = vec4(sampled.r, sampled.r, sampled.r, 1.0);

   // sampled = texture(text, texCoord_transfer);
   //  outputColor = colour_transfer;

   outputColor = sampled;

}
