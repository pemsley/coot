/*
 * src/rama-plot-phi-psi-markers.shader
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

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 texCoord;
layout (location = 2) in vec2 positions; // instanced (e.g. phi, psi positions)

out vec2 texCoord_transfer;

uniform vec2 position;
uniform vec2 scales;

uniform vec2 window_resize_position_correction;
uniform vec2 window_resize_scales_correction;

// This shader is for textures, uses HUDTextureMesh with the
// draw_instances() draw method.
//
// It is used for Rama-plot points

void main() {

   vec2 scaled_vertices = (vertex + positions) * scales;
   vec2 p1 = scaled_vertices + position;
   vec2 p2 = p1 * window_resize_scales_correction;
   vec2 p3 = p2 + window_resize_position_correction;
   // vec2 p3 = p1;
   gl_Position = vec4(p3, -1.0, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D text; // change this name to "image_texture"

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   vec4 sampled = texture(text, texCoord_transfer);
   outputColor = sampled;
   if (outputColor.a < 0.6)
      discard;
}
