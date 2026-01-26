/*
 * src/outline-of-active-residue.shader
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

layout(location = 0) in vec3 position; // origin-based cylinder
layout(location = 1) in vec3 normal;   // ditto
layout(location = 2) in vec4 colour;

uniform mat4 mvp;

out vec3 frag_pos_transfer;
out vec4 colour_transfer;

void main() {

   vec3 p2 = position;
   vec4 p3 = vec4(p2, 1.0);
   gl_Position = mvp * p3;
   gl_Position.z += 0.05;

   vec4 n1 = vec4(normal, 1.0);

   colour_transfer = colour;
}

#shader fragment

#version 330 core

in vec4 colour_transfer;

layout(location = 0) out vec4 outputColor;

uniform vec3 eye_position;
uniform vec4 background_colour;
uniform bool is_perspective_projection;
uniform bool do_depth_fog;


void main() {

   outputColor = vec4(0,1,0,1);

}
