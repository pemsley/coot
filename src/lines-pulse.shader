/*
 * src/lines-pulse.shader
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

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec3 atom_centre;

out vec4 colour_transfer;

void main() {

   mat4 t = transpose(view_rotation);
   vec4 p2 = vec4(position, 1.0);
   vec4 p3 = t * p2;
   vec4 p4 = p3 + vec4(atom_centre, 0.0); // 0.0 is important here
   gl_Position = mvp * p4;
   colour_transfer = colour;
}


#shader fragment

#version 330 core

in vec4 colour_transfer;
layout(location = 0) out vec4 out_col;

void main() {
   out_col = colour_transfer;
}
