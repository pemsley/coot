/*
 * src/draw-map-cap.shader
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
out vec4 colour_transfer;
out vec3 normal_transfer;

void main() {
   gl_Position = mvp * vec4(position, 1.0);
   colour_transfer = colour;
   normal_transfer = normal;
}

#shader fragment
#version 330 core

in vec4 colour_transfer;
in vec3 normal_transfer;
out vec4 output_colour;

void main() {

   vec3 lightdir = normalize(vec3(1,0,1));
   float dp = dot(normal_transfer, -lightdir);
   dp = clamp(dp, 0.0, 1.0);

   output_colour = colour_transfer + dp * colour_transfer;

}
