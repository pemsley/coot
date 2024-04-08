/*
 * src/background-image.shader
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

#shader vertex
// --- background-image.shader

#version 330 core

layout (location = 0) in vec2 vertex;   // -1 to +1
layout (location = 1) in vec2 texCoords; //  0 to 1
out vec2 texCoords_transfer;

void main() {
   vec4 p1 = vec4(vertex, 0.099, 1.0); // check the depth
   gl_Position = p1;
   texCoords_transfer = texCoords;
}


#shader fragment
// --- background-image.shader

#version 330 core

uniform sampler2D image_texture;
in vec2 texCoords_transfer;

out vec4 outputColor;

void main() {
   vec4 t = texture(image_texture, texCoords_transfer);
   outputColor = t;
}
