/*
 * src/ligand-view.shader
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

layout(location = 0) in vec2 position;

uniform float aspect_ratio;

void main() {

   // shader offsets
   float x_off =  0.02;
   float y_off = -0.20;

   vec2 scaled_pos = 0.05 * position;
   scaled_pos.x /= aspect_ratio;
   vec2 offset_pos = scaled_pos + vec2(-0.6, -0.6) + vec2(x_off, y_off);
   gl_Position = vec4(offset_pos, -1.0, 1.0);

}

#shader fragment

#version 330 core

// uniform vec4 colour;

out vec4 out_col;

void main() {

   // pass the colour as a uniform
   out_col = vec4(0.5, 0.5, 0.5, 1.0);
}
