/*
 * src/instanced-meshes-for-shadow-map.shader
 *
 * Copyright 2023 by Medical Research Council
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
// meshes-shadow-map.shader

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

layout(location = 3) in vec4 model_rotation_translation_scale_0; // instanced
layout(location = 4) in vec4 model_rotation_translation_scale_1;
layout(location = 5) in vec4 model_rotation_translation_scale_2;
layout(location = 6) in vec4 model_rotation_translation_scale_3;
layout(location = 7) in vec4 colour_instanced;

uniform mat4 mvp;
uniform mat4 view_rotation;

void main() {

   // gl_Position = mvp * vec4(position, 1.0); // simple

   mat4 model_rotation_translation_scale = mat4(model_rotation_translation_scale_0,
                                                model_rotation_translation_scale_1,
                                                model_rotation_translation_scale_2,
                                                model_rotation_translation_scale_3);
   vec3 t_pos = position;
   vec4 p4 = vec4(t_pos, 1.0);
   vec4 frag_pos = model_rotation_translation_scale * p4;
   gl_Position = mvp * frag_pos;

}

#shader fragment
// meshes-shadow-map.shader

#version 330 core

void main() {
}
