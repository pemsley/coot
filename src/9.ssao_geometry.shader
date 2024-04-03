/*
 * src/9.ssao_geometry.shader
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

#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;

out vec3 frag_pos_transfer;
out vec3 normal_transfer;

uniform bool invertedNormals;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    vec4 viewPos = view * model * vec4(aPos, 1.0);
    frag_pos_transfer = viewPos.xyz;
    
    mat3 normalMatrix = transpose(inverse(mat3(view * model)));
    normal_transfer = normalize(normalMatrix * (invertedNormals ? -aNormal : aNormal));

    gl_Position = projection * viewPos;
}


#shader fragment

#version 330 core

// actually, only gPosition and gNormal are used.
// so scrap gAlbedo
layout (location = 0) out vec3 gPosition;
layout (location = 1) out vec3 gNormal;

in vec3 frag_pos_transfer;
in vec3 normal_transfer;

void main() {
   // store the fragment position vector in the first gbuffer texture
   gPosition = frag_pos_transfer;
   // also store the per-fragment normals into the gbuffer
   gNormal = normal_transfer;

   gPosition = vec3(1,1,0);
   gNormal   = vec3(0,0,1);
   
}
