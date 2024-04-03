/*
 * src/9.ssao_blur.shader
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

// --------------- 9.ssao_blur.shader ----------------

#shader vertex

#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoords;

out vec2 TexCoords;

void main()
{
    TexCoords = aTexCoords;
    gl_Position = vec4(aPos, 1.0);
}

#shader fragment

#version 330 core
out float FragColor;

in vec2 TexCoords;

uniform sampler2D ssaoInput;

uniform int blur_size;

void main() {
    vec2 texelSize = 1.0 / vec2(textureSize(ssaoInput, 0));
    float result = 0.0;
    float n_blur = 0.0;
    for (int x = -blur_size; x <= blur_size; ++x) {
        for (int y = -blur_size; y <= blur_size; ++y) {
            vec2 offset = vec2(float(x), float(y)) * texelSize;
            float ts = texture(ssaoInput, TexCoords + offset).r;
            // This test removes the SSAO outline.
            // Perhaps use a weight w = 1 = ts ?
            if (ts != 1.0) { // don't blur with the background
               result += ts;
               n_blur += 1.0; // counting int to float
            }
        }
    }
    if (n_blur > 0)
       FragColor = result / n_blur;
    else
       FragColor = texture(ssaoInput, TexCoords).r;

    // FragColor = texture(ssaoInput, TexCoords).r;

}  
