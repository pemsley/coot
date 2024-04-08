/*
 * src/hud-labels.shader
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

// rename this "hud-image-texture.shader"
// because it is a shader for HUD images,
// not text.

#shader vertex

#version 330 core

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 texCoord;

out vec2 texCoord_transfer;

uniform vec2 position;
uniform vec2 scales;

// This shader is for textures

void main() {

   // Note to self: the text of the tooltip needs to go over the
   // background of the tooltip

   vec2 scaled_vertices = vertex * scales; // vec2(0.1, 0.05);
   gl_Position = vec4(scaled_vertices + position , -1.0, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D text; // change this confusing name - "image_texture"

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   bool this_is_the_hud_bar_labels = false; // pass this as a uniform

   vec4 sampled = texture(text, texCoord_transfer);
   // sampled = vec4(text_colour.r, text_colour.r, text_colour.r, sampled.r);
   outputColor = sampled;

   // 20231003-PE this is a quick fix. It should really use the background colour
   // and be a box with nice margins around the letters - but that would mean another
   // shader. This will do for now.
   //
   // if (outputColor.a < 0.5) discard;
   if (outputColor.a < 0.5) outputColor = vec4(0.0, 0.0, 0.0, 1.0);

   // outputColor.a = 0.9; // why did I have this? For the rama underlying distribution? Hmm.
                           // OK I guess I need a uniform for that if I'm going to use this shader
                           // for that.
}
