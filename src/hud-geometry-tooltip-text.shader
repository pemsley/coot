/*
 * src/hud-geometry-tooltip-text.shader
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

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 texCoord;

out vec2 texCoord_transfer;

uniform vec2 position;
uniform vec2 scales;

// This shader is for textures of text, the output colour
// depends on text_colour and the red value sampled for the
// text texture. Not a good fit for just showing a HUD image

void main() {

   // Note to self: the text of the tooltip needs to go over the
   // background of the tooltip
   // 20220214-PE actually for the moment, I have removed the background, so that we see ust "raw" text.

   vec2 scaled = vertex * scales;
   gl_Position = vec4(scaled + position, -1.0, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D text;
uniform vec4 text_colour; // only use the first 3 of course.

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   vec4 sampled = texture(text, texCoord_transfer);
   sampled = vec4(text_colour.r, text_colour.g, text_colour.b, sampled.r);
   // sampled = vec4(0.2, 0.7, 0.9, sampled.r);

   outputColor = sampled;
}
