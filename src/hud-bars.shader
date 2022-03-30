
#shader vertex

#version 330 core

// this shader uses instancing: The basic quad has the vertices
// and each block has colour and position/displacement (which places the
// bar on the screen)
//
// This shader is also used for HUD buttons

layout (location = 0) in vec2 vertex;
layout (location = 1) in float shade;
layout (location = 2) in vec4 colour;  // this and below are instanced
layout (location = 3) in vec2 position_offset;
layout (location = 4) in float scale_x;
layout (location = 5) in float scale_y;

out vec4 colour_transfer;

// culture clash! should the bars/buttons each know their offset position (and scales)
// (and they will stretch if the window is widened)
uniform vec2 position_offset_as_uniform;
uniform vec2 scales_as_uniform;

uniform vec2 window_resize_position_correction;
uniform vec2 window_resize_scales_correction;

uniform bool relative_to_right;
uniform bool relative_to_top;

void main()
{

   vec2 scaled_vertices = vec2(vertex.x * scale_x + position_offset.x,
                               vertex.y * scale_y + position_offset.y);
   vec2 p1 = scaled_vertices;

   // 20220329-PE I can't get this to work. I think that it works differently to hud-image-texture.shader
   // so comment out the corrections for now.
   vec2 p2 = p1; // * window_resize_scales_correction;
   vec2 p3 = p2; // + window_resize_position_correction;
   // if (relative_to_right) p3.x += 1.0;
   // if (relative_to_top)   p3.y += 1.0;
   vec2 p4 = p3;

   // we want the tooltip to appear "over" this bar. Maybe
   // I could turn off depth test for that? Hmm. Not done at the moment

   gl_Position = vec4(p4, -0.999, 1.0);

   // was:
   // gl_Position = vec4(scale_x * vertex.x + position_offset.x,
   // scale_y * vertex.y + position_offset.y,
   //                       -0.999, 1.0);

   // we adds colour adjust so that the buttons look smooth shaded
   //
   float c = 0.24;
   float s3 = shade * shade * shade * shade * shade;
   vec4 colour_adjust = vec4(s3 * c, s3 * c, s3 * c, 0.1);
   colour_transfer = colour + colour_adjust;
}

#shader fragment

#version 330 core

in vec4 colour_transfer;
out vec4 colour;

void main()
{
   colour = colour_transfer;
}
