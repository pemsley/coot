
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

// these act on the whole mesh, not per-button/bar
uniform vec2 scales; // because the window has been widened, say.
uniform vec2 offset_position; // not the same as position_offset!

out vec4 colour_transfer;

void main()
{
   // we want the tooltip to appear "over" this bar. Maybe
   // I could turn off depth test for that? Hmm. Not done at the moment
   //

   // These 2 adjustments are now uniforms.
   // float plot_scale = 0.5;
   // vec2 plot_offset = vec2(-0.85, -0.85);

   // vec2 p1 = vec2(vertex.x * scales.x * scales.x, vertex.y * scale_y * scales.y);

   // this gives a bottom right box that scales as the window is reshaped.
   // Not too bad
   // vec2 p1 = vec2(vertex.x * scale_x, vertex.y * scale_y);
   // vec2 p2 = p1 + position_offset; // + offset_position;
   // gl_Position = vec4(p2, -1.0, 1.0);

   vec2 p1 = vec2(vertex.x * scale_x, vertex.y * scale_y);
   vec2 p2 = p1 + position_offset;
   vec2 p3 = p2 * scales;
   vec2 p4 = p3 + offset_position;

   gl_Position = vec4(p4, -1.0, 1.0);
   colour_transfer = colour;

}

#shader fragment

#version 330 core

in vec4 colour_transfer;
out vec4 out_colour;

void main()
{
   out_colour = colour_transfer;
}

