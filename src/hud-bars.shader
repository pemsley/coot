
#shader vertex

#version 330 core

// this shader uses instancing: The basic quad has the vertices
// and each block has colour and position/displacement (which places the
// bar on the screen)

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec4 colour;
layout (location = 2) in vec2 position_offset;
layout (location = 3) in float scale;

out vec4 colour_transfer;

void main()
{
   // we want the tooltip to appear "over" this bar. Maybe
   // I could turn of depth cueing for that? Hmm. Not done at the moment
   gl_Position = vec4(scale * vertex.x + position_offset.x,
                      vertex.y + position_offset.y, -0.999, 1.0);
   colour_transfer = colour;
}

#shader fragment

#version 330 core

in vec4 colour_transfer;
out vec4 colour;

void main()
{
   colour = colour_transfer;
}


