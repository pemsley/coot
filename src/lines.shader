

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;

out vec4 colour_transfer;

void main() {

   gl_Position = mvp * vec4(position, 1.0);
   colour_transfer = colour;
}


#shader fragment

#version 330 core

in vec4 colour_transfer;
layout(location = 0) out vec4 out_col;

void main() {
   out_col = colour_transfer;
}
