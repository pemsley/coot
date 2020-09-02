

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec3 atom_centre;

out vec4 colour_transfer;

void main() {

   vec3 p1 = position - atom_centre;
   mat4 t = transpose(view_rotation);
   vec4 p2 = vec4(p1, 1.0);
   vec4 p3 = t * p2;
   vec4 p4 = p3 + vec4(atom_centre, 0.0); // 0.0 is important here
   gl_Position = mvp * p4;
   colour_transfer = colour;
}


#shader fragment

#version 330 core

in vec4 colour_transfer;
layout(location = 0) out vec4 out_col;

void main() {
   out_col = colour_transfer;
}
