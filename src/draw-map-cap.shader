
#shader vertex
#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
out vec4 colour_transfer;
out vec3 normal_transfer;

void main() {
   gl_Position = mvp * vec4(position, 1.0);
   colour_transfer = colour;
   normal_transfer = normal;
}

#shader fragment
#version 330 core

in vec4 colour_transfer;
in vec3 normal_transfer;
out vec4 output_colour;

void main() {

   vec3 lightdir = normalize(vec3(1,0,1));
   float dp = dot(normal_transfer, -lightdir);
   dp = clamp(dp, 0.0, 1.0);

   output_colour = colour_transfer + dp * colour_transfer;

}
