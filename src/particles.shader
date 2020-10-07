
#shader vertex

#version 330 core

layout(location = 0) in vec3 position; // per vertex for the quad/hex
layout(location = 1) in vec4 colour;   // per vertex for the quad/hex
layout(location = 2) in vec4 instance_colour;
layout(location = 3) in vec3 instance_translation;

uniform mat4 mvp;
out vec4 colour_transfer;

void main() {

   gl_Position = mvp * vec4(position + instance_translation, 1.0);
   colour_transfer = colour * instance_colour;

}


#shader fragment

#version 330 core

in  vec4 colour_transfer;
out vec4 outputColor;

void main() {

    outputColor = colour_transfer;

}
