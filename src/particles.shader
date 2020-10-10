
#shader vertex

#version 330 core

layout(location = 0) in vec3 position; // per vertex for the quad/hex
layout(location = 1) in vec3 normal;   // per vertex for the quad/hex
layout(location = 2) in vec4 colour;   // per vertex for the quad/hex
layout(location = 3) in vec3 instance_translation;
layout(location = 4) in vec4 instance_colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform float rotation_angle;
out vec4 colour_transfer;

void main() {

   vec3 p1 = position;
   float a = -rotation_angle; // negative so that they spin clockwise
   vec4 p2 = vec4(p1.x * cos(a) - p1.y * sin(a),
                  p1.x * sin(a) + p1.y * cos(a), 0.0, 1.0);
   mat4 trans = transpose(view_rotation);
   vec4 p3 = trans * p2;
   vec4 p4 = p3 + vec4(instance_translation, 0.0);
   gl_Position = mvp * p4;
   colour_transfer = instance_colour; // * instance_colour;

}


#shader fragment

#version 330 core

in  vec4 colour_transfer;
out vec4 outputColor;

void main() {

    outputColor = colour_transfer;

}
