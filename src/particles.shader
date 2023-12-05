
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

   float scale = 0.3;
   float sc_inst_trans = 1.0;
   // scale = 0.3 // for crows
   // sc_inst_trans = 0.1 // for crows

   // vec3 position_scaled = 0.0000000000000000002 * position;
   vec3 position_scaled = 0.2 * position;
   vec3 p1 = position_scaled * scale;
   //   float a = -rotation_angle; // negative so that they spin clockwise
   //   vec4 p2 = vec4(p1.x * cos(a) - p1.y * sin(a),
   //                  p1.x * sin(a) + p1.y * cos(a), 0.0, 1.0);
   vec4 p2 = vec4(p1, 1.0);
   mat4 trans = transpose(view_rotation);
   vec4 p3 = 0.1 * trans * p2;
   vec4 p4 = p3 + vec4(sc_inst_trans * instance_translation, 0.0);

   gl_Position = mvp * p4;
   gl_Position = vec4(position, 1.0);
   gl_Position = mvp * vec4((position + instance_translation), 1.0);

   vec4 p5 = vec4(0.0000000000000000001 * position, 0.0);
   p5 = vec4(position, 1.0);
   vec4 p6 = p5 + vec4(instance_translation, 0.0);
   gl_Position = mvp * p6;

   // gl_Position = vec4(0.0000000000000000002 * position, 0.0);

   colour_transfer = instance_colour; // * instance_colour;

   vec3 pp = position;
   if (position.x == 0.0) pp.x = 0.1;
   if (position.y == 0.0) pp.y = 0.1;
   colour_transfer = vec4(pp, 1.0);
   if (instance_translation.x != 0) colour_transfer = vec4(0,0,1,1);
   if (instance_translation.y != 0) colour_transfer = vec4(0,0,1,1);
   if (instance_translation.z != 0) colour_transfer = vec4(0,0,1,1);
}


#shader fragment

#version 330 core

in  vec4 colour_transfer;
out vec4 outputColor;

void main() {

    outputColor = colour_transfer;
    // outputColor = vec4(0.8, 0.1, 0.8, 1.0) + 0.001 * colour_transfer;
}
