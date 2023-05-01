
#shader vertex

#version 330 core

layout(location = 0) in vec3 position; // origin-based cylinder
layout(location = 1) in vec3 normal;   // ditto
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation; // the quaternion attached to what the mouse has done

out vec3 frag_pos_transfer;
out vec4 colour_transfer;

void main() {

   vec3 p2 = position;
   vec4 p3 = vec4(p2, 1.0);
   gl_Position = mvp * p3;
   gl_Position.z += 0.05;

   vec4 n1 = vec4(normal, 1.0);

   colour_transfer = colour;
}

#shader fragment

#version 330 core

in vec4 colour_transfer;

layout(location = 0) out vec4 outputColor;

uniform vec3 eye_position;
uniform vec4 background_colour;
uniform bool is_perspective_projection;
uniform bool do_depth_fog;


void main() {

   outputColor = vec4(0,1,0,1);

}
