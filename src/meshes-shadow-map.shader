

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

}

#shader fragment

#version 330 core

void main() {
}
