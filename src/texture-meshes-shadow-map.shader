

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec3 bitangent;
layout(location = 4) in vec4 colour;
layout(location = 5) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

}

#shader fragment

#version 330 core

void main() {
}
