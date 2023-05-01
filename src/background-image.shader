
#shader vertex
// --- background-image.shader

#version 330 core

layout (location = 0) in vec2 vertex;   // -1 to +1
layout (location = 1) in vec2 texCoords; //  0 to 1
out vec2 texCoords_transfer;

void main() {
   vec4 p1 = vec4(vertex, 0.099, 1.0); // check the depth
   gl_Position = p1;
   texCoords_transfer = texCoords;
}


#shader fragment
// --- background-image.shader

#version 330 core

uniform sampler2D image_texture;
in vec2 texCoords_transfer;

out vec4 outputColor;

void main() {
   vec4 t = texture(image_texture, texCoords_transfer);
   outputColor = t;
}
