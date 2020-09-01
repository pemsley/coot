
#shader vertex

#version 330 core

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 texCoord;

out vec2 texCoord_transfer;

uniform vec2 position;
uniform float scale;

void main() {
   float t_scale = 0.035; // .4 too big 0.3 too small
   vec2 t_position = vec2(-0.96, 0.89); // 0.9 too much, 0.88 too small
   gl_Position = vec4(scale * vertex + position, -1.0, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D text;

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {
   vec4 sampled = texture(text, texCoord_transfer);
   sampled = vec4(0.4, 0.7, 0.2, sampled.a);
   outputColor = sampled;
}
