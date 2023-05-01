
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;
layout(location = 3) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;
out vec4 colour_transfer;
out vec2 texCoord_transfer;

void main() {

   float scale = 0.0001; // for "real" atom label
   scale = 1.0;
   mat4 t = transpose(view_rotation);
   gl_Position = mvp * t * vec4(scale * position, 1.0);
   colour_transfer = colour;
   texCoord_transfer = texCoord;

}


#shader fragment

#version 330 core

uniform sampler2D text;

in vec4 colour_transfer;
in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   // This is for text in an image
   // vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);

   vec4 sampled = texture(text, texCoord_transfer);
   sampled = vec4(0.6, 0.6, 0.1, 1-sampled.r);

   // sampled = vec4(sampled.r, sampled.r, sampled.r, 1.0);

   // sampled = texture(text, texCoord_transfer);
   //  outputColor = colour_transfer;

   outputColor = sampled;

}
