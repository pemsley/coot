
#shader vertex

#version 330 core

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 texCoord;

out vec2 texCoord_transfer;

uniform vec2 position;
uniform vec2 scales;

// This shader is for textures

void main() {

   // Note to self: the text of the tooltip needs to go over the
   // background of the tooltip

   vec2 scaled_vertices = vertex * scales; // vec2(0.1, 0.05);
   gl_Position = vec4(scaled_vertices + position , -0.999, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D text;

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   bool this_is_the_hud_bar_labels = false; // pass this as a uniform

   vec4 sampled = texture(text, texCoord_transfer);
   sampled = sampled + sampled;

   if (this_is_the_hud_bar_labels) {
      sampled = vec4(0.4, 0.7, 0.2, sampled.a);
   } else {
      if ((sampled.r + sampled.g + sampled.b) < 0.01) {
         sampled.a = 0.0;
      } else {
         sampled.a = 1.0;
      }
   }
   outputColor = sampled;
}
