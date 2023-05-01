
// rename this "hud-image-texture.shader"
// because it is a shader for HUD images,
// not text.

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
   gl_Position = vec4(scaled_vertices + position , -1.0, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D text; // change this confusing name - "image_texture"

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   bool this_is_the_hud_bar_labels = false; // pass this as a uniform

   vec4 sampled = texture(text, texCoord_transfer);
   // sampled = vec4(text_colour.r, text_colour.r, text_colour.r, sampled.r);
   outputColor = sampled;

   if (outputColor.a < 0.5) discard;

   // outputColor.a = 0.9; // why did I have this? For the rama underlying distribution? Hmm.
                           // OK I guess I need a uniform for that if I'm going to use this shader
                           // for that.
}
