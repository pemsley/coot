
// rename this "hud-image-texture.shader"
// because it is a shader for HUD images,
// not text.

#shader vertex

#version 330 core

layout (location = 0) in vec2 vertex;   // -1 to +1
layout (location = 1) in vec2 texCoord; //  0 to 1

out vec2 texCoord_transfer;

// This shader is for textures - there is no instancing.
// 20210908-PE uses the same window resize scale and offset correction as is in rama-plot-phi-psi-markers.shader

uniform vec2 position;
uniform vec2 scales;

uniform vec2 window_resize_position_correction;
uniform vec2 window_resize_scales_correction;

uniform bool relative_to_right;
uniform bool relative_to_top;

void main() {

   vec2 scaled_vertices = vertex  * scales;
   vec2 p1 = scaled_vertices + position;
   vec2 p2 = p1 * window_resize_scales_correction;
   vec2 p3 = p2;
   if (relative_to_right) p3.x += 1.0;
   if (relative_to_top)   p3.y += 1.0;
   vec2 p4 = p3 + window_resize_position_correction;

   gl_Position = vec4(p4, 0.0, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D image_texture;

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   bool this_is_the_hud_bar_labels = false; // pass this as a uniform

   vec4 sampled = texture(image_texture, texCoord_transfer);

   outputColor = sampled;
   // outputColor = vec4(0.8 * sampled.r, 0.8 * sampled.r, 0.8 * sampled.r, 1.0);

   if (outputColor.a < 0.5) discard; // for text as textures rendering

   // outputColor.a = 0.9; // why did I have this? For the rama underlying distribution? Hmm.
                           // OK I guess I need a uniform for that if I'm going to use this shader
                           // for that.
}
