
#shader vertex

#version 330 core

// layout(location = 0) in vec3 position; // per vertex for the quad/hex
// layout(location = 1) in vec3 normal;   // per vertex for the quad/hex
// layout(location = 2) in vec4 colour;   // per vertex for the quad/hex
// layout(location = 3) in vec2 texCoord;
// layout(location = 4) in vec3 instance_translation;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec3 bitangent;
layout(location = 4) in vec4 colour;
layout(location = 5) in vec2 texCoord;
layout(location = 6) in vec3 instance_translation;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform float canvas_scale; // 0.8 for happy faces, 0.2 for anchored/fixed atoms
out vec2 texCoord_transfer;

void main() {

   // currently this code correctly places the vertices in
   // perspective mode, but not in orthographic. Baah!
   //
   mat4 t = transpose(view_rotation);
   vec4 p2 = vec4(canvas_scale * position, 1.0); // 1.0 is important here
   vec4 p3 = t * p2;
   vec4 p4 = p3 + vec4(instance_translation, 0.0); // 0.0 is important here

   gl_Position = mvp * p4;
   // gl_Position = vec4(0.000000000000000000005 * position, 0.0);
   texCoord_transfer = texCoord;

}


#shader fragment

#version 330 core

in vec2 texCoord_transfer;
uniform sampler2D face;
uniform float opacity;
uniform vec4 background_colour;
uniform bool is_perspective_projection;

out vec4 outputColor;

float get_fog_amount(float depth_in) {

   if (! is_perspective_projection) {
      return depth_in * depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d4;
   }
}


void main() {

   vec4 sampled = texture(face, texCoord_transfer);

   if (sampled.a < 0.1)
      discard;

   float fog_amount = get_fog_amount(gl_FragCoord.z);
   sampled.a *= opacity; // opacity not correct for anchored atoms marker?
   outputColor = sampled;

   outputColor = mix(outputColor, background_colour, fog_amount);

}

