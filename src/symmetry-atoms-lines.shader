

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec4 colour;

uniform mat4 mvp;

out vec4 line_colour;

void main() {
   gl_Position = mvp * vec4(position, 1.0);
   // gl_Position = vec4(0.5 * position, 1.0);
   line_colour = colour;
}

#shader fragment

#version 330 core

uniform bool do_depth_fog;
uniform vec4 background_colour;

in vec4 line_colour;
out vec4 outputColour;

float get_fog_amount(float depth_in) {

   bool is_perspective_projection = false; // this needs to be passed to draw_symmetry() and set as a uniform

   if (! is_perspective_projection) {
      return depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d * d;
   }

}

void main() {

   if (true) {
      float fog_amount = get_fog_amount(gl_FragCoord.z);
      outputColour = mix(line_colour, background_colour, fog_amount);
   } else {
      outputColour = vec4(0,0.5,0.5,1);
   }
}
