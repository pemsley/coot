
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;  // quad position (around the origin)
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;
layout(location = 3) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec3 label_position;

out vec4 colour_transfer;
out vec2 texCoord_transfer;

void main() {

   float scale = 0.00016;
   mat4 t = transpose(view_rotation);
   vec4 pos_down = scale * vec4(position, 1.0);
   vec4 p = pos_down + vec4(label_position, 1.0) * t;
   gl_Position = mvp * t * vec4(p);
   colour_transfer = colour;
   texCoord_transfer = texCoord;

}


#shader fragment

#version 330 core

uniform sampler2D text;
uniform bool do_depth_fog;
uniform bool is_perspective_projection;
uniform vec4 background_colour;

in vec4 colour_transfer;
in vec2 texCoord_transfer;

out vec4 outputColor;

float get_fog_amount(float depth_in) {

   // text is a bit wispy already.
   // Do don't fade it out with depth too heavily
   
   if (! is_perspective_projection) {
      return 0.6 * depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d4;
   }

}

void main() {

   // This is for text in an image
   // vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);

   vec4 sampled = texture(text, texCoord_transfer);
   sampled = vec4(1.0, 1.0, 1.0, sampled.r);
   float fog_amount = 0.0;
   if (do_depth_fog)
      fog_amount = get_fog_amount(gl_FragCoord.z);
   vec4 col = colour_transfer * sampled;
   outputColor = (1.0 - fog_amount) * col;

}
