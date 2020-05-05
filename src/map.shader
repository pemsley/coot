
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec4 background_colour;
uniform vec4 eye_position;
uniform int light_0_is_on;
uniform int light_1_is_on;
uniform vec4 light_0_position;
uniform vec4 light_1_position;
uniform vec4 light_0_diffuse_colour;
uniform vec4 light_1_diffuse_colour;
uniform bool is_perspective_projection;

out vec3 frag_pos;
out vec3 Normal;
out vec4 line_colour;
// out vec4 bg_colour;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   vec4 n = view_rotation * vec4(normal, 1.0);
   Normal = normalize(n).xyz;

   // bg_colour = background_colour;
   frag_pos = gl_Position.xyz;
   line_colour = colour;
}

#shader fragment

#version 330 core

in vec3 frag_pos;
in vec3 Normal;
in vec4 line_colour;

layout(location = 0) out vec4 out_col;

uniform vec4 eye_position;
uniform float map_opacity;
uniform vec4 background_colour;
uniform int light_0_is_on;
uniform int light_1_is_on;
uniform vec4 light_0_position;
uniform vec4 light_1_position;
uniform bool is_perspective_projection;


float get_fog_amount(float depth_in) {

   if (! is_perspective_projection) {
      return sqrt(depth_in);
   } else {
      // needs tweaking
      float d = depth_in * depth_in;
      return d * d;
   }

}


void main() {

  float specular_strength = 0.00004; // 1.5 is very shiny
  vec4 specular_light_colour = vec4(0.7, 0.7, 0.7, 1.0);

  // a light direction of 0,0,1 is good for fresnelly outlining (well, it used to be)
  // negative y means light from above
  vec3 lightdir = normalize(vec3(-2, -3, -2));

  // using the light_0_position given to the shader in a uniform
  // gives us rotating lights in perspective projection
  // but not in orthographic!
  lightdir = normalize(-light_0_position.xyz);
  float dp = dot(Normal, -lightdir);
  dp = max(dp, 0.0); // no negative dot products for diffuse for now, also, zero is avoided.

  float m = clamp(gl_FragCoord.z, 0.0f, 1.0f);
  float fog_amount = get_fog_amount(m);

  vec4 bg_col = background_colour; // needed?

  vec3 norm_2 = Normal;
  norm_2 = normalize(norm_2);
  vec3 reflect_dir = reflect(lightdir, norm_2);

  vec4 colour_local = line_colour;

  vec4 col_1 = colour_local;  // ambient
  float ambient_strength = 0.8;
  float diffuse_strength = 1.2;
  vec4 col_2 = diffuse_strength * colour_local * dp;
  vec4 col_3 = col_2 + col_1 * ambient_strength;
  vec4 col_4 = mix(col_3, bg_col, fog_amount);

  out_col = col_4;
  out_col.a = map_opacity;


}
