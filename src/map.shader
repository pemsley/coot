
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec4 background_colour;
uniform vec4 eye_position;
uniform float map_opacity;
uniform int light_0_is_on;
uniform int light_1_is_on;
uniform vec4 light_0_position;
uniform vec4 light_1_position;
uniform vec4 light_0_diffuse_colour;
uniform vec4 light_1_diffuse_colour;

out vec3 frag_pos;
out vec3 Normal;
out vec4 line_colour;
out vec4 bg_colour;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   vec4 n = view_rotation * vec4(normal, 1.0);
   Normal = normalize(n).xyz;

   bg_colour = background_colour;
   frag_pos = gl_Position.xyz;
   line_colour = colour;
}

#shader fragment

#version 330 core

in vec3 frag_pos;
in vec3 Normal;
in vec4 line_colour;
in vec4 bg_colour;

uniform vec4 eye_pos;
uniform float map_opacity;
uniform int light_0_is_on;
uniform int light_1_is_on;
uniform vec4 light_0_position;
uniform vec4 light_1_position;

layout(location = 0) out vec4 out_col;

float get_fog_amount(float depth_in) {

   // make this a uniform
   bool perspective_projection = false;
   if (! perspective_projection) {
      return depth_in;
   } else {
      // needs tweaking
      return depth_in * depth_in;
   }

}


void main() {

  float specular_strength = 0.3; // 1.5 is very shiny
  vec4 specular_light_colour = vec4(0.7, 0.7, 0.7, 1.0);

  // a light direction of 0,0,1 is good for fresnelly outlining
  vec3 lightdir = normalize(vec3(-2,-1, -1));

  // using the light_0_position give to the shader in a uniform
  // gives us rotating lights
  // lightdir = normalize(light_0_position.xyz);
  float dp = dot(Normal, -lightdir);
  dp = max(dp, 0.0); // no negative dot products for diffuse for now, also, zero is avoided.

  float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);

  float f_1 = 1.0 - m; // because glm::ortho() near and far are reversed?
  float fog_amount = get_fog_amount(gl_FragCoord.z);

  vec4 bg_col = bg_colour;

  // vec3 eye_pos =  vec3(0.0, 0.0, 1.0);
  vec3 eye_pos_3 =  eye_pos.xyz;

  vec3 view_dir = eye_pos_3 - frag_pos; // view_dir.z positive is a good idea.
  view_dir = normalize(view_dir);

  vec3 norm_2 = Normal;
  norm_2 = normalize(norm_2);
  vec3 reflect_dir = reflect(lightdir, norm_2);
  float dp_view_reflect = dot(view_dir, reflect_dir);
  dp_view_reflect = max(dp_view_reflect, 0.0);
  // when the exponent is low, the specular_strength needs to be reduced
  // a low exponent means lots of the map is specular (around the edges)
  float spec = pow(dp_view_reflect, 6.62);
  vec4 specular = specular_strength * spec * specular_light_colour;

  vec4 line_colour_local = line_colour;

  vec4 col_1 = line_colour_local;  // ambient
  float ambient_strength = 0.6;
  vec4 col_2 = line_colour_local * dp;
  // vec4 col_3 = 0.6 * col_2 + col_1 * ambient_strength;
  vec4 col_3 = col_2 + col_1 * ambient_strength + specular;
  vec4 col_4 = mix(col_3, bg_col, fog_amount);

  float op = map_opacity;
  out_col = col_4;
  // we need GL_BLEND to be enabled for this to have an effect.
  out_col.a = op;

}
