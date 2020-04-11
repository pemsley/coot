
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
out vec4 eye_pos_transfer;
out float map_opacity_transfer;
out int light_0_is_on_transfer; // can this be bool?
out int light_1_is_on_transfer;
out vec4 light_0_position_transfer;
out vec4 light_1_position_transfer;
out vec4 light_0_diffuse_colour_transfer;
out vec4 light_1_diffuse_colour_transfer;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   vec4 n = view_rotation * vec4(normal, 1.0);
   Normal = normalize(n).xyz;

   bg_colour = background_colour;
   frag_pos = gl_Position.xyz;
   line_colour = colour;
   eye_pos_transfer = eye_position;
   map_opacity_transfer = map_opacity;
   light_0_is_on_transfer = light_0_is_on;
   light_1_is_on_transfer = light_1_is_on;
   light_0_position_transfer = light_0_position;
   light_1_position_transfer = light_1_position;
}

#shader fragment

#version 330 core

in vec3 frag_pos;
in vec3 Normal;
in vec4 line_colour;
in vec4 bg_colour;
in vec4 eye_pos_transfer;
in float map_opacity_transfer;
in int light_0_is_on_transfer;
in int light_1_is_on_transfer;
in vec4 light_0_position_transfer;
in vec4 light_1_position_transfer;

layout(location = 0) out vec4 out_col;


void main() {

  float specular_strength = 0.7; // 1.5 is very shiny

  vec4 specular_light_colour = vec4(0.7, 0.7, 0.7, 1.0);
  vec3 lightdir = normalize(vec3(-2,-1, 5));
  lightdir = normalize(light_0_position_transfer.xyz);
  float dp = dot(Normal, -lightdir);
  dp = max(dp, 0.0); // no negative dot products for diffuse for now, also, zero is avoided.

  float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);

  float f_1 = 1.0 - m; // because glm::ortho() near and far are reversed?
  float f_2 = 1.0 - abs(f_1 - 0.7)/0.7;
  f_2 = f_1; // remove forward fogging for now

  vec4 bg_col = bg_colour;

  vec3 eye_pos =  vec3(0.0, 0.0, 1.0);
  // vec3 eye_pos =  eye_pos_transfer.xyz;

  vec3 view_dir = eye_pos - frag_pos; // view_dir.z positive is a good idea.
  view_dir = normalize(view_dir);

  vec3 norm_2 = Normal;
  norm_2 = normalize(norm_2);
  vec3 reflect_dir = reflect(lightdir, norm_2);
  float dp_view_reflect = dot(view_dir, reflect_dir);
  dp_view_reflect = max(dp_view_reflect, 0.0);
  // when the exponent is low, the specular_strength needs to be reduced
  // a low exponent means lots of the map is specular (around the edges)
  float spec = pow(dp_view_reflect, 6.62); // shows lots of shape of the density
  // float spec = pow(dp_view_reflect, 0.4);
  vec4 specular = specular_strength * spec * specular_light_colour;

  vec4 line_colour_local = line_colour;

  vec4 col_1 = line_colour_local;  // ambient
  float ambient_strength = 0.5;
  vec4 col_2 = line_colour_local * dp;
  vec4 col_3 = col_2 + col_1 * ambient_strength + specular;
  // col_3 = specular;
  vec4 col_4 = mix(bg_col, col_3, f_2);

  // the angle between the eye direction and the triangle normal)

  // we need blend to be enabled for this to have an effect, I think.
  if (false)
     col_4.a = 2.0 * spec;

  float op = map_opacity_transfer;
  out_col = col_4;
  out_col.a = op;

}
