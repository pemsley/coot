
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec4 background_colour;
uniform vec4 eye_position;

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

layout(location = 0) out vec4 out_col;


void main() {

  float specular_strength = 0.0000; // I don't like specular/phong map lights at the moment.

  vec4 light_colour = vec4(0.4, 0.4, 0.4, 1.0);
  vec3 lightdir = normalize(vec3(-2,-1,5));
  float dp = dot(Normal, -lightdir);

  float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);

  float f_1 = 1.0 - m; // because glm::ortho() near and far are reversed?
  float f_2 = 1.0 - abs(f_1 - 0.7)/0.7;
  f_2 = f_1; // remove forward fogging for now

  float flat_frac = 0.5;

  vec4 bg_col = bg_colour;

  vec3 eye_pos =  vec3(0.0, 0.0, 5.0);

  vec3 view_dir = eye_pos - frag_pos; // view_dir.z positive is a good idea.
  view_dir = normalize(view_dir);

  vec3 norm_2 = Normal;
  norm_2 += normalize(norm_2);
  vec3 reflect_dir = reflect(-lightdir, norm_2);
  float dp_view_reflect = dot(view_dir, reflect_dir);
  dp_view_reflect = max(dp_view_reflect, 0.0);
  float spec = pow(dp_view_reflect, 6.2);
  vec4 specular = specular_strength * spec * light_colour;

  vec4 line_colour_local = line_colour;
  line_colour_local.w = 0.3;

  vec4 col_1 = line_colour_local;
  vec4 col_2 = line_colour_local * dp;
  vec4 col_3 = mix(col_2, col_1, flat_frac) + specular;
  vec4 col_4 = mix(bg_col, col_3, f_2);

  out_col = col_4;

}
