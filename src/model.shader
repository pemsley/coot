
#shader vertex

#version 330 core

layout(location = 0) in vec3 model_rotation_matrix_0;
layout(location = 1) in vec3 model_rotation_matrix_1;
layout(location = 2) in vec3 model_rotation_matrix_2;
layout(location = 3) in vec3 model_translation;
layout(location = 4) in vec3 position; // origin-based cylinder
layout(location = 5) in vec3 normal;   // ditto
layout(location = 6) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation; // the quaternion attached to what the mouse has done

out vec3 frag_pos;
out vec3 Normal;
out vec4 tri_color;

void main() {

   mat3 model_rotation_matrix = mat3(model_rotation_matrix_0, model_rotation_matrix_1, model_rotation_matrix_2);

   vec3 p2 = position * model_rotation_matrix;
   vec4 p3 = vec4(p2 + model_translation, 1.0);
   gl_Position = mvp * p3;

   vec4 n1 = vec4(normal * model_rotation_matrix, 1.0);
   vec4 n2 = view_rotation * n1;

   Normal = normalize(n2.xyz);
   frag_pos =  p3.xyz;
   tri_color = colour;
}

#shader fragment

#version 330 core

in vec3 frag_pos;
in vec3 Normal;
in vec4 tri_color;

layout(location = 0) out vec4 out_col;

uniform vec4 eye_position;
uniform vec4 background_colour;
uniform bool is_perspective_projection;
uniform vec4 light_0_position;


float get_fog_amount(float depth_in) {

   // make this a uniform
   if (! is_perspective_projection) {
      return depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d * d;;
   }

}


void main() {

  float specular_strength = 0.1; // 1.5 is very shiny
  vec4 specular_light_colour = vec4(0.7, 0.7, 0.7, 1.0);

  // a light direction of 0,0,1 is good for fresnelly outlining (well, it used to be)
  vec3 lightdir = normalize(vec3(-2, 1, -2));

  // using the light_0_position give to the shader in a uniform
  // gives us rotating lights
  lightdir = normalize(-light_0_position.xyz);
  float dp = dot(Normal, -lightdir);
  dp = max(dp, 0.0); // no negative dot products for diffuse for now, also, zero is avoided.

  float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);

  float f_1 = 1.0 - m; // because glm::ortho() near and far are reversed?
  float fog_amount = get_fog_amount(gl_FragCoord.z);

  vec4 bg_col = background_colour; // needed?

  vec3 eye_pos_3 =  eye_position.xyz;

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

  vec4 colour_local = tri_color;

  vec4 col_1 = colour_local;  // ambient
  float ambient_strength = 0.3;
  vec4 col_2 = colour_local * dp;
  vec4 col_3 = col_1 * ambient_strength + col_2 + specular;
  vec4 col_4 = mix(col_3, bg_col, fog_amount);

  out_col = col_4;

}
