
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
uniform vec4 background_colour;

out vec3 Normal;
out vec4 tri_color;
out vec4 bg_colour;

void main() {

   mat3 model_rotation_matrix = mat3(model_rotation_matrix_0, model_rotation_matrix_1, model_rotation_matrix_2);

   vec3 p2 = position * model_rotation_matrix;
   vec4 p3 = vec4(p2 + model_translation, 1.0);
   gl_Position = mvp * p3;

   vec4 n1 = vec4(normal * model_rotation_matrix, 1.0);
   vec4 n2 = view_rotation * n1;

   Normal = normalize(n2).xyz;

   tri_color = colour;
   bg_colour = background_colour;
}

#shader fragment

#version 330 core

in vec3 Normal;
in vec4 tri_color;
in vec4 bg_colour;

layout(location = 0) out vec4 out_col;

void main() {

  vec3 lightdir = normalize(vec3(-2,-1,5));
  float dp = dot(Normal, -lightdir);

  float f_1 = 1.0 - gl_FragCoord.z; // because glm::ortho() near and far are reversed?
  float f_2 = 1.0 - abs(f_1 - 0.7)/0.7;
  vec4 col_1 = vec4(vec3(f_2), 1.0) * tri_color;
  vec4 col_2 = col_1 * dp;

  float flat_frac = 0.2;
  vec4 c_1 = col_2 * (1.0 - flat_frac) + col_1 * flat_frac;
  vec4 c_2 = mix(bg_colour, c_1, f_1);

  out_col = c_2;


}
