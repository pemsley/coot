
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec4 background_colour;

out vec3 Normal;
out vec4 line_colour;
out vec4 bg_colour;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   vec4 n = view_rotation * vec4(normal, 1.0);
   Normal = normalize(n).xyz;

   bg_colour = background_colour;

   line_colour = colour;
}

#shader fragment

#version 330 core

in vec3 Normal;
in vec4 line_colour;
in vec4 bg_colour;

void main() {

  vec3 lightdir = normalize(vec3(-2,-1,5));
  float dp = dot(Normal, -lightdir);

  float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);

  float f_1 = 1.0 - m; // because glm::ortho() near and far are reversed?
  float f_2 = 1.0 - abs(f_1 - 0.7)/0.7;

  float flat_frac = 0.4;

  // again

  vec4 bg_col = bg_colour;
  // bg_col = vec4(0.6, 0.2, 0.2, 1.0);

  vec4 col_1 = line_colour;
  vec4 col_2 = line_colour * dp;
  vec4 col_3 = mix(col_1, col_2, flat_frac);
  vec4 col_4 = mix(bg_col, col_3, f_2);
  gl_FragColor = col_4;

}
