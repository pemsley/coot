
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform mat4 view_rotation;

out vec3 Normal;
out vec4 line_color;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   vec4 n = view_rotation * vec4(normal, 1.0);
   Normal = normalize(n).xyz;

   line_color = colour;
}

#shader fragment

#version 330 core

in vec3 Normal;
in vec4 line_color;

void main() {

  vec4 background_color = vec4(0.0f, 0.0f, 0.0f, 0.0f);

  vec3 lightdir = normalize(vec3(-2,-1,3));
  float dp = dot(Normal, -lightdir);

  float m  = clamp(gl_FragCoord.x, 0.0f, 1.0f);
  // gl_FragCoord.z seems to have a range of 0-1, unlike gl_FragCoord.x, which seems to be -100,100 or so.
  vec4 col1 = vec4(vec3(1.0 - gl_FragCoord.z), 1.0) * line_color * dp;
  vec4 col2 = vec4(vec3(1.0 - gl_FragCoord.z), 1.0) * line_color;
  float flat_frac = 0.5;
  gl_FragColor = col1 * (1.0 - flat_frac) + col2 * flat_frac;


  // interesting effects of this kind: gl_FragColor += vec4(0.2, 0.2, 0.2, 1.0);
  // I want to make the back face be grey, not black.

}
