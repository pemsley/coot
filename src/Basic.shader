
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

uniform mat4 mvp;

out vec3 Normal;
out vec4 line_color;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   vec4 n = mvp * vec4(normal, 1.0);
   Normal = n.xyz;

   line_color = vec4(0.4, 0.5, 1.0, 0.5); // this can be a uniform
}

#shader fragment

#version 330 core

in vec3 Normal;
in vec4 line_color;

void main() {

  vec4 background_color = vec4(0.0f, 0.0f, 0.0f, 0.0f);

  vec3 lightdir = normalize(vec3(-2,-1,2));
  float dp = dot(Normal, -lightdir);

  float m  = clamp(gl_FragCoord.x, 0.0f, 1.0f);
  // gl_FragCoord.z seems to have a range of 0-1, unlike gl_FragCoord.x, which seems to be -100,100 or so.
  gl_FragColor = vec4(vec3(1.0 - gl_FragCoord.z), 1.0) * line_color * dp;

  

}
