
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

uniform mat4 mvp;

out vec3 Normal;

void main() {

   gl_Position = mvp * vec4(position, 1.0);
   Normal = normal;

}

#shader fragment

#version 330 core

in vec3 Normal;

void main() {

  vec4 line_color = vec4(0.3, 0.5, 0.7, 1.0);
  vec4 background_color = vec4(0.0f, 0.0f, 0.0f, 0.0f);

  vec3 lightdir = normalize(vec3(-2,-1,4));
  vec4 color = vec4(0.4, 0.6, 0.7, 1.0);
  float dp = dot(Normal, -lightdir);
  gl_FragColor = color * dp;
  // gl_FragColor = color; 
  
  // gl_FragColor = vec4(vec3(gl_FragCoord.z), 1.0);
  // float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);
  // gl_FragColor = mix(line_color, background_color, m);

  // depth cue
  //
  float m  = clamp(gl_FragCoord.x, 0.0f, 1.0f);
  // gl_FragCoord.z seems to have a range of 0-1, unlike gl_FragCoord.x, which seems to be -100,100 or so.
  m  = 0.9 * gl_FragCoord.z;
  gl_FragColor = mix(line_color, background_color, m);
  gl_FragColor = line_color;
  gl_FragColor = vec4(vec3(1.0 - gl_FragCoord.z), 1.0) * dp;

  

}
