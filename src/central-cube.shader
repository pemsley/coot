
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;

uniform mat4 mvp;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

}

#shader fragment

#version 330 core

void main() {

  vec4 line_color = vec4(0.8, 0.7, 0.5, 1.0);
  vec4 background_color = vec4(0.0f, 0.0f, 0.0f, 0.0f);

  // depth cue
  gl_FragColor = vec4(vec3(gl_FragCoord.z), 1.0);
  float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);

  // gl_FragColor = mix(line_color, background_color, m);
  gl_FragColor = mix(line_color, background_color, m);

}

