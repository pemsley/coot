
#shader vertex

#version 120

attribute vec3 position;

void main() {

   vec4 mvpp = gl_ModelViewProjectionMatrix * vec4(position, 1.0);
   gl_Position = mvpp;

   vec3 vVertex = vec3(gl_ModelViewProjectionMatrix * vec4(position, 1.0));


}

#shader fragment

#version 120

void main() {

  vec4 final_color = vec4(0.3, 0.5, 0.7, 1.0);
  vec4 background_color = vec4(0.0f, 0.0f, 0.0f, 0.0f);
  // gl_FragColor = final_color;

  // depth cue
  // gl_FragColor = vec4(vec3(gl_FragCoord.z), 1.0);
  float m  = clamp(gl_FragCoord.z, 0.0f, 1.0f);

  gl_FragColor = mix(final_color, background_color, m);

}
