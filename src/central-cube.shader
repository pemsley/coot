
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform vec4 line_colour;
uniform vec4 background_colour;

out vec4 tri_colour;
out vec4 bg_colour;

void main() {

   vec4 p2 = vec4(position, 1.0) * transpose(mvp);
   gl_Position = p2;

   tri_colour = vec4(0.5, 0.4, 0.4, 1.0);
   tri_colour = line_colour;
   bg_colour = background_colour;
}

#shader fragment

#version 330 core

in vec4 tri_colour;
in vec4 bg_colour;

layout(location = 0) out vec4 out_col;

void main() {

  float f_1 = 1.0 - gl_FragCoord.z; // because glm::ortho() near and far are reversed?
  vec4 col_1 = mix(bg_colour, tri_colour, f_1);

  out_col = col_1;

}
