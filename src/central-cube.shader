
#shader vertex

#version 330 core

layout(location = 0) in vec3 position; // origin-based cylinder

uniform mat4 mvp;
uniform mat4 view_rotation; // the quaternion attached to what the mouse has done
uniform vec4 background_colour;
uniform vec4 eye_position;

out vec4 tri_colour;
out vec4 bg_colour;

void main() {

   vec4 p2 = vec4(position, 1.0) * transpose(view_rotation);
   gl_Position = p2;

   tri_colour = vec4(0.8, 0.6, 0.6, 1.0);
   bg_colour = background_colour;
}

#shader fragment

#version 330 core

in vec4 tri_colour;
in vec4 bg_colour;

layout(location = 0) out vec4 out_col;

void main() {


  float f_1 = 1.0 - gl_FragCoord.z; // because glm::ortho() near and far are reversed?
  float f_2 = 1.0 - abs(f_1 - 0.7)/0.7;
  f_2 = f_1; // just testing

  vec4 col_1 = mix(bg_colour, tri_colour, f_2);

  out_col = col_1;

}
