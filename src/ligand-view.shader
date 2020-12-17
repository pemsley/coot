
#shader vertex

#version 330 core

layout(location = 0) in vec2 position;

uniform float aspect_ratio;

void main() {

   vec2 scaled_pos = 0.05 * position;
   scaled_pos.x /= aspect_ratio;
   vec2 offset_pos = scaled_pos + vec2(-0.6, -0.6);
   gl_Position = vec4(offset_pos, 0.999, 1.0);

}

#shader fragment

#version 330 core

out vec4 out_col;

void main() {

   out_col = vec4(0.7, 0.7, 0.7, 1.0);
}
