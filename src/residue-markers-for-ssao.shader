
#shader vertex

#version 330 core

// layout(location = 0) in vec3 position; // per vertex for the quad/hex
// layout(location = 1) in vec3 normal;   // per vertex for the quad/hex
// layout(location = 2) in vec4 colour;   // per vertex for the quad/hex
// layout(location = 3) in vec2 texCoord;
// layout(location = 4) in vec3 instance_translation;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec3 bitangent;
layout(location = 4) in vec4 colour;
layout(location = 5) in vec2 texCoord;
layout(location = 6) in vec3 instance_translation;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 frag_pos_transfer;
out vec3 normal_transfer;

// uniform mat4 mvp;
// nuniform mat4 view_rotation;
// uniform float canvas_scale; // 0.8 for happy faces, 0.2 for anchored/fixed atoms
// nout vec2 texCoord_transfer;
// out vec3 frag_pos_transfer;

void main() {

   mat4 mvp = projection * view * model;
   gl_Position = mvp * vec4(position, 1.0);

   mat3 normal_matrix = mat3(view);
   mat3 transpose_normal_matrix = transpose(normal_matrix);
   mat3 view_3_3 = mat3(view);

   normal_transfer = normal * view_3_3 * transpose_normal_matrix;
   vec4 frag_pos_transfer_v4 = view * model * vec4(position, 1.0);
   frag_pos_transfer = frag_pos_transfer_v4.xyz;
}


#shader fragment

#version 330 core

layout (location = 0) out vec3 gPosition;
layout (location = 1) out vec3 gNormal;

in vec3 frag_pos_transfer;
in vec3 normal_transfer;

void main() {

   vec3 normal = normal_transfer;

   //    vec3 normal_view_rotated = (vec4(normal, 1.0) * view_rotation).xyz;
   // normal_view_rotated = normal; // I don't understand how this works.

   gPosition = frag_pos_transfer;
   //    gNormal = normalize(normal_transfer); // do we need to normalize?
   gNormal = normal_transfer;

}
