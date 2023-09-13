
// this used to be obj.shader when it was in hmt.

#shader vertex

#version 330 core
// instanced-meshes-for-ssao.shader

// make it work like 9.ssao_geometry.shader

// the normal calculation is probably wrong

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

layout(location = 3) in vec4 model_rotation_translation_scale_0; // instanced
layout(location = 4) in vec4 model_rotation_translation_scale_1;
layout(location = 5) in vec4 model_rotation_translation_scale_2;
layout(location = 6) in vec4 model_rotation_translation_scale_3;
layout(location = 7) in vec4 colour_instanced;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

// These output values should be in ** screen space **
// so we need the view matrix.
out vec3 frag_pos_transfer;
out vec3 normal_transfer;

void main() {

   mat4 model_rotation_translation_scale = mat4(model_rotation_translation_scale_0,
                                                model_rotation_translation_scale_1,
                                                model_rotation_translation_scale_2,
                                                model_rotation_translation_scale_3);
   mat3 model_rotation = mat3(model_rotation_translation_scale_0.xyz,
                              model_rotation_translation_scale_1.xyz,
                              model_rotation_translation_scale_2.xyz);

   vec3 t_pos = position;
   vec3 n_dir = normal;

   vec4 p4 = vec4(t_pos, 1.0);
   vec4 frag_pos = model_rotation_translation_scale * p4;

   mat4 mvp = projection * view * model;
   gl_Position = mvp * frag_pos;

   // unlike in instanced-meshes-with-shaddws.shader, frag_pos_transfer is a vec3
   //
   frag_pos_transfer = (view * model * frag_pos).xyz;

   mat3 normal_matrix = transpose(inverse(mat3(view * model)));
   normal_transfer = normal_matrix * model_rotation * normal;

}

#shader fragment

#version 330 core
// instanced-meshes-for-ssao.shader

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
