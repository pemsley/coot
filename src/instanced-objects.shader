
// use the #shader "directive" to separate the shaders on parsing
// -----------------------------------------
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour; // instanced. For instanced objects, the colour is set
                                     // in the generic vertex (i.e. this), not the material.

layout(location = 3) in vec4 model_rotation_translation_scale_0; // instanced
layout(location = 4) in vec4 model_rotation_translation_scale_1;
layout(location = 5) in vec4 model_rotation_translation_scale_2;
layout(location = 6) in vec4 model_rotation_translation_scale_3;


uniform mat4 mvp;
uniform mat4 view_rotation;
uniform float time;

out vec4 colour_transfer;
out vec3 normal_transfer;
out vec4 frag_pos_transfer;

void main() {
   mat4 model_rotation_translation_scale = mat4(model_rotation_translation_scale_0,
                                                model_rotation_translation_scale_1,
                                                model_rotation_translation_scale_2,
                                                model_rotation_translation_scale_3);
   mat3 model_rotation = mat3(model_rotation_translation_scale_0.xyz,
                              model_rotation_translation_scale_1.xyz,
                              model_rotation_translation_scale_2.xyz);

   // vec3 t_pos = position + vec3(0,0, 0.5 * sin(0.01 * time));
   vec3 t_pos = position + vec3(0,0, 0.2 * sin(0.005 * time + 0.1 * gl_InstanceID));
   vec4 p4 = vec4(t_pos, 1.0);
   vec4 frag_pos = model_rotation_translation_scale * p4;

   gl_Position = mvp * frag_pos;

   normal_transfer = model_rotation * normal;
   colour_transfer = colour;
   frag_pos_transfer = frag_pos;
}


// -----------------------------------------
#shader fragment

#version 330 core

struct LightSource {
    bool is_on;
    bool directional;
    vec4 position;
    vec3 direction_in_molecule_coordinates_space;
    vec4 ambient;
    vec4 diffuse;
    vec4 specular;
    vec4 halfVector;
    vec3 spotDirection;
    float spotExponent;
    float spotCutoff;
    float spotCosCutoff;
    float constantAttenuation;
    float linearAttenuation;
    float quadraticAttenuation;
};
struct Material {
   float shininess;
   float specular_strength;
   vec4 specular;
};
uniform LightSource light_sources[2];
uniform vec3 eye_position;
uniform Material material;
uniform bool is_perspective_projection;
uniform vec4 background_colour;

in vec4 frag_pos_transfer;
in vec4 colour_transfer; // for instanced objects, the colour is set in the generic vertex, not the material
in vec3 normal_transfer;

out vec4 outputColor;

float get_fog_amount(float depth_in) {

   if (! is_perspective_projection) {
      return depth_in * depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d4;
   }
}

void main() {

   vec4 ct = colour_transfer;
   // ct = vec4(1,1,1,1);
   outputColor = vec4(0,0,0,0);

   float fog_amount = get_fog_amount(gl_FragCoord.z);

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {

         // ambient
         vec4 ambient = ct * light_sources[i].ambient * 0.2; // we are not using material here

         // diffuse
         vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space.xyz;
         vec3 norm_2 = normalize(normal_transfer); // not needed?

         // light_dir = vec3(0,0,1);

         float dp = dot(norm_2, light_dir);
         dp = max(dp, 0.0);
         vec4 diffuse = ct * light_sources[i].diffuse * dp * 0.8;

         // specular

         float shininess = material.shininess;
         vec3 eye_pos = eye_position;
         vec3 view_dir = normalize(eye_pos - frag_pos_transfer.xyz);
         vec3 light_dir_v3 = light_dir.xyz;
         vec3 reflect_dir = reflect(light_dir_v3, norm_2);
         reflect_dir = normalize(reflect_dir); // belt and braces
         float dp_view_reflect = dot(view_dir, -reflect_dir);
         dp_view_reflect = max(dp_view_reflect, 0.0);
         dp_view_reflect = min(dp_view_reflect, 1.0);

         float spec = material.specular_strength * pow(dp_view_reflect, shininess);
         // spec = 0;
         vec4 specular = spec * light_sources[i].specular;

         // final
         outputColor += ambient + diffuse + specular;

         // if (dp_view_reflect > 0.98) outputColor = vec4(1,0,0,1);
         if (dp > 1.0) outputColor = vec4(1,1,0,1);
         if (spec < 0.0) outputColor = vec4(0,1,1,1);
         if (dp_view_reflect < 0.0) outputColor = vec4(0,1,1,1);

         vec4 o = vec4(max(norm_2, 0.0), 1.0);
         o = vec4(0.95 * (light_dir.x),
                  0.95 * (light_dir.y),
                  0.95 * (light_dir.z),
                  1.0);
         // outputColor = o;

      }
   }

   outputColor = mix(outputColor, background_colour, fog_amount);


}
