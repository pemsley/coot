
// use the #shader "directive" to separate the shaders on parsing
// -----------------------------------------
#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;

out vec4 colour_transfer;
out vec3 normal_transfer;
out vec3 frag_pos_transfer;

void main() {
   vec4 n = vec4(normal, 1.0);
   gl_Position = mvp * vec4(position, 1.0);

   normal_transfer = normalize(normal); // probably the normalize() here is not needed
   colour_transfer = colour;
   frag_pos_transfer = position;
}


// -----------------------------------------
#shader fragment

#version 330 core

struct LightSource {
   bool is_on;
   bool directional;
   vec3 position;
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

uniform LightSource light_sources[2];
uniform vec3 eye_position;
uniform vec4 background_colour;
uniform bool do_depth_fog;

in vec4 colour_transfer;
in vec3 normal_transfer;
in vec3 frag_pos_transfer;

out vec4 outputColor;

struct Material {
   float shininess;
   float specular_strength;
   vec4 specular;
};
uniform Material material;

float get_fog_amount(float depth_in) {

   const bool is_perspective_projection = true;

   if (! is_perspective_projection) {
      return depth_in;
   } else {
      float d = depth_in;
      return 0.9 * d * d * d * d;
   }
}

void main() {

   vec4 bg_col = background_colour;
   vec4 sum_col = vec4(0,0,0,0);

   float specular_strength = material.specular_strength;
   float shininess = material.shininess;

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {
         vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space;
         float dp = dot(normal_transfer, light_dir);
         // we can't have specular lights where there is no diffuse light
         if (dp <= 0.0)
            specular_strength = 0.0;
         dp = clamp(dp, 0.1, 1.0); // no negative dot products for diffuse

         // why don't I use the light source?
         vec4 lsa = vec4(0.4, 0.4, 0.4, 1.0);
         vec4 lsd = vec4(0.6, 0.6, 0.6, 1.0);
         vec4 ambient = colour_transfer * lsa * 0.1;
         vec4 diffuse = colour_transfer * lsd * dp * 0.8;

         // specular
         vec3 eye_pos = eye_position;
         vec3 norm_2 = normalize(normal_transfer); // not needed, I think
         vec3 view_dir = normalize(eye_pos - frag_pos_transfer);
         vec3 reflect_dir = reflect(-light_dir, norm_2);
         reflect_dir = normalize(reflect_dir); // belt and braces
         float dp_view_reflect = dot(view_dir, reflect_dir);
         dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
         float spec = pow(dp_view_reflect, shininess);
         vec4 specular = specular_strength * spec * light_sources[i].specular;

         sum_col += ambient + diffuse + specular;

         // vec3 light_to_eye = normalize(eye_pos - 100.0 * light_dir);
         // sum_col = vec4(0.5 * view_dir + vec3(0.5,0.5,0.5), 1.0);
         // sum_col = vec4(0.5 * light_to_eye + vec3(0.5,0.5,0.5), 1.0);
      }
   }

   float fog_amount = 0.0;
   if (do_depth_fog)
      fog_amount = get_fog_amount(gl_FragCoord.z);
   outputColor += mix(sum_col, bg_col, fog_amount);

}

