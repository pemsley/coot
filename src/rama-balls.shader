

#shader vertex

#version 330 core

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in float displacement;

layout (location = 2) in vec4 colour;
layout (location = 3) in vec3 position;
layout (location = 4) in float size;
layout (location = 5) in float specular_strength;
layout (location = 6) in float shininess;

uniform mat4 mvp;
uniform mat4 view_rotation;

out vec4 colour_transfer;
out vec3 frag_pos_transfer;
out vec3 normal_transfer;
out float specular_strength_transfer;
out float shininess_transfer;

void main() {

   // I think that the configuration/placement of "displacement" is wrong at the moment.
   // displacement is a small number around 0, that gets multiplied by a displacement_scale
   // which is big for big baddies.

   vec4 p3 = vec4(size * (1.0 + displacement) * vertex_position + position, 1.0);
   gl_Position = mvp * p3;

   vec4 n1 = vec4(vertex_position, 1.0);

   frag_pos_transfer = p3.xyz;
   normal_transfer = n1.xyz;
   colour_transfer = colour;
   specular_strength_transfer = specular_strength;
   shininess_transfer = shininess;
}

#shader fragment

#version 330 core

in vec4 colour_transfer;
in vec3 frag_pos_transfer;
in vec3 normal_transfer;
in float specular_strength_transfer;
in float shininess_transfer;

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
uniform bool is_perspective_projection;
uniform bool do_depth_fog;
uniform bool do_diffuse_lighting; // really: all lighting vs ambient only (for demo, really)

layout(location = 0) out vec4 outputColor;

float get_fog_amount(float depth_in) {

   if (! is_perspective_projection) {
      return depth_in;
   } else {
      // needs tweaking
      float d = depth_in;
      float d4 = d * d * d * d;
      return d * d;
   }

}

void main() {

   vec4 bg_col = background_colour;
   vec4 sum_col = vec4(0,0,0,0);

   float shininess = shininess_transfer;
   float specular_strength = specular_strength_transfer;

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {
         vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space;
         float dp = dot(normal_transfer, light_dir);
         // we can't have specular lights where there is no diffuse light
         if (dp <= 0.0)
            specular_strength = 0.0;
         dp = clamp(dp, 0.0, 1.0); // no negative dot products for diffuse

         vec4 lsa = vec4(0.4, 0.4, 0.4, 1.0); // fix these
         vec4 lsd = vec4(0.6, 0.6, 0.6, 1.0);
         vec4 ambient  = colour_transfer * lsa * 0.15;
         vec4 diffuse  = colour_transfer * lsd * dp * 0.9;

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

         // sum_col += ambient + diffuse; // + specular;

         sum_col += ambient + diffuse + specular;

      }
   }
   float fog_amount = 0.0;
   if (do_depth_fog)
      fog_amount = get_fog_amount(gl_FragCoord.z);
   outputColor += mix(sum_col, bg_col, fog_amount);

}

