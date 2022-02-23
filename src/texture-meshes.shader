
// this used to be obj.shader when it was in hmt.

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec3 bitangent;
layout(location = 4) in vec4 colour;
layout(location = 5) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;

out vec3 frag_pos_transfer;
out vec4 colour_transfer;
out vec3 normal_transfer;
out vec2 texCoord_transfer;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   // gl_Position = vec4(2 * position - vec3(1.0, 1.0, 0.0), 1.0);

   mat3 normal_matrix = mat3(view_rotation);
   mat3 transpose_normal_matrix = transpose(normal_matrix);

   colour_transfer = colour;
   normal_transfer = normal * transpose_normal_matrix;
   texCoord_transfer = texCoord;
   frag_pos_transfer = position;

}

#shader fragment

#version 330 core

struct LightSource {
   bool is_on;
   bool directional;
   vec3 position;
   vec3 direction;
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

uniform sampler2D base_texture;  // diffuse_map

uniform vec3 eye_position;
uniform vec3 eye_position_in_molecule_coordinates_space;
uniform mat4 view_rotation;
uniform bool do_depth_fog;

in vec3 frag_pos_transfer;
in vec4 colour_transfer;
in vec3 normal_transfer;
in vec2 texCoord_transfer;

out vec4 output_colour;

void main() {

   vec3 normal = normal_transfer;
   vec3 normal_view_rotated = (vec4(normal, 1.0) * view_rotation).xyz;
   normal_view_rotated = normal; // I don't understand how this works.
   output_colour = vec4(0,0,0,0);

   bool use_normal_map = false;
   bool invert_normal_map = true;

   float scale_factor_n_lights = 0.95;

   vec4 sampled = texture(base_texture, texCoord_transfer);
   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {

         float shininess = 128;
         float specular_strength = 0.5;

         vec4 ambient  = 0.3 * sampled * light_sources[i].ambient;
         // vec3 light_dir = (vec4(light_sources[i].direction, 1.0) * view_rotation).xyz;
         mat4 ivr = transpose(view_rotation);
         vec3 light_dir = (vec4(light_sources[i].direction_in_molecule_coordinates_space, 1.0) * ivr).xyz;
         // light_dir = light_sources[i].direction;
         float dp = dot(normal_transfer, light_dir);
         if (dp < 0)
            specular_strength = 0.0; // no shiny insides
         dp = clamp(dp, 0.0, 1.0);

         // by hand "Material" information.
         // in future we could either pass a Material (which would have to be be added
         // in the TextureMesh constructor)
         // or
         // get "Material" information (shininess) from other textures (sampler2Ds)

         vec3 view_dir = normalize(eye_position_in_molecule_coordinates_space - frag_pos_transfer);
         vec3 reflect_dir = reflect(-light_dir, normal);
         reflect_dir = normalize(reflect_dir); // belt and braces
         float dp_view_reflect = dot(view_dir, reflect_dir);
         dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
         float spec = pow(dp_view_reflect, shininess);
         vec4 specular = specular_strength * spec * vec4(0.8, 0.8, 0.8, 1.0);

         output_colour += scale_factor_n_lights * ambient;
         output_colour += scale_factor_n_lights * sampled * dp;
         output_colour += scale_factor_n_lights * specular;
      }
   }

   // output_colour = vec4(0.95, 0.95, 0.95, 1.0);

}
