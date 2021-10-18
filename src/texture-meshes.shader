
// this used to be obj.shader when it was in hmt.

#shader vertex

#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;
layout(location = 3) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;
// uniform sampler2D roughness_map; // i.e. displacement map

out vec3 frag_pos_transfer;
out vec4 colour_transfer;
out vec3 normal_transfer;
out vec2 texCoord_transfer;

void main() {

   gl_Position = mvp * vec4(position, 1.0);
   // bool do_displacement_mapping = false; // meh
   // if (do_displacement_mapping) {
      // vec4 rm_sampled = texture(roughness_map, texCoord);
      // vec4 displacement = 0.5 * rm_sampled.r * vec4(normal, 0.0) * view_rotation;
   // gl_Position += displacement;
   // }
   colour_transfer = colour;
   normal_transfer = normal;
   texCoord_transfer = texCoord;
   frag_pos_transfer = position;

}

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

uniform sampler2D base_texture;  // diffuse_map
uniform sampler2D normal_map;

// uniform sampler2D occlusion_map; // ambient occlusion
// uniform sampler2D roughness_map; // i.e. displacement map
// uniform sampler2D specular_map;
uniform vec3 eye_position;
uniform mat4 view_rotation;
uniform bool do_depth_fog;

in vec3 frag_pos_transfer;
in vec3 normal_transfer;
in vec4 colour_transfer;
in vec2 texCoord_transfer;

out vec4 output_colour;

void main() {

   vec3 normal = normal_transfer;
   vec3 normal_view_rotated = (vec4(normal, 1.0) * view_rotation).xyz;
   normal_view_rotated = normal; // I don't understand how this works.
   output_colour = vec4(0,0,0,0);

   bool use_normal_map = false;
   bool invert_normal_map = true;

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {

         if (use_normal_map) {

            vec4 ambient  = 0.1 * colour_transfer * light_sources[i].ambient;

            // vec3 lightdir = normalize(vec3(-2,-5,0)); // different monkey
            vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space;

            vec4 nm_sampled = texture(normal_map, texCoord_transfer); // why is this a vec4? is this right?
            vec4 normal_map_normal_tmp = 2.0 * nm_sampled;
            vec4 normal_map_normal = vec4(normal_map_normal_tmp.x -1.0,
                                          normal_map_normal_tmp.y -1.0,
                                          normal_map_normal_tmp.z -1.0, 1.0);

            // if (invert_normal_map)
            // normal_map_normal = vec4(-normal_map_normal.xyz, 1.0);

            vec3 normal_map_normal_rotated = normalize(normal_map_normal * view_rotation).xyz;
            float normal_factor = 0.4;
            vec3 normal_combined = normal_factor * normal_view_rotated + normal_factor * normal_map_normal_rotated;
            float dp = dot(normal_combined, light_dir);
            dp = clamp(dp, 0.0, 1.0);
            float dp_nm = dot(normal, normal_map_normal_rotated); // testing

            // specular
            float shininess = 128;
            float specular_strength = 0.04;

            if (dp == 0.0) specular_strength = 0.0;
            vec3 eye_pos = eye_position;
            vec3 view_dir = normalize(eye_pos - frag_pos_transfer);
            vec3 reflect_dir = reflect(-light_dir, normal_combined);
            reflect_dir = normalize(reflect_dir); // belt and braces
            float dp_view_reflect = dot(view_dir, reflect_dir);
            dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
            float spec = pow(dp_view_reflect, shininess);
            vec4 specular = specular_strength * spec * dp * dp_view_reflect * light_sources[i].specular;

            vec4 sampled = texture(base_texture, texCoord_transfer);

            // output_colour += ambient + sampled * dp * 2.0 + specular;

            output_colour += ambient + sampled * dp_nm * 2.0 + specular;

            // output_colour = vec4(1,1,0,1);
         } else {

            // smooth surface

            vec4 ambient  = 0.01 * colour_transfer * light_sources[i].ambient;

            vec3 light_dir = normalize(vec3(1, 0, 4));
            light_dir = (vec4(light_dir, 1.0) * view_rotation).xyz;
            float dp = dot(normal_transfer, light_dir);
            dp = 0.5 * clamp(dp, 0.0, 1.0);

            // specular
            float shininess = 64;
            float specular_strength = 0.02;
            vec3 eye_pos = eye_position;
            vec3 view_dir = normalize(eye_pos - frag_pos_transfer);
            vec3 reflect_dir = reflect(-light_dir, normal);
            reflect_dir = normalize(reflect_dir); // belt and braces
            float dp_view_reflect = dot(view_dir, reflect_dir);
            dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
            float spec = pow(dp_view_reflect, shininess);
            vec4 specular = specular_strength * spec * vec4(0.8, 0.8, 0.8, 1.0);

            vec4 sampled = texture(base_texture, texCoord_transfer);
            output_colour += ambient + sampled * dp + specular;
            // output_colour = vec4(0,1,0,1);
         }
      }
   }

}
