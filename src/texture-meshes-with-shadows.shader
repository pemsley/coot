
// this used to be obj.shader when it was in hmt.

#shader vertex

// texture-meshes-with-shadows.shader
#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 tangent;
layout(location = 3) in vec3 bitangent;
layout(location = 4) in vec4 colour;
layout(location = 5) in vec2 texCoord;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform mat4 light_space_mvp;

out vec3 frag_pos_transfer;
out vec4 colour_transfer;
out vec3 normal_transfer;
out vec2 texCoord_transfer;
out vec4 frag_pos_light_space_transfer;
out mat3 TBN_transfer;

void main() {

   gl_Position = mvp * vec4(position, 1.0);

   mat3 normal_matrix = mat3(view_rotation);
   mat3 transpose_normal_matrix = transpose(normal_matrix);

   colour_transfer = colour;
   normal_transfer = normal * transpose_normal_matrix;
   texCoord_transfer = texCoord;
   frag_pos_transfer = position;
   frag_pos_light_space_transfer = light_space_mvp * vec4(position, 1.0); // 15, like gl_Position

   vec3 T = normal_matrix * tangent;
   vec3 B = normal_matrix * bitangent;
   vec3 N = normal_matrix * normal;

   TBN_transfer = mat3(T,B,N);

}

#shader fragment

// texture-meshes-with-shadows.shader
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
uniform sampler2D specular_map;
uniform sampler2D normal_map;
uniform sampler2D shadow_map;

// spacenerd3_14: shininess = 100/(500 * roughness + 0.01)

uniform bool specular_map_included;
uniform bool normal_map_included;
uniform bool reversed_normals;

uniform vec3 eye_position;
uniform vec3 eye_position_in_molecule_coordinates_space;
uniform mat4 view_rotation;
uniform bool do_depth_fog;
uniform float shadow_strength;
uniform int shadow_softness;

in vec3 frag_pos_transfer;
in vec4 colour_transfer;
in vec3 normal_transfer;
in vec2 texCoord_transfer;
in vec4 frag_pos_light_space_transfer;
in mat3 TBN_transfer;

out vec4 output_colour;

// useful for debugging and demonstration, this flag mean show *just* shadows
uniform bool show_shadows;


float calc_shadow(float dp_light_to_fragment) {

   // calc_shadow() returns 0.0 for full shadow

   vec3 pos = frag_pos_light_space_transfer.xyz * 0.5 + 0.5;
   // points outside the shadow map frustrum can have pos.z
   // values greater than 1.0, treat them as if they
   // are not shadowed, i.e.
   // (depth + bias) > pos.z
   //
   if (pos.z > 1.0)
      pos.z = 1.0;
   // if the dp is small then
   float bias = (1.0 + dp_light_to_fragment) * 0.01;
   // if (bias < 0.001) bias = 0.001;
   // bias = 0.0006;
   // setting bias to 0, makes the shadow box (ortho projection box)
   // easy to see
   bias = 0.005;
   bias = 0.002;
   float shadow = 0.0;
   vec2 texelSize = 1.0 / textureSize(shadow_map, 0);
   int pixel_count = 0;
   int ss = shadow_softness; // unsigned int to int conversion
   for (int ix = -ss; ix<=ss; ix++) {
      for (int iy = -ss; iy<=ss; iy++) {
         float depth = texture(shadow_map, pos.xy + vec2(ix,iy) * texelSize).r;
         shadow += ((depth + bias) < pos.z) ? 0.0 : 1.0;
         pixel_count++;
      }
   }
   // return ((depth + bias) < pos.z) ? 0.0 : 1.0;
   return (1.0-shadow_strength) + (shadow * shadow_strength) / pixel_count;

   // try again (simple)
   // float depth = texture(shadow_map, pos.xy).r;
   // return ((depth + bias) < pos.z) ? 0.0 : 1.0;
}

void main() {

   vec3 normal = normal_transfer;
   vec3 normal_view_rotated = (vec4(normal, 1.0) * view_rotation).xyz;
   normal_view_rotated = normal; // I don't understand how this works.
   output_colour = vec4(0,0,0,0);

   bool use_normal_map = false;
   bool invert_normal_map = true;

   float scale_factor_n_lights = 0.85;

   vec4 sampled          = texture(base_texture, texCoord_transfer);
   vec4 sampled_specular = texture(specular_map, texCoord_transfer);

   // if (show_shadows) { // show *just* shadows
   if (true) { // show *just* shadows

      mat4 ivr = transpose(view_rotation);
      vec3 light_dir = (vec4(light_sources[0].direction_in_molecule_coordinates_space, 1.0) * ivr).xyz;
      float dp_raw = dot(normal_transfer, light_dir);
      float shadow = calc_shadow(dp_raw); // calc_shadow() returns 0.0 for full shadow
      output_colour = vec4(vec3(shadow), 1.0);

   } else {

      for (int i=0; i<2; i++) {

         if (light_sources[i].is_on) {

            float shininess = 8;
            float specular_strength = 0.1;
            if (specular_map_included) {
               specular_strength = 3.0 * sampled_specular.r;
               shininess = 16;
            }

            vec4 ambient = 0.3 * sampled * light_sources[i].ambient;

            mat4 ivr = transpose(view_rotation);
            vec3 light_dir = (vec4(light_sources[i].direction_in_molecule_coordinates_space, 1.0) * ivr).xyz;

            float dp_raw = dot(normal_transfer, light_dir);
            float dp = dp_raw;
            if (dp < 0)
               specular_strength = 0.0; // no shiny insides
            dp = clamp(dp_raw, 0.0, 1.0);

            // by hand "Material" information.
            // in future we could either pass a Material (which would have to be be added
            // in the TextureMesh constructor)
            // or
            // get "Material" information (shininess) from other textures (sampler2Ds)

            bool flag_raised = false;

            if (normal_map_included) {

               vec4 nm_sampled = texture(normal_map, texCoord_transfer);
               vec3 nm_sampled_xyz = nm_sampled.xyz;
               vec3 normal_map_normal = 2.0 * nm_sampled_xyz - 1.0;

               vec3 normal_rotated = TBN_transfer * normal_map_normal;

               dp = dot(normal_rotated, light_dir);
               dp = clamp(dp, 0.0, 1.0);
               // if (reversed_normals) flag_rasised = true;
               // float dp_nm = dot(normal, normal_map_normal_rotated); // testing
               // ambient = 0.65 * nm_sampled;
            } else {
               // ambient = vec4(0,0,1,1); fun
            }

            vec3 view_dir = normalize(eye_position_in_molecule_coordinates_space - frag_pos_transfer);
            vec3 reflect_dir = reflect(-light_dir, normal);
            reflect_dir = normalize(reflect_dir); // belt and braces
            float dp_view_reflect = dot(view_dir, reflect_dir);
            dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
            float spec = pow(dp_view_reflect, shininess);
            vec4 specular = specular_strength * spec * vec4(0.8, 0.8, 0.8, 1.0);

            float shadow = calc_shadow(dp_raw); // calc_shadow() returns 0.0 for full shadow

            if (false) {
               float ambient_sum = ambient.r + ambient.g + ambient.b;
               float a = 0.3333 * ambient_sum;
               ambient.rgb = vec3(a,a,a);
               float sampled_sum = sampled.r + sampled.g + sampled.b;
               a = 0.3333 * sampled_sum;
               sampled.rgb = vec3(a,a,a);
            }

            output_colour += scale_factor_n_lights * ambient;
            output_colour += scale_factor_n_lights * shadow * sampled * dp;
            output_colour += scale_factor_n_lights * shadow * specular;

            // hack/test
            // output_colour = vec4(vec3(shadow), 1.0);

            if (flag_raised)
               output_colour = vec4(1,0,0,1);
         }
      }
   }

}
