/*
 * src/instanced-meshes-with-shadows.shader
 *
 * Copyright 2023 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */



// -----------------------------------------
#shader vertex

#version 330 core

// This is instanced-meshes-with-shadows.shader

// keep these because they are s_generic_vertex:
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

layout(location = 3) in vec4 model_rotation_translation_scale_0; // instanced
layout(location = 4) in vec4 model_rotation_translation_scale_1;
layout(location = 5) in vec4 model_rotation_translation_scale_2;
layout(location = 6) in vec4 model_rotation_translation_scale_3;
layout(location = 7) in vec4 colour_instanced;

uniform mat4 mvp;
uniform mat4 light_space_mvp;
uniform mat4 view_rotation;
uniform float time;
uniform bool transferred_colour_is_instanced;
uniform bool do_pulse;
uniform bool do_rotate_z;
uniform float pulsing_amplitude; // = 0.25;  // make these uniforms?
uniform float pulsing_frequency; // = 5.0;
uniform float pulsing_phase_distribution; // = 0.2;
uniform float z_rotation_angle;

out vec4 colour_transfer;
out vec3 normal_transfer;
out vec4 frag_pos_transfer;
out vec4 frag_pos_light_space_transfer;

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

   gl_Position = mvp * frag_pos;

   normal_transfer = model_rotation * n_dir;
   colour_transfer = colour;
   if (transferred_colour_is_instanced) colour_transfer = colour_instanced;
   frag_pos_transfer = frag_pos;
   frag_pos_light_space_transfer = light_space_mvp * vec4(frag_pos_transfer);

   // hack test
   // colour_transfer = vec4(model_rotation_translation_scale_0.xyz, 1.0);

}


// -----------------------------------------
#shader fragment

// This is instanced-meshes-with-shadows.shader

#version 330 core

in vec4 colour_transfer;
in vec3 normal_transfer;
in vec4 frag_pos_transfer;
in vec4 frag_pos_light_space_transfer;

out vec4 outputColor;

struct LightSource {
   bool is_on;
   bool directional;
   vec4 position; // is this actually sent?
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
   vec4 ambient;
   vec4 diffuse;
   vec4 specular;
};
uniform Material material;

uniform LightSource light_sources[2];
uniform mat4 view_rotation;
uniform vec3 eye_position;
uniform vec3 eye_position_in_molecule_coordinates_space;
uniform vec4 background_colour;
uniform bool is_perspective_projection;
uniform bool do_depth_fog;
uniform float opacity; // was map_opacity
uniform bool do_specular;
uniform float shadow_strength;
uniform int shadow_softness;


uniform sampler2D shadow_map;
// useful for debugging and demonstration, this flag mean show *just* shadows
uniform bool show_shadows;

float get_fog_amount(float depth_in) {

   if (do_depth_fog) {
      if (! is_perspective_projection) {
         float d = depth_in;
         return d * d;
      } else {
         // needs tweaking
         float d = depth_in;
         float d4 = d * d * d * d;
         // d4 = 0.25 * d; // 20220202-PE crow m,v,p matrices
         d4 = d * d * d * d;
         return 0.75 * d4;
      }
   } else {
      return 0.0;
   }
}

float calc_shadow(float dp_light_to_fragment) {

   // uses shadow_strength

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
   bias = 0.003;
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

   return (1.0-shadow_strength) + (shadow * shadow_strength) / pixel_count;
   // return texture(shadow_map, pos.xy).r;

}

void main() {

   vec4 ct = colour_transfer;
   // ct = vec4(1,1,1,1);
   vec4 running_col = vec4(0,0,0,0);

   float fog_amount = get_fog_amount(gl_FragCoord.z);

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {

         // ambient
         vec4 ambient = ct * light_sources[i].ambient * 0.2; // we are not using material here

         // diffuse
         vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space.xyz;
         vec3 norm_2 = normalize(normal_transfer); // not needed?

         // light_dir = vec3(0,0,1);

         float dp_raw = dot(norm_2, light_dir);
         float dp = max(dp_raw, 0.0);
         vec4 diffuse = ct * light_sources[i].diffuse * dp * 0.8;

         // specular

         float shininess = material.shininess;
         float specular_strength = material.specular_strength;
         if (dp_raw < 0.0)
            specular_strength = 0.0; // no shiny interiors

         float shadow_diffuse = calc_shadow(dp); // calc_shadow() returns 0.0 for full shadow
         float shadow_specular = shadow_diffuse;

         // if (shadow_specular < 0.95)
         //    shadow_specular = 0;

         //
         vec3 eye_pos_in_view = vec3(vec4(eye_position, 1.0) * view_rotation);
         vec3 view_dir = normalize(eye_pos_in_view - frag_pos_transfer.xyz);
         vec3 light_dir_v3 = light_dir.xyz;
         vec3 reflect_dir = reflect(light_dir_v3, norm_2);
         reflect_dir = normalize(reflect_dir); // belt and braces
         float dp_view_reflect = dot(view_dir, -reflect_dir);
         dp_view_reflect = max(dp_view_reflect, 0.0);
         dp_view_reflect = min(dp_view_reflect, 1.0);

         float spec = specular_strength * pow(dp_view_reflect, shininess);

         vec4 specular = 2.0 * spec * light_sources[i].specular;

         // final

         // running_col += ambient + diffuse * shadow_diffuse + specular * shadow_specular * specular_strength;
         vec4 test_col = vec4(0.4, 0.0, 0.4, 1.0);

         vec4 addition_col = vec4(0,0,0,0);

         if (shadow_diffuse > 1.99)
            addition_col = test_col;
         else
            addition_col += ambient + diffuse * shadow_diffuse + specular * shadow_specular * specular_strength;

         // addition_col = vec4(vec3(shadow_diffuse), 1.0);
         // addition_col = vec4(0.3, 0.3, 0.3, 1.0);

         running_col += addition_col;

         // if (dp_view_reflect > 0.98) outputColor = vec4(1,0,0,1);
         if (dp > 1.0) running_col = vec4(1,1,0,1);
         if (spec < 0.0) running_col = vec4(0,1,1,1);
         if (dp_view_reflect < 0.0) running_col = vec4(0,1,1,1);

         vec4 o = vec4(max(norm_2, 0.0), 1.0);
         o = vec4(0.95 * (light_dir.x),
                  0.95 * (light_dir.y),
                  0.95 * (light_dir.z),
                  1.0);
         // outputColor = o; // testing
      }
   }

   outputColor = mix(running_col, background_colour, fog_amount);

   // outputColor = vec4(0.8, 0.2, 0.8, 1.0);

}

