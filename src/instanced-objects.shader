/*
 * src/instanced-objects.shader
 *
 * Copyright 2020 by Medical Research Council
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

// use the #shader "directive" to separate the shaders on parsing
// -----------------------------------------
#shader vertex

#version 330 core

// This is instanced-objects.shader

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
uniform mat4 view_rotation;
uniform float time;
uniform bool transferred_colour_is_instanced;
uniform bool do_pulse;
uniform bool do_rotate_z;
uniform float pulsing_amplitude; // = 0.25;  // make these uniforms?
uniform float pulsing_frequency; // = 5.0;
uniform float pulsing_phase_distribution; // = 0.2;
uniform float z_rotation_angle;
uniform float stereo_x_scale;
uniform float stereo_x_offset;

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

   vec3 t_pos = position;
   vec3 n_dir = normal;
   // the normal doesn't change for translation along z (pulsing)
   if (do_pulse)
      t_pos = position + vec3(0,0, pulsing_amplitude * sin(0.001 * pulsing_frequency * time + pulsing_phase_distribution * gl_InstanceID));
   if (do_rotate_z) {
      float cos_theta = cos(time * 0.014 * z_rotation_angle);
      float sin_theta = sin(time * 0.014 * z_rotation_angle);
      t_pos = vec3(position.x * cos_theta - position.y * sin_theta, position.x * sin_theta + position.y * cos_theta, position.z);
      n_dir = vec3(  normal.x * cos_theta -   normal.y * sin_theta,   normal.x * sin_theta +   normal.y * cos_theta, normal.z);
   }

   vec4 p4 = vec4(t_pos, 1.0);
   vec4 frag_pos = model_rotation_translation_scale * p4;

   gl_Position = mvp * frag_pos * vec4(stereo_x_scale, 1.0f, 1.0f, 1.0f) + vec4(stereo_x_offset, 0.0f, 0.0f, 0.0f);

   normal_transfer = model_rotation * n_dir;
   colour_transfer = colour;
   if (transferred_colour_is_instanced) colour_transfer = colour_instanced;
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
   vec4 ambient;
   vec4 diffuse;
   vec4 specular;
};

uniform mat4 view_rotation;
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
   vec4 running_col = vec4(0,0,0,0);

   float fog_amount = get_fog_amount(gl_FragCoord.z);

   for (int i=0; i<2; i++) {
      if (light_sources[i].is_on) {

         // ambient
         vec4 ambient = ct * light_sources[i].ambient * material.ambient;

         // diffuse
         vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space.xyz;
         vec3 norm_2 = normalize(normal_transfer); // not needed?

         // light_dir = vec3(0,0,1);

         float dp_raw = dot(norm_2, light_dir);
         float dp = max(dp_raw, 0.0);
         vec4 diffuse = ct * light_sources[i].diffuse * dp * 5.0 * material.diffuse;

         // specular

         float shininess = material.shininess;
         float specular_strength = material.specular_strength;
         if (dp_raw < 0.0)
            specular_strength = 0.0; // no shiny interiors
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
         // spec = 0;
         vec4 specular = 3.0 * spec * light_sources[i].specular;

         // final
         running_col += ambient + diffuse + specular;

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
