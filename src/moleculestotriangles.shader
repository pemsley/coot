/*
 * src/moleculestotriangles.shader
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

// 20220213-PE question to self: what does this do that meshes.shader does not do?

#shader vertex

#version 330 core
// moleculestotriangles.shader

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

uniform mat4 mvp;
uniform float stereo_x_scale;
uniform float stereo_x_offset;

out vec4 colour_transfer;
out vec3 normal_transfer;
out vec3 frag_pos_transfer;
out mat3 model_rotation_transfer;

void main() {
   vec4 n = vec4(normal, 1.0);
   // gl_Position = mvp * vec4(position, 1.0);
   gl_Position = mvp * vec4(position, 1.0) * vec4(stereo_x_scale, 1.0f, 1.0f, 1.0f) + vec4(stereo_x_offset, 0.0f, 0.0f, 0.0f);
   normal_transfer = normalize(normal); // probably the normalize() here is not needed
   colour_transfer = colour;
   frag_pos_transfer = position;
   model_rotation_transfer = mat3(mvp);
}


// -----------------------------------------
#shader fragment

#version 330 core
// moleculestotriangles.shader

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
uniform bool do_specular;
uniform float opacity;
uniform bool show_shadows;
uniform sampler2D shadow_map;

struct Material {
   vec4 ambient;
   vec4 diffuse;
   vec4 specular;
   float shininess;
   float specular_strength;
};
uniform Material material;

in vec4 colour_transfer;
in vec3 normal_transfer;
in vec3 frag_pos_transfer;
in mat3 model_rotation_transfer;

out vec4 outputColor;


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
   outputColor = vec4(0,0,0,0);

   float specular_strength = material.specular_strength;
   float shininess = material.shininess;

   for (int i=0; i<2; i++) {
      // if (i > 0) continue;
      if (light_sources[i].is_on) {
         vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space;
         float dp = dot(normal_transfer, light_dir);
         // we can't have specular lights where there is no diffuse light
         if (dp <= 0.0)
            specular_strength = 0.0;
         dp = clamp(dp, 0.1, 1.0); // no negative dot products for diffuse

         vec4 lsa = light_sources[i].ambient;
         vec4 lsd = light_sources[i].diffuse;
         vec4 ambient = colour_transfer * lsa * 0.1; // I scale this down here - it's bright. why?
         vec4 diffuse = colour_transfer * lsd * dp;

         // specular
         vec3 eye_pos = eye_position * transpose(model_rotation_transfer);
         vec3 norm_2 = normal_transfer; // normalize(normal_transfer); // not needed, I think
         vec3 view_dir = normalize(eye_pos - frag_pos_transfer);
         vec3 reflect_dir = reflect(light_dir, norm_2);
         // reflect_dir = normalize(reflect_dir); // belt and braces
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

   // outputColor.a = 0.1;

   if (show_shadows)
      outputColor = vec4(0.7, 0.2, 0.1, 0.2);

}
