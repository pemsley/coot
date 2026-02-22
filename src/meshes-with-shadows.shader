/*
 * src/meshes-with-shadows.shader
 *
 * Copyright 2022 by Medical Research Council
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
// This is meshes-with-shadows.shader

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

out vec3 frag_pos_transfer;
out vec4 frag_pos_light_space_transfer;
out vec3 normal_transfer;
out vec4 colour_transfer;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform mat4 light_space_mvp;

void main() {
   vec4 n = vec4(normal, 1.0);
   gl_Position = mvp * vec4(position, 1.0);
   normal_transfer = normal;  // this does not match texture-meshes-with-shadows.shader
   mat3 view_rotation_3 = mat3(view_rotation);
   mat3 transpose_view_matrix = transpose(view_rotation_3);
   normal_transfer = normal;  // * transpose_view_matrix; // 20220326-PE these normals change the colour of the surface
                                                          // dependent on the view - so that blue is always facing the  viewer
   normal_transfer = normal; // this is like meshes.shader - the dp is calculated using eye in molecule space
   frag_pos_transfer = position;
   frag_pos_light_space_transfer = light_space_mvp * vec4(position, 1.0);
   colour_transfer = colour;

   // hack test
   // colour_transfer = vec4(normal_transfer, 1.0);
}


// -----------------------------------------
#shader fragment

// This is meshes-with-shadows.shader

#version 330 core

in vec3 frag_pos_transfer;
in vec4 frag_pos_light_space_transfer;
in vec3 normal_transfer;
in vec4 colour_transfer;

out vec4 output_colour;

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
uniform bool do_fresnel;
uniform float shadow_strength;
uniform int shadow_softness;
uniform float fresnel_bias;
uniform float fresnel_scale;
uniform float fresnel_power;
uniform vec4  fresnel_colour;


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

}

void main() {

   vec4 sum_col = vec4(0,0,0,0);

   if (show_shadows) { // show *just* shadows

      mat4 ivr = transpose(view_rotation);
      vec3 light_dir = (vec4(light_sources[0].direction_in_molecule_coordinates_space, 1.0) * ivr).xyz;
      float dp_raw = dot(normal_transfer, light_dir);
      float shadow = calc_shadow(dp_raw); // calc_shadow() returns 0.0 for full shadow

      vec3 pos_light_space = frag_pos_light_space_transfer.xyz * 0.5 + 0.5; // 0 to 1
      sum_col = vec4(vec3(shadow), 1.0);

      if (true) {
	 if (pos_light_space.z == 0.0 ) sum_col = vec4(1,1,0,1);
	 if (pos_light_space.z >  0.75) sum_col = vec4(1,0,0,1);
	 if (pos_light_space.z <  0.0 ) sum_col = vec4(0,0,1,1);
      }

   } else {

      float scale_factor_n_lights = 1.2; // was 0.75. why is the dwarf blacksmith dark?

      output_colour = vec4(0,0,0,0);
      vec4 bg_col = background_colour;

      float specular_strength = material.specular_strength;
      if (do_specular) {
         // no change
      } else {
         specular_strength = 0.0;
      }
      float shininess = material.shininess;

      // specular_strength *= 32;
      // shininess *= 8;

      // specular_strength = 0.0;

      float fog_amount = get_fog_amount(gl_FragCoord.z);

      for (int i=0; i<2; i++) {
	 if (light_sources[i].is_on) {
	    vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space;
	    float dp_raw = dot(normal_transfer, light_dir);
            float dp = dp_raw;
            if (dp_raw <= 0.0)
               specular_strength = 0.0; // we can't have specular lights where there is no diffuse light.
                                        // no weird shiny ball interiors

            // testing
            // for "flat" old-coot look, set dp = 0.0 and specular_strength = 0.0 here.
            // dp = 0.0;
            // specular_strength = 0.0;

	    dp = clamp(dp, 0.1, 1.0); // no negative dot products for diffuse

            // specular_strength = specular_strength * 0.0;

            // there is a mishmash of colouring schemes here. The vertex has a colour
            // but we also have access to the ambient and diffuse Material colours
            // (and they can be changed via the API, whereas the vertex colour cannot).
            // Let's try to handle both. A bit messy.

	    vec4 ambient = light_sources[i].ambient * scale_factor_n_lights * colour_transfer * material.ambient;
	    vec4 diffuse = light_sources[i].diffuse * scale_factor_n_lights * colour_transfer * material.diffuse * dp;
	    // diffuse = dp * vec4(0.5, 0.5, 0.5, 1.0);

#if 1 // largely from moleculestotriangles.shader, but using a different eye position space.
      // Also the same as meshes.shader.

            // specular
            vec3 eye_pos = eye_position; // eye position is (say) (0,0,40) - it's in view space, but the
                                         // frag_pos_transfer (position) is in model (molecular coordinates) space

            eye_pos = eye_position_in_molecule_coordinates_space;

            vec3 norm_2 = normal_transfer;
            vec3 view_dir = normalize(eye_pos - frag_pos_transfer);
            vec3 reflect_dir = reflect(-light_dir, norm_2);
            reflect_dir = normalize(reflect_dir); // belt and braces
            float dp_view_reflect = dot(view_dir, reflect_dir);
            dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
            float spec = pow(dp_view_reflect, shininess);
            vec4 specular = specular_strength * spec * light_sources[i].specular;
            if (! do_specular)
               specular = vec4(0,0,0,0);
#endif

            // Don't allow specularity if the pixel is in any amount of shadow.
            //
	    float shadow_diffuse = calc_shadow(dp); // calc_shadow() returns 0.0 for full shadow
            float shadow_specular = shadow_diffuse;
            if (shadow_specular < 0.95)
               shadow_specular = 0;

            vec4 add_col = ambient + diffuse * shadow_diffuse + specular * shadow_specular * specular_strength;
	    // vec4 add_col = ambient + diffuse * shadow_diffuse;

            // is this right? It looks a bit weird - but maybe it's the effects shader...

            add_col = mix(add_col, background_colour, 0.0000001 * fog_amount);
            add_col.a = opacity;

            // add_col = vec4(vec3(shadow_diffuse), 1.0);

	    sum_col += add_col;

            // if (shadow_diffuse == 1.0)
            //        sum_col = vec4(1,0,1,1);

	    // shadow is not being sampled correctly. Is the uniform for the sampler2d
	    // being set? (does it need to be?)

	    // if (shadow == 1.0)
	    // sum_col = vec4(1,1,0,1);

	    // vec3 light_to_eye = normalize(eye_pos - 100.0 * light_dir);
	    // sum_col = vec4(0.5 * view_dir + vec3(0.5,0.5,0.5), 1.0);
	    // sum_col = vec4(0.5 * light_to_eye + vec3(0.5,0.5,0.5), 1.0);
	 }
      }

      if (do_fresnel) {
         // R can be calculated in the vertex shader
         vec3 eye_to_frag_pos_uv = normalize(frag_pos_transfer - eye_position_in_molecule_coordinates_space);
         float dp_eye = dot(normalize(normal_transfer), eye_to_frag_pos_uv);
         // (I.N) should be 1.0 if we are looking staight on to (perpendicular to) the surface.
         float I_dot_N = dp_eye;
         if (I_dot_N <= 0.0) {
            float bias  = 0.01;
            float scale = 0.5;
            float power = 2.2;
            // I_dot_N = clamp(I_dot_N, 0.0, 1.0); // should not be needed.
            float R0 = fresnel_bias + fresnel_scale * pow((1.0 + I_dot_N), fresnel_power);
            float R = clamp(R0, 0.0, 1.0);
            float gr = 0.8; // grey
            // vec4 fresnel_light_colour = vec4(gr, gr, gr, 1.0);
            vec4 fresnel_light_colour = fresnel_colour;

            vec4 colour_from_fresnel = R * fresnel_light_colour;

            sum_col += colour_from_fresnel;
         }
      }

      sum_col = mix(sum_col, background_colour, fog_amount);

      if (fog_amount > 1.0) sum_col = vec4(1,0,1,1);
      if (fog_amount < 0.0) sum_col = vec4(0,1,0,1);

   }

   output_colour = sum_col;

   // 20220326-PE testing the normals:
   // output_colour = vec4(normal_transfer, 1.0);

   // output_colour = colour_transfer;
   // output_colour = vec4(1,0,1,1);

}

