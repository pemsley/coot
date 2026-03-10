/*
 * src/meshes.shader
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
// meshes.shader

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 colour;

// rama-balls material change would need a float here for the shininess manipulater - hmm!
// maybe it needs its own layout and shader. Oh. They do! rama-balls.shader. Hmm!
// So how do we use that for a Mesh for intermediate atoms?

out vec3 frag_pos_transfer;
out vec3 normal_transfer;
out vec4 colour_transfer;

uniform mat4 mvp;
uniform mat4 view_rotation;
uniform float stereo_x_scale;
uniform float stereo_x_offset;

void main() {
   vec4 n = vec4(normal, 1.0);
   // gl_Position = mvp * vec4(position, 1.0);
   gl_Position = mvp * vec4(position, 1.0) * vec4(stereo_x_scale, 1.0f, 1.0f, 1.0f) + vec4(stereo_x_offset, 0.0f, 0.0f, 0.0f);
   normal_transfer = normal; // Hmmm! 20220209-PE normals are in "molecule" space (as are the light positions)
   colour_transfer = colour;
   frag_pos_transfer = position;
}


// -----------------------------------------
#shader fragment

#version 330 core
// meshes.shader

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
   vec4 ambient;
   vec4 diffuse;
   float shininess;
   float specular_strength;
   vec4 specular;
};
uniform LightSource light_sources[2];
uniform vec3 eye_position;
uniform vec3 eye_position_in_molecule_coordinates_space; // sent by Shader::setup_light()
uniform mat4 view_rotation;
uniform Material material;
uniform bool is_perspective_projection;
uniform vec4 background_colour;
uniform float opacity; // was map_opacity
uniform bool do_specular;
uniform bool do_depth_fog;
uniform bool do_fresnel;
uniform float fresnel_bias;
uniform float fresnel_scale;
uniform float fresnel_power;
uniform vec4  fresnel_colour;

in vec3 frag_pos_transfer;
in vec3 normal_transfer;
in vec4 colour_transfer;

out vec4 output_colour;

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
         return d4;
      }
   } else {
      return 0.0;
   }
}


void main() {

   // make this true for Coot 0.x-like shading
   bool flat_shading = false;

   if (flat_shading) {
      output_colour = colour_transfer;
   } else {

      float fog_amount = get_fog_amount(gl_FragCoord.z);

      output_colour = vec4(0,0,0,0); // or do I mean 0,0,0,1 ?

      float specular_strength = material.specular_strength;
      float shininess = material.shininess;

      for (int i=0; i<2; i++) {
         // if (i > 0) continue; // debugging lights
         if (light_sources[i].is_on) {
            vec3 light_dir = light_sources[i].direction_in_molecule_coordinates_space;
            float dp = dot(normal_transfer, light_dir);

            if (dp < 0.0)
               specular_strength = 0.0;

            // We can't have specular lights where there is no diffuse light
            //
            // 20220115-PE That's true, but it might help debugging the specularity
            // if this is left on for the moment. Question to self: does the
            // specular shading in texture-meshes.shader work correctly? it calculates
            // light_dir differently, like this:
            //         mat4 ivr = transpose(view_rotation);
            //         vec3 light_dir = (vec4(light_sources[i].direction_in_molecule_coordinates_space, 1.0) * ivr).xyz;
            //  I tried that here but it doesn't work right.
            //
            // if (dp <= 0.0)
            // specular_strength = 0.0;

            dp = clamp(dp, 0.0, 1.0); // no negative dot products for diffuse

            // ambient
            vec4 ambient = light_sources[i].ambient * material.ambient * colour_transfer;

            // diffuse
            vec4 diffuse = light_sources[i].diffuse * dp * material.diffuse * colour_transfer;

            // specular - this is wrong - can result in stretched specular spots.
            //            It looks like the code in moleculestotriangles.shader - does that
            //            do the specular shading correctly?
            //
#if 0 // 20220115-PE as was
            vec3 eye_pos = eye_position;
            vec3 norm_2 = normalize(normal_transfer); // not needed?
            vec3 view_dir = normalize(eye_pos - frag_pos_transfer.xyz); // frag_pos_transfer is a vec3 in model.shader
            vec3 reflect_dir = reflect(-light_dir, norm_2);
            reflect_dir = normalize(reflect_dir); // belt and braces
            float dp_view_reflect = dot(view_dir, reflect_dir);
            dp_view_reflect = clamp(dp_view_reflect, 0.0, 1.0);
            float spec = pow(dp_view_reflect, shininess);
            vec4 specular = specular_strength * spec * light_sources[i].specular;
#endif

#if 1 // largely from moleculestotriangles.shader, but using a different eye position space.

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

            output_colour += ambient + diffuse + specular;

         }
      }

      vec4 bgc = background_colour;
      output_colour = mix(output_colour, bgc, fog_amount);
      output_colour.a = opacity;

      // also I think it's a good idea to clamp the output_colour here to max of 0.99
      // so that wild specularity is suppressed.

   }


   if (do_fresnel) {
      // R can be calculated in the vertex shader
      vec3 eye_to_frag_pos_uv = normalize(frag_pos_transfer - eye_position_in_molecule_coordinates_space);
      float dp_eye = dot(normalize(normal_transfer), eye_to_frag_pos_uv);
      // (I.N) should be 1.0 if we are looking staight on to (perpendicular to) the surface.
      float I_dot_N = dp_eye;
      if (I_dot_N <= 0.0) {
         // I_dot_N = clamp(I_dot_N, 0.0, 1.0); // should not be needed.
         float R0 = fresnel_bias + fresnel_scale * pow((1.0 + I_dot_N), fresnel_power);
         float R = clamp(R0, 0.0, 1.0);
         float gr = 0.8; // grey
         // vec4 fresnel_light_colour = vec4(gr, gr, gr, 1.0);
         vec4 fresnel_light_colour = fresnel_colour;

         vec4 colour_from_fresnel = R * fresnel_light_colour;

         output_colour += colour_from_fresnel;
      }
   }


   // testing hacks:

   // output_colour = vec4(normal_transfer, 1.0);

   // material.specular_strength = 0.4
   // if (material.specular_strength > 0.4)  // great than 0.2, 0.3  less than 0.4
   // output_colour = vec4(1,0,1,1);

   // output_colour = vec4(1,0,1,1);

 
}

