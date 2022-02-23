
#shader vertex
// effects.shader

#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoords;

out vec2 TexCoords;

void main() {

   TexCoords = aTexCoords;
   gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);

}


#shader fragment
// effects.shader

#version 330 core

in vec2 TexCoords;

uniform sampler2D screenTexture;
uniform sampler2D screenDepth;
uniform sampler2D ssao;

uniform bool use_ssao;
uniform float ssao_strength;
uniform vec4 background_colour;

// show *just* the ssao, that is
uniform bool show_ssao;
uniform int effects_output_type;

layout(location = 0) out vec4 out_color;


float get_fog_amount(float depth_in) {

   bool do_depth_fog = true;
   bool is_perspective_projection = true;
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
         d4 = d;
         return d4;
      }
   } else {
      return 0.0;
   }
}


void main() {

   float depth = texture(screenDepth,   TexCoords).x;

   // this line make a difference. If it is not there, we are
   // editing the whole texture. If it is there, we are editing
   // only the pixels of the object (not the background).
   //

   // gl_FragDepth = depth; // needed for depth for next shader // test remove

   float ao_mapped = 1.0;
   out_color = vec4(0,0,0,1);

   if (use_ssao) {
      float ao = texture(ssao, TexCoords).r;
      float inv_ao = 1.0 - ao; // now 1.0 is strong AO, 0.0 is no AO
      float str_inv_ao = 1.0 * ssao_strength * inv_ao;
      ao_mapped = 1.0 - str_inv_ao;
   }

   if (show_ssao) {
      out_color = vec4(vec3(ao_mapped), 1.0);
   } else {
      vec3 t = texture(screenTexture, TexCoords).rgb;
      vec4 t4 = vec4(t, 1.0);
#if 0
      float aafa = ao_mapped * get_fog_amount(gl_FragCoord.z);
      vec4 fogged_ssao = vec4(vec3(aafa, aafa, aafa), 1.0);
      out_color = 2.4 * t4 * fogged_ssao;
#endif
#if 1
      float compensating_scale = 1.0 + 0.1 * ssao_strength; // use more?

      // this doesn't work becuse if I add it there is a big colour descrepancy
      // betwween pixels at the back (0.99) and background pixels at 1.0
      // if (depth == 1.0)
      // compensating_scale = 1.0;

      out_color = compensating_scale * ao_mapped * t4;
      // Yikes! gl_FragCoord.z == 0.5 for all pixels.
      // if (gl_FragCoord.z == 0.5)
      // out_color = vec4(1,0,0,1);
      // maybe I can use depth, texture(screenDepth, TexCoords)?
#endif
      // out_color = texture(screenTexture, TexCoords);
      // out_color = vec4(vec3(ao_mapped), 1.0);
   }

   // see graphics-info.h effects_shader_output_type

   if (effects_output_type == 0) {
      // standard
   }
   if (effects_output_type == 1) {
      out_color = texture(screenTexture, TexCoords);
   }
   if (effects_output_type == 2) {
      out_color = vec4(vec3(texture(ssao, TexCoords).r), 1.0);
   }
   if (effects_output_type == 3) {
      out_color = vec4(vec3(texture(screenDepth, TexCoords).r), 1.0);
   }

   // test hack!
   // out_color = texture(screenTexture, TexCoords);
   // out_color = texture(screenDepth, TexCoords);


   // out_color = vec4(1,1,0,1);

   // out_color = vec4(texture(ssao, TexCoords).r);
   // out_color = 0.001 * gl_FragCoord;
   // out_color = texture(screenTexture, TexCoords);

   // out_color = texture(screenDepth, TexCoords);


}

