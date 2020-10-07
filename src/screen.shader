
#shader vertex

#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoords;

out vec2 TexCoords;

void main() {

   TexCoords = aTexCoords;
   gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);

}


#shader fragment

#version 330 core

in vec2 TexCoords;

uniform bool do_ambient_occlusion;
uniform sampler2D screenTexture;
uniform sampler2D screenDepth;

layout(location = 0) out vec4 out_color;

float rand(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

vec3 occlude() {

   vec3 r = vec3(0.0, 0.0, 0.2); // return r;
   bool do_long_scale = true;
   int n_pixels_max = 36;
   // most of the image:
   int n_sampled = 0;
   int n_closer_neighbours = 0; // ambient occlusion (testing)
   vec2 tex_scale = 1.0/textureSize(screenTexture, 0); // the size of single texel
   float depth_centre = texture(screenDepth, TexCoords).x;
   float closer_sum = 0.0;
   float bias = 0.01;
   for (int ix=0; ix<n_pixels_max; ix++) {
      for (int iy=0; iy<n_pixels_max; iy++) {
         if ((abs(ix)+abs(iy))> n_pixels_max) continue;
         // deal with double counting at some stage
         if (true) {
            {
               vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, tex_scale.y * iy);
               float depth_this = texture(screenDepth, offset_coords).x;
               n_sampled++;
               if ((depth_this + bias) < depth_centre) {
                  n_closer_neighbours++;
                  closer_sum += depth_this - depth_centre;
               }
            }
            {
               vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, -tex_scale.y * iy);
               float depth_this = texture(screenDepth, offset_coords).x;
               n_sampled++;
               if ((depth_this + bias) < depth_centre) {
                  n_closer_neighbours++;
                  closer_sum += depth_this - depth_centre;
               }
            }
            {
               vec2 offset_coords = TexCoords + vec2(-tex_scale.x * ix, tex_scale.y * iy);
               float depth_this = texture(screenDepth, offset_coords).x;
               n_sampled++;
               if ((depth_this + bias) < depth_centre) {
                  n_closer_neighbours++;
                  closer_sum += depth_this - depth_centre;
               }
            }
            {
               vec2 offset_coords = TexCoords + vec2(-tex_scale.x * ix, -tex_scale * iy);
               float depth_this = texture(screenDepth, offset_coords).x;
               n_sampled++;
               if ((depth_this + bias) < depth_centre) {
                  n_closer_neighbours++;
                  closer_sum += depth_this - depth_centre;
               }
            }
         }
      }
   }

   float scale_from_long = 1.0;
   int n_long_scale_samples = 100;
   int n_long_scale = 0;
   if (do_long_scale) {
      // 0.028 gives a nice/interesting/strange "outliney" look to zoomed out density
      float long_scale = 0.48; // what is the coordinates system?
      bool found_one = false;
      for(int i=0; i<n_long_scale_samples; i++) {
         // note to self - are these really random numbers - how can I test it?
         // test the mean and sd, perhaps.
         vec2 xy_1 = gl_FragCoord.xy + vec2(1.1 * float(i),  1.1 * float(i));
         vec2 xy_2 = gl_FragCoord.xy + vec2(1.1 * float(i), -1.1 * float(i));
         float r_1 = rand(xy_1) * 2.0 - 1.0;
         float r_2 = rand(xy_2) * 2.0 - 1.0;
         vec2 offset_coords = vec2(long_scale * tex_scale.x * r_1, long_scale * tex_scale.y * r_2);
         vec2 offseted_coords = TexCoords + offset_coords;
         float depth_sampled = texture(screenDepth, offseted_coords).x;
         if ((depth_sampled+bias) < depth_centre)
            n_long_scale++;
      }
      scale_from_long = 1.0 - 0.95 * float(n_long_scale)/float(n_long_scale_samples); // was 0.5
      scale_from_long = clamp(scale_from_long, 0.0f, 1.0f);
   }
   if (n_sampled > 0) {
      r = texture(screenTexture, TexCoords).rgb;
      if ((2 * n_closer_neighbours) >= n_sampled) {
         float aos = float(n_closer_neighbours)/float(n_sampled);
         aos = float(n_closer_neighbours)/float(n_sampled); // 0.5 to 1
         float f = 2.0 * aos - 1.0;  // 0.0 to 1.0 (very occluded to no occluded)
         float ff = 1.0 - f * 0.95; // was 0.7
         r *= ff;
         r *= scale_from_long;
         if (scale_from_long == 0.0) r = vec3(1.0, 0.0, 0.0);
         // if (float(n_long_scale)/float(n_long_scale_samples) > 0.9) r = vec3(1.0, 0.0, 0.0);
      }
   } else{
      r = vec3(0.0, 1.0, 0.0);
   }
   return r;
}


void main() {

   vec3 result = vec3(0,0,0);

   // result      = texture(screenTexture, TexCoords).rgb;
   float depth = texture(screenDepth,   TexCoords).x;
   gl_FragDepth = depth; // needed for depth for next shader

   // bool do_occlude = true;
   if (do_ambient_occlusion) {
      result = occlude();
   } else {
      result = texture(screenTexture, TexCoords).rgb;
   }

   out_color = vec4(result, 1.0);

}

