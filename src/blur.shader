
#shader vertex

#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoords;

uniform float zoom;

out vec2 TexCoords;
out float zoom_transfer;

void main() {

   TexCoords = aTexCoords;
   zoom_transfer = zoom;
   gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);

}


#shader fragment

#version 330 core

in vec2 TexCoords;
in float zoom_transfer;

uniform sampler2D screenTexture;
uniform sampler2D screenDepth;

layout(location = 0) out vec4 out_color;

float depth_scale(float depth_in_centre, float depth_in_ij) {

   // -1 means there is no blurring to be done.

   // xx_ij is the point being blurred. xx_centre is the point
   // being blurred into. Keep a clear head.

   // the deeper into the image we go, (depth_in_centre approaches 1.0)
   // then the more we want to the colour to be influenced by the surroundings.

   float lim = 0.5;
   if (depth_in_centre < lim) {
      return -1.0;
   } else {
      if (depth_in_ij < lim) {
         return -1.0;
      } else {
         return (depth_in_ij-lim)/(1.0-lim);
      }
   }
}

vec3 sampling_blur(int n_pixels_max) {

   float depth_centre = texture(screenDepth, TexCoords).x;
   vec3 result = vec3(1.0, 1.0, 0.0);

   // centre is the point being blurred *into*

   vec3 orig_colour = texture(screenTexture, TexCoords).rgb; // don't blur
   float lim = 0.5;
   if (depth_centre < lim) {
      return orig_colour;
   } else {
      vec3 sum_outer = vec3(0,0,0);
      vec3 sum_inner = vec3(0,0,0);
      vec2 tex_scale = 1.0/textureSize(screenTexture, 0);
      // the 0.03 here depends on how much of the origin colour we add back
      // that is currently 0.5
      float centre_points_scale = 0.03; // 0.225;
      int n_inner_neighbs = 0;
      int n_outer_neighbs = 0;
      for (int ix= -n_pixels_max; ix<=n_pixels_max; ix++) {
         for (int iy= -n_pixels_max; iy<=n_pixels_max; iy++) {
            float r_sqrd = float(ix*ix + iy*iy) / float(n_pixels_max * n_pixels_max);
            if (r_sqrd > 1.0) continue;
            vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, tex_scale.y * iy);
            float depth_ij = texture(screenDepth, offset_coords).x;
            if (depth_ij == 1.0) continue;
            vec3 colour_ij = texture(screenTexture, offset_coords).rgb;
            // depth_scale() return -1 for no blurring.
            float dbrs = depth_scale(depth_centre, depth_ij);
            if (ix == 0 && iy == 0) {
               // sum_inner += colour_ij;
            } else {
               if (abs(ix) < 3 && abs(iy) < 3) {
                  sum_inner += colour_ij;
                  n_inner_neighbs++;
               }
               if (dbrs < 0.0) {
                  // nothing
               } else {
                  float blur_radius = (depth_ij - lim) / (1.0 - lim);
                  // this scaling needs to be zoom dependent.
                  float k = 0.001; // gaussian scale
                  dbrs = 1.0;
                  float gauss = 1.0/dbrs * exp(-0.5 * k * r_sqrd/(blur_radius*blur_radius));
                  float depth_factor = 2.0 * (1.0 - depth_centre); // 0 -> 1
                  depth_factor = 1.0;
                  sum_outer += colour_ij * gauss * (1.0 - depth_ij) * depth_factor;
                  n_outer_neighbs++;
               }
            }
         }
      }
      float Sc = 1.5/float(n_pixels_max*n_pixels_max);
      // result = 0.5 * orig_colour + 0.1 * centre_points_scale * sum_inner + Sc * sum_outer;
      // result = 0.1 * orig_colour + 0.1 * centre_points_scale * sum_inner + Sc * sum_outer;
      result = 0.2 * orig_colour + 0.6 * centre_points_scale * sum_inner + Sc * sum_outer;
      if (n_inner_neighbs == 0 && n_outer_neighbs == 0)
         result = orig_colour;
   }
   return result;
}

vec3 make_outline() {
   float depth_centre = texture(screenDepth, TexCoords).x;
   vec3 orig_colour = texture(screenTexture, TexCoords).rgb; // don't blur
   vec3 result = orig_colour; // update this as needed.
   vec2 tex_scale = 1.0/textureSize(screenTexture, 0);
   int n_deep_neighbs = 0;
   for (int ix= -1; ix<=1; ix++) {
      for (int iy= -1; iy<=1; iy++) {
         vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, tex_scale.y * iy);
         float depth_ij = texture(screenDepth, offset_coords).x;
         if ((depth_ij - depth_centre) > 0.1) {
            n_deep_neighbs++;
         }
      }
   }
   if (n_deep_neighbs > 1) {
      result = vec3(0.1, 0.1, 0.1);
   }
   return result;
}


void main() {

   vec3 result = vec3(0,0,0);

   bool do_depth_blur = false;
   bool do_outline    = true;

   if (do_depth_blur) {
      result = sampling_blur(14);
   } else {
      if (do_outline) {
         result = make_outline();
      } else {
         float depth_centre = texture(screenDepth, TexCoords).x;
         result = texture(screenTexture, TexCoords).rgb; // don't blur
         if (zoom_transfer < 10.0)
            if (depth_centre > 0.99)
               result = vec3(0.2, 0.3, 0.1);
      }
   }

   out_color = vec4(result, 1.0);

}
