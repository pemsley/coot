
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

uniform sampler2D screenTexture;
uniform sampler2D screenDepth;

layout(location = 0) out vec4 out_color;


vec3 occlude() {

       vec3 r = vec3(0.0, 0.0, 0.2); // return r;
       int n_pixels_max = 8;
       // most of the image:
       int n_sampled = 0;
       int n_closer_neighbours = 0; // ambient occlusion (testing)
       vec2 tex_scale = 1.0/textureSize(screenTexture, 0); // the size of single texel
       float depth_centre = texture(screenDepth, TexCoords).x;
       float closer_sum = 0.0;
       for (int ix=0; ix<n_pixels_max; ix++) {
          for (int iy=0; iy<n_pixels_max; iy++) {
             // deal with double counting at some stage
             if (true) {
                {
                   vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, tex_scale.y * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   n_sampled++;
                   if (depth_this < depth_centre) {
                      n_closer_neighbours++;
                      closer_sum += depth_this - depth_centre;
                   }
                }
                {
                   vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, -tex_scale.y * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   n_sampled++;
                   if (depth_this < depth_centre) {
                      n_closer_neighbours++;
                      closer_sum += depth_this - depth_centre;
                   }
                }
                {
                   vec2 offset_coords = TexCoords + vec2(-tex_scale.x * ix, tex_scale.y * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   n_sampled++;
                   if (depth_this < depth_centre) {
                      n_closer_neighbours++;
                      closer_sum += depth_this - depth_centre;
                   }
                }
                {
                   vec2 offset_coords = TexCoords + vec2(-tex_scale.x * ix, -tex_scale * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   n_sampled++;
                   if (depth_this < depth_centre) {
                      n_closer_neighbours++;
                      closer_sum += depth_this - depth_centre;
                   }
                }
             }
          }
       }
       // Note to self, ambient occlusion looks good when the background is
       // black - when it is grey, the ambient occlusion needs to be less
       if (n_sampled > 0) {
          r = texture(screenTexture, TexCoords).rgb;
          if ((2 * n_closer_neighbours) >= n_sampled) {
             float aos = float(n_closer_neighbours)/float(n_sampled);
             aos = closer_sum;
             float f = 1.0 + 0.05 * aos;
             f = clamp(f, 0.0, 1.0);
             r *=  f;
             // r = vec3(0.0, 1.0, 0.0);
          }
       } else{
          r = vec3(1.0, 0.0, 0.0);
       }
       return r;
}


void main() {

   vec3 result = vec3(0,0,0);

   // result      = texture(screenTexture, TexCoords).rgb;
   float depth = texture(screenDepth,   TexCoords).x;
   gl_FragDepth = depth;
   if (depth > 1.0) {
      result = vec3(0,0,0);
   } else {
      result = occlude();
   }

   // result = vec3(depth, depth, depth);
   // result = vec3(0.3, 0.0, 0.3); // dark purple

   out_color = vec4(result, 1.0);

}

