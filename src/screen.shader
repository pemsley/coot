
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

// when we are close to the rotation centre, tightness is nearly 1.0
// when we are far (at the back), tightness is nearly 0.0
//
float get_weight_xy(float r_sqrd, float tightness) {

   float r = sqrt(r_sqrd); // between 0.0 and 1.0
   float t = clamp(tightness, 0.0, 1.0);
   float f_1 = 1.0 - 0.24 * t * r; // 0.2 needs tweaking
   float f_2 = clamp(f_1, 0.0, 1.0);
   return f_2;

}

// Note: the weight here imply (or are based on) linear gradient
// for the depth. If the projection matrix become perspective,
// then the depth will be non-linear and will need to be
// transformed.

// blur things at the back more than things at the front
// - how much should this pixel get blurred?
float get_weight_z_1(float depth_centre) {

    float w = abs(depth_centre - 0.158);
    return w;

}

// "this" pixel is contributing to the colour of another pixel.
//  What weight should this pixel have to colour that pixel?
float get_weight_z_2(float depth_this, float depth_centre) {

    float c = 0.52; // 0.66 looks good, but we,
                    // need a "protected" z region - say from  0.05 to 0.3
    float depth_delta = depth_this - depth_centre;
    float w = 0.0; // too far, by default.
    if (depth_delta >= -0.1) {
       w = (c - depth_delta)/c;
    } else {
       w = 0.0;
    }
    float weight = clamp(w, 0.0, 1.0);
    return weight;
}

vec3 sampling_blur() {

    // at some stage, fix the double addition when x =0 and when y = 0

    vec2 tex_scale = 1.0/textureSize(screenTexture, 0); // the size of single texel
    // the bigger the z value (it ranges from 0 (front) to 1.0 (back)), the more I want
    // to blur this pixel. Bigger z means more sampling
    vec3 sum = vec3(0,0,0); // return this
    int n_pixels_max = 8;
    float depth_centre = texture(screenDepth, TexCoords).x;
    float w_1 = get_weight_z_1(depth_centre); // w_1 is 0.0 at the front/no-blur
    float tightness = 1.0 - w_1;
    vec3 r = vec3(0.0, 0.125, 0.0);
    if (w_1 > 0.05 && w_1 < 0.3) { // the front (in sceen z) doesn't get blurred
       return texture(screenTexture, TexCoords).rgb;
    } else {
       // most of the image:
       int n_sampled = 0;
       int n_closer_neighbours = 0; // ambient occlusion (testing)
       float sum_weight = 0.0;
       for (int ix=0; ix<n_pixels_max; ix++) {
          for (int iy=0; iy<n_pixels_max; iy++) {
             float r_sqrd = (ix*ix + iy*iy) / (n_pixels_max * n_pixels_max);
             float weight_xy = get_weight_xy(r_sqrd, tightness);
             if (weight_xy > 0.0) {
                {
                   vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, tex_scale.y * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   float weight_z_2 = get_weight_z_2(depth_this, depth_centre);
                   float weight = weight_xy * weight_z_2;
                   if (weight > 0.0) {
                      sum += texture(screenTexture, offset_coords).rgb * weight;
                      n_sampled += 1;
                      sum_weight += weight;
                   }
                   if (depth_this < depth_centre) n_closer_neighbours++;
                }
                {
                   vec2 offset_coords = TexCoords + vec2(tex_scale.x * ix, -tex_scale.y * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   float weight_z_2 = get_weight_z_2(depth_this, depth_centre);
                   float weight = weight_xy * weight_z_2;
                   if (weight > 0.0) {
                      sum += texture(screenTexture, offset_coords).rgb * weight;
                      n_sampled += 1;
                      sum_weight += weight;
                   }
                   if (depth_this < depth_centre) n_closer_neighbours++;
                }
                {
                   vec2 offset_coords = TexCoords + vec2(-tex_scale.x * ix, tex_scale.y * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   float weight_z_2 = get_weight_z_2(depth_this, depth_centre);
                   float weight = weight_xy * weight_z_2;
                   if (weight > 0.0) {
                      sum += texture(screenTexture, offset_coords).rgb * weight;
                      n_sampled += 1;
                      sum_weight += weight;
                   }
                   if (depth_this < depth_centre) n_closer_neighbours++;
                }
                {
                   vec2 offset_coords = TexCoords + vec2(-tex_scale.x * ix, -tex_scale * iy);
                   float depth_this = texture(screenDepth, offset_coords).x;
                   float weight_z_2 = get_weight_z_2(depth_this, depth_centre);
                   float weight = weight_xy * weight_z_2;
                   if (weight > 0.0) {
                      sum += texture(screenTexture, offset_coords).rgb * weight;
                      n_sampled += 1;
                      sum_weight += weight;
                   }
                   if (depth_this < depth_centre) n_closer_neighbours++;
                }
             }
          }
       }
       if (n_sampled > 0) {
          r = texture(screenTexture, TexCoords).rgb;
          float f = 1.0;
          if ((2 * n_closer_neighbours) >= n_sampled) {
             float aos = float(n_closer_neighbours)/float(n_sampled);
             float f = 2.0 - 2.0 * aos;
             r  *=  f;
          }
          r = 0.25 * (r + 3 * sum/sum_weight);
          // r = sum/sum_weight;
       } else{
          r = vec3(1.0, 0.0, 0.0);
       }
       return r;
    }
}

void main() {

   vec3 result = vec3(0,0,0);

   if (false) {
      result = sampling_blur();
   } else {
      result = texture(screenTexture, TexCoords).rgb; // don't blur
   }

   if (false) {
      // test the depth - I want shades of grey for the moment.
      // result = vec3(z,z,z);
      // if (z == 1.0) result.r = 0.9;
      // if (z < 0.158) result = vec3(0.9, 0.0, 0.0);
   }

   out_color = vec4(result, 1.0);

}

