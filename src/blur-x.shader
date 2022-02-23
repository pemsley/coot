
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

uniform sampler2D colourTexture;
uniform sampler2D  depthTexture;

layout(location = 0) out vec4 out_colour;

void main() {

   vec3  sampled       = texture(colourTexture, TexCoords).rgb;
   float depth_sampled = texture( depthTexture, TexCoords).r;
   float kern[11];
   kern[ 0] = 1.0;
   kern[ 1] = 0.9607894391523232;
   kern[ 2] = 0.8521437889662113;
   kern[ 3] = 0.697676326071031;
   kern[ 4] = 0.5272924240430485;
   kern[ 5] = 0.36787944117144233;
   kern[ 6] = 0.23692775868212176;
   kern[ 7] = 0.14085842092104503;
   kern[ 8] = 0.07730474044329971;
   kern[ 9] = 0.039163895098987066;
   kern[10] = 0.01831563888873418;
   vec2 tex_scale = 1.0/textureSize(colourTexture, 0); // the size of single texel

   int n_pixels_max = 10;
   vec3 sum = vec3(0.0, 0.0, 0.0);
   float weight_sum = 0.0;
   for (int ix=-n_pixels_max; ix<=n_pixels_max; ix++) {
      float k = kern[abs(ix)];
      float weight = 1.0;
      vec2 offset_coords = TexCoords + vec2(1.0 * tex_scale.x * ix, 0.0);
      vec3 t = texture(colourTexture, offset_coords).rgb;
      sum += t * vec3(k,k,k) * weight;
      weight_sum += weight;
   }
   vec3 result = 2.0 * sum/weight_sum;

   gl_FragDepth = depth_sampled;

   // result = texture(colourTexture, TexCoords).rgb;
   out_colour = vec4(result, 1.0);

   // out_colour = vec4(1,0,1,1);

   // if (depth_sampled > 0.49)
   // out_colour = vec4(1,0,0,1);
   // else
   // out_colour = vec4(0,0,1,1);

   // out_colour = vec4(sampled, 1.0);
   float ds = depth_sampled;
   // out_colour = vec4(ds, ds, ds, 1.0);
   // if (ds < 0.8 ) out_colour = vec4( 0, ds, ds,  1.0);
   // if (ds > 0.99) out_colour = vec4(ds, ds, 1.0, 1.0);

   // testing, does the blur-y shader do anything? Let's by-pass the x-shader for now.
   // out_colour = vec4(sampled, 1.0);

}

