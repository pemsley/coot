
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

layout(location = 0) out vec4 out_colour;

void main() {

   vec3 result = texture(screenTexture, TexCoords).rgb;
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
   vec2 tex_scale = 1.0/textureSize(screenTexture, 0); // the size of single texel

   int n_pixels_max = 10;
   vec3 sum = vec3(0.0, 0.0, 0.0);
   float weight_sum = 0.0;
   for (int iy=-n_pixels_max; iy<=n_pixels_max; iy++) {
      float k = kern[abs(iy)];
      vec2 offset_coords = TexCoords + vec2(0.0, tex_scale.y * iy);
      float weight = 1.0;
      vec3 t = texture(screenTexture, offset_coords).rgb;
      sum += t * vec3(k,k,k) * weight;
      weight_sum += weight;
   }

   result = 2.4 * sum/weight_sum;
   // result = texture(screenTexture, TexCoords).rgb;

   out_colour = vec4(result, 1.0);

   // out_colour = vec4(0,1,1,0);

}

