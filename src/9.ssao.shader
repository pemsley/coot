
#shader vertex
// This is 9.ssao.shader
// ---------------- 9.ssao.shader shaderGeometryPass ---------------------

#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoords;

out vec2 TexCoords;

void main()
{
    TexCoords = aTexCoords;
    gl_Position = vec4(aPos, 1.0);
}

#shader fragment
// This is 9.ssao.shader

// ---------------- 9.ssao.shader shaderGeometryPass ---------------------

#version 330 core
out float FragColor;

in vec2 TexCoords;

uniform sampler2D gPosition;
uniform sampler2D gNormal;
uniform sampler2D texNoise;
uniform int n_ssao_kernel_samples;

uniform mat4 projection;
uniform vec3 samples[512]; // 512 is max, only up to n_ssao_kernel_samples is acutally set.
uniform float radius;
uniform float bias;

// tile noise texture over screen based on screen dimensions divided by noise size
float x_size = textureSize(gPosition, 0).x;
float y_size = textureSize(gPosition, 0).y;
vec2 noiseScale = vec2(x_size/4.0, y_size/4.0);

void main() {

   // get input for SSAO algorithm
   vec3 fragPos = texture(gPosition, TexCoords).xyz;
   vec3 normal = normalize(texture(gNormal, TexCoords).rgb);
   vec3 randomVec = normalize(texture(texNoise, TexCoords * noiseScale).xyz);
   // create TBN change-of-basis matrix: from tangent-space to view-space
   vec3 tangent = normalize(randomVec - normal * dot(randomVec, normal));
   vec3 bitangent = cross(normal, tangent);
   mat3 TBN = mat3(tangent, bitangent, normal);
   // iterate over the sample kernel and calculate occlusion factor
   float occlusion = 0.0;
   int n_hits = 0; // non-background
   for(int i = 0; i < n_ssao_kernel_samples; ++i) {
      // get sample position
      vec3 samplePos = TBN * samples[i]; // from tangent to view-space
      samplePos = fragPos + samplePos * 3.0 * radius;

      // project sample position (to sample texture) (to get position on screen/texture)
      vec4 offset = vec4(samplePos, 1.0);
      offset = projection * offset; // from view to clip-space
      offset.xyz /= offset.w; // perspective divide
      offset.xyz = offset.xyz * 0.5 + 0.5; // transform to range 0.0 - 1.0
        
      // get sample depth
      float sampleDepth = texture(gPosition, offset.xy).z; // get depth value of kernel sample
      if (sampleDepth == 1.0) {
         // n_hits++;
      } else {
         // // range check & accumulate
         float rangeCheck = smoothstep(0.0, 1.0, radius / abs(fragPos.z - sampleDepth));
         occlusion += (sampleDepth >= (samplePos.z + bias) ? 1.0 : 0.0) * rangeCheck;           
         n_hits++;
      }
   }

   if (n_hits == 0) n_hits = 1;
   float occ = 1.0 - (occlusion / float(n_hits));

   FragColor = occ;

}
