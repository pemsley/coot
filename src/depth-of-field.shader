
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

uniform sampler2D screenTexture1;
uniform sampler2D screenTexture2;
uniform sampler2D screenDepth;

layout(location = 0) out vec4 out_colour;

void main() {

   vec3 t1 = texture(screenTexture1, TexCoords).rgb;
   vec3 t2 = texture(screenTexture2, TexCoords).rgb;
   float d = texture(screenDepth, TexCoords).r;

   // out_colour = vec4(result, 1.0);

   float mid = 0.5; // the real rotation centre changes depending on
                    // the clipping planes. mid should be calculated
                    // and passed as a uniform
   float dd = d - mid;
   float add = abs(dd * dd);

   float mf = abs(16.4 * (dd) * dd);
   mf = clamp(mf, 0.0, 1.0);

   out_colour = vec4(mix(t2, t1, mf), 1.0);

   // out_colour = vec4(add, add, add, 1.0);

}

