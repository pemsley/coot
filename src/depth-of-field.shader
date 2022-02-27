
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

uniform float focus_blur_z_depth;
uniform float focus_blur_strength;

uniform sampler2D screenTexture1;
uniform sampler2D screenTexture2;
uniform sampler2D screenDepth;

layout(location = 0) out vec4 out_colour;

void main() {

   vec3 t1 = texture(screenTexture2, TexCoords).rgb; // starting
   vec3 t2 = texture(screenTexture1, TexCoords).rgb; // blurred
   float d = texture(screenDepth, TexCoords).r;

   vec4 t1a = vec4(t1, 1.0); // make these dependent on d
   vec4 t2a = vec4(t2, 1.0);

   // out_colour = vec4(result, 1.0);

   // perspective depth can be calcuated in the draw() function
   // and passed here

   // https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/building-basic-perspective-projection-matrix
   // Remapping the Z-Coordinate

   float rc_z = focus_blur_z_depth; // the real rotation centre changes depending on
                      // the clipping planes. mid should be calculated
                      // and passed as a uniform
   float dd = d - rc_z;
   float mf = abs(focus_blur_strength * dd);
   mf = clamp(mf, 0.0, 1.0);

   out_colour = vec4(mix(t1a, t2a, mf));

   // out_colour = t1a;

   //    out_colour = vec4(vec3(d), 1.0);

   // if (d == 1.0) out_colour = vec4(1,0,1,1);

}

