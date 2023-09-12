
// --------------- 9.ssao_blur.shader ----------------

#shader vertex

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

#version 330 core
out float FragColor;

in vec2 TexCoords;

uniform sampler2D ssaoInput;

uniform int blur_size;

void main() {
    vec2 texelSize = 1.0 / vec2(textureSize(ssaoInput, 0));
    float result = 0.0;
    float n_blur = 0.0;
    for (int x = -blur_size; x <= blur_size; ++x) {
        for (int y = -blur_size; y <= blur_size; ++y) {
            vec2 offset = vec2(float(x), float(y)) * texelSize;
            float ts = texture(ssaoInput, TexCoords + offset).r;
            // This test removes the SSAO outline.
            // Perhaps use a weight w = 1 = ts ?
            if (ts != 1.0) { // don't blur with the background
               result += ts;
               n_blur += 1.0; // counting int to float
            }
        }
    }
    if (n_blur > 0)
       FragColor = result / n_blur;
    else
       FragColor = texture(ssaoInput, TexCoords).r;

    // FragColor = texture(ssaoInput, TexCoords).r;

}  
