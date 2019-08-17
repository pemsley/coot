
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
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D screenTexture;
uniform sampler2D screenDepth;

void main() {

    vec3 result = texture(screenTexture, TexCoords).rgb;
    ivec2 texcoord = ivec2(floor(gl_FragCoord.xy));
    vec4 depth = texelFetch(screenDepth, texcoord, 0);
    vec3 depth_v2 = texture(screenTexture, TexCoords).rgb;
    float z = depth.r;

    vec4 c = vec4(result, 1.0);

    if (z > 0.4) c = vec4(1.0, 0.0, 0.0, 1.0);

    FragColor = c;

}

