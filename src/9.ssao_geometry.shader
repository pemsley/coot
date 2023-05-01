
#shader vertex

#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;

out vec3 frag_pos_transfer;
out vec3 normal_transfer;

uniform bool invertedNormals;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    vec4 viewPos = view * model * vec4(aPos, 1.0);
    frag_pos_transfer = viewPos.xyz;
    
    mat3 normalMatrix = transpose(inverse(mat3(view * model)));
    normal_transfer = normalize(normalMatrix * (invertedNormals ? -aNormal : aNormal));

    gl_Position = projection * viewPos;
}


#shader fragment

#version 330 core

// actually, only gPosition and gNormal are used.
// so scrap gAlbedo
layout (location = 0) out vec3 gPosition;
layout (location = 1) out vec3 gNormal;

in vec3 frag_pos_transfer;
in vec3 normal_transfer;

void main() {
   // store the fragment position vector in the first gbuffer texture
   gPosition = frag_pos_transfer;
   // also store the per-fragment normals into the gbuffer
   gNormal = normal_transfer;

}
