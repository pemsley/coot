
uniform float Time;
varying highp vec3 normal, eyeVec;
varying highp vec3 lightDir0;
varying highp vec4 diffuse0;
varying highp vec3 lightDir1;
varying highp vec4 diffuse1;
varying float fogFactor;

uniform mat4 mygl_ModelViewMatrix;
uniform mat4 mygl_ProjectionMatrix;
uniform mat3 mygl_NormalMatrix;

attribute vec4 mygl_Vertex;
attribute vec4 mygl_Color;
attribute vec3 mygl_Normal;

uniform int mygl_ColorMaterial;
uniform mediump int mygl_UseLight0;
uniform mediump int mygl_UseLight1;

struct mygl_LightSourceParameters {
    vec4 ambient;
    vec4 diffuse;
    vec4 specular;
    vec4 position;
    vec4 halfVector;
    vec3 spotDirection;
    float spotExponent;
    float spotCutoff;
    float spotCosCutoff;
    float constantAttenuation;
    float linearAttenuation;
    float quadraticAttenuation;
};

uniform mygl_LightSourceParameters mygl_LightSource[2];

uniform int mygl_UseColorArray;

struct mygl_MaterialParameters {
    vec4 emission;
    vec4 ambient;
    vec4 diffuse;
    vec4 specular;
    float shininess;
};
uniform mygl_MaterialParameters mygl_FrontMaterial;
uniform mygl_MaterialParameters mygl_BackMaterial;

void main()
{
    vec4 modified_gl_Vertex = mygl_Vertex;
    
    vec3 vVertex = vec3(gl_ModelViewMatrix * modified_gl_Vertex );
    
    eyeVec = -vVertex;

    vec4 frontColor = mygl_FrontMaterial.diffuse;
    if (mygl_UseColorArray == 1){
        frontColor = mygl_Color;
    }
    
    normal = gl_NormalMatrix * normalize(mygl_Normal);
    gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * modified_gl_Vertex;
    // lightDir0 = (gl_LightSource[0].position).xyz - vVertex;
    lightDir0 = (gl_LightSource[0].position).xyz;
    diffuse0 = (frontColor*gl_LightSource[0].diffuse) / 255.;
    // lightDir1 = (gl_LightSource[1].position).xyz - vVertex;
    lightDir1 = (gl_LightSource[1].position).xyz;
    diffuse1 = (frontColor*gl_LightSource[1].diffuse) / 255.;
    
    //gl_FogFragCoord = length(vVertex);
    fogFactor = clamp((gl_Fog.end - vVertex[2]) / (gl_Fog.end - gl_Fog.start), 0., 1.);
    
}
