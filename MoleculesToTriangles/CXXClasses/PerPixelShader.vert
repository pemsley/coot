uniform float Time;
uniform vec3 offset;
uniform float scaleIn;
uniform float scaleOut;
varying vec3 normal, eyeVec;
varying vec3 lightDir0, lightDir1, lightDir2;
varying vec4 diffuse0, diffuse1, diffuse2;

//noisemaker 
vec3 mynoise3(vec3 co)
{
	vec3 n1 = normalize(co);
	float noise1 = (fract(sin(dot(n1.xy ,vec2(12.9898,78.233))) * 43758.5453));
	float noise2 = (fract(sin(dot(n1.yz ,vec2(12.9898,78.233))) * 43758.5453));
	float noise3 = (fract(sin(dot(n1.zx ,vec2(12.9898,78.233))) * 43758.5453));
	return clamp(vec3(noise1,noise2,noise3),0.0,1.0);
}

void main()
{	
	normal = gl_NormalMatrix * gl_Normal;
    
    float timeFactor = 0.;
	//1.4 * ((1.+sin((Time/10.)*2.*3.141592654))/2.);
    vec4 modified_gl_Vertex = gl_Vertex + vec4((timeFactor * gl_Normal), 0.);
    
    //Here to add a wobbling, using the noise3 function 
    
	//vec3 vertex = gl_Vertex.xyz + mynoise3(gl_Vertex.xyz * Time) * scaleOut;
	//vec3 vertex = gl_Vertex.xyz + noise3(offset + gl_Vertex.xyz * scaleIn) * scaleOut;
	//vec4 modified_gl_Vertex = vec4(vertex, 1.);
    
    
	vec3 vVertex = vec3(gl_ModelViewMatrix * modified_gl_Vertex );
	eyeVec = -vVertex;
    
	lightDir0 = vec3(gl_LightSource[0].position.xyz - vVertex);
    diffuse0 = gl_Color*gl_LightSource[0].diffuse;
	lightDir1 = vec3(gl_LightSource[1].position.xyz - vVertex);
    diffuse1 = gl_Color*gl_LightSource[1].diffuse;
	lightDir2 = vec3(gl_LightSource[2].position.xyz - vVertex);
    diffuse2 = gl_Color*gl_LightSource[2].diffuse;
    
	gl_Position = (gl_ProjectionMatrix * gl_ModelViewMatrix * modified_gl_Vertex);
}
