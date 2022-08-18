
varying highp vec3 normal, eyeVec;
varying highp vec3 lightDir0;
varying highp vec4 diffuse0;
varying highp vec3 lightDir1;
varying highp vec4 diffuse1;
varying float fogFactor;

uniform int mygl_UseLight0;
uniform int mygl_UseLight1;

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

struct mygl_MaterialParameters {
	vec4 emission;   
	vec4 ambient;    
	vec4 diffuse;    
	vec4 specular;   
	float shininess; 
};

uniform mygl_MaterialParameters mygl_FrontMaterial;
uniform mygl_MaterialParameters mygl_BackMaterial;

struct mygl_FrontLightModelProductParameters {
	vec4 sceneColor;
};

uniform mygl_FrontLightModelProductParameters mygl_FrontLightModelProduct;

void main (void)
{
    vec4 final_color = mygl_FrontMaterial.emission + 
                       (mygl_FrontLightModelProduct.sceneColor * mygl_FrontMaterial.ambient);
    highp vec3 N = normalize(normal);

    if (mygl_UseLight0 == 1)
    {
        float dist0 = 1.0*length(lightDir0);
        float att0 = 1.0 / (gl_LightSource[0].constantAttenuation +
                            gl_LightSource[0].linearAttenuation * dist0 +
                            gl_LightSource[0].quadraticAttenuation * dist0 * dist0);
        final_color += (gl_LightSource[0].ambient * mygl_FrontMaterial.ambient)*att0;

        highp vec3 L = normalize(lightDir0);

        float lambertTerm = dot(N,L); // both normalized

        if (lambertTerm > 0.0) {
            // the surface normal and the light direction are approximately the same 
	    // so we should see a shiny surface.  This is mostly true for the triangles
	    // when vVertex is removed from lightDir0.

            final_color += att0*diffuse0*lambertTerm;

            highp vec3 E = normalize(eyeVec);
	    E = vec3(0,0,1);
            highp vec3 R = reflect(-L, N);
            highp float specular = pow( max(dot(R, E), 0.0), mygl_FrontMaterial.shininess );
            final_color += 0.4 * gl_LightSource[0].specular * (specular*att0);
        }
    }

    if (mygl_UseLight1 == 1)
    {
        float dist1 = 1.0*length(lightDir1);
        float att1 = 1.0 / (gl_LightSource[1].constantAttenuation +
                            gl_LightSource[1].linearAttenuation * dist1 +
                            gl_LightSource[1].quadraticAttenuation * dist1 * dist1);
        final_color += 0.4 * (gl_LightSource[1].ambient * mygl_FrontMaterial.ambient)*att1;
	// final_color = vec4(0.0, 1.0, 0.0, 0.0);

        highp vec3 L = normalize(lightDir1);

        float lambertTerm = dot(N,L);

        if (lambertTerm > 0.0) {
            final_color += att1*diffuse1*lambertTerm;	
		
            // highp vec3 E = normalize(eyeVec);
	    highp vec3 Elocal = vec3(0,1,0.3);
            highp vec3 E = normalize(Elocal);
            highp vec3 R = reflect(-L, N);
            highp float specular = pow( max(dot(R, E), 0.0), mygl_FrontMaterial.shininess);
            final_color += 0.4 * gl_LightSource[1].specular * (specular*att1);
	    // final_color = vec4(1.0, 0.0, 0.0, 0.0);
        }
    }

    gl_FragColor = mix(gl_Fog.color, 0.7 * final_color, fogFactor);
    //gl_FragColor = gl_Fog.color;
}
