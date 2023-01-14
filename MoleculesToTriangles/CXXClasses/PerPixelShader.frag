varying vec3 normal, eyeVec;
varying vec3 lightDir0, lightDir1, lightDir2;
varying vec4 diffuse0, diffuse1, diffuse2;

void main (void)
{
	vec4 final_color = gl_FrontMaterial.emission + (gl_FrontLightModelProduct.sceneColor * gl_FrontMaterial.ambient);
    
	float dist0 = length(lightDir0);
	float att0 = 1.0 / (gl_LightSource[0].constantAttenuation +
                        gl_LightSource[0].linearAttenuation * dist0 +
                        gl_LightSource[0].quadraticAttenuation * dist0 * dist0);
	final_color += att0*(gl_LightSource[0].ambient * gl_FrontMaterial.ambient);
    
	vec3 N = normalize(normal);
	vec3 L = normalize(lightDir0);
	
	float lambertTerm = dot(N,L);
	
	if(lambertTerm > 0.0)
	{
		final_color += att0*diffuse0 * 
        lambertTerm;	
		
		vec3 E = normalize(eyeVec);
		vec3 R = reflect(-L, N);
		float specular = pow( max(dot(R, E), 0.0), 
                             gl_FrontMaterial.shininess );
		final_color += att0 * gl_LightSource[0].specular * 
        gl_FrontMaterial.specular * 
        specular;	
	}
    
	float dist1 = length(lightDir1);
	float att1 = 1.0 / (gl_LightSource[1].constantAttenuation +
                        gl_LightSource[1].linearAttenuation * dist1 +
                        gl_LightSource[1].quadraticAttenuation * dist1 * dist1);
    final_color += att1*(gl_LightSource[1].ambient * gl_FrontMaterial.ambient);
    
	L = normalize(lightDir1);
	
	lambertTerm = dot(N,L);
	
	if(lambertTerm > 0.0)
	{
		final_color += att1*diffuse1 * 
        lambertTerm;	
		
		vec3 E = normalize(eyeVec);
		vec3 R = reflect(-L, N);
		float specular = pow( max(dot(R, E), 0.0), 
                             gl_FrontMaterial.shininess );
		final_color += att1 * gl_LightSource[1].specular * 
        gl_FrontMaterial.specular * 
        specular;	
	}
    
	float dist2 = length(lightDir2);
	float att2 = 1.0 / (gl_LightSource[2].constantAttenuation +
                        gl_LightSource[2].linearAttenuation * dist2 +
                        gl_LightSource[2].quadraticAttenuation * dist2 * dist2);
    final_color += att2*(gl_LightSource[2].ambient * gl_FrontMaterial.ambient);
    
	L = normalize(lightDir2);
	
	lambertTerm = dot(N,L);
	
	if(lambertTerm > 0.0)
	{
		final_color += att2*diffuse2 * 
        lambertTerm;	
		
		vec3 E = normalize(eyeVec);
		vec3 R = reflect(-L, N);
		float specular = pow( max(dot(R, E), 0.0), 
                             gl_FrontMaterial.shininess );
		final_color += att2 * gl_LightSource[2].specular * 
        gl_FrontMaterial.specular * 
        specular;	
	}
    
    float z = gl_FragCoord.z / gl_FragCoord.w;
    //float fog = (gl_Fog.end - z) * gl_Fog.scale;

    float fog = (gl_Fog.end - eyeVec.z) * gl_Fog.scale;

    float fogFactor = clamp(fog, 0., 1.);
    
	gl_FragColor = mix(gl_Fog.color, final_color, fogFactor);			
	//gl_FragColor = final_color;
}
