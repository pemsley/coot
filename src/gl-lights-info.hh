
#ifndef GL_LIGHTS_INFO_T
#define GL_LIGHTS_INFO_T

class gl_lights_info_t {
public:
   gl_lights_info_t() {
      is_on = true;
      directional = true;
      position = glm::vec4(0,0,1,1);
      ambient  = glm::vec4(1,1,1,1);
      diffuse  = glm::vec4(1,1,1,1);
      specular = glm::vec4(1,1,1,1);
      shininess = 10;
      constant_attenuation = 1.0;
      linear_attenuation = 1.0;
      quadratic_attenuation = 1.0;
   }
   bool is_on;
   bool directional;
   glm::vec4 position;
   glm::vec4 ambient;
   glm::vec4 diffuse;
   glm::vec4 specular;
   float shininess;
   float constant_attenuation;
   float linear_attenuation;
   float quadratic_attenuation;
   // also spot-light parameters?
   // later.

};

#endif // GL_LIGHTS_INFO_T
