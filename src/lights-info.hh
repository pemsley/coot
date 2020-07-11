

#ifndef LIGHTS_INFO_T_HH
#define LIGHTS_INFO_T_HH

#include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>

class lights_info_t {
public:
   lights_info_t() {
      is_on = true;
      directional = true;
      position = glm::vec4(0,0,1,1);
      direction = glm::vec3(0,0,1);
      ambient  = glm::vec4(1,1,1,1);
      diffuse  = glm::vec4(1,1,1,1);
      specular = glm::vec4(1,1,1,1);
      shininess = 10.0;
      constant_attenuation  = 1.0;
      linear_attenuation    = 1.0;
      quadratic_attenuation = 1.0;
   }
   bool is_on;
   bool directional;
   glm::vec4 position;
   glm::vec3 direction;
   glm::vec4 ambient;
   glm::vec4 diffuse;
   glm::vec4 specular;
   glm::vec3 spot_direction;
   glm::vec4 half_vector;
   float spot_exponent;
   float spot_cutoff;
   float speot_cos_cutoff;
   float shininess;
   float constant_attenuation;
   float linear_attenuation;
   float quadratic_attenuation;
   void rotate_light_direction(float angle_rad, const glm::vec3 &axis) {
      direction = glm::rotate(direction, angle_rad, axis);
   }

};

#endif // GL_LIGHTS_INFO_T_HH
