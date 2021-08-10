
#ifndef COOT_MATERIAL_HH
#define COOT_MATERIAL_HH

#include <glm/ext.hpp>

class Material {
public:
   glm::vec4 ambient;
   glm::vec4 diffuse;
   glm::vec4 specular;
   bool do_specularity;
   float shininess;
   float specular_strength;
   Material(const glm::vec4 &ambient,
            const glm::vec4 &diffuse,
            const glm::vec4 &specular,
            const float &s) : ambient(ambient), diffuse(diffuse), specular(specular), shininess(s) {
      specular_strength = 1.0;
      do_specularity = false;
   }
   Material() {
      ambient  = glm::vec4(0.5, 0.5, 0.5, 1.0);
      diffuse  = glm::vec4(0.5, 0.5, 0.5, 1.0);
      specular = glm::vec4(0.5, 0.5, 0.5, 1.0);
      do_specularity = false;
      shininess = 66.0;
      specular_strength = 1.0;
   }
   void turn_specularity_on(bool state) {
      do_specularity = state;
   }
};


#endif // COOT_MATERIAL_HH
