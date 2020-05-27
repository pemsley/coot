

#include <glm/ext.hpp>

class Material {
public:
   glm::vec4 ambient;
   glm::vec4 diffuse;
   glm::vec4 specular;
   float shininess;
   Material(const glm::vec4 &ambient,
            const glm::vec4 &diffuse,
            const glm::vec4 &specular,
            const float &s) : ambient(ambient), diffuse(diffuse), specular(specular), shininess(s) {}
   Material() {}
};
