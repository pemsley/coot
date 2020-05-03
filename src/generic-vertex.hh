
class generic_vertex {
public:
   glm::mat3 model_rotation_matrix; // orientation
   glm::vec3 model_translation; // the coordinates of the first atom of the bond
   glm::vec3 pos;
   glm::vec3 normal; // normalized when set
   glm::vec4 colour;
};

