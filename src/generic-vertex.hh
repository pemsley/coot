
#ifndef GENERIC_VERTEX_HH
#define GENERIC_VERTEX_HH

#include <glm/glm.hpp>

// for standard objects at the origin - typically used in instancing
// 20230109-PE is this used? Having this here may be confusing.
//
class vn_vertex {
public:
   glm::vec3 pos;
   glm::vec3 normal; // normalized on input
   vn_vertex(const glm::vec3 &pos_in,
             const glm::vec3 &norm_in) :
             pos(pos_in), normal(norm_in) {}
   vn_vertex() {}
};

// simple vertex
class s_generic_vertex {
public:
   glm::vec3 pos;
   glm::vec3 normal; // normalized on input
   glm::vec4 color;  // make this "colour"
   s_generic_vertex(const glm::vec3 &pos_in,
                    const glm::vec3 &norm_in,
                    const glm::vec4 &col_in) : pos(pos_in), normal(norm_in), color(col_in) {}
   explicit s_generic_vertex(const vn_vertex &vn) :
      pos(vn.pos), normal(vn.normal), color(glm::vec4(0.5, 0.5, 0.5, 1.0)) {}
   s_generic_vertex() {}
};


// 20220222-PE
// same as below, but below as written first and I don't want to change the class name
// at the moment. Maybe later. This is for the simple-lines representation of the molecule.
// Actually the vector is not stored in Mesh. The vector gets clears up when the
// buffer construction function ends
class simple_atoms_line_vertex {
public:
   glm::vec3 pos;
   glm::vec4 colour;
   simple_atoms_line_vertex(const glm::vec3 &pos_in,
                            const glm::vec4 &col_in) : pos(pos_in), colour(col_in) {}
   simple_atoms_line_vertex() {}
};

class symmetry_atoms_line_vertex {
public:
   glm::vec3 pos;
   glm::vec4 colour;
   symmetry_atoms_line_vertex(const glm::vec3 &pos_in,
                              const glm::vec4 &col_in) : pos(pos_in), colour(col_in) {}
   symmetry_atoms_line_vertex() {}
};

// for instanced objects (those animated/moving) the position in molecular space
// and the colour is dictacted by the instancing matrices (and colours)
// We don't need to create them here.
class position_normal_vertex {
public:
   glm::vec3 pos; // based around the origin
   glm::vec3 normal; // normalized on input
   position_normal_vertex(const glm::vec3 &pos_in, const glm::vec3 &norm_in) : pos(pos_in), normal(norm_in) {}
   position_normal_vertex() {}
};


// moved into coot utils.

// class vertex_with_rotation_translation {
// public:
//    glm::mat3 model_rotation_matrix; // orientation
//    glm::vec3 model_translation; // the coordinates of the first atom of the bond
//    glm::vec3 pos;
//    glm::vec3 normal; // normalized when set
//    glm::vec4 colour;
//    vertex_with_rotation_translation(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &c) : pos(p), normal(n), colour(c) {}
//    vertex_with_rotation_translation(const s_generic_vertex &v, const glm::vec3 &atom_position, float scale) :
//       model_rotation_matrix(glm::mat3(1.0f)), model_translation(atom_position),
//       pos(v.pos * scale), normal(v.normal), colour(v.color) {}
//    vertex_with_rotation_translation() {}
// };



#endif

