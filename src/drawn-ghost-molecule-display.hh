
#include "api/ghost-molecule-display.hh"
#include "stereo-eye.hh"
#include "Shader.hh"
#include "Mesh.hh"

class drawn_ghost_molecule_display_t : public coot::ghost_molecule_display_t {

public:
   drawn_ghost_molecule_display_t() : mesh("ghost-name-here-A") {}
   drawn_ghost_molecule_display_t(const clipper::RTop_orth &rtop_in,
                                  int SelHnd_in,
                                  const std::string &name_in) :
      coot::ghost_molecule_display_t(rtop_in, SelHnd_in, name_in),
      mesh("ghost-name-here-B") {}
   drawn_ghost_molecule_display_t(const coot::ghost_molecule_display_t &g) :
      ghost_molecule_display_t(g), mesh("ghost-name-here-C") {}

   Mesh mesh;

   void draw(Shader *shader,
             stereo_eye_t eye,
             const glm::mat4 &mvp,
             const glm::mat4 &view_rotation_matrix,
             const std::map<unsigned int, lights_info_t> &lights,
             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
             const glm::vec4 &background_colour);
   // override virtual function
   void update_bonds(mmdb::Manager *mol);
};
