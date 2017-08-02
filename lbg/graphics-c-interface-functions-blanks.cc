
// #include "coot-utils/coot-coord-utils.hh"

extern "C" {
   
   /* \brief display/undisplay the given additional representation  */
   void set_show_additional_representation(int imol, int representation_number, int on_off_flag) {}
   
   /* \brief display/undisplay all the additional representations for the given molecule  */
   void set_show_all_additional_representations(int imol, int on_off_flag) {}

   void all_additional_representations_off_except(int imol, int representation_number,
						  short int ball_and_sticks_off_too_flag) {}

}

namespace coot {
   class residue_spec_t;
}
namespace clipper {
   class Coord_orth;
}

void orient_view(int imol,
 		 const coot::residue_spec_t &central_residue_spec,
 		 const coot::residue_spec_t &neighbour_residue_spec) {}

void set_rotation_centre(const clipper::Coord_orth &x) {}
