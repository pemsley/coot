
#include "graphics-info.h"

void
graphics_info_t::check_if_in_range_defines() {

   bool iaof = false; // use intermediates atom for picking?
   pick_info naii = atom_pick_gtk3(iaof);

   // rotamers are done by picking the residue close to the centre of the screen
   // (i.e. on button press, not atom pick)
   //
   // check_if_in_rotamer_define_gtk4(naii);
   
}

void
graphics_info_t::check_if_in_rotamer_define_gtk4(const pick_info &naii) {

   if (in_rotamer_define) {
      if (naii.success == GL_TRUE) {
	 do_rotamers(naii.atom_index, naii.imol);
	 in_rotamer_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_rotamer_togglebutton");
	 if (false) {
	    if (moving_atoms_asc) {
	       std::cout << "debug moving atoms to moving-atoms.pdb" << std::endl;
	       moving_atoms_asc->mol->WritePDBASCII("moving-atoms.pdb");
	    } else {
	       std::cout << "debug no moving atoms object" << std::endl;
	    }
	 }
      }
   }
}
