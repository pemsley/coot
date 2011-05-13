
// I don't know how to handle this yet.
//
// These functions are needed by the canvas item callbacks
// (e.g. residue circles), yet they are in src (and we don't want to
// include c-interface.h and all its dependencies here).
//
// 
// so for current convenience, I simply repeat the functions declarations here.

extern "C" {
   
   /* \brief display/undisplay the given additional representation  */
   void set_show_additional_representation(int imol, int representation_number, int on_off_flag);
   
   /* \brief display/undisplay all the additional representations for the given molecule  */
   void set_show_all_additional_representations(int imol, int on_off_flag);

   void all_additional_representations_off_except(int imol, int representation_number,
						  short int ball_and_sticks_off_too_flag);


}

// c++ functions
// 
void orient_view(int imol,
		 const coot::residue_spec_t &central_residue_spec,
		 const coot::residue_spec_t &neighbour_residue_spec);

void add_animated_ligand_interation(int imol,
				    const coot::atom_spec_t &spec_1,
				    const coot::atom_spec_t &spec_2
				    // type etc.
				    );

