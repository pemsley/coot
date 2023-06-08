
#include "graphics-info.h"

void
graphics_info_t::check_if_in_range_defines() {

   bool iaof = false; // use intermediates atom for picking?
   pick_info naii = atom_pick_gtk3(iaof);

   // rotamers are done by picking the residue close to the centre of the screen
   // (i.e. on button press, not atom pick)
   //
   // check_if_in_rotamer_define_gtk4(naii);


   // ------------------------- distance ---------------------------

   if (in_distance_define) {

      pick_info nearest_atom_index_info;
      nearest_atom_index_info = atom_pick_gtk3(false);

      if ( nearest_atom_index_info.success == GL_TRUE ) {

         int im = nearest_atom_index_info.imol;
         std::cout << "geometry: on molecule number: " << im << std::endl;
         // some visual feedback, label the atom:
         molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);

         if (in_distance_define == 1) {
            in_distance_define = 2; // flag for next atom pick
            geometry_atom_index_1 = nearest_atom_index_info.atom_index;
            geometry_atom_index_1_mol_no = nearest_atom_index_info.imol;
            mmdb::Atom *atom1 = molecules[im].atom_sel.atom_selection[geometry_atom_index_1];
            distance_pos_1 = coot::Cartesian(atom1->x, atom1->y, atom1->z);
            std::cout << "click on a second atom" << std::endl;
            graphics_draw();
         } else {

            // in_distance_define == 2
            geometry_atom_index_2 = nearest_atom_index_info.atom_index;
            geometry_atom_index_2_mol_no = nearest_atom_index_info.imol;

            mmdb::Atom *atom2 = molecules[im].atom_sel.atom_selection[geometry_atom_index_2];
            coot::Cartesian pos2 = coot::Cartesian(atom2->x, atom2->y, atom2->z);

            // 20190104-PE Why were we using the symmetry function?
            //             display_geometry_distance_symm(geometry_atom_index_1_mol_no, distance_pos_1,
            //                                            geometry_atom_index_2_mol_no, pos2);

            add_measure_distance(distance_pos_1, pos2); // calls graphics_draw()

            unset_geometry_dialog_distance_togglebutton();
            in_distance_define = 0;  // clear flag
            pick_pending_flag = 0;
            normal_cursor();
         }
      }
   }

   // ------------------------- angle ---------------------------

   if (in_angle_define) {
      // We need a cleaner way to know if this was an atom pick or a symm atom pick.
      // Let's sort it out in the beginning:
      short int picked = 0;
      pick_info nearest_atom_index_info = atom_pick_gtk3(false);
      coot::Cartesian pos;

      if (nearest_atom_index_info.success == GL_TRUE) {
         picked = 1;
         int im = nearest_atom_index_info.imol;
         mmdb::Atom *atom = molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index];
         molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
         pos = coot::Cartesian(atom->x, atom->y, atom->z);
      } else {
         coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick();
         if (symm_nearest_atom_index_info.success == GL_TRUE) {
            picked = 1;
            int im = symm_nearest_atom_index_info.imol;
            // some visual feedback, label the atom:
            molecules[im].add_atom_to_labelled_symm_atom_list(symm_nearest_atom_index_info.atom_index,
                                                              symm_nearest_atom_index_info.symm_trans,
                                                              symm_nearest_atom_index_info.pre_shift_to_origin);
            pos = symm_nearest_atom_index_info.hybrid_atom.pos;
         }
      }

      if (picked) {

         if (in_angle_define == 1) {
            in_angle_define = 2; // flag for next atom pick
            angle_tor_pos_1 = pos;
            graphics_draw();

         } else {

            if (in_angle_define == 2) {
               in_angle_define = 3; // flag for next atom pick
               angle_tor_pos_2 = pos;
               graphics_draw();

            } else {
               // in_angle_define == 3
               angle_tor_pos_3 = pos;
               graphics_draw();

               add_measure_angle(); // uses class members that we have just set

               in_angle_define = 0;  // clear flag
               pick_pending_flag = 0;
               normal_cursor();
               unset_geometry_dialog_angle_togglebutton();
            }
         }
         graphics_draw();
      }
   }


   // ------------------------- torsion ---------------------------

   //
   if (in_torsion_define) {

      // We need a cleaner way to know if this was an atom pick or a symm atom pick.
      // Let's sort it out in the beginning:
      bool picked = 0;
      pick_info nearest_atom_index_info = atom_pick_gtk3(false);
      mmdb::Atom *atom = 0;
      coot::Cartesian pos;

      if (nearest_atom_index_info.success == GL_TRUE) {
         picked = true;
         int im = nearest_atom_index_info.imol;
         atom = molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index];
         molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
         pos = coot::Cartesian(atom->x, atom->y, atom->z);
      } else {
         coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick();
         if (symm_nearest_atom_index_info.success == GL_TRUE) {
            picked = 1;
            int im = symm_nearest_atom_index_info.imol;
            // some visual feedback, label the atom:
            molecules[im].add_atom_to_labelled_symm_atom_list(symm_nearest_atom_index_info.atom_index,
                                                              symm_nearest_atom_index_info.symm_trans,
                                                              symm_nearest_atom_index_info.pre_shift_to_origin);
            pos = symm_nearest_atom_index_info.hybrid_atom.pos;
         }
      }

      if (picked) {

         if (in_torsion_define == 1) {
            angle_tor_pos_1 = pos;
            in_torsion_define = 2; // flag for next atom pick
            graphics_draw();
         } else {
            if (in_torsion_define == 2) {
               angle_tor_pos_2 = pos;
               in_torsion_define = 3; // flag for next atom pick
               graphics_draw();
            } else {
               if (in_torsion_define == 3) {
                  angle_tor_pos_3 = pos;
                  in_torsion_define = 4; // flag for next atom pick
                  graphics_draw();
               } else {
                  // in_torsion_define == 4
                  angle_tor_pos_4 = pos;
                  display_geometry_torsion(); // does a draw
                  in_torsion_define = 0; // clear up.
                  pick_pending_flag = 0;
                  normal_cursor();
                  unset_geometry_dialog_torsion_togglebutton();
               }
            }
         }
      }
   }


}

#if 0 // 20230526-PE - it doesn't work like this

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

#endif
