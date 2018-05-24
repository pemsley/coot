
#include "Python.h"
#include "graphics-info.h"

// static 
gint
graphics_info_t::drag_refine_refine_intermediate_atoms() {

   int retprog = -1;
#ifdef HAVE_GSL

   graphics_info_t g;

   if (! g.last_restraints) {
      std::cout << "Null last restraints " << std::endl;
      return retprog;
   }

   // While the results of the refinement are a conventional result
   // (unrefined), let's continue.  However, there are return values
   // that we will stop refining and remove the idle function is on a
   // GSL_ENOPROG(RESS) and GSL_SUCCESS.... actually, we will remove
   // it on anything other than a GSL_CONTINUE
   //

   // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
   // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
   // coot::restraint_usage_Flags flags = coot::BONDS_AND_PLANES;
   // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_PARALLEL_PLANES;
   coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;

   if (do_torsion_restraints) {
      if (use_only_extra_torsion_restraints_for_torsions_flag) { 
	 // flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
	 flags = coot::TYPICAL_RESTRAINTS;
      } else {
	 flags = coot::TYPICAL_RESTRAINTS_WITH_TORSIONS;
      }
   }

   if (do_rama_restraints)
      // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
      // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_RAMA_AND_PARALLEL_PLANES;
      flags = coot::ALL_RESTRAINTS;
   
   if (do_torsion_restraints && do_rama_restraints) {

      // Do we really need this fine control?
      
//       if (use_only_extra_torsion_restraints_for_torsions_flag) { 
// 	 // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
// 	 // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_RAMA_AND_PARALLEL_PLANES;
//       } else {
// 	 // This changes the function to using torsions (for non-peptide torsions)
// 	 // flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
// 	 flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_RAMA_AND_PARALLEL_PLANES;
//       }

   }

   if (g.auto_clear_atom_pull_restraint_flag) {
      // returns true when the restraint was turned off.
      std::vector<coot::atom_spec_t> specs_for_removed_restraints =
	 g.last_restraints->turn_off_atom_pull_restraints_when_close_to_target_position();
      if (specs_for_removed_restraints.size()) {
	 atom_pulls_off(specs_for_removed_restraints);
	 g.clear_atom_pull_restraints(specs_for_removed_restraints, true);
      }
   }

   // print_initial_chi_squareds_flag is 1 the first time then we turn it off.
   int steps_per_frame = dragged_refinement_steps_per_frame;
   if (! g.last_restraints->include_map_terms())
      steps_per_frame *= 6;

   // flags = coot::BONDS_AND_NON_BONDED;
   // flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS; // seems OK in parallel
   // flags = coot::BONDS_ANGLES_TORSIONS_NON_BONDED_CHIRALS_AND_TRANS_PEPTIDE_RESTRAINTS;
   // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_GEMAN_MCCLURE_DISTANCES; // crashes
   // flags = coot::TYPICAL_RESTRAINTS; // crashes
   // flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS; // seems OK in parallel

   // It is inconvenient that we can't do this:
   // flags += coot::restraint_usage_Flags(coot::JUST_RAMAS);
   // so turn restraint_usage_Flags into a class

   // so, 1087 is crashy,  59 is fine

   if (false)
      std::cout << "debug:: in drag_refine_refine_intermediate_atoms() calling minimize() with "
		<< flags << std::endl;

   g.last_restraints->set_lennard_jones_epsilon(graphics_info_t::lennard_jones_epsilon);

   graphics_info_t::saved_dragged_refinement_results =
      g.last_restraints->minimize(flags, steps_per_frame, print_initial_chi_squareds_flag);

   retprog = graphics_info_t::saved_dragged_refinement_results.progress;
   print_initial_chi_squareds_flag = 0;
   int do_disulphide_flag = 0;
   int draw_hydrogens_flag = 0;
   if (molecules[imol_moving_atoms].draw_hydrogens())
      draw_hydrogens_flag = 1;
   bool do_rama_markup = graphics_info_t::do_intermediate_atoms_rama_markup;
   bool do_rota_markup = graphics_info_t::do_intermediate_atoms_rota_markup;

   // wrap the filling of the rotamer probability tables
   //
   coot::rotamer_probability_tables *tables_pointer = NULL;

   if (do_rota_markup) {
      if (! rot_prob_tables.tried_and_failed()) {
	 if (rot_prob_tables.is_well_formatted()) {
	    tables_pointer = &rot_prob_tables;
	 } else {
	    rot_prob_tables.fill_tables();
	    if (rot_prob_tables.is_well_formatted()) {
	       tables_pointer = &rot_prob_tables;
	    }
	 }
      } else {
	 do_rota_markup = false;
      }
   }

   Bond_lines_container bonds(*g.moving_atoms_asc, imol_moving_atoms,
			      do_disulphide_flag, draw_hydrogens_flag,
			      do_rama_markup, do_rota_markup,
			      tables_pointer
			      );

   g.regularize_object_bonds_box.clear_up();
   g.regularize_object_bonds_box = bonds.make_graphical_bonds(g.ramachandrans_container,
							      do_rama_markup, do_rota_markup);

   // debug
   // coot::geometry_distortion_info_container_t gdic = last_restraints.geometric_distortions(flags);

   char *env = getenv("COOT_DEBUG_REFINEMENT");
   if (env)
      g.tabulate_geometric_distortions(*last_restraints);

   // Update the Accept/Reject Dialog if it exists (and it should do,
   // if we are doing dragged refinement).
   if (accept_reject_dialog) {
      if (saved_dragged_refinement_results.lights.size() > 0) {

	 if (false)
	    std::cout << "debug:: here in drag_refine_refine_intermediate_atoms() calling "
		      << "update_accept_reject_dialog_with_results() with rr.info \""
		      << saved_dragged_refinement_results.info_text << "\"" << std::endl;
	 update_accept_reject_dialog_with_results(accept_reject_dialog,
						  coot::CHI_SQUAREDS,
						  saved_dragged_refinement_results);
      }
   }

   if (true)
      if (do_coot_probe_dots_during_refine_flag)
	 g.do_interactive_coot_probe();


   g.run_post_intermediate_atoms_moved_hook_maybe();
   
#endif // HAVE_GSL

   return retprog;
}


// return true if flip moving_atoms_asc was found
bool graphics_info_t::pepflip_intermediate_atoms() {

   bool status = false;
   if (moving_atoms_asc->mol) {
      status = true;

      mmdb::Atom *at_close = NULL;
      float min_dist_sqrd = 4.0;

      coot::Cartesian pt(graphics_info_t::RotationCentre());

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) { 
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	 coot::Cartesian atom_pos(at->x, at->y, at->z);
	 coot::Cartesian diff = atom_pos - pt;
	 if (diff.amplitude_squared() < min_dist_sqrd) {
	    min_dist_sqrd = diff.amplitude_squared();
	    at_close = at;
	 }
      }

      if (at_close) {

	 mmdb::Residue *res_this = at_close->residue;
	 std::string atom_name = at_close->name;

	 // if N is at active atom then we want prev,this
	 // otherwise we want this,next
	 //
	 mmdb::Residue *res_1 = NULL;
	 mmdb::Residue *res_2 = NULL;
	 if (atom_name == " N  ") {
	    res_1 = moving_atoms_asc->get_previous(res_this);
	    res_2 = res_this;
	 } else {
	    res_1 = res_this;
	    res_2 = moving_atoms_asc->get_next(res_this);
	 }

	 if (res_1 && res_2) {
	    mmdb::Atom *at_1_ca = res_1->GetAtom(" CA ");
	    mmdb::Atom *at_1_c  = res_1->GetAtom(" C  ");
	    mmdb::Atom *at_1_o  = res_1->GetAtom(" O  ");
	    mmdb::Atom *at_2_ca = res_2->GetAtom(" CA ");
	    mmdb::Atom *at_2_n  = res_2->GetAtom(" N  ");
	    mmdb::Atom *at_2_h  = res_2->GetAtom(" H  ");

	    if (at_1_ca && at_2_ca) {
	       clipper::Coord_orth base(at_1_ca->x, at_1_ca->y, at_1_ca->z);
	       clipper::Coord_orth  top(at_2_ca->x, at_2_ca->y, at_2_ca->z);
	       clipper::Coord_orth dir = top - base;
	       coot::util::rotate_atom_about(dir, base, M_PI, at_1_c);
	       coot::util::rotate_atom_about(dir, base, M_PI, at_1_o);
	       coot::util::rotate_atom_about(dir, base, M_PI, at_2_n);
	       coot::util::rotate_atom_about(dir, base, M_PI, at_2_h); // does null check

	       drag_refine_refine_intermediate_atoms();

	    }
	 }
      }
   }
   graphics_draw();
   return status;
}

// return true if the isomerisation was made
//
bool
graphics_info_t::cis_trans_conversion_intermediate_atoms() {

   bool state = false;
   if (moving_atoms_asc->mol) {

      mmdb::Atom *at_close = NULL;
      float min_dist_sqrd = 4.0;

      coot::Cartesian pt(graphics_info_t::RotationCentre());

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	 coot::Cartesian atom_pos(at->x, at->y, at->z);
	 coot::Cartesian diff = atom_pos - pt;
	 if (diff.amplitude_squared() < min_dist_sqrd) {
	    min_dist_sqrd = diff.amplitude_squared();
	    at_close = at;
	 }
      }

      if (at_close) {

	 mmdb::Residue *res_this = at_close->residue;
	 std::string atom_name = at_close->name;

	 // if N is at active atom then we want prev,this
	 // otherwise we want this,next
	 //
	 mmdb::Residue *res_1 = NULL;
	 mmdb::Residue *res_2 = NULL;
	 if (atom_name == " N  ") {
	    res_1 = moving_atoms_asc->get_previous(res_this);
	    res_2 = res_this;
	 } else {
	    res_1 = res_this;
	    res_2 = moving_atoms_asc->get_next(res_this);
	 }

	 if (res_1 && res_2) {
	    mmdb::Atom *at_1_ca = res_1->GetAtom(" CA ");
	    mmdb::Atom *at_1_c  = res_1->GetAtom(" C  ");
	    mmdb::Atom *at_1_o  = res_1->GetAtom(" O  ");
	    mmdb::Atom *at_2_ca = res_2->GetAtom(" CA ");
	    mmdb::Atom *at_2_n  = res_2->GetAtom(" N  ");
	    mmdb::Atom *at_2_h  = res_2->GetAtom(" H  ");

	    bool is_N_flag = false;
	    if (atom_name == " N  ")
	       is_N_flag = true;

	    if (at_1_ca && at_2_ca) {

	       mmdb::Manager *standard_residues_mol = standard_residues_asc.mol;
	       mmdb::Manager *mol = moving_atoms_asc->mol;
	       int s = coot::util::cis_trans_conversion(at_close, is_N_flag, mol, standard_residues_mol);
	       state = true; // we found intermedate atoms and tried to convert them

	       // the atoms have moved, we need to continue the refinement.

	       drag_refine_refine_intermediate_atoms();

	       graphics_draw(); // needed

	    }
	 }
      }
   }
   return state;
}


#ifdef USE_PYTHON
void graphics_info_t::register_post_intermediate_atoms_moved_hook(PyObject *function) {

   std::cout << "::::::::::: set post_intermediate_atoms_moved_hook to " << function << std::endl;
   post_intermediate_atoms_moved_hook = function;

}
#endif



#ifdef USE_PYTHON
void graphics_info_t::run_post_intermediate_atoms_moved_hook_maybe() {

   if (post_intermediate_atoms_moved_hook) {

      graphics_info_t g;
      PyObject *o = g.get_intermediate_atoms_bonds_representation();

      if (PyBool_Check(o)) {
	 // no useful representation, (must be False)
      } else {
	 PyObject *py_main = PyImport_AddModule("__main__");
	 if (py_main) {
	    // wraps PyDict_SetItemString (String because the key is a string (the variable name))
	    PyModule_AddObject(py_main, "intermediate_atoms_representation_internal", o);
	 }
      }
   }
}
#endif
