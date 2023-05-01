

  

/* testing/debugging - fails - sigh. */
/*   tmp_list = g_list_nth (GTK_MENU_SHELL (render_optionmenu_1_menu)->children, 1); */
/*   if (tmp_list) { */
/*      printf("got a tmp list\n"); */
/*      child = tmp_list->data; */
/*      if (GTK_BIN (child)->child) { */
/* 	printf("got a  (child)->child\n"); */
/* 	if (GTK_MENU(menu)->old_active_menu_item) */
/* 	   gtk_widget_unref (GTK_MENU(menu)->old_active_menu_item); */
/* 	GTK_MENU(menu)->old_active_menu_item = child; */
/* 	gtk_widget_ref (GTK_MENU(menu)->old_active_menu_item); */
/*      } */
/*   } else {  */
/*      printf("failed to get a tmp list\n"); */
/*   }  */



/*! \brief toggle the display of map in molecule number imol

  @param imol is the molecule number
  @param imap is ignored */
int toggle_display_map(int imol, int imap); 

/*! \brief toggle the display of coordinates molecule number imol */
int toggle_display_mol(int imol); 

/*! \brief toggle the active state (clickable) of coordinates molecule
  number imol */
int toggle_active_mol (int imol); 

// widget (toggle button) call-backs
// 
int toggle_display_map(int imol, int imap) { 
   
   graphics_info_t g;
   int i = g.molecules[imol].toggle_display_map(imap); 
   graphics_draw();
   std::vector<std::string> command_strings;
   command_strings.push_back("toggle-display-map");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   command_strings.push_back(graphics_info_t::int_to_string(imap));
   add_to_history(command_strings);
   return i;
} 


int toggle_display_mol(int imol) { 

   graphics_info_t g; 
   
   int i = g.molecules[imol].toggle_display_mol();
   graphics_draw();
   std::vector<std::string> command_strings;
   command_strings.push_back("toggle-display-mol");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
   return i; 
}

// Return the new pickable? state.
// 
int toggle_active_mol(int imol) { 

   graphics_info_t g; 

   int i = g.molecules[imol].toggle_active_mol(); 

   return i; 

} 


// molecule class functions:
   int toggle_display_map(int i_map) { 
      
      // int i_map corresponds to the multi dimensionality of drawit_for_map
      
      drawit_for_map = 1 - drawit_for_map; 
      
      return drawit_for_map; 

   }

   int toggle_display_mol() { 
//       std::cout << "changing display status from " << drawit
// 		<< " to "<< 1 - drawit << std::endl;
      drawit = 1 - drawit; 
      return drawit; 
   }

   int toggle_active_mol() { 
//       std::cout << "changing pickable_atom_selection status from " 
// 		<< pickable_atom_selection
// 		<< " to "<< 1 - pickable_atom_selection << std::endl;
      pickable_atom_selection = 1 - pickable_atom_selection; 
      return pickable_atom_selection; 
   } 



         if (mutate_auto_fit_do_post_refine_flag) {
	    // Run refine zone with autoaccept, autorange on
	    // the "clicked" atom:
	    mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	    std::string chain_id = at->GetChainID();
	    refine_auto_range(naii.imol, chain_id.c_str(), at->GetSeqNum(),
			      at->altLoc);
	 }

mmdb::realtype u, v, w;
      // get cell_trans_for_symm
      coot::Cartesian centre_point = centre_of_molecule();
      atom_sel.mol->Orth2Frac(centre_point.x(), centre_point.y(), centre_point.z(), u, v, w);

      // doub ud = u
      int iu = int (rint(u));
      int iv = int (rint(v));
      int iw = int (rint(w));  // 0.6 -> 1 and -0.6 -> -1
      
      std::cout << "DEBUG:: on reading coords cell_trans_for_symm: "
		<< cell_trans_for_symm.us << " " 
		<< cell_trans_for_symm.vs << " " 
		<< cell_trans_for_symm.ws << std::endl;


void rotate_intermediate_atoms_round_vector(const clipper::Coord_orth &direction,
					       const clipper::Coord_orth &position,
					       const clipper::Coord_orth &origin_shift,
					       double angle);

   std::vector<std::string> command_strings;
   command_strings.push_back("");

   add_to_history(command_strings);


// Old style buttons for rotate/translate
// 
// rottrans_buttons class calls back this function on button pressed mouse motion
// 
void
graphics_info_t::rot_trans_obj(int x_diff, const std::string &button_label) { 

   Cartesian centre = unproject_xyz(0, 0, 0.5);
   Cartesian front  = unproject_xyz(0, 0, 0.0);
   Cartesian right  = unproject_xyz(1, 0, 0.5);
   Cartesian top    = unproject_xyz(0, 1, 0.5);

   Cartesian screen_x = (right - centre);
   Cartesian screen_y = (top   - centre);
   Cartesian screen_z = (front - centre);

   screen_x.unit_vector_yourself();
   screen_y.unit_vector_yourself();
   screen_z.unit_vector_yourself();

   float x_add = 0.0;
   float y_add = 0.0;
   float z_add = 0.0;

   if (button_label == "rotate_translate_obj_xtrans_button") { 
      x_add = screen_x.x() * x_diff * 0.002 * zoom;
      y_add = screen_x.y() * x_diff * 0.002 * zoom;
      z_add = screen_x.z() * x_diff * 0.002 * zoom;
   }
   if (button_label == "rotate_translate_obj_ytrans_button") { 
      x_add = screen_y.x() * x_diff * 0.002 * zoom;
      y_add = screen_y.y() * x_diff * 0.002 * zoom;
      z_add = screen_y.z() * x_diff * 0.002 * zoom;
   }
   if (button_label == "rotate_translate_obj_ztrans_button") {
      x_add = screen_z.x() * x_diff * 0.002 * zoom;
      y_add = screen_z.y() * x_diff * 0.002 * zoom;
      z_add = screen_z.z() * x_diff * 0.002 * zoom;
   }

   short int do_rotation = 0; 
   clipper::Coord_orth screen_vector; 

   if (button_label == "rotate_translate_obj_xrot_button") {
      do_rotation = 1;
      screen_vector = clipper::Coord_orth(screen_x.x(), 
					   screen_x.y(), 
					   screen_x.z());
   }
   if (button_label == "rotate_translate_obj_yrot_button") {
      do_rotation = 1;
      screen_vector = clipper::Coord_orth(screen_y.x(), 
					  screen_y.y(), 
					  screen_y.z());
   }
   if (button_label == "rotate_translate_obj_zrot_button") {
      do_rotation = 1;
      screen_vector = clipper::Coord_orth(screen_z.x(), 
					  screen_z.y(), 
					  screen_z.z());
   }
   
   if (do_rotation) { 
      if (rot_trans_rotation_origin_atom) {
	 
	 // int indx = rot_trans_atom_index_rotation_origin_atom;
	 mmdb::Atom *rot_centre = rot_trans_rotation_origin_atom;
	 clipper::Coord_orth rotation_centre(rot_centre->x, 
					     rot_centre->y, 
					     rot_centre->z);
	 
	 for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	    clipper::Coord_orth co(moving_atoms_asc->atom_selection[i]->x,
				   moving_atoms_asc->atom_selection[i]->y,
				   moving_atoms_asc->atom_selection[i]->z);
	    clipper::Coord_orth new_pos = 
	       rotate_round_vector(screen_vector, co, rotation_centre, x_diff * 0.01);
	    moving_atoms_asc->atom_selection[i]->x = new_pos.x();
	    moving_atoms_asc->atom_selection[i]->y = new_pos.y();
	    moving_atoms_asc->atom_selection[i]->z = new_pos.z();
	 }
      } else { 
	 // this should never happen.
	 std::cout << "sorry - rotation atom not found" << std::endl;
      }
   }

   for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
      moving_atoms_asc->atom_selection[i]->x += x_add;
      moving_atoms_asc->atom_selection[i]->y += y_add;
      moving_atoms_asc->atom_selection[i]->z += z_add;
   }
   int do_disulphide_flag = 0;
   Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
   regularize_object_bonds_box.clear_up();
   regularize_object_bonds_box = bonds.make_graphical_bonds();
   gtk_widget_draw(glarea, NULL);
}
// Old unused code
// 
void 
autobuild_ca_on() { 

   graphics_info_t g; 

   g.autobuild_flag = 1;

   // we use break so that we get only one skeleton. 
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    if (g.molecules[imol].xskel_is_filled == 0) { 
	       
	       cout << "----------------------------------------" << endl; 
	       cout << "must Calculate the Map Skeleton first..." << endl; 
	       cout << "----------------------------------------" << endl; 
	       
	    } else {

	       float map_cutoff  = g.skeleton_level;
	       //
	       // save a pointer to this map btw.
	       // 
	       BuildCas bc(g.molecules[imol].xmap_list[imap], map_cutoff); 


	       // mark segments by connectivity
	       // 
	       int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan, 
							  g.molecules[imol].xmap_list[imap],
							  map_cutoff); 

	       cout << "INFO:: There were " << nsegments << " different segments" << endl; 
	       
	       // --------------------------------------------------------------
	       // --------------------------------------------------------------
	       cout << "cuckoo code:" << endl; 

	       bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
	       g.molecules[imol].update_clipper_skeleton();

	       // --------------------------------------------------------------
	       // --------------------------------------------------------------


	       // bc.depth_search_skeleton_testing_2(); 

	       // add branch points
	       //
	       cout << "INFO:: finding branch points..." << endl; 
	       //
	       vector<coot::Cartesian> branch_pts =
		  bc.find_branch_points(g.molecules[imol].xmap_list[imap],
					g.molecules[imol].xskel_cowtan,
					map_cutoff);
	       // 
	       cout << "INFO:: converting branch points to asc..." << endl; 
	       // 
	       atom_selection_container_t branch_pts_as_asc =
		  bc.convert_to_atoms(g.molecules[imol].xmap_list[imap],
				      branch_pts, "branch points");
 
	       // slight tangle here (internally to BuildCas,
	       // branch_points_symm_expanded is a vector<coot::Cartesian>, but
	       // we want an asc, so need to convert, so we need cell and
	       // symm (so we pass a const reference to the map). 
	       // 
	       atom_selection_container_t s_e_branch_pts_as_asc = 
		  bc.symmetry_expanded_branch_points(g.molecules[imol].xmap_list[imap]); 

	       cout << "INFO:: c-interface branch points conversion to atoms done!" << endl; 

	       cout << "INFO:: c-interface converting skeleton points to atoms..." << endl; 

// 	       asc_and_grids all_skels_pts_in_asu =
// 		  bc.all_skel_pts_in_asu(g.molecules[imol].xmap_list[imap],
// 					 g.molecules[imol].xskel_cowtan,
// 					 map_cutoff); // was 0.2

	       asc_and_grids all_skels_pts_in_asu =
		  bc.toplevel_skel_pts_in_asu(); // use internal segment_map

	       cout << "INFO:: c-interface expanding skeleton points by symmetry..." << endl; 

	       atom_selection_container_t big_ball =
		  bc.build_big_ball(g.molecules[imol].xmap_list[imap],
				    all_skels_pts_in_asu.asc, 
				    all_skels_pts_in_asu.grid_points); 

	       GraphicalSkel cowtan; 

	       int n_tips = cowtan.N_tips(g.molecules[imol].xmap_list[imap],
					  g.molecules[imol].xskel_cowtan,
					  map_cutoff);

 	       bc.interconnectedness(n_tips);

	       // now make that atom_selection_container for branch points:
	       // 

	       // Turn this back on when we have filled in molecule and map control.
	       // As it stood this code makes the molecule "active" (i.e. clickable)
	       // but the atoms were not displayed. 
	       // 
//  	       g.molecules[imol_new].initialize_coordinate_things_on_read_molecule(label);
//  	       g.molecules[imol_new].atom_sel = s_e_branch_pts_as_asc; 
//  	       g.molecules[imol_new].makebonds(0.1); // we don't want to join branch points
// 	       g.n_molecules++;

 	       std::string label = "branch points (symm expanded)";
	       asc_to_graphics(s_e_branch_pts_as_asc, label, ATOM_BONDS, 0.1); 

	       // debug_atom_selection_container(g.molecules[imol].atom_sel); 



	       // now make the atom_selection_container for the big ball
	       // viewable:
	       // 

// 	       g.molecules[imol_new+1].initialize_coordinate_things_on_read_molecule("big ball");
// 	       g.molecules[imol_new+1].atom_sel = big_ball; 
// 	       g.molecules[imol_new+1].make_ca_bonds(3.72, 3.85); //uses atom_sel
// 	       g.n_molecules++;
	       
	       // CA_BONDS !?
	       // 
	       asc_to_graphics(big_ball, "big ball", CA_BONDS, 3.72, 3.85); 


	       
	       // bc.ca_grow(g.molecules[imol].xmap_list[imap]); 
	       bc.ca_grow_recursive(); 

	       // show the cas
// 	       g.molecules[imol_new+2].initialize_coordinate_things_on_read_molecule("Auto-built C-alphas"); 
// 	       g.molecules[imol_new+2].atom_sel = bc.grown_Cas(); 
// 	       g.molecules[imol_new+2].make_ca_bonds(2.6, 4.3);  // Yikes! :-)
// 	       g.n_molecules++;

	       //
	       asc_to_graphics(bc.grown_Cas(), "Auto-built C-alphas", CA_BONDS, 2.6,4.3); 

	       bc.grown_Cas().mol->WritePDBASCII("autobuilt.pdb"); 

	    }
	 }
      }
   }
   graphics_draw();
}


// return -1 on failure
int
graphics_info_t::find_atom_index_in_moving_atoms(char *chain_id, int resno, char *atom_name) const { 

   int iret = -1;
   if (moving_atoms_asc->mol != NULL) { 
      int SelHnd = moving_atoms_asc->mol->NewSelection();
      moving_atoms_asc->mol->SelectAtoms (SelHnd, 0, chain_id,
					  resno, // starting resno, an int
					  "*", // any insertion code
					  resno, // ending resno
					  "*", // ending insertion code
					  "*", // any residue name
					  atom_name, // atom name
					  "*", // elements
					  "*"  // alt loc.
					  );

      int nSelAtoms; 
      mmdb::PPAtom local_SelAtom = NULL; 
      moving_atoms_asc->mol->GetSelIndex(SelHnd, local_SelAtom, nSelAtoms);   
      if (nSelAtoms > 0) {
	 for(int i=0; i<moving_atoms_asc->n_selected_atoms; i++) { 
	    // pointer comparison
	    if (local_SelAtom[0] == moving_atoms_asc->atom_selection[i]) {
	       iret = i;
	       break;
	    }
	 }
      }
   } else {
      std::cout << "INFO:: moving atoms mol is empty (in find_atom_index_in_moving_atoms)!\n"
		<< std::endl;
   }
   return iret;
} 





      std::vector<int> ins((max_resno - min_resno + 1), 0);
      for (int i=0; i<single_insertions.size(); i++) {
	 std::cout << " adjusting for insertion at " << single_insertions[i].resno
		   << " from " << single_insertions[i].resno+1 << " to "
		   << max_resno << std::endl;
	 for(int idx=single_insertions[i].resno+1; idx<=max_resno; idx++) {
	    ins[idx-min_resno]++;
	 }
      }

      // Let's see our insertion offsets as a function of residue
      // number:
      int tandem_insertion_start_resno;
      for(int idx=0; idx<ins.size(); idx++) {
	 std::cout << idx + min_resno + 1 << " " << ins[idx];
	 if (idx > 0)
	    if (ins[idx] != ins[idx-1]+1) {
	       std::cout << " interesting ";
	    }
	 std::cout << std::endl;
      }

      for(int idx=0; idx<ins.size(); idx++) {
	 
      }

   }




		      
		      // 		      ------------------------------------
   std::cout << "INFO:: There are " << residue_info_edits->size()
	     << " edits to apply" << std::endl;
   for (int i=0; i<residue_info_edits->size(); i++) {
      int imol = (*residue_info_edits)[i].molecule_number;
      if (imol < n_molecules) {
	 if (imol >= 0) { 
	    // a molecule-info-class-other function
	    molecules[imol].apply_atom_edit((*residue_info_edits)[i]);
	 } else {
	    std::cout << "ERROR::  Trapped error in apply_residue_info_changes"
		      << std::endl << "imol is " << imol << std::endl;
	 }
      } else {
	 std::cout << "ERROR::  Trapped error in apply_residue_info_changes"
		   << std::endl << "imol is " << imol << std::endl;
      }
   }

   if (residue_info_edits->size() > 0) {
      int imol = (*residue_info_edits)[0].molecule_number;
      molecules[imol].make_bonds_type_checked();
      gtk_widget_draw(glarea, NULL);
   }
		      

// ----
// pair.second = 0 for failure
// pair.first  = 1 for success
// 
std::pair<clipper::RTop_orth, short int>
molecule_class_info_t::get_ori_to_this_res(mmdb::Residue *res) const {

   std::pair<clipper::RTop_orth, short int>  pair;
   // clipper::RTop_orth op;

   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   res->GetAtomTable(residue_atoms, nResidueAtoms);
   if (nResidueAtoms == 0) {
      std::cout << " something broken in atom residue selection in ";
      std::cout << "write_residue_as_mol, got 0 atoms" << std::endl;
   } else {

      clipper::Coord_orth ca, c, n;
      int n_found = 0;
      
      for(int iat=0; iat<nResidueAtoms; iat++) {
	 std::string atom_name = residue_atoms[iat]->name;
	 if (atom_name == " CA ") {
	    n_found++;
	    ca = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
	 if (atom_name == " C  ") {
	    n_found++;
	    c  = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
	 if (atom_name == " N  ") {
	    n_found++;
	    n  = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
      }

      if (n_found != 3) {
	 std::cout << "DISASTER! Not all necessary atoms found in residue\n";
	 pair.second = 0; // failure
	 return pair;
      }

      // testing
//       ca = clipper::Coord_orth(0,0,0);
//       n  = clipper::Coord_orth(0.87, 0, 1.23);
//       c  = clipper::Coord_orth(0.83, 0, -1.18);
//       clipper::Coord_orth cb(-1.03, -1.11, 0);

      
      // now get the vectors of the orientation:
      //
      clipper::Coord_orth can_unit = clipper::Coord_orth((n - ca).unit());
      clipper::Coord_orth cac_unit = clipper::Coord_orth((c - ca).unit());

      clipper::Coord_orth bisector ((can_unit + cac_unit).unit());
      clipper::Coord_orth diff_unit((can_unit - cac_unit).unit());

      clipper::Coord_orth cross_prod(clipper::Coord_orth::cross(diff_unit,bisector));
      clipper::Coord_orth cpu = clipper::Coord_orth(cross_prod.unit());

//       std::cout << "bisector   " << bisector.format() << std::endl;
//       std::cout << "diff_unit  " << diff_unit.format() << std::endl;
//       std::cout << "cross prod " << cross_prod.format() << std::endl;

      // bisector   -> new x axis
      // diff_unit  -> new z axis
      // cross_prod -> new y axis

      clipper::Mat33<double> m(  bisector.x(),   bisector.y(),   bisector.z(),
			              cpu.x(),        cpu.y(),        cpu.z(),
			        diff_unit.x(),  diff_unit.y(),  diff_unit.z());

      pair.first = clipper::RTop_orth(m.transpose(), ca);
      pair.second = 1;

   }
   return pair;
} 

void do_create_mutate_sequence_gui() {

   GtkWidget *w = wrapped_create_mutate_sequence_dialog();

   GtkWidget *molecule_option_menu = lookup_widget(w, "align_and_mutate_molecule_optionmenu");
   GtkWidget *chain_option_menu = lookup_widget(w, "align_and_mutate_chain_optionmenu");
   GtkWidget *textwindow = lookup_widget(w, "align_and_mutate_sequence_text");
   GtkSignalFunc callback = GTK_SIGNAL_FUNC(align_and_mutate_molecule_menu_item_activate);
   GtkSignalFunc chain_callback = GTK_SIGNAL_FUNC(align_and_mutate_chain_option_menu_item_activate);

   graphics_info_t g;
   // get the active molecule:
   int imol = graphics_info_t::align_and_mutate_imol;
   if (imol == -1) { 
      for (int i=0; i<g.n_molecules; i++) {
	 if (g.molecules[i].has_model()) {
	    imol = i;
	    break;
	 }
      }
   }
   if (imol >= 0) { 
   
      g.fill_option_menu_with_coordinates_options(molecule_option_menu, callback, imol);
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu, imol,
								      chain_callback);
      graphics_info_t::align_and_mutate_chain_from_optionmenu = set_chain;
   }
   
   gtk_widget_show(w);
}






void add_to_history(std::vector<std::string> ls) {

   std::string command;
   if (ls.size() > 0) {
      graphics_info_t g;
      if (g.python_history) {
	 command = pythonize_command_name(ls[0]);
	 command += " (";
	 for (int i=1; i<ls.size()-1; i++) {
	    command += ls[i];
	    command += ", ";
	 }
	 command += ls.back();
	 command += ")";

	 // std::cout << "python history command: " << command << std::endl;
	 add_to_python_history(command);
      }
      if (g.guile_history) {
	 command = "(";
	 for (int i=0; i<ls.size()-1; i++) {
	    command += ls[i];
	    command += " ";
	 }
	 command += ls.back();
	 command += ")";
	 // std::cout << "guile history command: " << command << std::endl;
	 add_to_guile_history(command);
      }
   }
}


      if () {
	 // update the graph then 
	 if (vr.size() > 0) {
	    mmdb::PResidue *selresidues[vr.size()];
	    for (int i=0; i<vr.size(); i++) {
	       selresidues[i] = vr[i];
	    }
	 }
      }



   if (g.go_to_atom_molecule() >= 0) {

      std::cout << "DEBUG:: setting active position to " << g.go_to_atom_molecule()
		<< std::endl;
      gtk_menu_set_active(GTK_MENU(menu), g.go_to_atom_molecule());
      
   } else { 




// We are passed an GtkOptionMenu *option_menu
//
void fill_option_menu_with_coordinates_options(GtkWidget *option_menu,
					       GtkSignalFunc callback_func) { 

   graphics_info_t g;
   cout << "DEBUG:: fill_option_menu_with_coordinates_options called" << endl;

   // like the column labels from an mtz file, similarly fill this
   // option_menu with items that correspond to molecules that have
   // coordinates.
   //

   // Get the menu of the optionmenu (which was set in interface.c:
   // gtk_option_menu_set_menu (GTK_OPTION_MENU (go_to_atom_molecule_optionmenu),
   //                           go_to_atom_molecule_optionmenu_menu);
   //
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));

   // for the strangeness of destroying and re adding the menu to the
   // option menu, set the comments in the
   // fill_close_option_menu_with_all_molecule_options function

   gtk_widget_destroy(menu);
   menu = gtk_menu_new();

//       lookup_widget(option_menu,
// 		    "go_to_atom_molecule_optionmenu_menu");

   /* create a menu for the optionmenu button.  The various molecule
    numbers will be added to this menu as menuitems*/

   // GtkWidget *optionmenu_menu = gtk_menu_new();
   GtkWidget *menuitem;
   int menu_item_count = 0;

   for (int imol=0; imol<graphics_n_molecules(); imol++) {

//       std::cout << "in fill_option_menu_with_coordinates_options, "
// 		<< "g.molecules[" << imol << "].atom_sel.n_selected_atoms is "
// 		<< g.molecules[imol].atom_sel.n_selected_atoms << std::endl;
      
      if (is_valid_model_molecule(imol)) { 

	 std::string ss = g.molecules[imol].dotted_chopped_name();
	 menuitem = gtk_menu_item_new_with_label (ss.c_str());

	 // GtkSignalFunc callback_func = 
	 
	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     callback_func,
			     GINT_TO_POINTER(imol)); 

	 // Note that we probably don't need to do the following
	 // because we already pass a GINT_TO_POINTER(imol) in the
	 // signal connect.
	 //
	 // But on reflection.. perhaps we do because we do a
	 // menu_get_active in save_go_to_atom_mol_menu_ative_position
	 // 
	 // we set user data on the menu item, so that when this goto
	 // Atom widget is cancelled, we can whatever was the molecule
	 // number corresponding to the active position of the menu
	 //
	 // Should be freed in on_go_to_atom_cancel_button_clicked
	 // (callbacks.c)
	 // 
	 int *n = (int *) g_malloc(sizeof(int));
	 *n = imol;
	  
	 gtk_object_set_user_data(GTK_OBJECT(menuitem), (char *) n);

	 if (g.go_to_atom_molecule() == imol) {
	    gtk_menu_set_active(GTK_MENU(menu), g.go_to_atom_molecule());
	    std::cout << "DEBUG:: active menu item set to " << g.go_to_atom_molecule()
		      << " for go_to_molecule number " << g.go_to_atom_molecule()
		      << std::endl;
	 }

	 // we do need this bit of course:
	 gtk_menu_append(GTK_MENU(menu), menuitem); 
	 gtk_widget_show(menuitem);
	 menu_item_count++;
      }
   }
   
   std::cout << "DEBUG:: fill_option_menu_with_coordinates_options: "
	     << "g.go_to_atom_mol_menu_active_position "
	     << g.go_to_atom_mol_menu_active_position << std::endl;
   
   std::cout << "DEBUG:: fill_option_menu_with_coordinates_options: "
	     << "g.go_to_atom_molecule() "
	     << g.go_to_atom_molecule() << std::endl;

   // set any previously saved active position:
   // as it used to be until 20050526
//    if (g.go_to_atom_mol_menu_active_position >= 0) {
//       gtk_menu_set_active(GTK_MENU(menu),g.go_to_atom_mol_menu_active_position);

   if (g.go_to_atom_molecule() < 0) {
      // set active to the first displayed molecule:
      int ipos = 0;
      int iset = 0;
      for (int imol=0; imol<graphics_n_molecules(); imol++) {
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    if (graphics_info_t::molecules[imol].drawit) { 
 	       std::cout << "Setting active molecule position to " << ipos
 			 << std::endl;
	       gtk_menu_set_active(GTK_MENU(menu),ipos);  
	       iset = 1;
	       break;
	    }
	    ipos++;
	 }
      }
      // std::cout << "DEBUG:: iset: " << iset << std::endl;
   }

   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);
}
      



 	 // add in sigmaA calculation c.f. make_map_from_cif_generic()
	 
	 std::cout << "sigmaa and scaling..." << std::endl; 
	 
	 clipper::HKL_data<clipper::datatypes::E_sigE<float> > eo(mydata); 
	 clipper::HKL_data<clipper::datatypes::E_sigE<float> > ec(mydata);
	 // sigmaa foms (and phases) go here
	 clipper::HKL_data<clipper::datatypes::Phi_fom<float> > phifom(mydata);
	 
	 // make a list of working data for the E's
	 for ( clipper::HKL_info::HKL_reference_index ih=mydata.first();
	       !ih.last(); ih.next() ) {
	    if (!myfsigf[ih].missing()) {
	       eo[ih].E() =  myfsigf[ih].f()/sqrt(ih.hkl_class().epsilon()); 
	       ec[ih].E() = fphidata[ih].f()/sqrt(ih.hkl_class().epsilon());
	       eo[ih].sigE() = 1.0;
	       ec[ih].sigE() = 1.0;
	    }
	 }
	 // Scales for the Es.
	 int nprm = 10;
	 std::vector<clipper::ftype> params_init( nprm, 1.0 );
	 std::cout << "aquiring scales..." << std::endl; 
	 clipper::BasisFn_binner basis_fo(mydata, nprm, 2.0);
	 clipper::TargetFn_scaleEsq<clipper::datatypes::E_sigE<float> > target_fo(eo); // Fo**2
	 clipper::ResolutionFn meanfosq(mydata, basis_fo, target_fo, params_init);
	 
	 clipper::BasisFn_binner basis_fc(mydata, nprm, 2.0);
	 clipper::TargetFn_scaleEsq<clipper::datatypes::E_sigE<float> > target_fc(ec);
	 clipper::ResolutionFn meanfcsq(mydata, basis_fc, target_fc, params_init);
	 
	 //  get |Eo|, |Ec|
	 for ( clipper::HKL_info::HKL_reference_index ih=mydata.first();
	       !ih.last(); ih.next() ) {
	    if (!myfsigf[ih].missing()) {
	       eo[ih].E() *= sqrt(meanfosq.f(ih));
	       ec[ih].E() *= sqrt(meanfcsq.f(ih));
	       // 	std::cout << "meanfosq and meanfcsq: " << meanfosq.f(ih)
	       // 		  << "   " << meanfcsq.f(ih) << std::endl;
	    }
	 }
	 
	 std::cout << "refining scales..." << std::endl; 
	 // Refine those values:
	 clipper::BasisFn_binner basis_sigmaa(mydata, nprm, 2.0);
	 clipper::TargetFn_sigmaa_omegaa<clipper::datatypes::E_sigE<float> > target_sigmaa(eo,ec); // get sigmaa
	 clipper::ResolutionFn_nonlinear sigmaa(mydata,basis_sigmaa, target_sigmaa,
						params_init, 3.0);


	 std::cout << "constructing coefficients..." << std::endl; 
	 clipper::ftype x, sigma_a;
	 for ( clipper::HKL_info::HKL_reference_index ih=mydata.first();
	       !ih.last(); ih.next() ) {
	    if (!myfsigf[ih].missing()) {
	       phifom[ih].phi() = fphidata[ih].phi();
	       sigma_a = target_sigmaa.sigmaa( sigmaa.f(ih));
	       // 	std::cout <<  "sigmaa.f(ih) " <<  sigmaa.f(ih)
	       // 		  << ", Eo E: " << eo[ih].E()
	       // 		  << ", fo : " <<  myfsigf[ih].f()
	       // 		  << ", ec[ih].E(): " << ec[ih].E()
	       // 		  << ", sigma_a: " << sigma_a << std::endl; 
	       x = 2.0*eo[ih].E() * ec[ih].E()*sigma_a/(1.0-sigma_a*sigma_a);
	       if (ih.hkl_class().centric()) {
		  phifom[ih].fom() = clipper::Util::sim(x);
	       } else {
		  phifom[ih].fom() = tanh(0.5*x);
	       }
	    }
	 }
  
	 map_fphidata.update(); // do we need to do this?
	 map_fphidata.compute(myfsigf, phifom,
			      clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());


   if (widget_type == COOT_ACCEPT_REJECT_WINDOW) {
      if (graphics_info_t::accept_reject_dialog_x_position > -100) {
	 if (graphics_info_t::accept_reject_dialog_y_position > -100) {
	    gtk_widget_set_uposition(window,
				     graphics_info_t::accept_reject_dialog_x_position,
				     graphics_info_t::accept_reject_dialog_y_position);
	 }
      }
   }

   
	 



   /* We try as .phs and .cif files.

   Called after we select the mtz filename, sets up the column label
   widget and displays it. */
void
manage_column_selector(const char *filename) { 
   
   struct f_phi_columns_t a;
   int i;
   GtkWidget *column_label_window;
   GtkWidget *optionmenu1_menu, *optionmenu2_menu, *optionmenu3_menu;
   GtkWidget *menuitem;
   GtkWidget *resolution_frame;
   
   GtkWidget *optionmenu_f, *optionmenu_phi, *optionmenu_weight;
   GtkCheckButton *check_weights;
   int is_phs;
   int is_cif;
   int is_cns_data;
   int is_expert_mode_flag = 0;


   a = get_f_phi_columns(filename); 

   if (a.read_success == 0 ) { /*  not a valid mtz file */
      printf("%s is not a valid mtz file\n", filename);
      is_phs = try_read_phs_file(filename); /* Try reading the data file 
						   as an XtalView .phs file */ 
      if (! is_phs)
	 is_cif = try_read_cif_file(filename); 

      if (! is_cif) 
	 is_cns_data = try_read_cns_data_file(filename);

      return; 
   } 

   /* Else filename was OK */

   /* Save the global.  */
   save_f_phi_columns = a;

/*    printf("The F columns: \n"); */
/*    for (i=0; i< a.n_f_cols; i++) {  */
/*       printf("%s\n", a.f_cols[i].column_label); */
/*    } */

/*    printf("The P columns: \n"); */
/*    for (i=0; i< a.n_phi_cols; i++) {  */
/*       printf("%s\n", a.phi_cols[i].column_label); */
/*    } */


/* Recall that save_f_phi_columns is declared global in mtz-extras.h */

   column_label_window = create_column_label_window(); 

   /* compiler warning: we need to pass a char * not const char * for
      filename  */
   gtk_object_set_user_data(GTK_OBJECT(column_label_window), 
			    (char *) filename);

   filename = (const char *) gtk_object_get_user_data(GTK_OBJECT(column_label_window));


   /* The default column labels are at the top of the list.  The
      selcted_{f,phi}_cols get changed by the menubutton function
      (callbacks I suppose, but they are not listed there, they are in
      this file). */
   save_f_phi_columns.selected_f_col      = 0;
   save_f_phi_columns.selected_phi_col    = -1; /* unset */
   save_f_phi_columns.selected_weight_col = 0;
   save_f_phi_columns.selected_refmac_fobs_col = 0;
   save_f_phi_columns.selected_refmac_sigfobs_col = 0;
   save_f_phi_columns.selected_refmac_r_free_col = 0;


   optionmenu_f = lookup_widget(column_label_window, "optionmenu1"); 

   fill_f_optionmenu(optionmenu_f, is_expert_mode_flag); /* flag: 0 */


   /* And now do the same for the P (phase) columns */

   optionmenu2_menu = gtk_menu_new(); 
   
   /* set the default to the top for phases (if we have phase columns) */
   if (a.n_phi_cols > 0) 
     save_f_phi_columns.selected_phi_col    = 0; /* set to top */

   for (i=0; i< a.n_phi_cols; i++) { 
      menuitem = make_menu_item(a.phi_cols[i].column_label,
				GTK_SIGNAL_FUNC(phase_button_select),
				GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(optionmenu2_menu), menuitem);
      gtk_widget_show(menuitem);
      if (! (strncmp(a.phi_cols[i].column_label, "PHWT", 4))) { 
				/* was PHWT */
	gtk_menu_set_active(GTK_MENU(optionmenu2_menu), i);
	save_f_phi_columns.selected_phi_col = i;
      }
   }

   
   optionmenu_phi = lookup_widget(column_label_window, 
				    "optionmenu2"); 

   /* Link the menu  to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu_phi), optionmenu2_menu);


   /* And now do the same for the W (weight) columns */


   optionmenu3_menu = gtk_menu_new(); 
   
   for (i=0; i< a.n_weight_cols; i++) { 
      menuitem = make_menu_item(a.weight_cols[i].column_label,
				GTK_SIGNAL_FUNC(weight_button_select),
				GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(optionmenu3_menu), menuitem);
      gtk_widget_show(menuitem);
   }
   
   optionmenu_weight = lookup_widget(column_label_window, 
				    "optionmenu3"); 
 
   /* Link the menu  to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu_weight), 
			    optionmenu3_menu);


   /* By default, we want the use weights checkbutton to be off */

   check_weights = GTK_CHECK_BUTTON(lookup_widget(column_label_window, 
						  "use_weights_checkbutton"));
   
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_weights), FALSE);


   /* New addition: the refmac buttons  */
   setup_refmac_parameters(column_label_window, a);

  /* and now show the widget finally. */
   
   gtk_widget_show (column_label_window);


}


 
coot::mtz_column_types_info_t
coot::get_f_phi_columns(const std::string &filename) {

   coot::mtz_column_types_info_t a;
   a.read_success = 0;

   CMtz::MTZ *m = CMtz::MtzGet(filename.c_str(), 0);
   if (m) { 
      CMtz::MTZXTAL **xtals = CMtz::MtzXtals(m);

      int nxtals = 0;
      while (*xtals) {
	 nxtals++;
	 xtals++;
      }
      std::cout << "INFO:: There were " << nxtals << " crystals in " << filename
		<< std::endl;
   }
   return a;
}




 
	    // I am confused about the indexing of Ca1 and Ca2: Hmmm..
	    //
	    // Question:  what is the range of acceptable vaules for Ca2 indexing?
	    //            For atom_selection2 it's n_selected_atoms_2, isn't it?
	    //            But for 1kna.pdb onto 1pdq.pdb SSMAlign->Ca2[52] = -1074790400
	    //

	    // nres2 is for reference structure and nres1 is for
	    // moving structure.
	    std::cout << "DEBUG:: SSMAlign->nres1 (moving) " << SSMAlign->nres1
		      << " SSMAlign->nres2 (reference) " << SSMAlign->nres2 << std::endl;
	    for (int i=0; i<SSMAlign->nsel1 || i<SSMAlign->nsel2; i++) {
	       res_type1 = "---";
	       res_no_str_1 = "--";
	       slc_1 = "-";
	       if (i<SSMAlign->nres1)
		  if (SSMAlign->Ca1[i] != -1)
		     if (SSMAlign->Ca1[i] < n_selected_atoms_1) {
			res_type1 = atom_selection1[SSMAlign->Ca1[i]]->residue->name;
			rn1       = atom_selection1[SSMAlign->Ca1[i]]->GetSeqNum();
			res_no_str_1 = graphics_info_t::int_to_string(rn1);
			slc_1 = coot::util::three_letter_to_one_letter(res_type1);
		     }
	       reference_seq += slc_1;
	    }
	    
	    for (int i=0; i<SSMAlign->nsel1 || i<SSMAlign->nsel2; i++) {
	       res_type2 = "---";
	       res_no_str_2 = "--";
	       slc_2 = "-";
	       if (i<SSMAlign->nres2)
		  if (SSMAlign->Ca2[i] != -1)
		     if (SSMAlign->Ca2[i] < n_selected_atoms_2) {
// 			std::cout << "DEBUG:: geting atom selection atom "
// 				  << SSMAlign->Ca2[i] << " of " << n_selected_atoms_2 << std::endl;
			res_type2 = atom_selection2[SSMAlign->Ca2[i]]->residue->name;
			rn2       = atom_selection2[SSMAlign->Ca2[i]]->GetSeqNum();
			res_no_str_2 = graphics_info_t::int_to_string(rn2);
			slc_2 = coot::util::three_letter_to_one_letter(res_type2);
		     }
	       moving_seq    += slc_2;
	    }

	    
	    std::cout << "\n> " << name2 << "\n";
	    std::cout << reference_seq << "\n";
	    std::cout << "> " << name1 << "\n";
	    std::cout <<    moving_seq << "\n\n";






	       if (vs.size() != dots_type_data.size()) {
		  if (vs.size() > 1) { 
		     std::string new_object_name = vs[vs.size()-2];
		     new_object_name += " ";
		     new_object_name += vs[vs.size()-1];
		     if (generic_object_has_objects_p(obj_no)) { 
			obj_no = new_generic_object_number(new_object_name.c_str());
			set_display_generic_object(obj_no, 1);
		     }
		     dots_type_data = vs;
		  }
	       } else {
		  // are all the bits of vs and dots_type_data the same?
		  short int same=1;
		  if (n_bits > 1) { 
		     for (int ibit=n_bits-2; ibit<n_bits; ibit++) {
			if (vs[ibit] != dots_type_data[ibit]) {
			   std::string new_object_name = vs[vs.size()-2];
			   new_object_name += " ";
			   new_object_name += vs[vs.size()-1];
			   if (generic_object_has_objects_p(obj_no)) { 
			      obj_no = new_generic_object_number(new_object_name.c_str());
			      set_display_generic_object(obj_no, 1);
			   }
			   dots_type_data = vs;
			}
		     }
		  }
	       } 
	    



	       
      int trial_passes = 0;
      std::string trial_value = "A";

      while (trial_passes == 0) {

// 	 std::cout << trial_value << "  :";
// 	 for (int j=0; j<rv.size(); j++)
// 	    std::cout << rv[j] << " ";
// 	 std::cout << ": " << r << std::endl;

	 if (coot::is_member_p(this_model_chains, trial_value)) {
	    
	    unsigned int idx = r.find(trial_value);
	    if (idx != std::string::npos) {
	       if (int(idx) < int(int((r.length())-1))) { 
		  r = r.substr(0, idx) + r.substr(idx+1);
	       } else {
		  r = r.substr(0, idx);
	       }
	    } else {
	       if (r.length() > 1)
		  r = r.substr(1);
	       else
		  r = "A";
	    }
	 } else {
	    if (coot::is_member_p(rv, trial_value)) {
	       unsigned int idx = r.find(trial_value);
	       if (idx != std::string::npos) {
		  std::cout << "here: r is :" << r << ": and idx is " << idx << " on finding :"
			    << trial_value << ": and r length " << r.length() << std::endl;
		  if (int(idx) < int(int(r.length())-1)) {
		     std::cout << "path 1" << std::endl;
		     r = r.substr(0, idx) + r.substr(idx+1);
		  } else {
		     std::cout << "path 2" << std::endl;
		     r = r.substr(0, idx);
		  }
		  std::cout << "r is now :" << r << ":" << std::endl;
	       } else {
		  r = "A";
	       }
	    } else {
	       trial_passes = 1;
	    }
	 }
	 if (trial_passes == 0) 
	    trial_value = r.substr(0, 1);
      }
//        std::cout << "accepting mapped chain " << trial_value << " for chain "
//  		<< i << " was: " << adding_model_chains[i] << std::endl;
      rv.push_back(trial_value);
   }


      Bond_lines_container bonds(atom_sel);
      bonds_box.clear_up();
      bonds_box = bonds.make_graphical_bonds();
      bonds_box_type = coot::NORMAL_BONDS;
      gtk_widget_draw(graphics_info_t::glarea, NULL);



// --------------

	    residue_range_atom_index_1 = naii.atom_index;
	    residue_range_mol_no = naii.imol;
	    in_range_define_for_refine = 2;
	    if (a_is_pressed) {
	       rot_trans_rotation_origin_atom = 0; // flag for Ctrl left
						   // mouse behaviour (we
						   // don't want to rotate
						   // the atoms)
	       watch_cursor();
	       int auto_range_flag = 1;
	       refine(residue_range_mol_no,
		      auto_range_flag,
		      residue_range_atom_index_1,
		      residue_range_atom_index_1);
	       in_range_define_for_refine = 0;
	       normal_cursor();
	       pick_pending_flag = 0;
	       model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");


	       //
		  std::vector<std::vector<int> >
		  std::pair<std::string, std::string> atom_names =
		     c.atom_names_of_bond(i);
		  std::cout << "Highlight bond between :"
			    << atom_names.first << ": :"
			    << atom_names.second << ":" << std::endl;
	       

	 // I don't think we need these if the total window (top
	 // level) size is correct and the pane position is correct.
// 	 if (graphics_info_t::display_manager_maps_vbox_x_size != -1) {
// 	    GtkWidget *w_maps = lookup_widget(widget, "display_map_vbox");
// 	    gtk_widget_set_usize(w_maps,
// 				 graphics_info_t::display_manager_maps_vbox_x_size,
// 				 graphics_info_t::display_manager_maps_vbox_y_size);
// 	    GtkWidget *w_mols = lookup_widget(widget, "display_molecule_vbox");
// 	    gtk_widget_set_usize(w_mols,
// 				 graphics_info_t::display_manager_molecules_vbox_x_size,
// 				 graphics_info_t::display_manager_molecules_vbox_y_size);
// 	 }

		  

   int done_copy = 0;
   if (atom_sel.n_selected_atoms > 0) {

      if (residue_range_1 > residue_range_2) { 
	 int t = residue_range_2;
	 residue_range_2 = residue_range_1;
	 residue_range_1 = t;
      }
	 
      int imod = 1;
      make_backup();
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
      
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (to_chain == chain_p) { 
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = to_chain->GetResidue(ires);
	       int resno = residue_p->GetSeqNum();
	       if (resno <= residue_range_2 && resno >= residue_range_1) { 
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  mmdb::Atom *at;
		  for (int iat=0; iat<n_atoms; iat++) {
		     at = residue_p->GetAtom(iat);
		     clipper::Coord_orth p(at->x, at->y, at->z);
		     clipper::Coord_orth tp = p.transform(a_to_b_transform);
		     at->x = tp.x();
		     at->y = tp.y();
		     at->z = tp.z();
		  }
	       }
	    }
	 }
      }
   }
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   have_unsaved_changes_flag = 1;
   make_bonds_type_checked();
   trim_atom_label_table();
   update_symmetry();
   return done_copy;

}
		  
		  hydrogen_torsion_count = 0;
		  if (monomer_torsions.size() > 0) { 
		     for(unsigned int i=0; i<monomer_torsions.size(); i++) {

			if ((bond + hydrogen_torsion_count) == i) { 
			   
			   std::cout << "DEBUG:: testing indexes " << ibond_user << " vs "
				     << i << std::endl;
			   
			   std::string atom1 = monomer_torsions[i].atom_id_1();
			   std::string atom2 = monomer_torsions[i].atom_id_4();
			   
			   if ( (!r.second.is_hydrogen(atom1) && !r.second.is_hydrogen(atom2))
				|| find_hydrogen_torsions) {

			      std::string atom2 = monomer_torsions[i].atom_id_2();
			      std::string atom3 = monomer_torsions[i].atom_id_3();
			      atom_names = std::pair<std::string, std::string> (atom2, atom3);
			      std::cout << "atom names :" << atom_names.first << ": :"
					<< atom_names.second << ":" << std::endl;
			   } else {
			      hydrogen_torsion_count ++;
			   } 
			}
		     }
		  }


// Bad version of: replace chain-ids with seg-ids algorithm

      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
      
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);

	    // So what shall we use as the segid that gets transfered
	    // to the chain?  Let's use the segid of the first atom in
	    // the first residue of the chain (segids are a property
	    // of atoms, not chains).

	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    mmdb::Atom *at = 0;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int n_atoms = residue_p->GetNumberOfAtoms();
	 
	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  if (at)
		     break;
	       }
	       if (at)
		  break;
	    }

	    if (at) {
	       std::cout << "got starting atom " << at << std::endl;
	       std::string segid(at->segID);
	       std::string new_chain_id(chain_p->GetChainID());
	       new_chain_id += ":";
	       new_chain_id += segid;
	       if (new_chain_id.length() > 7)
		  new_chain_id = new_chain_id.substr(0,7);
	       chain_p->SetChainID((char *)new_chain_id.c_str());
	       changed_flag++;
	    }
	 }
      }
   }

      if (imol_active != -1) {
	 for (int imol=0; imol<n_molecules; imol++) {
	    if (imol==imol_active) {
	       gtk_menu_set_active(GTK_MENU(menu), go_to_atom_mol_menu_active_position); 
	    }
	 }
      } else {
	 // the old (pre-April 2007) way
	 if (go_to_atom_mol_menu_active_position >= 0) {
	    gtk_menu_set_active(GTK_MENU(menu), go_to_atom_mol_menu_active_position); 
	 }
      }


// This is the callback when the OK button of the residue info was pressed.
// 
void
graphics_info_t::apply_residue_info_changes(GtkWidget *dialog) {

   // This is where we accumulate the residue edits:
   std::vector<coot::select_atom_info> local_atom_edits;

   GtkWidget *residue_info_occ_vbox =
      lookup_widget(dialog, "residue_info_occ_vbox");
   GtkWidget *residue_info_tempfactor_vbox =
      lookup_widget(dialog, "residue_info_tempfactor_vbox");

   // ---------------------------------------------------------
   // OK, let's do it another way, not using residue_info_edits
   // ---------------------------------------------------------

//    coot::residue_spec_t *res_spec =
//       (coot::residue_spec_t *) gtk_object_get_user_data(GTK_OBJECT(dialog));

   GList *ls;

   ls = gtk_container_children(GTK_CONTAINER(residue_info_occ_vbox));

   int atom_count = 0;
   int imol = -1; 
   while (ls) {
      coot::select_atom_info *ai =
	 (coot::select_atom_info *) gtk_object_get_user_data(GTK_OBJECT(ls->data));
//       std::cout << "    " << atom_count << " " << ai->chain_id << " "
// 		<< ai->residue_number << " " << ai->atom_name << " "
// 		<< ai->altconf << std::endl;

      imol = ai->molecule_number;  // hehe
      mmdb::Atom *at = ai->get_atom(graphics_info_t::molecules[imol].atom_sel.mol);
      std::pair<short int, float>  occ_entry = graphics_info_t::float_from_entry(GTK_WIDGET(ls->data));
      if (occ_entry.first) {
	 if (at) { 
	    // std::cout << "    occ comparison " << occ_entry.second << " "
	    // << at->occupancy << std::endl;
	    if (abs(occ_entry.second - at->occupancy) > 0.009) {
	       coot::select_atom_info local_at = *ai;
	       local_at.add_occ_edit(occ_entry.second);
	       local_atom_edits.push_back(local_at);
	    }
	 }
      }
      
      // next atom
      ls = ls->next;
      atom_count++;
   }

   // do that block again for b factor:
   // 
   ls = gtk_container_children(GTK_CONTAINER(residue_info_tempfactor_vbox));

   atom_count = 0;
   while (ls) {
      coot::select_atom_info *ai =
	 (coot::select_atom_info *) gtk_object_get_user_data(GTK_OBJECT(ls->data));
//       std::cout << "    " << atom_count << " " << ai->chain_id << " "
// 		<< ai->residue_number << " " << ai->atom_name << " "
// 		<< ai->altconf << std::endl;
      imol = ai->molecule_number;  // hehe
      mmdb::Atom *at = ai->get_atom(graphics_info_t::molecules[imol].atom_sel.mol);
      std::pair<short int, float>  temp_entry = graphics_info_t::float_from_entry(GTK_WIDGET(ls->data));
      if (temp_entry.first) {
	 if (at) { 
	    // std::cout << "    temp comparison " << temp_entry.second
	    // << " " << at->tempFactor << std::endl;
	    if (abs(temp_entry.second - at->tempFactor) > 0.009) {
	       coot::select_atom_info local_at = *ai;
	       local_at.add_b_factor_edit(temp_entry.second);
	       local_atom_edits.push_back(local_at);
	    }
	 }
      }
      
      // next atom
      ls = ls->next;
      atom_count++;
   }

   if (local_atom_edits.size() >0)
      if (imol >= 0)
	 graphics_info_t::molecules[imol].apply_atom_edits(local_atom_edits);

   residue_info_edits->resize(0);
   // delete res_spec; // can't do this: the user may press the button twice
}

// ------------------------------------
   static void graphics_draw() {
     if (glarea_2) { 
       do_expose_swap_buffers_flag = 0;
       if (glarea) { 
	 gtk_widget_queue_draw(glarea);
	 if (make_movie_flag)
	   dump_a_movie_image();
       }
       gtk_widget_queue_draw(glarea_2);

       // now swap the buffers then:
       if (display_mode == coot::HARDWARE_STEREO_MODE) {
	 graphics_info_t::coot_swap_buffers(glarea, 1); // flag for hardware stereo
       } else {
	 graphics_info_t::coot_swap_buffers(glarea, 0);
       }
       if (glarea_2)
	 graphics_info_t::coot_swap_buffers(glarea_2, 0); 

       do_expose_swap_buffers_flag = 1;
     } else { 
       // there was only one gl context, so normal thing
       do_expose_swap_buffers_flag = 1;
       gtk_widget_queue_draw(glarea);
       // graphics_info_t::coot_swap_buffers(glarea, 0);
     }
   }

   static void coot_swap_buffers(GtkWidget *widget, short int in_stereo_flag) { 
     if (! in_stereo_flag) { 
       /* Swap backbuffer to front */
       GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);
       if (gdk_gl_drawable_is_double_buffered (gldrawable)) { 
	 gdk_gl_drawable_swap_buffers (gldrawable);
       } else { 
	 glFlush ();
       }
       graphics_info_t::Increment_Frames();
     }
   }

// ------------------------------------------

	    std::cout << "DEBUG::::::::::::::::::" << std::endl;
	    for (int idebug=0; idebug<bfa_chain_info.size(); idebug++)
	       std::cout << "DEBUG::::: " << idebug << " "
			 << bfa_chain_info[idebug].residue_properties.size()
			 << std::endl;
	    

      std::cout << "::::::::::::::::::::::::::::::::::::: mol (post-filtering): "
		<< mol.count_atoms() << std::endl;
      for(int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
	 for(int ires=mol[ifrag].min_res_no(); ires<=mol[ifrag].max_residue_number(); ires++) {
	    for (int iat=0; iat<mol[ifrag][ires].atoms.size(); iat++) {
	       std::cout << "::::::::::::::::::::::: mol atom: "
			 << " " << mol[ifrag].fragment_id
			 << " " << mol[ifrag][ires].seqnum
			 << " " << mol[ifrag][ires][iat].name
			 << " :" << mol[ifrag][ires][iat].altLoc
			 << ": " << mol[ifrag][ires][iat].pos.format() << std::endl;
	    }
	 }
      }

      std::cout << "::::::::::::::::::::::::::::::::::::: range_mol: "
		<< range_mol.count_atoms() << std::endl;
      for(int ifrag=0; ifrag<range_mol.fragments.size(); ifrag++) {
	 for(int ires=range_mol[ifrag].min_res_no(); ires<=range_mol[ifrag].max_residue_number(); ires++) {
	    for (int iat=0; iat<range_mol[ifrag][ires].atoms.size(); iat++) {
	       std::cout << "::::::::::::::::::::::: range mol atom: "
			 << " " << range_mol[ifrag].fragment_id
			 << " " << range_mol[ifrag][ires].seqnum
			 << " " << range_mol[ifrag][ires][iat].name
			 << " :" << range_mol[ifrag][ires][iat].altLoc
			 << ": " << range_mol[ifrag][ires][iat].pos.format() << std::endl;
	    }
	 }
      }
      
      std::cout << "::::::::::::::::::::::::::::::::::::: moved_mol: "
		<< moved_mol.count_atoms() << std::endl;
      for(int ifrag=0; ifrag<moved_mol.fragments.size(); ifrag++) {
	 for(int ires=moved_mol[ifrag].min_res_no(); ires<=moved_mol[ifrag].max_residue_number(); ires++) {
	    for (int iat=0; iat<moved_mol[ifrag][ires].atoms.size(); iat++) {
	       std::cout << "::::::::::::::::::::::: moved mol atom: "
			 << " " << moved_mol[ifrag].fragment_id
			 << " " << moved_mol[ifrag][ires].seqnum
			 << " " << moved_mol[ifrag][ires][iat].name
			 << " :" << moved_mol[ifrag][ires][iat].altLoc
			 << ": " << moved_mol[ifrag][ires][iat].pos.format() << std::endl;
	    }
	 }
      }

      std::cout << "::::::::::::::::::::::::::::::::::::: mol (pre-filtering): "
		<< mol.count_atoms() << std::endl;
      for(int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
	 for(int ires=mol[ifrag].min_res_no(); ires<=mol[ifrag].max_residue_number(); ires++) {
	    for (int iat=0; iat<mol[ifrag][ires].atoms.size(); iat++) {
	       std::cout << "::::::::::::::::::::::: mol atom: "
			 << " " << mol[ifrag].fragment_id
			 << " " << mol[ifrag][ires].seqnum
			 << " " << mol[ifrag][ires][iat].name
			 << " :" << mol[ifrag][ires][iat].altLoc
			 << ": " << mol[ifrag][ires][iat].pos.format() << std::endl;
	    }
	 }
      }
       


			      for (it=fixed_atom_specs.begin();
				   it != fixed_atom_specs.end();
				   it++) {
				 if (atom_spec == *it) {
				    fixed_atom_specs.erase(it);
				    break;
				 }
			      }

			      for (unsigned int ispec=0; ispec<fixed_atom_specs.size(); ispec++) {
				 if (atom_spec == fixed_atom_specs[ispec]) {
				    // delete this!
				    std::cout << "delete this fixed atom" << std::endl;
				    break;
				 }
			      }


#ifdef USE_PYTHON
PyObject *mark_intermediate_atom_as_fixed_py(int imol, PyObject *atom_spec, int state) {
   PyObject *retval = Py_False;
   std::pair<bool, coot::atom_spec_t> p = make_atom_spec_py(atom_spec);
   if (p.first) {
      graphics_info_t::mark_atom_as_fixed(imol, p.second, state);
   }
   return retval;
}
#endif // USE_PYTHON 





   } else {
       mmdb::mmcif::PData   mmCIF;
       mmCIF = new CMMCIFData();
       int rc = mmCIF->ReadMMCIFData ((char *)cif_dictionary_filename.c_str());
       if (rc==CIFRC_CantOpenFile)  {
 	 std::cout << "Can't read cif file " << cif_dictionary_filename << std::endl;
       } else {
 	 std::cout << "Dictionary cif file " << cif_dictionary_filename
 		   << " ReadMMCIFData OK" << std::endl;
       } 

//       mmdb::mmcif::File ciffile;
//       int ierr = ciffile.ReadMMCIFFile((char *)cif_dictionary_filename.c_str());
//       if (ierr!=mmdb::mmcif::CIFRC_Ok) {
// 	 std::cout << "dirty mmCIF file? " << cif_dictionary_filename << std::endl;
// 	 std::cout << "    Bad mmdb::mmcif::CIFRC_Ok on ReadMMCIFFile" << std::endl;
//       } else {
// 	 std::cout << "Dictionary cif file " << cif_dictionary_filename
// 		   << " ReadMMCIFFile OK" << std::endl;
//       } 

      Pmmdb::mmcif::File ciffile;
      ciffile = new mmdb::mmcif::File();
      int ierr = ciffile->ReadMMCIFFile((char *)cif_dictionary_filename.c_str());
      if (ierr!=mmdb::mmcif::CIFRC_Ok) {
	 std::cout << "dirty mmCIF file? " << cif_dictionary_filename << std::endl;
	 std::cout << "    Bad mmdb::mmcif::CIFRC_Ok on ReadMMCIFFile" << std::endl;
      } else {
	 std::cout << "Dictionary cif file " << cif_dictionary_filename
		   << " ReadMMCIFFile OK" << std::endl;
      } 
   } 

	 std::string chain_id = asi.chain_id;
	 int resno_start = asi.resno_start;
	 int resno_end   = asi.resno_end;
	 std::string altconf = asi.altconf;
	 molecules[imol].atom_sel.mol->SelectAtoms(SelHnd, 0, (char *) chain_id.c_str(),
						   resno_start, // starting resno, an int
						   "*", // any insertion code
						   resno_start, // ending resno
						   "*", // ending insertion code
						   "*", // any residue name
						   "*", // atom name
						   "*", // elements
						   altconf.c_str()  // alt loc.
						   );

int test_torsion_derivs() {

   int r = 0;
   std::string file_name = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(file_name);
   std::string chain_id = "A";
   char *chn = (char *) chain_id.c_str(); // mmdb thing.  Needs updating on new mmdb?
   int resno = 59;
   coot::protein_geometry geom;
   geom.init_standard();
   int selHnd = atom_sel.mol->NewSelection();
   int nSelResidues; 
   mmdb::PResidue *SelResidues = NULL;
   atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
			chn,
			resno-1, "",
			resno+1, "",
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"",   // altLocs
			mmdb::SKEY_NEW // selection key
			);
   atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
   
   int enable_rama_refinement = 0;
   int side_step = 0;
   bool use_flanking_residues = 1;
   // bool output_numerical_gradients = 1;
   bool output_numerical_gradients = 0;
   residue_selection_t refined_res_sel =
      testing_func_probabilities_refine_fragment(atom_sel, SelResidues, nSelResidues,
						 chain_id, resno, geom,
						 enable_rama_refinement,
						 side_step,
						 use_flanking_residues,
						 output_numerical_gradients);
   return r;
}

{ 
#if defined USE_GUILE && !defined WINDOWS_MINGW

   std::string s("(molecule-chooser-gui ");
   s += "\"Choose Molecule for Undo operations\" ";
   s += "set-undo-molecule)";
   // std::cout << s << std::endl;

   safe_scheme_command(s);

#else
#ifdef USE_PYGTK
 
   std::string s("molecule_chooser_gui(");
   s += "\"Choose Molecule for Undo operations\",";
   s += "lambda imol: set_undo_molecule(imol))";
   // std::cout << s << std::endl;

   safe_python_command(s);
 
#endif // PYGTK
#endif // USE_GUILE   

   add_to_history_simple("show-set-undo-molecule-chooser");
}



#ifdef USE_PYTHON
void show_set_undo_molecule_chooser_py() {
 
#ifdef USE_PYGTK
   std::string s("molecule_chooser_gui(");
   s += "\"Choose Molecule for Undo operations\",";
   s += "lambda imol: set_undo_molecule(imol))";
   // std::cout << s << std::endl;

   safe_python_command(s);
 
   add_to_history_simple("show-set-undo-molecule-chooser");
#else
   std::cout << "BL INFO:: wont work since there is no pygtk!" << std::endl;
#endif // PYGTK
} 
#endif // USE_PYTHON




pick_info
pick_intermediate_atom(const atom_selection_container_t &SelAtom) {
   coot::Cartesian front = unproject(0.0);
   coot::Cartesian back  = unproject(1.0);
   short int pick_mode = PICK_ATOM_ALL_ATOM;
   pick_info pi = pick_atom(SelAtom, -1, front, back, pick_mode, 0);
   if (pi.success == GL_TRUE) {
      std::cout << "flashing picked intermediate atom" << std::endl;
      mmdb::Atom *at = graphics_info_t::molecules[pi.imol].atom_sel.atom_selection[pi.atom_index];
      clipper::Coord_orth co(at->x, at->y, at->z);
      graphics_info_t::flash_position(co);
   } 
   return pi;
}





#else 

   std::cout << "============== gtk2 path =================" << std::endl;

// How do you get a menu of a menu item?  (Not like this, this gets
// the submenu, (if there is a submenu attached to this menu item).

   GtkWidget *menu = gtk_menu_item_get_submenu(GTK_MENU_ITEM(item));
   const char *t = gtk_menu_get_title(GTK_MENU(menu));
   if (t) 
     std::cout << "         got a t: " << t << std::endl;
   else 
     std::cout << "         got a NULL t: " << std::endl;

   if (t) 
      graphics_info_t::change_chain_id_from_chain = t; // a static std::string

#endif // GTK_MAJOR_VERSION




   m = GL_matrix(mat);

   glPopMatrix();
   glLoadIdentity();
   m.print_matrix();
   glMultMatrixf(m.get());



   for (unsigned int i=0; i<contacting_pairs_vec.size(); i++) {
      coot::restraints_container_t rc;
      std::pair<std::string, bool> link_info =
	 rc.find_link_type_rigourous(contacting_pairs_vec[i].first, contacting_pairs_vec[i].second, geom);
      std::cout << "DEUBG:: Covalent test found link :"  << link_info.first << ": " << link_info.second
		<< std::endl;
   }


	 
	 mmdb::Atom *at = new mmdb::Atom();
	 std::string h_name = " H";
	 if (i<100) { // protection (unlikely needed of course)
	    h_name += coot::util::int_to_string(i);
	    if (i<10)
	       h_name += " ";
	 }
	 at->SetAtomName(h_name.c_str());
	 at->SetElementName(" H");
	 at->SetCoordinates(pt.x(), pt.y(), pt.z(), 1.0, 20.0);
	 ligand_res_copy->AddAtom(at);


// ================================ EM map ======================================

      // Was this an EM map, or some other such "not crystallographic"
      // entity? (identified by P1 and angles 90, 90, 90).
      //
      // if so, fill nx_map.
      // 
      if (xmap_list[0].spacegroup().num_symops() == 1) { // P1
	 if (((xmap_list[0].cell().descr().alpha() - M_PI/2) <  0.0001) && 
	     ((xmap_list[0].cell().descr().alpha() - M_PI/2) > -0.0001) &&
	     ((xmap_list[0].cell().descr().beta()  - M_PI/2) > -0.0001) &&
	     ((xmap_list[0].cell().descr().beta()  - M_PI/2) <  0.0001) &&
	     ((xmap_list[0].cell().descr().gamma() - M_PI/2) > -0.0001) &&
	     ((xmap_list[0].cell().descr().gamma() - M_PI/2) <  0.0001)) { 

	    std::cout << "=================== EM Map ====================== " << std::endl;

	    clipper::Grid_range gr = xmap_list[0].grid_asu();  // for P1, this is right.
	    nx_map.init(xmap_list[0].cell(),
			xmap_list[0].grid_sampling(),
			gr);

	    std::cout << "INFO:: created NX Map with grid " << nx_map.grid().format() << std::endl;

	    nx_map = 0.0;

	    if (1) { 

	       // populate the map
	       clipper::Xmap_base::Map_reference_index ix(xmap_list[0]);
	       for (clipper::NXmap_base::Map_reference_index inx = nx_map.first();
		    ! inx.last();
		    inx.next() ) {
		  ix.set_coord( inx.coord() + gr.min() );
		  nx_map[inx] = xmap_list[0][ix];
	       }
	    }

	    if (0) { 
	       clipper::NXmap_base::Map_reference_index inx;
	       for (inx = nx_map.first(); !inx.last(); inx.next()) {
		  std::cout << inx.coord().format() << "  " <<  nx_map[inx] << std::endl;
	       }
	    }

	    
	 }
      }


std::vector<std::vector<std::string> >
coot::pi_stacking_container_t::get_ligand_aromatic_ring_list(const coot::dictionary_residue_restraints_t &monomer_restraints) const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::pair<std::string, std::string> > bonds;
   for (unsigned int irest=0; irest<monomer_restraints.bond_restraint.size(); irest++) {
      if (monomer_restraints.bond_restraint[irest].type() == "aromatic") {
	 std::pair<std::string, std::string> p(monomer_restraints.bond_restraint[irest].atom_id_1_4c(),
					       monomer_restraints.bond_restraint[irest].atom_id_2_4c());
	 bonds.push_back(p);
      }
   }
   
   coot::aromatic_graph_t arom(bonds);
   std::vector<std::vector<std::string> > ring_list = arom.ring_list();

   if (0) {
      std::cout << "----------- " << ring_list.size() << " rings ---------- " << std::endl;
      for (unsigned int i=0; i<ring_list.size(); i++) {
	 std::cout << "ring " << i << "\n   ";
	 for (unsigned int j=0; j<ring_list[i].size(); j++) { 
	    std::cout << ring_list[i][j] << "  ";
	 }
	 std::cout << std::endl;
      }
   }
   return ring_list;
}


	    GList *ls_2 = gtk_container_children(GTK_CONTAINER(item));

	    while (ls_2) {

	       GtkMenuItem *item_2 = GTK_MENU_ITEM(ls_2->data);
	       std::string t_2 = gtk_menu_item_get_label(GTK_MENU_ITEM(item_2));
	       std::cout << " calculate menu child  " << ls_2->data << " " << t_2 << std::endl;
	       ls_2 = ls_2->next;
	    }

--with-enterprise-tools with optional complicated/experimental ligand-related dependencies

	       // the mapping between atom name and the atom index in
	       // the residue/rdkit-molecule.
	       // 
	       // std::map<std::string, int> name_map =
	       // coot::make_flat_ligand_name_map(flat_res);

// ------------- helper function to orient_view() --------------------
// can throw an std::runtime exception;
clipper::Coord_orth
molecule_class_info_t::get_vector(const coot::residue_spec_t &central_residue_spec, // ligand typically
				  const coot::residue_spec_t &neighbour_residue_spec) const {

   clipper::Coord_orth r(0,0,0);

   mmdb::Residue *central_residue = get_residue(central_residue_spec);
   mmdb::Residue *neighbour_residue = get_residue(neighbour_residue_spec);

   if (! central_residue) {
      std::string message = "Missing residue ";
      message += central_residue_spec.chain;
      message += " "; 
      message += central_residue_spec.resno;
      throw std::runtime_error(message);
   } else { 
      if (! neighbour_residue) {
	 std::string message = "Missing residue ";
	 message += neighbour_residue_spec.chain;
	 message += " ";
	 message += neighbour_residue_spec.resno;
	 throw std::runtime_error(message);
      } else {
	 // OK! go.

	 double min_dist = 42e32;
	 clipper::Coord_orth shortest_dist(0,0,0); // "best"
	 mmdb::PPAtom c_residue_atoms = 0;
	 int c_n_residue_atoms;
	 central_residue->GetAtomTable(c_residue_atoms, c_n_residue_atoms);
	 mmdb::PPAtom n_residue_atoms = 0;
	 int n_n_residue_atoms;
	 neighbour_residue->GetAtomTable(n_residue_atoms, n_n_residue_atoms);
	 for (unsigned int iat=0; iat<c_n_residue_atoms; iat++) {
	    if (! c_residue_atoms[iat]->isTer()) { 
	       clipper::Coord_orth pt_1(c_residue_atoms[iat]->x,
					c_residue_atoms[iat]->y,
					c_residue_atoms[iat]->z);
	       for (unsigned int jat=0; jat<n_n_residue_atoms; jat++) {
		  if (! c_residue_atoms[jat]->isTer()) { 
		     clipper::Coord_orth pt_2(n_residue_atoms[jat]->x,
					      n_residue_atoms[jat]->y,
					      n_residue_atoms[jat]->z);
		     double d = (pt_1 - pt_2).lengthsq();
		     if (d < min_dist) {
			d = min_dist;
			shortest_dist = pt_2 - pt_1;
			r = shortest_dist;
		     } 
		  }
	       }
	    }
	 }

	 if (shortest_dist.lengthsq() < 0.0001) {
	    std::string message = "bad inter-residue vector: No atoms or overlapping atoms?";
	    throw std::runtime_error(message);
	 }
      }
   }

   return r;
} 

// Tinker with the atom positions of residue
// Return 1 on success.
// We need to pass the asc for the mol because we need it for seekcontacts()
// Of course the asc that is passed is the moving atoms asc.
// 
short int 
graphics_info_t::update_residue_by_chi_change_old(mmdb::Residue *residue,
						  atom_selection_container_t &asc,
						  int nth_chi, double diff) {

   // add phi/try rotamers? no.
   coot::chi_angles chi_ang(residue, 0);
   
   // Contact indices:
   //
   coot::contact_info contact = coot::getcontacts(asc);

   std::vector<std::vector<int> > contact_indices;

   // Change the coordinates in the residue to which this object contains a
   // pointer:
   // 
   // Return non-zero on failure:
   // change_by() generates the contact indices (trying regular_residue = 1)
   std::pair<short int, float> istat = chi_ang.change_by(nth_chi, diff, geom_p);

   if (istat.first) { // failure

      if (0) {
	 mmdb::PPAtom residue_atoms;
	 int n_residue_atoms;
	 
	 residue->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (unsigned int i=0; i<n_residue_atoms; i++) {
	    std::cout << "before " << i << " " << residue_atoms[i] << std::endl;
	 }
      }

      // std::cout << "DEBUG:: Simple chi change failed" << std::endl;

      // Hack, hack.  See comments on get_contact_indices_from_restraints.
      //
      bool add_reverse_contacts = 0;
      contact_indices = coot::util::get_contact_indices_from_restraints(residue, geom_p, 0,
									add_reverse_contacts);
      // try looking up the residue in

      istat = chi_ang.change_by(nth_chi, diff, contact_indices, geom_p,
				chi_angles_clicked_atom_spec,
				find_hydrogen_torsions_flag);
      if (0) {
	 mmdb::PPAtom residue_atoms;
	 int n_residue_atoms;
	 
	 residue->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (unsigned int i=0; i<n_residue_atoms; i++) {
	    std::cout << "after  " << i << " " << residue_atoms[i] << std::endl;
	 }
      }
   }

   setup_flash_bond_internal(nth_chi);

   if (istat.first == 0) { // success.
      display_density_level_this_image = 1;
      display_density_level_screen_string = "  Chi ";
      display_density_level_screen_string += int_to_string(nth_chi);
      display_density_level_screen_string += "  =  ";
      display_density_level_screen_string += float_to_string( (180.8 /M_PI) * istat.second);
      statusbar_text(display_density_level_screen_string);
   } 

   return istat.first;
} 


	 
	 std::cout << "===== set the f for " << count << " reflections" << std::endl;
	 count = 0;
	 for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) {
	    if (count > 100)
	       break;
	    clipper::data32::F_phi f = fphi[ih.hkl()];
	    std::cout << ih.hkl().format() << " " << f.f() << " " << f.phi() << std::endl;
	    count++;
	 }
	 



   std::vector<fle_ligand_bond_t> get_fle_ligand_bonds(mmdb::Residue *res_ref,
						       const std::vector<mmdb::Residue *> &residues,
						       const std::map<std::string, std::string> &name_map);
   

// This function is currently not used.
// 
// if there is a name map, use it, otherwise just bond to the found
// atom.  Using the name map allows bond to hydrogens hanging off of
// ligand atoms (e.g. the H on an O or an H on an N). The Hs are not
// in res_ref (often/typically), but are in the flat residue.
// 
std::vector<coot::fle_ligand_bond_t>
coot::get_fle_ligand_bonds(mmdb::Residue *ligand_res,
			   const std::vector<mmdb::Residue *> &residues,
			   const std::map<std::string, std::string> &name_map) { 
   std::vector<coot::fle_ligand_bond_t> v;

   // do it by hand, rather than create a molecule and use SeekContacts();
   double bond_length = 3.3;

   double bl_2 = bond_length * bond_length;
   mmdb::PPAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_res->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);
   for (unsigned int iat=0; iat<n_ligand_residue_atoms; iat++) {
      std::string ele_ligand = ligand_residue_atoms[iat]->element;
      if (ele_ligand != " C") { 
	 clipper::Coord_orth ref_pt(ligand_residue_atoms[iat]->x,
				    ligand_residue_atoms[iat]->y,
				    ligand_residue_atoms[iat]->z);
	 for (unsigned int ires=0; ires<residues.size(); ires++) {
	    mmdb::PPAtom residue_atoms = 0;
	    int n_residue_atoms;
	    residues[ires]->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (unsigned int jat=0; jat<n_residue_atoms; jat++) {
	       std::string ele = residue_atoms[jat]->element;
	       if (ele != " C") { 
		  clipper::Coord_orth pt(residue_atoms[jat]->x,
					 residue_atoms[jat]->y,
					 residue_atoms[jat]->z);
		  double diff_2 = (ref_pt - pt).lengthsq();
		  if (diff_2 < bl_2) {
		     std::string ligand_atom_name = ligand_residue_atoms[iat]->name;
		     double bl = sqrt(diff_2);
		     coot::residue_spec_t res_spec(residues[ires]);

		     // donor/acceptor from the residue to the ligand
		     // 
		     int bond_type = fle_ligand_bond_t::H_BOND_DONOR_SIDECHAIN;

		     // it's not a donor from a mainchain oxygen
		     //
		     std::string residue_atom_name(residue_atoms[jat]->name);
		     if (ele == " O")
			bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_SIDECHAIN;
		     if (ele == " N")
			bond_type = fle_ligand_bond_t::H_BOND_DONOR_SIDECHAIN;
		     if (residue_atom_name == " O  ")
			bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_MAINCHAIN;
		     if (residue_atom_name == " N  ")
			bond_type = fle_ligand_bond_t::H_BOND_DONOR_MAINCHAIN;
		     if (residue_atom_name == " H  ")
			bond_type = fle_ligand_bond_t::H_BOND_DONOR_MAINCHAIN;

		     // We don't want to map between hydroxyl oxygens -> it attached H
		     // when the residue to which it has a bond is a metal!  (and in
		     // that case, we should probably remove the Hs from the ligand
		     // anyway - but not today).
		     //
		     if (is_a_metal(residues[ires])) {
			bond_type = fle_ligand_bond_t::METAL_CONTACT_BOND;
		     } else {
			// note, can't use [] because name_map is const
			std::map<std::string, std::string>::const_iterator it =
			   name_map.find(ligand_atom_name);
			if (it != name_map.end()) {
			   // If the map happens, that's presumably because we found a H
			   // attached to an N (or an H attached to an O), either way, we
			   // are sitting now on an H.
			   ligand_atom_name = it->second;

			   // is this OK?
			   // 
			   ele = " H";
			   
			   bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_SIDECHAIN;
			   if (coot::is_main_chain_p(residue_atoms[jat]))
			      bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_MAINCHAIN;
			}
		     }

		     // that was not enough.  We still need to do some
		     // more rejections.  Here's a start:
		     if (! ((ele_ligand == " O") && (ele == " O"))) {
			
			coot::fle_ligand_bond_t bond(ligand_atom_name, bond_type, bl, res_spec);
			v.push_back(bond);
			break;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return v;
}
