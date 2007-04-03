	 if (mutate_auto_fit_do_post_refine_flag) {
	    // Run refine zone with autoaccept, autorange on
	    // the "clicked" atom:
	    CAtom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	    std::string chain_id = at->GetChainID();
	    refine_auto_range(naii.imol, chain_id.c_str(), at->GetSeqNum(),
			      at->altLoc);
	 }

realtype u, v, w;
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
	 CAtom *rot_centre = rot_trans_rotation_origin_atom;
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
      PPCAtom local_SelAtom = NULL; 
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
molecule_class_info_t::get_ori_to_this_res(CResidue *res) const {

   std::pair<clipper::RTop_orth, short int>  pair;
   // clipper::RTop_orth op;

   PPCAtom residue_atoms;
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
	    PCResidue *selresidues[vr.size()];
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
      
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (to_chain == chain_p) { 
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = to_chain->GetResidue(ires);
	       int resno = residue_p->GetSeqNum();
	       if (resno <= residue_range_2 && resno >= residue_range_1) { 
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  CAtom *at;
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
      
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 CChain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);

	    // So what shall we use as the segid that gets transfered
	    // to the chain?  Let's use the segid of the first atom in
	    // the first residue of the chain (segids are a property
	    // of atoms, not chains).

	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    CAtom *at = 0;
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
