

// This seems to be in the wrong place.  lig_build functions should
// not depend on lbg_info_t.  Bleugh.  Fix at some stage.
// 
void
widgeted_molecule_t::highlight_maybe(int x, int y, bool delete_mode) const {

   bool highlighted_bond = 0;

   lbg_info_t l;
   l.remove_bond_and_atom_highlighting();

   double dx(x);
   double dy(y);
   for (unsigned int i=0; i<bonds.size(); i++) {
      
      if (bonds[i].over_bond(dx, dy,
			     atoms[bonds[i].get_atom_1_index()],
			     atoms[bonds[i].get_atom_2_index()])) {
	 l.highlight_bond(bonds[i], delete_mode);
	 highlighted_bond = 1;
	 break; // don't hightlight 2 bonds
      }
   }

   // check for highlighting atoms
   if (! highlighted_bond) {
      for (unsigned int i=0; i<atoms.size(); i++) {
	 if (atoms[i].over_atom(dx, dy)) {
	    l.highlight_atom(atoms[i], i, delete_mode);
	    break; // only highlight one atom
	 }
      }
   }
}


template class std::vector<lig_build::bond_t> lig_build::molecule_t<lig_build::atom_t, lig_build::bond_t>::bonds_with_vertex(const atom_position_t &pos);

		  std::cout << "bond number " << ib << "   " 
			    << "atom_index: " << atom_index << "   index_1: " << index_1
			    << "   index_2: " << index_2 << std::endl;

   if (1) {
      GooCanvasItem *rect_item_1  =
	 goo_canvas_rect_new (root,
			      candidate_pos_1.x - 1.0, 
			      candidate_pos_1.y - 1.0,
			      2.0, 2.0,
			      "line-width", 2.0,
			      "stroke-color", "blue",
			      NULL);

      GooCanvasItem *rect_item_2  =
	 goo_canvas_rect_new (root,
			      candidate_pos_2.x - 1.0, 
			      candidate_pos_2.y - 1.0,
			      2.0, 2.0,
			      "line-width", 2.0,
			      "stroke-color", "red",
			      NULL);


   }


      if (1) {
	 std::cout << "   central_atom_pos: " << central_atom_pos << std::endl;
	 std::cout << "   p1 :              " << p1 << std::endl;
	 GooCanvasItem *rect_item_1  =
	    goo_canvas_rect_new (root,
				 central_atom_pos.x - 2.0, 
				 central_atom_pos.y - 2.0,
				 4.0, 4.0,
				 "line-width", 2.0,
				 "stroke-color", "blue",
				 NULL);
	 GooCanvasItem *rect_item_2  =
	    goo_canvas_rect_new (root,
				 p1.x - 2.0, 
				 p1.y - 2.0,
				 4.0, 4.0,
				 "line-width", 2.0,
				 "stroke-color", "red",
				 NULL);
      }


class molecule_widget_item {
   GooCanvasItem *ci;
protected:
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      }
      ci = NULL;
   }
};



	 // draw that bond
	 GooCanvasItem *bci =  goo_canvas_polyline_new_line(root,
							    atom_pos_index[i].first.x,
							    atom_pos_index[i].first.y,
							    atom_pos_index[j].first.x,
							    atom_pos_index[j].first.y,
							    "stroke-color", stroke_colour.c_str(),
							    NULL);


		     // and add a line to the canvas:
		     double x1 = highlight_data.get_atom_1_pos().x;
		     double y1 = highlight_data.get_atom_1_pos().y;
		     double x2 = mol.atoms[new_atoms[0]].atom_position.x;
		     double y2 = mol.atoms[new_atoms[0]].atom_position.y;

		     GooCanvasItem *bci = goo_canvas_polyline_new_line(root, x1, y1, x2, y2,
								       "stroke-color", dark,
								       NULL);
		     widgeted_bond_t bond(highlighted_atom_index, new_atoms[0],
					  lig_build::bond_t::SINGLE_BOND, bci);


	    std::vector<int> bonds = mol.bonds_having_atom_with_atom_index(ind_1);
	    std::string ele = mol.atoms[ind_1].element;
	    std::string atom_name = mol.make_atom_name_by_using_bonds(ele, bonds);
	    mol.atoms[ind_1].update_name_maybe(atom_name, root);

	    bonds = mol.bonds_having_atom_with_atom_index(ind_2);
	    ele = mol.atoms[ind_2].element;
	    atom_name = mol.make_atom_name_by_using_bonds(ele, bonds);
	    mol.atoms[ind_2].update_name_maybe(atom_name, root);




std::vector<std::string>
lbg_info_t::compare_vs_sbase(CGraph *graph1) const {

   std::vector<std::string> v;
   if (! SBase) {
      std::cout << "SBase not initialized" << std::endl;
   } else {
      int nStructures = SBase->GetNofStructures();
      std::cout << "there are " << nStructures << " in the sbase " << std::endl;
      int H_flag = 0;
      int minMatch = 4;
      
      for (int is=0; is<nStructures; is++) {
	 CGraphMatch match;
	 CGraph *G = new CGraph;
	 int status = SBase->GetGraph(is, G, H_flag);
	 if (status !=  SBASE_Ok) {
	    std::cout << "bad status on get graph " << is << std::endl;
	 } else {
	    // good
	    match.MatchGraphs(graph1, G, minMatch, 1);
	    int n_match = match.GetNofMatches();
	    std::cout << "INFO:: match NumberofMatches (potentially similar graphs) "
		      << n_match << " matching vs structure " << is << std::endl;
	    if (n_match > 0) {
	       match.PrintMatches();
	    }
	 }
	 delete G;
      }
   }
   return v;
}



      int nStructures = SBase->GetNofStructures();
      std::cout << "there are " << nStructures << " structures in the sbase " << std::endl;
      int H_flag = 0;
      int minMatch = 4;
      
      for (int is=0; is<nStructures; is++) {
	 std::cout << " testing structure " << is << std::endl;
	 CGraph *G = new CGraph;
	 int status = SBase->GetGraph(is, G, H_flag);
	 if (status != SBASE_Ok) {
	    std::cout << "bad status on get graph " << is << std::endl;
	 } else {
	    // good
	    CGraphMatch match;
	    match.MatchGraphs(graph1, G, minMatch, 1);
	    int n_match = match.GetNofMatches();
	    std::cout << "INFO:: match NumberofMatches (potentially similar graphs) "
		      << n_match << std::endl;
	    if (n_match > 0) 
	       match.PrintMatches();
	 }
	 delete G;
      }


lbg_info_t::match_results_t
lbg_info_t::residue_from_best_match(CGraph &graph1, CGraph &graph2,
				    CGraphMatch &match, int n_match,
				    CSBStructure *SBS) const {

   lbg_info_t::match_results_t r("", "", NULL);
   
   int best_match = UNASSIGNED_INDEX;
   bool apply_rtop_flag = 1;

   bool success = 0;
   clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
   clipper::Coord_orth pt_dum(0,0,0);
   clipper::RTop_orth rtop(m_dum, pt_dum);
   clipper::RTop_orth best_rtop(m_dum, pt_dum);
   double best_match_sum = 1e20;
   int best_n_match = -99;
   std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > best_matching_atoms;

   for (int imatch=0; imatch<n_match; imatch++) {
      std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > matching_atoms; 
      int n;
      realtype p1, p2;
      ivector FV1, FV2;
      match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
      // 	    For understanding only.  
      // 	    std::cout << "Match number: " << imatch << "  " << p1*100 << "% "
      // 		      << p2*100 << "% "<< std::endl;
      std::vector<clipper::Coord_orth> coords_1_local;
      std::vector<clipper::Coord_orth> coords_2_local;
      for (int ipair=1; ipair<=n; ipair++) {
	 PCVertex V1 = graph1.GetVertex ( FV1[ipair] );
	 PCVertex V2 = graph2.GetVertex ( FV2[ipair] );
	 if ((!V1) || (!V2))  {
	    std::cout << "Can't get vertices for match "
		      << ipair << std::endl;
	 } else  {

	    CAtom *at1 = NULL; // = cleaned_res_moving->atom[V1->GetUserID()];
	    CAtom *at2 = NULL; // = cleaned_res_reference->atom[V2->GetUserID()];
	    coords_1_local.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
	    coords_2_local.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
	    std::pair<std::string, std::string> atom_info_1(at1->name, at1->altLoc);
	    std::pair<std::string, std::string> atom_info_2(at2->name, at2->altLoc);
	    std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > atom_pair(atom_info_1, atom_info_2);
	    matching_atoms.push_back(atom_pair);
	 }
      }
	    
      double dist_sum = 0.0;
      clipper::RTop_orth rtop_local(clipper::Mat33<double>(0,0,0,0,0,0,0,0,0),
				    clipper::Coord_orth(0,0,0)); // unset
      if (apply_rtop_flag) {

	 // 	       for (unsigned int iat=0; iat<4; iat++) {
	 // 		  std::cout << "debug:: getting rtop " << coords_1_local[iat].format() << " vs "
	 // 			    << coords_2_local[iat].format() << std::endl;
	 // 	       }
	       
	 rtop_local = clipper::RTop_orth(coords_1_local, coords_2_local);
	 for (unsigned int i=0; i<coords_1_local.size(); i++) {
	    dist_sum += clipper::Coord_orth::length(coords_2_local[i],
						    coords_1_local[i].transform(rtop_local));
	 }
      } else {
	 for (unsigned int i=0; i<coords_1_local.size(); i++) {
	    dist_sum += clipper::Coord_orth::length(coords_2_local[i], coords_1_local[i]);
	 }
      } 
      if (dist_sum < best_match_sum) {
	 // Debugging
	 std::cout << "DEBUG:: better dist_sum: " << dist_sum << std::endl;
	 best_rtop = rtop_local;
	 best_match_sum = dist_sum;
	 best_match = imatch;
	 best_matching_atoms = matching_atoms;
	 best_n_match = coords_1_local.size();
      }
   } // imatch loop

   if (best_match != -1) {
      rtop = best_rtop;
      success = 1;
   }
   return r;
}




   std::vector<int> bonds_with_atom_index_start;
   for (unsigned int i=0; i<bonds.size(); i++) {
      if (bonds[i].get_atom_1_index() == atom_index_start)
	 bonds_with_atom_index_start.push_back(i);
      if (bonds[i].get_atom_2_index() == atom_index_start)
	 bonds_with_atom_index_start.push_back(i);
   }


bool
widgeted_molecule_t::can_find_self_p(int atom_index_start, int depth) const {

   if (depth == 0) {
      return 0;
   } else {
      
   }
}


bool
widgeted_molecule_t::find_bonded_atoms_with_no_pass(int start_atom_index,
						    int this_atom_index,
						    const std::vector<int> &no_pass_atoms,
						    int depth) const {

   std::vector<int> atoms_bonded_to_atom_index_start;
   std::vector<int> local_no_pass_atoms = no_pass_atoms;
   if (depth == 0) {
      return 0;
   } else {

      // get a list of all the atoms that are bonded to this atom
      
      for (unsigned int i=0; i<bonds.size(); i++) { 
	 if (bonds[i].get_atom_1_index() == this_atom_index) {
	    if (! member(bonds[i].get_atom_2_index(), no_pass_atoms)) {
	       int idx = bonds[i].get_atom_2_index();
	       if (idx == start_atom_index) { 
		  return local_no_pass_atoms;
	       } else { 
		  atoms_bonded_to_atom_index_start.push_back(idx);
		  local_no_pass_atoms.push_back(this_atom_index);
	       }
	    } 
	 }
	 if (bonds[i].get_atom_2_index() == this_atom_index) {
	    if (! member(bonds[i].get_atom_1_index(), no_pass_atoms)) {
	       int idx = bonds[i].get_atom_1_index();
	       if (idx == start_atom_index) {
		  return local_no_pass_atoms;
	       } else { 
		  atoms_bonded_to_atom_index_start.push_back(idx);
		  local_no_pass_atoms.push_back(this_atom_index);
	       }
	    }
	 }
      }

      for (unsigned int iat=0; iat<atoms_bonded_to_atom_index_start.size(); iat++) { 
	 std::vector<int> v =
	    find_bonded_atoms_with_no_pass(start_atom_index,
					   atoms_bonded_to_atom_index_start[iat],
					   local_no_pass_atoms, depth-1);
      }
   }
   
   return atoms_bonded_to_atom_index_start; // compile, not correct value
}


      gtk_tree_store_append(tree_store_similarities, &toplevel, NULL);
      gtk_tree_store_set(tree_store_similarities, &toplevel, 0, f, -1);

   GtkCellRenderer *cell_renderer = gtk_cell_renderer_text_new();
   g_object_set_data (G_OBJECT (cell_renderer), "column", GINT_TO_POINTER (0));


	    // debug
	    std::cout << "words: ";
	    for (unsigned int iw=0; iw<words.size(); iw++) { 
	       std::cout << "  :" << words[iw] << ":";
	    }
	    std::cout << std::endl;

                  std::cout << "... parsing for x :" << words[2] << ":" << std::endl;
		  double pos_x = lig_build::string_to_float(words[2]);
		  std::cout << "got x " << pos_x << std::endl;
		  std::cout << "... parsing for y :" << words[3] << ":" << std::endl;
		  double pos_y = lig_build::string_to_float(words[3]);
		  std::cout << "got y " << pos_y << std::endl;
		  std::cout << "... parsing for z :" << words[4] << ":" << std::endl;
		  double pos_z = lig_build::string_to_float(words[4]);
		  std::cout << "got z " << pos_z << std::endl;
		  std::cout << "... parsing for sa :" << words[5] << ":" << std::endl;
		  double sa    = lig_build::string_to_float(words[5]);
		  std::cout << "got sa " << sa << std::endl;
		  clipper::Coord_orth pt(pos_x, pos_y, pos_z);
		  std::cout << "got atom name :" << atom_name << ": and pos "
			    << pt.format() << " and accessibility: " << sa
			    << std::endl;


   GooCanvasItem *group = goo_canvas_group_new (root,
						"line_width", 1.0,
						// "fill-color-rgba", 0x5555bb20,
						"fill-color", "blue",
						NULL);

   int n_circles = int(sa*20) + 1;
   if (n_circles> 5 ) n_circles = 5;

   for (unsigned int i=0; i<n_circles; i++) { 
      double rad = 3.0 * double(i);
       GooCanvasItem *cirle = goo_canvas_ellipse_new(group,
						     pos.x, pos.y,
						     rad, rad,
						     "fill-color", "blue",
						     NULL);
   }

      clipper::Coord_orth cp(residue_circles[i].pos_x,
			     residue_circles[i].pos_y,
			     residue_circles[i].pos_z);
      lig_build::pos_t pos = mol.input_coords_to_canvas_coords(cp);



	    // fill solvent_accessible_atoms vector before reading the
	    // molecule.  When the molecule is read, the solvent
	    // accessibilites are grafted onto the various different
	    // molecules.
	    //
	    l->read_solvent_accessibilities(sa_file);
	    
   std::string res_info_file = coot_dir + "coot-tmp-fle-view-residue-info.txt";
