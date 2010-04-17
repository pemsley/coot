// Unfortunately, the bonds and atoms must also have canvas items as
// part of their make-up.  Why?
//
// Consider changing a carbon to a nitrogen (or vica vera).  We want
// the molecular description of the *bond* to stay the same, we want
// the atom to change its element (and the widget representing it to
// be changed or generated or deleted (when going back to Carbon)) and
// the representation of the bond should be changed from a line to the
// atom point to a line that approaches (but does not touch) the atom
// point.  Which means that the bond and its representation are at
// different places.
//
// Now try to delete that bond. Which is the widget that needs to be
// removed from the canvas?
//
// We can only know that if the canvas item is part of the bond
// description.
//


class solvent_accessible_atom_t {
public:
   std::string atom_name;
   clipper::Coord_orth pt;
   double solvent_accessibility;
   solvent_accessible_atom_t(const std::string &at,
			     const clipper::Coord_orth &pt_in,
			     double sa) {
      atom_name = at;
      pt = pt_in;
      solvent_accessibility = sa;
   }
};



// ====================================================================
//                     widgeted_atom_t
// ====================================================================

class widgeted_atom_t : public lig_build::atom_t {
   std::string font_colour;
   double solvent_accessibility;
   GooCanvasItem *ci;
   std::string atom_name; // typically names from a PDB file.
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      }
      ci = NULL;
   }

   // then general form, use this not the below 2 (which are used by
   // this function)
   // 
   GooCanvasItem *make_canvas_text_item(const std::string atom_id_in,
					const std::string &fc,
					GooCanvasItem *root) {
      GooCanvasItem *text_item = NULL;
      if (atom_id_in != "C") {
	 if (atom_id_in.length() > 2) {
	    text_item = make_subscripted_canvas_item(atom_id_in, fc, root);
	 } else { 
	    text_item = make_convention_canvas_item(atom_id_in, fc, root);
	 }
      }
      return  text_item;
   }
   
   GooCanvasItem *make_convention_canvas_item(const std::string atom_id_in,
					       const std::string &fc,
					       GooCanvasItem *root) const {

      GooCanvasItem *text_item = goo_canvas_text_new(root, atom_id_in.c_str(),
						     atom_position.x, atom_position.y, -1,
						     GTK_ANCHOR_CENTER,
						     "font", "Sans 10",
						     "fill_color", fc.c_str(),
						     NULL);
      return text_item;
   }
   
   GooCanvasItem *make_subscripted_canvas_item(const std::string atom_id_in,
					       const std::string &fc,
					       GooCanvasItem *root) const {

      GooCanvasItem *group = goo_canvas_group_new (root,
						   "fill_color", fc.c_str(),
						   NULL);
      std::string p1 = atom_id_in.substr(0,2);
      std::string p2 = atom_id_in.substr(2);
      
      GooCanvasItem *item_1 = goo_canvas_text_new(group, p1.c_str(),
						  atom_position.x, atom_position.y, -1,
						  GTK_ANCHOR_CENTER,
						  "font", "Sans 10",
						  "fill_color", fc.c_str(),
						  NULL);

      double pos_p2_x = atom_position.x + 13;
      double pos_p2_y = atom_position.y + 2;
      
      GooCanvasItem *item_2 = goo_canvas_text_new(group, p2.c_str(),
						  pos_p2_x, pos_p2_y, -1,
						  GTK_ANCHOR_CENTER,
						  "font", "Sans 8",
						  "fill_color", fc.c_str(),
						  NULL);
      return group;
   }

public:
   widgeted_atom_t(lig_build::atom_t &at_in, GooCanvasItem *ci_in) : lig_build::atom_t(at_in) {
      ci = ci_in;
      font_colour = "yellow";
      solvent_accessibility = -1;
   }
   widgeted_atom_t(lig_build::pos_t pos_in,
		   std::string ele_in,
		   int charge_in,
		   GooCanvasItem *ci_in) :
      lig_build::atom_t(pos_in, ele_in, charge_in) {
      ci = ci_in;
      font_colour = "hotpink";
      solvent_accessibility = -1;
   }
   GooCanvasItem *get_canvas_item() const { return ci; }
   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      clear(root);
      ci = new_item;
   }
   bool update_atom_id_maybe(const std::string &atom_id_in,
			     GooCanvasItem *root) {
      return update_atom_id_maybe(atom_id_in, font_colour, root);
   }
   bool update_atom_id_maybe(const std::string &atom_id_in,
			     const std::string &fc,
			     GooCanvasItem *root) {
      bool changed_status = 0;
      font_colour = fc;
      GooCanvasItem *text_item = NULL;
      if (! is_closed()) { 
	 if (atom_id_in != get_atom_id()) {
	    changed_status = set_atom_id(atom_id_in);
	    if (changed_status) {
	       text_item = make_canvas_text_item(atom_id_in, fc, root);
	       update_canvas_item(text_item, root);
	    }
	 }
      } else {
	 update_canvas_item(text_item, root); // close atom, replace with null.
      } 
      return changed_status;
   }
   void update_atom_id_forced(const std::string &atom_id_in,
			      const std::string &fc, 
			      GooCanvasItem *root) {
      set_atom_id(atom_id_in);
      GooCanvasItem *text_item = make_canvas_text_item(atom_id_in, fc, root);
      update_canvas_item(text_item, root);
   }
   
   void update_atom_id_forced(const std::string &atom_id_in,
			      GooCanvasItem *root) {
      update_atom_id_forced(atom_id_in, font_colour, root);
   }
   void add_solvent_accessibility(double sa) {
      solvent_accessibility = sa;
   }
   double get_solvent_accessibility() const { return solvent_accessibility; } // negative for none.
   void close(GooCanvasItem *root) {
      // std::cout << " closing subclass atom" << std::endl;
      lig_build::atom_t::close();
      update_canvas_item(NULL, root);
   }
   void set_atom_name(const std::string atom_name_in) {
      atom_name = atom_name_in;
   }
   std::string get_atom_name() const {
      return atom_name;
   } 
};

// ====================================================================
//                     widgeted_bond_t
// ====================================================================

class widgeted_bond_t : public lig_build::bond_t {
   GooCanvasItem *ci;
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      }
      ci = NULL;
   }

   void construct_internal(const lig_build::atom_t &atom_first,
			   const lig_build::atom_t &atom_second,
			   bond_type_t bt, GooCanvasItem *root) {
      bool shorten_first = 0;
      bool shorten_second = 0;
      if (atom_first.atom_id != "C") { 
	 shorten_first = 1;
      } 
      if (atom_second.atom_id != "C") { 
	 shorten_second = 1;
      }
      lig_build::pos_t pos_1 =  atom_first.atom_position;
      lig_build::pos_t pos_2 = atom_second.atom_position;
      ci = canvas_item_for_bond(pos_1, pos_2, shorten_first, shorten_second, bt, root);
   }

   // all bonds are made this way...
   // 
   GooCanvasItem *canvas_item_for_bond(const lig_build::pos_t &pos_1,
				       const lig_build::pos_t &pos_2,
				       bool shorten_first,
				       bool shorten_second,
				       bond_type_t bt,
				       GooCanvasItem *root) const;

   // We need to make a shorter bond canvas line because we have (say)
   // changed a carbon to a N (bond canvas line now does not
   // completely extend to the atom position).
   // 
   void make_new_canvas_item_given_type(const lig_build::atom_t &atom_changed,
					const lig_build::atom_t &atom_other,
					lig_build::bond_t::bond_type_t bt,
					GooCanvasItem *root) {

      lig_build::pos_t A = atom_changed.atom_position;
      lig_build::pos_t B =   atom_other.atom_position;

      bool shorten_first = 0;
      bool shorten_second = 0;
      if (atom_changed.atom_id != "C")
	 shorten_first = 1;
      if (atom_other.atom_id != "C")
	 shorten_second = 1;
      GooCanvasItem *new_line = canvas_item_for_bond(A, B, shorten_first, shorten_second,
						     bt, root);
      update_canvas_item(new_line, root);
   }

   GooCanvasItem * canvas_item_double_bond(const lig_build::pos_t &pos_1,
					   const lig_build::pos_t &pos_2,
					   GooCanvasItem *root) const;
					
   GooCanvasItem * canvas_item_double_aromatic_bond(const lig_build::pos_t &pos_1,
						    const lig_build::pos_t &pos_2,
						    GooCanvasItem *root) const;
					
   GooCanvasItem * make_wedge_bond_item(const lig_build::pos_t &pos_1,
					const lig_build::pos_t &pos_2,
					const lig_build::bond_t::bond_type_t &bt,
					GooCanvasItem *root) const;
   GooCanvasItem * make_wedge_out_bond_item(const lig_build::pos_t &pos_1,
					    const lig_build::pos_t &pos_2,
					    GooCanvasItem *root) const;
   GooCanvasItem * make_wedge_in_bond_item(const lig_build::pos_t &pos_1,
					   const lig_build::pos_t &pos_2,
					   GooCanvasItem *root) const;


   // -------------------- widgeted_bond_t public --------------------------

public:

   // this is for widgeted_bond_t that are invalid (to be assigned later).
   widgeted_bond_t() : lig_build::bond_t() {
      ci = NULL;
   } 
   
   // Now we use a constructor that does the creation of the canvas item too
   //
   widgeted_bond_t(int first, int second, 
		   const lig_build::atom_t &atom_first, const lig_build::atom_t &atom_second,
		   bond_type_t bt, GooCanvasItem *root) :
      lig_build::bond_t(first, second, bt) {
      construct_internal(atom_first, atom_second, bt, root);
   }
   // as above, but we give the centre of the ring too.
   widgeted_bond_t(int first, int second, 
		   const lig_build::atom_t &atom_first, const lig_build::atom_t &atom_second,
		   lig_build::pos_t centre_pos_in,
		   bond_type_t bt, GooCanvasItem *root) :
      bond_t(first, second, centre_pos_in, bt) {
      construct_internal(atom_first, atom_second, bt, root);
   }
   
   // old constructor.  Consider deleting.
   // 
   // widgeted_bond_t(int first, int second, bond_type_t bt, GooCanvasItem *ci_in) :
   // lig_build::bond_t(first, second, bt) {
   // ci = ci_in;
   // is_closed_ = 0;
   // }

   // old constructor.  Consider deleting.
   // 
   // widgeted_bond_t(int first, int second,
   // lig_build::pos_t centre_pos_in,
   // lig_build::bond_t::bond_type_t bt,
   // GooCanvasItem *ci_in) :
   // lig_build::bond_t(first, second, centre_pos_in, bt) {
   // ci = ci_in;
   // is_closed_ = 0;
   // }


   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      clear(root);
      ci = new_item;
   }

   // We need to make a shorter bond canvas line because we have (say)
   // changed a carbon to a N (bond canvas line now does not
   // completely extend to the atom position).
   void make_new_canvas_item(const lig_build::atom_t &atom_changed,
			     const lig_build::atom_t &atom_other,
			     GooCanvasItem *root) {

      lig_build::bond_t::bond_type_t bt = get_bond_type();
      make_new_canvas_item_given_type(atom_changed, atom_other, bt, root);
   }
   void change_bond_order(const lig_build::atom_t &atom_changed,
			  const lig_build::atom_t &atom_other,
			  GooCanvasItem *root) {
      change_bond_order(atom_changed, atom_other, 0, root);
   } 
   void change_bond_order(const lig_build::atom_t &atom_changed,
			  const lig_build::atom_t &atom_other,
			  bool allow_triple_toggle,
			  GooCanvasItem *root) {
      lig_build:: atom_t at_1 = atom_changed;
      lig_build:: atom_t at_2 = atom_other;
      lig_build::bond_t::bond_type_t bt = get_bond_type();
      if (bt == lig_build::bond_t::SINGLE_BOND) { 
	 if (allow_triple_toggle)
	    bt = lig_build::bond_t::TRIPLE_BOND;
	 else 
	    bt = lig_build::bond_t::DOUBLE_BOND;
      } else { 
	 if (bt == lig_build::bond_t::DOUBLE_BOND)
	    if (allow_triple_toggle)
	       bt = lig_build::bond_t::TRIPLE_BOND;
	    else
	       bt = lig_build::bond_t::SINGLE_BOND;
	 else 
	    if (bt == lig_build::bond_t::TRIPLE_BOND)
	       bt = lig_build::bond_t::DOUBLE_BOND;
	    else
	       if (bt == lig_build::bond_t::IN_BOND) {
		  // std::cout << " add a reverse direction here " << std::endl;
		  std::swap(atom_1, atom_2);
		  std::swap(at_1, at_2);
		  bt = lig_build::bond_t::OUT_BOND;
	       } else {
		  if (bt == lig_build::bond_t::OUT_BOND) { 
		     bt = lig_build::bond_t::IN_BOND;
		  }
	       }
      }
      
//       std::cout << "changing bond type from " << get_bond_type() << " to "
// 		<< bt << std::endl;
      set_bond_type(bt);
      make_new_canvas_item_given_type(at_1, at_2, bt, root);
   }
   void close(GooCanvasItem *root) {
      // std::cout << " closing sub-class bond" << std::endl;
      lig_build::bond_t::close();
      update_canvas_item(NULL, root);
   }
   int mmdb_bond_type() const {
      int mmdb_bt = 1;
      switch (get_bond_type()) { 
      case SINGLE_BOND:
	 mmdb_bt = 1;
	 break;
      case IN_BOND:
	 mmdb_bt = 1;
	 break;
      case OUT_BOND:
	 mmdb_bt = 1;
	 break;
      case DOUBLE_BOND:
	 mmdb_bt = 2;
	 break;
      case TRIPLE_BOND:
	 mmdb_bt = 3;
	 break;
      case BOND_UNDEFINED:
	 mmdb_bt = UNASSIGNED_INDEX;
	 break;
      }
      return mmdb_bt;
   }
   void add_centre(const lig_build::pos_t &centre_in) {
      set_centre_pos(centre_in);
   }

   int get_other_index(int atom_index) const {
      int idx = get_atom_1_index();
      if (idx == atom_index)
	 idx = get_atom_2_index();
      return idx;
   }

};


// ====================================================================
//                     widgeted_molecule_t
// ====================================================================

#define MAX_SEARCH_DEPTH 9

class widgeted_molecule_t : public lig_build::molecule_t<widgeted_atom_t, widgeted_bond_t> {
   
private:
   std::string group;
   void init() {
      mol_in_max_y = 0;
      mol_in_min_y = 0;
      scale_correction.first = 0;
      scale_correction.second = 1;
   }
   bool member(const int &ind, const std::vector<int> &no_pass_atoms) const {
      bool found = 0;
      for (unsigned int i=0; i<no_pass_atoms.size(); i++) { 
	 if (no_pass_atoms[i] == ind) {
	    found = 1;
	    break;
	 }
      }
      return found;
   } 
   // Return a vector of bonds.  If empty, then it didn't find self.
   // 
   std::pair<bool, std::vector<int> >
   found_self_through_bonds(int atom_index_start, int atom_index_other) const;
   std::pair<bool, std::vector<int> >
   find_bonded_atoms_with_no_pass(int atom_index_start,
				  int atom_index_other, // must pass through this
				  int this_atom_index,
				  const std::vector<int> &no_pass_atoms,
				  int depth) const;
   void debug_pass_atoms(int atom_index, int this_atom_index, 
			 int depth,  const std::vector<int> &local_no_pass_atoms) const;
   std::pair<bool, double>
   get_scale_correction(const lig_build::molfile_molecule_t &mol_in) const;
   int get_number_of_atom_including_hydrogens() const;
   // return negative if not solvent accessibility available.
   double get_solvent_accessibility(const clipper::Coord_orth &pt,
				    const std::vector<solvent_accessible_atom_t> &sa) const;
   std::string get_atom_name(const clipper::Coord_orth &pt, CMMDBManager *mol) const;

public:
   widgeted_molecule_t() { init(); }
   widgeted_molecule_t(const lig_build::molfile_molecule_t &mol_in, CMMDBManager *pdb_mol);

   // return 0 as first if not highlighting a bond
   std::pair<bool, widgeted_bond_t> highlighted_bond_p(int x, int y) const;

   // return -1 as the atom index if not highlighting an atom.
   std::pair<int, widgeted_atom_t> highlighted_atom_p(int x, int y) const;
   bool write_mdl_molfile(const std::string &file_name) const;
   bool write_minimal_cif_file(const std::string &file_name) const;
   bool close_bond(int ib, GooCanvasItem *root, bool handle_post_delete_stray_atoms_flag);
   bool close_atom(int iat, GooCanvasItem *root);
   std::vector<int> get_unconnected_atoms() const;

   // don't count closed bonds.
   std::vector<int> bonds_having_atom_with_atom_index(int test_atom_index) const;
   bool operator==(const widgeted_molecule_t &mol_other) const;
   int n_stray_atoms() const; // unbonded atoms
   std::vector<int> stray_atoms() const;
   void translate(const lig_build::pos_t &delta); // move the atoms
   lig_build::pos_t input_coords_to_canvas_coords(const clipper::Coord_orth &in) const;
      

   // make private when bug is fixed.
   lig_build::pos_t centre_correction;
   std::pair<bool, double> scale_correction;
   double mol_in_min_y;
   double mol_in_max_y;

   void map_solvent_accessibilities_to_atoms(std::vector<solvent_accessible_atom_t> solvent_accessible_atoms);

   // can throw an exception
   // 
   lig_build::pos_t get_atom_canvas_position(const std::string &atom_name);
   
};
