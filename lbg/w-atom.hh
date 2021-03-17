// ====================================================================
//                     widgeted_atom_t
// ====================================================================

#ifndef W_ATOM_HH
#define W_ATOM_HH

#include <goocanvas.h>
#include "lidia-core/lig-build.hh"


class widgeted_atom_t : public lig_build::atom_t , ligand_layout_graphic_primitives {
   std::string font_colour;
   double solvent_accessibility;
   GooCanvasItem *ci;
   // std::string atom_name; // typically names from a PDB file. 20111229 base class now
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      }
      ci = NULL;
   }
   void blank_text(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 // set text to ""
	 g_object_set(ci, "visibility", GOO_CANVAS_ITEM_INVISIBLE, NULL);
	 // gtk_widget_hide(GTK_WIDGET(ci)); invalid cast
      } else {
	 // std::cout << "debug::       blank_text() ci not found in children of root "<< std::endl;
      } 
   }

   // then general form, use this not the below 2 (which are used by
   // this function)
   // 
   GooCanvasItem *make_canvas_text_item(const lig_build::atom_id_info_t &atom_id_info_in,
					const std::string &fc,
					GooCanvasItem *root);

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
   ~widgeted_atom_t() {}
   GooCanvasItem *get_canvas_item() const { return ci; }
   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      blank_text(root);
      ci = new_item;
   }
   bool update_atom_id_maybe(const lig_build::atom_id_info_t &atom_id_info_in,
			     GooCanvasItem *root) {
      return update_atom_id_maybe(atom_id_info_in, font_colour, root);
   }
   bool update_atom_id_maybe(const lig_build::atom_id_info_t &atom_id_info_in,
			     const std::string &fc,
			     GooCanvasItem *root) {
      bool changed_status = 0;
      font_colour = fc;
      GooCanvasItem *text_item = NULL;
      if (! is_closed()) {
	 if (atom_id_info_in.atom_id != get_atom_id()) {
	    changed_status = set_atom_id(atom_id_info_in.atom_id);
	    if (changed_status) {
	       text_item = make_canvas_text_item(atom_id_info_in, fc, root);
	       update_canvas_item(text_item, root);
	    }
	 }
      } else {
	 update_canvas_item(text_item, root); // close atom, replace with null.
      } 
      return changed_status;
   }
   void update_atom_id_forced(const lig_build::atom_id_info_t &atom_id_info_in,
			      const std::string &fc, 
			      GooCanvasItem *root) {

      set_atom_id(atom_id_info_in.atom_id);
      GooCanvasItem *text_item = make_canvas_text_item(atom_id_info_in, fc, root);
      update_canvas_item(text_item, root);
   }
   
   void update_atom_id_forced(const lig_build::atom_id_info_t &atom_id_info_in,
			      GooCanvasItem *root) {
      update_atom_id_forced(atom_id_info_in, font_colour, root);
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
//    std::string get_atom_name() const {
//       return atom_name;
//    } 
   std::vector<coot::bash_distance_t> bash_distances;
};



#endif // W_ATOM_HH

