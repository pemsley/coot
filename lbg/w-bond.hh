// ====================================================================
//                     widgeted_bond_t
// ====================================================================

#ifndef W_BOND_HH
#define W_BOND_HH

#include <gtk/gtk.h>
#include <goocanvas.h>
#include "lidia-core/lig-build.hh"


static gboolean
on_wmolecule_key_press_event (GooCanvasItem *item,
			      GooCanvasItem *target,
			      GdkEventKey *event,
			      gpointer data)
{
  gchar *id = 0;
  // id = g_object_get_data (G_OBJECT (item), "id");
  g_print ("%s received key-press event\n", id ? id : "unknown");
  return FALSE;
}


class widgeted_bond_t : public lig_build::bond_t, ligand_layout_graphic_primitives {
   GooCanvasItem *ci;
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      } else {
	 std::cout << "WARNING:: widgeted_bond_t::clear() failed to remove child "
		   << ci << " for root " << root << std::endl;
      }
      ci = NULL;
   }

   void construct_internal(const lig_build::atom_t &atom_first,
			   const lig_build::atom_t &atom_second,
			   bool shorten_first,
			   bool shorten_second,
			   bond_type_t bt,
			   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
			   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
			   GooCanvasItem *root) {
      
// Now CH3 are superatoms that require shortening.      
//       bool shorten_first = 0;
//       bool shorten_second = 0;
//       if (atom_first.atom_id != "C") { 
// 	 shorten_first = 1;
//       } 
//       if (atom_second.atom_id != "C") { 
// 	 shorten_second = 1;
//       }

      ci = canvas_item_for_bond(atom_first, atom_second, shorten_first, shorten_second, bt,
				other_connections_to_first_atom,
				other_connections_to_second_atom,
				root);

      if (false)
	 std::cout << "construct_internal() shorten first: " << shorten_first
		   << " shorten_second " << shorten_second << " made ci " << ci << std::endl;
      
      g_signal_connect (ci, "key_press_event",
			G_CALLBACK (on_wmolecule_key_press_event), NULL);
      
   }

   // all bonds are made this way...
   // 
   GooCanvasItem *canvas_item_for_bond(const lig_build::atom_t &at_1,
				       const lig_build::atom_t &at_2,
				       bool shorten_first,
				       bool shorten_second,
				       bond_type_t bt,
				       const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
				       const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
				       GooCanvasItem *root) const;

   // We need to make a shorter bond canvas line because we have (say)
   // changed a carbon to a N (bond canvas line now does not
   // completely extend to the atom position).
   // 
   void make_new_canvas_item_given_type(const lig_build::atom_t &atom_changed,
					const lig_build::atom_t &atom_other,
					bool shorten_first,
					bool shorten_second,
					lig_build::bond_t::bond_type_t bt,
					const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
					const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
					GooCanvasItem *root) {

      lig_build::pos_t A = atom_changed.atom_position;
      lig_build::pos_t B =   atom_other.atom_position;

//       bool shorten_first = 0;
//       bool shorten_second = 0;
//       if (atom_changed.atom_id != "C")
// 	 shorten_first = 1;
//       if (atom_other.atom_id != "C")
// 	 shorten_second = 1;
      
      GooCanvasItem *new_line = canvas_item_for_bond(atom_changed, atom_other,
						     shorten_first, shorten_second,
						     bt,
						     other_connections_to_first_atom,
						     other_connections_to_second_atom,
						     root);
      update_canvas_item(new_line, root);
   }

   GooCanvasItem * canvas_item_double_bond_simple(const lig_build::pos_t &pos_1,
						  const lig_build::pos_t &pos_2,
						  GooCanvasItem *root) const;
					
   GooCanvasItem * canvas_item_double_bond(const lig_build::pos_t &pos_1,
					   const lig_build::pos_t &pos_2,
					   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
					   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
					   GooCanvasItem *root) const;
					
   GooCanvasItem * canvas_item_double_with_shortened_side_bond(const lig_build::pos_t &pos_1,
							       const lig_build::pos_t &pos_2,
							       GooCanvasItem *root) const;
					
   GooCanvasItem * make_wedge_bond_item(const lig_build::pos_t &pos_1,
					const lig_build::pos_t &pos_2,
					const lig_build::bond_t::bond_type_t &bt,
					const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
					GooCanvasItem *root) const;
   GooCanvasItem * make_wedge_out_bond_item(const lig_build::pos_t &pos_1,
					    const lig_build::pos_t &pos_2,
					    const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
					    GooCanvasItem *root) const;
   // which uses
   GooCanvasItem *make_sheared_or_darted_wedge_bond(const lig_build::pos_t &pos_1,
					    const lig_build::pos_t &pos_2,
					    const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
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

   ~widgeted_bond_t() {}

   // Now we use a constructor that does the creation of the canvas item too
   //
   widgeted_bond_t(int first, int second, 
		   const lig_build::atom_t &atom_first,
		   const lig_build::atom_t &atom_second,
		   bool shorten_first, bool shorten_second,
		   bond_type_t bt,
		   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
		   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
		   GooCanvasItem *root) :
      lig_build::bond_t(first, second, bt) {
      construct_internal(atom_first, atom_second, shorten_first, shorten_second, bt,
			 other_connections_to_first_atom,
			 other_connections_to_second_atom, root);
   }
   // as above, but we give the centre of the ring too.
   widgeted_bond_t(int first, int second, 
		   const lig_build::atom_t &atom_first,
		   const lig_build::atom_t &atom_second,
		   bool shorten_first, bool shorten_second,
		   lig_build::pos_t centre_pos_in,
		   bond_type_t bt,
		   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
		   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
		   GooCanvasItem *root) :
      bond_t(first, second, centre_pos_in, bt) {
      
      // because there was a ring, these atoms don't need CH3 shortening (so we
      // can rely on the atom type)
      // 20160909-PE that is no longer true now that we have atom names - so all
      // bonds are shortened.
      // 
      // bool shorten_first = false;
      // bool shorten_second = false;
      
      construct_internal(atom_first, atom_second, shorten_first, shorten_second, bt,
			 other_connections_to_first_atom,
			 other_connections_to_second_atom, root);
   }
   
   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      clear(root);
      ci = new_item;
   }

   void update(const lig_build::atom_t &at_1, const lig_build::atom_t &at_2,
	       bool shorten_first, bool shorten_second,
	       const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
	       const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
	       GooCanvasItem *root) {
      
      clear(root);
      ci = canvas_item_for_bond(at_1, at_2, shorten_first, shorten_second, get_bond_type(),
				other_connections_to_first_atom,
				other_connections_to_second_atom, root);
   }

   // for cases when we know that we don't have a wedge bond - or don't care if we do.
   void update(const lig_build::atom_t &at_1, const lig_build::atom_t &at_2,
	       bool shorten_first, bool shorten_second, GooCanvasItem *root) {
      clear(root);
      std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_first_atom;
      std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > other_connections_to_second_atom;
      ci = canvas_item_for_bond(at_1, at_2, shorten_first, shorten_second, get_bond_type(),
				other_connections_to_first_atom,
				other_connections_to_second_atom,
				root);
   }

   void rotate_canvas_item(gdouble cx, gdouble cy, gdouble degrees) {
      if (ci)
	 wrap_goo_canvas_item_rotate(ci, degrees, cx, cy);
      else
	 std::cout << "ERROR: traped null ci in wmolecule rotate_canvas_item() "
		   << std::endl;
   }

   // We need to make a shorter bond canvas line because we have (say)
   // changed a carbon to a N (bond canvas line now does not
   // completely extend to the atom position).
   void make_new_canvas_item(const lig_build::atom_t &atom_changed,
			     const lig_build::atom_t &atom_other,
			     bool shorten_first,
			     bool shorten_second,
			     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
			     const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
			     GooCanvasItem *root) {

      lig_build::bond_t::bond_type_t bt = get_bond_type();
      make_new_canvas_item_given_type(atom_changed, atom_other,
				      shorten_first, shorten_second, bt,
				      other_connections_to_first_atom,
				      other_connections_to_second_atom, root);
   }
   void change_bond_order(const lig_build::atom_t &atom_changed,
			  const lig_build::atom_t &atom_other,
			  bool shorten_first, bool shorten_second,
			  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
			  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
			  GooCanvasItem *root) {
      change_bond_order(atom_changed, atom_other, shorten_first, shorten_second, 0,
			other_connections_to_first_atom,
			other_connections_to_second_atom, root);
   }
   void change_bond_order(const lig_build::atom_t &atom_changed,
			  const lig_build::atom_t &atom_other,
			  bool shorten_first, bool shorten_second,
			  bool allow_triple_toggle,
			  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,			  
			  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,			  
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
	 if (bt == lig_build::bond_t::DOUBLE_BOND) { 
	    if (allow_triple_toggle)
	       bt = lig_build::bond_t::TRIPLE_BOND;
	    else
	       bt = lig_build::bond_t::SINGLE_BOND;
	 } else {
	    if (bt == lig_build::bond_t::TRIPLE_BOND) { 
	       bt = lig_build::bond_t::DOUBLE_BOND;
	    } else { 
	       if (bt == lig_build::bond_t::IN_BOND) {
		  std::swap(atom_1, atom_2);
		  std::swap(at_1, at_2);
		  bt = lig_build::bond_t::OUT_BOND;
	       } else {
		  if (bt == lig_build::bond_t::OUT_BOND) { 
		     bt = lig_build::bond_t::IN_BOND;
		  }
	       }
	    }
	 }
      }

      if (0)
	 std::cout << "changing bond type from " << get_bond_type() << " to "
		   << bt << std::endl;
      
      set_bond_type(bt);
      make_new_canvas_item_given_type(at_1, at_2, shorten_first, shorten_second, bt,
				      other_connections_to_first_atom,
				      other_connections_to_second_atom, root);
   }
   void close(GooCanvasItem *root) {
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
      case AROMATIC_BOND:
      case DELOC_ONE_AND_A_HALF:
      case SINGLE_OR_DOUBLE:
      case SINGLE_OR_AROMATIC:
      case DOUBLE_OR_AROMATIC:
      case BOND_ANY:
      case BOND_UNDEFINED:
	 mmdb_bt = UNASSIGNED_INDEX;
	 break;
      }
      return mmdb_bt;
   }

   // debugging
   GooCanvasItem *get_ci() const { return ci;}


}; // end of widgeted_bond_t



#endif // W_BOND_HH

