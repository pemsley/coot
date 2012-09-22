 /* src/main.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2010 by The University of Oxford.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef RAMA_PLOT_HH
#define RAMA_PLOT_HH

#include <map>
#include <string>
#include "coot-coord-utils.hh"

// ------------------------ Canvas stuff -----------------------------------
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
#ifdef HAVE_GTK_CANVAS

#include <gtk/gtk.h>
#include <gdk_imlib.h>
#include <gtk-canvas.h>

#else 
#ifdef HAVE_GNOME_CANVAS	// usually defined for WINDOWS_MINGW
#include <libgnomecanvas/libgnomecanvas.h>
// We use GtkCanvas in the code, so for GNOME Canvas that's GnomeCanvas
#include <gtk/gtk.h>
typedef GnomeCanvas     GtkCanvas;
typedef GnomeCanvasItem GtkCanvasItem;
typedef GnomeCanvasPoints GtkCanvasPoints;
#endif
#endif // HAVE_GNOME_CANVAS
// ------------------------ End of Canvas stuff ------------------------------


#include <mmdb/mmdb_manager.h>
#include "clipper/core/ramachandran.h"
#include "coot-rama.hh"

namespace coot { 

   class canvas_tick_t {
   public:
      double x, y;
      short int axis;
      // i==0 for x axis, i==1 for y axis.
      canvas_tick_t(int i, double a, double b) {
	 x = a;
	 y = b;
	 axis = i;
      }
      double start_x() {return x;};
      double start_y() {return y;};
      double end_x() { if (axis == 0) return x; else return x - 10;};
      double end_y() { if (axis == 1) return y; else return y + 10;}; 
   }; 
   
   class mouse_util_t {
   public:
      enum {MODEL_NUMBER_UNSET = -1}; 
      mouse_util_t() {
	 model_number = -1; // unset
	 mouse_over_secondary_set = 0; } // use residue_spec_t constructor to signal unset.
      residue_spec_t spec;
      int model_number;
      bool mouse_over_secondary_set;  // used in mouse_point_check()
                                      // and the routine which pops up
                                      // the mouseover label.
   };
   

// class rama_matcher_score {

// public:
//    int i;
//    double dist;
//    rama_matcher_score(int in, double distn) {
//       i = in;
//       dist = distn; 
//    } 

// };



   // the container class with the residue spec -> phi-psi map (associative array)
   //
   // generating from a CMMDBManager makes a vector of these (one for each model).
   // 
   class phi_psis_for_model_t {

   public:
      int model_number;
      std::map<residue_spec_t, util::phi_psi_t> phi_psi;
      phi_psis_for_model_t(int model_number_in) {
	 model_number = model_number_in;
      }
      void add_phi_psi(const residue_spec_t &spec, const util::phi_psi_t &phi_psi_in) {
	 phi_psi[spec] = phi_psi_in;
      }
      util::phi_psi_t operator[](const residue_spec_t &spec) {
 	 return phi_psi[spec];
      }
      unsigned int size() { return phi_psi.size(); } 
   };

   class diff_sq_t {
      
      double v_; 

      util::phi_psi_t pp_1;
      util::phi_psi_t pp_2;
      residue_spec_t res_1;
      residue_spec_t res_2;
      
   public:
      diff_sq_t() {}; 
      diff_sq_t(const util::phi_psi_t &pp_1_in, const util::phi_psi_t &pp_2_in,
		const residue_spec_t &r1, const residue_spec_t r2,
		double d) {
	 pp_1 = pp_1_in;
	 pp_2 = pp_2_in;
	 res_1 = r1;
	 res_2 = r2;
	 v_ = d;
      }
      util::phi_psi_t phi_psi_1() const { return pp_1; } 
      util::phi_psi_t phi_psi_2() const { return pp_2; }
      double v() const { return v_;}
   };

   bool compare_phi_psi_diffs(const diff_sq_t &d1, const diff_sq_t d2);

   // A container for information for testing whether a kleywegt
   // phi/psi pair goes over the board of the phi/phi plot (if it
   // should go over the border () we say is_wrapped = 1.
   //
   // If it doesn't go over the border, we draw the conventional
   // (potentially long) arrows.
   // 
   class rama_kleywegt_wrap_info {
   public:
      std::pair<float, float> primary_border_point;
      std::pair<float, float> secondary_border_point;
      short int is_wrapped;
      rama_kleywegt_wrap_info() {
	 is_wrapped = 0;
      }
   };

   class rama_stats_container_t {
   public:
      int n_ramas;
      int n_preferred;
      int n_allowed;
      rama_stats_container_t() {
	 n_ramas = 0;
	 n_preferred = 0;
	 n_allowed = 0;
      }
      void operator+=(const rama_stats_container_t &sc) {
	 n_ramas += sc.n_ramas;
	 n_preferred += sc.n_preferred;
	 n_allowed += sc.n_allowed;
      }
   };

class rama_plot {

   int imol;  // which molecule in mapview did this come from?
   clipper::Ramachandran rama;
   clipper::Ramachandran r_gly, r_pro, r_non_gly_pro;
   GtkCanvas *canvas;
   GtkWidget *dialog; // the container for the canvas
   std::vector<GtkCanvasItem *> canvas_item_vec; // we save them so that
					         // we can destroy them.
   GtkCanvasItem *big_box_item; 
   float step; // the "angular" size of the background blocks
   std::vector<int> ifirst_res; // offset between actual residue number and
		                // position in the phi_psi vector

   // a set of phi-psis for each model.
   std::vector<phi_psis_for_model_t> phi_psi_model_sets;
   std::vector<phi_psis_for_model_t> secondary_phi_psi_model_sets;
   GtkTooltips *tooltips;
   std::string fixed_font_str;

   void setup_internal(float level_prefered, float level_allowed);
   // gint rama_button_press(GtkWidget *widget, GdkEventButton *event);
   
   GtkCanvasItem *tooltip_item;
   GtkCanvasItem *tooltip_item_text;
   clipper::Ramachandran::TYPE displayed_rama_type;
   short int drawing_differences;
   int n_diffs; // number of differences in the 2 molecule ramachandran plot
                // typically 10.
   std::vector<diff_sq_t> diff_sq;
   double zoom;  // initial 0.8
   short int have_sticky_labels; // when we mouse away from a residue
				 // that has a tooltip-like label,
				 // should the label disappear?  If
				 // no, then we have sticky labels.

   double rama_threshold_allowed;   // 0.05 
   double rama_threshold_preferred; // 0.002
   void init_internal(const std::string &mol_name,
		      float level_prefered, float level_allowed,
		      float block_size,
		      short int hide_butttons = 0,
		      short int is_kleywegt_plot = 0); // called by init(int imol)
   void draw_green_box(double phi, double psi);
   short int phipsi_edit_flag;   // for active canvas (can move phi/psi point)
   short int backbone_edit_flag; // for passive canvas
   void clear_canvas_items();
   void clear_last_canvas_item();
   void clear_last_canvas_items(int n);
   std::pair<int, int> molecule_numbers_; // needed for undating kleywegt plots
   std::pair<std::string, std::string> chain_ids_; // ditto.
   int dialog_position_x; 
   int dialog_position_y;
   bool kleywegt_plot_uses_chain_ids;
   void hide_stats_frame();
   void counts_to_stats_frame(const rama_stats_container_t &sc);

   
   bool allow_seqnum_offset_flag; // was from a shelx molecule with A 1->100 and B 201->300
   int seqnum_offset; // for shelx molecule as above, what do we need to add to seqnum_1 to get the
                      // corresponding residue in the B chain (in the above example it is 100).
   int get_seqnum_2(int seqnum_1) const;
   void set_seqnum_offset(int imol1, int imol2,
			  CMMDBManager *mol1,
			  CMMDBManager *mol2,
			  const std::string &chain_id_1,
			  const std::string &chain_id_2);
   util::phi_psi_t green_box;
   void draw_green_box(); // use above stored position to draw green square
   bool green_box_is_sensible(util::phi_psi_t gb) const; // have the phi and psi been set to
                                                         // something sensible?
   void recentre_graphics_maybe(mouse_util_t t);
   mouse_util_t mouse_point_check_differences(double worldx, double worldy) const;

   void find_phi_psi_differences();
   void find_phi_psi_differences_internal(const std::string &chain_id1,
					  const std::string &chain_id2,
					  bool use_chain_ids);
   void find_phi_psi_differences(const std::string &chain_id1,
				 const std::string &chain_id2);

   mouse_util_t mouse_point_check_internal(const coot::phi_psis_for_model_t &phi_psi_set,
					   int imod, 
					   double worldx, double worldy,
					   bool is_secondary) const;

   
public:

   enum rama_position_t {RAMA_OUTLIER, RAMA_ALLOWED, RAMA_PREFERRED, RAMA_UNKNOWN};
   void resize_rama_canvas_internal(GtkWidget *widget, GdkEventConfigure *event);
   rama_plot() {
      big_box_item = 0; dialog = 0;
      dialog_position_x = -100; dialog_position_y = -100; };

   // consider destructor where we should
   // gtk_object_destroy(big_box_item) if it is non-zero.
   void init(const std::string &type);
   // typically level_prefered = 0.02, level_allowed is 0.002, block_size is 10.0;
   void init(int imol_no, const std::string &mol_name, float level_prefered, float level_allowed,
	     float block_size_for_background, short int is_kleywegt_plot_flag);

   void allow_seqnum_offset();
   void set_n_diffs(int nd);


   // The graphics interface, given that you have a CMMDBManager. 
   // 
   void draw_it(CMMDBManager *mol);
   void draw_it(int imol1, int imol2, CMMDBManager *mol1, CMMDBManager *mol2); // no chain ids.
   void draw_it(int imol1, int imol2,
		CMMDBManager *mol1, CMMDBManager *mol2,
		const std::string &chain_id_1, const std::string &chain_id_2);
   
   void draw_it(const util::phi_psi_t &phipsi);
   void draw_it(const std::vector<util::phi_psi_t> &phipsi);

   int molecule_number() { return imol; } // mapview interface

   void basic_white_underlay(); // Not const because we modify canvas_item_vec. 
   void display_background();   // Likewise.
   void setup_canvas(); 
   void black_border();
   void cell_border(int i, int j, int step);

   void generate_phi_psis(CMMDBManager *mol);
   void generate_phi_psis(CMMDBManager *mol, bool primary_flag);
   void generate_secondary_phi_psis(CMMDBManager *mol);
   void generate_phi_psis_by_selection(CMMDBManager *mol,
				       bool is_primary_flag,
				       int SelectionHandle);
				       

   
   void generate_phi_psis(CMMDBManager *mol, int SelectionHandle);
   
   void generate_phi_psis_debug();

   void draw_axes();
   void draw_zero_lines();

   int draw_phi_psi_point(const util::phi_psi_t &phi_psi, bool as_white_flag);
   void draw_kleywegt_arrow(const util::phi_psi_t &phi_psi_primary,
			    const util::phi_psi_t &phi_psi_secondary,
			    GtkCanvasPoints *points);
   rama_kleywegt_wrap_info test_kleywegt_wrap(const util::phi_psi_t &phi_psi_primary,
					      const util::phi_psi_t &phi_psi_secondary)
      const;
						    
   int draw_phi_psi_point_internal(const util::phi_psi_t &phi_psi,
				   bool as_white_flag, int box_size);

   
   rama_stats_container_t draw_phi_psi_points();
   rama_stats_container_t draw_phi_psi_points_for_model(const coot::phi_psis_for_model_t &pp_set); 


   // sanitise the arguments
   void big_square(const std::string &chainid,
		   int ires,
		   const std::string &ins_code);
   void big_square(int model_number, 
		   const std::string &chainid,
		   int ires,
		   const std::string &ins_code);


   // Not sure about this yet
//    void add_phi_psi(std::vector <util::phi_psi_t> *phi_psi_vec,
// 		    CMMDBManager *mol, const char *chain_id,
// 		    CResidue *prev, CResidue *this_res, CResidue *next_res);

   util::phi_psi_pair_helper_t make_phi_psi_pair(CMMDBManager *mol1,
						 CMMDBManager *mol2,
						 const std::string &chain_id1,
						 const std::string &chain_id2,
						 int resno,
						 const std::string &ins_code) const;

   // SelResidue is guaranteed to have 3 residues (there is no protection
   // for that in this function).
   std::pair<bool, util::phi_psi_t> get_phi_psi(PCResidue *SelResidue) const;

   void map_mouse_pos(double x, double y);
   void mouse_motion_notify(GdkEventMotion *event, double x, double y);
   void mouse_motion_notify_editphipsi(GdkEventMotion *event, double x, double y);
   gint button_press (GtkWidget *widget, GdkEventButton *event); 
   gint button_press_conventional (GtkWidget *widget, GdkEventButton *event); 
   gint button_press_editphipsi (GtkWidget *widget, GdkEventButton *event); 
   gint button_press_backbone_edit (GtkWidget *widget, GdkEventButton *event); 

   void draw_phi_psis_on_canvas(char *filename);
   // void update_canvas(CMMDBManager *mol);  // mol was updated(?)

   void draw_2_phi_psi_sets_on_canvas(char *file1, char *file2);
   void draw_2_phi_psi_sets_on_canvas(CMMDBManager *mol1, 
				      CMMDBManager *mol2);
   void draw_2_phi_psi_sets_on_canvas(CMMDBManager *mol1, 
				      CMMDBManager *mol2,
				      std::string chain_id1, std::string chain_id2);

   void draw_phi_psi_differences(); 

   CMMDBManager *rama_get_mmdb_manager(std::string pdb_name);
   // void tooltip_like_box(int i, short int is_seconary);
   void tooltip_like_box(const mouse_util_t &t); 
   int draw_phi_psi_as_gly(const util::phi_psi_t &phi_psi); 
   void residue_type_background_as(std::string res);
   void all_plot(clipper::Ramachandran::TYPE type);

   // provides a object that contains the residue number and the chain
   // number (for phi_psi_set indexing) of the mouse position
   // 
   mouse_util_t mouse_point_check(double x, double y) const;

   void all_plot_background_and_bits(clipper::Ramachandran::TYPE type); 

   gint key_release_event(GtkWidget *widget, GdkEventKey *event);

   // ramachandran difference plot
   short int is_kleywegt_plot() const {
      return drawing_differences;
   }

   // where we also say where to put the ramachandran dialog:
   // Call this function before init() :-)
   void set_position(int x_position, int y_position) {
      dialog_position_x = x_position;
      dialog_position_y = y_position;
   }

   std::pair<int, int> molecule_numbers() const { // needed for undating kleywegt plots
      return molecule_numbers_;
   }
   std::pair<std::string, std::string> chain_ids() const { 
      return chain_ids_;
   }

   void zoom_in(); 
   void zoom_out();

   void set_kleywegt_plot_uses_chain_ids() { kleywegt_plot_uses_chain_ids = 1; }
   bool kleywegt_plot_uses_chain_ids_p() { return kleywegt_plot_uses_chain_ids; }

   unsigned int n_phi_psi_model_sets() const { return phi_psi_model_sets.size(); }

   phi_psis_for_model_t get_phi_psis_for_model(unsigned int model_no) const {
      phi_psis_for_model_t r(-1);
      if (model_no < phi_psi_model_sets.size())
	 r = phi_psi_model_sets[model_no];
      return r;
   } 

   void debug() const; 

   void destroy_yourself();

}; 

} // namespace coot

#endif // HAVE_GTK_CANVAS

#endif //  RAMA_PLOT_HH
