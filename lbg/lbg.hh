/* lbg/lbg.hh
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
 * Copyright 2013, 2015, 2016 by Medical Research Council
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
n * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef LBG_HH
#define LBG_HH

// Don't compile any of this if we don't have the canvas.


#ifdef USE_PYTHON
#   ifndef HAVE_INCLUDED_PYTHON
#      define HAVE_INCLUDED_PYTHON
#      include <Python.h>
#   endif
#endif 


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#endif 

#include <iostream>
#include <map>
#include <queue>

#include <gtk/gtk.h>
#include <goocanvas.h>

#include "gsl/gsl_multimin.h"

#include <mmdb2/mmdb_manager.h>
#ifndef MONOMER_DIR_STR
#define MONOMER_DIR_STR "COOT_CCP4SRS_DIR"
#endif 

#include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-molfile.hh"
#include "utils/coot-utils.hh"

#include "wmolecule.hh"

#include "solvent-exposure-difference.hh"
#include "pli/flev-annotations.hh"
#include "pli/pi-stacking.hh"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "graphics-c-interface-functions.hh"
#endif 

// static double LIGAND_TO_CANVAS_SCALE_FACTOR = 23;
static double LIGAND_TO_CANVAS_SCALE_FACTOR = 17;
static double SINGLE_BOND_CANVAS_LENGTH= LIGAND_TO_CANVAS_SCALE_FACTOR * 1.54;


bool save_togglebutton_widgets(GtkBuilder *builder);

void lbg_handle_toggle_button(GtkToggleToolButton *tb, GtkWidget *canvas, int mode);
GtkWidget *get_canvas_from_scrolled_win(GtkWidget *scrolled_window);


void lbg_scale_adj_changed( GtkWidget *widget, GtkSpinButton *spinbutton);


// extern "C" { 
// static gboolean on_residue_circle_clicked(GooCanvasItem  *item,
// 					  GooCanvasItem  *target_item,
// 					  GdkEventButton *event,
// 					  gpointer        user_data);
// }

namespace coot {
   mmdb::Residue *get_first_residue_helper_fn(mmdb::Manager *mol);
}

// ====================================================================
//                     lbg_info_t
// ====================================================================

class lbg_info_t {

public:

   enum { UNASSIGNED_INDEX = -1 };
      
   class highlight_data_t {
      int n_atoms_;
      lig_build::pos_t pos_1_;
      lig_build::pos_t pos_2_;
      std::pair<int,int> bond_indices; // the atom indices of the bond.
      int atom_index; // the index of the atom in the molecule (single atom highlighting)
      bool has_ring_centre_flag;
      lig_build::pos_t ring_centre;
      GooCanvasItem *highlight_widget;
      lig_build::polygon_position_info_t get_new_polygon_centre_using_2_atoms(int n_edges,
								const double &radius) const;
      lig_build::polygon_position_info_t
      get_new_polygon_centre_using_1_atom(int n_edges,
					  bool spiro,
					  const double &radius_std,
					  const double &radius_corr,
					  const widgeted_molecule_t &mol) const;
   public:
      highlight_data_t(GooCanvasItem *w_in,
		       std::pair<int, int> bond_indices_in,
		       const lig_build::pos_t &p1,
		       const lig_build::pos_t &p2) {
	 n_atoms_ = 2;
	 pos_1_ = p1;
	 pos_2_ = p2;
	 highlight_widget = w_in;
	 atom_index = UNASSIGNED_INDEX; // no value
	 bond_indices = bond_indices_in;
	 has_ring_centre_flag = false;
      }
      highlight_data_t(GooCanvasItem *w_in,
		       const lig_build::pos_t &p, int index_in) {
	 n_atoms_ = 1;
	 pos_1_ = p;
	 highlight_widget = w_in;
	 atom_index = index_in;
	 bond_indices = std::pair<int, int> (UNASSIGNED_INDEX, UNASSIGNED_INDEX);
	 has_ring_centre_flag = false;
      }

      highlight_data_t() {
	 highlight_widget = NULL;
	 n_atoms_ = 0;
	 atom_index = UNASSIGNED_INDEX; // unset
	 bond_indices = std::pair<int, int> (UNASSIGNED_INDEX, UNASSIGNED_INDEX);
	 has_ring_centre_flag = false;
      }
      bool has_contents() const {
	 if (highlight_widget)
	    return 1;
	 else
	    return 0;
      }
      bool single_atom() const {
	 return (n_atoms_ == 1);
      }
      int get_atom_index() const { return atom_index; }
      lig_build::pos_t get_atom_1_pos() const {
	 return pos_1_;
      } 
      void clear(GooCanvasItem *root) {
	 n_atoms_ = 0;
	 if (highlight_widget) {
	    gint child_index = goo_canvas_item_find_child(root, highlight_widget);
	    if (child_index != -1) {
	       goo_canvas_item_remove_child(root, child_index);
	       highlight_widget = NULL;
	    }
	 } else {
	    std::cout << "in clear() NULL highlight_widget" << std::endl;
	 }
      }

      bool has_ring_centre() const { return has_ring_centre_flag; }
      
      // The old way was to use clear() where we lose the button up if
      // the button up happens on the highlight item. So instead, we
      // simply undisplay the widget.
      //
      // Memory loss, I think but better than missing a button
      // release.
      //
      void undisplay() {
	 n_atoms_ = 0;
	 if (highlight_widget) {
	    g_object_set(highlight_widget, "visibility", GOO_CANVAS_ITEM_INVISIBLE, NULL);
	    highlight_widget = NULL; // eek.
	 } else {
	    std::cout << "in undisplay() NULL highlight_widget" << std::endl;
	 }
      }

      void set_ring_centre(const lig_build::pos_t &pos) {
	 has_ring_centre_flag = true;
	 ring_centre = pos;
      } 

      
      lig_build::polygon_position_info_t
      get_new_polygon_centre(int n_edges,
			     bool spiro_flag,
			     const double &radius_standard,
			     const double &radius_corrected,
			     const widgeted_molecule_t &mol) const;
      std::pair<int, int> get_bond_indices() const { return bond_indices; }
      // behind this class's back the atom indices of the bond can get
      // swapped.  Here is a means to tell that to this class.
      void swap_bond_indices() { std::swap(bond_indices.first, bond_indices.second);}
   }; // finish highlight_data_t class


   // contains a coordinate on the flat ligand (input mol coords) that
   // corresponds to an atom on the ligand - when we add read the
   // ligand mol file, we'll attach to an atom that it needs a bond to
   // the residue (that encapsulates this bond description).
   // 
   class bond_to_ligand_t {
   private:
      bool is_set_; 
   public:
      // sync to c-interface-ligands.hh
      enum { H_BOND_DONOR_MAINCHAIN,
	     H_BOND_DONOR_SIDECHAIN,
	     H_BOND_ACCEPTOR_MAINCHAIN, 
	     H_BOND_ACCEPTOR_SIDECHAIN,
	     METAL_CONTACT_BOND,
	     BOND_COVALENT,	     
	     BOND_OTHER };  
      std::string ligand_atom_name;
      double bond_length;
      int bond_type; // acceptor or donor
      bond_to_ligand_t(const std::string &n, double b) {
	 ligand_atom_name = n;
	 bond_length = b;
	 is_set_ = 1;
      }
      bond_to_ligand_t() { is_set_ = 0;}
      bool is_set() const { return is_set_; } 
   };

   class residue_circle_t {
      bool   se_diff_set_;
      double se_holo; 
      double se_apo;
      int stacking_type;
      std::vector<std::string> ligand_ring_atom_names;
      std::string ligand_cation_atom_name;
      bool is_a_primary_residue_; // primary residues are stacking or
				  // bond and are placed first looking
				  // for non-overlaping
				  // interaction-bond lines.
   public:
      // CATION_PI_STACKING sets ligand_cationic_atom_name, not the
      // ligand_ring_atom_names vector.
      //
      enum { PI_PI_STACKING,
	     PI_CATION_STACKING, // for cations on the protein residues (ligand pi)
	     CATION_PI_STACKING, // for cations on the ligand (protein TRY, PRO, TRP)
      };

      clipper::Coord_orth trans_rel_pos_3d;
      clipper::Coord_orth residue_centre_real;
      
      lig_build::pos_t pos; // coordinate system of the ligand atoms
      coot::residue_spec_t spec;
      std::string residue_type;
      std::string residue_label;
      std::vector<bond_to_ligand_t> bonds_to_ligand;
      double water_dist_to_protein; 
      residue_circle_t(const clipper::Coord_orth &pos_in,
		       const clipper::Coord_orth &click_pos_in,
		       coot::residue_spec_t spec_in,
		       const std::string &type_in,
		       const std::string &label_in) {
	 trans_rel_pos_3d = pos_in;
	 residue_centre_real = click_pos_in;
	 spec = spec_in;
	 residue_type = type_in;
	 residue_label = label_in;
	 se_holo = 0.0;
	 se_apo  = 0.0;
	 se_diff_set_ = 0;
	 stacking_type = -1; // should have enumerated value
	 is_a_primary_residue_ = 0;
	 water_dist_to_protein = 100;
      }
      void set_canvas_pos(const lig_build::pos_t &pos_in) {
	 pos = pos_in;
      }
      void add_bond_to_ligand(const bond_to_ligand_t &bl) {
	 bonds_to_ligand.push_back(bl);
	 is_a_primary_residue_ = 1;
      }
      void set_solvent_exposure_diff(double se_holo_in, double se_apo_in) {
	 se_holo = se_holo_in;
	 se_apo = se_apo_in;
	 se_diff_set_ = 1;
      }
      void set_stacking(const std::string &type_string,
			const std::vector<std::string> &ligand_ring_atom_names_in,
			const std::string &ligand_cation_atom_name_in) {
	 is_a_primary_residue_ = 1;
	 if (type_string == "pi-pi")
	    stacking_type = PI_PI_STACKING;
	 if (type_string == "pi-cation")
	    stacking_type = PI_CATION_STACKING;
	 if (type_string == "cation-pi")
	    stacking_type = CATION_PI_STACKING;
	 if (stacking_type == PI_PI_STACKING)
	    ligand_ring_atom_names = ligand_ring_atom_names_in;
	 if (stacking_type == PI_CATION_STACKING)
	    ligand_ring_atom_names = ligand_ring_atom_names_in;
	 if (stacking_type == CATION_PI_STACKING)
	    ligand_cation_atom_name = ligand_cation_atom_name_in;
      }
      int get_stacking_type() const {
	 return stacking_type;
      }
      std::string get_ligand_cation_atom_name() const {
	 return ligand_cation_atom_name;
      }
      bool se_diff_set() const {
	 return se_diff_set_;
      }
      std::pair<double, double> solvent_exposures() const {
	 return std::pair<double, double> (se_holo, se_apo);
      }
      void set_water_dist_to_protein(double d) {
	 water_dist_to_protein = d;
      } 
      std::vector<std::string> get_ligand_ring_atom_names() const {
	 return ligand_ring_atom_names;
      }
      bool has_ring_stacking_interaction() const {
	 return (ligand_ring_atom_names.size() > 0);
      }
      bool is_a_primary_residue() const {
	 return is_a_primary_residue_;
      }
      std::vector<std::pair<lig_build::pos_t, double> >
      get_attachment_points(const widgeted_molecule_t &mol) const;
      // friend std::ostream& operator<<(std::ostream &s, residue_circle_t r);
   };
   // std::ostream& operator<<(std::ostream &s, residue_circle_t r);


   class grid_index_t {

      int ii_;
      int jj_;
   public: 
      enum { INVALID_INDEX = -1 };
      grid_index_t(int i, int j) {
	 ii_ = i;
	 jj_ = j;
      }
      // hhmmmm!  needed to compile contour_fragment
      // constructor, which takes a const reference
      // to a grid_index_t as an argument. I don't
      // understand why this is needed - or if it
      // works.      
      grid_index_t() {
	 ii_ = INVALID_INDEX;
	 jj_ = INVALID_INDEX;
      }
      bool is_valid_p() const {
	 return (ii_ != INVALID_INDEX);
      }
      int i() const { return ii_;}
      int j() const { return jj_;}
      bool operator==(const grid_index_t &grid_in) const {
	 if (grid_in.i() != ii_) { 
	    return 0;
	 } else { 
	    if (grid_in.j() != jj_) { 
	       return 0;
	    } else { 
	       return 1; // they match
	    }
	 }
      } 
   };

   class ligand_grid {
      double scale_fac;
      lig_build::pos_t top_left;
      lig_build::pos_t bottom_right;
      std::vector<std::vector<double> > grid_;
      int x_size_;
      int y_size_;
      void normalize(); // scale peak value to 1.0
      std::pair<int, int> canvas_pos_to_grid_pos(const lig_build::pos_t &atom_pos) const;
      int square_type(int ii, int jj, float contour_level) const;
      std::vector<std::vector<lig_build::pos_t> > make_contour_lines(const std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > &line_fragments) const;
      double substitution_value(double r_squared, double bash_dist) const;
      // can throw a std::runtime_error (if result is out of grid)
      grid_index_t grid_pos_nearest(const lig_build::pos_t &pos) const;
      
      
   public:
      // (low means low numbers, not low on the canvas)
      // 
      ligand_grid(const lig_build::pos_t &low_x_and_y,
		  const lig_build::pos_t &high_x_and_y);

      void plot_contour_lines(const std::vector<std::vector<lig_build::pos_t> > &contour_lines, GooCanvasItem *root) const;
      enum { MS_NO_CROSSING = -2,
	     MS_NO_SQUARE = -1, 
	     MS_UP_0_0,
	     MS_UP_0_1,
	     MS_UP_1_0,
	     MS_UP_1_1,
	     MS_UP_0_0_and_0_1,
	     MS_UP_0_0_and_1_0,
	     MS_UP_0_0_and_1_1, // hideous valley
	     MS_UP_0_1_and_1_0, // hideous valley
	     MS_UP_0_1_and_1_1,
	     MS_UP_1_0_and_1_1,
	     MS_UP_0_0_and_0_1_and_1_0,
	     MS_UP_0_0_and_0_1_and_1_1,
	     MS_UP_0_0_and_1_0_and_1_1,
	     MS_UP_0_1_and_1_0_and_1_1,
	     };
      
      // lig_build::pos_t to_canvas_pos(const int &ii, const int &jj) const;
      lig_build::pos_t to_canvas_pos(const double &ix, const double &iy) const;
      
      // Actually, not exactly zero but something small.
      // Don't return a grid-point/position that matches anything in
      // already_positioned.
      // 
      std::pair<lbg_info_t::grid_index_t, lig_build::pos_t>
      find_nearest_zero(const lig_build::pos_t &pos,
			const std::vector<lbg_info_t::grid_index_t> &already_positioned) const;

      // arg is not const reference because get_ring_centres() caches
      // the return value inside mol.
      // 
      void fill(widgeted_molecule_t mol);
      double get(int i, int j) const {
	 return grid_[i][j];
      }
      int x_size() const {
	 return x_size_;
      }
      int y_size() const {
	 return y_size_;
      }
      void add_quadratic(const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points);
      lig_build::pos_t find_minimum_position() const;

      void avoid_ring_centres(std::vector<std::vector<std::string> > ring_atoms_list,
			      const widgeted_molecule_t &mol);

      void add_for_accessibility(double bash_dist, const lig_build::pos_t &atom_pos);

      // fudge in a smoothly varing function (so that the contouring
      // behaves smoothly, rather that the jaggies that we'd get if we
      // added a top hat function values to 1.0A.
      // 
      void add_for_accessibility_no_bash_dist_atom(double scale, const lig_build::pos_t &atom_pos);

      void show_contour(GooCanvasItem *root, float contour_level) const;
      // the "cutting" of the contour behaves differently if the
      // unlimited atom is a member of a ring (compared to if it is
      // not).
      void show_contour(GooCanvasItem *root, float contour_level,
			const std::vector<widgeted_atom_ring_centre_info_t> &unlimited_atoms,
			const std::vector<std::vector<std::string> > &ring_atoms_list) const;

   };

   // The lines constituting the fragment and the indices of the next
   // square for the contour line that we are chasing (the indices are
   // not necessarility valid).
   // 
   class contour_fragment {

   public:
      enum { X_AXIS_LOW, X_AXIS_HIGH, Y_AXIS_LOW, Y_AXIS_HIGH };
      class coordinates {
	 float frac_x;
	 float frac_y;
	 int i_ax;
	 bool x_y_axis;
	 
      public:
	 coordinates() { frac_x = 0; frac_y = 0; i_ax = 0; }
	 coordinates(float f, int i) {
	    if (f>1.0)
	       std::cout << "-----> Bad frac " << f << std::endl;
	    if (f<0.0)
	       std::cout << "-----> Bad frac " << f << std::endl;
	    frac_y = 11; // should not be used
	    frac_x = f;
	    i_ax = i;
	    if (i == X_AXIS_LOW)
	       frac_y = 0.0;
	    else
	       frac_y = 1.0;
	    if (i_ax != X_AXIS_LOW)
	       if (i_ax != X_AXIS_HIGH)
		  std::cout << "Bad axis to coordinates(f, i) f: "
			    << f << "  i: " << i << std::endl;
	 } 
	 coordinates(int i, float f) {
	    if (f>1.0)
	       std::cout << "----->  Bad frac " << f << std::endl;
	    if (f<0.0)
	       std::cout << "----->  Bad frac " << f << std::endl;
	    frac_x = 11; // should not be used
	    frac_y = f;
	    i_ax = i;
	    if (i == Y_AXIS_LOW)
	       frac_x = 0.0;
	    else
	       frac_x = 1.0;
	    if (i_ax != Y_AXIS_LOW)
	       if (i_ax != Y_AXIS_HIGH)
		  std::cout << "Bad axis to coordinates(i, f) i: "
			    << i << "  f: " << f << std::endl;
	 }
	 float get_frac_x() { return frac_x; } 
	 float get_frac_y() { return frac_y; } 
      };
      
      grid_index_t grid_index_next;
      lig_build::pos_t start_point; // on either the x or y axis
      lig_build::pos_t end_point;
      contour_fragment(int ms_type,
		       const float &contour_level,
		       const grid_index_t &grid_index_prev,
		       const grid_index_t &grid_index,
		       const ligand_grid &grid);

      typedef std::pair<coordinates, coordinates> cp_t;
      std::vector<cp_t> coords;
      std::pair<double, double> get_coords(int ii, int jj, int coord_indx) {
	 coordinates c;
	 if (coord_indx == 0)
	    c = coords[0].first;
	 if (coord_indx == 1)
	    c = coords[0].second;
	 
	 // these are for hideous value (two crossing vectors)
	 if (coord_indx == 2)
	    c = coords[1].first;
	 if (coord_indx == 3)
	    c = coords[1].second;

	 return std::pair<double, double> (ii+c.get_frac_x(), jj+c.get_frac_y());
      } 

   };




   class optimise_residue_circles {
   private:

      class angle {
      public:
	 angle(int i_atom_index_in, int ires_1_in, int ires_2_in) {
	    i_atom_index = i_atom_index_in;
	    ires_1_index = ires_1_in;
	    ires_2_index = ires_2_in;
	 }
	 int i_atom_index;
	 int ires_1_index;
	 int ires_2_index;
      };
      int status; 
      static double f(const gsl_vector *v, void *params);
      static void  df(const gsl_vector *v, void *params, gsl_vector *df); 
      static void fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
      
      std::vector<residue_circle_t> starting_circles;
      std::vector<residue_circle_t>  current_circles;

      widgeted_molecule_t mol;
      
      void numerical_gradients(gsl_vector *x, gsl_vector *df, void *params) const;
      std::vector<angle> angles;
      void setup_angles(); // uses current_circles and mol and fills angles

      // these can be unstaticed by passing an lbg_info_t to f(), df()
      // as part of params.  Hmmm... it is already.  Should be easy to
      // fix then.  But not now.
      // 
      static bool score_vs_ligand_atoms;
      static bool score_vs_ring_centres;
      static bool score_vs_other_residues;
      static bool score_vs_other_residues_for_angles;
      static bool score_vs_original_positions;
      static bool score_vs_ligand_atom_bond_length;

   public:
      // we pass two vectors here because (for trajectory-view) we
      // don't want to restart the minimisation with the current
      // positions - we want to minimise against the (constant)
      // original positions.
      optimise_residue_circles(const std::vector<residue_circle_t> &r, // starting points
			       const std::vector<residue_circle_t> &c, // current points
			       const widgeted_molecule_t &mol,
			       const std::vector<int> &primary_indices);
      std::pair<int, std::vector<residue_circle_t> > solution() const;
      // return GSL minimisation status;
      int get_gsl_min_status() const { return status; }
      std::vector<int> primary_indices;  // indices in residues_circles
                                         // of the primary residues.
   };

#ifdef MAKE_ENHANCED_LIGAND_TOOLS   
   class alert_info_t {
   public:
      std::string smarts;
      std::string smarts_name;
      RDKit::MatchVectType matches;
      alert_info_t(const std::string &smarts_in,
		   const std::string &smarts_name_in,
		   const RDKit::MatchVectType &matches_in) {
	 smarts = smarts_in;
	 smarts_name = smarts_name_in;
	 matches = matches_in;
      } 
   };
#endif

private:
   bool stand_alone_flag;

   // sometime, of course the ligand_spec doesn't mean anything (not
   // set on construction), hence a flag for that.
   // 
   std::pair<bool, coot::residue_spec_t> ligand_spec_pair;

   bool try_stamp_polygon(int n_edges, int x_pos_centre, int y_pos_centre,
			  bool is_spiro, bool is_aromatic);
   void stamp_polygon_anywhere(int n_edges, int x_pos_centre, int y_pos_centre,
			       bool is_aromatic, GooCanvasItem *root);
   std::vector<int> stamp_polygon(int n_edges, lig_build::polygon_position_info_t ppi,
				  bool aromatic_flag, GooCanvasItem *root);
   std::vector<int> try_stamp_polygon_using_highlighted_data(int n_edges,
							     bool spiro_flag,
							     bool aromomatic_flag,
							     GooCanvasItem *root);
   bool try_stamp_hexagon_aromatic(int x_pos_centre, int y_pos_centre, bool shift_is_pressed);
   bool handle_charge_change();
   std::vector<widgeted_molecule_t> previous_molecules;
   int save_molecule_index;
   bool in_delete_mode_;
   bool draw_flev_annotations_flag; // residue_attributes (circles and
				    // bonds), solvent accessiblity and
				    // substitution contour.
   highlight_data_t highlight_data;
   // we need to store the atom that was the rotation centre when we are
   // rotating a bond. Sometimes we may need to delete a bond as it is rotated
   // We need to be careful that we rotate atoms as we rotate canvas items.
   // most_recent data are relevant if we still have left mouse pressed.
   // 
   highlight_data_t most_recent_drag_atom;
   bool most_recent_bond_made_new_atom_flag;
   bool most_recent_bond_is_being_dragged_flag;
   int atom_index_of_atom_to_which_the_latest_bond_was_added;
   
   bool is_atom_element(int addition_mode) const;
   bool is_bond(int addition_mode) const;

   bool button_down_bond_addition; // set on bond addition by button down, unset by
				   // button up.  Used to distinguish between canvas drag
				   // and mouse-based rotation of the bond.
   
   lig_build::pos_t penultimate_atom_pos; // used with the above as the anchor point of
					  // the rotation of the last added atom.
   int penultimate_atom_index;
   int ultimate_atom_index;
   
   GooCanvasItem *latest_bond_canvas_item; // For (above mentioned) bond rotation, we
					   // need to rotate the bond.  This is the
					   // canvas item of the bond.
   
   bool latest_bond_was_extended; // Again, for the above system, when we are rotating a
				  // bond normally, the rotation_bond() is fine.
				  // However, if, when rotating the bond, the bond
				  // snaped/extended to a highlighted atom, then simply
				  // rotating the atom position and the widget won't do
				  // (when the bond is snapped/extended the new atom (at
				  // the end of the rotating stick) is removed and the
				  // bond is adjusted to use the atom index of the
				  // highlighted atom.  In that case, we don't want to
				  // simply rotate the bond, that would be bad. OK,
   
   bool try_change_to_element(int addition_element_mode); // check for highlighted atom;
   bool try_add_or_modify_bond(int canvas_addition_mode, int x, int y,
			       bool button_1_is_pressed); //  ditto.
   // return "was changed" and new-atom-was-created status pair
   std::pair<bool, bool> add_bond_to_atom(unsigned int atom_index, int canvas_addition_mode);
   // return new-atom-was-created status
   bool add_bond_to_atom_with_0_neighbours(unsigned int atom_index, int canvas_addition_mode);
   bool add_bond_to_atom_with_1_neighbour(unsigned int atom_index, int canvas_addition_mode,
					  unsigned int bond_index);
   bool add_bond_to_atom_with_2_neighbours(unsigned int atom_index, int canvas_addition_mode,
					   const std::vector<unsigned int> &bond_indices);
   bool add_bond_to_atom_with_3_neighbours(unsigned int atom_index, int canvas_addition_mode,
					   const std::vector<unsigned int> &bond_indices);
   std::string to_element(int addition_mode) const;
   std::string font_colour(int addition_element_mode) const;
   std::string font_colour(const std::string &ele) const;
   lig_build::bond_t::bond_type_t addition_mode_to_bond_type(int canvas_addition_mode) const;
   void try_stamp_bond_anywhere(int canvas_addition_mode, int x_mouse, int y_mouse); // always modifies.
   bool change_atom_element(unsigned int atom_index, std::string new_element, std::string fc);
   void change_atom_id_maybe(unsigned int atom_index);
   lig_build::pos_t mouse_at_click;
   std::string mdl_file_name; // for save function.
   void add_search_combobox_text() const;
   bool make_saves_mutex;
   double canvas_scale;
   bool use_graphics_interface_flag;
   void init_internal() {
      in_delete_mode_ = 0;
      key_group = NULL;      
      save_molecule_index = UNASSIGNED_INDEX;
      make_saves_mutex = 1; // allow saves initially
      search_similarity = 0.95;
      coot_mdl_ready_time = 0;
      canvas_scale = 1.0;
      canvas_drag_offset =  lig_build::pos_t(0,0);
      top_left_correction = lig_build::pos_t(0,0);
      most_recent_bond_made_new_atom_flag = false;
      standard_residue_circle_radius = 19;
      button_down_bond_addition = false;
      latest_bond_canvas_item = 0;
      penultimate_atom_index = UNASSIGNED_INDEX;
      ultimate_atom_index = UNASSIGNED_INDEX;
      latest_bond_was_extended = 0;
      atom_index_of_atom_to_which_the_latest_bond_was_added = UNASSIGNED_INDEX;
      stand_alone_flag = 0;
      ligand_spec_pair.first = 0; // unset ligand_spec
      use_graphics_interface_flag = 1; // default: show gui windows and widgets.
      mdl_file_name = "coot.mol";
      atom_X = "H";
      comp_id = "LIG";
      lbg_atom_x_dialog = NULL;
      lbg_atom_x_entry = NULL;
      get_url_func_ptr_flag = false;
      prodrg_import_func_ptr = NULL;
      sbase_import_func_ptr = NULL;
      get_drug_mdl_file_function_pointer = NULL;

      orient_view_func                               = NULL;
      set_rotation_centre_func                       = NULL;
      set_show_additional_representation_func        = NULL;
      all_additional_representations_off_except_func = NULL;
      
      draw_flev_annotations_flag = false;
//       lbg_nitrogen_toggle_toolbutton = NULL;
//       lbg_carbon_toggle_toolbutton   = NULL;
//       lbg_oxygen_toggle_toolbutton   = NULL;
//       lbg_sulfur_toggle_toolbutton   = NULL;
//       lbg_phos_toggle_toolbutton     = NULL;
//       lbg_fluorine_toggle_toolbutton = NULL;
//       lbg_chlorine_toggle_toolbutton = NULL;
//       lbg_bromine_toggle_toolbutton  = NULL;
//       lbg_single_bond_toggle_toolbutton  = NULL;
//       lbg_double_bond_toggle_toolbutton  = NULL;
//       lbg_ring_8_toggle_toolbutton = NULL;
//       lbg_ring_7_toggle_toolbutton = NULL;
//       lbg_ring_6_toggle_toolbutton = NULL;
//       lbg_ring_6_arom_toggle_toolbutton = NULL;
//       lbg_ring_5_toggle_toolbutton = NULL;
//       lbg_ring_4_toggle_toolbutton = NULL;
//       lbg_ring_3_toggle_toolbutton = NULL;
      lbg_show_alerts_checkbutton    = NULL;
      lbg_alert_hbox_outer = NULL;
      alert_group = NULL; // group for alert annotations
      show_alerts_user_control = false; // no pattern matching available
      geom_p = NULL; // no (const) geometry passed/set
      display_atom_names   = false;
      display_atom_numbers = false;
#ifdef MAKE_ENHANCED_LIGAND_TOOLS   
      show_alerts_user_control = false;
      bond_pick_pending = false;
      atom_pick_pending = false;
#ifdef USE_PYTHON
      user_defined_alerts_smarts_py = NULL;
      setup_silicos_it_qed_default_func();
      setup_user_defined_alert_smarts();
#endif      
#endif      
   }
   
   // return a status and a vector of atoms (bonded to atom_index) having
   // only one bond.
   // 
   std::pair<bool, std::vector<unsigned int> > 
   have_2_stubs_attached_to_atom(unsigned int atom_index,
				 const std::vector<unsigned int> &bond_indices) const;
   void squeeze_in_a_4th_bond(unsigned int atom_index, int canvas_addition_mode,
			      const std::vector<unsigned int> &bond_indices);
   std::vector<double>
   get_angles(unsigned int atom_index, const std::vector<unsigned int> &bond_indices) const;
   lig_build::pos_t  new_pos_by_bisection(unsigned int atom_index,
					  const std::vector<unsigned int> &bond_indices,
					  const std::vector<double> &angles,
					  GooCanvasItem *root) const;
   bool all_closed_rings(unsigned int atom_index, const std::vector<unsigned int> &bond_indices) const;
   std::vector<lig_build::pos_t>
   get_centres_from_bond_indices(const std::vector<unsigned int> &bond_indices) const;
   lig_build::pos_t get_new_pos_not_towards_ring_centres(unsigned int atom_index,
							 const std::vector<unsigned int> &bond_indices) const;
   mmdb::Manager *get_cmmdbmanager(const std::string &filename) const;

   // sbase functions
   float search_similarity;
   
#ifdef HAVE_CCP4SRS   
   coot::match_results_t residue_from_best_match(mmdb::math::Graph &graph_1, mmdb::math::Graph &graph_2,
						 mmdb::math::GraphMatch &match, int n_match, 
						 ccp4srs::Monomer *monomer_p) const;
#endif   
   // not const because try_dynamic_add() can be called (to make images):
   void display_search_results(const std::vector<coot::match_results_t> &v);
   void rotate_latest_bond(int x, int y);
   void rotate_or_extend_latest_bond(int x, int y);
   bool extend_latest_bond_maybe(); // use hilight_data

   // return the bond between the ligand "core" and the atom_index
   // (one of bond_indices)
   // 
   widgeted_bond_t orthogonalise_2_bonds(unsigned int atom_index,
					 const std::vector<unsigned int> &attached_bonds,
					 const std::vector<unsigned int> &bond_indices);
   std::vector<residue_circle_t> residue_circles;
   std::pair<bool,lig_build::pos_t> get_residue_circles_top_left() const;
   lig_build::pos_t top_left_correction; // 0,0 by default
   
   // a set of handles (returned from
   // additional_representation_by_attributes()) that correspond to
   // the residues in residue_circles.  If there are no additional
   // representations, then an empty vector is handled.
   std::vector<int> additional_representation_handles;

   widgeted_molecule_t translate_molecule(const lig_build::pos_t &delta);
   void translate_residue_circles(const lig_build::pos_t &delta);
   void clear_canvas();
   
   // if you don't have an add_rep_handle, then pass -1 (something negative)
   // 
   void draw_residue_circle_top_layer(const residue_circle_t &rc, const lig_build::pos_t &p,
				      int add_rep_handle);
   double standard_residue_circle_radius;
   

   std::pair<std::string, std::string>
   get_residue_circle_colour(const std::string &residue_type) const;

   
   std::vector<int> get_primary_indices() const;
   // twiddle residue_circles, taking in to account the residues that
   // bond to the ligand (including stacking residues).
   // 
   void initial_residues_circles_layout();

   // fiddle with the position of the residue_circles[primary_index].
   //
   void initial_primary_residue_circles_layout(const lbg_info_t::ligand_grid &grid,
					       int primary_index,
					       const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points);

   // untrap residues as needed.
   void position_non_primaries(const lbg_info_t::ligand_grid &grid,
			       const std::vector<int> &primary_indices);

   
   // return 1 on solution having problems, 0 for no problems, also
   // return a list of the residue circles with problems.
   // 
   std::pair<bool, std::vector<int> > solution_has_problems_p() const;
      


   // minimise layout energy
   std::pair<int, std::vector<residue_circle_t> >
   optimise_residue_circle_positions(const std::vector<residue_circle_t> &start,
				     const std::vector<residue_circle_t> &current,
				     const std::vector<int> &primary_indices) const;
   lig_build::pos_t canvas_drag_offset;
   bool is_sane_drag(const lig_build::pos_t &delta) const;

   void draw_solvent_accessibility_of_atom(const lig_build::pos_t &pos, double sa,
					   GooCanvasItem *root);
   void draw_solvent_accessibility_of_atoms();
   
   void draw_substitution_contour();
   void draw_bonds_to_ligand();
   void draw_solvent_exposure_circle(const residue_circle_t &residue_circle,
				     const lig_build::pos_t &ligand_centre,
				     GooCanvasItem *group);
   std::string get_residue_solvent_exposure_fill_colour(double radius_extra) const;
   
   // click_pos is where we recentre in 3D graphics when the annotation
   // (line) is clicked.
   void draw_annotated_stacking_line(const lig_build::pos_t &ligand_ring_centre,
				     const lig_build::pos_t &residue_pos,
				     int stacking_type,
				     const clipper::Coord_orth &click_pos);
   void handle_read_draw_coords_mol_and_solv_acc(const std::string &coot_pdb_file,
						 const std::string &coot_mdl_file,
						 const std::string &sa_file);
   // for debugging the minimization vs original positions
   std::vector<residue_circle_t> offset_residues_from_orig_positions(); 
   void show_grid(const lbg_info_t::ligand_grid &grid);
   void show_mol_ring_centres(); // not const because mol.get_ring_centres() caches
   void show_unlimited_atoms(const std::vector<widgeted_atom_ring_centre_info_t> &ua);
   void show_ring_centres(std::vector<std::vector<std::string> > ring_atoms_list,
			  const widgeted_molecule_t &mol);
   // this can cache ring centres in mol if they are not there already
   void show_ring_centres(widgeted_molecule_t &mol);
   

   std::string grid_intensity_to_colour(int val) const;
   std::string sixteen_to_hex_let(int v) const;
   void reposition_problematics_and_reoptimise(const std::vector<int> &problematics,
					       const std::vector<int> &primary_indices);
   void recentre_considering_residue_centres();  // move atoms and residues
   std::vector<residue_circle_t>
   filter_residue_waters(const std::vector<residue_circle_t> &r_in,
			 double max_dist_water_to_ligand_atom,
			 double max_dist_water_to_protein_atom) const;

   void refine_residue_circle_positions(); // changes the positions of in residue_circles

   std::vector<solvent_accessible_atom_t>
   convert(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
	   const coot::flev_attached_hydrogens_t &ah) const;
   
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   RDKit::Bond::BondType convert_bond_type(const lig_build::bond_t::bond_type_t &t) const;
   RDKit::Bond::BondDir  convert_bond_dir(const lig_build::bond_t::bond_type_t &t) const;
   
   std::string get_smiles_string_from_mol_rdkit() const;
   std::vector<alert_info_t> alerts(const RDKit::ROMol &mol) const;
   void rdkit_mol_post_read_handling(RDKit::RWMol *m, const std::string &file_name, unsigned int iconf=0);
#ifdef USE_PYTHON   
   PyObject *silicos_it_qed_default_func;
   PyObject *silicos_it_qed_properties_func;
   PyObject *silicos_it_qed_pads;
   PyObject * get_callable_python_func(const std::string &module_name,
				       const std::string &function_name) const;
   PyObject *user_defined_alerts_smarts_py;
   void setup_user_defined_alert_smarts();
   void setup_silicos_it_qed_default_func(); // try to get the python function, or set it to null.
#endif    
#endif
   std::string get_smiles_string_from_mol_openbabel() const;

   // this is the value generated in annotate(), containing aromatic
   // and non-aromatic rings.
   // 
   std::vector<std::vector<std::string> > ring_atoms_list;

   // allow access of the Search button callback to the search
   // similarity combox box text.
   double get_search_similarity() const;
   std::string atom_X; // initially "H"

   std::vector<std::pair<std::string, std::string> > alert_smarts() const;
   std::vector<std::pair<std::string, std::string> > user_defined_alert_smarts() const;
   GooCanvasItem *alert_group;

   // get bottom of flev items, used in positioning the key
   double bottom_of_flev_items();

   // geometry
   // (not const so that we can call try_dynamic_add())
   coot::protein_geometry *geom_p;

public:
   lbg_info_t(GtkWidget *canvas_in, coot::protein_geometry *geom_p_in) {
      canvas = canvas_in;
      init_internal();
      geom_p = geom_p_in;
   }
   lbg_info_t() { init_internal(); }
   lbg_info_t(int imol_in, coot::protein_geometry *geom_p_in) {
      init_internal();
      geom_p = geom_p_in;
      imol = imol_in;
      if (0)
	 std::cout << "in lbg_info_t(imol) mdl_file_name is now :"
		   << mdl_file_name << ":" << std::endl;
   }
      
   // toggle button modes, mutually exclusive
   enum { NONE, TRIANGLE, SQUARE, PENTAGON, HEXAGON, HEXAGON_AROMATIC, HEPTAGON, OCTAGON,
	  ATOM_C, ATOM_N, ATOM_O, ATOM_H, ATOM_S, ATOM_P, ATOM_F, ATOM_CL, ATOM_I, ATOM_BR, ATOM_X,
	  CHARGE, ADD_SINGLE_BOND, ADD_DOUBLE_BOND, ADD_TRIPLE_BOND, ADD_STEREO_OUT_BOND,
	  DELETE_MODE};
#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
   bool init(GtkBuilder *builder); // return success status (true is good).
   void setup_lbg_drag_and_drop(GtkWidget *lbg_window);
#endif // GTK_VERSION
   int imol; // the coot molecule number from which this plot was
	     // generated (quite possibly -1, i.e. no coot molecule)
   GtkWidget *lbg_window;
   GtkWidget *about_dialog; 
   GtkWidget *lbg_apply_button;
   GtkWidget *lbg_search_combobox;
   GtkWidget *open_dialog;
   GtkWidget *save_as_dialog;
   GtkWidget *lbg_export_as_pdf_dialog;
   GtkWidget *lbg_export_as_png_dialog;
   GtkWidget *lbg_export_as_svg_dialog;
   GtkWidget *lbg_sbase_search_results_dialog;
   GtkWidget *lbg_sbase_search_results_vbox;
   GtkWidget *lbg_smiles_dialog;
   GtkWidget *lbg_smiles_entry;
   GtkWidget *lbg_import_from_smiles_dialog;
   GtkWidget *lbg_import_from_smiles_entry;
   GtkWidget *lbg_import_from_comp_id_dialog;
   GtkWidget *lbg_import_from_comp_id_entry;
   GtkWidget *lbg_import_from_comp_id_hydrogens_checkbutton;
   GtkWidget *lbg_statusbar;
   GtkWidget *lbg_toolbar_layout_info_label;
   GtkWidget *lbg_atom_x_dialog;
   GtkWidget *lbg_atom_x_entry;
   GtkWidget *lbg_qed_hbox;
   GtkWidget *lbg_qed_text_label;
   GtkWidget *lbg_qed_progressbar;
   GtkWidget *lbg_alert_hbox; // controlled by alerts in structure
   GtkWidget *lbg_alert_hbox_outer; // controled by user wanting to see alerts in structure
   GtkWidget *lbg_alert_name_label;
   GtkWidget *lbg_show_alerts_checkbutton; 
   GtkWidget *lbg_get_drug_dialog;
   GtkWidget *lbg_get_drug_entry;
   GtkWidget *lbg_get_drug_menuitem;
   GtkWidget *lbg_flip_rotate_hbox;
   GtkWidget *lbg_clean_up_2d_toolbutton;
   GtkWidget *lbg_search_database_frame;
   GtkWidget *lbg_scale_spinbutton;
   GtkWidget *lbg_view_rotate_entry;
   GtkWidget *lbg_qed_properties_vbox; // hide if not enhanced-ligand
   GtkWidget *lbg_qed_properties_progressbars[8];
   GtkWidget *lbg_srs_search_results_scrolledwindow;
   GtkWidget *lbg_srs_search_results_vbox;
//    GtkWidget *lbg_nitrogen_toggle_toolbutton;
//    GtkWidget *lbg_carbon_toggle_toolbutton;
//    GtkWidget *lbg_oxygen_toggle_toolbutton;
//    GtkWidget *lbg_sulfur_toggle_toolbutton;
//    GtkWidget *lbg_phos_toggle_toolbutton;
//    GtkWidget *lbg_fluorine_toggle_toolbutton;
//    GtkWidget *lbg_chlorine_toggle_toolbutton;
//    GtkWidget *lbg_bromine_toggle_toolbutton;
//    GtkWidget *lbg_single_bond_toggle_toolbutton;
//    GtkWidget *lbg_double_bond_toggle_toolbutton;
//    GtkWidget *lbg_ring_8_toggle_toolbutton;
//    GtkWidget *lbg_ring_7_toggle_toolbutton;
//    GtkWidget *lbg_ring_6_toggle_toolbutton;
//    GtkWidget *lbg_ring_6_arom_toggle_toolbutton;
//    GtkWidget *lbg_ring_5_toggle_toolbutton;
//    GtkWidget *lbg_ring_4_toggle_toolbutton;
//    GtkWidget *lbg_ring_3_toggle_toolbutton;
   GtkWidget *canvas;
   GooCanvasItem *key_group;
   std::map<std::string, GtkToggleToolButton *> widget_names;

   // when we have multiple molecule, these things should go together
   widgeted_molecule_t mol;
   bool display_atom_names;
   bool display_atom_numbers; // e.g. N:11, C:12
   
   int canvas_addition_mode;
   void lbg_toggle_button_my_toggle(GtkToggleToolButton *tb);
   void save_molecule(); // moved so that function lbg() can save on start.
#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
   bool save_togglebutton_widgets(GtkBuilder *builder);
#endif // GTK_VERSION
   void handle_item_add(GdkEventButton *event);
   void handle_item_delete(GdkEventButton *event);
   void untoggle_others_except(GtkToggleToolButton *button_toggled_on);
   bool item_highlight_maybe(int x, int y);
   void highlight_bond(const lig_build::bond_t &bond, bool delete_mode);
   void highlight_atom(const lig_build::atom_t &atom, int atom_index, bool delete_mode);
   void remove_bond_and_atom_highlighting();
   void clear_button_down_bond_addition();
   bool button_down_bond_addition_state() { return button_down_bond_addition; } 
   void handle_drag(GdkModifierType state, int x, int y); // could be rotate bond or drag canvas
   void set_in_delete_mode(bool v) {
      in_delete_mode_ = v;
   }
   void set_draw_flev_annotations(bool v) {
      draw_flev_annotations_flag = v;
   } 
   bool in_delete_mode_p() const { return in_delete_mode_; }
   double radius(int n_edges) const; // depends on zoom? (for future).
   void clear(bool do_descriptor_updates);
   std::string get_stroke_colour(int i, int n) const;
   void drag_canvas(int mouse_x, int mouse_y);
   void write_pdf(const std::string &file_name) const;
   void write_ps(const std::string &file_name) const;
   void write_png(const std::string &file_name);
   void write_svg(const std::string &file_name) const;
   void set_mouse_pos_at_click(int xpos, int ypos) {
      mouse_at_click = lig_build::pos_t(double(xpos), double(ypos));
   }
   void render(); // uses internal data member mol
   void update_descriptor_attributes(); // this is not in render_from_molecule() because it can/might be slow.
   void update_apply_button_sensitivity(); // turn off the "Apply" button is the molecule is not sane.
   void delete_hydrogens();
   void undo();
#ifdef HAVE_CCP4SRS
   // not const because try_dynamic_add() can be called.
   void search();
#endif

   // update the internal class variable widgeted_molecule_t mol from mol_in
   void import_from_widgeted_molecule(const widgeted_molecule_t &mol_in);
   void import_molecule_from_file(const std::string &file_name); // mol or cif
   void import_molecule_from_cif_file(const std::string &file_name); // cif
   // 20111021 try to read file_name as a MDL mol or a mol2 file.
   void import_mol_from_file(const std::string &file_name);
   // read an MDL mol file.
   widgeted_molecule_t  import_mol_file(const lig_build::molfile_molecule_t &mol_in,
					const std::string &filename,
					mmdb::Manager *pdb_mol);
   void import_via_rdkit_from_restraints_dictionary(const coot::dictionary_residue_restraints_t &dict, bool show_hydrogens_status);

   void import_mol_from_smiles_file(const std::string &file_name);
   void import_mol_from_smiles_string(const std::string &smiles);
   void import_mol_from_comp_id(const std::string &comp_id,
				bool show_hydrogens_flag);
   
   static gboolean on_highlight_key_press_event (GooCanvasItem *item,
						 GooCanvasItem *target,
						 GdkEventKey *event,
						 gpointer data);
   void handle_key_press_button_toggle(int key, bool ctrl_is_pressed);
   
   // and the version of that not going via an intermediate molfile_molecule_t
#ifdef MAKE_ENHANCED_LIGAND_TOOLS   
   widgeted_molecule_t import_rdkit_mol(RDKit::ROMol *rdkm, int iconf) const;
#endif // MAKE_ENHANCED_LIGAND_TOOLS   

   static void on_sbase_search_result_button_clicked(GtkButton *button, gpointer user_data);
   static gboolean watch_for_mdl_from_coot(gpointer user_data);
   time_t coot_mdl_ready_time;
   std::pair<lig_build::pos_t, lig_build::pos_t> flev_residues_extents() const; // for canvas sizing
   void read_draw_residues(const std::string &file_name);
   void draw_all_flev_annotations(); // which calls draw_all_residue_attribs();
   void draw_all_flev_residue_attribs();
   void draw_all_flev_ligand_annotations();
   
   
   std::vector<residue_circle_t> read_residues(const std::string &file_name) const;
   // if you don't have add_rep_handles, pass an empty vector.
   void draw_residue_circles(const std::vector<residue_circle_t> &v,
			     const std::vector<int> &add_rep_handles);
   void draw_stacking_interactions(const std::vector<residue_circle_t> &v);
   std::vector<solvent_accessible_atom_t> 
   read_solvent_accessibilities(const std::string &file_name) const;
   void read_files_from_coot();
   std::string file_to_string(const std::string &filename) const;
   std::string get_flev_analysis_files_dir() const;
   void show_key();
   void hide_key();

   // non-enterprise path
   void update_statusbar_smiles_string() const;
   void update_statusbar_smiles_string(const std::string &smiles_string) const;
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   RDKit::RWMol rdkit_mol(const widgeted_molecule_t &mol) const;
   // do these need to be RWMols?
   void update_statusbar_smiles_string(const RDKit::ROMol &rdkm) const;
   std::string get_smiles_string(const RDKit::ROMol &mol) const;

   void update_qed(const RDKit::RWMol &rdkm);
   void update_qed_properties(const std::vector<std::pair<double, double> > &d);
   void reset_qed_properties_progress_bars(); // on exception on molecule editing and clear()
   
   void update_alerts(const RDKit::RWMol &rdkm);
   std::string get_smiles_string_from_mol(const RDKit::RWMol &mol) const;
   bool bond_pick_pending;
   bool atom_pick_pending;
   // and bond pick is stored here:
   int pending_action_data_picked_bond_index;
   bool handle_bond_picking_maybe();
#endif
   // although these depend on/manipulate rdkit-based entities - they
   // only clear them up, not generate them, so, in order that we
   // don't complicated lbg-callbacks.cc with MAKE_ENHANCED_LIGAND_TOOLS
   // dependencies, let's put those functions outside the
   // MAKE_ENHANCED_LIGAND_TOOLS dependency here.
   void clear_canvas_alerts();
   bool show_alerts_user_control; 
   
   // can throw an exception
   std::string get_smiles_string_from_mol() const;
   void set_stand_alone() { stand_alone_flag = 1; }
   bool is_stand_alone() { return stand_alone_flag; }

   bool annotate(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
		 const std::vector<coot::fle_residues_helper_t> &centres,
		 const std::vector<int> &addition_representation_handles,
		 const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
		 const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
		 const coot::flev_attached_hydrogens_t &ah,
		 const coot::pi_stacking_container_t &pi_stack_info,
		 const coot::dictionary_residue_restraints_t &restraints);

   void set_ligand_spec(const coot::residue_spec_t &spec) {
      ligand_spec_pair.first = 1;
      ligand_spec_pair.second = spec;
   }

   std::pair<bool, coot::residue_spec_t> get_ligand_spec() const {
      return ligand_spec_pair;
   }

   void no_graphics_mode() {
      use_graphics_interface_flag = false;
   }

   void show_atom_x_dialog();
   void set_atom_x_string(const std::string &s);

   void write_mdl_molfile_using_default_file_name() const;

   void set_default_mdl_file_name(const std::string &file_name) {
      mdl_file_name = file_name;
   }

   void clear_and_redraw();
   void clear_and_redraw(const lig_build::pos_t &delta);

   // drag and drop callbacks
   int handle_lbg_drag_and_drop_string(const std::string &uri);
   int handle_lbg_drag_and_drop_single_item(const std::string &uri);
   int handle_lbg_drag_and_drop_chemspider_image(const std::string &uri);
   int handle_lbg_drag_and_drop_chemspider_structure(const std::string &uri);
   int handle_lbg_drag_and_drop_pubchem_image(const std::string &uri);
   int handle_lbg_drag_and_drop_filesystem_file(const std::string &uri);
   int handle_lbg_drag_and_drop_drugbank(const std::string &uri,
					 const std::string &url_file_name_file);
   int handle_lbg_drag_and_drop_mol_file(const std::string &uri_clean,
					 const std::string &url_file_name_file);
   int handle_lbg_drag_and_drop_smiles(const std::string &smiles);

   std::string get_id_string(const std::string &s, int prefix_len, int max_len) const;
   int get_chemspider_mol(const std::string &id_string);
   int get_pubchem_sid_mol(const std::string &id_string);
   int get_pubchem_cid_mol(const std::string &id_string);
   int get_pubchem_mol_generic(const std::string &id_string, const std::string &id_type);


   // we want to use curl functions, but they are declared - and
   // stored in the src functions.  We can't move them down the
   // hierarchy because curl function handlers are stored in
   // graphics_info_t.  So pass lbg a pointer to a function
   // (coot_get_url).
   //
   int (*get_url_func_ptr) (const char *s1, const char *s2);
   bool get_url_func_ptr_flag; 
   void set_curl_function(int (*f) (const char *s1, const char *s2)) {
      get_url_func_ptr = f;
      get_url_func_ptr_flag = true;
   }

   // handle PRODRG output
   // 
   void (*prodrg_import_func_ptr) (std::string file_name_in, std::string comp_id);
   void set_prodrg_import_function(void (*f) (std::string, std::string)) {
      prodrg_import_func_ptr = f;
   }
   // handle SBase input, i.e. when a user clicks on a sbase-search
   // button result, then call an action which affects the graphics.
   // the string is the comp_id to import.
   void (*sbase_import_func_ptr) (std::string file_name);
   void set_sbase_import_function(void (*f) (std::string)) {
      sbase_import_func_ptr = f;
   }
   // Let's have a wrapper around that so that lbg-search doesn't need to ask if sbase_import_func_ptr
   // is valid or not.
   void import_srs_monomer(const std::string &comp_id);

   void import_prodrg_output(const std::string &prodrg_mdl_file_name, const std::string &comp_id) {
      if (prodrg_import_func_ptr) {
	 prodrg_import_func_ptr(prodrg_mdl_file_name, comp_id);
      } else {

	 // all we do is write the file.
	 // update the status bar.
	 //
	 std::string status_string = "  Wrote file " + prodrg_mdl_file_name;
	 guint statusbar_context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(lbg_statusbar),
								   status_string.c_str());
	 gtk_statusbar_push(GTK_STATUSBAR(lbg_statusbar),
			    statusbar_context_id,
			    status_string.c_str());
      
	 
      }
   }

   // instead of hardwiring "DRG" into on_lbg_apply_button_clicked(), allow the user to set the
   // three-letter-code (we call that variable comp_id)
   std::string comp_id; // make private?
   std::string get_comp_id() const { return comp_id; }

   // handle the net transfer of drug (to mdl file)
   //
   std::string (*get_drug_mdl_file_function_pointer) (std::string drug_name);
   void set_get_drug_mdl_file_function(std::string (*get_drug_mdl_file_function_pointer_in) (std::string drug_name)) {
      get_drug_mdl_file_function_pointer = get_drug_mdl_file_function_pointer_in;
   }

   void get_drug_using_entry_text(); // uses lbg_get_drug_entry
   void get_drug(const std::string &drug_name); // get mol file and load it

   void new_lbg_window();
   void clean_up_2d_representation(); // using rdkit

#ifdef HAVE_CCP4SRS
   // not const because we can call geom_p->try_dynamic_add()
   GtkWidget *get_image_widget_for_comp_id(const std::string &comp_id, int imol, ccp4srs::Manager *srs_manager);
#endif // HAVE_CCP4SRS

   void pe_test_function();

   // ----------------------------------------------------- functions from src ----------------------
private:

   void (*orient_view_func) (int imol,
			     const coot::residue_spec_t &central_residue_spec,
			     const coot::residue_spec_t &neighbour_residue_spec);
   void (*set_rotation_centre_func) (const clipper::Coord_orth &pos);
   void (*set_show_additional_representation_func) (int imol, int representation_number, int on_off_flag);
   void (*all_additional_representations_off_except_func) (int imol, int representation_number,
							   short int ball_and_sticks_off_too_flag);
   
public:
   void set_orient_view_func(void (*f)(int imol,
				       const coot::residue_spec_t &central_residue_spec,
				       const coot::residue_spec_t &neighbour_residue_spec)) {
      orient_view_func = f;
   }
   void set_set_rotation_centre_func(void (*f) (const clipper::Coord_orth &pos)) {
      set_rotation_centre_func = f;
   }
   void set_set_show_additional_representation_func(void (*f) (int imol, int representation_number, int on_off_flag)) {
      set_show_additional_representation_func = f;
   }
   void set_all_additional_representations_off_except_func(void (*f) (int imol,
								      int representation_number,
								      short int ball_and_sticks_off_too_flag)) {
      all_additional_representations_off_except_func = f;
   }

   // flipping
   void flip_molecule(int axis);
   void rotate_z_molecule(double angle); // in degrees
   void rotate_z_molecule(const std::string &angle); // in degrees (used in on_lbg_view_rotate_apply_button_clicked
                                                     // callback).
   void scale_canvas(double sf);

   // -- actually run the functions if they were set:
   void orient_view(int imol_in,
		    const coot::residue_spec_t &central_residue_spec,
		    const coot::residue_spec_t &neighbour_residue_spec) {
      if (orient_view_func) {
	 (*orient_view_func)(imol_in, central_residue_spec, neighbour_residue_spec);
      } 
   } 
   void set_rotation_centre(const clipper::Coord_orth &pos) {
      if (set_rotation_centre_func) {
	 (*set_rotation_centre_func)(pos);
      } 
   } 
   void set_show_additional_representation(int imol_in, int representation_number, int on_off_flag) {
      if (set_show_additional_representation_func) {
	 (*set_show_additional_representation_func)(imol_in, representation_number, on_off_flag);
      }
   } 
   void all_additional_representations_off_except(int imol_in, int representation_number,
						  short int ball_and_sticks_off_too_flag) {
      if (all_additional_representations_off_except_func) {
	 (*all_additional_representations_off_except_func) (imol_in, representation_number, ball_and_sticks_off_too_flag);
      }
   }

   void set_display_atom_names(bool state) {
      display_atom_names = state;
   }
   void set_display_atom_numbers(bool state) {
      display_atom_numbers = state;
   }
   
   
};

// return pointer to an lbg_info_t.  Caller deletes.
//
// If there is no molecule (typically, mol is NULL) then pass -1 for
// the imol.
// 
lbg_info_t *lbg(lig_build::molfile_molecule_t mm,
		std::pair<bool, coot::residue_spec_t> ligand_spec_pair,
		mmdb::Manager *mol,
		const std::string &view_name, // annotate the decoration
		const std::string &molecule_file_name,
		int imol, // molecule number of the molecule of the
		coot::protein_geometry *geom_p_in,
			  // layed-out residue
		bool use_graphics_interface_flag,
		bool stand_alone_flag_in,
		int (*get_url_func_pointer) (const char *s1, const char *s2),
		void (*prodrg_import_function_pointer) (std::string file_name_in, std::string comp_id),
		void (*sbase_import_function_pointer) (std::string comp_id),
		std::string (*get_drug_mdl_file_function_pointer) (std::string drug_name)
		);

#endif // LBG_HH

