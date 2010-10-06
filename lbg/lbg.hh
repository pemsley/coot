/* lbg/lbg.hh
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
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

#ifndef LBG_HH
#define LBG_HH

#include <iostream>
#include <map>
#include <queue>


#ifdef MAKE_ENTERPRISE_TOOLS
#include "rdkit-interface.hh"
#endif 

#include <gtk/gtk.h>
#include <goocanvas.h>

#include "gsl/gsl_multimin.h"

#include "mmdb_manager.h"
#include "mmdb_sbase.h"
#define MONOMER_DIR_STR "COOT_SBASE_DIR"

#include "lig-build.hh"
#include "lbg-molfile.hh"

#include "some-coot-utils.hh"

#include "wmolecule.hh"

#define dark "#111111"

// Sorry Bernie, I don't see how this is useful - tell me why.  It
// causes a compilation issue for me (no compat defined in includes).
// That can be addressed another way of course, but I don't want to do
// that if this is not useful.
// 
// #include "coot-sysdep.h"

static double LIGAND_TO_CANVAS_SCALE_FACTOR = 23;
static double SINGLE_BOND_CANVAS_LENGTH= LIGAND_TO_CANVAS_SCALE_FACTOR * 1.54;


bool lbg(lig_build::molfile_molecule_t mm, CMMDBManager *mol, const std::string &molecule_file_name,
	 bool stand_alone_flag_in = 0);

bool save_togglebutton_widgets(GtkBuilder *builder);

void lbg_handle_toggle_button(GtkToggleToolButton *tb, GtkWidget *canvas, int mode);
GtkWidget *get_canvas_from_scrolled_win(GtkWidget *scrolled_window);



// ====================================================================
//                     lbg_info_t
// ====================================================================

class lbg_info_t {

public:
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
      }
      highlight_data_t(GooCanvasItem *w_in,
		       const lig_build::pos_t &p, int index_in) {
	 n_atoms_ = 1;
	 pos_1_ = p;
	 highlight_widget = w_in;
	 atom_index = index_in;
	 bond_indices = std::pair<int, int> (UNASSIGNED_INDEX, UNASSIGNED_INDEX);
      }

      highlight_data_t() {
	 highlight_widget = NULL;
	 n_atoms_ = 0;
	 atom_index = UNASSIGNED_INDEX; // unset
	 bond_indices = std::pair<int, int> (UNASSIGNED_INDEX, UNASSIGNED_INDEX);
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
      lig_build::polygon_position_info_t
      get_new_polygon_centre(int n_edges,
			     bool spiro_flag,
			     const double &radius_standard,
			     const double &radius_corrected,
			     const widgeted_molecule_t &mol) const;
      std::pair<int, int> get_bond_indices() const { return bond_indices; }
   }; // finish highlight_data_t class

   // a container for the results of the comparison vs SBase graph matching.
   //
   class match_results_t {
   public:
      bool success;
      std::string name;
      std::string comp_id;
      CResidue *res;
      // clipper::RTop_orth
      match_results_t(const std::string &comp_id_in, const std::string &name_in, CResidue *res_in) {
	 name = name_in;
	 comp_id = comp_id_in;
	 res = res_in;
	 if (res_in)
	    success = 1;
	 else
	    success = 0;
      }
   };

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
      double pos_x; // input coordinate reference frame
      double pos_y;
      double pos_z;
      lig_build::pos_t pos; // coordinate system of the ligand atoms
      std::string residue_type;
      std::string residue_label;
      std::vector<bond_to_ligand_t> bonds_to_ligand;
      double water_dist_to_protein; 
      residue_circle_t(const double &x_in, const double &y_in, const double &z_in,
		       const std::string &type_in,
		       const std::string &label_in) {
	 pos_x = x_in;
	 pos_y = y_in;
	 pos_z = z_in;
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
      grid_index_t(int i, int j) {
	 ii_ = i;
	 jj_ = j;
      }
      grid_index_t() {} // hhmmmm!  needed to compile contour_fragment
			// constructor, which takes a const reference
			// to a grid_index_t as an argument. I don't
			// understand why this is needed - or if it
			// works.
      enum { INVALID_INDEX = -1 };
      
      bool is_valid_p() const {
	 return (ii_ != INVALID_INDEX);
      }
      int i() const { return ii_;}
      int j() const { return jj_;}
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
      void plot_contour_lines(const std::vector<std::vector<lig_build::pos_t> > &contour_lines, GooCanvasItem *root) const;
      double substitution_value(double r_squared, double bash_dist) const;


      
   public:
      // (low means low numbers, not low on the canvas)
      // 
      ligand_grid(const lig_build::pos_t &low_x_and_y,
		  const lig_build::pos_t &high_x_and_y);

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

      void add_for_accessibility(double bash_dist, const lig_build::pos_t &atom_pos);

      // fudge in a smoothly varing function (so that the contouring
      // behaves smoothly, rather that the jaggies that we'd get if we
      // added a top hat function values to 1.0A.
      // 
      void add_for_accessibility_no_bash_dist_atom(double scale, const lig_build::pos_t &atom_pos);

      void show_contour(GooCanvasItem *root, float contour_level) const;
      void show_contour(GooCanvasItem *root, float contour_level,
			const std::vector<widgeted_atom_ring_centre_info_t> &unlimited_atoms) const;

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


private:
   bool stand_alone_flag;

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
   std::vector<widgeted_molecule_t> previous_molecules;
   int save_molecule_index;
   bool in_delete_mode_;
   highlight_data_t highlight_data;
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
   bool try_add_or_modify_bond(int canvas_addition_mode, int x, int y); //  ditto.
   bool add_bond_to_atom(int atom_index, int canvas_addition_mode);
   void add_bond_to_atom_with_0_neighbours(int atom_index, int canvas_addition_mode);
   void add_bond_to_atom_with_1_neighbour(int atom_index, int canvas_addition_mode,
					  int bond_index);
   void add_bond_to_atom_with_2_neighbours(int atom_index, int canvas_addition_mode,
					   const std::vector<int> &bond_indices);
   void add_bond_to_atom_with_3_neighbours(int atom_index, int canvas_addition_mode,
					   const std::vector<int> &bond_indices);
   std::string to_element(int addition_mode) const;
   std::string font_colour(int addition_element_mode) const;
   std::string font_colour(const std::string &ele) const;
   lig_build::bond_t::bond_type_t addition_mode_to_bond_type(int canvas_addition_mode) const;
   void try_stamp_bond_anywhere(int canvas_addition_mode, int x_mouse, int y_mouse); // always modifies.
   bool change_atom_element(int atom_index, std::string new_element, std::string fc);
   void change_atom_id_maybe(int atom_index);
   lig_build::pos_t mouse_at_click;
   void save_molecule();
   std::string mdl_file_name; // for save function.
   void add_search_combobox_text() const;
   bool make_saves_mutex;
   double canvas_scale; 
   void init_internal() {
      in_delete_mode_ = 0;
      save_molecule_index = UNASSIGNED_INDEX;
      make_saves_mutex = 1; // allow saves initially
      search_similarity = 0.93;
      coot_mdl_ready_time = 0;
      canvas_scale = 1.0;
      canvas_drag_offset = lig_build::pos_t(0,0);
      standard_residue_circle_radius = 19;
      button_down_bond_addition = 0;
      latest_bond_canvas_item = 0;
      penultimate_atom_index = -1;
      ultimate_atom_index = -1;
      latest_bond_was_extended = 0;
      stand_alone_flag = 0;
   }
   
   // return a status and a vector of atoms (bonded to atom_index) having
   // only one bond.
   // 
   std::pair<bool, std::vector<int> > 
   have_2_stubs_attached_to_atom(int atom_index, const std::vector<int> &bond_indices) const;
   void squeeze_in_a_4th_bond(int atom_index, int canvas_addition_mode,
			      const std::vector<int> &bond_indices);
   std::vector<double>
   get_angles(int atom_index, const std::vector<int> &bond_indices) const;
   lig_build::pos_t  new_pos_by_bisection(int atom_index,
					  const std::vector<int> &bond_indices,
					  const std::vector<double> &angles,
					  GooCanvasItem *root) const;
   bool all_closed_rings(int atom_index, const std::vector<int> &bond_indices) const;
   std::vector<lig_build::pos_t>
   get_centres_from_bond_indices(const std::vector<int> &bond_indices) const;
   lig_build::pos_t get_new_pos_not_towards_ring_centres(int atom_index,
							 const std::vector<int> &bond_indices) const;
   CMMDBManager *get_cmmdbmanager(const std::string &filename) const;

   // sbase functions
   CSBase *SBase;
   float search_similarity;
   int init_sbase(const std::string &sbase_monomer_dir_in);
   std::vector<match_results_t> compare_vs_sbase(CGraph *graph1,
						 float similarity,
						 int n_vertices) const;
   // return 0 on strangeness, to pass in search.
   // 

   int get_min_match(const int &n1) const {
      int most_1 = int (search_similarity * float(n1));
      return most_1;
      // return n1;
   }

   // not used.
   int get_min_match(const int &n1, const int &n2) const {
      int r = 0;
      int most_1 = int (search_similarity * float(n1));
      int most_2 = int (search_similarity * float(n2));
      if ((n2>=most_1) && (n1>=most_2))  {
	 r = (most_2 > most_1) ? most_2 : most_1;
      }
      return r;
   }
   match_results_t residue_from_best_match(CGraph &graph_1, CGraph &graph_2,
					   CGraphMatch &match, int n_match, 
					   CSBStructure *SBS) const;
   void display_search_results(const std::vector<lbg_info_t::match_results_t> v) const;
   void rotate_latest_bond(int x, int y);
   void rotate_or_extend_latest_bond(int x, int y);
   void extend_latest_bond(); // use hilight_data

   // return the bond between the ligand "core" and the atom_index
   // (one of bond_indices)
   // 
   widgeted_bond_t orthogonalise_2_bonds(int atom_index,
					 const std::vector<int> &attached_bonds,
					 const std::vector<int> &bond_indices);
   std::vector<residue_circle_t> residue_circles;

   widgeted_molecule_t translate_molecule(const lig_build::pos_t &delta);
   void translate_residue_circles(const lig_build::pos_t &delta);
   void clear_canvas();
   void draw_residue_circle_top_layer(const residue_circle_t &rc, const lig_build::pos_t &p);
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
   void draw_substitution_contour();
   void draw_bonds_to_ligand();
   void draw_solvent_exposure_circle(const residue_circle_t &residue_circle,
				     const lig_build::pos_t &ligand_centre,
				     GooCanvasItem *group);
   std::string get_residue_solvent_exposure_fill_colour(double radius_extra) const;
   void draw_annotated_stacking_line(const lig_build::pos_t &ligand_ring_centre,
				     const lig_build::pos_t &residue_pos,
				     int stacking_type);
   void draw_all_residue_attribs();
   void handle_read_draw_coords_mol_and_solv_acc(const std::string &coot_pdb_file,
						 const std::string &coot_mdl_file,
						 const std::string &sa_file);
   // for debugging the minimization vs original positions
   std::vector<residue_circle_t> offset_residues_from_orig_positions(); 
   void show_grid(const lbg_info_t::ligand_grid &grid);
   void show_mol_ring_centres(); // not const because mol.get_ring_centres() caches
   std::string grid_intensity_to_colour(int val) const;
   std::string sixteen_to_hex_let(int v) const;
   void reposition_problematics_and_reoptimise(const std::vector<int> &problematics,
					       const std::vector<int> &primary_indices);
   std::vector<residue_circle_t>
   filter_residue_waters(const std::vector<residue_circle_t> &r_in,
			 double max_dist_water_to_ligand_atom,
			 double max_dist_water_to_protein_atom) const;

   PCGraph makeTestQueryGraph() const;  // debugging

#ifdef MAKE_ENTERPRISE_TOOLS
   RDKit::RWMol rdkit_mol(const widgeted_molecule_t &mol) const;
   RDKit::Bond::BondType convert_bond_type(const lig_build::bond_t::bond_type_t &t) const;
   std::string get_smiles_string_from_mol_rdkit() const;
#endif
   std::string get_smiles_string_from_mol_openbabel() const;




public:
   lbg_info_t(GtkWidget *canvas_in) {
      canvas = canvas_in;
      init_internal();
   }
   lbg_info_t() { init_internal(); }
   // toggle button modes, mutually exclusive
   enum { NONE, TRIANLE, SQUARE, PENTAGON, HEXAGON, HEXAGON_AROMATIC, HEPTAGON, OCTAGON,
	  ATOM_C, ATOM_N, ATOM_O, ATOM_S, ATOM_P, ATOM_F, ATOM_CL, ATOM_I, ATOM_BR, ATOM_X,
	  CHARGE, ADD_SINGLE_BOND, ADD_DOUBLE_BOND, ADD_TRIPLE_BOND, ADD_STEREO_OUT_BOND,
	  DELETE_MODE};
   void init(GtkBuilder *builder);
   GtkWidget *lbg_window;
   GtkWidget *about_dialog; 
   GtkWidget *search_combobox;
   GtkWidget *open_dialog;
   GtkWidget *save_as_dialog;
   GtkWidget *lbg_export_as_pdf_dialog;
   GtkWidget *lbg_export_as_png_dialog;
   GtkWidget *lbg_sbase_search_results_dialog;
   GtkWidget *lbg_sbase_search_results_vbox;
   GtkWidget *lbg_smiles_dialog;
   GtkWidget *lbg_smiles_entry;
   GtkWidget *lbg_statusbar;
   GtkWidget *lbg_toolbar_layout_info_label;
   GtkWidget *canvas;
   std::map<std::string, GtkToggleToolButton *> widget_names;
   widgeted_molecule_t mol;
   int canvas_addition_mode;
   bool save_togglebutton_widgets(GtkBuilder *builder);
   void handle_item_add(GdkEventButton *event);
   void handle_item_delete(GdkEventButton *event);
   void untoggle_others_except(GtkToggleToolButton *button_toggled_on);
   bool item_highlight_maybe(int x, int y);
   void highlight_bond(const lig_build::bond_t &bond, bool delete_mode);
   void highlight_atom(const lig_build::atom_t &atom, int atom_index, bool delete_mode);
   void remove_bond_and_atom_highlighting();
   void clear_button_down_bond_addition();
   void handle_drag(GdkModifierType state, int x, int y); // could be rotate bond or drag canvas
   void set_in_delete_mode(bool v) {
      in_delete_mode_ = v;
   }
   bool in_delete_mode_p() const { return in_delete_mode_; }
   double radius(int n_edges) const; // depends on zoom? (for future).
   void clear();
   std::string get_stroke_colour(int i, int n) const;
   void drag_canvas(int mouse_x, int mouse_y);
   void write_pdf(const std::string &file_name) const;
   void write_png(const std::string &file_name) const;
   void set_mouse_pos_at_click(int xpos, int ypos) {
      mouse_at_click = lig_build::pos_t(double(xpos), double(ypos));
   }
   void render_from_molecule(const widgeted_molecule_t &mol_in);
   void undo();
   void search() const;
   widgeted_molecule_t  import(const lig_build::molfile_molecule_t &mol_in,
			       const std::string &filename,
			       CMMDBManager *pdb_mol);
   static void on_sbase_search_result_button_clicked(GtkButton *button, gpointer user_data);
   static gboolean watch_for_mdl_from_coot(gpointer user_data);
   time_t coot_mdl_ready_time;
   void read_draw_residues(const std::string &file_name);
   std::vector<residue_circle_t> read_residues(const std::string &file_name) const;
   void draw_residue_circles(const std::vector<residue_circle_t> &v);
   void draw_stacking_interactions(const std::vector<residue_circle_t> &v);
   std::vector<solvent_accessible_atom_t> 
   read_solvent_accessibilities(const std::string &file_name) const;
   void read_files_from_coot();
   std::string file_to_string(const std::string &filename) const;
   std::string get_flev_analysis_files_dir() const;
   void show_key();
   void hide_key();
   void update_statusbar_smiles_string() const;
   // can throw an exception
   std::string get_smiles_string_from_mol() const;
   void set_stand_alone() { stand_alone_flag = 1; }
   bool is_stand_alone() { return stand_alone_flag; } 
};

#endif // LBG_HH
