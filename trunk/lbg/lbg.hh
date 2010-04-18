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

static double SINGLE_BOND_CANVAS_LENGTH=35.27;

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
   public:
      double pos_x; // input coordinate reference frame
      double pos_y;
      double pos_z;
      lig_build::pos_t pos; // coordinate system of the ligand atoms
      std::string residue_type;
      std::string residue_label;
      std::vector<bond_to_ligand_t> bonds_to_ligand;
      residue_circle_t(const double &x_in, const double &y_in, const double &z_in,
		       const std::string &type_in,
		       const std::string &label_in) {
	 pos_x = x_in;
	 pos_y = y_in;
	 pos_z = z_in;
	 residue_type = type_in;
	 residue_label = label_in;
      }
      void set_canvas_pos(const lig_build::pos_t &pos_in) {
	 pos = pos_in;
      }
      void add_bond_to_ligand(const bond_to_ligand_t &bl) {
	 bonds_to_ligand.push_back(bl);
      }
      // friend std::ostream& operator<<(std::ostream &s, residue_circle_t r);
   };
   // std::ostream& operator<<(std::ostream &s, residue_circle_t r);

   class optimise_residue_circles {
   private:
      int status; 
      static double f(const gsl_vector *v, void *params);
      static void  df(const gsl_vector *v, void *params, gsl_vector *df); 
      static void fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
      
      std::vector<residue_circle_t> starting_circles;
      std::vector<residue_circle_t>  current_circles;
      widgeted_molecule_t mol;
      void numerical_gradients(gsl_vector *x, gsl_vector *df, void *params) const;
      static bool score_vs_ligand_atoms;
      static bool score_vs_ring_centres;
      static bool score_vs_other_residues;
      static bool score_vs_original_positions;
      
   public:
      // we pass two vectors here because (for trajectory-view) we
      // don't want to restart the minimisation with the current
      // positions - we want to minimise against the (constant)
      // original positions.
      optimise_residue_circles(const std::vector<residue_circle_t> &r, // starting points
			       const std::vector<residue_circle_t> &c, // current points
			       const widgeted_molecule_t &mol);
      std::pair<int, std::vector<residue_circle_t> > solution() const;
      int get_status() const { return status; }
   };


private:
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
   void orthogonalise_2_bonds(int atom_index,
			      const std::vector<int> &attached_bonds,
			      const std::vector<int> &bond_indices);
   std::vector<residue_circle_t> residue_circles;

   widgeted_molecule_t translate_molecule(const lig_build::pos_t &delta);
   void translate_residue_circles(const lig_build::pos_t &delta);
   void clear_canvas();
   void add_residue_circle(const residue_circle_t &rc, const lig_build::pos_t &p);

   std::pair<std::string, std::string>
   get_residue_circle_colour(const std::string &residue_type) const;
   // minimise layout energy
   std::pair<int, std::vector<residue_circle_t> >
   optimise_residue_circle_positions(const std::vector<residue_circle_t> &start,
				     const std::vector<residue_circle_t> &current) const;
   lig_build::pos_t canvas_drag_offset;
   bool is_sane_drag(const lig_build::pos_t &delta) const;

   void draw_solvent_accessibility_of_atom(const lig_build::pos_t &pos, double sa,
					   GooCanvasItem *root);
   void draw_bonds_to_ligand();
   
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
   GtkWidget *about_dialog; 
   GtkWidget *search_combobox;
   GtkWidget *open_dialog;
   GtkWidget *save_as_dialog;
   GtkWidget *lbg_export_as_pdf_dialog;
   GtkWidget *lbg_export_as_png_dialog;
   GtkWidget *lbg_sbase_search_results_dialog;
   GtkWidget *lbg_sbase_search_results_vbox;
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
   std::vector<solvent_accessible_atom_t> 
   read_solvent_accessibilities(const std::string &file_name) const;
};

#endif // LBG_HH
