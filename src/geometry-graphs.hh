/* Sec/geometry-graphs.hh
 * 
 * Copyright 2004, 2005 by The University of York
 * Author: Paul Emsley
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#ifndef HAVE_GEOMETRY_GRAPHS_HH
#define HAVE_GEOMETRY_GRAPHS_HH

#include "goograph/goograph.hh"

#include "ideal/simple-restraint.hh"

namespace coot {

   enum geometry_graph_type {GEOMETRY_GRAPH_GEOMETRY,
			     GEOMETRY_GRAPH_B_FACTOR,
			     GEOMETRY_GRAPH_DENSITY_FIT,
			     GEOMETRY_GRAPH_OMEGA_DISTORTION,
			     GEOMETRY_GRAPH_ROTAMER,
			     GEOMETRY_GRAPH_NCS_DIFFS,
			     SEQUENCE_VIEW,
			     RAMACHANDRAN_PLOT,
////B
			     GEOMETRY_GRAPH_CALC_B_FACTOR
////E
   };

   gint on_geometry_graph_block_clicked(GooCanvasItem *item,
                                        GooCanvasItem *target,
					GdkEvent *event, 
					gpointer data);

   class geometry_graph_block_info_generic { 
   public:
      int imol;
      int resno;
      atom_spec_t atom_spec;
      double distortion;
      std::string distortion_info_string;
      geometry_graph_block_info_generic(int imol_in, 
					int resno_in, 
					const atom_spec_t &atom_spec_in, 
					double dist_in, 
					const std::string &str) {
	 resno = resno_in;
	 atom_spec = atom_spec_in;
	 distortion = dist_in;
	 distortion_info_string = str;
	 imol = imol_in;
      }
   };
   
   class geometry_graph_block_info : public geometry_graph_block_info_generic { 
   public:
      GooCanvas *canvas;  // for button press callback
      geometry_graph_block_info(int imol_in, 
				int resno_in, 
				const atom_spec_t &atom_spec_in, 
				double dist_in, 
				const std::string &str, 
				GooCanvas *w) :
	 geometry_graph_block_info_generic(imol_in, resno_in, atom_spec_in, dist_in, str)
      {
	 canvas = w;
      }
   };
   
   class b_factor_block_info_t {
   public:
      int resno;
      double b_factor_var;
      std::string atom_name;
      std::string info_string;
      b_factor_block_info_t() {}  // vector
      b_factor_block_info_t(int resno_in, double b_factor_var_in, const std::string &s) {
	 resno = resno_in;
	 b_factor_var = b_factor_var_in;
	 info_string = s;
      }
   };
   

   class geometry_graphs {
      int imol; // for block click callback (and unset widget from
		// graphics info on destroy)
      geometry_graph_type graph_type;
      int n_chains;
      int max_chain_length;
      std::string graph_label;
      GooCanvas *canvas;
      void clear_tooltip_box();
      void tooltip_like_box(const geometry_graph_block_info &bi,
			    const GdkEvent *event);
      GooCanvasItem *tooltip_item;
      GooCanvasItem *tooltip_item_text;
      
      std::vector<std::string> chain_index; // could be a map
      std::vector<int> offsets; // tricky indexing.  If we have a
				// residue number of 1 in the chain,
				// offset is 0.  If we have 10, it is 9.

      std::vector<std::vector<GooCanvasItem *> > blocks;
      void setup_internal();
      void setup_canvas(int n_chains, int max_chain_length);
      void plot_block(const geometry_graph_block_info &block_info,
		      int offset, int chain_number);
      void plot_blocks(const std::map<coot::residue_spec_t, std::pair<double, std::string> > &residue_distortions,
		       int chain_number);
      std::string make_distortion_string(const coot::geometry_distortion_info_t &extra_distortion,
					 const coot::geometry_distortion_info_container_t &dc) const;

      void label_chain(const std::string &label, int ichain) const;
      void draw_chain_axis(int nres, int ichain) const;
      void draw_chain_axis_tick_and_tick_labels(int min_resno, int max_resno, int chain_number) const;
      void delete_block(int chain_number, int resno);
      void delete_block(int chain_number, const coot::residue_spec_t &rs);
      std::vector<std::string> colour_list;
      std::string distortion_to_colour(const double &distortion) const;
      double distortion_max;
      std::string fixed_font_str;
      int chain_id_to_chain_index(const std::string &chain_id) const;

      void render_geometry_distortion_blocks_internal(const geometry_distortion_info_container_t &dc);
      void render_geometry_distortion_blocks_internal_linear(const geometry_distortion_info_container_t &dc,
				  int min_res, int max_resno);
      std::vector<geometry_distortion_info_container_t>
      geometric_distortions_from_mol(const atom_selection_container_t &asc) const;
      GtkBuilder *builder;

      double sane_occupancy(const double &occ_in) {
	 // SHELX atoms have occupancies that depend on FVARS and can
	 // be things like -81 and 61.  Often they are 11.0, which
	 // means fixed at 1.0.  So, if they are depending on FVARS,
	 // presume dual-occupancy and return 0.5.
	 double occ = occ_in;
	 if ((occ_in < 11.1) && (occ_in > 10.9)) { // like SHELX
	    occ = 1.0;
	 } else { 
	    if (occ_in > 1.0) { 
	       occ = 1.0;
	    } else {
	       if (occ_in < -1) {// SHELX negative
		  occ = 0.5;
	       } else {
		  if (occ_in < 0.0) {
		     occ = 0.0;
		  }
	       }
	    }
	 }
	 return occ;
      }

      GtkWidget *widget_from_builder(const std::string &wn) const {
         if (builder)
            return GTK_WIDGET(gtk_builder_get_object(builder, wn.c_str()));
         else
            return nullptr;
      }

      GtkWidget *create_geometry_graphs_dialog_gtk3();

   public:
      geometry_graphs(geometry_graph_type type, int imol, std::string graph_label,
		      int nchains, int max_chain_length);

      int get_imol() const { return imol; }
      void render_to_canvas(const geometry_distortion_info_container_t &dv, 
			    int chain_number);
      void render_to_canvas(const geometry_distortion_info_container_t &dv, 
			    const std::string &chain_id_in) {
	 int ch =  chain_id_to_chain_index(chain_id_in);
	 if (ch >= 0) 
	    render_to_canvas(dv, ch);
      }
      void render_to_canvas(const std::vector<geometry_graph_block_info_generic> &gbi,
			    int ichain,
			    const std::string &chain_id,
			    int max_resno,
			    int min_resno,
			    int offset);
      void render_b_factor_blocks(int imol,
				  int ichain,
				  const std::string  &chain_id,
				  int offset, 
				  const std::vector<b_factor_block_info_t> &biv);

      void render_omega_blocks(const omega_distortion_info_container_t &om_dist,
			       int chain_number,
			       const std::string &chain_id,
			       int offset_in);
      void update_omega_blocks(const omega_distortion_info_container_t &om_dist,
			       int chain_number,
			       const std::string &chain_id);
      
      void update_residue_blocks(const geometry_distortion_info_container_t &dv);
      // 20100211 old style linear
      void update_residue_blocks_linear(const geometry_distortion_info_container_t &dv); 
      void update_residue_blocks(const std::vector<geometry_graph_block_info_generic> &dv);
      
      void set_sensitivity(double d); // becomes the distortion_max.
      void button_press(GooCanvasItem *item, GdkEvent *event, 
			const geometry_graph_block_info &binfo);
      void mouse_over(GooCanvasItem *item, GdkEvent *event, 
		      const geometry_graph_block_info &binfo);
      void delete_block(const std::string &chain_id, int resno);
      GtkWidget *dialog() const;
      int Imol() const { return imol; }
      geometry_graph_type Graph_Type() const { return graph_type; }
      void close_yourself();

      static void density_fit_rescale_button_callback(GtkButton *button, gpointer user_data);

   };
}


#endif // HAVE_GEOMETRY_GRAPHS_HH

