/* src/goograph.cc
 * 
 * Copyright 2011, 2012 by The University of Oxford
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

#ifndef HAVE_GOOGRAPH_HH
#define HAVE_GOOGRAPH_HH

#ifdef HAVE_GOOCANVAS

#include <vector>
#include <string>
#include <iostream>

#include <gtk/gtk.h>
#include <goocanvas.h>
#include "lig-build.hh"

namespace coot {

   class graph_trace_info_t {
   public:
      enum {PLOT_TYPE_LINE, PLOT_TYPE_BAR, PLOT_TYPE_SMOOTHED_LINE};
      int plot_type;
      std::string colour;
      double line_width;
      bool dashed;
      graph_trace_info_t () {
	 dashed = false;
	 plot_type = PLOT_TYPE_LINE;
      }
      std::vector<std::pair<double, double> > data;
   };

   class goograph {
      class annotation_line_t {
      public:
	 annotation_line_t(lig_build::pos_t pos_1_in,
			   lig_build::pos_t pos_2_in,
			   std::string colour_in,
			   double line_width_in,
			   bool dashed_in,
			   bool start_arrow_in,
			   bool end_arrow_in) {
	    pos_1 = pos_1_in;
	    pos_2 = pos_2_in;
	    colour = colour_in;
	    line_width = line_width_in;
	    dashed = dashed_in;
	    start_arrow = start_arrow_in;
	    end_arrow = end_arrow_in;
	 }
	 lig_build::pos_t pos_1;
	 lig_build::pos_t pos_2;
	 std::string colour;
	 double line_width;
	 bool dashed;
	 bool start_arrow;
	 bool end_arrow;
      };
      class annotation_text_t {
      public:
	 annotation_text_t(std::string text_in,
			   lig_build::pos_t pos_in,
			   std::string colour_in,
			   std::string font_in) {
	    text = text_in;
	    pos = pos_in;
	    colour = colour_in;
	    font = font_in;
	 }
	 std::string text;
	 lig_build::pos_t pos;
	 std::string colour;
	 std::string font;
      };

      GooCanvas *canvas;
      GtkWidget *dialog;
      int dialog_width_orig;
      int dialog_height_orig;
      int dialog_width; // variable
      int dialog_height;// variable
      std::vector<GooCanvasItem *> items;
      std::vector<annotation_line_t> annotation_lines;
      std::vector<annotation_text_t> annotation_texts;
      std::string title_string;
      double extents_min_x;
      double extents_max_x;
      double extents_min_y;
      double extents_max_y;
      double tick_major_x;
      double tick_major_y;
      double tick_minor_x;
      double tick_minor_y;
      double canvas_offset_x; 
      double canvas_offset_y; 
      std::string x_axis_label;
      std::string y_axis_label;
      std::string plot_title;
      std::vector<graph_trace_info_t> traces;
      std::string dark; // axis colour
      double data_scale_x;
      double data_scale_y;
      void init();
      void init_widgets();
      bool is_valid_trace(int trace_id) {
	 if (trace_id < 0) { 
	    return false;
	 } else {
	    if (trace_id >= traces.size()) {
	       return false;
	    } else {
	       return true;
	    }
	 }
      }
      void plot_bar_graph(int trace_id);
      void plot_line_graph(int trace_id);
      void plot_smoothed_line_graph(int trace_id);
      static void goograph_close_callback(GtkWidget *button,
					  GtkWidget *dialog);
      static gint reshape(GtkWidget *widget, GdkEventConfigure *event);

      void set_bounds(double right, double bottom);
      void draw_axes();
      void draw_ticks();
      void draw_data();
      void draw_y_ticks(int tick_type, double tick_major_step, double tick_step, double tick_length_multiplier);
      void draw_x_ticks(int tick_type, double tick_major_step, double tick_step, double tick_length_multiplier);
      void draw_ticks_generic(int axis, int tick_type,
			      double tick_major_step,
			      double tick_step,
			      double tick_length_multiplier);
      lig_build::pos_t world_to_canvas(const lig_build::pos_t &p) const {
// 	 return lig_build::pos_t(canvas_offset_x + (p.x-extents_min_x)*data_scale_x,
// 				 canvas_offset_y - (p.y-extents_min_y)*data_scale_y);
	 return lig_build::pos_t((canvas_offset_x + (p.x-extents_min_x)*data_scale_x) * dialog_width / dialog_width_orig,
				 (canvas_offset_y - (p.y-extents_min_y)*data_scale_y) * dialog_height/ dialog_height_orig);
      }
      void draw_title();
      double y_range() const { return extents_max_y - extents_min_y; }
      double x_range() const { return extents_max_x - extents_min_x; }
      bool close_float_p(const double &f1, const double &f2) {
	    return (fabs(f1-f2) < 0.0001);
      }
      double median_bin_width(int trace_id) const;
      double calc_tick(double range) const;
      void draw_axis_label(int axis);
      void draw_annotation_lines();
      void draw_annotation_texts();
      void set_data_scales();
   public:
      enum {X_AXIS, Y_AXIS};
      enum {MAJOR_TICK, MINOR_TICK};
      goograph() {
	 init();
	 init_widgets();
	 extents_min_x =  9999999990.0;
	 extents_min_y =  9999999990.0;
	 extents_max_x = -9999999990.0;
	 extents_max_y = -9999999990.0;
	 tick_major_x = 0.1;
	 tick_minor_x = 0.05;
	 tick_major_y = 0.1;
	 tick_minor_y = 0.05;
	 dark = "#111111";
	 canvas_offset_x = 70.0; // how much is the Y axis displaced
				 // from the left-hand edge of the
				 // canvas?
	 canvas_offset_y = 370.0; // how much is the size of the //
	                          // canvas - and include an offset of
	                          // // the axis from the bottom edge.
	                          // (the smaller the number the
	                          // greater the displacement (without
	                          // resizing).
	 data_scale_x = 1.0;
	 data_scale_y = 1.0;
      }
      void clear();
      void draw_graph();
      void set_trace_type(int trace_id, int plot_type, bool dashed=false);
      void set_trace_colour(int trace, const std::string colour);
      GtkWidget *get_canvas() const; // for embedding in other windows
      GtkWidget *show_dialog();            // for graph in dialog, return the close button so that we
                                           // can add a callback that NULLs the pointer to a goograph
      void set_extents(int axis, double min, double max); 
      void set_ticks(int axis, double tick_major, double tick_minor);
      void set_axis_label(int axis, const std::string &label);
      void set_plot_title(const std::string &title);
      void set_data(int trace_id, const std::vector<std::pair<double, double> > &data);
      int trace_new();
      void plot_trace(int trace_id);
      void add_annotation_line(const lig_build::pos_t &pos_1,
			       const lig_build::pos_t &pos_2,
			       const std::string &colour,
			       double line_width,
			       bool dashed_flag,
			       bool start_arrow,
			       bool end_arrow);
      // font can be "" meaning default ("Sans 9");
      // colour can be "" meaning standard dark colour
      void add_annotation_text(const std::string &text,
			       const lig_build::pos_t &pos_1,
			       const std::string &colour,
			       const std::string &font);
      void clear_traces_and_annotations();
   };
}

#endif // HAVE_GOOCANVAS
#endif // HAVE_GOOGRAPH_HH

