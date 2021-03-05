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

#include <vector>
#include <string>
#include <iostream>

#include <gtk/gtk.h>
#include <goocanvas.h>
#include "lidia-core/lig-build.hh"

namespace coot {

   class graph_trace_info_t {
      std::vector<std::pair<double, double> > data;
   public:
      enum {PLOT_TYPE_LINE, PLOT_TYPE_BAR, PLOT_TYPE_SMOOTHED_LINE, PLOT_TYPE_SCATTER};
      int plot_type;
      std::string colour;
      double line_width;
      bool dashed;
      bool y_data_are_ints;
      graph_trace_info_t () {
	 dashed = false;
	 y_data_are_ints = false;
	 plot_type = PLOT_TYPE_LINE;
      }
      void set_data(const std::vector<std::pair<double, double> > &data_in);
      const std::vector<std::pair<double, double> > &get_data() const { return data; }
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
      class annotation_box_t {
      public:
	 annotation_box_t() {} // needed because we derive a class from this
	                       // and initial that class has nothing. Maybe
	                       // is_set could go here.
	 annotation_box_t(const lig_build::pos_t &top_left_in,
			  const lig_build::pos_t &bottom_right_in,
			  const std::string &outline_colour_in,
			  double line_width_in,
			  const std::string &fill_colour_in) {
	    top_left = top_left_in;
	    bottom_right = bottom_right_in;
	    line_width = line_width_in;
	    fill_colour = fill_colour_in;
	    outline_colour = outline_colour_in;
	 }
	 lig_build::pos_t top_left;
	 lig_build::pos_t bottom_right;
	 std::string outline_colour;
	 double line_width;
	 std::string fill_colour;
      };

      GooCanvas *canvas;
      GtkWidget *dialog;
      bool dark_mode; // means light text
      bool use_dialog_for_graph_flag;
      int dialog_width_orig;
      int dialog_height_orig;
      int dialog_width; // variable
      int dialog_height;// variable
      std::vector<GooCanvasItem *> items;
      std::vector<annotation_line_t> annotation_lines;
      std::vector<annotation_text_t> annotation_texts;
      std::vector<annotation_box_t>  annotation_boxes;
      std::string title_string;
      bool   extents_x_are_set;
      bool   extents_y_are_set;
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
      bool draw_x_axis_flag;
      bool draw_y_axis_flag;
      bool draw_x_ticks_flag;
      bool draw_y_ticks_flag;
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
	    if (trace_id >= int(traces.size())) {
	       return false;
	    } else {
	       return true;
	    }
	 }
      }
      void plot_bar_graph(int trace_id);
      void plot_line_graph(int trace_id);
      void plot_smoothed_line_graph(int trace_id);
      void plot_scatter_plot(int trace_id);
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

	 lig_build::pos_t pos((canvas_offset_x + (p.x-extents_min_x)*data_scale_x) * dialog_width / dialog_width_orig,
			      (canvas_offset_y - (p.y-extents_min_y)*data_scale_y) * dialog_height/ dialog_height_orig);
	 if (false)
	    std::cout << "world_to_canvas() input: " << p << " returning " << pos << std::endl;
	 return pos;
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
      void draw_annotation_boxes();
      void draw_contour_level_bar();
      void set_data_scales(int axis);
   public:
      enum {X_AXIS, Y_AXIS};
      enum {MAJOR_TICK, MINOR_TICK};
      class annotation_box_info_t {
      public:
	 annotation_box_info_t() { is_set = false; }
	 annotation_box_info_t(float x, float y) : x_mouse(x), y_mouse(y) {
	    current_contour_level = -1;
	    func = 0;
	    is_set = false;
	 }
	 annotation_box_info_t(int imol_in, float rmsd_in, float cc_in, void (*func_in)(int, float)) {
	    imol = imol_in;
	    rmsd = rmsd_in;
	    current_contour_level = cc_in;
	    func = func_in;
	    is_set = true;
	 }
	 bool is_set;
	 float x_mouse;
	 float y_mouse;
	 int imol;
	 float current_contour_level;
	 float rmsd;
	 void (*func)(int, float); // to set the contour level
      };
      class annotated_box_info_t : public annotation_box_t {
      public:
	 annotated_box_info_t() : abi() {}
	 annotated_box_info_t(const lig_build::pos_t &pos_top_left,
			      const lig_build::pos_t &pos_bottom_right,
			      const std::string &outline_colour,
			      double line_width,
			      const std::string &fill_colour,
			      int imol,
			      float rmsd,
			      float current_contour_level,
			      void (*func)(int, float)) : annotation_box_t(pos_top_left,
									   pos_bottom_right,
									   outline_colour,
									   line_width,
									   fill_colour),
							  abi(imol, rmsd, current_contour_level, func) {}
	 annotation_box_info_t abi;
      };
      goograph() {
	 init();
	 init_widgets();
      }
      goograph(int width, int height) {
	 init();
	 init_widgets();
	 dialog_width = width;
	 dialog_height = height;
	 dialog_width_orig = width;
	 dialog_height_orig = height;

	 canvas_offset_x = 10.0; // how much is the Y axis displaced
                                 // from the left-hand edge of the canvas?

	 // if drawing of the x axis has been turned off, then we want
	 // less space at the bottom, but for now, we will
	 // make it 8% offset up.
	 int x_off = 0.08 * height;
	 canvas_offset_y = height - x_off;; // how much is the size of the canvas
      }
      void clear();
      void draw_graph();
      void set_trace_type(int trace_id, int plot_type, bool dashed=false);
      void set_trace_colour(int trace, const std::string &colour);
      GtkWidget *get_canvas() const; // for embedding in other windows
      GtkWidget *show_dialog();            // for graph in dialog, return the close button so that we
                                           // can add a callback that NULLs the pointer to a goograph
      void set_extents(int axis, double min, double max); 
      void set_ticks(int axis, double tick_major, double tick_minor);
      void set_axis_label(int axis, const std::string &label);
      void set_draw_axis(int axis, bool draw_state);
      void set_draw_ticks(int axis, bool draw_state);
      void set_plot_title(const std::string &title);

      void set_use_dialog_for_graph(bool flag); // default true currently
      void enable_dark_mode(bool flag);

      // if using auto-extents, then the first use of set_data() sets the extents
      // (they are not over-ridden by subsequent set_data() usage)
      //
      void set_data(int trace_id, const std::vector<std::pair<double, double> > &data);
      int trace_new();
      void plot_trace(int trace_id);
      std::pair<double, double> min_max_x() const {
	 return std::pair<double, double> (extents_min_x, extents_max_x);
      }
      std::pair<double, double> min_max_y() const {
	 return std::pair<double, double> (extents_min_y, extents_max_y);
      }
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
      void add_annotation_box(const lig_build::pos_t &pos_top_left,
			      const lig_build::pos_t &pos_bottom_right,
			      const std::string &outline_colour,
			      double line_width,
			      const std::string &fill_colour);
      void add_contour_level_box(float contour_level,
				 const std::string &outline_colour,
				 double line_width,
				 const std::string &fill_colour,
				 int imol, float rmsd,
				 void (*func)(int, float));

      void clear_traces_and_annotations();

      annotated_box_info_t contour_level_bar; // a special item that has a callback

      static
      bool on_goograph_active_button_box_press_event(GooCanvasItem  *item,
						     GooCanvasItem  *target_item,
						     GdkEventButton *event,
						     gpointer        user_data);

      static
      bool on_goograph_active_button_box_release_event(GooCanvasItem  *item,
						       GooCanvasItem  *target_item,
						       GdkEventButton *event,
						       gpointer        user_data);

   };
}

#endif // HAVE_GOOGRAPH_HH

