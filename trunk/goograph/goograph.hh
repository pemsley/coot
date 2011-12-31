

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
      graph_trace_info_t () {
      }
      std::vector<std::pair<double, double> > data;
   };

   class goograph {
      GtkWidget *dialog;
      GooCanvas *canvas;
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
      void draw_graph();
      void draw_axes();
      void draw_ticks();
      void draw_data();
      void draw_y_ticks(int tick_type, double tick_step, double tick_length_multiplier);
      void draw_x_ticks(int tick_type, double tick_step, double tick_length_multiplier);
      void draw_ticks_generic(int axis, int tick_type,
			      double tick_step, double tick_length_multiplier);
      lig_build::pos_t world_to_canvas(const lig_build::pos_t &p) const {
	 return lig_build::pos_t(canvas_offset_x + (p.x-extents_min_x)*data_scale_x,
				 canvas_offset_y - (p.y-extents_min_y)*data_scale_y);
      }
      double y_range() const { return extents_max_y - extents_min_y; }
      double x_range() const { return extents_max_x - extents_min_x; }
      bool close_float_p(const double &f1, const double &f2) {
	    return (fabs(f1-f2) < 0.0001);
      }
      double median_bin_width(int trace_id) const;

   public:
      enum {X_AXIS, Y_AXIS};
      enum {PLOT_TYPE_LINE, PLOT_TYPE_BAR, PLOT_TYPE_SMOOTHED_LINE};
      enum {MAJOR_TICK, MINOR_TICK};
      goograph() {
	 init_widgets();
	 extents_min_x = 9999999990.0;
	 extents_min_y = 9999999990.0;
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
      void show_dialog();
      void set_extents(int axis, double min, double max); 
      void set_ticks(int axis, double tick_major, double tick_minor);
      void set_axis_label(int axis, const std::string &label);
      void set_plot_title(const std::string &title);
      void set_data(int trace_id, const std::vector<std::pair<double, double> > &data);
      int trace_new();
      void plot(int trace_id, int plot_type);
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
		    
   };
}

#endif // HAVE_GOOCANVAS
#endif // HAVE_GOOGRAPH_HH

