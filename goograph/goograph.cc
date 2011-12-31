
#ifdef HAVE_GOOCANVAS
#include <iostream>
#include <algorithm>
#include "goograph.hh"
#include "coot-utils.hh"


void
coot::goograph::show_dialog() {
   draw_graph();
   gtk_widget_show(dialog);

}

void
coot::goograph::init_widgets() {

   canvas = GOO_CANVAS(goo_canvas_new());

   dialog = gtk_dialog_new();
   gtk_window_set_default_size(GTK_WINDOW(dialog), 600, 500);
   gtk_object_set_data(GTK_OBJECT(dialog), "goograph", dialog);
   gtk_window_set_default_size(GTK_WINDOW(dialog), 600, 500);
   GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;
   GtkWidget *vbox_inner = gtk_vbox_new(FALSE, 2);
   GtkWidget *scrolled_window = gtk_scrolled_window_new (NULL, NULL);
   gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
					 GTK_WIDGET(vbox_inner));
   gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scrolled_window), TRUE, TRUE, 2);
   gtk_widget_show(scrolled_window);
   gtk_widget_show(vbox_inner);
   gtk_container_add(GTK_CONTAINER(vbox_inner), GTK_WIDGET(canvas));
   gtk_widget_show(GTK_WIDGET(canvas));
   GtkWidget *close_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Close", 2);
   gtk_widget_show(close_button);
   g_signal_connect(G_OBJECT(close_button), "clicked",
		    G_CALLBACK(goograph_close_callback),
		    (gpointer) dialog);

   double left = 0;
   double top = 0;
   double right = 1000;;
   double bottom = 1000;;
   goo_canvas_set_bounds(canvas, left, top, right, bottom); // we can't "invert" these.
   
}

// static
void
coot::goograph::goograph_close_callback(GtkWidget *button,
					GtkWidget *dialog) {
   gtk_widget_destroy(dialog);
   dialog = NULL;
} 


void
coot::goograph::draw_graph() {

   draw_axes();
   draw_ticks();
}

void
coot::goograph::draw_axes() {

   std::string col = "blue";
   std::string colour = "red";
   GooCanvasItem *root = goo_canvas_get_root_item(canvas);

   // debug item
//    GooCanvasItem *rect =
//       goo_canvas_rect_new(root, 150, 150, 30.5, 30.5,
// 			  "line-width", 1.0,
// 			  "fill_color", colour.c_str(),
// 			  "stroke-color", colour.c_str(),
// 			  NULL);
//    rect =
//       goo_canvas_rect_new(root, -5, -5, 50.5, 50.5,
// 			  "line-width", 1.0,
// 			  // "fill_color", "orange",
// 			  "fill-color-rgba", 0x5555cc30,
// 			  "stroke-color", "orange",
// 			  NULL);
   
//    double y = 5.5;
//    GooCanvasItem *h_line =
//       goo_canvas_polyline_new_line(root,
// 				   -20.0, y,
// 				   1000.4, y,
// 				   "line-width", 7.0,
// 				   "stroke-color", col.c_str(),
// 				   NULL);

   GooCanvasLineDash *dash=goo_canvas_line_dash_new (2, 5.7, 0.0);
   lig_build::pos_t A_yaxis(extents_min_x, extents_min_y);
   lig_build::pos_t B_yaxis(extents_min_x, extents_min_y + y_range()*1.1);
   lig_build::pos_t wAy = world_to_canvas(A_yaxis);
   lig_build::pos_t wBy = world_to_canvas(B_yaxis);
   lig_build::pos_t A_xaxis(extents_min_x, extents_min_y);
   lig_build::pos_t B_xaxis(extents_min_x + x_range()*1.1, extents_min_y);
   std::cout << "in draw_axes() X: extents " << extents_min_x << " "
	     << extents_max_x << " data_scale_x " << data_scale_x << std::endl;
   std::cout << "in draw_axes() Y: extents " << extents_min_y << " "
	     << extents_max_y << " data_scale_y " << data_scale_y << std::endl;
   
   lig_build::pos_t wAx = world_to_canvas(A_xaxis);
   lig_build::pos_t wBx = world_to_canvas(B_xaxis);
   gboolean start_arrow = 0;
   gboolean   end_arrow = 1;
   GooCanvasItem *item_1 = goo_canvas_polyline_new_line(root,
							wAy.x, wAy.y,
							wBy.x, wBy.y,
							"line-width", 2.5,
							"start_arrow", start_arrow,
							"end_arrow",   end_arrow,
							"stroke-color", dark.c_str(),
							NULL);
   GooCanvasItem *item_2 = goo_canvas_polyline_new_line(root,
							wAx.x, wAx.y,
							wBx.x, wBx.y,
							"line-width", 2.5,
							"start_arrow", start_arrow,
							"end_arrow",   end_arrow,
							"stroke-color", dark.c_str(),
							NULL);
} 

void
coot::goograph::draw_ticks() {

   draw_y_ticks(MINOR_TICK, tick_minor_y, 0.5);
   draw_y_ticks(MAJOR_TICK, tick_major_y, 1.0);

   draw_x_ticks(MINOR_TICK, tick_minor_x, 0.5);
   draw_x_ticks(MAJOR_TICK, tick_major_x, 1.0);
   
}

void
coot::goograph::draw_y_ticks(int tick_type, double tick_step, double tick_length_multiplier) {
   draw_ticks_generic(Y_AXIS, tick_type, tick_step, tick_length_multiplier);
}

void
coot::goograph::draw_x_ticks(int tick_type, double tick_step, double tick_length_multiplier) {
   draw_ticks_generic(X_AXIS, tick_type, tick_step, tick_length_multiplier);
}

void
coot::goograph::draw_ticks_generic(int axis, int tick_type,
				   double tick_step, double tick_length_multiplier) {

   GooCanvasItem *root = goo_canvas_get_root_item(canvas);
   double extents_min = 0.0;
   double extents_max = 0.0;
   double tick_major_step = 0;

   if (axis == Y_AXIS) {
      extents_min = extents_min_y;
      extents_max = extents_max_y;
      tick_major_step = tick_major_y;
   } 
   if (axis == X_AXIS) {
      extents_min = extents_min_x;
      extents_max = extents_max_x;
      tick_major_step = tick_major_x;
   } 
   // 
   for (double tick_pos = extents_min; tick_pos <= extents_max; tick_pos+=tick_step) {

      bool do_tick = true;
      if (tick_type == MINOR_TICK) {
	 // if the minor tick overlaps a major tick then don't display
	 // it (set do_tick to false).
	 for (double tick_pos_inner = extents_min;
	      tick_pos_inner <= extents_max;
	      tick_pos_inner+=tick_major_step) {
	    if (close_float_p(tick_pos_inner, tick_pos)) {
	       do_tick = false;
	       break;
	    }
	 }
      }
      if (do_tick) {
	 lig_build::pos_t A_axis;
	 lig_build::pos_t B_axis;
	 double tick_label_x_off = 0;
	 double tick_label_y_off = 0;
	 if (axis == Y_AXIS) { 
	    A_axis = lig_build::pos_t(extents_min_x, tick_pos);
	    B_axis = lig_build::pos_t(extents_min_x -x_range()*0.04*tick_length_multiplier, tick_pos);
	    tick_label_x_off = 10;
	    tick_label_y_off = 0;
	 }
	 if (axis == X_AXIS) { 
	    A_axis = lig_build::pos_t(tick_pos, extents_min_y);
	    B_axis = lig_build::pos_t(tick_pos, extents_min_y - y_range()*0.06*tick_length_multiplier);
	    tick_label_x_off = 7;
	    tick_label_y_off = -13;
	 }
	 lig_build::pos_t wA = world_to_canvas(A_axis);
	 lig_build::pos_t wB = world_to_canvas(B_axis);
      
	 GooCanvasItem *tick =
	    goo_canvas_polyline_new_line(root,
					 wA.x, wA.y,
					 wB.x, wB.y,
					 "line-width", 1.0,
					 "stroke-color", dark.c_str(),
					 NULL);

	 GtkAnchorType anchor_type = GTK_ANCHOR_CENTER;
	 if (axis == X_AXIS)
	    anchor_type = GTK_ANCHOR_SOUTH_WEST;

	 if (tick_type == MAJOR_TICK) {
	    int n_dec_pl = 0;
	    if (axis == Y_AXIS) { 
	       if (y_range() < 5)
		  n_dec_pl = 1;
	       if (y_range() < 2)
		  n_dec_pl = 2;
	    }
	    std::string txt =
	       coot::util::float_to_unspaced_string_using_dec_pl(tick_pos, n_dec_pl);
	    GooCanvasItem *text = goo_canvas_text_new(root, txt.c_str(),
						      wB.x - tick_label_x_off,
						      wB.y - tick_label_y_off,
						      -1,
						      anchor_type,
						      "font", "Sans 7",
						      "fill_color", dark.c_str(),
						      NULL);
	 }
      }
   }
}




void
coot::goograph::set_extents(int axis, double min, double max) {

   if (axis == X_AXIS) {
      extents_min_x = min;
      extents_max_x = max;
      double delta = extents_max_x - extents_min_x;
      if (delta > 0.1) 
	 data_scale_x = 400.0/delta;
      std::cout << "data_scale_x " << data_scale_x << std::endl;
   }
   if (axis == Y_AXIS) {
      extents_min_y = min;
      extents_max_y = max;
      double delta = extents_max_y - extents_min_y;
      if (delta > 0.1)
	 data_scale_y = 300/delta;
      std::cout << "data_scale_y " << data_scale_y << std::endl;
   }
}

void
coot::goograph::set_ticks(int axis, double tick_major, double tick_minor) {

   if (axis == X_AXIS) {
      tick_major_x = tick_major;
      tick_minor_x = tick_minor;
   } 
   if (axis == Y_AXIS) {
      tick_major_y = tick_major;
      tick_minor_y = tick_minor;
   } 

}


void
coot::goograph::set_axis_label(int axis, const std::string &label) {

   GooCanvasItem *root = goo_canvas_get_root_item(canvas);
   lig_build::pos_t A;
   if (axis == X_AXIS) {
      A = lig_build::pos_t(extents_min_x+x_range()*0.8, extents_min_y -y_range()*0.12);
   }
   if (axis == Y_AXIS) {
      A = lig_build::pos_t(extents_min_x-x_range()*0.14, y_range()*1.15);
   }
   lig_build::pos_t wA = world_to_canvas(A);
   GtkAnchorType anchor_type = GTK_ANCHOR_NORTH_WEST;
   GooCanvasItem *text =
      goo_canvas_text_new(root, label.c_str(),
			  wA.x, wA.y,
			  -1,
			  anchor_type,
			  "font", "Sans 9",
			  "fill_color", dark.c_str(),
			  NULL);
} 

void
coot::goograph::set_plot_title(const std::string &title) {
   if (dialog) 
      gtk_window_set_title (GTK_WINDOW(dialog), title.c_str());
} 


int
coot::goograph::trace_new() {

   coot::graph_trace_info_t tr;
   traces.push_back(tr);
   int n = traces.size();
   return n-1;
}

void
coot::goograph::set_data(int trace_id, const std::vector<std::pair<double, double> > &data_in) {

   if (is_valid_trace(trace_id)) {
      std::vector<std::pair<double, double> > sorted_data = data_in;
      std::sort(sorted_data.begin(), sorted_data.end());
      traces[trace_id].data = sorted_data;

      if (data_in.size()) { 
	 double max_x = -10000000;
	 double min_x =  10000000;
	 double max_y = -10000000;
	 double min_y =  10000000;
	 for (unsigned int i=0; i<data_in.size(); i++) { 
	    if (data_in[i].first < min_x)
	       min_x = data_in[i].first;
	    if (data_in[i].first > max_x)
	       max_x = data_in[i].first;
	    if (data_in[i].second < min_y)
	       min_y = data_in[i].second;
	    if (data_in[i].second > max_y)
	       max_y = data_in[i].second;
	 }
	 set_extents(X_AXIS, min_x, max_x);
	 set_extents(Y_AXIS, min_y, max_y);

	 std::cout << "in set_data: x_range() is " << x_range() << std::endl;
	 set_ticks(X_AXIS, x_range()*0.1, x_range()*0.02);
	 set_ticks(Y_AXIS, y_range()*0.1, y_range()*0.02);
      }
   }
}

void
coot::goograph::plot(int trace_id, int plot_type) {

   if (is_valid_trace(trace_id)) {
      if (plot_type == coot::goograph::PLOT_TYPE_BAR) {
	 plot_bar_graph(trace_id);
      }
      if (plot_type == coot::goograph::PLOT_TYPE_LINE) {
	 plot_line_graph(trace_id);
      }
      if (plot_type == coot::goograph::PLOT_TYPE_SMOOTHED_LINE) {
	 plot_smoothed_line_graph(trace_id);
      }
   } 
}

void
coot::goograph::plot_bar_graph(int trace_id) {

   if (is_valid_trace(trace_id)) { 
      const std::vector<std::pair<double, double> > &data = traces[trace_id].data;
      GooCanvasItem *root = goo_canvas_get_root_item(canvas);
      std::string colour = "#70e070";
      double mbw = median_bin_width(trace_id);

      for (unsigned int i=0; i<data.size(); i++) {
	 double width  = mbw * data_scale_x;
	 double height = -(data[i].second - extents_min_y) * data_scale_y;
	 lig_build::pos_t A(data[i].first-mbw*0.5, extents_min_y);
	 lig_build::pos_t wA = world_to_canvas(A);

	 std::cout << "  plot_bar_graph() " << A << " -> " << wA << " "
		   << width << " " << height << " " 
		   << " data_scale_x " << data_scale_x
		   << " data_scale_y " << data_scale_y
		   << std::endl;
	 
	 GooCanvasItem *rect =
	    goo_canvas_rect_new(root,
				wA.x, wA.y,
				width, height, 
				"line-width", 1.0,
				"fill_color", colour.c_str(),
				"stroke-color", "#333333",
				NULL);
      }
   }
}

void
coot::goograph::plot_line_graph(int trace_id) {

   if (is_valid_trace(trace_id)) { 
      const std::vector<std::pair<double, double> > &data = traces[trace_id].data;
      if (data.size()) { 
	 GooCanvasItem *root = goo_canvas_get_root_item(canvas);
	 int n_points = data.size();
	 GooCanvasPoints *points  = goo_canvas_points_new(n_points);
	 for (unsigned int i=0; i<data.size(); i++) {
	    lig_build::pos_t p(data[i].first, data[i].second);
	    lig_build::pos_t wp = world_to_canvas(p);
	    points->coords[2*i  ] = wp.x;
	    points->coords[2*i+1] = wp.y;
	 }

	 GooCanvasItem *line =
	    goo_canvas_polyline_new(root, 0, 0,
				    "line-width", 2.0,
				    "stroke-color", dark.c_str(),
				    "points", points,
				    NULL);
	 goo_canvas_points_unref(points); 
      }
   }
}

void
coot::goograph::plot_smoothed_line_graph(int trace_id) {

   if (is_valid_trace(trace_id)) { 
      const std::vector<std::pair<double, double> > &data = traces[trace_id].data;
      if (data.size()) { 
	 GooCanvasItem *root = goo_canvas_get_root_item(canvas);
	 int n_points = data.size();
	 std::string path_data;
	 for (unsigned int i=0; i<data.size(); i++) {
	    lig_build::pos_t p(data[i].first, data[i].second);
	    lig_build::pos_t wp = world_to_canvas(p);
	    std::string l;
	    if (i == 0)
	       l = "M";
	    if (i == 1)
	       l = "C";
 	    if (i == (n_points-2))
	       l = "S";
	    path_data += l;
	    path_data += " ";
	    path_data += coot::util::float_to_string(wp.x);
	    path_data += ",";
	    path_data += coot::util::float_to_string(wp.y);
	    path_data += " ";
	 }

	 GooCanvasItem *line =
	    goo_canvas_path_new(root, path_data.c_str(), 
				"line-width", 3.0,
				"stroke-color", dark.c_str(),
				NULL);
      }
   }
} 



double
coot::goograph::median_bin_width(int trace_id) const {

   double bw = 0.1;
   const std::vector<std::pair<double, double> > &data = traces[trace_id].data;
   std::vector<double> x(data.size(), 0);
   std::sort(x.begin(), x.end());
   if (data.size() > 1) {
      std::vector<double> bin_widths;
      for (unsigned int i=0; i<data.size(); i++)
	 x[i] = data[i].first;
      
      for (unsigned int i=0; i<x.size()-1; i++) { 
	 double this_bw = x[i+1] - x[i];
	 bin_widths.push_back(this_bw);
      }

      int idx = bin_widths.size() / 2;
      bw = bin_widths[idx];
   }
   return bw;
} 

void
coot::goograph::add_annotation_line(const lig_build::pos_t &pos_1,
				    const lig_build::pos_t &pos_2,
				    const std::string &colour,
				    double line_width,
				    bool dashed_flag,
				    bool start_arrow_in,
				    bool end_arrow_in) {
   
   lig_build::pos_t wA = world_to_canvas(pos_1);
   lig_build::pos_t wB = world_to_canvas(pos_2);
   gboolean start_arrow = start_arrow_in;
   gboolean   end_arrow =   end_arrow_in;
   GooCanvasItem *root = goo_canvas_get_root_item(canvas);
   GooCanvasItem *line =
      goo_canvas_polyline_new_line(root,
				   wA.x, wA.y,
				   wB.x, wB.y,
				   "line-width", line_width,
				   "stroke-color", colour.c_str(),
				   "start_arrow", start_arrow,
				   "end_arrow",   end_arrow,
				   NULL);

} 

void
coot::goograph::add_annotation_text(const std::string &text,
				    const lig_build::pos_t &pos_1,
				    const std::string &colour_in,
				    const std::string &font_in) {

   std::string font = font_in;
   std::string colour = colour_in;
   lig_build::pos_t wA = world_to_canvas(pos_1);
   if (font.empty())
      font = "Sans 9";
   if (colour.empty())
      colour = dark;
   GtkAnchorType anchor_type = GTK_ANCHOR_CENTER;
   GooCanvasItem *root = goo_canvas_get_root_item(canvas);
   GooCanvasItem *text_item =
      goo_canvas_text_new(root, text.c_str(),
			  wA.x, wA.y,
			  -1,
			  anchor_type,
			  "font", font.c_str(),
			  "fill_color", colour.c_str(),
			  NULL);

}


#endif


