/* src/goograph.cc
 * 
 * Copyright 2011, 2012 by The University of Oxford
 * Copyright 2012, 2014, 2015 by Medical Research Council
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

#include <iostream>
#include <algorithm>
#include "utils/coot-utils.hh"
#include "goograph.hh"

// returns the close button the first time it is called, otherwise returns null.
GtkWidget *
coot::goograph::show_dialog() {

   GtkWidget *close_button = NULL;

   draw_graph();

   if (! dialog)
      if (use_dialog_for_graph_flag)
         dialog = gtk_dialog_new();

   gtk_window_set_default_size(GTK_WINDOW(dialog), dialog_width, dialog_height);
   gtk_window_set_title (GTK_WINDOW(dialog), title_string.c_str());
   g_object_set_data(G_OBJECT(dialog), "goograph_dialog", dialog);
   g_object_set_data(G_OBJECT(dialog), "goograph", this);
   // GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;
   GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
   GtkWidget *vbox_inner = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
   GtkWidget *scrolled_window = gtk_scrolled_window_new (NULL, NULL);
   // gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(vbox_inner)); // gtk2
   gtk_container_add(GTK_CONTAINER(scrolled_window), GTK_WIDGET(vbox_inner)); // gtk3
   gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scrolled_window), TRUE, TRUE, 2);
   gtk_widget_show(scrolled_window);
   gtk_widget_show(vbox_inner);
   gtk_container_add(GTK_CONTAINER(vbox_inner), GTK_WIDGET(canvas));
   gtk_widget_show(GTK_WIDGET(canvas));
   close_button = gtk_dialog_add_button(GTK_DIALOG(dialog), "Close", 2);
   gtk_widget_show(close_button);
   g_signal_connect(G_OBJECT(close_button), "clicked",
                    G_CALLBACK(goograph_close_callback),
                    (gpointer) dialog);
   g_signal_connect(G_OBJECT(dialog), "configure_event",
                    G_CALLBACK(reshape),
                    (gpointer) dialog);

   if (use_dialog_for_graph_flag)
      gtk_widget_show(dialog);

   return close_button;
}

void
coot::goograph::init() {

   canvas = NULL;
   dialog = NULL;

   dark_mode = false;
   use_dialog_for_graph_flag = true;
   dialog_width_orig = 600;
   dialog_height_orig = 500;

   dialog_width = dialog_width_orig;
   dialog_height = dialog_height_orig;

   extents_x_are_set = false;
   extents_y_are_set = false;
   extents_min_x =  9999999990.0;
   extents_min_y =  9999999990.0;
   extents_max_x = -9999999990.0;
   extents_max_y = -9999999990.0;
   tick_major_x = 0.1;
   tick_minor_x = 0.05;
   tick_major_y = 0.1;
   tick_minor_y = 0.05;
   draw_x_axis_flag = true;
   draw_y_axis_flag = true;
   draw_x_ticks_flag = true;
   draw_y_ticks_flag = true;
   dark = "#111111";
   canvas_offset_x = 70.0; // how much is the Y axis displaced
   // from the left-hand edge of the canvas?
   canvas_offset_y = 370.0; // how much is the size of the
                            // canvas - and include an offset of
                            // the axis from the bottom edge.
                            // (the smaller the number the
                            // greater the displacement (without resizing).
   data_scale_x = 1.0;
   data_scale_y = 1.0;
}

void
coot::goograph::init_widgets() {

   dialog = NULL;
   if (! canvas) {
      canvas = GOO_CANVAS(goo_canvas_new());
   }
   // set_bounds(1000, 1000); // does this do any good?
}

void
coot::goograph::set_bounds(double right, double bottom) {

   double left = 0;
   double top = 0;
   // goo_canvas_set_bounds(canvas, left, top, right, bottom); // we can't "invert" these.
}

GtkWidget *
coot::goograph::get_canvas() const {
   return GTK_WIDGET(canvas);
}

void
coot::goograph::set_use_dialog_for_graph(bool flag) {
   use_dialog_for_graph_flag = flag;
}

void
coot::goograph::enable_dark_mode(bool flag) {
   dark_mode = flag;
}


// static
void
coot::goograph::goograph_close_callback(GtkWidget *button,
					GtkWidget *goograph_dialog) {
   gtk_widget_hide(goograph_dialog);
}

// static
gint
coot::goograph::reshape(GtkWidget *widget, GdkEventConfigure *event) {

   gint status = 0;
   coot::goograph *g_p = static_cast<coot::goograph *> (g_object_get_data(G_OBJECT(widget), "goograph"));

   if (g_p) {
      bool do_redraw = false;
      GtkAllocation allocation;
      gtk_widget_get_allocation(widget, &allocation);
      if (g_p->dialog_width != allocation.width)
 	 do_redraw = true;
      if (g_p->dialog_height != allocation.height)
 	 do_redraw = true;
      if (do_redraw) { 
 	 g_p->dialog_width  = allocation.width;
 	 g_p->dialog_height = allocation.height;
	 // g_p->set_data_scales();
 	 g_p->draw_graph();
	 // status = 1; // this slows things down dramatically!
      }
   }
   return status;
}



void
coot::goograph::draw_graph() {
   clear();
   draw_axes();
   draw_ticks();
   draw_title();
   draw_axis_label(X_AXIS);
   draw_axis_label(Y_AXIS);
   for (unsigned int i=0; i<traces.size(); i++) { 
      plot_trace(i);
   }
   draw_annotation_lines();
   draw_annotation_texts();
   draw_annotation_boxes();
   draw_contour_level_bar(); // for map histograms, handles callback
}

void
coot::goograph::set_draw_axis(int axis, bool draw_state) {

   if (axis == X_AXIS)
      draw_x_axis_flag = draw_state;
   if (axis == Y_AXIS)
      draw_y_axis_flag = draw_state;

}

void
coot::goograph::set_draw_ticks(int axis, bool draw_state) {

   if (axis == X_AXIS)
      draw_x_ticks_flag = draw_state;
   if (axis == Y_AXIS)
      draw_y_ticks_flag = draw_state;

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
//    GooCanvasLineDash *dash=goo_canvas_line_dash_new (2, 5.7, 0.0);

   lig_build::pos_t A_yaxis(extents_min_x, extents_min_y);
   lig_build::pos_t B_yaxis(extents_min_x, extents_min_y + y_range()*1.1);
   lig_build::pos_t wAy = world_to_canvas(A_yaxis);
   lig_build::pos_t wBy = world_to_canvas(B_yaxis);
   lig_build::pos_t A_xaxis(extents_min_x, extents_min_y);
   lig_build::pos_t B_xaxis(extents_min_x + x_range()*1.1, extents_min_y);
   if (false) {
      std::cout << "in draw_axes() X: extents " << extents_min_x << " "
		<< extents_max_x << " data_scale_x " << data_scale_x << std::endl;
      std::cout << "in draw_axes() Y: extents " << extents_min_y << " "
		<< extents_max_y << " data_scale_y " << data_scale_y << std::endl;
   }
   std::string grey = "#333333";
   
   lig_build::pos_t wAx = world_to_canvas(A_xaxis);
   lig_build::pos_t wBx = world_to_canvas(B_xaxis);
   gboolean start_arrow = 0;
   gboolean   end_arrow = 1;
   if (draw_y_axis_flag) {
      GooCanvasItem *item_1 = goo_canvas_polyline_new_line(root,
							   wAy.x, wAy.y,
							   wBy.x, wBy.y,
							   "line-width", 2.5,
							   "start_arrow", start_arrow,
							   "end_arrow",   end_arrow,
							   "stroke-color", grey.c_str(),
							   NULL);
      items.push_back(item_1);
   }
   if (draw_x_axis_flag) {
      GooCanvasItem *item_2 = goo_canvas_polyline_new_line(root,
							   wAx.x, wAx.y,
							   wBx.x, wBx.y,
							   "line-width", 2.5,
							   "start_arrow", start_arrow,
							   "end_arrow",   end_arrow,
							   "stroke-color", grey.c_str(),
							   NULL);
      items.push_back(item_2);
   }
}

void
coot::goograph::clear() {

   // If you are looking at this.... first check that your goograph is
   // still in scope! i.e. use a new pointer.
   // 
   if (canvas) {
      if (GOO_IS_CANVAS(canvas)) {
	 GooCanvasItem *root = goo_canvas_get_root_item(canvas);
	 if (root) { 
	    for (unsigned int i=0; i<items.size(); i++) {
	       gint child_index = goo_canvas_item_find_child(root, items[i]);
	       if (child_index != -1) {
		  goo_canvas_item_remove_child(root, child_index);
	       }
	    }
	 }
      } else {
	 std::cout << "ERROR:: canvas is not a GOOCANVAS " << canvas << std::endl;
      }
   }
   items.clear();
} 

void
coot::goograph::draw_ticks() {

   if (draw_y_ticks_flag) {
      draw_y_ticks(MINOR_TICK, tick_major_y, tick_minor_y, 0.5);
      draw_y_ticks(MAJOR_TICK, tick_major_y, tick_major_y, 1.0);
   }

   if (draw_x_ticks_flag) {
      draw_x_ticks(MINOR_TICK, tick_major_x, tick_minor_x, 0.5);
      draw_x_ticks(MAJOR_TICK, tick_major_x, tick_major_x, 1.0);
   }
}

void
coot::goograph::draw_y_ticks(int tick_type, double tick_major_step, double tick_step, double tick_length_multiplier) {

   draw_ticks_generic(Y_AXIS, tick_type, tick_major_step, tick_step, tick_length_multiplier);
}

void
coot::goograph::draw_x_ticks(int tick_type, double tick_major_step, double tick_step, double tick_length_multiplier) {
   draw_ticks_generic(X_AXIS, tick_type, tick_major_step, tick_step, tick_length_multiplier);
}

void
coot::goograph::draw_ticks_generic(int axis, int tick_type,
				   double tick_major_step,
				   double tick_step,
				   double tick_length_multiplier) {

   if (! canvas) {
      std::cout << "ERROR:: in draw_ticks_generic() canvas is null " << std::endl;
      return;
   }

   if (! GOO_IS_CANVAS(canvas)) {
      std::cout << "ERROR:: in draw_ticks_generic() canvas is not a goocanvas " << canvas << std::endl;
   }

   GooCanvasItem *root = goo_canvas_get_root_item(canvas);

   if (! root) {
      std::cout << "ERROR:: in draw_ticks_generic() root is null " << std::endl;
      return;
   }
      
   
   double extents_min = 0.0;
   double data_extents_min = 0.0;
   double extents_max = 0.0;

   bool tick_y_integers_only = false;

   if (axis == Y_AXIS) {
      tick_y_integers_only = true;
      for (unsigned int i=0; i<traces.size(); i++) { 
	 if (! traces[i].y_data_are_ints) {
	    tick_y_integers_only = false;
	    break;
	 } 
      }
   }

   // extents not set, no data yet.
   if (extents_min_x > extents_max_x) {
      std::cout << "no extents path - returning" << std::endl;
      return;
   } 

   if (axis == Y_AXIS) {
      data_extents_min = extents_min_y;
      extents_min = extents_min_y;
      extents_max = extents_max_y * 1.05;
   }
   
   if (axis == X_AXIS) {

      data_extents_min = extents_min_x;
      extents_min = extents_min_x;
      // extents_max = extents_max_x * 1.05;  // why did I have this?
      extents_max = extents_max_x; //  * 1.005;
      if (false) 
	 std::cout << "X_AXIS path"
		   << " extents_min_x"  << extents_min_x
		   << " extents_max_x " << extents_max_x
		   << " extents_min "   << extents_min
		   << " extents_max "   << extents_max
		   << std::endl;
   }

   // what is the tick that is a multiple of tick_major_step, that is
   // greater than or equal to extents_min?
   // 
   double tick_major_start = 0.0;
   // sanity check first
   if (tick_major_step > 0) { 
      if (tick_step > 0) { 
	 int n = int(extents_min/tick_major_step);
	 double diff = double(n)*tick_major_step - extents_min;
	 tick_major_start = extents_min + diff;
	 if (false) // debug
	    std::cout << "   extents_min: " << extents_min
		      << "  tick_major_step " << tick_major_step
		      << " n " << n
		      << " diff " << diff 
		      << " tick_major_start " << tick_major_start
		      << std::endl;
	 if (tick_major_start < extents_min) 
	    extents_min = tick_major_start + tick_major_step;
	 if (tick_type == MINOR_TICK) { 
	    extents_min -= tick_major_step;
	    while (extents_min < data_extents_min) { 
	       extents_min += tick_step;
	    }
	 }
      }
   }
   

   // 
   for (double tick_pos = extents_min; tick_pos <= extents_max; tick_pos+=tick_step) {

      bool do_tick = true;
      if (tick_type == MINOR_TICK) {
	 if (tick_y_integers_only)
	    do_tick = false;
	 // if the minor tick overlaps a major tick then don't display
	 // it (set do_tick to false).
	 for (double tick_pos_inner = tick_major_start;
	      tick_pos_inner <= extents_max;
	      tick_pos_inner+=tick_major_step) {
	    if (close_float_p(tick_pos_inner, tick_pos)) {
	       do_tick = false;
	       break;
	    }
	 }
      }
      if (tick_type == MAJOR_TICK) {
	 if (tick_y_integers_only) {
	    if (fabs(nearbyint(tick_pos)-tick_pos)>0.1)
	       do_tick = false;
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
      
         std::string grey = "#333333";
         if (dark_mode) grey = "#aaaaaa";

	 GooCanvasItem *tick =
	    goo_canvas_polyline_new_line(root,
					 wA.x, wA.y,
					 wB.x, wB.y,
					 "line-width", 1.0,
					 "stroke-color", grey.c_str(),
					 NULL);
	 items.push_back(tick);

	 GooCanvasAnchorType anchor_type = GOO_CANVAS_ANCHOR_CENTER;
	 if (axis == X_AXIS)
	    anchor_type = GOO_CANVAS_ANCHOR_SW;

	 // if it's a MAJOR_TICK, then we want a text label too
	 // 
	 if (tick_type == MAJOR_TICK) {
	    int n_dec_pl = 0;
	    double range = y_range();
	    if (axis == X_AXIS)
	       range = x_range();

	    if (range < 5)
	       n_dec_pl = 1;
	    if (range < 2)
	       n_dec_pl = 2;

	    std::string txt =
	       coot::util::float_to_unspaced_string_using_dec_pl(tick_pos, n_dec_pl);
	    if (tick_y_integers_only)
	       txt = coot::util::int_to_string(int(nearbyint(tick_pos)));
		  
	    GooCanvasItem *text = goo_canvas_text_new(root, txt.c_str(),
						      wB.x - tick_label_x_off,
						      wB.y - tick_label_y_off,
						      -1,
						      anchor_type,
						      "font", "Sans 7",
						      "fill_color", grey.c_str(),
						      NULL);
	    items.push_back(text);
	 }
      }
   }
}

void
coot::goograph::set_extents(int axis, double min, double max) {

   if (! extents_x_are_set) {
      if (axis == X_AXIS) {
	 extents_min_x = min;
	 extents_max_x = max;
	 double x_major_tick = calc_tick(x_range());
	 set_ticks(X_AXIS, x_major_tick, x_major_tick*0.2);
	 extents_x_are_set = true;
      }
   }
   if (! extents_y_are_set) {
      if (axis == Y_AXIS) {
	 extents_min_y = min;
	 extents_max_y = max;
	 double y_major_tick = calc_tick(y_range());
	 set_ticks(Y_AXIS, y_major_tick, y_major_tick*0.2);
	 extents_y_are_set = true;
      }
   }

   set_data_scales(axis);
}

void
coot::goograph::set_data_scales(int axis) {

   double dsx = data_scale_x;
   double dsy = data_scale_y;

   if (axis == X_AXIS) { 
      double delta_x = extents_max_x - extents_min_x;
      if (delta_x > 0) {
	 // was this to cope with resizing graphs?
	 data_scale_x = double(dialog_width)/double(dialog_width_orig) * 400.0/delta_x;
	 // it's more important to get small graphs to work properly
	 data_scale_x = 0.6 * double(dialog_width)/delta_x;
      } else {
	 data_scale_x = 0.1; // who knows what a sane fall-back is.
      }
      
      if (false)
	 std::cout << "DEBUG:: ------------------- set_data_scales() "
		   << "delta_x " << delta_x << " "
		   << "data_scale_x was " << dsx << " now  " << data_scale_x << std::endl;
   }
   
   if (axis == Y_AXIS) { 
      double delta_y = extents_max_y - extents_min_y;
      if (delta_y > 0.0) {
	 data_scale_y = double(dialog_height)/double(dialog_height_orig) * 300.0/delta_y;
	 data_scale_y = 0.6 * double(dialog_height)/delta_y;
      } else {
	 data_scale_y = 1;
      }

      if (false)
	 std::cout << "DEBUG:: ------------------- set_data_scales() "
		   << "delta_y " << delta_y << " "
		   << "data_scale_y was " << dsy << " now  " << data_scale_y << std::endl;
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
   
   if (axis == X_AXIS) {
      x_axis_label = label;
   }
   if (axis == Y_AXIS) {
      y_axis_label = label;
   }
} 

void
coot::goograph::draw_axis_label(int axis) {

   GooCanvasItem *root = goo_canvas_get_root_item(canvas);
   lig_build::pos_t A;
   std::string label;
   bool do_it = true;
   if (axis == X_AXIS) {
      A = lig_build::pos_t(extents_min_x+x_range()*0.8, extents_min_y -y_range()*0.12);
      label = x_axis_label;
      if (draw_x_axis_flag == false)
	 do_it = false;
   }
   if (axis == Y_AXIS) {
      A = lig_build::pos_t(extents_min_x-x_range()*0.14, extents_min_y + y_range()*1.15);
      label = y_axis_label;
      if (draw_y_axis_flag == false)
	 do_it = false;
   }
   if (do_it) {
      std::string grey = "#333333";
      if (dark_mode) grey = "#aaaaaa";
      // std::cout << "draw_axis_label() " << label << std::endl;
      lig_build::pos_t wA = world_to_canvas(A);
      GooCanvasAnchorType anchor_type = GOO_CANVAS_ANCHOR_NORTH_WEST;
      GooCanvasItem *text =
	 goo_canvas_text_new(root, label.c_str(),
			     wA.x, wA.y,
			     -1,
			     anchor_type,
			     "font", "Sans 9",
			     "fill_color", grey.c_str(),
			     NULL);
      items.push_back(text);
   }
}

void
coot::goograph::set_plot_title(const std::string &title) {
   title_string = title;
}

void
coot::goograph::draw_title() {

   if (! title_string.empty()) {
      lig_build::pos_t A(extents_min_x + 0.5 * x_range(),
			 extents_min_y + 1.15 * y_range());
      lig_build::pos_t wA = world_to_canvas(A);
      GooCanvasItem *root = goo_canvas_get_root_item(canvas);
      GooCanvasAnchorType anchor_type = GOO_CANVAS_ANCHOR_CENTER;
      std::string grey = "#333333";
      if (dark_mode) grey = "#aaaaaa";
      GooCanvasItem *text =
	 goo_canvas_text_new(root, title_string.c_str(),
			     wA.x, wA.y,
			     -1,
			     anchor_type,
			     "font", "Sans 11",
			     "fill_color", grey.c_str(),
			     NULL);
      items.push_back(text);
   }
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
      traces[trace_id].set_data(sorted_data);

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
	 if (false) {
	    std::cout << "   in set_data() setting X exents "
		      << min_x << " " << max_x << std::endl;
	    std::cout << "   in set_data() setting Y exents "
		      << min_y << " " << max_y << std::endl;
	 }
	 set_extents(X_AXIS, min_x, max_x);
	 set_extents(Y_AXIS, min_y, max_y);

	 // std::cout << "in set_data: x_range() is " << x_range() << std::endl;
	 double x_major_tick = calc_tick(x_range());
	 double y_major_tick = calc_tick(y_range());
	 set_ticks(X_AXIS, x_major_tick, x_major_tick*0.2);
	 set_ticks(Y_AXIS, y_major_tick, y_major_tick*0.2);
      }
   }
}

void
coot::graph_trace_info_t::set_data(const std::vector<std::pair<double, double> > &data_in) {

   data = data_in;
   
   // check y data for ints
   y_data_are_ints = true;
   for (unsigned int i=0; i<data.size(); i++) { 
      if (fabs(data[i].second-nearbyint(data[i].second)) > 0.000001) {
	 y_data_are_ints = false;
	 break;
      }
   }
}


// return closes multiple of 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100...
double
coot::goograph::calc_tick(double range) const { 

   double r = 1;
   if (range > 0) {
      double m = 1.0/log(10);
      double r_1 = log(range*0.1)*m;
      double rr_1 = fabs(r_1 - nearbyint(r_1));

      double r_2 = log(range/2)*m;
      double rr_2 = fabs(r_2 - nearbyint(r_2));
      
      double r_5 = log(range/5)*m;
      double rr_5 = fabs(r_5 - nearbyint(r_5));
      
      double nbi = nearbyint(r_1);
      r = pow(10, nbi);
      if (rr_2 < rr_1) { 
 	 nbi = nearbyint(r_2);
 	 r = 0.2*pow(10, nbi);
 	 if (rr_5 < rr_2) {
 	    nbi = nearbyint(r_5);
 	    r = 0.5*pow(10, nbi);
 	 }
      }
      if (0) 
	 std::cout << "  calc_tick() range " << range << "    "
		   << " r_1 " << r_1
		   << " rr_1 " << rr_1
		   << " r_2 " << r_2
		   << " rr_2 " << rr_2
		   << " r_5 " << r_5
		   << " rr_5 " << rr_5
		   << " nbi: " << nbi
		   << " r: " << r << std::endl;
   } 
   return r;
}

void
coot::goograph::set_trace_type(int trace_id, int plot_type, bool dashed) {
   
   if (is_valid_trace(trace_id)) {
      if (plot_type == coot::graph_trace_info_t::PLOT_TYPE_BAR)
	 traces[trace_id].plot_type = plot_type;
      if (plot_type == coot::graph_trace_info_t::PLOT_TYPE_LINE)
	 traces[trace_id].plot_type = plot_type;
      if (plot_type == coot::graph_trace_info_t::PLOT_TYPE_SMOOTHED_LINE)
	 traces[trace_id].plot_type = plot_type;
      if (plot_type == coot::graph_trace_info_t::PLOT_TYPE_SCATTER)
	 traces[trace_id].plot_type = plot_type;
      traces[trace_id].dashed = dashed;
   } 
}

void
coot::goograph::set_trace_colour(int trace_id, const std::string &colour) {

   if (is_valid_trace(trace_id))
      traces[trace_id].colour = colour;
} 


void
coot::goograph::plot_trace(int trace_id) {

   if (is_valid_trace(trace_id)) {
      if (traces[trace_id].plot_type == graph_trace_info_t::PLOT_TYPE_BAR) {
	 plot_bar_graph(trace_id);
      }
      if (traces[trace_id].plot_type == graph_trace_info_t::PLOT_TYPE_LINE) {
	 plot_line_graph(trace_id);
      }
      if (traces[trace_id].plot_type == graph_trace_info_t::PLOT_TYPE_SMOOTHED_LINE) {
	 plot_smoothed_line_graph(trace_id);
      }
      if (traces[trace_id].plot_type == graph_trace_info_t::PLOT_TYPE_SCATTER) {
	 plot_scatter_plot(trace_id);
      }
   } 
}

void
coot::goograph::plot_scatter_plot(int trace_id) {

   GooCanvasItem *root = goo_canvas_get_root_item(canvas);
   if (is_valid_trace(trace_id)) { 
      const std::vector<std::pair<double, double> > &data = traces[trace_id].get_data();
      std::string colour = traces[trace_id].colour;
      if (colour.empty())
	 colour = "#70e070";
      double mbw = median_bin_width(trace_id);
      double line_width = 1.0;
      double radius = 3; 

      for (unsigned int i=0; i<data.size(); i++) {
	 
	 // double width  = mbw * data_scale_x;
	 // double height = -(data[i].second - extents_min_y) * data_scale_y;
	 // lig_build::pos_t A(data[i].first-mbw*0.5, extents_min_y);
	 // 	    goo_canvas_ci(root,
	 // 				wA.x, wA.y,
	 // 				width, height, 
	 // 				"line-width", 1.0,
	 // 				"fill_color", colour.c_str(),
	 // 				"stroke-color", "#333333",
	 // 			  NULL);

	 
	 lig_build::pos_t p(data[i].first, data[i].second);
	 lig_build::pos_t wA = world_to_canvas(p);

	 GooCanvasItem *ring =
	 goo_canvas_ellipse_new(root, wA.x, wA.y,
				radius, radius,
				"line_width", line_width,
				// "fill-color-rgba", 0xffbb3350,
				NULL);
	 
	 items.push_back(ring);
    	 
      }
   }
}


void
coot::goograph::plot_bar_graph(int trace_id) {

   GooCanvasItem *root = goo_canvas_get_root_item(canvas);
   if (is_valid_trace(trace_id)) { 
      const std::vector<std::pair<double, double> > &data = traces[trace_id].get_data();
      std::string colour = traces[trace_id].colour;
      if (colour.empty())
	 colour = "#70e070";
      double mbw = median_bin_width(trace_id);

      if (false)
         std::cout << "debug  in plot_bar_graph() mbw is "
		   << mbw << " data_scale_x " << data_scale_x << std::endl;

      for (unsigned int i=0; i<data.size(); i++) {
	 double width  = mbw * data_scale_x;
	 double height = -(data[i].second - extents_min_y) * data_scale_y;
	 lig_build::pos_t A(data[i].first-mbw*0.5, extents_min_y);

	 if (data[i].second > 0) {
	    lig_build::pos_t wA = world_to_canvas(A);
	    if (false)
	       std::cout << "plot_bar_graph() " << i << " "
			 << A << " -> " << wA << " "
			 << width << " " << height << " "
			 << " " << data[i].first << " " << data[i].second
			 << " data-scales: " << data_scale_x << " " << data_scale_y
			 << std::endl;

	    GooCanvasItem *rect =
	       goo_canvas_rect_new(root,
				   wA.x, wA.y,
				   width, height, 
				   "line-width", 1.0,
				   "fill_color", colour.c_str(),
				   "stroke-color", "#333333",
				   NULL);
	    items.push_back(rect);
	 }
      }
   }
}


void
coot::goograph::plot_line_graph(int trace_id) {

   if (is_valid_trace(trace_id)) { 
      const std::vector<std::pair<double, double> > &data = traces[trace_id].get_data();
      std::string colour = traces[trace_id].colour;
      bool dashed = traces[trace_id].dashed;
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
	 if (colour.empty())
	    colour = dark;

	 GooCanvasItem *line =
	    goo_canvas_polyline_new(root, 0, 0,
				    "line-width", 2.0,
				    "stroke-color", colour.c_str(),
				    "points", points,
				    NULL);
	 goo_canvas_points_unref(points); 
	 items.push_back(line);
      }
   }
}

void
coot::goograph::plot_smoothed_line_graph(int trace_id) {

   if (is_valid_trace(trace_id)) {
      const std::vector<std::pair<double, double> > &data = traces[trace_id].get_data();
      std::string colour = traces[trace_id].colour;
      bool dashed = traces[trace_id].dashed;
      if (data.size()) { 
	 GooCanvasItem *root = goo_canvas_get_root_item(canvas);
	 unsigned int n_points = data.size();
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

	 if (colour.empty())
	    colour = dark;

	 GooCanvasLineDash *dash;
	 if (dashed)
	    dash = goo_canvas_line_dash_new (2, 5.7, 3.0);
	 else
	    dash = goo_canvas_line_dash_new (2, 4.0, 0.0);
	 
	 GooCanvasItem *line =
	    goo_canvas_path_new(root, path_data.c_str(), 
				"line-width", 3.0,
				"line-dash", dash,
				"stroke-color", colour.c_str(),
				NULL);
	 items.push_back(line);
      }
   }
}



double
coot::goograph::median_bin_width(int trace_id) const {

   double bw = 0.1;
   const std::vector<std::pair<double, double> > &data = traces[trace_id].get_data();
   std::vector<double> x(data.size(), 0);
   std::sort(x.begin(), x.end());
   if (data.size() > 1) {
      std::vector<double> bin_widths;
      for (unsigned int i=0; i<data.size(); i++)
	 x[i] = data[i].first;
      
      for (unsigned int i=0; i<x.size()-1; i++) { 
	 double this_bw = x[i+1] - x[i];
	 if (this_bw > 0)
	    bin_widths.push_back(this_bw);
      }

      if (bin_widths.size() > 0) {
	 int idx = bin_widths.size() / 2;
	 bw = bin_widths[idx];
      }
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

   annotation_line_t line(pos_1, pos_2, colour, line_width, dashed_flag, start_arrow_in, end_arrow_in);
   annotation_lines.push_back(line);

} 

void
coot::goograph::draw_annotation_lines() {

   GooCanvasItem *root = goo_canvas_get_root_item(canvas);

   for (unsigned int i=0; i<annotation_lines.size(); i++) { 
      lig_build::pos_t wA = world_to_canvas(annotation_lines[i].pos_1);
      lig_build::pos_t wB = world_to_canvas(annotation_lines[i].pos_2);
      gboolean start_arrow = annotation_lines[i].start_arrow;
      gboolean   end_arrow = annotation_lines[i].end_arrow;
      double   line_width  = annotation_lines[i].line_width;
      std::string colour   = annotation_lines[i].colour;
      GooCanvasItem *line =
	 goo_canvas_polyline_new_line(root,
				      wA.x, wA.y,
				      wB.x, wB.y,
				      "line-width", line_width,
				      "stroke-color", colour.c_str(),
				      "start_arrow", start_arrow,
				      "end_arrow",   end_arrow,
				      NULL);
      items.push_back(line);
   }

} 

void
coot::goograph::add_annotation_text(const std::string &text,
				    const lig_build::pos_t &pos_1,
				    const std::string &colour_in,
				    const std::string &font_in) {
   annotation_text_t at(text, pos_1, colour_in, font_in);
   annotation_texts.push_back(at);

}

void
coot::goograph::add_annotation_box(const lig_build::pos_t &pos_top_left,
				   const lig_build::pos_t &pos_bottom_right,
				   const std::string &outline_colour,
				   double line_width,
				   const std::string &fill_colour) {

   annotation_box_t b(pos_top_left, pos_bottom_right, outline_colour, line_width, fill_colour);
   annotation_boxes.push_back(b);

}

void
coot::goograph::add_contour_level_box(float current_contour_level,
				      const std::string &outline_colour,
				      double line_width,
				      const std::string &fill_colour,
				      int imol,
				      float rmsd,
				      void (*func)(int, float)) {

   lig_build::pos_t dum;

   annotated_box_info_t b(dum, dum, outline_colour, line_width, fill_colour,
			  imol, rmsd, current_contour_level, func);
   contour_level_bar = b;

}



void
coot::goograph::draw_annotation_texts() { 

   for (unsigned int i=0; i<annotation_texts.size(); i++) { 

      std::string text = annotation_texts[i].text;
      std::string font = annotation_texts[i].font;
      std::string colour = annotation_texts[i].colour;
      lig_build::pos_t wA = world_to_canvas(annotation_texts[i].pos);
      if (font.empty())
	 font = "Sans 9";
      if (colour.empty())
	 colour = dark;
      GooCanvasAnchorType anchor_type = GOO_CANVAS_ANCHOR_CENTER;
      GooCanvasItem *root = goo_canvas_get_root_item(canvas);
      GooCanvasItem *text_item =
	 goo_canvas_text_new(root, text.c_str(),
			     wA.x, wA.y,
			     -1,
			     anchor_type,
			     "font", font.c_str(),
			     "fill_color", colour.c_str(),
			     NULL);
      items.push_back(text_item);
   }
}

bool
coot::goograph::on_goograph_active_button_box_press_event(GooCanvasItem  *item,
					  GooCanvasItem  *target_item,
					  GdkEventButton *event,
					  gpointer        user_data) {

   annotation_box_info_t *bi =
      static_cast<annotation_box_info_t *> (g_object_get_data(G_OBJECT (item), "box-info"));
   if (bi) {
      bi->x_mouse = event->x;
      bi->y_mouse = event->y;
   } else {
      std::cout << "error:: null bi in on_goograph_active_button_box_press_event" << std::endl;
   }
   return TRUE;
}

// static
bool
coot::goograph::on_goograph_active_button_box_release_event(GooCanvasItem  *item,
							    GooCanvasItem  *target_item,
							    GdkEventButton *event,
							    gpointer        user_data) {

   // Let's extract where we started:
   annotation_box_info_t *bi =
      static_cast<annotation_box_info_t *> (g_object_get_data(G_OBJECT (item), "box-info"));
   if (bi) {
      double current_x = event->x;
      double delta_x = current_x - bi->x_mouse;

      double delta_contour = delta_x * 7.0*bi->rmsd/300.0;
      double new_contour_level = bi->current_contour_level + delta_contour;

      if (bi->func) {
	 goograph *self = static_cast<goograph *> (user_data);
	 bi->func(bi->imol, new_contour_level);
	 self->contour_level_bar.abi.current_contour_level += delta_contour;
	 delete bi; // draw_graph() makes a new bi
	 // redraw self
	 self->draw_graph();
      } else {
	 std::cout << "null func" << std::endl;
      }
   }

   return TRUE;
}

void
coot::goograph::draw_annotation_boxes() {

   for (unsigned int i=0; i<annotation_boxes.size(); i++) {

      std::string fill_col = annotation_boxes[i].fill_colour;
      std::string stroke_col = annotation_boxes[i].outline_colour;

      lig_build::pos_t wA = world_to_canvas(annotation_boxes[i].top_left);
      lig_build::pos_t wB = world_to_canvas(annotation_boxes[i].bottom_right);
      double width  = wB.x - wA.x;
      double height = wB.y - wA.y;
      double lw = annotation_boxes[i].line_width;

      GooCanvasItem *root = goo_canvas_get_root_item(canvas);
      GooCanvasItem *rect =
	 goo_canvas_rect_new(root,
			     wA.x, wA.y,
			     width, height,
			     "line-width", lw,
			     "fill_color", fill_col.c_str(),
			     "stroke-color", stroke_col.c_str(),
			     NULL);
      items.push_back(rect);
   }
}

// special function for the map density histogram

void
coot::goograph::draw_contour_level_bar() {

   if (contour_level_bar.abi.is_set) {
      std::string fill_col = contour_level_bar.fill_colour;
      std::string stroke_col = contour_level_bar.outline_colour;

      double y_max_graph = extents_max_y * 1.25; // the graph/canvas has size above extents_max_y

      double contour_level_bar_width = contour_level_bar.abi.rmsd * 0.14;

      double ccl = contour_level_bar.abi.current_contour_level;

      lig_build::pos_t tl(ccl-0.5*contour_level_bar_width, y_max_graph);
      lig_build::pos_t br(ccl+0.5*contour_level_bar_width, -y_max_graph*0.05);

      lig_build::pos_t wA = world_to_canvas(tl);
      lig_build::pos_t wB = world_to_canvas(br);
      double width  = wB.x - wA.x;
      double height = wB.y - wA.y;
      double lw = contour_level_bar.line_width;

      GooCanvasItem *root = goo_canvas_get_root_item(canvas);
      GooCanvasItem *rect =
	 goo_canvas_rect_new(root,
			     wA.x, wA.y,
			     width, height,
			     "line-width", lw,
			     "fill_color", fill_col.c_str(),
			     "stroke-color", stroke_col.c_str(),
			     NULL);
      items.push_back(rect);

      annotation_box_info_t *bi = new annotation_box_info_t(contour_level_bar.abi);

      g_object_set_data(G_OBJECT(rect), "box-info", bi);

      g_signal_connect(rect, "button_press_event",
		       G_CALLBACK(on_goograph_active_button_box_press_event), NULL);
      g_signal_connect(rect, "button_release_event",
		       G_CALLBACK(on_goograph_active_button_box_release_event), this);

   }
}

void
coot::goograph::clear_traces_and_annotations() {

   traces.clear();
   annotation_texts.clear();
   annotation_lines.clear();
   extents_min_x =  9999999990.0;
   extents_min_y =  9999999990.0;
   extents_max_x = -9999999990.0;
   extents_max_y = -9999999990.0;

}


