// GTK4/Cairo sketch replacement for coot::goograph.

#include "graph2d.hh"
#include <algorithm>
#include <cmath>
#include <limits>

namespace coot {

graph2d::graph2d(int w, int h) : width_req(w), height_req(h) {
   drawing_area = gtk_drawing_area_new();
   gtk_widget_set_size_request(drawing_area, w, h);
   gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area),
                                  &graph2d::draw_cb, this, nullptr);
}

int graph2d::trace_new() {
   traces.emplace_back();
   return int(traces.size()) - 1;
}

void graph2d::set_trace_data(int id,
                             const std::vector<std::pair<double,double>> &d) {
   if (id < 0 || id >= int(traces.size())) return;
   traces[id].data = d;
   if (!extents_x_set || !extents_y_set) update_auto_extents();
   queue_redraw();
}

void graph2d::set_trace_type(int id, graph_trace_t::plot_type_t t) {
   if (id < 0 || id >= int(traces.size())) return;
   traces[id].plot_type = t;
   queue_redraw();
}

void graph2d::set_trace_colour(int id, double r, double g, double b) {
   if (id < 0 || id >= int(traces.size())) return;
   traces[id].r = r; traces[id].g = g; traces[id].b = b;
   queue_redraw();
}

void graph2d::set_extents(axis_t ax, double lo, double hi) {
   if (ax == X_AXIS) { x_lo = lo; x_hi = hi; extents_x_set = true; }
   else              { y_lo = lo; y_hi = hi; extents_y_set = true; }
   queue_redraw();
}

void graph2d::set_axis_label(axis_t ax, const std::string &s) {
   if (ax == X_AXIS) x_label = s; else y_label = s;
   queue_redraw();
}

void graph2d::set_draw_axis (axis_t ax, bool st) {
   if (ax == X_AXIS) draw_x_axis = st; else draw_y_axis = st;
   queue_redraw();
}

void graph2d::set_draw_ticks(axis_t ax, bool st) {
   if (ax == X_AXIS) draw_x_tick = st; else draw_y_tick = st;
   queue_redraw();
}

void graph2d::queue_redraw() {
   if (drawing_area) gtk_widget_queue_draw(drawing_area);
}

GtkWidget *graph2d::show_dialog(GtkWindow *parent) {
   GtkWidget *win = gtk_window_new();
   if (parent) gtk_window_set_transient_for(GTK_WINDOW(win), parent);
   gtk_window_set_title(GTK_WINDOW(win),
                        plot_title.empty() ? "Graph" : plot_title.c_str());
   gtk_window_set_default_size(GTK_WINDOW(win), width_req, height_req);
   gtk_window_set_child(GTK_WINDOW(win), drawing_area);
   gtk_window_present(GTK_WINDOW(win));
   return win;
}

void graph2d::update_auto_extents() {
   if (traces.empty()) return;
   double xmin =  std::numeric_limits<double>::infinity();
   double xmax = -std::numeric_limits<double>::infinity();
   double ymin =  std::numeric_limits<double>::infinity();
   double ymax = -std::numeric_limits<double>::infinity();
   bool any = false;
   for (const auto &t : traces) {
      for (const auto &p : t.data) {
         xmin = std::min(xmin, p.first);  xmax = std::max(xmax, p.first);
         ymin = std::min(ymin, p.second); ymax = std::max(ymax, p.second);
         any = true;
      }
   }
   if (!any) return;
   if (xmin == xmax) { xmin -= 0.5; xmax += 0.5; }
   if (ymin == ymax) { ymin -= 0.5; ymax += 0.5; }
   if (!extents_x_set) { x_lo = xmin; x_hi = xmax; }
   if (!extents_y_set) { y_lo = ymin; y_hi = ymax; }
}

double graph2d::nice_tick(double range) {
   if (range <= 0.0) return 1.0;
   double exp10 = std::pow(10.0, std::floor(std::log10(range)));
   double f = range / exp10;
   double n = (f < 1.5) ? 0.2 : (f < 3.0) ? 0.5 : (f < 7.0) ? 1.0 : 2.0;
   return n * exp10;
}

graph2d::plot_box_t graph2d::plot_box(int w, int h) const {
   return { left_margin, top_margin,
            std::max(1.0, w - left_margin - right_margin),
            std::max(1.0, h - top_margin  - bottom_margin) };
}

void graph2d::draw_cb(GtkDrawingArea *, cairo_t *cr,
                      int w, int h, gpointer user_data) {
   static_cast<graph2d*>(user_data)->render(cr, w, h);
}

void graph2d::render(cairo_t *cr, int w, int h) {
   // background
   cairo_set_source_rgb(cr, 1, 1, 1);
   cairo_paint(cr);

   plot_box_t pb = plot_box(w, h);

   if (draw_x_axis || draw_y_axis) draw_axes (cr, pb);
   if (draw_x_tick || draw_y_tick) draw_ticks(cr, pb);

   for (const auto &t : traces) draw_trace(cr, pb, t);

   draw_labels(cr, w, h, pb);
   draw_title (cr, w);
}

void graph2d::draw_axes(cairo_t *cr, const plot_box_t &pb) {
   cairo_set_source_rgb(cr, 0, 0, 0);
   cairo_set_line_width(cr, 1.0);
   if (draw_y_axis) {
      cairo_move_to(cr, pb.x, pb.y);
      cairo_line_to(cr, pb.x, pb.y + pb.h);
   }
   if (draw_x_axis) {
      cairo_move_to(cr, pb.x,          pb.y + pb.h);
      cairo_line_to(cr, pb.x + pb.w,   pb.y + pb.h);
   }
   cairo_stroke(cr);
}

void graph2d::draw_ticks(cairo_t *cr, const plot_box_t &pb) {
   cairo_set_source_rgb(cr, 0, 0, 0);
   cairo_set_line_width(cr, 1.0);
   cairo_set_font_size(cr, 10.0);

   const double tl = 4.0; // tick length

   if (draw_x_tick) {
      double step = nice_tick(x_hi - x_lo);
      double start = std::ceil(x_lo / step) * step;
      for (double v = start; v <= x_hi + 1e-9; v += step) {
         double px = w2p_x(v, pb);
         cairo_move_to(cr, px, pb.y + pb.h);
         cairo_line_to(cr, px, pb.y + pb.h + tl);
         char buf[32];
         std::snprintf(buf, sizeof buf, "%g", v);
         cairo_text_extents_t te; cairo_text_extents(cr, buf, &te);
         cairo_move_to(cr, px - te.width * 0.5, pb.y + pb.h + tl + te.height + 2);
         cairo_show_text(cr, buf);
      }
   }
   if (draw_y_tick) {
      double step = nice_tick(y_hi - y_lo);
      double start = std::ceil(y_lo / step) * step;
      for (double v = start; v <= y_hi + 1e-9; v += step) {
         double py = w2p_y(v, pb);
         cairo_move_to(cr, pb.x - tl, py);
         cairo_line_to(cr, pb.x,      py);
         char buf[32];
         std::snprintf(buf, sizeof buf, "%g", v);
         cairo_text_extents_t te; cairo_text_extents(cr, buf, &te);
         cairo_move_to(cr, pb.x - tl - te.width - 3, py + te.height * 0.4);
         cairo_show_text(cr, buf);
      }
   }
   cairo_stroke(cr);
}

void graph2d::draw_title(cairo_t *cr, int w) {
   if (plot_title.empty()) return;
   cairo_set_source_rgb(cr, 0, 0, 0);
   cairo_set_font_size(cr, 13.0);
   cairo_text_extents_t te;
   cairo_text_extents(cr, plot_title.c_str(), &te);
   cairo_move_to(cr, (w - te.width) * 0.5, top_margin * 0.7);
   cairo_show_text(cr, plot_title.c_str());
}

void graph2d::draw_labels(cairo_t *cr, int w, int h, const plot_box_t &pb) {
   cairo_set_source_rgb(cr, 0, 0, 0);
   cairo_set_font_size(cr, 11.0);
   if (!x_label.empty()) {
      cairo_text_extents_t te; cairo_text_extents(cr, x_label.c_str(), &te);
      cairo_move_to(cr, pb.x + (pb.w - te.width) * 0.5, h - 6.0);
      cairo_show_text(cr, x_label.c_str());
   }
   if (!y_label.empty()) {
      cairo_text_extents_t te; cairo_text_extents(cr, y_label.c_str(), &te);
      cairo_save(cr);
      cairo_move_to(cr, 12.0, pb.y + (pb.h + te.width) * 0.5);
      cairo_rotate(cr, -M_PI * 0.5);
      cairo_show_text(cr, y_label.c_str());
      cairo_restore(cr);
   }
}

void graph2d::draw_trace(cairo_t *cr, const plot_box_t &pb,
                         const graph_trace_t &t) {
   if (t.data.empty()) return;
   cairo_set_source_rgb(cr, t.r, t.g, t.b);
   cairo_set_line_width(cr, t.line_width);
   if (t.dashed) {
      double d[] = { 4.0, 3.0 };
      cairo_set_dash(cr, d, 2, 0);
   } else {
      cairo_set_dash(cr, nullptr, 0, 0);
   }

   switch (t.plot_type) {
   case graph_trace_t::LINE:
   case graph_trace_t::SMOOTHED_LINE: // TODO spline
      for (size_t i = 0; i < t.data.size(); ++i) {
         double px = w2p_x(t.data[i].first,  pb);
         double py = w2p_y(t.data[i].second, pb);
         if (i == 0) cairo_move_to(cr, px, py);
         else        cairo_line_to(cr, px, py);
      }
      cairo_stroke(cr);
      break;
   case graph_trace_t::BAR: {
      double half = 0.4 * pb.w / std::max<size_t>(1, t.data.size());
      for (const auto &p : t.data) {
         double px = w2p_x(p.first, pb);
         double py = w2p_y(p.second, pb);
         double y0 = w2p_y(std::max(0.0, y_lo), pb);
         cairo_rectangle(cr, px - half, std::min(py, y0),
                         2 * half, std::fabs(py - y0));
      }
      cairo_fill(cr);
      break;
   }
   case graph_trace_t::SCATTER:
      for (const auto &p : t.data) {
         double px = w2p_x(p.first, pb);
         double py = w2p_y(p.second, pb);
         cairo_arc(cr, px, py, 2.5, 0, 2 * M_PI);
         cairo_fill(cr);
      }
      break;
   }
}

}
