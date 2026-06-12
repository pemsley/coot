#ifndef COOT_GRAPH2D_HH
#define COOT_GRAPH2D_HH

// GTK4/Cairo replacement sketch for coot::goograph.
// Goocanvas has no GTK4 port, so this abandons the retained-mode
// GooCanvasItem model and redraws from the data each frame.

#include <string>
#include <utility>
#include <vector>
#include <gtk/gtk.h>

namespace coot {

   class graph_trace_t {
   public:
      enum plot_type_t { LINE, BAR, SMOOTHED_LINE, SCATTER };
      std::vector<std::pair<double,double>> data;
      plot_type_t plot_type = LINE;
      double r = 0.2, g = 0.4, b = 0.8;   // RGB in [0,1]
      double line_width = 1.5;
      bool   dashed = false;
   };

   class graph2d {
   public:
      enum axis_t { X_AXIS, Y_AXIS };

      graph2d(int width = 500, int height = 360);

      // Data
      int  trace_new();
      void set_trace_data(int id, const std::vector<std::pair<double,double>> &d);
      void set_trace_type (int id, graph_trace_t::plot_type_t t);
      void set_trace_colour(int id, double r, double g, double b);

      // Layout
      void set_extents(axis_t ax, double lo, double hi);
      void set_axis_label(axis_t ax, const std::string &s);
      void set_plot_title(const std::string &s)  { plot_title = s; queue_redraw(); }
      void set_draw_axis (axis_t ax, bool state);
      void set_draw_ticks(axis_t ax, bool state);

      // Widget access / presentation
      GtkWidget *get_widget()   const { return drawing_area; }     // embed
      GtkWidget *show_dialog(GtkWindow *parent = nullptr);         // standalone

      // Force a repaint (public so callers can nudge after batch edits)
      void queue_redraw();

   private:
      // Data
      std::vector<graph_trace_t> traces;
      std::string plot_title;
      std::string x_label, y_label;

      bool   extents_x_set = false, extents_y_set = false;
      double x_lo = 0.0, x_hi = 1.0;
      double y_lo = 0.0, y_hi = 1.0;

      bool draw_x_axis = true,  draw_y_axis = true;
      bool draw_x_tick = true,  draw_y_tick = true;

      // Widget
      GtkWidget *drawing_area = nullptr;
      int width_req, height_req;
      // margins used to position the plot inside the widget
      double left_margin = 50.0, right_margin = 18.0;
      double top_margin  = 28.0, bottom_margin = 40.0;

      // GTK4 draw callback
      static void draw_cb(GtkDrawingArea *da, cairo_t *cr,
                          int w, int h, gpointer user_data);
      void render(cairo_t *cr, int w, int h);

      // Rendering helpers
      struct plot_box_t { double x, y, w, h; };
      plot_box_t plot_box(int w, int h) const;
      void update_auto_extents();
      void draw_axes   (cairo_t *cr, const plot_box_t &pb);
      void draw_ticks  (cairo_t *cr, const plot_box_t &pb);
      void draw_title  (cairo_t *cr, int w);
      void draw_labels (cairo_t *cr, int w, int h, const plot_box_t &pb);
      void draw_trace  (cairo_t *cr, const plot_box_t &pb,
                        const graph_trace_t &t);

      // world (data) -> widget pixels
      inline double w2p_x(double x, const plot_box_t &pb) const {
         return pb.x + (x - x_lo) * pb.w / (x_hi - x_lo);
      }
      inline double w2p_y(double y, const plot_box_t &pb) const {
         return pb.y + pb.h - (y - y_lo) * pb.h / (y_hi - y_lo);
      }
      static double nice_tick(double range);
   };
}

#endif
