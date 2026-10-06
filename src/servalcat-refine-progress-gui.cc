/*
 * src/servalcat-refine-progress-gui.cc
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include <gtk/gtk.h>
#include <cairo/cairo.h>

#include "graphics-info.h"
#include "coot-utils/servalcat-refine-progress.hh"
#include "widget-from-builder.hh"
#include "servalcat-refine-progress-gui.hh"

namespace {

   // One labelled data series within a plot panel.
   struct series_t {
      std::vector<std::pair<double, double> > pts; // (x, y)
      double r, g, b;
      std::string label;
   };

   enum class x_axis_mode_t { cycle, inv_res_sq };

   // the drawing area in the pop-out dialog (nullptr until the dialog is made)
   GtkWidget *popout_drawing_area = nullptr;
   GtkWidget *popout_window = nullptr;

   GdkRGBA widget_foreground(GtkWidget *w) {
      GdkRGBA fg = {0.0f, 0.0f, 0.0f, 1.0f}; // default: a dark fg (light background)
#if GTK_MINOR_VERSION >= 10
      gtk_widget_get_color(w, &fg);
#endif
      return fg;
   }

   // A light foreground colour implies a dark background.
   bool background_is_dark(const GdkRGBA &fg) {
      double lum = 0.2126 * fg.red + 0.7152 * fg.green + 0.0722 * fg.blue;
      return lum > 0.5;
   }

   std::vector<std::pair<double, double> >
   to_double_series(const std::vector<std::pair<int, double> > &in) {
      std::vector<std::pair<double, double> > out;
      out.reserve(in.size());
      for (const auto &p : in)
         out.push_back(std::make_pair(static_cast<double>(p.first), p.second));
      return out;
   }

   // Draw one plot panel (axes, autoscaled y, title, legend and series) into the
   // rectangle (ox, oy, w, h). y is autoscaled to the data; x is either the cycle
   // number (integer ticks from 0) or 1/d^2 (labelled as resolution in Angstrom).
   void draw_panel(cairo_t *cr, double ox, double oy, double w, double h,
                   const std::string &title,
                   const std::vector<series_t> &series,
                   const GdkRGBA &fg, bool dark_bg,
                   x_axis_mode_t x_mode) {

      const double ml = 38.0, mr = 10.0, mt = 16.0, mb = 20.0;
      double pw = w - ml - mr;
      double ph = h - mt - mb;
      if (pw <= 2.0 || ph <= 2.0) return;
      double plx = ox + ml;        // plot left
      double ply = oy + mt;        // plot top
      double pbx = plx + pw;       // plot right
      double pby = ply + ph;       // plot bottom

      // title
      cairo_set_font_size(cr, 10.0);
      cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.95);
      cairo_move_to(cr, ox + 2.0, oy + 11.0);
      cairo_show_text(cr, title.c_str());

      // x and y data ranges
      bool have_pts = false;
      double x_min = 0.0, x_max = 1.0, y_min = 0.0, y_max = 1.0;
      for (const auto &s : series) {
         for (const auto &p : s.pts) {
            if (! have_pts) { x_min = x_max = p.first; y_min = y_max = p.second; have_pts = true; }
            x_min = std::min(x_min, p.first); x_max = std::max(x_max, p.first);
            y_min = std::min(y_min, p.second); y_max = std::max(y_max, p.second);
         }
      }

      if (x_mode == x_axis_mode_t::cycle) { x_min = 0.0; if (x_max < 1.0) x_max = 1.0; }
      if (x_max <= x_min) x_max = x_min + 1.0;

      // pad y a little; guard a flat series
      if (y_max <= y_min) { y_max = y_min + 1.0; y_min -= 1.0; }
      double y_pad = 0.08 * (y_max - y_min);
      y_min -= y_pad; y_max += y_pad;

      auto x_to_px = [&](double x){ return plx + (x - x_min) / (x_max - x_min) * pw; };
      auto y_to_px = [&](double y){ return pby - (y - y_min) / (y_max - y_min) * ph; };

      // axes (left and bottom)
      cairo_set_line_width(cr, 1.0);
      cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.6);
      cairo_move_to(cr, plx, ply); cairo_line_to(cr, plx, pby);
      cairo_move_to(cr, plx, pby); cairo_line_to(cr, pbx, pby);
      cairo_stroke(cr);

      // y tick labels (min and max)
      cairo_set_font_size(cr, 8.0);
      cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.75);
      char lab[48];
      g_snprintf(lab, sizeof lab, "%.3g", y_max);
      cairo_move_to(cr, ox + 2.0, ply + 7.0);  cairo_show_text(cr, lab);
      g_snprintf(lab, sizeof lab, "%.3g", y_min);
      cairo_move_to(cr, ox + 2.0, pby);        cairo_show_text(cr, lab);

      // x tick labels
      if (x_mode == x_axis_mode_t::cycle) {
         int n_cyc = static_cast<int>(std::lround(x_max));
         int step = 1; while ((n_cyc / step) > 6) step++;
         for (int c = 0; c <= n_cyc; c += step) {
            double px = x_to_px(static_cast<double>(c));
            cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.5);
            cairo_move_to(cr, px, pby); cairo_line_to(cr, px, pby + 3.0); cairo_stroke(cr);
            g_snprintf(lab, sizeof lab, "%d", c);
            cairo_text_extents_t ext; cairo_text_extents(cr, lab, &ext);
            cairo_move_to(cr, px - ext.width/2.0 - ext.x_bearing, pby + 12.0);
            cairo_show_text(cr, lab);
         }
      } else {
         // inverse resolution squared: label as resolution d (A) at a few ticks
         for (int i = 0; i <= 4; i++) {
            double frac = i / 4.0;
            double s = x_min + frac * (x_max - x_min);
            double px = x_to_px(s);
            cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.5);
            cairo_move_to(cr, px, pby); cairo_line_to(cr, px, pby + 3.0); cairo_stroke(cr);
            if (s > 1e-6) {
               double d = 1.0 / std::sqrt(s);
               g_snprintf(lab, sizeof lab, "%.1f", d);
               cairo_text_extents_t ext; cairo_text_extents(cr, lab, &ext);
               double tx = px - ext.width/2.0 - ext.x_bearing;
               if (tx < ox) tx = ox;
               if (tx + ext.width > ox + w) tx = ox + w - ext.width;
               cairo_move_to(cr, tx, pby + 12.0);
               cairo_show_text(cr, lab);
            }
         }
      }

      if (! have_pts) {
         cairo_set_font_size(cr, 9.0);
         cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.5);
         cairo_move_to(cr, plx + 6.0, ply + ph/2.0);
         cairo_show_text(cr, "waiting for data...");
         return;
      }

      // series: polyline + markers
      for (const auto &s : series) {
         if (s.pts.empty()) continue;
         cairo_set_line_width(cr, 1.6);
         cairo_set_source_rgb(cr, s.r, s.g, s.b);
         for (std::size_t i=0; i<s.pts.size(); i++) {
            double px = x_to_px(s.pts[i].first), py = y_to_px(s.pts[i].second);
            if (i == 0) cairo_move_to(cr, px, py); else cairo_line_to(cr, px, py);
         }
         cairo_stroke(cr);
         for (const auto &p : s.pts) {
            cairo_arc(cr, x_to_px(p.first), y_to_px(p.second), 1.8, 0.0, 2.0 * M_PI);
            cairo_fill(cr);
         }
      }

      // legend (top-right, inside the plot)
      cairo_set_font_size(cr, 8.0);
      double ly = ply + 8.0;
      for (const auto &s : series) {
         if (s.label.empty()) continue;
         cairo_text_extents_t ext; cairo_text_extents(cr, s.label.c_str(), &ext);
         double lx = pbx - ext.width - 16.0;
         cairo_set_line_width(cr, 2.0);
         cairo_set_source_rgb(cr, s.r, s.g, s.b);
         cairo_move_to(cr, lx, ly - 3.0); cairo_line_to(cr, lx + 11.0, ly - 3.0); cairo_stroke(cr);
         cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.85);
         cairo_move_to(cr, lx + 14.0, ly);
         cairo_show_text(cr, s.label.c_str());
         ly += 11.0;
      }
      (void) dark_bg;
   }

   // The series colours, adapted to a light or dark background.
   struct palette_t {
      double rwork[3], rfree[3], ccwork[3], ccfree[3], mll[3];
      double geom_bond[3], geom_angle[3], geom_chiral[3], geom_planar[3];
   };

   palette_t make_palette(bool dark) {
      palette_t p;
      auto set = [](double *c, double r, double g, double b){ c[0]=r; c[1]=g; c[2]=b; };
      if (dark) {
         set(p.rwork,      0.45, 0.70, 1.00); set(p.rfree,      1.00, 0.55, 0.40);
         set(p.ccwork,     0.40, 0.85, 0.65); set(p.ccfree,     0.95, 0.75, 0.35);
         set(p.mll,        0.80, 0.60, 1.00);
         set(p.geom_bond,  0.45, 0.80, 0.55); set(p.geom_angle, 0.50, 0.70, 1.00);
         set(p.geom_chiral,0.95, 0.70, 0.40); set(p.geom_planar,1.00, 0.55, 0.55);
      } else {
         set(p.rwork,      0.15, 0.35, 0.80); set(p.rfree,      0.85, 0.30, 0.12);
         set(p.ccwork,     0.10, 0.55, 0.35); set(p.ccfree,     0.80, 0.55, 0.10);
         set(p.mll,        0.50, 0.25, 0.70);
         set(p.geom_bond,  0.15, 0.55, 0.25); set(p.geom_angle, 0.20, 0.35, 0.75);
         set(p.geom_chiral,0.75, 0.45, 0.10); set(p.geom_planar,0.80, 0.20, 0.20);
      }
      return p;
   }

   // the shared, parsed progress data (owned by graphics_info_t)
   coot::servalcat_refine_progress_t *progress_data() {
      return graphics_info_t::servalcat_refine_progress_p;
   }

   series_t summary_series(const char *key, const double *col, const std::string &label) {
      series_t s;
      s.r = col[0]; s.g = col[1]; s.b = col[2]; s.label = label;
      if (progress_data())
         s.pts = to_double_series(progress_data()->summary_series(key));
      return s;
   }

   series_t geom_series(const char *restraint, const double *col, const std::string &label) {
      series_t s;
      s.r = col[0]; s.g = col[1]; s.b = col[2]; s.label = label;
      if (progress_data())
         s.pts = to_double_series(progress_data()->geom_series("r.m.s.Z", restraint));
      return s;
   }

   series_t binned_series(const char *column, const double *col, const std::string &label) {
      series_t s;
      s.r = col[0]; s.g = col[1]; s.b = col[2]; s.label = label;
      if (progress_data() && ! progress_data()->empty())
         s.pts = progress_data()->binned_series(progress_data()->last_cycle_index(), column);
      return s;
   }

   // Compact overlay: R-factors (top) and geometry r.m.s.Z (bottom), vs cycle.
   void draw_progress_overlay(GtkDrawingArea *area, cairo_t *cr,
                              int width, int height, G_GNUC_UNUSED gpointer user_data) {
      GdkRGBA fg = widget_foreground(GTK_WIDGET(area));
      bool dark = background_is_dark(fg);
      palette_t pal = make_palette(dark);
      double h2 = height / 2.0;

      std::vector<series_t> r_factors = {
         summary_series("Rwork", pal.rwork, "Rwork"),
         summary_series("Rfree", pal.rfree, "Rfree") };
      draw_panel(cr, 0, 0, width, h2, "R factors", r_factors, fg, dark, x_axis_mode_t::cycle);

      std::vector<series_t> geom = {
         geom_series("Bond distances, non H", pal.geom_bond,  "bond"),
         geom_series("Bond angles, non H",    pal.geom_angle, "angle") };
      draw_panel(cr, 0, h2, width, h2, "Geometry r.m.s.Z", geom, fg, dark, x_axis_mode_t::cycle);
   }

   // Full pop-out dialog: a 2-column grid of panels.
   void draw_progress_full(GtkDrawingArea *area, cairo_t *cr,
                           int width, int height, G_GNUC_UNUSED gpointer user_data) {
      GdkRGBA fg = widget_foreground(GTK_WIDGET(area));
      bool dark = background_is_dark(fg);
      palette_t pal = make_palette(dark);

      const int n_rows = 3;
      double cw = width / 2.0;
      double rh = height / static_cast<double>(n_rows);

      // row 0: R factors vs cycle | binned R factors vs resolution
      std::vector<series_t> r_factors = {
         summary_series("Rwork", pal.rwork, "Rwork"),
         summary_series("Rfree", pal.rfree, "Rfree") };
      draw_panel(cr, 0, 0, cw, rh, "R factors vs cycle", r_factors, fg, dark, x_axis_mode_t::cycle);

      std::vector<series_t> binned_r = {
         binned_series("Rwork", pal.rwork, "Rwork"),
         binned_series("Rfree", pal.rfree, "Rfree") };
      draw_panel(cr, cw, 0, cw, rh, "R factors vs resolution (A)", binned_r, fg, dark, x_axis_mode_t::inv_res_sq);

      // row 1: CCF vs cycle | binned CCF vs resolution
      std::vector<series_t> ccf = {
         summary_series("CCFworkavg", pal.ccwork, "CCFwork"),
         summary_series("CCFfreeavg", pal.ccfree, "CCFfree") };
      draw_panel(cr, 0, rh, cw, rh, "CCF vs cycle", ccf, fg, dark, x_axis_mode_t::cycle);

      std::vector<series_t> binned_cc = {
         binned_series("CCFwork", pal.ccwork, "CCFwork"),
         binned_series("CCFfree", pal.ccfree, "CCFfree") };
      draw_panel(cr, cw, rh, cw, rh, "CCF vs resolution (A)", binned_cc, fg, dark, x_axis_mode_t::inv_res_sq);

      // row 2: -LL vs cycle | geometry r.m.s.Z vs cycle
      std::vector<series_t> mll = { summary_series("-LL", pal.mll, "-LL") };
      draw_panel(cr, 0, 2*rh, cw, rh, "-log likelihood vs cycle", mll, fg, dark, x_axis_mode_t::cycle);

      std::vector<series_t> geom = {
         geom_series("Bond distances, non H", pal.geom_bond,   "bond"),
         geom_series("Bond angles, non H",    pal.geom_angle,  "angle"),
         geom_series("Chiral centres",        pal.geom_chiral, "chiral"),
         geom_series("Planar groups",         pal.geom_planar, "planar") };
      draw_panel(cr, cw, 2*rh, cw, rh, "Geometry r.m.s.Z vs cycle", geom, fg, dark, x_axis_mode_t::cycle);
   }

   void redraw_graphs() {
      GtkWidget *da = widget_from_builder("servalcat-refine-progress-drawing-area");
      if (da) gtk_widget_queue_draw(da);
      if (popout_drawing_area) gtk_widget_queue_draw(popout_drawing_area);
   }

   // ---- buttons -----------------------------------------------------------

   void on_popout_clicked(G_GNUC_UNUSED GtkButton *button, G_GNUC_UNUSED gpointer user_data) {
      if (! popout_window) {
         popout_window = gtk_window_new();
         gtk_window_set_title(GTK_WINDOW(popout_window), "Servalcat Refinement Progress");
         gtk_window_set_default_size(GTK_WINDOW(popout_window), 760, 560);
         popout_drawing_area = gtk_drawing_area_new();
         gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(popout_drawing_area),
                                        draw_progress_full, nullptr, nullptr);
         gtk_window_set_child(GTK_WINDOW(popout_window), popout_drawing_area);
         // keep the pointers valid: just hide on close rather than destroy
         g_signal_connect(popout_window, "close-request",
                          G_CALLBACK(+[](GtkWindow *w, G_GNUC_UNUSED gpointer d) -> gboolean {
                             gtk_widget_set_visible(GTK_WIDGET(w), FALSE);
                             return TRUE; // handled - do not destroy
                          }), nullptr);
      }
      gtk_widget_set_visible(popout_window, TRUE);
      gtk_window_present(GTK_WINDOW(popout_window));
      redraw_graphs();
   }

   void on_close_clicked(G_GNUC_UNUSED GtkButton *button, G_GNUC_UNUSED gpointer user_data) {
      servalcat_refine_progress_hide();
   }

   // Set the overlay drawing-area draw function and connect the buttons - once.
   void ensure_overlay_wired() {
      static bool done = false;
      if (done) return;

      GtkWidget *da = widget_from_builder("servalcat-refine-progress-drawing-area");
      if (da)
         gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(da), draw_progress_overlay, nullptr, nullptr);

      GtkWidget *popout_button = widget_from_builder("servalcat-refine-progress-popout-button");
      if (popout_button)
         g_signal_connect(popout_button, "clicked", G_CALLBACK(on_popout_clicked), nullptr);

      GtkWidget *close_button = widget_from_builder("servalcat-refine-progress-close-button");
      if (close_button)
         g_signal_connect(close_button, "clicked", G_CALLBACK(on_close_clicked), nullptr);

      done = true;
   }
}

// ---- public interface (servalcat-refine-progress-gui.hh) ------------------

void servalcat_refine_progress_show(const std::string &stats_file_name) {

   if (! graphics_info_t::servalcat_refine_progress_p)
      graphics_info_t::servalcat_refine_progress_p = new coot::servalcat_refine_progress_t;
   else
      graphics_info_t::servalcat_refine_progress_p->cycles.clear();

   graphics_info_t::servalcat_refine_progress_stats_file_name = stats_file_name;

   ensure_overlay_wired();

   GtkWidget *frame = widget_from_builder("servalcat-refine-progress-frame");
   if (frame) gtk_widget_set_visible(frame, TRUE);

   redraw_graphs();
}

void servalcat_refine_progress_update() {

   if (! graphics_info_t::servalcat_refine_progress_p) return;
   const std::string &fn = graphics_info_t::servalcat_refine_progress_stats_file_name;
   if (fn.empty()) return;

   // a failed/partial read leaves the previous cycles untouched (see parse_json_string)
   graphics_info_t::servalcat_refine_progress_p->parse_json_file(fn);
   redraw_graphs();
}

void servalcat_refine_progress_hide() {

   GtkWidget *frame = widget_from_builder("servalcat-refine-progress-frame");
   if (frame) gtk_widget_set_visible(frame, FALSE);
}
