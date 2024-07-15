
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>
#include <filesystem>

#include <signal.h>

#include "json.hpp"
using json = nlohmann::json;

class point_t {
public:
   point_t(double x, double y) : x(x), y(y) {}
   explicit point_t(std::pair<double, double> &xy) : x(xy.first), y(xy.second) {}
   explicit point_t(std::pair<double, double> &xy, double x_scale, double y_scale, double x_offset, double y_offset) :
      x(x_offset + x_scale * xy.first), y(y_offset + y_scale * xy.second) {}
   double x;
   double y;
   point_t canvas_convert() const {
      double new_x = x;
      double new_y = -y + 100;
      return point_t(new_x, new_y);
   }
   point_t operator+(const point_t &p) const {
      return point_t(x + p.x, y + p.y);
   }
   void operator +=(const point_t &p) {
      x += p.x;
      y += p.y;
   }
};

std::string double_to_string(const double &d, unsigned int n_dec_pl) {
   std::stringstream s;
   s << std::right << std::fixed;
   s << std::setprecision(n_dec_pl);
   s << d;
   std::string ss = s.str();
   return ss;
}


std::string make_line(const point_t &p1, const point_t &p2, const std::string &colour, double width,
                      bool dashed, bool do_arrow_head) {

   point_t p1c = p1.canvas_convert();
   point_t p2c = p2.canvas_convert();

   do_arrow_head = false;

   std::string line;
   line += "   <line x1=\"";
   line += std::to_string(p1c.x);
   line += "\" y1=\"";
   line += std::to_string(p1c.y);
   line += "\" x2=\"";
   line += std::to_string(p2c.x);
   line += "\" y2=\"";
   line += std::to_string(p2c.y);
   line += "\" " ;
   line += "style=\"stroke:";
   line += colour;
   line += ";stroke-width:";
   line += std::to_string(width);
   line += "\"";
   if (dashed)
      line += " stroke-dasharray=\"1\"";
   if (do_arrow_head)
      // line += " marker-end=\"url(#arrow)\"";
      line += " marker-end=url(#arrow)";
   line += " />\n";

   return line;
}

std::vector<std::string> make_grid_lines(unsigned int n_cycles,
                                         double x_scale, double y_scale,
                                         double x_offset, double y_offset) {

   unsigned int n_lines = 5;
   if (n_cycles > n_lines) n_lines = n_cycles;

   std::vector<std::string> v;
   for (unsigned int i=1; i<=n_lines; i++) { // not at x = 0;
      std::pair<double, double> xy_1(i, 0.0);
      std::pair<double, double> xy_2(i, 1.2);
      point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
      point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
      std::string l = make_line(pt_1, pt_2, "#aaaaaa", 0.4, true, false);
      v.push_back(l);
   }
   return v;
}

class summary_data_t {
public:
   summary_data_t() : cycle_number(0), fsc_average(0), mll(0) {}
   summary_data_t(int c, double f, double m) : cycle_number(c), fsc_average(f), mll(m) {}
   int cycle_number;
   double fsc_average;
   double mll;
};


class summary_data_container_t {
public:
   summary_data_container_t() { init(); }
   std::vector<summary_data_t> data;
   double x_scale;
   double y_scale;
   double x_offset;
   double y_offset;

   void init() {
      x_scale = 10.0;
      y_scale = 80.0;
      x_offset = 22.0;
      y_offset = 20.0; // move it down the page
   }
   void add(const summary_data_t &d) { data.push_back(d); }
   void output() {

      std::string output_file_name = "summary.svg";
      std::string s = make_complete_svg();
      std::ofstream f(output_file_name);
      f << s;
      f.close();
   }

   std::string make_complete_svg() {

      std::string s;
      std::string svg_header_1 = "<svg xmlns=\"http://www.w3.org/2000/svg\"\n   xmlns:xlink=\"http://www.w3.org/1999/xlink\" ";
      std::string svg_header_2 = ">\n";

      std::string svg_footer = "</svg>\n";

      double min_x = -10;
      double min_y = -10;
      double max_x = 400;
      double max_y = 200;
      double x_offset_orig = x_offset;

      std::string viewBox_string = "viewBox=" + std::string("\"") +
         std::to_string(min_x) + std::string(" ") +
         std::to_string(min_y) + std::string(" ") +
         std::to_string(max_x) + std::string(" ") +
         std::to_string(max_y) + std::string("\"");

      s += svg_header_1;
      s += viewBox_string;
      s += svg_header_2;

      s += graph_internals("fsc_average");

      x_offset += 140;
      s += graph_internals("mll");
      s += svg_footer;
      x_offset = x_offset_orig;
      return s;
   }

   std::string make_svg() {
      double x_offset_orig = x_offset;
      std::string s;
      s += graph_internals("fsc_average");
      x_offset += 150;
      s += graph_internals("mll");
      x_offset = x_offset_orig;
      return s;
   }

   std::vector<std::string> make_grid_lines(unsigned int n_cycles) const {
      return ::make_grid_lines(n_cycles, x_scale, y_scale, x_offset, y_offset);
   }

   std::vector<std::string> make_graph_lines(const std::string &type) const {

      // convert this to polylines

      std::cout << "make_graph_lines() " << type << std::endl;

      std::string colour = "#662222";
      std::vector<std::string> v;

      if (type == "fsc_average") {
         double this_graph_scale = 1.2;
         std::vector<std::pair<double, double> > data_points;
         for (unsigned int i=0; i<data.size(); i++)
            data_points.push_back(std::make_pair(static_cast<double>(i), data[i].fsc_average));
         for (unsigned int i=0; i<data.size(); i++) {
            if (i < (data.size() - 1)) {
               point_t p1(data_points[i  ], x_scale, this_graph_scale * y_scale, x_offset, y_offset);
               point_t p2(data_points[i+1], x_scale, this_graph_scale * y_scale, x_offset, y_offset);
               std::string line = make_line(p1, p2, colour, 0.5, false, false);
               // std::cout << "fsc " << i << " " << data_points[i].first << " " << data_points[i].second << " " << line << std::endl;
               v.push_back(line);
            }
         }
      }

      if (type == "mll") {
         colour = "#222299";
         double y_scale_mll = 0.00000024;
         std::vector<std::pair<double, double> > data_points;
         for (unsigned int i=0; i<data.size(); i++)
            data_points.push_back(std::make_pair(static_cast<double>(i), data[i].mll));
         for (unsigned int i=0; i<data.size(); i++) {
            if (i < (data.size() - 1)) {
               // std::cout << "mll data_points: " << data_points[i].first << "  " << data_points[i].second << std::endl;
               point_t p1(data_points[i  ], x_scale, y_scale * y_scale_mll, x_offset, y_offset);
               point_t p2(data_points[i+1], x_scale, y_scale * y_scale_mll, x_offset, y_offset);
               std::string line = make_line(p1, p2, colour, 0.5, false, false);
               v.push_back(line);
            }
         }
      }
      return v;
   }

   std::string make_x_axis_label(const std::string &x_axis_text) const {

      point_t p(40 + x_offset, y_offset - 14);
      point_t pc = p.canvas_convert();
      std::string s = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"5\" x=\"";
      s += std::to_string(pc.x);
      s += "\" y=\"";
      s += std::to_string(pc.y);
      s += "\">";
      s += x_axis_text;
      s += "</text>\n";
      return s;
   }

   std::string make_y_axis_label(const std::string &y_axis_text) const {
      point_t p(-27 + x_offset, 50 + y_offset);
      point_t pc = p.canvas_convert();
      std::string s = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"5\" x=\"";
      s += std::to_string(pc.x);
      s += "\" y=\"";
      s += std::to_string(pc.y);
      s += "\">";
      s += y_axis_text;
      s += "</text>\n";
      return s;
   }

   std::string make_tick_marks_x_axis(unsigned int n_cycles) const {
      std::string s;
      unsigned int n_lines = 5;
      if (n_cycles > n_lines) n_lines = n_cycles;
      for (unsigned int i=1; i<=n_lines; i++) { // not at x = 0;
         std::pair<double, double> xy_1(i, 0.0);
         std::pair<double, double> xy_2(i, -0.03);
         point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
         point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
         std::string l = make_line(pt_1, pt_2, "#222222", 0.4, false, false);
         s += l;
      }
      return s;
   }

   std::string make_tick_marks_y_axis() const {
      std::string s = "   \n";
      unsigned int n_ticks = 5;
      for (unsigned int i=0; i<=n_ticks; i++) {
         double f = static_cast<double>(i) / static_cast<double>(n_ticks);
         double y = 1.2 * f;
         std::pair<double, double> xy_1(0.0,  y);
         std::pair<double, double> xy_2(-0.3, y);
         point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
         point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
         std::string l = make_line(pt_1, pt_2, "#222222", 0.4, false, false);
         s += l;
      }
      return s;
   }

   std::string make_x_axis_tick_labels(unsigned int n_cycles, double x_scale, double x_offset) const {
      std::string s;
      unsigned int n_lines = 5;
      if (n_cycles > n_lines) n_lines = n_cycles;
      for (unsigned int i=1; i<=n_lines; i++) { // not at x = 0;
         point_t p( x_scale * i + x_offset - 1.0, y_offset - 6);
         point_t pc = p.canvas_convert();
         std::string l = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"3\" ";
         l += "x=\"";
         l += std::to_string(pc.x);
         l += "\" y=\"";
         l += std::to_string(pc.y);
         l += "\">";
         l += std::to_string(i);
         l += "</text>\n";
         s += l;
      }
      return s;
   }

   std::string make_y_axis_tick_labels(const std::string &graph_type, double x_scale, double x_offset) const {
      std::string s;
      unsigned int n_lines = 5;
      for (unsigned int i=0; i<=n_lines; i++) { // not at x = 0;
         double f = static_cast<double>(i) / static_cast<double>(n_lines);
         double v = 0.0;
         double y = 0.0;
         double tw_offset = 12;
         std::string value_as_string;
         if (graph_type == "fsc_average") {
            tw_offset = 12;
            v = f;
            y = 96.0 * f; // 100 is too much
            value_as_string = double_to_string(v, 1);
         }
         if (graph_type == "mll") {
            tw_offset = 10; // 24;
            v = 5.0 * f;
            y = 96.0 * f; // 100 is too much
            value_as_string = double_to_string(v, 0);
         }
         point_t p(x_offset - tw_offset, y + y_offset - 2);
         point_t pc = p.canvas_convert();
         std::string l = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"6\" ";
         l += "x=\"";
         l += std::to_string(pc.x);
         l += "\" y=\"";
         l += std::to_string(pc.y);
         l += "\">";
         l += value_as_string;
         l += "</text>\n";
         s += l;
      }
      return s;
   }

   std::string make_main_title(const std::string &title_string) const {
      point_t p(12 + x_offset, 110 + y_offset);
      point_t pc = p.canvas_convert();
      std::string s =  "   <text font-family=\"Helvetica, sans-serif\" font-size=\"9\" ";
      s += "x=\"";
      s += std::to_string(pc.x);
      s += "\" y=\"";
      s += std::to_string(pc.y);
      s += "\"";
      s += ">";
      s += title_string;
      s += "</text>\n";
      return s;
   }

   std::string graph_internals(const std::string &graph_type) {

      point_t origin(x_offset, y_offset);
      point_t x_axis_max(x_offset + 100 + 10, y_offset);
      point_t y_axis_max(0 + x_offset, y_offset + 100 + 5);
      std::string line_x_axis = make_line(origin, x_axis_max, "#111111", 1, false, true);
      std::string line_y_axis = make_line(origin, y_axis_max, "#111111", 1, false, true);
      std::string parts = line_x_axis + line_y_axis;
      int n_cycles = 10;
      std::vector<std::string> grid_lines = make_grid_lines(n_cycles);
      for (const auto &line : grid_lines)
         parts += line;
      if (graph_type == "fsc_average") {
         std::vector<std::string> graph_lines = make_graph_lines("fsc_average");
         for (const auto &line : graph_lines)
            parts += line;
      }
      if (graph_type == "mll") {
         std::vector<std::string> graph_lines = make_graph_lines("mll");
         for (const auto &line : graph_lines)
            parts += line;
      }

      std::string x_axis_text = "Cycle Number";
      std::string y_axis_text;
      if (graph_type == "fsc_average") y_axis_text = "FSC";
      if (graph_type == "mll")         y_axis_text = "-LL/M";

      std::string x_axis_label = make_x_axis_label(x_axis_text);
      std::string y_axis_label = make_y_axis_label(y_axis_text);

      parts += x_axis_label;
      parts += y_axis_label;

      std::string tick_marks_x = make_tick_marks_x_axis(n_cycles);
      std::string tick_marks_y = make_tick_marks_y_axis();

      parts += tick_marks_x;
      parts += tick_marks_y;

      std::string tick_labels_x = make_x_axis_tick_labels(n_cycles, x_scale, x_offset);
      std::string tick_labels_y = make_y_axis_tick_labels(graph_type, x_scale, x_offset);

      parts += tick_labels_x;
      parts += tick_labels_y;

      std::string graph_title;
      if (graph_type == "fsc_average") graph_title = "Summary Data FSC";
      if (graph_type == "mll")         graph_title = "Summary Data -LL";

      std::string title = make_main_title(graph_title);

      parts += title;

      return parts;
   }
};

class binned_data_t {
public:
   binned_data_t(const double &reso_d, const double &fsc_FC_full, const double &Rcmplx_FC_full,
                 const double &cc_FC_full, const double &mcos_FC_full,
                 const double &power_FP, const double &power_FC) : reso(reso_d) {
                  data["fsc_FC_full"] = fsc_FC_full;
                  data["Rcmplx_FC_full"] = Rcmplx_FC_full;
                  data["cc_FC_full"] = cc_FC_full;
                  data["mcos_FC_full"] = mcos_FC_full;
                  data["power_FP"] = power_FP;
                  data["power_FC"] = power_FC;
                 }
   double reso;
   std::map<std::string, double> data;

};

class binned_data_container_t {
public:
   binned_data_container_t() { init(); }
   std::map<unsigned int, std::vector<binned_data_t> > data;
   void add(int cycle, const binned_data_t &b) { data[cycle].push_back(b); }
   void output() {

      std::string output_file_name = "binned-data.svg";
      std::string s = make_complete_svg();
      std::ofstream f(output_file_name);
      f << s;
      f.close();
   }

   double x_scale;
   double y_scale;
   double x_offset;
   double y_offset;

   void init() {
      x_scale = 10.0;
      y_scale = 80.0;
      x_offset = 22.0;
      y_offset = -180.0; // move it down the page
   }

   // in Angstroms
   double get_max_resolution() const {
      double r = 0.1; // A
      if (data.begin() != data.end()) {
         r = 10000;
         const std::vector<binned_data_t> &v = data.begin()->second;
         for (const auto &bin : v) {
            if (bin.reso < r) {
               r = bin.reso;
            }
         }
      }
      return r;
   }

   std::string make_graph_lines(const std::string &graph_type, unsigned int cycle_number, const std::string &colour) const {
      std::string s;
      std::vector<std::pair<double, double> > data_points;
      double this_graph_scale = 100.0;
      if (graph_type == "power_FP") this_graph_scale = 0.000001;
      if (graph_type == "power_FC") this_graph_scale = 0.000001;

      auto data_vec = data.at(cycle_number);
      if (true) {
         for (unsigned int i=0; i<data_vec.size(); i++) {
            const binned_data_t &binned_data = data_vec.at(i);
            double inv_res = 1.0/binned_data.reso; // reso is in A
            double inv_res_sq = inv_res * inv_res;
            double xx = x_offset + static_cast<double>(i) * 2.4;
            try {
               point_t p(xx, y_offset + this_graph_scale * binned_data.data.at(graph_type));
               point_t pc = p.canvas_convert();
               data_points.push_back(std::make_pair(pc.x, pc.y));
            }
            catch(const std::exception &e) {
               std::cerr << " make_graph_lines() A " << e.what() << " key: " << graph_type << '\n';
            }
         }
      }

      if (! data_points.empty()) {
         std::string polyline = "<polyline points=\"";

         for (unsigned int i=0; i<data_points.size(); i++) {
            polyline += std::to_string(data_points[i].first);
            polyline += ",";
            polyline += std::to_string(data_points[i].second);
            if (i < (data_points.size() -1))
               polyline += ",";
         }
         polyline += "\" fill=\"none\" stroke=\"";
         polyline += colour;
         polyline += "\" />\n";
         s += polyline;
      }
      return s;
   }

   std::string graph_internals(const std::string &graph_type, unsigned int cycle_number, const std::string &colour) const {
      std::string s;
      s += make_graph_lines(graph_type, cycle_number, colour);
      return s;
   }

   std::string make_tick_marks_x_axis() const {
      std::string s;
      unsigned int n_ticks = 5;
      for (unsigned int i=0; i<=n_ticks; i++) {
         double f = static_cast<double>(i) / static_cast<double>(n_ticks);
         double x = f * 10.5; // not 10.
         std::pair<double, double> xy_1(x,   0);
         std::pair<double, double> xy_2(x, -0.05);
         point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
         point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
         std::string l = make_line(pt_1, pt_2, "#222222", 0.4, false, false);
         s += l;
      }
      return s;
   }

   std::string make_tick_marks_y_axis() const {
      std::string s;
      unsigned int n_ticks = 5;
      for (unsigned int i=0; i<=n_ticks; i++) {
         double f = static_cast<double>(i) / static_cast<double>(n_ticks);
         double y = f * 1.2;
         std::pair<double, double> xy_1(0.0,  y);
         std::pair<double, double> xy_2(-0.3, y);
         point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
         point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
         std::string l = make_line(pt_1, pt_2, "#222222", 0.4, false, false);
         s += l;
      }
      return s;
   }

   std::string make_x_axis_tick_labels(double max_resolution) const {

      std::string s;
      unsigned int n_lines = 5;
      double inv_resol_sq = 1.0/(max_resolution * max_resolution);
      for (unsigned int i=0; i<=n_lines; i++) { // not at x = 0;
         double f = static_cast<double>(i) / static_cast<double>(n_lines);
         double v0 = f * inv_resol_sq;
         double v = std::sqrt(1.0/v0);
         double x = f * 100.5;
         double y = 0.0;
         double tw_offset = 0;
         std::string value_as_string = double_to_string(v,1);
         point_t p(x + x_offset - tw_offset, y_offset - 12);
         point_t pc = p.canvas_convert();
         std::string l = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"6\" ";
         l += "x=\"";
         l += std::to_string(pc.x);
         l += "\" y=\"";
         l += std::to_string(pc.y);
         l += "\">";
         l += value_as_string;
         l += "</text>\n";
         s += l;
      }
      return s;
   }

   std::string make_y_axis_tick_labels() const {
      std::string s;
      unsigned int n_lines = 5;
      for (unsigned int i=0; i<=n_lines; i++) { // not at x = 0;
         double f = static_cast<double>(i) / static_cast<double>(n_lines);
         double v = f;
         double y = 96.0 * f;
         double tw_offset = 12;
         std::string value_as_string = double_to_string(v,1);
         point_t p(x_offset - tw_offset, y + y_offset - 2);
         point_t pc = p.canvas_convert();
         std::string l = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"6\" ";
         l += "x=\"";
         l += std::to_string(pc.x);
         l += "\" y=\"";
         l += std::to_string(pc.y);
         l += "\">";
         l += value_as_string;
         l += "</text>\n";
         s += l;
      }
      return s;
   }

   std::string make_main_title(const std::string &title_string) const {

      point_t p(12 + x_offset, 110 + y_offset);
      point_t pc = p.canvas_convert();
      std::string s =  "   <text font-family=\"Helvetica, sans-serif\" font-size=\"9\" ";
      s += "x=\"";
      s += std::to_string(pc.x);
      s += "\" y=\"";
      s += std::to_string(pc.y);
      s += "\"";
      s += ">";
      s += title_string;
      s += "</text>\n";
      return s;
   }

   std::string make_complete_svg() {

      std::string s;
      std::string svg_header_1 = "<svg xmlns=\"http://www.w3.org/2000/svg\"\n   xmlns:xlink=\"http://www.w3.org/1999/xlink\" ";
      std::string svg_header_2 = ">\n";

      std::string svg_footer = "</svg>\n";

      double min_x = -10;
      double min_y = -30;
      double max_x = 400;
      double max_y = 300;

      std::string viewBox_string = "viewBox=" + std::string("\"") +
         std::to_string(min_x) + std::string(" ") +
         std::to_string(min_y) + std::string(" ") +
         std::to_string(max_x) + std::string(" ") +
         std::to_string(max_y) + std::string("\"");

      s += svg_header_1;
      s += viewBox_string;
      s += svg_header_2;

      // add graphs here
      s += make_svg();

      s += svg_footer;
      return s;
   }

   std::string make_svg() {

      std::string s;

      std::vector<std::pair<point_t, std::string> > graph_names;
      graph_names.push_back(std::make_pair(point_t( 0,     0), "fsc_FC_full"));
      graph_names.push_back(std::make_pair(point_t(140,    0), "Rcmplx_FC_full"));
      graph_names.push_back(std::make_pair(point_t(  0, -140), "cc_FC_full"));
      graph_names.push_back(std::make_pair(point_t(140, -140), "mcos_FC_full"));
      // graph_names.push_back(std::make_pair(point_t(  0, -240), "power_FP"));
      // graph_names.push_back(std::make_pair(point_t(140, -240), "power_FC"));

      double x_offset_orig = x_offset;
      double y_offset_orig = y_offset;
      double max_resolution = get_max_resolution();

      for (unsigned int i_graph=0; i_graph<graph_names.size(); i_graph++) {

         const std::string &graph_name = graph_names[i_graph].second;
         const point_t &offset = graph_names[i_graph].first;
         x_offset = x_offset_orig + offset.x;
         y_offset = y_offset_orig + offset.y;

         point_t origin(x_offset, y_offset);
         point_t x_axis_max(x_offset + 100 + 10, y_offset);
         point_t y_axis_max(0 + x_offset, y_offset + 100 + 5);
         std::string line_x_axis = make_line(origin, x_axis_max, "#111111", 1, false, true);
         std::string line_y_axis = make_line(origin, y_axis_max, "#111111", 1, false, true);
         s += line_x_axis;
         s += line_y_axis;

         std::string tick_marks_x = make_tick_marks_x_axis();
         std::string tick_marks_y = make_tick_marks_y_axis();
         s += tick_marks_x;
         s += tick_marks_y;

         std::string y_axis_tick_labels = make_y_axis_tick_labels();
         std::string x_axis_tick_labels = make_x_axis_tick_labels(max_resolution);
         s += y_axis_tick_labels;
         s += x_axis_tick_labels;

         std::string title_string = make_main_title(graph_name);
         s += title_string;

         unsigned int last_cycle_number = data.size() -1;

         std::map<unsigned int, std::vector<binned_data_t> >::const_iterator it;
         for (it=data.begin(); it!=data.end(); ++it) {
            unsigned int cycle_number = it->first;
            std::string colour = "grey";
            if (cycle_number == last_cycle_number)
               colour = "#111111";
            s += graph_internals(graph_name, cycle_number, colour);
         }
      }
      x_offset = x_offset_orig;
      y_offset = y_offset_orig;
      return s;
   }
};

class geom_data_t {
   public:
   geom_data_t() { }
   // outer keys "r.m.s.d." and "r.m.s.Z"
   std::map<std::string, std::map<std::string, double> > geom_data;
};

class geom_data_container_t {
   public:
   geom_data_container_t() { init(); }
   std::map<unsigned int, geom_data_t> data;
   std::vector<std::string> keys;
   void init() {
      keys = {
         "Bond distances, non H",
         "Bond distances, H",
         "Bond angles, non H",
         "Bond angles, H",
         "Torsion angles, period 1",
         "Torsion angles, period 2",
         "Torsion angles, period 3",
         "Torsion angles, period 6",
         "Chiral centres",
         "Planar groups",
         "VDW nonbonded",
         "VDW torsion",
         "VDW hbond",
         "B values (bond)",
         "B values (angle)",
         "B values (others)" };

      x_scale = 10.0;
      y_scale = 80.0;
      x_offset = 22.0;
      y_offset = 20.0; // move it down the page
      outer_keys = { "r.m.s.d.", "r.m.s.Z"};
      outer_keys.erase(outer_keys.begin()); // only r.m.s.Z

   }
   double x_scale;
   double y_scale;
   double x_offset;
   double y_offset;
   std::vector<std::string> outer_keys;


   std::string make_tick_marks_x_axis(unsigned int n_cycles) const {
      std::string s;
      unsigned int n_lines = 5;
      if (n_cycles > n_lines) n_lines = n_cycles;
      for (unsigned int i=1; i<=n_lines; i++) { // not at x = 0;
         std::pair<double, double> xy_1(i, 0.0);
         std::pair<double, double> xy_2(i, -0.03);
         point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
         point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
         std::string l = make_line(pt_1, pt_2, "#222222", 0.4, false, false);
         s += l;
      }
      return s;
   }

   std::string make_tick_marks_y_axis() const {
      std::string s = "   \n";
      unsigned int n_ticks = 5;
      for (unsigned int i=0; i<=n_ticks; i++) {
         double f = static_cast<double>(i) / static_cast<double>(n_ticks);
         double y = 1.2 * f;
         std::pair<double, double> xy_1(0.0,  y);
         std::pair<double, double> xy_2(-0.3, y);
         point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
         point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
         std::string l = make_line(pt_1, pt_2, "#222222", 0.4, false, false);
         s += l;
      }
      return s;
   }

   std::string make_x_axis_tick_labels(unsigned int n_cycles) const {
      std::string s;
      unsigned int n_lines = 5;
      if (n_cycles > n_lines) n_lines = n_cycles;
      for (unsigned int i=1; i<=n_lines; i++) { // not at x = 0;
         point_t p( x_scale * i + x_offset - 1.0, y_offset - 6);
         point_t pc = p.canvas_convert();
         std::string l = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"3\" ";
         l += "x=\"";
         l += std::to_string(pc.x);
         l += "\" y=\"";
         l += std::to_string(pc.y);
         l += "\">";
         l += std::to_string(i);
         l += "</text>\n";
         s += l;
      }
      return s;
   }

   std::string make_y_axis_tick_labels(const std::string &graph_type) const {
      std::string s;
      unsigned int n_lines = 5;
      for (unsigned int i=0; i<=n_lines; i++) { // not at x = 0;
         double f = static_cast<double>(i) / static_cast<double>(n_lines);
         double v = 0.0;
         double y = 0.0;
         double tw_offset = 12;
         std::string value_as_string;
         point_t p(x_offset - tw_offset, y + y_offset - 2);
         point_t pc = p.canvas_convert();
         std::string l = "   <text font-family=\"Helvetica, sans-serif\" font-size=\"6\" ";
         l += "x=\"";
         l += std::to_string(pc.x);
         l += "\" y=\"";
         l += std::to_string(pc.y);
         l += "\">";
         l += value_as_string;
         l += "</text>\n";
         s += l;
      }
      return s;
   }


   std::string make_y_axis_tick_labels() {
      std::string s;
      return s;
   }

   std::string make_main_title(const std::string &title_string) const {

      point_t p(12 + x_offset, 110 + y_offset);
      point_t pc = p.canvas_convert();
      std::string s =  "   <text font-family=\"Helvetica, sans-serif\" font-size=\"8\" ";
      s += "x=\"";
      s += std::to_string(pc.x);
      s += "\" y=\"";
      s += std::to_string(pc.y);
      s += "\"";
      s += ">";
      s += title_string;
      s += "</text>\n";
      return s;
   }

   std::vector<std::string> make_grid_lines(unsigned int n_cycles) {

      unsigned int n_lines = 5;
      if (n_cycles > n_lines) n_lines = n_cycles;

      std::vector<std::string> v;
      for (unsigned int i=1; i<=n_lines; i++) { // not at x = 0;
         std::pair<double, double> xy_1(i, 0.0);
         std::pair<double, double> xy_2(i, 1.2);
         point_t pt_1(xy_1, x_scale, y_scale, x_offset, y_offset);
         point_t pt_2(xy_2, x_scale, y_scale, x_offset, y_offset);
         std::string l = make_line(pt_1, pt_2, "#999999", 0.3, true, false);
         v.push_back(l);
      }
      return v;
   }

   std::string make_graph_lines(const std::string &outer_key,
                                 const std::string &graph_name,
                                 const std::string &colour) const {
      std::string s;
      std::map<unsigned int, geom_data_t>::const_iterator it;
      std::vector<point_t> data_points;
      for (it=data.begin(); it!=data.end(); ++it) {
         unsigned int cycle_number = it->first;
         const geom_data_t &geom_data = it->second;
         double this_graph_scale = 3.0;
         if (graph_name == "Bond distances, non H") this_graph_scale = 2000.0;
         if (graph_name == "Bond distances, H")     this_graph_scale = 2000.0;
         if (graph_name == "Bond angles, non H")    this_graph_scale = 30.0;
         if (graph_name == "Bond angles, H")        this_graph_scale = 30.0;
         if (graph_name == "Chiral centres")        this_graph_scale =  500.0;
         if (graph_name == "Planar groups")         this_graph_scale = 1000.0;
         if (graph_name == "VDW nonbonded")         this_graph_scale =  100.0;
         if (graph_name == "VDW torsion")           this_graph_scale =  100.0;
         if (graph_name == "VDW hbond")             this_graph_scale =  100.0;
         if (outer_key == "r.m.s.Z") this_graph_scale = 40.0;
         try {
            double y = geom_data.geom_data.at(outer_key).at(graph_name);
            double xx = x_offset + cycle_number * 10.0;
            double yy = y_offset + y * this_graph_scale;
            point_t p(xx, yy);
            point_t pc = p.canvas_convert();
            data_points.push_back(pc);
         }
         catch (const std::exception &e) {
            std::cout << "make_graph_lines() B " << e.what() << " key " << outer_key << " graph_name " << graph_name
               << std::endl;
         }
      }

      if (! data_points.empty()) {
         std::string polyline = "<polyline points=\"";

         for (unsigned int i=0; i<data_points.size(); i++) {
            polyline += std::to_string(data_points[i].x);
            polyline += ",";
            polyline += std::to_string(data_points[i].y);
            if (i < (data_points.size() -1))
               polyline += ",";
         }
         polyline += "\" fill=\"none\" stroke=\"";
         polyline += colour;
         polyline += "\" />\n";
         s += polyline;
      }
      return s;
   }


   std::string make_svg() {
      std::string s;

      std::vector<std::pair<point_t, std::string> > graph_names = {
         std::make_pair(point_t(  0, 0), "Bond distances, non H"),
         std::make_pair(point_t(120, 0), "Bond distances, H"),
         std::make_pair(point_t(240, 0), "Bond angles, non H"),
         std::make_pair(point_t(360, 0), "Bond angles, H"),
         std::make_pair(point_t(  0, -140), "Torsion angles, period 1"),
         std::make_pair(point_t(120, -140), "Torsion angles, period 2"),
         std::make_pair(point_t(240, -140), "Torsion angles, period 3"),
         std::make_pair(point_t(360, -140), "Torsion angles, period 6"),
         std::make_pair(point_t(  0, -280), "Chiral centres"),
         std::make_pair(point_t(120, -280), "Planar groups"),
         std::make_pair(point_t(240, -280), "VDW nonbonded"),
         std::make_pair(point_t(360, -280), "VDW torsion"),
         std::make_pair(point_t(  0, -420), "VDW hbond"),
         std::make_pair(point_t(120, -420), "B values (bond)"),
         std::make_pair(point_t(240, -420), "B values (angle)"),
         std::make_pair(point_t(360, -420), "B values (others)"),
      };

      double x_offset_orig = x_offset;
      double y_offset_orig = y_offset;

      unsigned int n_cycles = data.rbegin()->first;

      unsigned int key_count = 0;
      for (const auto &outer_key : outer_keys) {
         key_count++;

         for (unsigned int i_graph=0; i_graph<graph_names.size(); i_graph++) {

            const std::string &graph_name = graph_names[i_graph].second;
            point_t offset = graph_names[i_graph].first;
            if (true) {
               offset += point_t(300, 0);
            }
            x_offset = x_offset_orig + offset.x;
            y_offset = y_offset_orig + offset.y;

            point_t origin(x_offset, y_offset);
            point_t x_axis_max(x_offset + 100 + 10, y_offset);
            point_t y_axis_max(0 + x_offset, y_offset + 100 + 5);
            std::string line_x_axis = make_line(origin, x_axis_max, "#111111", 1, false, true);
            std::string line_y_axis = make_line(origin, y_axis_max, "#111111", 1, false, true);
            s += line_x_axis;
            s += line_y_axis;

            std::string tick_marks_x = make_tick_marks_x_axis(n_cycles);
            std::string tick_marks_y = make_tick_marks_y_axis();
            s += tick_marks_x;
            s += tick_marks_y;

            std::string x_axis_tick_labels = make_x_axis_tick_labels(n_cycles);
            std::string y_axis_tick_labels = make_y_axis_tick_labels();
            s += y_axis_tick_labels;
            s += x_axis_tick_labels;

            std::cout << "calling make_graph_lines with n_cycles " << n_cycles << std::endl;
            std::vector<std::string> grid_lines = make_grid_lines(n_cycles);
            for (const auto &gl : grid_lines) s += gl;

            std::string title = outer_key + " " + graph_name;
            std::string title_string = make_main_title(title);
            s += title_string;

            // outer keys "r.m.s.d." and "r.m.s.Z"
            // in geom_data_t: std::map<std::string, std::map<std::string, double> > geom_data;
            // std::map<unsigned int, geom_data_t> data;

            std::string colour = "#333333";
            s += make_graph_lines(outer_key, graph_name, colour);

            x_offset = x_offset_orig;
            y_offset = y_offset_orig;
         }
      }
      return s;
   }

   std::string make_complete_svg() {

      std::string s;
      std::string svg_header_1 = "<svg xmlns=\"http://www.w3.org/2000/svg\"\n   xmlns:xlink=\"http://www.w3.org/1999/xlink\" ";
      std::string svg_header_2 = ">\n";

      std::string svg_footer = "</svg>\n";

      double min_x =  100; // -10
      double min_y = -130; // -30
      double max_x =  800; // 400
      double max_y =  800;

      std::string viewBox_string = "viewBox=" + std::string("\"") +
         std::to_string(min_x) + std::string(" ") +
         std::to_string(min_y) + std::string(" ") +
         std::to_string(max_x) + std::string(" ") +
         std::to_string(max_y) + std::string("\"");

      s += svg_header_1;
      s += viewBox_string;
      s += svg_header_2;

      // add graphs here
      s += make_svg();

      s += svg_footer;
      return s;
   }

   void output() {
      std::string output_file_name = "geom-data.svg";
      std::string s = make_complete_svg();
      std::ofstream f(output_file_name);
      f << s;
      f.close();
   }
};

void
parse_cycle(json j,
            summary_data_container_t *summary_data_container_p,
            binned_data_container_t *bdc_p,
            geom_data_container_t *geom_data_container_p) {

   int nth_cycle = j["Ncyc"];
   json j_data = j["data"];
   json j_summary = j_data["summary"];
   float fsc_average = j_summary["FSCaverage"];
   float mll = j_summary["-LL"];
   std::cout << "-------------- nth_cycle: " << nth_cycle << " fsc: " << fsc_average << " -ll: " << mll << std::endl;
   summary_data_container_p->add(summary_data_t(nth_cycle, fsc_average, mll));

   // "binned" data inside "data"
   json j_binned = j_data["binned"];
   unsigned int s = j_binned.size();
   for (unsigned int i=0; i<s; i++) {
      json binned_data = j_binned[i];
      double d_min          = binned_data["d_min"];
      double d_max          = binned_data["d_max"];
      double fsc_FC_full    = binned_data["fsc_FC_full"];
      double Rcmplx_FC_full = binned_data["Rcmplx_FC_full"];
      double cc_FC_full     = binned_data["cc_FC_full"];
      double mcos_FC_full   = binned_data["mcos_FC_full"];
      double power_FP       = binned_data["power_FP"];
      double power_FC       = binned_data["power_FC"];
      unsigned int n_coeffs = binned_data["ncoeffs"];
      double d_mid = 2.0/(1.0/d_min + 1.0/d_max);
      binned_data_t bd(d_mid, fsc_FC_full, Rcmplx_FC_full, cc_FC_full, mcos_FC_full, power_FP, power_FC);
      bdc_p->add(nth_cycle, bd);
   }

   geom_data_container_t &geom_data_container = *geom_data_container_p;
   std::vector<std::string> data_types = {"r.m.s.d.", "r.m.s.Z"};
   data_types.erase(data_types.begin()); // only use r.m.s.Z
   json j_geom = j["geom"];
   json j_geom_summary = j_geom["summary"];

   for (const auto &data_type : data_types) {
      json j_data_type = j_geom_summary[data_type];
      if (j_data_type.is_null()) {
         std::cout << "========= failed to find " << data_type << std::endl;
      } else {
         for (const auto &key : geom_data_container.keys) {
            try {
               double x = j_data_type[key];
               geom_data_container.data[nth_cycle].geom_data[data_type][key] = x;
            }
            catch (const std::exception &e) {
               // the value hasn't been set yet
               std::cout << "soft-error: parse_cycle(): " << e.what() << std::endl;
            }
         }
      }
   }
}

void make_consolidated_graph_set(summary_data_container_t *summary_data_container,
                                 binned_data_container_t *binned_data_container,
                                 geom_data_container_t *geom_data_container,
                                 const std::string &svg_file_name,
                                 bool wrap_in_refresh_html) {

   auto make_svg = [summary_data_container, binned_data_container, geom_data_container] {
      std::string s;
      s += summary_data_container->make_svg();
      s += binned_data_container->make_svg();
      s += geom_data_container->make_svg();
      return s;
   };

   auto make_complete_svg = [make_svg] {

      std::string s;
      std::string svg_header_1 = "<svg xmlns=\"http://www.w3.org/2000/svg\"\n   xmlns:xlink=\"http://www.w3.org/1999/xlink\" ";
      std::string svg_header_2 = ">\n";

      std::string svg_footer = "</svg>\n";

      double min_x =  -10;
      double min_y =  -60;
      double max_x =  860;
      double max_y =  800;

      std::string viewBox_string = "viewBox=" + std::string("\"") +
         std::to_string(min_x) + std::string(" ") +
         std::to_string(min_y) + std::string(" ") +
         std::to_string(max_x) + std::string(" ") +
         std::to_string(max_y) + std::string("\"");

      s += svg_header_1;
      s += viewBox_string;
      s += svg_header_2;

      // add graphs here
      s += make_svg();

      s += svg_footer;
      return s;
   };

   auto output = [] (const std::string &s, const std::string &output_file_name) {
      std::ofstream f(output_file_name);
      f << s;
      f.close();
   };

   std::string s;
   s = make_complete_svg();

   if (wrap_in_refresh_html) {
      std::string html_top = "<!DOCTYPE html><html><head><meta http-equiv=\"refresh\" content=\"2\"></head><body>\n";
      std::string html_foot = "</body><html>\n";
      std::string ss = html_top + s + html_foot;
      s = ss;
   }

   output(s, svg_file_name);
   std::cout << "debug:: in make_consolidated_graph_set() wrote " << svg_file_name << std::endl;
}

// throw an exception on unable to convert
int
string_to_int(const std::string &s) {

   int i;
   std::istringstream sstream(s);

   if (sstream>>i) {
      return i;
   } else {
      std::string mess = "Cannot convert \"";
      mess += s;
      mess += "\" to an integer";
      throw std::runtime_error(mess);
   }
}

std::string file_to_string(const std::string &file_name) {
   std::fstream f(file_name);
   std::string s;
   f.seekg(0, std::ios::end);
   s.reserve(f.tellg());
   f.seekg(0, std::ios::beg);
   s.assign((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
   return s;
}

int get_cycle_number_max(json j) {
   int size = j.size();
   int cycle_number = -1;
   for (int i=0; i<size; i++) {
      json j_c = j[i];
      int c = j_c["Ncyc"];
      // std::cout << " c "  << c << std::endl;
      if (c > cycle_number)
         cycle_number = c;
   }

   return cycle_number;
}

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if defined _MSC_VER
#define snprintf _snprintf
#else
#include <unistd.h>
#endif

void track_process_make_graphs_svg(const std::string &json_file_name, int pid) {

   bool continue_status = true;
   time_t latest_modified_time = 0;
   summary_data_container_t summary_data_container;
   binned_data_container_t binned_data_container;
   geom_data_container_t geom_data_container;
   bool found_begining = false;
   while (continue_status) {
      int kill_status = kill(pid,0);
      if (kill_status == 0) {
         if (std::filesystem::exists(json_file_name)) {

            struct stat buf;
            int stat_status = stat(json_file_name.c_str(), &buf);
            if (stat_status == 0) {
               time_t mtime = buf.st_mtime;
               if (mtime > latest_modified_time) {
                  latest_modified_time = mtime;
                  // std::cout << "updated " << mtime << std::endl;
                  std::string s = file_to_string(json_file_name);
                  json j = json::parse(s);
                  int cycle_number = get_cycle_number_max(j);
                  if (cycle_number == 0) found_begining = true;
                  std::cout << "cycle_number: " << cycle_number << " found_begining " << found_begining << std::endl;
                  if (found_begining) {
                     std::cout << "... parsing..." << std::endl;
                     json item = j.at(cycle_number);
                     parse_cycle(item, &summary_data_container, &binned_data_container, &geom_data_container);

                     // std::string svg_file_name = "servalcat-graphs.svg";
                     std::string svg_file_name = "servalcat-graphs.html";
                     make_consolidated_graph_set(&summary_data_container, &binned_data_container,
                                                 &geom_data_container, svg_file_name, true);
                  }
               }
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(500));

         } else {
            std::cout << "No such file " << json_file_name << std::endl;
            // continue_status = false;
         }
      } else {
         std::cout << "process is dead" << std::endl;
         continue_status = false;
      }
   }
}

int main(int argc, char **argv) {

   int status = 0;
   bool done = false;

   if (argc == 2) {

      std::string json_file_name = argv[1];
      try {

         summary_data_container_t summary_data_container;
         binned_data_container_t binned_data_container;
         geom_data_container_t geom_data_container;

         std::string s = file_to_string(json_file_name);
         json j = json::parse(s);
         int cycle_index = 0;
         while (true) {

            try {
               json item = j.at(cycle_index);
               parse_cycle(item, &summary_data_container, &binned_data_container, &geom_data_container);

               // next time
               cycle_index++;
            }
            catch (const nlohmann::detail::out_of_range &oor) {
               // std::cout << "out of range " << cycle_index << std::endl;
               break;
            }
         }

         summary_data_container.output();
         binned_data_container.output();
         geom_data_container.output();

         std::string svg_file_name = "servalcat-graphs.svg";
         make_consolidated_graph_set(&summary_data_container, &binned_data_container,
                                     &geom_data_container, svg_file_name, false);
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING:: main(): " << e.what() << std::endl;
         status = 2; // file not found
      }

      done = true;
   }
   if (argc == 3) {
      std::string json_file_name = argv[1];
      std::string pid_string = argv[2];
      try {
         int pid = string_to_int(pid_string);
         track_process_make_graphs_svg(json_file_name, pid);
      }
      catch (const std::exception &e) {
         std::cout << e.what() << " Failed  to convert to int" << std::endl;
      }


      done = true;
   }

   if (! done) {
      std::cout << "Usage: servalcat-tracker json-file-name" << std::endl;
   }

   return status;
}
