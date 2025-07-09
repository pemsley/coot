#ifndef COOT_LIDIA_CORE_SVG_CONTAINER_HH
#define COOT_LIDIA_CORE_SVG_CONTAINER_HH

#include <string>
#include <iostream> // for debugging - remove later.

#include "lig-build.hh"

class svg_container_t {

public:
   svg_container_t() { init(); }
   explicit svg_container_t(const std::string &s) : svg(s) { init(); }
   std::string svg;
   std::string svg_header_1;
   std::string svg_header_2;
   std::string svg_footer;
   float min_x;
   float min_y;
   float max_x;
   float max_y;

   void init() {
      svg_header_1 = "<svg xmlns=\"http://www.w3.org/2000/svg\"\n    xmlns:xlink=\"http://www.w3.org/1999/xlink\" ";
      svg_header_2 = ">\n";
      svg_footer = "</svg>\n";
      svg.reserve(2048);
      min_x =  100000.0;
      min_y =  100000.0;
      max_x = -100000.0;
      max_y = -100000.0;
   }

   void set_bounds(float min_xx, float min_yy, float max_xx, float max_yy) {
      min_x = min_xx;
      min_y = min_yy;
      max_x = max_xx;
      max_y = max_yy;
   }

   bool update_bounds(float min_xx, float min_yy, float max_xx, float max_yy) {
      bool bounding_box_updated = false;
#if 0
      if (min_xx < min_x) { std::cout << "debug:: min_x was " << min_x << " now " << min_xx << std::endl; }
      if (min_yy < min_y) { std::cout << "debug:: min_y was " << min_y << " now " << min_yy << std::endl; }
      if (max_xx > max_x) { std::cout << "debug:: max_x was " << max_x << " now " << max_xx << std::endl; }
      if (max_yy > max_y) { std::cout << "debug:: max_y was " << max_y << " now " << max_yy << std::endl; }
#endif

      if (min_xx < min_x) { min_x = min_xx; bounding_box_updated = true; }
      if (min_yy < min_y) { min_y = min_yy; bounding_box_updated = true; }
      if (max_xx > max_x) { max_x = max_xx; bounding_box_updated = true; }
      if (max_yy > max_y) { max_y = max_yy; bounding_box_updated = true; }
      return bounding_box_updated;
   }

   // used to add a key to flev
   void add_to_y_bounds(double yy) {
      max_y += yy;
   }

   std::string make_viewbox_string() const {

      float width  = max_x - min_x;
      float height = max_y - min_y;
      std::string viewBox_string = "viewBox=" + std::string("\"") +
         std::to_string(min_x) + std::string(" ") +
         std::to_string(min_y) + std::string(" ") +
         std::to_string(width) + std::string(" ") +
         std::to_string(height) + std::string("\"");
      return viewBox_string;
   }

   void add(const std::string &s) {
      svg += s;
   }

   void prepend(const svg_container_t &svgc_in) {
      svg = svgc_in.svg + svg;
      update_bounds(svgc_in.min_x, svgc_in.min_y, svgc_in.max_x, svgc_in.max_y);
   }

   void add(const svg_container_t &svgc_in) {
      svg += svgc_in.svg;
      update_bounds(svgc_in.min_x, svgc_in.min_y, svgc_in.max_x, svgc_in.max_y);
   }

   void add_comment(const std::string &comment) {
      svg += "<!-- ";
      svg += comment;
      svg += " -->\n";
   }

   void add_line(const lig_build::pos_t &p1, const lig_build::pos_t &p2,
                 double line_width, const std::string &stroke_colour, bool dashed) {
      std::string s;
      s += "   <line x1=\"" + std::to_string(p1.x) + "\" y1=\"" + std::to_string(-p1.y) + "\" ";
      s += "x2=\"" + std::to_string(p2.x) + "\" y2=\"" + std::to_string(-p2.y) + "\" ";
      s += "style=\"stroke:" + stroke_colour + ";stroke-width:" + std::to_string(line_width) + ";";
      s += "stroke-linecap:round;";
      if (dashed)
         s += "stroke-dasharray:0.1,0.2;";
      s += "\" />\n";
      svg += s;
   };



   std::string compose(bool add_background_rect) const {

      auto quoted = [] (float v) {
         std::string s = std::to_string(v);
         return std::string("'" + s + "'");
      };

      auto make_background_rect = [quoted] (float min_x, float min_y, float max_x, float max_y) {

         float w = max_x - min_x;
         float h = max_y - min_y;
         std::string s = "<!-- background-rectangle -->\n";
         s += "   <rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width=";
         // s += "400.0' height='400.0' x='0.0' y='0.0'";
         s += quoted(w);
         s += " ";
         s += "height=";
         s += quoted(h);
         s += " ";
         s += "x=";
         s += quoted(min_x);
         s += " ";
         s += "y=";
         s += quoted(min_y);
         s += "> </rect>\n";
         return s;
      };

      if (max_x > min_x) {
         std::string s = svg_header_1;
         s += make_viewbox_string();
         s += svg_header_2;
         if (add_background_rect)
            s += make_background_rect(min_x, min_y, max_x, max_y);
         s += svg;
         s += svg_footer;
         return s;
      } else {
         return "";
      }
   }

};



#endif // SVG_COOT_LIDIA_CORE_CONTAINER_HH
