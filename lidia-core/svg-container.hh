#ifndef COOT_LIDIA_CORE_SVG_CONTAINER_HH
#define COOT_LIDIA_CORE_SVG_CONTAINER_HH

#include <string>
#include <iostream> // for debugging - remove later.

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
#if 1
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

   std::string make_viewbox_string() const {

      std::string viewBox_string = "viewBox=" + std::string("\"") +
         std::to_string(min_x) + std::string(" ") +
         std::to_string(min_y) + std::string(" ") +
         std::to_string(max_x) + std::string(" ") +
         std::to_string(max_y) + std::string("\"");
      return viewBox_string;
   }

   void add(const std::string &s) {
      svg += s;
   }

   void add(const svg_container_t &svgc_in) {
      svg += svgc_in.svg;
      update_bounds(svgc_in.min_x, svgc_in.min_y, svgc_in.max_x, svgc_in.max_y);
   }

   std::string compose() const {

      std::string s = svg_header_1;
      s += make_viewbox_string();
      s += svg_header_2;
      s += svg;
      s += svg_footer;
      return s;
   }

};



#endif // SVG_COOT_LIDIA_CORE_CONTAINER_HH
