
#ifndef HAVE_GENERIC_DISPLAY_OBJECT_HH
#define HAVE_GENERIC_DISPLAY_OBJECT_HH

#include "coot-utils.hh"  // for colour_holder

namespace coot { 

   class generic_display_line_t { 
   public:
     std::pair<clipper::Coord_orth, clipper::Coord_orth> coords;
     generic_display_line_t(const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords_in) { 
       coords = coords_in;
     }
   };

   // a container for lines of the same colour and width
   class generic_display_line_set_t { 
   public:
     coot::colour_holder colour;
     std::string colour_name;
     int width;
     std::vector<generic_display_line_t> lines;
     generic_display_line_set_t(const coot::colour_holder &colour_in, 
				const std::string &colour_name_in,
				const int &width_in) { 
       colour = colour_in;
       colour_name = colour_name_in;
       width = width_in;
     }
     void add_line(const generic_display_line_t &line) { 
       lines.push_back(line);
     } 
   };

   // a container for points of the same colour and size
   class generic_display_point_set_t { 
   public:
     coot::colour_holder colour;
     std::string colour_name;
     int size;
     std::vector<clipper::Coord_orth> points;
     generic_display_point_set_t(const colour_holder &colour_in, 
				 const std::string &colour_name_in,
				 const int &size_in) { 
       colour = colour_in;
       colour_name = colour_name_in;
       size = size_in;
     }
     void add_point(const clipper::Coord_orth &point) { 
       points.push_back(point);
     } 
   };

   class generic_display_object_t {

   public:
      short int is_displayed_flag;
      short int is_closed_flag; // don't make buttons for closed display objects
      std::string name; // name the object, use for indexing, perhaps
      std::vector<generic_display_line_set_t> lines_set;
      std::vector<generic_display_point_set_t> points_set;
      std::vector<int> GL_display_list_handles;
      generic_display_object_t(const std::string &n) { 
	 name = n;
	 is_displayed_flag = 0;
	 is_closed_flag = 0;
      }
      void add_line(const coot::colour_holder &colour_in,
		    const std::string &colour_name,
		    const int &width_in, 
		    const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords_in);
      void add_point(const coot::colour_holder &colour_in,
		     const std::string &colour_name,
		     const int &size_in, 
		     const clipper::Coord_orth &coords_in);
      void raster3d(std::ofstream &render_stream) const;
      void clear() {
	 lines_set.clear();
	 points_set.clear();
      } 
      void close_yourself() { 
	 std::string name = "Closed Generic Display Object";
	 clear();
	 is_closed_flag = 1;
      } 
      static coot::colour_holder colour_values_from_colour_name(const std::string &colour_name);
   };

   class generic_text_object_t { 
   public:
     int handle;
     std::string s;
     float x, y, z;
     generic_text_object_t(const std::string &s_in, int handle_in, 
			   float x_in, float y_in, float z_in) { 
       handle = handle_in;
       s = s_in;
       x = x_in;
       y = y_in;
       z = z_in;
     }
   };

}

#endif // HAVE_GENERIC_DISPLAY_OBJECT_HH
