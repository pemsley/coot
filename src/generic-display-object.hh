
#ifndef HAVE_GENERIC_DISPLAY_OBJECT_HH
#define HAVE_GENERIC_DISPLAY_OBJECT_HH

#include "coot-utils.hh"  // for colour_holder

#include "coot-colour.hh"

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

      class arrow_t {
      public:
	 arrow_t() { fract_head_size = 0.3;}
	 arrow_t(const clipper::Coord_orth &pt1, const clipper::Coord_orth &pt2) {
	    start_point = pt1;
	    end_point = pt2;
	    fract_head_size = 0.3;
	 }
	 clipper::Coord_orth start_point;
	 clipper::Coord_orth end_point;
	 coot::colour_t col;
	 float fract_head_size;
      };
      class sphere_t {
      public:
	 sphere_t() {}
	 sphere_t(const clipper::Coord_orth &centre_in, float r) {
	    centre = centre_in;
	    radius = r;
	 }
	 clipper::Coord_orth centre;
	 coot::colour_t col;
	 float radius;
      };
      class torus_t {
      public:
	 torus_t() {}
	 torus_t(const clipper::Coord_orth &pt1,
		 const clipper::Coord_orth &pt2,
		 float r1, float r2) {
	    start_point = pt1;
	    end_point   = pt2;
	    radius_1 = r1;
	    radius_2 = r2;
	    n_ring_atoms = 6;
	 }
	 clipper::Coord_orth start_point;
	 clipper::Coord_orth end_point;
	 coot::colour_t col;
	 float radius_1;
	 float radius_2;
	 int n_ring_atoms;
      };
      
      bool is_displayed_flag;
      bool is_closed_flag; // don't make buttons for closed display objects
      bool is_solid_flag;
      bool is_transparent_flag;
      float opacity; // 0.0 -> 1.0
      std::string name; // name the object, use for indexing, perhaps
      std::vector<generic_display_line_set_t> lines_set;
      std::vector<generic_display_point_set_t> points_set;
      std::vector<arrow_t> arrows;
      std::vector<sphere_t> spheres;
      std::vector<torus_t> tori;
      std::vector<int> GL_display_list_handles;
      generic_display_object_t(const std::string &n) { 
	 name = n;
	 is_displayed_flag = false;
	 is_closed_flag = false;
	 opacity = 1.0;
	 is_transparent_flag = false;
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
	 is_closed_flag = true;
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
