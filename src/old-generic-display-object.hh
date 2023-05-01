
#ifndef HAVE_OLD_GENERIC_DISPLAY_OBJECT_HH
#define HAVE_OLD_GENERIC_DISPLAY_OBJECT_HH

#include "utils/colour-holder.hh"
#include "utils/dodec.hh"

#include "coot-colour.hh"

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name);

namespace coot {

   class old_generic_display_line_t {
   public:
      std::pair<clipper::Coord_orth, clipper::Coord_orth> coords;
      explicit old_generic_display_line_t(const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords_in) : coords(coords_in) {}
   };

   // a container for lines of the same colour and width
   class old_generic_display_line_set_t {
   public:
      colour_holder colour;
      std::string colour_name;
      int width;
      std::vector<old_generic_display_line_t> lines;
      old_generic_display_line_set_t(const colour_holder &colour_in,
                                     const std::string &colour_name_in,
                                     const int &width_in) : colour(colour_in), colour_name(colour_name_in) {
         width = width_in;
      }
      void add_line(const old_generic_display_line_t &line) {
         lines.push_back(line);
      }
   };

   // a container for points of the same colour and size
   class old_generic_display_point_set_t {
   public:
      colour_holder colour;
      std::string colour_name;
      int size;
      std::vector<clipper::Coord_orth> points;
      old_generic_display_point_set_t(const colour_holder &colour_in,
                                      const std::string &colour_name_in,
                                      const int &size_in) : colour(colour_in), colour_name(colour_name_in) {
         size = size_in;
      }
      void add_point(const clipper::Coord_orth &point) {
         points.push_back(point);
      }
   };

   class old_generic_display_object_t {

   public:

      class arrow_t {
      public:
         arrow_t() { fract_head_size = 0.3;}
         arrow_t(const clipper::Coord_orth &pt1, const clipper::Coord_orth &pt2) : start_point(pt1), end_point(pt2) {
            fract_head_size = 0.3;
         }
         clipper::Coord_orth start_point;
         clipper::Coord_orth end_point;
         colour_t col;
         float fract_head_size;
      };
      class sphere_t {
      public:
         sphere_t() {}
         sphere_t(const clipper::Coord_orth &centre_in, float r) : centre(centre_in) {
            radius = r;
         }
         clipper::Coord_orth centre;
         colour_t col;
         float radius;
      };
      class torus_t {
      public:
         torus_t() {}
         torus_t(const clipper::Coord_orth &pt1,
                 const clipper::Coord_orth &pt2,
                 float r1, float r2) : start_point(pt1), end_point(pt2) {
            radius_1 = r1;
            radius_2 = r2;
            n_ring_atoms = 6;
         }
         clipper::Coord_orth start_point;
         clipper::Coord_orth end_point;
         colour_t col;
         float radius_1;
         float radius_2;
         int n_ring_atoms;
      };
      // arc is part of a torus
      class arc_t {
      public:
         arc_t(float start_angle_in, float end_angle_in,
               const clipper::Coord_orth &start_point_in,
               const clipper::Coord_orth &start_dir_in,
               const clipper::Coord_orth &normal_in,
               float radius_in, float radius_inner_in) :
            start_point(start_point_in), normal(normal_in), start_dir(start_dir_in) {
            start_angle = start_angle_in;
            end_angle = end_angle_in;
            radius = radius_in;
            radius_inner = radius_inner_in;
         }
         clipper::Coord_orth start_point;
         clipper::Coord_orth normal;
         clipper::Coord_orth start_dir;
         float start_angle;
         float end_angle;
         colour_t col;
         float radius;
         float radius_inner;
      };

      class dodec_t {
      public:
         dodec_t(const dodec &d_in, double size_in, const clipper::Coord_orth &pos_in) : d(d_in), position(pos_in) {
            size = size_in;
         }
         dodec d;
         double size;
         clipper::Coord_orth position;
         colour_holder col;
      };

      class pentakis_dodec_t { // perhaps this should inherit from above
      public:
         pentakis_dodec_t(const pentakis_dodec &pkdd_in, double size_in,
                          const clipper::Coord_orth &pos_in) :pkdd(pkdd_in), position(pos_in) {
            size = size_in;
         }
         pentakis_dodec pkdd;
         double size;
         clipper::Coord_orth position;
         colour_holder col;
      };

      enum {UNDEFINED = -1, INTERMEDIATE_ATOMS=-9};
      int  imol;
      bool is_displayed_flag;
      bool is_closed_flag; // don't make buttons for closed display objects
      bool is_solid_flag;
      bool is_transparent_flag;
      float opacity; // 0.0 -> 1.0
      std::string name; // name the object, use for indexing, perhaps
      std::vector<old_generic_display_line_set_t> lines_set;
      std::vector<old_generic_display_point_set_t> points_set;
      std::vector<arrow_t> arrows;
      std::vector<sphere_t> spheres;
      std::vector<torus_t> tori;
      std::vector<arc_t> arcs;
      std::vector<dodec_t> dodecs;
      std::vector<pentakis_dodec_t> pentakis_dodecs;
      std::vector<int> GL_display_list_handles;
      explicit old_generic_display_object_t(const std::string &n) : name(n) {
         is_displayed_flag = false;
         is_closed_flag = false;
         opacity = 1.0;
         is_transparent_flag = false;
         imol = UNDEFINED;
         is_solid_flag = false;
      }
      old_generic_display_object_t(int imol_in, const std::string &n) : name(n) {
         is_displayed_flag = false;
         is_closed_flag = false;
         opacity = 1.0;
         is_transparent_flag = false;
         imol = imol_in;
         is_solid_flag = false;
      }
      void add_line(const colour_holder &colour_in,
                    const std::string &colour_name,
                    const int &width_in,
                    const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords_in);
      void add_point(const colour_holder &colour_in,
                     const std::string &colour_name,
                     const int &size_in,
                     const clipper::Coord_orth &coords_in);
      void add_dodecahedron(const colour_holder &colour_in,
                            const std::string &colour_name,
                            double radius, const clipper::Coord_orth &pos);
      void add_pentakis_dodecahedron(const colour_holder &colour_in,
                                     const std::string &colour_name,
                                     double stellation_factor,
                                     double radius,
                                     const clipper::Coord_orth &pos);
      void raster3d(std::ofstream &render_stream) const;
      void clear() {
         lines_set.clear();
         points_set.clear();
      }
      void close_yourself() {
         name = "Closed Generic Display Object";
         clear();
         is_closed_flag = true;
      }
      static colour_holder colour_values_from_colour_name(const std::string &colour_name);
      int get_imol() const { return imol; }
      bool is_valid_imol() { return (imol != UNDEFINED && imol != INTERMEDIATE_ATOMS); }
      bool is_intermediate_atoms_object() const { return (imol == INTERMEDIATE_ATOMS); }
      void attach_to_molecule(int imol_in) { imol = imol_in; }
      void attach_to_intermediate_atoms() { imol = INTERMEDIATE_ATOMS; }
   };

   class old_generic_text_object_t {
   public:
      int handle;
      std::string s;
      float x, y, z;
      old_generic_text_object_t(const std::string &s_in, int handle_in,
                                float x_in, float y_in, float z_in) : s(s_in) {
         handle = handle_in;
         x = x_in;
         y = y_in;
         z = z_in;
      }
   };

}

#endif // HAVE_OLD_GENERIC_DISPLAY_OBJECT_HH
