
#ifndef MESHED_GENERIC_DISPLAY_OBJECT_HH
#define MESHED_GENERIC_DISPLAY_OBJECT_HH

#include "utils/coot-utils.hh"  // for colour_holder
#include "utils/dodec.hh"

#include "coords/graphical-bonds-container.hh"

// How many ways of specifying a colour does an application need?
// At least 3, it turns out.
#include "Mesh.hh"
#include "utils/colour-holder.hh"
#include "coot-utils/arc-info.hh"
#include "coot-colour.hh"

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name);

coot::colour_holder colour_values_from_colour_name(const std::string &colour_name);

glm::vec3 coord_orth_to_glm(const clipper::Coord_orth &co);

class meshed_generic_display_object {
public:
   class arrow_t {
   public:
      arrow_t() { fract_head_size = 0.3; radius = 0.15; }
      arrow_t(const clipper::Coord_orth &pt1, const clipper::Coord_orth &pt2) :
         start_point(pt1), end_point(pt2) {
         fract_head_size = 0.3;
         radius = 0.15;
      }
      clipper::Coord_orth start_point;
      clipper::Coord_orth end_point;
      float fract_head_size;
      coot::colour_holder col; // use this type of colour because we use cylinder
      float radius;
   };
   class sphere_t {
   public:
      sphere_t() {}
      sphere_t(const clipper::Coord_orth &centre_in, float r) : centre(centre_in) {
         radius = r;
      }
      clipper::Coord_orth centre;
      glm::vec4 col;
      float radius;
   };
   class torus_t {
   public:
      torus_t() {}
      torus_t(const clipper::Coord_orth &p,
              const clipper::Coord_orth &n,
              float r1, float r2) : position(p), normal(n), radius_1(r1), radius_2(r2) {
         n_ring_atoms = 6;
         height_scale = 1.0;
      }
      clipper::Coord_orth position;
      clipper::Coord_orth normal;
      coot::colour_holder col;
      float radius_1;
      float radius_2;
      float height_scale;
      int n_ring_atoms;
   };
   // arc is part of a torus
   class arc_t {
   public:
      arc_t(float delta_angle_in,
            const clipper::Coord_orth &start_point_in,
            const clipper::Coord_orth &start_dir_in,
            const clipper::Coord_orth &normal_in,
            float radius_in, float radius_inner_in) :
         normal(normal_in),
         start_point(start_point_in),
         start_dir(start_dir_in),
         delta_angle(delta_angle_in),
         radius(radius_in),
         radius_inner(radius_inner_in) {}
      arc_t(coot::arc_info_type &ai, float radius_in, float radius_inner_in,
            const coot::colour_holder &ch) :
         normal(ai.normal), start_point(ai.start_point), start_dir(ai.start_dir),
         orientation_matrix(ai.orientation_matrix),
         delta_angle(ai.delta),
         col(ch),
         radius(radius_in),
         radius_inner(radius_inner_in) {}
      clipper::Coord_orth normal;
      clipper::Coord_orth start_point;
      clipper::Coord_orth start_dir;
      clipper::Mat33<double> orientation_matrix;
      float delta_angle;
      coot::colour_holder col;
      float radius;
      float radius_inner;
   };
   class dodec_t {
   public:
      dodec_t(const dodec &d_in, double size_in, const clipper::Coord_orth &pos_in) :
         d(d_in), size(size_in), position(pos_in) { }
      dodec d;
      double size;
      clipper::Coord_orth position;
      coot::colour_holder col;
   };

   class pentakis_dodec_t { // perhaps this should inherit from above
   public:
      pentakis_dodec_t(const pentakis_dodec &pkdd_in, double size_in,
                       const clipper::Coord_orth &pos_in) : pkdd(pkdd_in), size(size_in), position(pos_in) { }
      pentakis_dodec pkdd;
      double size;
      clipper::Coord_orth position;
      coot::colour_holder col;
   };
   enum {UNDEFINED = -1, INTERMEDIATE_ATOMS=-9};
   meshed_generic_display_object() { imol = UNDEFINED; }
   explicit meshed_generic_display_object(const Mesh &mesh_in) : mesh(mesh_in) { imol = UNDEFINED; }
   std::map<unsigned int, std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > > origin_octasphere_map;
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
      wrapped_make_octasphere(unsigned int num_subdivisions, const glm::vec3 &position,
                              float radius, const glm::vec4 &col);
   void init(const graphical_bonds_container &bonds_box, bool background_is_black_flag);
   int imol;
   int get_imol() const { return imol; }
   bool is_valid_imol() { return imol != INTERMEDIATE_ATOMS && imol != UNDEFINED; }
   bool is_intermediate_atoms_object() const { return imol == INTERMEDIATE_ATOMS; }
   Mesh mesh;
   void attach_to_intermediate_atoms() { imol = INTERMEDIATE_ATOMS; }
   void attach_to_molecule(int imol_in) { imol = imol_in; }
   void clear() { mesh.clear(); imol = UNDEFINED; }
   void close_yourself() { clear(); mesh.close(); }
   void add(const sphere_t &sphere) {
      std::cout << "FIXME:: meshed-generic-display-object add a sphere here " << sphere.centre.format() << std::endl;
   }
   void add_line(const coot::colour_holder &colour, const std::string &colour_name, int line_width,
                 const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords);
   void add_arrow(const arrow_t &arrow);

   // I need to say which cap type, flat or rounded.
   enum cap_type { FLAT_CAP, ROUNDED_CAP };
   void add_cylinder(const std::pair<glm::vec3, glm::vec3> &start_end,
                     const coot::colour_holder &col, float radius,
                     unsigned int n_slices,
                     bool cap_start, bool cap_end,
                     cap_type start_cap_type, cap_type end_cap_type, bool do_faces=true,
                     float unstubby_cap_factor=1.0);
   void add_cone(const std::pair<glm::vec3, glm::vec3> &start_end,
                 const coot::colour_holder &col, float base_radius, float top_radius,
                 unsigned int n_slices,
                 bool cap_start, bool cap_end,
                 cap_type start_cap_type, cap_type end_cap_type);
   void add_point(const coot::colour_holder &colour_in,
                  const std::string &colour_name,
                  const int &size_in,
                  const clipper::Coord_orth &coords_in,
                  unsigned int num_subdivisions);
   void add_dodecahedron(const coot::colour_holder &colour_in,
                         const std::string &colour_name,
                         double radius, const clipper::Coord_orth &pos);
   void add_pentakis_dodecahedron(const coot::colour_holder &colour_in,
                                  const std::string &colour_name,
                                  double stellation_factor,
                                  double radius,
                                  const clipper::Coord_orth &pos);
   void add_arc(const arc_t &arc);
   void add_torus(const torus_t &torus);
   void raster3d(std::ofstream &render_stream) const;

};


#endif // MESHED_GENERIC_DISPLAY_OBJECT_HH
