/*
 * src/meshed-generic-display-object.hh
 *
 * Copyright 2020 by Medical Research Council
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
 *
 */

#ifndef MESHED_GENERIC_DISPLAY_OBJECT_HH
#define MESHED_GENERIC_DISPLAY_OBJECT_HH

#include "utils/dodec.hh"

#include "coords/graphical-bonds-container.hh"

// How many ways of specifying a colour does an application need?
// At least 3, it turns out.

#include "generic-vertex.hh"
#include "coot-utils/g_triangle.hh"
#include "Material.hh"

#include "Mesh.hh"
#include "utils/colour-holder.hh"
#include "coot-utils/arc-info.hh"

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name);

coot::colour_holder colour_values_from_colour_name(const std::string &colour_name);

// useful utils
glm::vec4 colour_holder_to_glm(const coot::colour_holder &ch);
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
   class cone_t {
   public:
      cone_t() { radius = 0.15; }
      cone_t(const clipper::Coord_orth &pt1, const clipper::Coord_orth &pt2, float r = 0.15f) :
         start_point(pt1), end_point(pt2) {
         radius = r;
      }
      clipper::Coord_orth start_point;
      clipper::Coord_orth end_point;
      coot::colour_holder col;
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

   class object_info_t {
   public:
      coot::colour_holder colour;
      clipper::Coord_orth position;
   };

   // for multi-add points
   class point_info_t {
   public:
      point_info_t(const coot::colour_holder &ch, const clipper::Coord_orth &p, int w) : colour(ch), position(p), width(w) {}
      coot::colour_holder colour;
      clipper::Coord_orth position;
      int width;
   };

   // for multi-add lines
   class line_info_t {
   public:
      line_info_t(const coot::colour_holder &ch, const clipper::Coord_orth &p1, const clipper::Coord_orth &p2, float w) :
         colour(ch), position_start(p1), position_end(p2), radius(w) {}
      coot::colour_holder colour;
      clipper::Coord_orth position_start;
      clipper::Coord_orth position_end;
      float radius;
   };

   enum {UNDEFINED = -1, INTERMEDIATE_ATOMS=-9};
   meshed_generic_display_object() : mesh(Mesh("init_meshed_generic_display_object-A"))
      { imol = UNDEFINED; wireframe_mode = false; }
   explicit meshed_generic_display_object(const std::string &name_in) : mesh(Mesh(name_in)) {
      imol = UNDEFINED;
      mesh.name = name_in;
      wireframe_mode = false; }
   explicit meshed_generic_display_object(const Mesh &mesh_in) : mesh(mesh_in) {
      imol = UNDEFINED; wireframe_mode = false; }
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
   std::vector<object_info_t> info; // a place to store the positions and colours
                                    // so that they can be retrieved in a python function
                                    // get_generic_object_info()

   bool wireframe_mode;
   void attach_to_intermediate_atoms() { imol = INTERMEDIATE_ATOMS; }
   void attach_to_molecule(int imol_in) { imol = imol_in; }
   void clear() {
      mesh.clear();
   }
   void close_yourself() { clear();
      mesh.close();
   }
   void add_sphere(const sphere_t &sphere);

   void add_line(const coot::colour_holder &colour, const std::string &colour_name, int line_width,
                 const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords);
   void add_dashed_line(const coot::colour_holder &colour, const std::string &colour_name,
                        const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords,
                        const Material &material, float line_width_scale = 1.0, unsigned int n_segments = 5);
   void add_arrow(const arrow_t &arrow);

   // I need to say which cap type, flat or rounded.
   enum cap_type { FLAT_CAP, ROUNDED_CAP };
   void add_cylinder(const std::pair<glm::vec3, glm::vec3> &start_end,
                     const coot::colour_holder &col, float radius,
                     unsigned int n_slices,
                     bool cap_start, bool cap_end,
                     cap_type start_cap_type, cap_type end_cap_type, bool do_faces=false,
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
   void add_points(std::vector<point_info_t> &pos, unsigned int num_subdivisions);
   void add_lines(std::vector<line_info_t> &liv);
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

   void translate(const coot::Cartesian &t);

   void remove_last_object(); // remove from info vector and remove 182 triangles from the mesh (that's a bit of a hack)

   void raster3d(std::ofstream &render_stream) const;

};


#endif // MESHED_GENERIC_DISPLAY_OBJECT_HH
