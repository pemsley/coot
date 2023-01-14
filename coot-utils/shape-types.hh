#ifndef SHAPE_TYPES_HH
#define SHAPE_TYPES_HH

#include <clipper/core/coords.h>
#include <glm/glm.hpp>
#include "utils/colour-holder.hh"
#include "arc-info.hh"

namespace shapes {

   class arrow_t {
   public:
      arrow_t() { fract_head_size = 0.3; radius = 0.15; cone_radius = 0.3f; }
      arrow_t(const clipper::Coord_orth &pt1, const clipper::Coord_orth &pt2) :
         start_point(pt1), end_point(pt2) {
         fract_head_size = 0.3;
         radius = 0.15;
         cone_radius = 0.3;
      }
      clipper::Coord_orth start_point;
      clipper::Coord_orth end_point;
      float fract_head_size;
      coot::colour_holder col; // use this type of colour because we use cylinder
      float radius;
      float cone_radius;
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

}

#endif // SHAPE_TYPES_HH
