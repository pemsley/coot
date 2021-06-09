#ifndef TUBE_FINDER_HH
#define TUBE_FINDER_HH

#include <clipper/core/xmap.h>
#include "gsl/gsl_multimin.h"

namespace coot {

   class tube_finder_t {

      class simplex_param_t {
      public:
         clipper::Coord_orth test_point_centre_orig;
         std::vector<clipper::Coord_orth> original_positions;
         const clipper::Xmap<float> *xmap;
      };
      
      // return the fitted rtop_orth
      std::pair<bool, clipper::Coord_orth> 
      fit_stick_in_map(const clipper::Coord_orth &test_point,
                       const clipper::Mat33<double> &orientation,
                       const clipper::Xmap<float> &xmap,
                       float density_level_crit) const;
      std::pair<bool, clipper::Coord_orth>
      fit_to_map_by_simplex_rigid(const clipper::Coord_orth &test_point,
                                  const clipper::Mat33<double> &orientation,
                                  const clipper::Xmap<float> &xmap,
                                  float density_level_crit) const;
      static double my_f_simplex_rigid_internal(const gsl_vector *v, void *params);
      clipper::Coord_orth apply_shifts_rigid_internal(gsl_vector *x, const simplex_param_t &par) const;
      static clipper::RTop_orth construct_matrix(const gsl_vector *v);
      // return negative on failure
      static float sphere_variance(const clipper::Coord_orth &centre_point,
                                   const std::vector<clipper::Coord_orth> &sphere_points,
                                   float radius,
                                   const clipper::Xmap<float> &xmap);

   public:
      explicit tube_finder_t(const clipper::Xmap<float> &xmap);
      std::vector<clipper::Coord_orth> get_positions() const;
      std::vector<clipper::Coord_orth> positions; // may be separate in the strands and helices

   };
}

#endif // TUBE_FINDER_HH
