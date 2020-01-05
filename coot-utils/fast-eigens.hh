
#include <tuple>
#include <clipper/core/map_interp.h>

namespace coot {

   // maybe should be called "sloppy eigens"?

   std::tuple<double, double, double> fast_eigens(clipper::Matrix<double> &m, bool sort_eigenvalues);

}
