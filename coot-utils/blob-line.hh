#ifndef BLOB_LINE_HH
#define BLOB_LINE_HH

#include <utility>
#include <clipper/core/xmap.h>

namespace coot {

std::pair<bool, clipper::Coord_orth>
find_peak_along_line_favour_front(const clipper::Coord_orth &p1,
                                  const clipper::Coord_orth &p2,
                                  float contour_level,
                                  const clipper::Xmap<float> &xmap);

}


#endif // BLOB_LINE_HH
