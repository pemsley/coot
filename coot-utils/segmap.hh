
#include "clipper/core/coords.h"
#include "clipper/core/xmap.h"
#include "clipper/ccp4/ccp4_map_io.h"

// Another go at segmenting maps with symmetry considerations

namespace coot {

   // I want to give the symmetry as a n-fold a vector and centre
   //
   // or: as n-fold and presume that the axis is the z axis and goes through the middle
   // 
   // or: as n-fold and presume that the axis is the z axis and goes through the origin

   class segmap {
      std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > > find_peaks(float c) const;
      clipper::Xmap<float> flood_from_peaks(const std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > > &peaks,
					    float cut_off_for_flooding);
   public:
      segmap(const clipper::Xmap<float> &xmap_in) : xmap(xmap_in) {}
      const clipper::Xmap<float> &xmap;
      void proc(bool do_write_flag, const std::string &file_name);
      // remove "dust" that is smaller than 5.0A across (by default)
      void dedust(float vol=5.0);

   };

}
