
#include <clipper/ccp4/ccp4_map_io.h>

class clipper_map_file_wrapper : public clipper::CCP4MAPfile {
public:
   clipper_map_file_wrapper() : clipper::CCP4MAPfile() { }
   void wrap_open_read(const clipper::String &filename_in) {
      open_read(filename_in);
   }
   bool starts_at_zero() const {
      if (grid_map_.min() == clipper::Coord_grid(0,0,0)) {
	 return true;
      } else {
	 return false;
      }
   }
};

