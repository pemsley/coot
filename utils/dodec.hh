
#ifndef DODEC_HH
#define DODEC_HH

#include <vector>
#include <string>
#include "clipper/core/coords.h"

class dodec {
   std::vector<clipper::Coord_orth> points;
   std::vector<std::vector<unsigned int> > rings; // point indices
   bool known_ring(std::vector<unsigned int> &ring,
		   const std::vector<std::vector<unsigned int> > &rings) const;
   void assign_face_rings();
public:
   dodec();
   std::vector<clipper::Coord_orth> coords() const { return points; }
   // ordered set of points
   std::vector<unsigned int> face(unsigned int face_number) const; // only call with 0->11 inclusive
   void test(const std::string &file_name) const;
};


#endif // DODEC_HH

