/*
 * utils/dodec.hh
 *
 * Copyright 2013 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


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
   const clipper::Coord_orth &get_point(const unsigned int &idx) const { return points[idx]; }
   void test(const std::string &file_name) const;
};

class pentakis_dodec {
   void init();
public:
   explicit pentakis_dodec(double height_in);
   pentakis_dodec() : prism_vertex_height(sqrt(3)) { init(); }
   dodec d;
   double prism_vertex_height; // for pentagonal prism from the origin.
   std::vector<clipper::Coord_orth> pyrimid_vertices;
   // for triangles, see make-a-dodec.cc
};


#endif // DODEC_HH

