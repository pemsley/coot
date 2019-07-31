/* utils/dodec.cc
 *
 * Copyright 2013 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#include <fstream>
#include <algorithm> // for reverse()
#include "dodec.hh"

dodec::dodec() {

   double phi  = 0.5 * (1 + sqrt(5));
   double r3   = sqrt(3);
   double rr3  = 1.0/r3;
   double r3p  = 1.0/(r3 * phi);
   double pdr3 = phi/r3;

   clipper::Coord_orth p0( rr3,  rr3,  rr3);
   clipper::Coord_orth p1(-rr3,  rr3,  rr3);
   clipper::Coord_orth p2( rr3, -rr3,  rr3);
   clipper::Coord_orth p3( rr3,  rr3, -rr3);

   clipper::Coord_orth p4(-rr3, -rr3,  rr3);
   clipper::Coord_orth p5(-rr3,  rr3, -rr3);
   clipper::Coord_orth p6( rr3, -rr3, -rr3);
   clipper::Coord_orth p7(-rr3, -rr3, -rr3);

   clipper::Coord_orth p8 (  0,  r3p,  pdr3);
   clipper::Coord_orth p9 (  0, -r3p,  pdr3);
   clipper::Coord_orth p10(  0,  r3p, -pdr3);
   clipper::Coord_orth p11(  0, -r3p, -pdr3);

   clipper::Coord_orth p12( r3p,  pdr3,  0);
   clipper::Coord_orth p13( r3p, -pdr3,  0);
   clipper::Coord_orth p14(-r3p,  pdr3,  0);
   clipper::Coord_orth p15(-r3p, -pdr3,  0);

   clipper::Coord_orth p16( pdr3, 0,  r3p);
   clipper::Coord_orth p17(-pdr3, 0,  r3p);
   clipper::Coord_orth p18( pdr3, 0, -r3p);
   clipper::Coord_orth p19(-pdr3, 0, -r3p);

   points.push_back(p0 ); points.push_back(p1 ); points.push_back(p2 );
   points.push_back(p3 ); points.push_back(p4 ); points.push_back(p5 );
   points.push_back(p6 ); points.push_back(p7 ); points.push_back(p8 );
   points.push_back(p9 ); points.push_back(p10); points.push_back(p11);
   points.push_back(p12); points.push_back(p13); points.push_back(p14);
   points.push_back(p15); points.push_back(p16); points.push_back(p17);
   points.push_back(p18); points.push_back(p19);

   assign_face_rings();

}

void
dodec::test(const std::string &file_name) const {

   std::ofstream f(file_name.c_str());

   if (f) {
      for (unsigned int i=0; i<points.size(); i++) {
	 f << "  "
	   << points[i].x() << " "
	   << points[i].y() << " "
	   << points[i].z() << "\n";
      }
   }
   f.close();
}

std::vector<unsigned int>
dodec::face(unsigned int face_number) const {
   return rings[face_number];
}


void
dodec::assign_face_rings() {

   std::vector<std::pair<unsigned int, unsigned int> > pairs;

   for (unsigned int i=0; i<points.size(); i++) {
      for (unsigned int j=0; j<points.size(); j++) {
	 if (i != j) {
	    double d = (points[i]-points[j]).lengthsq();
	    if (d < 0.51) {
	       // std::cout << "  " << i << " " << j << " " << d << std::endl;
	       std::pair<unsigned int, unsigned int> p(i,j);
	       pairs.push_back(p);
	    }
	 }
      }
   }

   for (unsigned int i1=0; i1<pairs.size(); i1++) {
      for (unsigned int i2=0; i2<pairs.size(); i2++) {
	 if (pairs[i2].first == pairs[i1].second) {
	    for (unsigned int i3=0; i3<pairs.size(); i3++) {
	       if (pairs[i3].first == pairs[i2].second) {
		  for (unsigned int i4=0; i4<pairs.size(); i4++) {
		     if (pairs[i4].first == pairs[i3].second) {
			for (unsigned int i5=0; i5<pairs.size(); i5++) {
			   if (pairs[i5].first == pairs[i4].second) {
			      if (pairs[i5].second == pairs[i1].first) {
				 std::vector<unsigned int> ring(5);
				 ring[0] = pairs[i1].first;
				 ring[1] = pairs[i2].first;
				 ring[2] = pairs[i3].first;
				 ring[3] = pairs[i4].first;
				 ring[4] = pairs[i5].first;
				 if (!known_ring(ring, rings)) {
				    // we need the rings to go round in a direction (relative to the origin) so that
				    // the normals are calculated pointing out consistently. These have been hand
				    // tested one by one - terrible hack!
				    if ((rings.size() == 1) || (rings.size() == 4) || (rings.size() == 5) ||
					(rings.size() == 6) || (rings.size() == 7) || (rings.size() == 10)) {
				       std::reverse(ring.begin(), ring.end());
				       rings.push_back(ring);
				    } else {
				       // don't be perverse
				       rings.push_back(ring);
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (false) {
      for (unsigned int ii=0; ii<rings.size(); ii++) {
	 std::cout << "ring: ";
	 for (unsigned int jj=0; jj<rings[ii].size(); jj++) {
	    std::cout << " " << rings[ii][jj];
	 }
	 std::cout << std::endl;
      }
   }
}

bool
dodec::known_ring(std::vector<unsigned int> &ring,
		    const std::vector<std::vector<unsigned int> > &rings) const {

   bool r = false;
   for (unsigned int j=0; j<rings.size(); j++) {
      unsigned int n_match = 0;
      // are the all of the contents of ring in any ring in rings?
      for (unsigned int i=0; i<ring.size(); i++) {
	 for (unsigned int k=0; k<rings[j].size(); k++) {
	    if (ring[i] == rings[j][k])
	       n_match++;
	 }
      }
      if (n_match == 5) {
	 r = true;
	 break;
      }
   }
   return r;
}


pentakis_dodec::pentakis_dodec(double height_in) {

   prism_vertex_height = height_in;
   init();
}

void
pentakis_dodec::init() {


   // the radius of d is sqrt(3) perhaps?  So if height_in is that, the we get ??? (I don't know the name) "buckyball"

   std::vector<clipper::Coord_orth> coords = d.coords();
   pyrimid_vertices.resize(12);

   for (unsigned int i=0; i<12; i++) {

      // these are the indices of the vertices that make up the given face
      std::vector<unsigned int> v = d.face(i);

      // what is the coordinates of the point along the line from the
      // origin, through the middle of the plane and somewhat above
      // it?

      clipper::Coord_orth face_centre_sum(0,0,0);

      for (unsigned int j=0; j<5; j++) {
         face_centre_sum += coords[v[j]];
      }

      clipper::Coord_orth face_centre(0.2 * face_centre_sum);
      clipper::Coord_orth face_centre_unit(face_centre.unit());

      clipper::Coord_orth pyrimid_vertex(prism_vertex_height * face_centre_unit);
      pyrimid_vertices[i] = pyrimid_vertex;
   }

}
