
// header-here

#include <fstream>
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
