
#include <iostream>
#include "clipper/core/ramachandran.h"


int main(int argc, char **argv) {

   clipper::Ramachandran r_non_gly_pro;

   double rama_threshold_preferred = 0.02; 
   double rama_threshold_allowed = 0.0012;
   
   r_non_gly_pro.init(clipper::Ramachandran::NonGlyPro);
   r_non_gly_pro.set_thresholds(rama_threshold_preferred,rama_threshold_allowed);

   double phi = -57;
   double psi = -48;
   double p;

   double step = 5; 
   for (double phi=0.0; phi<360.0; phi += step) { 
      for (double psi=0.0; psi<360.0; psi += step) {
	 p = r_non_gly_pro.probability(clipper::Util::d2rad(phi),
				       clipper::Util::d2rad(psi));
	 std::cout << phi << " " << psi << " " << 10000*p << "\n";
	 // std::cout << 10000*p <<  " "; 
      }
      std::cout << "\n";
   }
   
   return 0; 
} 
