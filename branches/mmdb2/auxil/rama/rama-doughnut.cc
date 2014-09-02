
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

   double r1 = 13;
   double r2 = 20;

   std::cout << "(new-generic-object-number \"doughnut\")" << std::endl;
   double step = 1.5;
   short int even_odd = 1; 
   for (double phi=0.0; phi<360.0; phi += step) { 
      for (double psi=0.0; psi<360.0; psi += step) {
	 p = r_non_gly_pro.probability(clipper::Util::d2rad(phi),
				       clipper::Util::d2rad(psi));
	 std::cout << "; " << phi << " " << psi << " " << 1000*p << "\n";
	 std::string colour = "white";

	 if (1) { 

	    // The position on a doughnut, given phi, psi and r1 and r2 is:

	    // First, put it at the origin and rotate about the x axis:
	 
	    double xringbase = 0;
	    double yringbase = (r1 + 4.0 * p) * sin(clipper::Util::d2rad(180 + phi));
	    double zringbase = (r1 + 4.0 * p) * cos(clipper::Util::d2rad(180 + phi));

	    // Now push the ring up the y axis and rotate about the z
	    // axis by psi:
	    // 
	    double x = (xringbase * sin(clipper::Util::d2rad(psi)) + (yringbase + r2) * cos(clipper::Util::d2rad(psi)));
	    double y = (xringbase * cos(clipper::Util::d2rad(psi)) - (yringbase + r2) * sin(clipper::Util::d2rad(psi)));
	    double z = zringbase;
	 
	    // 	 std::cout << "(to-generic-object-add-point 0 \""
	    //  		   << colour << "\" 2 " << xringbase << " " << yringbase << " " << zringbase
	    //  		   << ")\n";

	    double pstep = 0.06;
	    if (p > pstep)
	       colour = "blue";
//  	    if (p > 2.0 * pstep)
//  	       colour = "sky";
	    if (p > 3.0 * pstep)
	       colour = "cyan";
	    if (p > 4.0 * pstep)
	       colour = "sea";
	    if (p > 5.0 * pstep)
	       colour = "aquamarine";
	    if (p > 6.0 * pstep)
	       colour = "green";
// 	    if (p > 7.0 * pstep)
// 	       colour = "greentint";
	    if (p > 8.0 * pstep)
	       colour = "forestgreen";
	    if (p > 9.0 * pstep)
	       colour = "yellowgreen";
	    if (p > 10.0 * pstep)
	       colour = "yellow";
	    if (p > 11.0 * pstep)
	       colour = "goldenrod";
	    if (p > 12.0 * pstep)
	       colour = "orange";
	    if (p > 13.0 * pstep)
	       colour = "orangered";
	    if (p > 14.0 * pstep)
	       colour = "red";
	    std::cout << "(to-generic-object-add-point 0 \""
		      << colour << "\" 5 " << x << " " << y << " " << z
		      << ")\n";
	 }
	 even_odd = 1 - even_odd;

      }
      std::cout << "\n";
   }
   std::cout << "(set-display-generic-object 0 1)" << std::endl;
   
   return 0; 
} 
