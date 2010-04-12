
#include <iostream>
#include "lig-build.hh"

bool close_double(double a, double b) {

   double f = fabs(a-b);
   return (f<0.001);
}

int main(int argc, char **argv) {

   lig_build::atom_position_t p1(0,2);
   lig_build::atom_position_t p2(2,0);
   lig_build::atom_position_t p3(sqrt(3),1);
   lig_build::atom_position_t p4(1,sqrt(3));
   lig_build::atom_position_t p5(1,-2);
   lig_build::atom_position_t p6(-1,-2);

   lig_build::atom_position_t u = p1.unit_vector();

   if (close_double(u.length(), 1.0)) {
      std::cout << "PASS  unit length" << std::endl;
   } else { 
      std::cout << "FAIL  unit length" << std::endl;
   }

   double length = p2.length();
   if (close_double(length, 2.0)) {
      std::cout << "PASS  length 1 " << std::endl;
   } else { 
      std::cout << "FAIL  length 1 " << length << std::endl;
   }

   length = p3.length();
   if (close_double(length, 2)) {
      std::cout << "PASS  length 2 " << std::endl;
   } else { 
      std::cout << "FAIL  length 2 " << length << std::endl;
   }

   double angle = lig_build::atom_position_t::angle(p1, p2);
   if (close_double(angle, 90.0)) {
      std::cout << "PASS  angle " << std::endl;
   } else { 
      std::cout << "FAIL  angle " << angle << std::endl;
   }

   angle = p1.axis_orientation();
   if (close_double(angle, 90.0)) {
      std::cout << "PASS  axis_orientation 1 " << std::endl;
   } else { 
      std::cout << "FAIL axis_orientation 1 " << angle << std::endl;
   }

   angle = p2.axis_orientation();
   if (close_double(angle, 0.0)) {
      std::cout << "PASS  axis_orientation 2 " << std::endl;
   } else { 
      std::cout << "FAIL axis_orientation 2 " << angle << std::endl;
   }
   angle = p3.axis_orientation();
   if (close_double(angle, 30.0)) {
      std::cout << "PASS  axis_orientation 3 " << std::endl;
   } else { 
      std::cout << "FAIL axis_orientation 3 " << angle << std::endl;
   }
   angle = p4.axis_orientation();
   if (close_double(angle, 60.0)) {
      std::cout << "PASS  axis_orientation 4 " << std::endl;
   } else { 
      std::cout << "FAIL axis_orientation 4 " << angle << std::endl;
   }

   lig_build::atom_position_t r = p5.rotate(90);
   if (! close_double(r.x, 2))
      std::cout << "FAIL  rotate - a "  << r << std::endl;
   else 
      if (! close_double(r.y, 1))
	 std::cout << "FAIL  rotate - b "  << r << std::endl;
      else 
	 std::cout << "PASS  rotate "<< std::endl;

   angle = lig_build::atom_position_t::angle(r, p5);
   if (! close_double(angle, 90.0))
      std::cout << "FAIL  rotate angle "  << r << " " << angle << std::endl;
   else 
      std::cout << "PASS  rotate angle "<< std::endl;

   double t = atan(2)/DEG_TO_RAD; // ~63 degrees.
   angle = p5.axis_orientation();
   if (! close_double(angle, -t))
      std::cout << "FAIL  neg axis orientation 1 "  << 180+t << " "
		<< angle << std::endl;
   else 
      std::cout << "PASS  neg axis orientation 1 "<< std::endl;

   angle = p6.axis_orientation();
   if (! close_double(angle, t-180))
      std::cout << "FAIL  neg axis orientation 2 "  << t-180 << " "
		<< angle << std::endl;
   else 
      std::cout << "PASS  neg axis orientation 2 "<< std::endl;

   return 0;
}
