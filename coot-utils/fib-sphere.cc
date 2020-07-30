
#include <math.h>
#include "fib-sphere.hh"

std::vector<clipper::Coord_orth>
coot::fibonacci_sphere(unsigned int n_samples) {

   // I watched the Coding Adventure on Boids
   double phi = M_PI * (3.0 - sqrt(5.0));

   std::vector<clipper::Coord_orth> points(n_samples);
   for (unsigned int i=0; i<n_samples; i++) {
      double d(n_samples - 1.0);
      double y = 1.0 - (static_cast<double>(i)/d) * 2.0;
      double radius = sqrt(1.0 - y * y);
      double theta = phi * static_cast<double>(i);
      double x = cos(theta) * radius;
      double z = sin(theta) * radius;
      points[i] = clipper::Coord_orth(x,y,z);
   }
   return points;
}

