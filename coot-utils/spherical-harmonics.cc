
// #include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_legendre.h>
#include "spherical-harmonics.hh"


void
coot::spherical_harmonics(clipper::NXmap<float> &nxmap) {

   /*  -- GSL 2 or so?
   const gsl_sf_legendre_t norm = GSL_SF_LEGENDRE_SCHMIDT;
   const size_t lmax = 3;

   double x = 0.5;

   size_t dim = gsl_sf_legendre_array_n(lmax);
   double *p = malloc(sizeof(double) * dim);

   gsl_sf_legendre_array(norm, lmax, x, p);
   */


   gsl_sf_result r;
   int i = gsl_sf_legendre_P1_e(-0.5, &r);
   std::cout << "gsl_sf " << i << " " << r.val << " " << r.err << std::endl;

}

