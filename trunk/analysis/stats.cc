
#include <algorithm> // for sorting

#include "stats.hh"

double
coot::stats::pnorm::erf(const double &z_in) const {

   double z = z_in;

   // https://oeis.org/A007680
   unsigned int n_tsec = 21;
   double tsec[] = { 1, 3, 10, 42, 216, 1320, 9360, 75600, 685440,
		     6894720, 76204800, 918086400, 11975040000,
		     168129561600.0, 2528170444800.0,
		     40537905408000.0, 690452066304000.0,
		     12449059983360000.0, 236887827111936000.0,
		     4744158915944448000.0, 99748982335242240000.0 };

   double z_lim = 2.545; // Taylor series expansion fails beyond this (3.6/sqrt(2))
   if (fabs(z) > z_lim) { 
      if (z < -z_lim) { 
	 z = -z_lim;
      } else { 
	 z = z_lim;
      }
   }
    
   double ts = 0;
   for (unsigned int i=0; i<n_tsec; i++) { 
      double exp_v = float(i*2+1);
      double m = 1; 
      if (i%2!=0) 
	 m = -1;
      ts += m * pow(z,exp_v)/tsec[i];
   }
   return (2/sqrt(M_PI)) * ts;
}



// 20150807-PE
// 
// Are the data in v_in drawn from a normal distribution with mean
// normal_mean and variance normal_variance?
//
// The KS value returned is the maximum value of the distance between
// the cumulaltive distributions.  i.e. the maxium of the fraction of
// the distribution that has passed by the time we get to a particular
// value
// 
double
coot::stats::get_kolmogorov_smirnov_vs_normal(const std::vector<double> &v_in,
					      const double &reference_mean,
					      const double &reference_sd) {

   std::vector<double> v = v_in;
   std::sort(v.begin(), v.end());
   double rdn = 1.0/double(v.size()); // for reciprocals.

   coot::stats::pnorm pn;
   double max_diff = 0;
   for (unsigned int i=0; i<v.size(); i++) {
      double f = v[i];
      // how many standard deviations from the mean is f?
      double z = (f - reference_mean)/reference_sd;
      double p = pn.get(z);
      double dist_frac = double(i) * rdn;
      double diff = fabs(dist_frac - p);
      if (diff > max_diff)
	 max_diff = diff;
   }
   return max_diff;
} 
