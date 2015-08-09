
#ifndef INCLUDE_KOLMOGOROV_HH
#define INCLUDE_KOLMOGOROV_HH

#include <vector>

namespace nicholls {

   double get_KS(const std::vector<double> &v1, const std::vector<double> &v2);

   // Assumption - data are positive.
   // Kullback-Liebler 
   std::pair<double, double> get_KL(const std::vector<double> &v1, const std::vector<double> &v2);

}

#endif // INCLUDE_KOLMOGOROV_HH
