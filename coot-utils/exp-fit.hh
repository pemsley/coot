#ifndef EXP_FIT_HH
#define EXP_FIT_HH

#include <cmath>
#include <vector>

namespace coot {

   class exponential_fit_with_offset {
   public:

      // for y = a + bexp(cx)
      double a; // the offset
      double b; // scale fact
      double c; // scale for x

      // A_data must be in increasing order of x (up to caller to sort that out)
      explicit exponential_fit_with_offset(const std::vector<std::pair<double, double> > &A_data);
      double at(const double &x) const {
         double v = a + b * std::exp(c * x);
         return v;
      }
      double average_deviation(const std::vector<std::pair<double, double> > &A_data) const {
         if (! A_data.empty()) {
            double sum_delta = 0.0;
            for (unsigned int i=0; i<A_data.size(); i++) {
               const double &x = A_data[i].first;
               const double &y = A_data[i].second;
               double y_hat = at(x);
               double delta = std::fabs(y_hat - y);
               sum_delta += delta;
            }
            double av_delta = sum_delta/static_cast<double>(A_data.size());
            return av_delta;
         } else {
            return 0.0;
         }
      }
   };
}

#endif // EXP_FIT_HH
