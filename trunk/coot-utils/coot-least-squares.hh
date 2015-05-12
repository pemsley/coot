
#include <vector>

namespace coot {

   class least_squares_fit { 
      double a_;
      double b_;
   public:
      least_squares_fit() {}
      least_squares_fit(const std::vector<std::pair<double, double> > &data);
      double m() const { return b_; }
      double c() const { return a_; }
   };

}
