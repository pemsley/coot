
#include <iostream>
#include "stats.hh"
#include "coot-utils/exp-fit.hh"

int main(int argc, char **argv) {

   int status = 0;

   if (true) {
      double A = 1.1;
      double B = 2.2;
      double C = 0.33;
      std::vector<std::pair<double, double> > data;
      unsigned int n = 20;
      for (unsigned int i=0; i<(n+1); i++) {
         double x = 0.1 * static_cast<double>(i);
         double y = A + B * std::exp(C * x);
         auto p = std::make_pair(x,y);
         data.push_back(p);
      }
      coot::exponential_fit_with_offset ef(data);
      for (unsigned int i=0; i<data.size(); i++) {
         const double &x = data[i].first;
         const double &y = data[i].second;
         double model_y = ef.at(x);
         std::cout << x << " " << y << " " << model_y << std::endl;
      }
   }

   if (false) {
      coot::stats::single s;
      s.add(1); s.add(2); s.add(3);
      s.add(4);
      s.add(5); s.add(6); s.add(7); s.add(8);
      s.add(9); s.add(10);

      std::pair<double, double> mi = s.median_and_iqr();

      std::cout << "median: " << mi.first << " irq " << mi.second << std::endl;
   }

   return status;
}
