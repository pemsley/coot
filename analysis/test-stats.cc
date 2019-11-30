
#include <iostream>
#include "stats.hh"

int main(int argc, char **argv) {

   int status = 0;

   coot::stats::single s;
   s.add(1); s.add(2); s.add(3);
   s.add(4);
   s.add(5); s.add(6); s.add(7); s.add(8);
   s.add(9); s.add(10);

   std::pair<double, double> mi = s.median_and_iqr();

   std::cout << "median: " << mi.first << " irq " << mi.second << std::endl;

   return status;
}
