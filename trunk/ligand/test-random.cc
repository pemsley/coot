

#include <iostream>
#include "wligand.hh"

int main(int argc, char **argv) {

   coot::wligand w;
   float v;

   int n = 1000000;

   // Let's say from -5 to +5
   //
   // 10 bins per 1 (sigma)
   // 


   // There is a problem with binning here. The int(v*10) is giving
   // bad numbers at v close to 0 (|v|<0.1).  Actually, for negative
   // v.  We need floor() or something.  Yes.  That's it.  Problem
   // solved.
   
   std::vector<int> bins(100, 0);
   int istat = 0;

//    std::cout << "testint " << int(0.1) << std::endl;
//    std::cout << "testint " << int(-0.1) << std::endl; // Hmmm...
//    std::cout << "test floor " << int(floor(-0.1)) << std::endl;

   for (int i=0; i<n; i++) { 
      v = w.get_random_normal_value();
      // std::cout << "iv= " << i << " " << v << std::endl; // yes,
                                          // these look normal, mean 0, std dev 1.
      int ibin = int(floor(v*10)) + 50;
      if (ibin >= bins.size()) { 
	 std::cout << "error in binning! " << v << " " << ibin << " >= "
		   << bins.size() << std::endl;
	 istat = 1; // error
      } else { 
	 bins[ibin]++;
      }
   }

   for (int i=0; i<bins.size(); i++) {
      std::cout << i-50 << " " << bins[i] << "\n";
   }
   return istat;
}
