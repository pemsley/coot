
#include <iostream>
#include "gl-matrix.h"

int main (int argc, char **argv) {


   std::vector<int> v;

   v.push_back(1);
   v.push_back(2);
   v.push_back(4);

   std::vector<int>::iterator it;
   for (it = v.begin(); it != v.end(); it++) {
      if ( *it == 2 )
	 v.erase(it);
   }
   for (unsigned int iv=0; iv<v.size(); iv++) {
      std::cout << "Test iterator " << iv << "   " << v[iv] << std::endl;
   }

   // ---------------------
   double m11 = 2;
   double m22 = 2;
   double m33 = 2;
   double m12 = 1;
   double m13 = 0;
   double m23 = 0;

   GL_matrix m(m11, m12, m13, m12, m22, m23, m13, m23, m33);

   std::pair<bool,GL_matrix> m2 = m.cholesky();

   std::cout << "Initial Matrix: " << std::endl;
   m.print_matrix();
   std::cout << "Cholesky decomposition of above: " << std::endl;
   m2.second.print_matrix();
   std::cout << "Multiplication simple of above: " << std::endl;
   GL_matrix sq1 = m2.second.mat_mult(m2.second);
   sq1.print_matrix();
   std::cout << "Multiplication of above using transpose: " << std::endl;
   GL_matrix sq2 = m2.second.mat_mult(m2.second.transpose());
   sq2.print_matrix();

   return 0;

}
