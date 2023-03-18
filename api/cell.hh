#ifndef COOT_API_CELL_HH
#define COOT_API_CELL_HH

namespace coot {
   //! Simple Cell class for transfering symmetry info.
   class Cell {
   public:
      //! the unit cell lengths and angles in radians.
      float a, b, c, alpha, beta, gamma;
      Cell(float a, float b, float c, float alpha, float beta, float gamma) :
         a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma) {}
      Cell() { a = -1; b = -1; c = -1; alpha = -1; beta = -1; gamma = -1; }
   };

}


#endif // COOT_API_CELL_HH
