#ifndef API_CELL_HH
#define API_CELL_HH

namespace api {

   //! POD cell type for moorhen.
   //!
   //! the angles are in radians.
   class cell_t {
   public:
      //! a
      float a;
      //! b
      float b;
      //! c
      float c;
      //! alpha
      float alpha;
      //! beta
      float beta;
      //! gamma
      float gamma;
      //! is_set
      bool is_set;
      cell_t() { is_set = false; }
      cell_t(float a, float b, float c, float alpha, float beta, float gamma) :
         a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma), is_set(true) {}
   };

}

#endif // API_CELL_HH
