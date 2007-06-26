
#ifndef HAVE_GL_MATRIX
#define HAVE_GL_MATRIX

enum { TRANSPOSE }; 

#ifdef HAVE_GSL
#include "gsl/gsl_linalg.h"
#endif

#include "Cartesian.h"

class GL_matrix { 

   float mat[16];

 public:

   GL_matrix();
   GL_matrix(float m11, float m12, float m13,
	     float m21, float m22, float m23,
	     float m31, float m32, float m33);
   
   GL_matrix( GL_matrix, int manip);  

   float* operator()(void) const;

   const float* get(void) const;

   void rotate_X(float angle);
   void rotate_Y(float angle);
   void rotate_Z(float angle);

   friend std::ostream& operator<<(std::ostream&, GL_matrix);

   // return a "this is a useful matrix" flag, so that we don't draw
   // elipsoids for atoms with non-positive definite U matrices.
   std::pair<bool,GL_matrix> cholesky() const; 

   GL_matrix transpose() const; 

   float cholesky_non_diag(const GL_matrix &l, int i, int j) const; 

   float cholesky_diag(const GL_matrix &l, int i) const; 


   float squared(float a) const { return a*a; }

   void from_quaternion(float quat[4]); // constructor really, if 
   // we could pass a quaternion, not just a float[4], I would fix
   // it now.

   GL_matrix mat_mult(GL_matrix in) const; 

   coot::Cartesian mult(const coot::Cartesian &in) const;

   void print_matrix() const;

   // irow and icol: 0 -> 2 inclusive
   float matrix_element(int icol, int irow) const;

};

#endif // HAVE_GL_MATRIX
