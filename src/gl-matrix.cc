
#include <iostream>
#include <math.h>

#include "gl-matrix.h"


GL_matrix::GL_matrix() {

   // we don't want to use an "if";
   //
   for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {

	 mat[i*4+j] = 0.0;
      }
   }
      
   for (int i=0; i<4; i++) {
      mat[i*5] = 1.0;
   }
}

GL_matrix::GL_matrix(const clipper::Mat33<double> &m) {

   // indexing bug fix from Ansgar Esztermann
   // 
   mat[ 0] = m(0,0);
   mat[ 1] = m(0,1);
   mat[ 2] = m(0,2);
   mat[ 3] = 0;
   mat[ 4] = m(1,0);
   mat[ 5] = m(1,1);
   mat[ 6] = m(1,2);
   mat[ 7] = 0;
   mat[ 8] = m(2,0);
   mat[ 9] = m(2,1);
   mat[10] = m(2,2);
   mat[11] = 0;
   mat[12] = 0; // don't skew
   mat[13] = 0;
   mat[14] = 0;
   mat[15] = 1;

}

// We want an orientation (doesn't matter which) that is along normal.
GL_matrix::GL_matrix(const clipper::Coord_orth &normal) {

   clipper::Coord_orth d_unit = normal;
   
   clipper::Coord_orth arb(0,0.1,0.9);
   if (d_unit.y() < d_unit.z())
      arb = clipper::Coord_orth(0.0, 0.9, 0.1);
   if (d_unit.x() < d_unit.y())
      arb = clipper::Coord_orth(0.9, 0.0, 0.1);
	    
   clipper::Coord_orth p1(clipper::Coord_orth::cross(arb, d_unit).unit());
   clipper::Coord_orth p2(clipper::Coord_orth::cross( p1, d_unit).unit());
   clipper::Coord_orth p3 = d_unit;
	    
   GL_matrix m(p1.x(), p1.y(), p1.z(),
	       p2.x(), p2.y(), p2.z(),
	       p3.x(), p3.y(), p3.z());

   mat[0] = p1.x(); mat[1] = p1.y(); mat[ 2] = p1.z();
   mat[4] = p2.x(); mat[5] = p2.y(); mat[ 6] = p2.z();
   mat[8] = p3.x(); mat[9] = p3.y(); mat[10] = p3.z();
   
   mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0.0;
   mat[15] = 1.0;
} 




GL_matrix::GL_matrix(float m11, float m12, float m13,
		     float m21, float m22, float m23,
		     float m31, float m32, float m33){

   //
   mat[0] = m11; mat[1] = m12; mat[ 2] = m13;
   mat[4] = m21; mat[5] = m22; mat[ 6] = m23;
   mat[8] = m31; mat[9] = m32; mat[10] = m33;

   mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0.0;
   mat[15] = 1.0;
}

std::ostream&
operator<<(std::ostream &s, const GL_matrix &m) {
   
   s <<    "(" << m.mat[ 0] << " " << m.mat[ 1] << " " << m.mat[ 2] << " " << m.mat[ 3] << ")\n";
   s <<    "(" << m.mat[ 4] << " " << m.mat[ 5] << " " << m.mat[ 6] << " " << m.mat[ 7] << ")\n";
   s <<    "(" << m.mat[ 8] << " " << m.mat[ 9] << " " << m.mat[10] << " " << m.mat[11] << ")\n";
   s <<    "(" << m.mat[12] << " " << m.mat[13] << " " << m.mat[14] << " " << m.mat[15] << ")\n";
   return s;
} 


void
GL_matrix::from_quaternion(float q[4]) { // quaternion q


   mat[0] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
   mat[1] = 2.0 * (q[0] * q[1] - q[2] * q[3]);
   mat[2] = 2.0 * (q[2] * q[0] + q[1] * q[3]);
   mat[3] = 0.0;
   
   mat[4] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
   mat[5]= 1.0 - 2.0 * (q[2] * q[2] + q[0] * q[0]);
   mat[6] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
   mat[7] = 0.0;
   
   mat[8] = 2.0 * (q[2] * q[0] - q[1] * q[3]);
   mat[9] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
   mat[10] = 1.0 - 2.0 * (q[1] * q[1] + q[0] * q[0]);
   mat[11] = 0.0;
   
   mat[12] = 0.0;
   mat[13] = 0.0;
   mat[14] = 0.0;
   mat[15] = 1.0;
}

GL_matrix::GL_matrix(float quat[4]) {
   from_quaternion(quat);
}

clipper::Mat33<double>
GL_matrix::to_clipper_mat() const {

   clipper::Mat33<double> m;

   m(0,0) = mat[0];
   m(0,1) = mat[1];
   m(0,2) = mat[2];

   m(1,0) = mat[4];
   m(1,1) = mat[5];
   m(1,2) = mat[6];

   m(2,0) = mat[8];
   m(2,1) = mat[9];
   m(2,2) = mat[10];

   return m;
}


void 
GL_matrix::rotate_X(float theta) {

   float cos_theta = cos(theta);
   float sin_theta = sin(theta);

   float t0, t1, t2, t3, u0, u1, u2, u3;

   t0 = mat[4]*cos_theta + mat[ 8]*sin_theta;
   t1 = mat[5]*cos_theta + mat[ 9]*sin_theta;
   t2 = mat[6]*cos_theta + mat[10]*sin_theta;
   t3 = mat[7]*cos_theta + mat[11]*sin_theta;

   u0 = -mat[4]*sin_theta + mat[ 8]*cos_theta;
   u1 = -mat[5]*sin_theta + mat[ 9]*cos_theta;
   u2 = -mat[6]*sin_theta + mat[10]*cos_theta;
   u3 = -mat[7]*sin_theta + mat[11]*cos_theta;

   mat[4] = t0;
   mat[5] = t1;
   mat[6] = t2;
   mat[7] = t3;

   mat[ 8] = u0;
   mat[ 9] = u1;
   mat[10] = u2;
   mat[11] = u3;
      
}

void
GL_matrix::rotate_Z(float angle) {


}

void
GL_matrix::rotate_Y(float theta) {

   // x' = x*cos(theta) - y*sin(theta);
   // y' = x*sin(theta) + y*cos(theta);
   
   float cos_theta = cos(theta);
   float sin_theta = sin(theta);

   float t0, t1, t2, t3, u0, u1, u2, u3;

   t0 = mat[0]*cos_theta - mat[ 8]*sin_theta;
   t1 = mat[1]*cos_theta - mat[ 9]*sin_theta;
   t2 = mat[2]*cos_theta - mat[10]*sin_theta;
   t3 = mat[3]*cos_theta - mat[11]*sin_theta;

   u0 = mat[0]*sin_theta + mat[ 8]*cos_theta; 
   u1 = mat[1]*sin_theta + mat[ 9]*cos_theta; 
   u2 = mat[2]*sin_theta + mat[10]*cos_theta; 
   u3 = mat[3]*sin_theta + mat[11]*cos_theta; 
   
   mat[0] = t0;
   mat[1] = t1;
   mat[2] = t2;
   mat[3] = t3;

   mat[ 8] = u0;
   mat[ 9] = u1;
   mat[10] = u2;
   mat[11] = u3;

}

const float*
GL_matrix::get() const {

   return mat;
}


#ifndef HAVE_GSL
// Given a symmetric positive definate matrix,
// return the lower triangular matrix by LU decomposition.
//
std::pair<bool, GL_matrix>
GL_matrix::cholesky() const {

   //
   GL_matrix l; // identity

   // I am not sure about this indexing - it may be the other
   // way round... yes, it looks like it is.. (eff). 
   // 
   // l_{11} = mat[0]   l_{12} = mat[1]   l_{13} = mat[2]
   // l_{21} = mat[4]   l_{22} = mat[5]   l_{23} = mat[6]
   // l_{31} = mat[8]   l_{32} = mat[9]   l_{33} = mat[10]
   //

   l.mat[0] = cholesky_diag    (l, 1);
   l.mat[4] = cholesky_non_diag(l, 2, 1);
   l.mat[8] = cholesky_non_diag(l, 3, 1);
   l.mat[1] = 0; // cholesky_non_diag(l, 1,2);
   l.mat[5] = cholesky_diag(l, 2);
   l.mat[9] = cholesky_non_diag(l, 3, 2);
   l.mat[2] = 0; // cholesky_non_diag(l, 1, 3);
   l.mat[6] = 0; // cholesky_non_diag(l, 2, 3);
   l.mat[10]= cholesky_diag(l, 3);

//    l.mat[0] = cholesky_diag    (l, 1);
//    l.mat[1] = cholesky_non_diag(l, 1,2);
//    l.mat[2] = cholesky_non_diag(l, 1,3);

//    l.mat[4] = cholesky_non_diag(l, 2,1);
//    l.mat[5] = cholesky_diag(l, 2);
//    l.mat[6] = cholesky_non_diag(l, 2,3);

//    l.mat[8] = cholesky_non_diag(l, 3,1);
//    l.mat[9] = cholesky_non_diag(l, 3,2);
//    l.mat[10]= cholesky_diag(l, 3);


   // fix up the other half of the triangular matrix.

   // and the 4th column:
   l.mat[ 3] = 0;
   l.mat[ 7] = 0;
   l.mat[11] = 0;
   l.mat[15] = 1; 

   return std::pair<bool,GL_matrix> (1, l); 
}

// Use the GSL for cholesky then:
#else 
std::pair<bool,GL_matrix>
GL_matrix::cholesky() const {

   double a_data[] = { mat[0], mat[1], mat[ 2],
		       mat[4], mat[5], mat[ 6], 
		       mat[8], mat[9], mat[10] }; 

   gsl_matrix_view m = gsl_matrix_view_array (a_data, 3, 3);

   // test the princal minors for being positive (required for
   // m.matrix being positive definite (which in turn is required
   // for m.matrix not to colapse gsl_linalg_cholesky_decomp in a
   // big heap)).

   float pm1 = (mat[5] * mat[10]) - (mat[6] * mat[ 9]);
   float pm2 = (mat[0] * mat[ 5]) - (mat[1] * mat[ 4]);
   float pm3 = (mat[0] * mat[10]) - (mat[2] * mat[ 8]);

   // std::cout << "PM: " << pm1 << " " << pm2 << " " << pm3 << std::endl;
   if ((pm1<0.0) || (pm2<0.0) || (pm3<0.0)) {
      // std::cout << "Ooops - negative principal minors!" << std::endl;

      return std::pair<bool, GL_matrix> (0, GL_matrix(gsl_matrix_get(&m.matrix, 0, 0),
						      gsl_matrix_get(&m.matrix, 0, 1),
						      gsl_matrix_get(&m.matrix, 0, 2),
						      gsl_matrix_get(&m.matrix, 1, 0),
						      gsl_matrix_get(&m.matrix, 1, 1),
						      gsl_matrix_get(&m.matrix, 1, 2),
						      gsl_matrix_get(&m.matrix, 2, 0),
						      gsl_matrix_get(&m.matrix, 2, 1),
						      gsl_matrix_get(&m.matrix, 2, 2)));

   } else {
      gsl_error_handler_t *old_handler;
      old_handler = gsl_set_error_handler(my_aniso_error_handler);
      int ic = gsl_linalg_cholesky_decomp (&m.matrix);
      gsl_set_error_handler(old_handler);
      
      return std::pair<bool, GL_matrix> (1, GL_matrix(gsl_matrix_get(&m.matrix, 0, 0),
						      gsl_matrix_get(&m.matrix, 0, 1),
						      gsl_matrix_get(&m.matrix, 0, 2),
						      gsl_matrix_get(&m.matrix, 1, 0),
						      gsl_matrix_get(&m.matrix, 1, 1),
						      gsl_matrix_get(&m.matrix, 1, 2),
						      gsl_matrix_get(&m.matrix, 2, 0),
						      gsl_matrix_get(&m.matrix, 2, 1),
						      gsl_matrix_get(&m.matrix, 2, 2)));
   }
} 
#endif // GSL for Cholesky


#ifdef HAVE_GSL
void my_aniso_error_handler (const char * reason,
			     const char * file,
			     int line,
			     int gsl_errno) {
   std::cout << "Non-positive definite anisotropic atom!" << std::endl;
}
#endif // HAVE_GSL


// move these functions into the header file so that they
// are inlined (when they work).
// 
float
GL_matrix::cholesky_diag(const GL_matrix &l, int i) const {

   //
   // l_{ii} = sqrt(a_{ii} - \sum^{i-1}_{k=1}{l_{ik}^2})
   //
   
   float a_ii = mat[(i-1)*5];

   float sum = 0.0;

   for(int k=1; k<=(i-1); k++) {
      sum += squared(l.mat[(i-1)*4 + k - 1]); // fix the indexing?
      // sum += squared(l.mat[(k-1)*4+i-1]); // fix the indexing?
   }

   if ( (a_ii - sum) < 0) {
      std::cout << "WARNING negative sqrt in cholesky_diag(" << i << ")" << std::endl;
      std::cout << "a_ii: " << a_ii << ", sum: " << sum << std::endl;
   } 
   return sqrt(a_ii - sum); 
}


// quesion: what is the relationship of i, j with the indexing of
// mat.
//

float
GL_matrix::cholesky_non_diag(const GL_matrix &l, int j, int i) const {

   //
   float l_ii = l.mat[(i-1)*5];  // the divisor

   float a_ji = mat[(j-1)*4+i-1]; 

   float sum = 0.0;

   for(int k=1; k<=i-1; k++) {
      sum += l.mat[(j-1)*4+k-1] * l.mat[(i-1)*4+k-1];
      // sum += l.mat[(k-1)*4+j-1] * l.mat[(k-1)*4+i-1]; // fix the indexing
   }

   if ( (a_ji - sum) < 0) {
      std::cout << "WARNING negative numerator in cholesky_diag("
	   << i << "," << j << ")" << std::endl;
      std::cout << "a_ji: " << a_ji << ", sum: " << sum << std::endl;
   } 
   return (a_ji - sum)/l_ii; 

} 

GL_matrix
GL_matrix::transpose() const {

  //
  return GL_matrix( (*this), TRANSPOSE ); 

} 

GL_matrix::GL_matrix( GL_matrix in, int manip) {
  //
  if (manip == TRANSPOSE) {

    mat[0] = in.mat[0];
    mat[1] = in.mat[4]; 
    mat[2] = in.mat[8];
    
    mat[4] = in.mat[1]; 
    mat[5] = in.mat[5]; 
    mat[6] = in.mat[9]; 

    mat[8]  = in.mat[2]; 
    mat[9]  = in.mat[6]; 
    mat[10] = in.mat[10];

    mat[ 3] = mat[ 7] = mat[11] = 0.0; 
    mat[12] = mat[13] = mat[14] = 0.0;
    mat[15] = 1.0;
  }
  
}

// premultipy current matrix by matrix in.
// return a GL_matrix (don't modify this GL_matrix)
// 
GL_matrix
GL_matrix::mat_mult(GL_matrix in) const {

  GL_matrix res;

  // top row
  res.mat[0] = in.mat[0]*mat[0] + in.mat[1]*mat[4] + in.mat[2]*mat[ 8];
  res.mat[1] = in.mat[0]*mat[1] + in.mat[1]*mat[5] + in.mat[2]*mat[ 9];
  res.mat[2] = in.mat[0]*mat[2] + in.mat[1]*mat[6] + in.mat[2]*mat[10];

  // middle row
  res.mat[4] = in.mat[4]*mat[0] + in.mat[5]*mat[4] + in.mat[6]*mat[ 8];
  res.mat[5] = in.mat[4]*mat[1] + in.mat[5]*mat[5] + in.mat[6]*mat[ 9];
  res.mat[6] = in.mat[4]*mat[2] + in.mat[5]*mat[6] + in.mat[6]*mat[10];

  // 3rd row
  res.mat[ 8] = in.mat[8]*mat[0] + in.mat[9]*mat[4] + in.mat[10]*mat[ 8];
  res.mat[ 9] = in.mat[8]*mat[1] + in.mat[9]*mat[5] + in.mat[10]*mat[ 9];
  res.mat[10] = in.mat[8]*mat[2] + in.mat[9]*mat[6] + in.mat[10]*mat[10];

  std::cout << "   check: "
       << in.mat[0] << "*" << mat[0] << " + "
       << in.mat[1] << "*" << mat[4] << " + "
       << in.mat[2] << "*" << mat[8] << " = " << res.mat[0] << std::endl;
  
  return res; 
} 



coot::Cartesian
GL_matrix::mult(const coot::Cartesian &in) const {

   return coot::Cartesian(mat[0]*in.x() + mat[1]*in.y() + mat[ 2]*in.z(),
			  mat[4]*in.x() + mat[5]*in.y() + mat[ 6]*in.z(),
			  mat[8]*in.x() + mat[9]*in.y() + mat[10]*in.z());
}


void 
GL_matrix::print_matrix() const { 

   std::cout << mat[0] << "  " << mat[1] << "  " << mat[2] << "  \n"
	     << mat[4] << "  " << mat[5] << "  " << mat[6] << "  \n"
	     << mat[8] << "  " << mat[9] << "  " << mat[10] << "\n";

}

// irow and icol: 0 -> 2 inclusive
float
GL_matrix::matrix_element(int icol, int irow) const {

   float f = 0;
   if (icol<3 && icol>=0 && irow<3 && irow>=0) { 
      return mat[4*irow + icol];
   } else {
      std::cout << "ERROR in GL_matrix index " << icol << " " << irow << std::endl;
      return f;
   }
   
}
