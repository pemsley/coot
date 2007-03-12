// mode: -*-c++-*-

// Written by Kevin Cowtan 2002
// 
// (debugged by Paul Emsley shortly after :)

#include "math.h" // for floor()
#include <vector>
#include "clipper/core/clipper_util.h"

template <class T> class MatrixSq2
{ 
private:
   std::vector<T> data;
   int n_; 

public:
   MatrixSq2() { n_ = 0; }
   MatrixSq2(int n) { resize(n); }
   MatrixSq2(int n, T val) { resize(n, val); }
   
   void resize(int n) {n_=n; data.resize(size()); }
   void resize(int n, T val) {n_=n; data.resize(size(), val); }
   int size() const { return n_*n_*n_; } 
   int order() const { return n_; }

   //! interpolate
   T interp(double x1, double x2) const {

      int i1, i2, j1, j2; 
      double w1, w2, v1, v2; 
      
      w1 = floor ( x1 ); 
      w2 = floor ( x2 ); 
      i1 = int (w1); 
      i2 = int (w2); 
      j1 = i1+1; 
      j2 = i2+1; 
//       if (j1 >= n_) { 
// 	j1 -= -n_;
//       }
//       if (j2 >= n_) { 
// 	j2 -= n_;
//       }
      i1 = clipper::Util::mod(i1, n_); 
      i2 = clipper::Util::mod(i2, n_); 
      j1 = clipper::Util::mod(j1, n_); 
      j2 = clipper::Util::mod(j2, n_); 
      w1 = x1 - w1;  // e.g 0.9 (for 99.9)
      w2 = x2 - w2; 
      v1 = 1.0 - w1; // e.g. 0.1
      v2 = 1.0 - w2; 
      return ( v2 * ( v1*(*this)(i1,i2) + w1*(*this)(j1,i2) ) + 
	       w2 * ( v1*(*this)(i1,j2) + w1*(*this)(j1,j2)) ); 
   }

   // read accessor:
   const T& operator() (int i1, int i2) const { 
      return data[i1*n_ + i2]; }
   // write accessor:
   T& operator() (int i1, int i2) { 
      return data[i1*n_ + i2]; }

}; 


template<class T> class MatrixSq3
{
   
private:
   std::vector<T> data;
   int n_;

public:
   //! null constructor
   MatrixSq3() { n_ = 0;  }
   //! constructor
   MatrixSq3( const int& n ) { resize( n ); }
   //! constructor
   MatrixSq3( const int& n, T val ) { resize( n, val ); }
   //! resize
   void resize( const int& n )
   { n_ = n; data.resize( size() ); }
   //! resize
   void resize( const int& n, T val)
   { n_ = n; data.resize( size(), val ); }
   int size() const { return n_*n_*n_; }    //!< size
   const int& order() const { return n_; }  //!< order

    //! interpolate
    T interp( const double& x1, const double& x2, const double& x3 ) const
    {
      const MatrixSq3& mat = *this;
      int i1, i2, i3, j1, j2, j3;
      double w1, w2, w3, v1, v2, v3;

      w1 = floor( x1 );
      w2 = floor( x2 );
      w3 = floor( x3 );
      i1 = int (w1);
      i2 = int (w2);
      i3 = int (w3);
      j1 = i1+1;
      j2 = i2+1;
      j3 = i3+1;
      w1 = x1 - w1;
      w2 = x2 - w2;
      w3 = x3 - w3;
      v1 = 1.0 - w1;
      v2 = 1.0 - w2;
      v3 = 1.0 - w3;
      return
        ( v3*( v2*( v1*mat(i1,i2,i3) + w1*mat(j1,i2,i3) ) +
               w2*( v1*mat(i1,j2,i3) + w1*mat(j1,j2,i3) ) ) +
          w3*( v2*( v1*mat(i1,i2,j3) + w1*mat(j1,i2,j3) ) +
               w2*( v1*mat(i1,j2,j3) + w1*mat(j1,j2,j3) ) ) );
    }

   //! read accessor
   const T& operator () ( const int& i1, const int& i2, const int& i3 ) const
   { return data[ ( i1 * n_ + i2 ) * n_ + i3 ]; }
   //! write accessor
   T& operator () ( const int& i1, const int& i2, const int& i3 )
   { return data[ ( i1 * n_ + i2 ) * n_ + i3 ]; }

};



template<class T> class MatrixSq4
{

private:
   std::vector<T> data;
   int n_;

public:
   //! null constructor
   MatrixSq4() { n_ = 0;  }
   //! constructor
   MatrixSq4( const int& n ) { resize( n ); }
   //! constructor
   MatrixSq4( const int& n, T val ) { resize( n, val ); }
   //! resize
   void resize( const int& n )
   { n_ = n; data.resize( size() ); }
   //! resize
   void resize( const int& n, T val )
   { n_ = n; data.resize( size(), val ); }
   int size() const { return n_*n_*n_*n_; }  //!< size
   const int& order() const { return n_; }   //!< order
   //! interpolate
   T interp( const double& x1, const double& x2, 
	     const double& x3, const double& x4 ) const
   {
      const MatrixSq4& mat = *this;
      int i1, i2, i3, i4, j1, j2, j3, j4;
      double w1, w2, w3, w4, v1, v2, v3, v4;

      w1 = floor( x1 );
      w2 = floor( x2 );
      w3 = floor( x3 );
      w4 = floor( x4 );
      i1 = int (w1);
      i2 = int (w2);
      i3 = int (w3);
      i4 = int (w4);
      j1 = i1+1;
      j2 = i2+1;
      j3 = i3+1;
      j4 = i4+1;
      w1 = x1 - w1;
      w2 = x2 - w2;
      w3 = x3 - w3;
      w4 = x4 - w4;
      v1 = 1.0 - w1;
      v2 = 1.0 - w2;
      v3 = 1.0 - w3;
      v4 = 1.0 - w4;
      return
	 ( v4*( v3*( v2*( v1*mat(i1,i2,i3,i4) + w1*mat(j1,i2,i3,i4) ) +
		     w2*( v1*mat(i1,j2,i3,i4) + w1*mat(j1,j2,i3,i4) ) ) +
		w3*( v2*( v1*mat(i1,i2,j3,i4) + w1*mat(j1,i2,j3,i4) ) +
		     w2*( v1*mat(i1,j2,j3,i4) + w1*mat(j1,j2,j3,i4) ) ) ) +
	   w4*( v3*( v2*( v1*mat(i1,i2,i3,j4) + w1*mat(j1,i2,i3,j4) ) +
		     w2*( v1*mat(i1,j2,i3,j4) + w1*mat(j1,j2,i3,j4) ) ) +
		w3*( v2*( v1*mat(i1,i2,j3,j4) + w1*mat(j1,i2,j3,j4) ) +
		     w2*( v1*mat(i1,j2,j3,j4) + w1*mat(j1,j2,j3,j4) ) ) ) );
   }
   //! read accessor
   const T& operator () ( const int& i1, const int& i2, const int& i3, const int& i4 ) const
   { return data[ ( ( i1 * n_ + i2 ) * n_ + i3 ) * n_ + i4 ]; }
   //! write accessor
   T& operator () ( const int& i1, const int& i2, const int& i3, const int& i4  )
   { return data[ ( ( i1 * n_ + i2 ) * n_ + i3 ) * n_ + i4 ]; }
};

