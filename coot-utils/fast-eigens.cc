/*! \file fast-eigens.cc
    Source code file for faster less accurate eigen vectors
*/
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

#include "fast-eigens.hh"

std::tuple<double, double, double>
coot::fast_eigens(clipper::Matrix<double> &m, bool sort_eigenvalues) {

   int n = 3;

   int cyc, j, q, p;
   float spp, spq, t, s, c, theta, tau, h, ap, aq, a_pq;

   clipper::Matrix<double> &mat = m;
   clipper::Matrix<double> evec(3,3);

   float eval[3];
   float b[3];
   float z[3];

   // Set evec to identity, eval & b to diagonal, z to 0.
   for (p = 0; p < n; p++)
      for (j = 0; q < n; q++)
         evec(p,q) = 0.0;
      
   for (p = 0; p < n; p++) {
      evec(p,p) = 1.0;
      eval[p] = b[p] = m(p,p);
   }

   for ( cyc = 1; cyc <= 10; cyc++ ) {

      // calc sum of diagonal, off-diagonal
      spp = spq = 0.0;
      for ( p=0; p<n-1; p++ ) {
	 for ( q=p+1; q<n; q++ )
	    spq += fabs(m(p,q));
	 spp += fabs(m(p,p));
      }
      // printf("cyc %d spq %f spp %f\n", cyc, spq, spp);
      if ( spq <= 1.0e-5 * spp ) break;

      // zero z
      for ( p = 0; p < n; p++ ) z[p] = 0.0;

      // now try and reduce each off-diagonal element in turn
      for( p=0; p<n-1; p++ ) {
	 for( q=p+1; q<n; q++ ) {
	    a_pq = m(p,q);
	    h = eval[q] - eval[p];
	    if ( fabs(a_pq) > 1.0e-12f*fabs(h) ) {
	       theta = 0.5f*h/a_pq;
	       t = 1.0f/(fabs(theta) + sqrt(1.0f + theta*theta));
	       if ( theta < 0.0f ) t = -t;
	    } else {
	       t = a_pq/h;
	    }

	    // calc trig properties
	    c   = 1.0f/sqrt(1.0f+t*t);
	    s   = t*c;
	    tau = s/(1.0f+c);
	    h   = t * a_pq;

	    // update eigenvalues
	    z[p] -= h;
	    z[q] += h;
	    eval[p] -= h;
	    eval[q] += h;

	    // rotate the upper diagonal of the matrix
	    m(p,q) = 0.0f;
	    for ( j = 0; j < p; j++ ) {
	       ap = m(j,p);
	       aq = m(j,q);
	       m(j,p) = ap - s * ( aq + ap * tau );
	       m(j,q) = aq + s * ( ap - aq * tau );
	    }
	    for ( j = p+1; j < q; j++ ) {
	       ap = m(p,j);
	       aq = m(j,q);
	       m(p,j) = ap - s * ( aq + ap * tau );
	       m(j,q) = aq + s * ( ap - aq * tau );
	    }
	    for ( j = q+1; j < n; j++ ) {
	       ap = m(p,j);
	       aq = m(q,j);
	       m(p,j) = ap - s * ( aq + ap * tau );
	       m(q,j) = aq + s * ( ap - aq * tau );
	    }
	    // apply corresponding rotation to result
	    for ( j = 0; j < n; j++ ) {
	       ap = evec(j,p);
	       aq = evec(j,q);
	       evec(j,p) = ap - s * ( aq + ap * tau );
	       evec(j,q) = aq + s * ( ap - aq * tau );
	    }
	 }
      }

      for ( p = 0; p < n; p++ ) {
	 b[p] += z[p];
	 eval[p] = b[p];
      }
   }

   // sort the eigenvalues
   if ( sort_eigenvalues ) {
      for ( p = 0; p < n; p++ ) {
	 j = p;        // set j to index of largest remaining eval
	 for ( q = p+1; q < n; q++ )
	    if ( eval[q] < eval[j] ) j = q;
	 // Util::swap( eval[p], eval[j] );  // now swap evals, evecs
	 float tmp = eval[p];
	 eval[p] = eval[j];
	 eval[j] = tmp;
	 for ( q = 0; q < n; q++ ) {
	    // Util::swap( evec( q, p ), evec( q, j ) );
	    float tmp = evec(q,p);
	    evec(q,p) = evec(q,j);
	    evec(q,j) = tmp;
	 }
      }
   }


   // adjust the input matrix reference, for the vectors

   mat = evec;

   // return eigenvalues

   return std::tuple<double, double, double> (eval[0], eval[1], eval[2]);

}
