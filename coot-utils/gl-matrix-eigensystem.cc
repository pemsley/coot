/* coot-utils/gl-matrix-eigensystem.cc
 *
 * Copyright 2024,  by Global Phasing Ltd.
 * Author: ClAuS Flensburg, Clemens Vonrhein, Gerard Bricogne
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <gsl/gsl_matrix.h>
#include "gl-matrix.h"

std::pair<bool,GL_matrix>
GL_matrix::eigensystem() const {

   double a_data[] = { mat[0], mat[1], mat[ 2],
                       mat[4], mat[5], mat[ 6],
                       mat[8], mat[9], mat[10] };

   gsl_matrix_view m = gsl_matrix_view_array (a_data, 3, 3);

   // test the princal minors for being positive (required for
   // m.matrix being positive definite.

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
      gsl_vector *eval = gsl_vector_alloc(3);
      gsl_matrix *evec = gsl_matrix_alloc(3, 3);
      gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);

      gsl_eigen_symmv(&m.matrix, eval, evec, w);
      gsl_eigen_symmv_free(w);
      gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

      gsl_set_error_handler(old_handler);

      float L[3];
      float V[9];
      float R[16];
      int l = 0;
      for (int i = 0; i < 3; i++) {
        double eval_i = gsl_vector_get(eval, i);
        //printf ("eigenvalue and eigenvalue^(1/2) = %g  %g ", eval_i, sqrt(eval_i));
        L[i] = (float) sqrt(eval_i);

        //printf ("eigenvector = ");
        for(int k=0; k<3; k++){
          double q = gsl_matrix_get(evec, k, i); // extract the columns
          //printf(" %g", q);
          V[l++] = (float) q;
        }
        //printf ("\n");
      }
      gsl_vector_free(eval);
      gsl_matrix_free(evec);

      // indexing like in GL_matrix of nine floats in gl-matrix.cc
      // Scale columns by the square-root of eigenvalues.
      R[ 0] = V[0]*L[0]; R[ 1] = V[1]*L[0]; R[ 2] = V[2]*L[0]; R[ 3] = 0.0f;
      R[ 4] = V[3]*L[1]; R[ 5] = V[4]*L[1]; R[ 6] = V[5]*L[1]; R[ 7] = 0.0f;
      R[ 8] = V[6]*L[2]; R[ 9] = V[7]*L[2]; R[10] = V[8]*L[2]; R[11] = 0.0f;
      R[12] = 0.0f     ; R[13] = 0.0f     ; R[14] = 0.0f     ; R[15] = 1.0f;

      return std::pair<bool, GL_matrix> (1, GL_matrix(
                                                      R[ 0], R[ 1], R[ 2],
                                                      R[ 4], R[ 5], R[ 6],
                                                      R[ 8], R[ 9], R[10]
                                                      ));
   }
}

