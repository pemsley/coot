/*
 * MoleculesToTriangles/CXXClasses/
 *
 * Copyright 2011 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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

//
//  lfit.h
//  AesopCD_macos
//

#ifndef lfit_h
#define lfit_h

void nrerror(char *string);
float *NRvector(int nl,int nh);
int *ivector_nr(int nl,int nh);
float **matrix(int nrl,int nrh,int ncl,int nch);
void free_vector (float *v, int nl, int nh);
void free_ivector_nr (int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void gaussj(float **a,int n,float **b, int m);
void lfit(float *x, float *y, float *sig, int ndata, 
          float *a, int ma, int *lista, int mfit, 
          float **covar, float *chisq, void (*funcs)(float x, float *afunc, int ma));
float Determinant(float *a, int n);
void ludcmp(float **a,int n,int *indx,float *d);
void ForcePolynomial (int mfit,
                      float *xes, float *yes, int *fixeds,
                      int count, int start, int finish);

#endif
