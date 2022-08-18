//
//  lfit.h
//  AesopCD_macos
//
//  Created by Martin Noble on 01/08/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
//

#ifndef lfit_h
#define lfit_h

//static void polynomials(float x,float *afunc,int ma);
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