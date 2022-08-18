//
//  lfit.c
//  AesopCD_macos
//
//  Created by Martin Noble on 01/08/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "lfit.h"

static void polynomials(float x,float *afunc,int ma)
{
    int i;
    float result;
    
    result = 1;
    for (i=1; i<= ma; i++){
        afunc[i] = result;
        result = result*x;
    }
}

void nrerror(char *string)
{
    fprintf(stderr,"Numerical recipes run time error\n");
}

float *NRvector(int nl,int nh)
{
    float *v;
    
    v = (float *) malloc ((unsigned) (nh-nl+1) * sizeof(float));
    if (!v) nrerror("vector allocation failed");
    return v-nl;
}

int *ivector_nr(int nl,int nh)
{
    int *v;
    
    v = (int *) malloc ((unsigned) (nh-nl+1) * sizeof(int));
    if (!v) nrerror("vector allocation failed");
    return v-nl;
}

float **matrix(int nrl,int nrh,int ncl,int nch)
{
    int i, j;
    float **m;
    
    m = (float **) malloc ((unsigned) (nrh-nrl+1) * sizeof(float *));
    if (!m) nrerror("failed allocating space for row pointers");
    m -= nrl;
    
    for (i = nrl; i <= nrh; i++){
        m[i] = (float *) malloc((unsigned) (nch-ncl + 1)*sizeof(float));
        if (!m[i]) nrerror ("failed allocating space for a row");
        m[i] -= ncl;
        for (j=ncl; j<=nch; j++) m[i][j] = 0.;
    }
    return m;
}

void free_vector (float *v, int nl, int nh)
{
    free ((char *) (v+nl));
}

void free_ivector_nr (int *v, int nl, int nh)
{
    free ((char *) (v+nl));
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
    int i;
    for (i=nrh; i>= nrl; i--) free ((char *) (m[i]+ncl));
    free ((char *) (m + nrl));
}

#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}

void gaussj(float **a, int n, float **b, int m)
{
    int *indxc, *indxr, *ipiv;
    int i, icol = 0, irow = 0, j, k = 0, l, ll;
    float big, dum, pivinv;
    
    indxc=ivector_nr(1,n);
    indxr=ivector_nr(1,n);
    ipiv=ivector_nr(1,n);
    
    for (j=1; j<=n; j++) ipiv[j]=0;
    
    for (i=1; i<=n; i++){
        
        big=0.0;
        for (j=1; j<=n; j++)
            if (ipiv[j] != 1)
                for (k=1; k<=n; k++){
                    if (ipiv[k] ==0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] >1) nrerror("GAUSSJ: Singular matrix-1");
                }
        ++(ipiv[icol]);
        
        if (irow != icol) {
            for (l=1; l<=n; l++) SWAP (a[irow][l], a[icol][l]);
            for (l=1; l<=m; l++) SWAP (b[irow][l], b[icol][l]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) nrerror("GaussJ: Singular matrix-2");
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l=1; l<=n; l++) a[icol][l] *= pivinv;
        for (l=1; l<=m; l++) b[icol][l] *= pivinv;
        for (ll=1; ll<=n; ll++)
            if (ll != icol) {
                dum=a[ll][icol];
                a[ll][icol] = 0.0;
                for (l=1; l<=n; l++) a[ll][l] -= a[icol][l]*dum;
                for (l=1; l<=m; l++) b[ll][l] -= b[icol][l]*dum;
            }
        
    }
    for (l=n; l>=1; l--) {
        if (indxr[l] != indxc[l])
            SWAP (a[k][indxr[l]], a[k][indxc[l]]);
    }
    free_ivector_nr(ipiv, 1, n);
    free_ivector_nr(indxr, 1, n);
    free_ivector_nr(indxc, 1, n);
}

float sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

void lfit(float *x, float *y, float *sig, int ndata, 
          float *a, int ma, int *lista, int mfit, 
          float **covar, float *chisq, void (*funcs)(float x, float *afunc, int ma))

/*   float x[], y[], sig[];     observations and estimated errors */
/*float a[];                    return coefficients, + initial values */
/*int ma;                       order of basis vector */
/*int lista[];                  identity of coefficients to be fit */
/*int mfit;                     length of lista vector */
/*float **covar;                covariance matrix */
/*float *chisq;                 chisq for parameters: length ma */
/*void (*funcs)();              */

{
    int k, kk, j, ihit, i;
    float ym, wt, sum, sig2i, **beta, *afunc;
    
    beta = matrix(1, ma, 1, 1);
    afunc = NRvector(1,ma);
    kk = mfit + 1;
    for (j=1; j<= ma; j++){
        ihit = 0;
        for (k=1; k<=mfit; k++)
            if (lista[k] == j) ihit++;
        if (ihit == 0)
            lista[kk++] = j;
        else if (ihit > 1) nrerror("Bad lista permutation in lfit-1");
    }
    if (kk != (ma+1)) nrerror("Bad lista permutation in lfit-2");
    for (j=1; j<= mfit; j++){
        beta[j][1] = 0.;
        for (k=1; k<=mfit; k++) covar[j][k] = 0.;
    }
    for (i=1; i<=ndata; i++) {
        (*funcs)(x[i], afunc, ma);
        ym = y[i];
        if (mfit < ma)
            for (j=(mfit+1); j<= ma; j++)
                ym -= a[lista[j]] * afunc[lista[j]];
        sig2i = 1.0/SQR(sig[i]);
        for (j=1; j<=mfit; j++){
            wt = afunc[lista[j]]*sig2i;
            for (k=1; k<=j; k++)
                covar[j][k] += wt*afunc[lista[k]];
            beta[j][1] += ym*wt;
        }
    }
    
    if (mfit > 1)
        for (j=2;  j<= mfit; j++)
            for (k=1; k<= j-1; k++)
                covar [k][j] = covar [j][k];
    gaussj(covar, mfit, beta, 1);
    for (j=1; j <= mfit; j++) 
        a[lista[j]] = beta[j][1];
    *chisq = 0.0;
    for (i=1; i<= ndata; i++){
        (*funcs) (x[i], afunc, ma);
        for (sum=0.0, j=1; j<= ma; j++)
            sum += a[j]*afunc[j];
        *chisq += SQR((y[i]-sum)/sig[i]);
    }
    free_vector(afunc, 1, ma);
    free_matrix(beta, 1, ma, 1, 1);
}

float Determinant(float *a, int n){
    float **NrMatrix;
    int *indx;
    float det;
    int i, j, k;
    
    NrMatrix = matrix(1, n, 1, n);
    indx = ivector_nr (1, n);
    
    for (i=0, k=0; i<n; i++){
        for (j=0; j<n; j++, k++){
            NrMatrix[i+1][j+1] = a[k];
        }
    }
    ludcmp (NrMatrix, n, indx, &det);
    for (i=1; i<=n; i++){
        det *= NrMatrix[i][i];
    }
    free_matrix (NrMatrix, 1, n, 1, n);
    free_ivector_nr (indx, 1, n);
    
    return(det);
}

#define TINY 1.0e-20;

void ludcmp(float **a,int n,int *indx,float *d)
{
    int i,imax = 0,j,k;
    float big,dum,sum,temp;
    float *vv;
    
    
    vv=NRvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) {
            printf("Singular matrix in routine LUDCMP\n");
            return;
        }
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_vector(vv,1,n);
    
    return;
}

#undef TINY

void ForcePolynomial (int mfit,
                      float *xes, float *yes, int *fixeds,
                      int count, int start, int finish)
{
    float *x, *y, *sig;
    int ndata;
    float *a;
    int ma;
    int *lista;
    float *afunc;
    float **covar;
    float chisq;
    
    int i, j;
    
    
    ma = mfit;
    ndata = count;
    
    x = NRvector(1,ndata);
    y = NRvector(1,ndata);
    sig = NRvector(1,ndata);
    a = NRvector(1,ma);
    afunc = NRvector(1,ma);
    covar = matrix(1,ma,1,ma);
    lista = ivector_nr(1,ma);
    
    for (i=1; i<=ma; i++)
        lista[i] = i;
    
    for (i = 1; i<= ndata; i++){
        x[i] = xes[i-1];
        y[i] = yes[i-1];
        sig[i] = 0.1;
    }
    
    lfit (x, y, sig, ndata, a, ma, lista, mfit, covar, &chisq, polynomials);

    /*
     *	Force "y"es to fall on the line unless they have been "fixed"
     */
    
    for (i=start; i<= finish; i++) {
        polynomials(x[i+1], afunc, ma);
        if (!fixeds[i]) {
            yes[i] = 0.;
            for (j=1; j<=ma; j++)
                yes[i] = yes[i] + a[j]*afunc[j];
        }
    }
    
    free_vector (x, 1, ndata);
    free_vector (y, 1, ndata);
    free_vector (sig, 1, ndata);
    free_vector (a, 1, ma);
    free_vector (afunc, 1, ma);
    free_matrix (covar, 1, ma, 1, ma);
    free_ivector_nr (lista, 1, ma);
}

