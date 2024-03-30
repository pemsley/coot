/*
 *  NRStuff.mm
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 17/07/2008.
 *  Copyright 2008 LMB, Oxford University. All rights reserved.
 *
 */

#include "NRStuff.h"

void NRSpline::calculateYDoublePrime(float ypLow, float ypHigh){
    int n = int(x.size());
    
    std::vector<float> u(n);
    y2.resize(n);
    if (ypLow > 0.99e30){
        y2[0]=u[0]=0.0;
    }
    else {
        y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-ypLow);
    }
    
    float sig, p;
    for (int i=1;i<n-1;i++) {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    float qn, un;
    if (ypHigh > 0.99e30){
        qn=un=0.0;
    }
    else {
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypHigh-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (int k=n-2;k>=0;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
    yDoublePrimeCalculated = true;
    
}

float NRSpline::yForXEquals(const float xVal)  {
    if (!yDoublePrimeCalculated) calculateYDoublePrime(1e30f, 1e30f);
    unsigned long n = x.size();
    //int n1 = y.size();
    //int n2 = y2.size();
    unsigned long klo,khi,k;
    float h,b,a;
    
    klo=0;
    khi=n-1;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (x[k] > xVal) khi=k;
        else klo=k;
    }
    h=x[khi]-x[klo];
    if (h == 0.0) {
        //std::cout << "Bad x value";
    }

    //MN Dealing with spline discontinuities by returning y element if x close to control point
    if (fabs(xVal-x[khi]) < 0.001) {
        return y[khi];
    }
    else if (fabs(xVal-x[klo]) < 0.001) {
        return y[klo];
    }

    a=(x[khi]-xVal)/h;
    b=(xVal-x[klo])/h;
    float c = y[khi];
    float d = y[klo];
    float e = y2[khi];
    float f = y2[klo];
    float answer=a*d+b*c+((a*a*a-a)*f+(b*b*b-b)*e)*(h*h)/6.0;
    return answer;
}
