//   CCP4 Molecular Graphics Program
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Library.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.



/****** lincrv.c ******/
/* Ken Shoemake, 1994 */

#include <vector>
#include "cartesian.h"

inline Cartesian lerp(const double &t, const double &a0, const double &a1, const Cartesian &p0, const Cartesian &p1)
{
     double t0=(a1-t)/(a1-a0), t1=1-t0;

    Cartesian p;
    p = t0*p0 + t1*p1;
    return p;
}

/* DialASpline(t,a,p,m,n,Cn,interp,val) computes a point val at parameter
    t on a spline with knot values a and control points p. The curve will have
    Cn continuity, and if interp is TRUE it will interpolate the control points.
    Possibilities include Langrange interpolants, Bezier curves, Catmull-Rom
    interpolating splines, and B-spline curves. Points have m coordinates, and
    n+1 of them are provided. The work array must have room for n+1 points.
 */

Cartesian DialASpline(double t, const std::vector<double> &a,  const std::vector<Cartesian> &p, int Cn, int interp)
{
    register int i, j, k, h, lo, hi;
    Cartesian val;
    std::vector<Cartesian> work(p.size());

    int n = (int)(p.size()) - 1;

    if (Cn>n-1) Cn = n-1;       /* Anything greater gives one polynomial */
    for (k=0; t> a[k]&&k<a.size(); k++);    /* Find enclosing knot interval */
    for (h=k; t==a[k]; k++);    /* May want to use fewer legs */
    if (k>n) {k = n; if (h>k) h = k;}
    h = 1+Cn - (k-h); k--;
    lo = k-Cn; hi = k+1+Cn;

    if (interp) {               /* Lagrange interpolation steps */
        int drop=0;
        if (lo<0) {lo = 0; drop += Cn-k;
                   if (hi-lo<Cn) {drop += Cn-hi; hi = Cn;}}
        if (hi>n) {hi = n; drop += k+1+Cn-n;
                   if (hi-lo<Cn) {drop += lo-(n-Cn); lo = n-Cn;}}
        for (i=lo; i<=hi; i++) work[i] = p[i];
        for (j=1; j<=Cn; j++) {
            for (i=lo; i<=hi-j; i++) {
                work[i] = lerp(t,a[i],a[i+j],work[i],work[i+1]);
            }
        }
        h = 1+Cn-drop;
    } else {                    /* Prepare for B-spline steps */
        if (lo<0) {h += lo; lo = 0;}
        for (i=lo; i<=lo+h; i++) work[i]= p[i];
        if (h<0) h = 0;
    }
    for (j=0; j<h; j++) {
        int tmp = 1+Cn-j;
        for (i=h-1; i>=j; i--) {
            work[lo+i+1] = lerp(t,a[lo+i],a[lo+i+tmp],work[lo+i],work[lo+i+1]);
        }
    }
    val = work[lo+h];
    return val;
}
