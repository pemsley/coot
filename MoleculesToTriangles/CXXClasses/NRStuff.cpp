/*
 * MoleculesToTriangles/CXXClasses/NRStuff.cpp
 *

https://github.com/erich666/GraphicsGems/blob/master/LICENSE.md

LICENSE

This code repository predates the concept of Open Source, and predates most licenses along such lines.
As such, the official license truly is:

EULA: The Graphics Gems code is copyright-protected. In other words, you cannot claim the text of the
code as your own and resell it. Using the code is permitted in any program, product, or library,
non-commercial or commercial. Giving credit is not required, though is a nice gesture. The code comes
as-is, and if there are any flaws or problems with any Gems code, nobody involved with Gems - authors,
editors, publishers, or webmasters - are to be held responsible. Basically, don't be a jerk, and
remember that anything free comes with no guarantee.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#define TRUE 1
#define FALSE 0
#define BIG (1.0e12)

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

#include "NRStuff.h"

/****** lincrv.c ******/
/* Ken Shoemake, 1994 */

/* DialASpline(t,a,p,m,n,Cn,interp,val) computes a point val at parameter
    t on a spline with knot values a and control points p. The curve will have
    Cn continuity, and if interp is TRUE it will interpolate the control points.
    Possibilities include Langrange interpolants, Bezier curves, Catmull-Rom
    interpolating splines, and B-spline curves. Points have m coordinates, and
    n+1 of them are provided. The work array must have room for n+1 points.
 */

void CoordSpline::DialASpline(double t, const std::vector<double> &a,  const std::vector<FCXXCoord> &p, int Cn, int interp, std::vector<FCXXCoord> &output, const int idx, std::vector<FCXXCoord> &work)
{
    int i, j, k, h, lo, hi;

    int n = (int)(p.size()) - 1;

    if (Cn>n-1) Cn = n-1;       /* Anything greater gives one polynomial */
    for (k=0; t> a[k]&&k<int(a.size()); k++);    /* Find enclosing knot interval */
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
        for (i=lo; i<=hi; i++){
          work[i][0] = p[i][0];
          work[i][1] = p[i][1];
          work[i][2] = p[i][2];
        }
        for (j=1; j<=Cn; j++) {
            for (i=lo; i<=hi-j; i++) {
                double t0=(a[i+j]-t)/(a[i+j]-a[i]), t1=1-t0;
                work[i][0] = t0*work[i][0] + t1*work[i+1][0];
                work[i][1] = t0*work[i][1] + t1*work[i+1][1];
                work[i][2] = t0*work[i][2] + t1*work[i+1][2];
            }
        }
        h = 1+Cn-drop;
    } else {                    /* Prepare for B-spline steps */
        if (lo<0) {h += lo; lo = 0;}
        for (i=lo; i<=lo+h; i++){
          work[i][0] = p[i][0];
          work[i][1] = p[i][1];
          work[i][2] = p[i][2];
        }
        if (h<0) h = 0;
    }
    for (j=0; j<h; j++) {
        int tmp = 1+Cn-j;
        for (i=h-1; i>=j; i--) {
            double t0=(a[lo+i+tmp]-t)/(a[lo+i+tmp]-a[lo+i]), t1=1-t0;
            work[lo+i+1][0] = t0*work[lo+i][0] + t1*work[lo+i+1][0];
            work[lo+i+1][1] = t0*work[lo+i][1] + t1*work[lo+i+1][1];
            work[lo+i+1][2] = t0*work[lo+i][2] + t1*work[lo+i+1][2];
        }
    }
    output[idx][0] = work[lo+h][0];
    output[idx][1] = work[lo+h][1];
    output[idx][2] = work[lo+h][2];
}

/*** lincrvtest.c ***/

std::vector <FCXXCoord> CoordSpline::SplineCurve(const std::vector<FCXXCoord> &ctlPts, int nsteps, int Cn, int iinterp){
   unsigned int i;
   std::vector<double> knots;
   double t;
   double maxx = -BIG;
   double minx = BIG;
   double maxy = -BIG;
   double miny = BIG;
   double maxz = -BIG;
   double minz = BIG;
   double maxt = -BIG;
   double mint = BIG;
   double tstep;
   double knotstep;

   FCXXCoord outputi;
   int interp;

   if (!iinterp){
     interp = FALSE;
   }else{
     interp =TRUE;
   }

   for(i=0;i<ctlPts.size();i++){
     minx = min(ctlPts[i][0],minx);
     maxx = max(ctlPts[i][0],maxx);
   }

   for(i=0;i<ctlPts.size();i++){
     miny = min(ctlPts[i][1],miny);
     maxy = max(ctlPts[i][1],maxy);
   }

   for(i=0;i<ctlPts.size();i++){
     minz = min(ctlPts[i][2],minz);
     maxz = max(ctlPts[i][2],maxz);
   }

   mint = minx;
   mint = min(mint,miny);
   mint = min(mint,minz);

   maxt = maxx;
   maxt = max(maxt,maxy);
   maxt = max(maxt,maxz);

   tstep = (maxt-mint)/double(nsteps-1);
   knotstep = (maxt-mint)/(ctlPts.size()-1);

   for(i=0;i<ctlPts.size();i++)
     knots.push_back(mint + (double)i*knotstep);
   knots.push_back(tstep);

   std::vector<FCXXCoord> output(nsteps);
   std::vector<FCXXCoord> work(ctlPts.size());
   for (int ii=0;ii<nsteps;ii++){
     t = mint + ii*tstep;
     DialASpline(t, knots, ctlPts, Cn, interp,output,ii,work);
   }

   return output;

}
