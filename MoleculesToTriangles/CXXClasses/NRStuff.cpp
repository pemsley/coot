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

void CoordSpline::DialASpline(float t, const std::vector<float> &a,  const std::vector<FCXXCoord> &p, int Cn, int interp, std::vector<FCXXCoord> &output, const int idx, std::vector<FCXXCoord> &work)
{
    int i, j, k, h, lo, hi;

    int n = (int)(p.size()) - 1;

    if (Cn>n-1) Cn = n-1;       /* Anything greater gives one polynomial */
    for (k=0; k<int(a.size()) && t> a[k]; k++);
    for (h=k; k<int(a.size()) && t==a[k]; k++);
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
          work[i][3] = p[i][3];
        }
        for (j=1; j<=Cn; j++) {
            for (i=lo; i<=hi-j; i++) {
                float t0=(a[i+j]-t)/(a[i+j]-a[i]), t1=1-t0;
                work[i][0] = t0*work[i][0] + t1*work[i+1][0];
                work[i][1] = t0*work[i][1] + t1*work[i+1][1];
                work[i][2] = t0*work[i][2] + t1*work[i+1][2];
                work[i][3] = t0*work[i][3] + t1*work[i+1][3];
            }
        }
        h = 1+Cn-drop;
    } else {                    /* Prepare for B-spline steps */
        if (lo<0) {h += lo; lo = 0;}
        for (i=lo; i<=lo+h; i++){
          work[i][0] = p[i][0];
          work[i][1] = p[i][1];
          work[i][2] = p[i][2];
          work[i][3] = p[i][3];
        }
        if (h<0) h = 0;
    }
    for (j=0; j<h; j++) {
        int tmp = 1+Cn-j;
        for (i=h-1; i>=j; i--) {
            float t0=(a[lo+i+tmp]-t)/(a[lo+i+tmp]-a[lo+i]), t1=1-t0;
            work[lo+i+1][0] = t0*work[lo+i][0] + t1*work[lo+i+1][0];
            work[lo+i+1][1] = t0*work[lo+i][1] + t1*work[lo+i+1][1];
            work[lo+i+1][2] = t0*work[lo+i][2] + t1*work[lo+i+1][2];
            work[lo+i+1][3] = t0*work[lo+i][3] + t1*work[lo+i+1][3];
        }
    }
    output[idx][0] = work[lo+h][0];
    output[idx][1] = work[lo+h][1];
    output[idx][2] = work[lo+h][2];
    output[idx][3] = work[lo+h][3];
}

/*** lincrvtest.c ***/

std::vector <FCXXCoord> CoordSpline::SplineCurve(const std::vector<FCXXCoord> &ctlPts, int nsteps, int Cn, int iinterp){

   if (false) {
      std::cout << "CoordSpline::SplineCurve() debug: ctlPts size() " << ctlPts.size() << std::endl;
      for (unsigned int ii=0; ii<ctlPts.size(); ii++)
         std::cout << "CoordSpline::SplineCurve() debug: ctlPts " << ii << " " << ctlPts[ii] << std::endl;
      std::cout << "CoordSpline::SplineCurve() debug: nsteps " << nsteps << std::endl;
      std::cout << "CoordSpline::SplineCurve() debug: Cn " << Cn << std::endl;
      std::cout << "CoordSpline::SplineCurve() debug: iinterp " << iinterp << std::endl;
   }

   std::vector<FCXXCoord> output;
   if(ctlPts.size()==1){
       output.push_back(ctlPts[0]);
       return output;
   }

   unsigned int i;
   std::vector<float> knots;
   float t;
   float maxx = -BIG;
   float minx = BIG;
   float maxy = -BIG;
   float miny = BIG;
   float maxz = -BIG;
   float minz = BIG;
   float maxt = -BIG;
   float mint = BIG;
   float tstep;
   float knotstep;

   FCXXCoord outputi;
   int interp;

   if (!iinterp){
     interp = FALSE;
   }else{
     interp =TRUE;
   }

   for(i=0;i<ctlPts.size();i++){
     minx = std::min(ctlPts[i][0],minx);
     maxx = std::max(ctlPts[i][0],maxx);
   }

   for(i=0;i<ctlPts.size();i++){
     miny = std::min(ctlPts[i][1],miny);
     maxy = std::max(ctlPts[i][1],maxy);
   }

   for(i=0;i<ctlPts.size();i++){
     minz = std::min(ctlPts[i][2],minz);
     maxz = std::max(ctlPts[i][2],maxz);
   }

   mint = minx;
   mint = std::min(mint,miny);
   mint = std::min(mint,minz);

   maxt = maxx;
   maxt = std::max(maxt,maxy);
   maxt = std::max(maxt,maxz);

   tstep = (maxt-mint)/float(nsteps-1);
   knotstep = (maxt-mint)/(ctlPts.size()-1);

   for(i=0;i<ctlPts.size();i++)
     knots.push_back(mint + (float)i*knotstep);
   // knots.push_back(tstep);
   knots.push_back(maxt);

   output.resize(nsteps);
   std::vector<FCXXCoord> work(ctlPts.size());
   for (int ii=0;ii<nsteps;ii++){
     t = mint + ii*tstep;
     DialASpline(t, knots, ctlPts, Cn, interp,output,ii,work);
   }

   return output;

}
