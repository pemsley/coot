/*
     util/catmull.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/


/*** lincrvtest.c ***/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "lincrv.h"
#include "catmull.h"
#include "cartesian.h"

#define TRUE 1
#define FALSE 0
#define BIG (1.0e12)

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

std::vector <Cartesian> SplineCurve(const std::vector<Cartesian> &ctlPts, int nsteps, int Cn, int iinterp){
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

   Cartesian outputi;
   int interp;

   if (!iinterp){
     interp = FALSE;
   }else{
     interp =TRUE;
   }

   for(i=0;i<ctlPts.size();i++){
     minx = min(ctlPts[i].get_x(),minx);
     maxx = max(ctlPts[i].get_x(),maxx);
   }

   for(i=0;i<ctlPts.size();i++){
     miny = min(ctlPts[i].get_y(),miny);
     maxy = max(ctlPts[i].get_y(),maxy);
   }

   for(i=0;i<ctlPts.size();i++){
     minz = min(ctlPts[i].get_z(),minz);
     maxz = max(ctlPts[i].get_z(),maxz);
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

   //std::cout << "mint: " << mint << "\n";
   //std::cout << "maxt: " << maxt << "\n";
   
   std::vector <Cartesian>output(nsteps);
   std::vector<CART3D> work(ctlPts.size());
   for (int ii=0;ii<nsteps;ii++){
     t = mint + ii*tstep;
     //outputi = DialASpline(t, knots, ctlPts, Cn, interp,output,ii);
     //output.push_back(outputi);
     DialASpline(t, knots, ctlPts, Cn, interp,output,ii,work);
   }

   return output;
   
}
