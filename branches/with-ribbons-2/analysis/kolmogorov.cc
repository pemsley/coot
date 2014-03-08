/* analysis/kolmogorov.cc
 * 
 * Copyright 2012 by The Medical Research Council
 * Author: Robert Nicholls
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the Lesser GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <algorithm>
#include <iostream>
#include "math.h"

#include "kolmogorov.hh"

// double
// nicholls::get_KS(const std::vector<double> &v1_in, const std::vector<double> &v2_in) {

//    std::vector<double> v1 = v1_in;
//    std::vector<double> v2 = v2_in;

//    std::sort(v1.begin(), v1.end());
//    std::sort(v2.begin(), v2.end());
   
//    unsigned int n1 = v1.size();
//    unsigned int n2 = v2.size();
//    if(n1==0 || n2==0) return -1.0;
    
//    unsigned int j=0;
//    double tmp = 0.0;
//    double result = 0.0;
//    for(unsigned int i=1; i<=n1; i++){
//       tmp = fabs( (double(i))/(double(n1)) - (double(j))/(double(n2)) );
//       if(tmp>result)
// 	 result = tmp;
//       if(i>=n1){                  //end of i array - increment j
// 	 if(j>=n2)               //end of j array - break
// 	    break;
// 	 j++;
// 	 i--;
//       } else if(j<n2){            //not end of i or j arrays
// 	 if(v1[i]>=v2[j]){       //if F[i]==F[j] then increment i and j
// 	    j++;
// 	    if(v1[i]>v2[j])     //if F[i]>F[j] then just increment j
// 	       i--;
// 	 }                       //if F[i]<F[j] then just increment i
//       }
//    }
//    return result;
// }



// double
// nicholls::get_KS(const std::vector<double> &v1_in, const std::vector<double> &v2_in) { 

//    std::vector<double> v1 = v1_in;
//    std::vector<double> v2 = v2_in;
//    std::sort(v1.begin(),v1.end());
//    std::sort(v2.begin(),v2.end());
//    unsigned int n1 = v1.size();
//    unsigned int n2 = v2.size();
//    if (n1==0 || n2==0) return -1.0;
   
//    double tmp = 0.0;
//    double result = 0.0;
//    unsigned int i = 0;
//    unsigned int j = 0;
//    while (1){      
//       tmp = fabs((double(i))/double(n1)-(double(j))/double(n2));
//       if(tmp>result){
// 	 result = tmp;
//       }
//       std::cout.precision(3);
//       std::cout << i << " " << j << "\t" << (double)i/n1 << "\t" << (double)j/n2 << "\t"
// 		<< tmp << "\t" << result << std::endl;
//       if(i<n1){
// 	 if(j<n2){
// 	    if(( (double (i+1))/double(n1)) < ((double(j+1))/double(n2))){
// 	       i++;
// 	    } else {
// 	       j++;
// 	    }
// 	 } else {
// 	    i++;
// 	 }
//       } else if(j<n2){
// 	 j++;
//       } else {
// 	 break;
//       }
//    }
//    return result;
// }



double
nicholls::get_KS(const std::vector<double> &v1_in, const std::vector<double> &v2_in) { 

   std::vector<double> v1 = v1_in;
   std::vector<double> v2 = v2_in;
   std::sort(v1.begin(),v1.end());
   std::sort(v2.begin(),v2.end());
   unsigned int n1 = v1.size();
   unsigned int n2 = v2.size();
   if(n1==0 || n2==0) return -1.0;
   double tmp = 0.0;
   double result = 0.0;
   unsigned int i = 1;
   unsigned int j = 1;
   while(1){      
      tmp = fabs((double)i/n1-(double)j/n2);
      if(tmp>result){
	 result = tmp;
      }
      // std::cout.precision(3);
      // std::cout << std::endl << v1[i-1] << " " << v2[j-1] << "\t" << (double)i/n1 << "\t" << (double)j/n2 << "\t" << tmp << "\t" << result;
      if(i<n1){
	 if(j<n2){
	    if(v1[i]==v2[j]){
	       i++;
	       j++;
	    } else if(v1[i]<v2[j]){
	       i++;
	    } else {
	       j++;
	    }
	 } else {
	    break;
	 }
      } else {
	 break;
      }
   }
   return result;
}




// Assumption - data are positive.
// 
std::pair<double, double>
nicholls::get_KL(const std::vector<double> &v1, const std::vector<double> &v2) {

   std::vector<double> result;
   unsigned int n1 = v1.size();
   unsigned int n2 = v2.size();
   if(n1==0 || n2==0) return std::pair<double, double> (-1, -1);
    
   double BIN_SIZE = 0.05;
   double SCALE1 = 1.0/(double)n1;
   double SCALE2 = 1.0/(double)n2;

   //get discrete probabilities (scaled density so that sum=1)
   std::vector<double> P;
   std::vector<double> Q;
   unsigned int bin_no;
   for(unsigned int i=0; i<n1; i++){
      bin_no = (unsigned int) floor(v1[i]/BIN_SIZE);
      while (P.size()<=bin_no) {
	 P.push_back(0.0);
      }
      P[bin_no] += SCALE1;
   }
   for(unsigned int i=0; i<n2; i++){
      bin_no = (unsigned int)floor(v2[i]/BIN_SIZE);
      while(Q.size()<=bin_no){
	 Q.push_back(0.0);
      }
      Q[bin_no] += SCALE2;
   }
   while(P.size()<Q.size())
      P.push_back(0.0);
   while(Q.size()<P.size())
      Q.push_back(0.0);


   //calculate contribution of intensity bins to KL divergence
   std::vector<double> KL1;
   std::vector<double> KL2;
   double logP;
   double logQ;
   for(unsigned int i=0; i<P.size(); i++){
      if(P[i]>0.0 && Q[i]>0.0){
	 logP = log(P[i]);
	 logQ = log(Q[i]);
	 //cout << endl << logP << "\t" << logQ << "\t" << P[i]*(logP-logQ) << "\t" << Q[i]*(logQ-logP);
	 KL1.push_back(P[i]*(logP-logQ));
	 KL2.push_back(Q[i]*(logQ-logP));
      } else {
	 KL1.push_back(0.0);
	 KL2.push_back(0.0);
      }
   }
    
   //calculate Kullback-Leibler divergences
   double KLdiv1 = 0.0;
   double KLdiv2 = 0.0;
   for(unsigned int i=0; i<KL1.size(); i++){
      KLdiv1 += KL1[i];
      KLdiv2 += KL2[i];
   }
   return std::pair<double, double> (KLdiv1, KLdiv2);
}
