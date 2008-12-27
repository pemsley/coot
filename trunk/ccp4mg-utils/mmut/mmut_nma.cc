/*
     mmut/mmut_nma.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

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
#include "mmut_nma.h"
#include <vector>
#include <matrix.h>
#include <cartesian.h>
#include <geomutil.h>

#include <iostream>
#include <iomanip>

#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define kBA 1.3806503e-03 // Boltzmann constant scaled for angstroms.
#define TEMP 298.155



NormalModeAnalysis::NormalModeAnalysis(const std::vector<Cartesian> &carts, const int type, const double cutoff){
  if(type==MMUT_ANM){
    if(cutoff<0.0)
      NMA_ANM(carts,13.0);
    else
      NMA_ANM(carts,cutoff);
    _type = MMUT_ANM;
  }else{
    if(cutoff<0.0)
      NMA_GNM(carts,7.0);
    else
      NMA_GNM(carts,cutoff);
    _type = MMUT_GNM;
  }
}

/*
NormalModeAnalysis::NormalModeAnalysis(const std::vector<Cartesian> &carts, const std::string &type, const double cutoff){
  if(type=="ANM"){
    if(cutoff<0.0)
      NMA_ANM(carts,13.0);
    else
      NMA_ANM(carts,cutoff);
    _type = MMUT_ANM;
  }else{
    if(cutoff<0.0)
      NMA_GNM(carts,7.0);
    else
      NMA_GNM(carts,cutoff);
    _type = MMUT_GNM;
  }
}
*/

void NormalModeAnalysis::NMA_ANM(const std::vector<Cartesian> &carts, const double cutoff){
  /* Generate the Kirchoff matrix */

  matrix k(carts.size(),carts.size());

  for(unsigned i=0;i<k.get_rows();i++){
    for(unsigned j=0;j<i;j++){
      k(i,j) = ((LineLength(carts[i],carts[j])<cutoff)?-1:0);
      k(j,i) = k(i,j);
    }
  }
  for(unsigned i=0;i<k.get_rows();i++){
    k(i,i) = 0.0;
    for(unsigned j=0;j<k.get_columns();j++){
      if(i!=j) k(i,i) -= k(i,j);
    }
  }

  /* Generate the Hessian for anisotropic network model */
  matrix a(k.get_rows()*3,k.get_columns()*3);

  for(unsigned i=0;i<k.get_rows();i++){
    for(unsigned j=0;j<k.get_columns();j++){
      double len = LineLength(carts[i],carts[j]);
      len = len * len;
      if(i!=j){
	double xij = (carts[i].get_x()-carts[j].get_x());
	double yij = (carts[i].get_y()-carts[j].get_y());
	double zij = (carts[i].get_z()-carts[j].get_z());
        a(i*3,  j*3)   = k(i,j)*xij*xij/len;
        a(i*3+1,j*3+1) = k(i,j)*yij*yij/len;
        a(i*3+2,j*3+2) = k(i,j)*zij*zij/len;
        a(i*3,  j*3+1) = a(i*3+1,j*3)   = k(i,j)*xij*yij/len;
        a(i*3,  j*3+2) = a(i*3+2,j*3)   = k(i,j)*xij*zij/len;
        a(i*3+1,j*3+2) = a(i*3+2,j*3+1) = k(i,j)*yij*zij/len;
      }
    }
  }

  for(unsigned i=0;i<k.get_rows();i++){
    for(unsigned j=0;j<k.get_columns();j++){
       if(i!=j){
         a(i*3,i*3) -= a(i*3,j*3);
         a(i*3+1,i*3+1) -= a(i*3+1,j*3+1);
         a(i*3+2,i*3+2) -= a(i*3+2,j*3+2);
         a(i*3+1,i*3)   -= a(i*3+1,j*3);
         a(i*3,i*3+1)   -= a(i*3,j*3+1);
         a(i*3+1,i*3+2) -= a(i*3+1,j*3+2);
         a(i*3+2,i*3+1) -= a(i*3+2,j*3+1);
         a(i*3+2,i*3)   -= a(i*3+2,j*3);
         a(i*3,i*3+2)   -= a(i*3,j*3+2);
       }
    }
  }
  //std::cout << gnm[0] << "\n";
  //std::cout << gnm[1] << "\n";
  //std::cout << a.get_columns() << "\n";
  //std::cout << a.Determinant() << "\n";
  //std::cout << a << "\n";
  // The number crunching
  hessian = a;
  eigen = a.SortEigenvalues(a.Eigen());
}

void NormalModeAnalysis::NMA_GNM(const std::vector<Cartesian> &carts, const double cutoff){

  /* Generate the Kirchoff matrix for isotropic Gaussian network model */
  matrix a(carts.size(),carts.size());
  for(unsigned i=0;i<a.get_rows();i++){
    for(unsigned j=0;j<i;j++){
      a(i,j) = ((LineLength(carts[i],carts[j])<cutoff)?-1:0);
      a(j,i) = a(i,j);
    }
  }
  for(unsigned i=0;i<a.get_rows();i++){
    a(i,i) = 0.0;
    for(unsigned j=0;j<a.get_columns();j++){
      if(i!=j) a(i,i) -= a(i,j);
    }
  }

  // The number crunching
  hessian = a;
  eigen = a.SortEigenvalues(a.Eigen());

}

std::vector<double>  NormalModeAnalysis::GetBValues() const {
  // These B-values have to be divided by the force constant.
  // Unfortunately we do not know that, so these values need 
  // to be scaled to agree with the experimental B-factors. 
  // This isn't quite self-defeating, we can now calculate the
  // force constant and use it for other purposes.
  std::vector<double> bvalues;

  if(_type==MMUT_ANM){
    std::cerr << "Warning, B-value calculation for anisotropic method not yet implemented!\n";
    return bvalues;
  }

  matrix my_hinv = GetDiagonalInverseHessian();
  //std::cout << "Inverse Hessian\n";
  //std::cout << my_hinv << "\n";

  double prefactor=8.0*M_PI*M_PI*kBA*TEMP;

  for(unsigned i=0;i<my_hinv.get_rows();i++){
     bvalues.push_back(prefactor*my_hinv(i,i));
  }

  return bvalues;
}

void NormalModeAnalysis::StoreInverseHessian() {
  hinv = GetInverseHessian();
}

std::vector<std::vector<double> > NormalModeAnalysis::GetModeShapes(const double gammainv) const {

  // Almost certainly only correct for GNM at the moment.
  // Returns N-1 (or N-6 for ANM) vectors of N-values describing mobility of nodes.
  // gammainv can be calculated by comparing experimental with calculated B-values.

  std::vector<std::vector<double> > mode_shapes;
  matrix eigenvecs = eigen[0];
  matrix eigenvals = eigen[1];
  int start = 1;
  if(GetType()==MMUT_ANM){
    start = 6;
  }

  double prefactor=8.0*M_PI*M_PI*kBA*TEMP*gammainv;

  for(unsigned i=start;i<eigenvecs.get_columns();i++){
    mode_shapes.push_back(std::vector<double>(eigenvecs.get_rows()));
    if(fabs(eigenvals(i,0))>1e-7){
      double einv = 1.0/eigenvals(i,0);
      for(unsigned ii=0;ii<eigenvecs.get_rows();ii++){
           mode_shapes.back()[ii] = prefactor*einv*eigenvecs(ii,i)*eigenvecs(ii,i);
      }
    }
  }
  return mode_shapes;
}

matrix NormalModeAnalysis::GetDiagonalInverseHessian() const {

  if(hinv.get_rows()==eigen[0].get_rows()) 
    return hinv;
     
  matrix eigenvecs = eigen[0];
  matrix eigenvals = eigen[1];
  matrix my_hinv = matrix(eigenvecs.get_rows(),eigenvecs.get_rows());

  int start = 1;
  if(GetType()==MMUT_ANM){
    start = 6;
  }

  for(unsigned i=start;i<eigenvecs.get_columns();i++){
    if(fabs(eigenvals(i,0))>1e-7){
      double einv = 1.0/eigenvals(i,0);
      for(unsigned ii=0;ii<my_hinv.get_rows();ii++){
           my_hinv(ii,ii) += einv*eigenvecs(ii,i)*eigenvecs(ii,i);
      }
    }
  }
  return my_hinv;
  
}

matrix NormalModeAnalysis::GetInverseHessian() const {

  if(hinv.get_rows()==eigen[0].get_rows()) 
    return hinv;
     
  matrix eigenvecs = eigen[0];
  matrix eigenvals = eigen[1];
  matrix my_hinv = matrix(eigenvecs.get_rows(),eigenvecs.get_rows());

  int start = 1;
  if(GetType()==MMUT_ANM){
    start = 6;
  }

  for(unsigned i=start;i<eigenvecs.get_columns();i++){
    if(fabs(eigenvals(i,0))>1e-7){
      double einv = 1.0/eigenvals(i,0);
      matrix u(eigenvecs.get_rows(),1);
      for(unsigned j=0;j<eigenvecs.get_rows();j++){
        u(j,0) = eigenvecs(j,i);
      }
      matrix uu = u*u.Transpose();
      for(unsigned ii=0;ii<my_hinv.get_rows();ii++){
        for(unsigned jj=0;jj<my_hinv.get_columns();jj++){
           my_hinv(ii,jj) += einv*uu(ii,jj);
        }
      }
    }
  }
  return my_hinv;
  
}

std::vector<std::vector<Cartesian> > NormalModeAnalysis::GetDisplacements(const int nsteps) const{

  std::vector<std::vector<Cartesian> > displacements;
  matrix evecs = eigen[0];
  matrix evals = eigen[1];

  if(_type==MMUT_GNM){
    std::cerr << "Cannot compute displacements for Gaussian Network Model, try Anisotropic method\n";
  }

  displacements = std::vector<std::vector<Cartesian> > (evecs.get_columns()-6);

  for(unsigned k=6;k<evecs.get_columns();k++){
    double t = 0.0; // And now we have to loop over time.
    double freq = sqrt(evals(k,0));
    double time_per_cycle = 1.0/freq;
    double tstep = time_per_cycle/nsteps;
    for(int step=0;step<nsteps;step++,t+=tstep){
      for(unsigned i=0;i<evecs.get_rows();i+=3){
         double x = evecs(i,k)*sin(2*M_PI*freq*t);
         double y = evecs(i+1,k)*sin(2*M_PI*freq*t);
         double z = evecs(i+2,k)*sin(2*M_PI*freq*t);
         displacements[k-6].push_back(Cartesian(x,y,z));
      }
    }
  }

  return displacements;

}
