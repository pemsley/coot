/*
     mmut/mmut_hbond.h: CCP4MG Molecular Graphics Program
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



#ifndef __MMUT_HBond__
#define __MMUT_HBond__

#include <utility>
#include <string>
#include <mmut_manager.h>
#include <mman_base.h>
#include <mmut_connectivity.h>

enum { MIN_D_A, MAX_D_A, MAX_H_A, MIN_DD_D_A, MIN_D_A_AA, MIN_H_A_AA, MIN_D_H_A };

DefineClass(CHBond);

class CHBond : public CMMANBase {

public :

  CHBond(PCMMUTManager molHndin, int selHndin);
  CHBond(PCMMUTManager molHndin);    
  CHBond(PCMMUTManager molHndin, int selHndin,
          PCMMUTManager molHndin2, int selHndin2);
 ~CHBond();
  void InitParams();
  void SetParams(int nv, double *value);
  int Calculate (bool separate_models=0 );
  int Calculate0 (int model=0);
  int LoadUDDHBType(PCMMUTManager molH); 

  std::string Print (bool geometry=1);
  Connectivity2 hbonds;

private:

  //the algorithm parameters parameters
  
  // Params
  realtype min_D_A;
  realtype max_D_A;
  realtype max_H_A;
  realtype min_DD_D_A;
  realtype min_D_A_AA;
  realtype min_H_A_AA;
  realtype min_D_H_A;

};

  
#endif
