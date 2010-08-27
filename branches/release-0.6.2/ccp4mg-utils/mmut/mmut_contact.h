/*
     mmut/mmut_contact.h: CCP4MG Molecular Graphics Program
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



#ifndef __MMUT_Contact__
#define __MMUT_Contact__

#include <utility>
#include <string>
#include <mmut_manager.h>
#include <mman_base.h>
#include <mmut_connectivity.h>


DefineClass(CContact);

class CContact : public CMMANBase {

public :
  
  CContact(PCMMUTManager molHndin, int selHndin);
  CContact(PCMMUTManager molHndin);    
  CContact(PCMMUTManager molHndin, int selHndin,
          PCMMUTManager molHndin2, int selHndin2);
 ~CContact();
  void InitParams();
  void SetParams(int v, double *value, int niv, int *ivalue);
  int Calculate (bool separate_models = 1);
  int Calculate0 (int model=0);

  std::string Print (bool geometry=1);
  Connectivity2 close_contacts;

private:

  //the algorithm parameters parameters
  
  // Params
  int test_VDW_radius;
  int label_VDW_radius;
  int exclude_hbondable;
  realtype simple_max_cutoff;
  realtype simple_min_cutoff;
  realtype VDW_fraction_min;
  realtype VDW_fraction_max;

};
#endif
