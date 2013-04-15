/*
     mmut/mmut_sasarea.h: CCP4MG Molecular Graphics Program
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



#ifndef __MMUT_SASArea__
#define __MMUT_SASArea__

#include <string>

#include "mmut_manager.h"
#include "mman_base.h"

enum SASMethods { SURFACEPOINTS , WODAKJANIN , LEERICHARDS };

DefineClass(CSASArea);

class CSASArea : public CMMANBase {

public :

  CSASArea(PCMMUTManager molHndin);    
  CSASArea(PCMMUTManager molHndin, int selHndin);
  CSASArea(PCMMUTManager molHndin, int selHndin,PCMMUTManager molHndin1, int selHndin1);
 ~CSASArea();
  void InitParams();
  void SetParams(int nv, double *value, int niv, int *ivalue);
  int SetMethod ( int meth );
  int Calculate_Contact ( void );
  int Calculate (int imodel=0,bool separate_models=0);
  int Calculate0 (int imodel);
  void  SetUDD(int atomUDD,int resUDD);
  std::string Print (int imodel, std::string selection1="", int set_selHnd1=-1,std::string selection2="", int set_selHnd2=-1 );

private:

  //the algorithm parameters parameters
  //PCMMANManager molHnd;
  int method;
  realtype HOHrad;
  realtype brick_margin;   
  realtype point_density;
  int exclude_solvent;
  AltLoc keep_altLoc;

  // the derived data
  int atomUDDHnd;
  int resUDDHnd;
  realtype total_area;

  int LeeAndRichards( realtype *radwithhoh,int imodel, int imode=0, int local_selHnd=-1 );
  int WodakAndJanin( realtype *radwithhoh,int imodel);
  int SurfacePoints( realtype *radwithhoh,int imodel);
  int ResSASArea (int imodel);
  //int GetAtomTypes(int natoms, PPCAtom selected_atoms, 
    //	int iatom_types[], int iatom_types_lookup[] );
  void  SortAB(realtype *arci, realtype *arcf, int &karc);
  void ArcLap(realtype *arci, realtype *arcf, int &karc);

};
#endif
