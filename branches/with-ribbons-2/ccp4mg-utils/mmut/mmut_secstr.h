/*
     mmut/mmut_secstr.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2005 University of York, CCLRC

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


#ifndef __MMUT_SecStr__
#define __MMUT_SecStr__

#include <mmdb/mmdb_manager.h>
#include "mmut_manager.h"
#include "mman_base.h"
#include <string>

// The recognised secondary structure types
enum enum_SecStr { NOSECSTR, BETA, BULGE, TURN3, TURN4, TURN5, ALPHA };


DefineClass(CSecStructure);
DefineStreamFunctions(CSecStructure);

class CSecStructure : public CMMANBase {

public :
    
  // Constructors
  CSecStructure( PCMMUTManager molHnd );
  CSecStructure(const PCMMUTManager molHnd, const int selHnd );


 // Destructor
 ~CSecStructure();
 
  int GetSecondaryStructure ( int &nres, ivector &secstrout,
            imatrix &hbondsout, int imodel=0 );
  void SetParams (int nv,double *value, int niv, int *ivalue);
  int **GetHBonds (int imodel = 0);
  int *GetSecStr (int imodel = 0); 
  int SetFlagBulge ( int flag );
  std::string Print(int imodel = 0);
  void ClearMemory();
  PPCAtom* GetHBondAtoms(int imodel=0);

 private:

  float NOmaxdist;   
  float NOmaxdist2;
  float NOCanglemin;
  int flagBulge;

  int nRes;
  imatrix hbonds;
  ivector secstr;
  PPCAtom *hbond_atoms;
     
  int hbondsN;  //The first element dimension of hbonds - required for freeing memory
  void InitParams();
  int CalculateSecondaryStructure(int imodel = 0);
  int InitMemory( int nRes );
  Boolean IsHBond ( PCResidue PCRes1, PCResidue PCRes );
};
#endif
