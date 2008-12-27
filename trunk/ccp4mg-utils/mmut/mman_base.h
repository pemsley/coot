/*
     mmut/mman_base.h: CCP4MG Molecular Graphics Program
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


#ifndef __MMAN_Base__
#define __MMAN_Base__

#include <mmut_manager.h>

DefineClass(CMMANBase);
DefineStreamFunctions(CMMANBase);
#define MAXSETS 5
enum enum_SelMode { SELECT_ALL , SELECT_SAME , SELECT_NOT };

class CMMANBase {

 public:
  CMMANBase (const PCMMUTManager molHndin, const int selHndin=-1,
             const PCMMUTManager molHndin2=NULL, const int selHndin2=-1  );
  ~CMMANBase ();

  int SetSelHandle ( const int selHndin);
  int SetSelHandle ( const int iset, const int selHndin, const PCMMUTManager molHndin=NULL);
  int GetSelection (PPCAtom &atomTable, int & nAtoms, const int model = 0 );
  int GetSelection (const int iset, PPCAtom &atomTable, int & nAtoms ,const int model = 0);
  int GetSelection (const int iset,PPCResidue &resTable, int & nRes, const int model = 0 );
  int GetSelection ( PPCResidue &resTable, int & nRes, const int model = 0 );
  void ClearSelection (const int iset=-1, const int clear_selHnd=1);
  int  GetOneModel(const int iset,const int selH, PPCAtom &atomTable ,
                        int &nAtoms, const int model);
  int  GetOneModel(const int iset,const int selH, PPCResidue &resTable ,
                         int &nRes, const int model);
  void SetExclusions(const int ex_solvent, const int ex_hydrogen, const int ex_alternate, const char* use_al="");

  PCMMUTManager GetMolHnd(const int iset=0);
  int  SetMolHnd(const int iset, const PCMMUTManager molHndin);

 protected:


  // Flags to exclude nasties that break analysis code
  int exclude_solvent;
  int exclude_hydrogen;
  int exclude_alternate;
  char* use_altLoc;
  Boolean own_selHnds;
  int selMode[MAXSETS];
  PCMMUTManager molHnds[MAXSETS];
  int selHnds[MAXSETS];   // Selection handles passed in from application
  int resSelHnds[MAXSETS]; // Selection handles managed by this class
  int nmrSelHnds[MAXSETS];
  int nmrResSelHnds[MAXSETS];
};

#endif
