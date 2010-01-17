/*
     mgapp/mgapp_base.h: CCP4MG Molecular Graphics Program
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


#ifndef __MgApp_Base__
#define __MgApp_Base__

#include <mman_manager.h>

DefineClass(CMgAppBase)
DefineStreamFunctions(CMgAppBase)

class CMgAppBase {

 public:
  CMgAppBase (const PCMMANManager molHndin );
  CMgAppBase (const PCMMANManager molHndin, const int selHndin );
  ~CMgAppBase ();

  void SetSelHandle ( const int selHndin);
  int GetSelection (PPCAtom &atomTable, int & nAtoms );
  int GetSelection (PPCResidue &resTable, int & nRes );
  void ClearSelection ();

 protected:
  PCMMANManager molHnd;
  int selHnd, resSelHnd;


};

#endif
