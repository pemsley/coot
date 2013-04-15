/*
     mgapp/mgapp_base.cc: CCP4MG Molecular Graphics Program
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


#include <string.h>
#include <mman_manager.h>
#include "mgapp_base.h"



//  ===================CMMgAppBase ============================================

//  A base class for all analysis classes

//------------------------------------------------------------------
CMgAppBase::CMgAppBase (const PCMMANManager molHndin) {
//------------------------------------------------------------------
  //Input data
  //printf ("CMgAppBase constructor\n");
  molHnd = molHndin;
  selHnd = -1;
  resSelHnd = -1;
}

//------------------------------------------------------------------
CMgAppBase::CMgAppBase (const PCMMANManager molHndin, const int selHndin ) {
//------------------------------------------------------------------
  //Input data
  molHnd = molHndin;
  selHnd = selHndin;
  resSelHnd = -1; 
}

//------------------------------------------------------------------
CMgAppBase::~CMgAppBase () {
//------------------------------------------------------------------
  if ( resSelHnd >= 0 ) molHnd->DeleteSelection ( resSelHnd );

}

//------------------------------------------------------------------
void CMgAppBase::SetSelHandle (const int selHndin ) {
//------------------------------------------------------------------
  //printf("SetSelHandle selHnd %i selHndin %i\n");
  selHnd = selHndin;
  if ( resSelHnd >= 0 ) {
    molHnd->DeleteSelection ( resSelHnd );
    resSelHnd = -1;
  }
}

//------------------------------------------------------------------
int CMgAppBase::GetSelection ( PPCAtom &atomTable , int &nAtoms ) {
//------------------------------------------------------------------

  //na = molHnd->GetNumberOfAtoms();
  if ( selHnd >= 0 )
    molHnd->GetSelIndex ( selHnd, atomTable, nAtoms );
  else {
    atomTable = NULL;
    molHnd->GetAtomTable ( atomTable, nAtoms );
  }
 //printf("nAtoms %i\n",nAtoms);
  return 0;
    
}
//------------------------------------------------------------------
int CMgAppBase::GetSelection ( PPCResidue &resTable , int &nRes ) {
//------------------------------------------------------------------
  //printf("selHnd,resSelHnd %i %i\n",selHnd,resSelHnd);
  nRes = 0;
  resTable = NULL;
  if ( resSelHnd >= 0 ) 
    molHnd->GetSelIndex ( resSelHnd, resTable, nRes );
  else if ( selHnd >= 0 ) {
    resSelHnd = molHnd->NewSelection();
    molHnd->Select(resSelHnd, STYPE_RESIDUE, selHnd, SKEY_NEW);
    molHnd->GetSelIndex ( resSelHnd, resTable, nRes );
  }
  else {
    //molHnd->GetResidueTable ( "/*", resTable, nRes );
    
    resSelHnd = molHnd->NewSelection();
    molHnd->Select(resSelHnd,STYPE_RESIDUE, 0, "*",ANY_RES,"*", ANY_RES,
                   "*","*","*","*","*",SKEY_OR);
    molHnd->GetSelIndex ( resSelHnd, resTable, nRes ); 
  }
  return 0;
    
}

//-----------------------------------------------------------------
void CMgAppBase::ClearSelection() {
//-----------------------------------------------------------------

  if ( resSelHnd >= 0 ) {
    molHnd->DeleteSelection ( resSelHnd );
    resSelHnd = -1;
  }
}
