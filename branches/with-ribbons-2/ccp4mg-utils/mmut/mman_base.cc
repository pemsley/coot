/*
     mmut/mman_base.cc: CCP4MG Molecular Graphics Program
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


#if defined(__sgi) || defined(sgi) || defined(__OSF1__) || defined(__osf__)
#include <string.h>
#else
#include <cstring>
#endif
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
#include <mmdb/mmdb_manager.h>
#include "mman_base.h"

using namespace std;


//  ===================CMMANBase ============================================

//  A base class for all analysis classes


//------------------------------------------------------------------
CMMANBase::CMMANBase (const PCMMUTManager molHndin, const int selHndin,
             const PCMMUTManager molHndin2, const int selHndin2 ) {
//------------------------------------------------------------------
  //Input data
  //cout << "molHndin " << molHndin << endl;
  for ( int i = 0; i < MAXSETS; i++) {
    selHnds[i] = -1;
    resSelHnds[i] = -1;
    nmrResSelHnds[i] = -1;
    nmrSelHnds[i] = -1;
    molHnds[i] = NULL;
  }

  exclude_solvent=0;
  exclude_hydrogen=0;
  exclude_alternate=0;
  own_selHnds = false;
  use_altLoc = "";

  molHnds[0] = molHndin;
  selHnds[0] = selHndin;
  if (molHndin2 != NULL ) molHnds[1] = molHndin2;
  selHnds[1] = selHndin2;
  //cout << "CMMANBase molHnd " << molHnds[0] << " " <<molHnds[1] << " " << selHnds[0] << "  " <<  selHnds[1] << endl;
}

//------------------------------------------------------------------
CMMANBase::~CMMANBase () {
//------------------------------------------------------------------
  // Clean up any sel handles created by this class
  //cout << "CMMANBase destructor";
  ClearSelection();
}

//-----------------------------------------------------------------------
void CMMANBase::SetExclusions(const int ex_solvent, const int ex_hydrogen,
  const int ex_alternate, const char* use_al) {
//-----------------------------------------------------------------------
  int save_selHnd;

  // If called with SetExclusions(1,1,1) this method apparently
  // does not work.  But it does work with the following diagnostic
  // line. I don't understand
  //cout << "SetExclusions " << ex_hydrogen << endl;

  exclude_solvent = ex_solvent;
  exclude_hydrogen =ex_hydrogen;
  exclude_alternate = ex_alternate;
 
  if ( exclude_solvent + exclude_hydrogen + exclude_alternate <= 0 ) {
    own_selHnds = false;
  }
  else {
    own_selHnds = true;
    // If we already have some selection handles set then 
    // reset them with the exclusions applied
    for (int is = 0; is< MAXSETS; is++) {
      if (selHnds[is]>0) {
        save_selHnd = selHnds[is];
        // Clear selections but BEWARE  do not want to clear
        // the original selHnd
        ClearSelection(is,0);
        SetSelHandle(is, save_selHnd );
      }
    }
  }   
}
//-----------------------------------------------------------------
int CMMANBase::SetMolHnd(const int iset,const PCMMUTManager molH) {
//-----------------------------------------------------------------
  if ( iset < 0 || iset >= MAXSETS ) {
    return 1;
  }
  else {
    molHnds[iset] = molH;
    return 0;
  }
}

//-----------------------------------------------------------------
PCMMUTManager CMMANBase::GetMolHnd(const int iset) {
//-----------------------------------------------------------------
  if (iset >= 0 && iset < MAXSETS) {
    return molHnds[iset];
  }
  else {
    return NULL;
  }
}

//------------------------------------------------------------------
int CMMANBase::SetSelHandle (const int iset, const int selHndin ,
                     const PCMMUTManager molHndin ) {
//------------------------------------------------------------------
  int is=0;
  int tmpHnd,tmpNAtoms;
  PPCAtom tmpAtoms;
  AltLoc aL;

  //cout << "CMMANBase::SetSelHandle " << own_selHnds << endl; 

  if ( iset < 0 || iset >= MAXSETS) return 1;

  // Delete existing resSelHnds etc - using original molHnds
  ClearSelection(iset);

  if (molHndin != NULL ) {
    SetMolHnd(iset,molHndin);
  }

  // If we do not have a molHnd for this selection set then
  // get the molHnd from a preceding set
  else if (molHnds[iset] == NULL) { 
    is = iset - 1;
    while ( is >= 0 ) {
      if ( GetMolHnd(is) != NULL ) { 
        SetMolHnd(iset, GetMolHnd(is));
      }
      else
        is = is-1;
    }
  }
  if  (is<0) return 2;

  // If own_selHnds is set then this class requires some
  // excusions from the input seelction set
  if ( own_selHnds ) {
    selHnds[iset] = molHnds[iset]->NewSelection();
    if (selHndin >= 0 ) 
      molHnds[iset]->Select(selHnds[iset],STYPE_ATOM,selHndin,SKEY_NEW);
    else
      // Select all atoms
      molHnds[iset]->SelectAtoms(selHnds[iset],0,0,SKEY_NEW);

    if (exclude_solvent>0) {
      molHnds[iset]->Select(selHnds[iset],STYPE_ATOM,
               "(HOH,H2O,WAT,SOL)",SKEY_CLR);
    }
    if ( exclude_hydrogen > 0 ) {
      //Exclude hydrogen 
      molHnds[iset]->Select(selHnds[iset],STYPE_ATOM,"/*/*/*/[H]",SKEY_CLR);
    }

    if ( exclude_alternate ) {
      // Alternate conformations - ??
      tmpHnd = molHnds[iset]->NewSelection();
      molHnds[iset]->Select(tmpHnd,STYPE_ATOM,selHnds[iset],SKEY_NEW);
      molHnds[iset]->Select(tmpHnd,STYPE_ATOM,"/*/*/*/*:",SKEY_CLR);
      molHnds[iset]->GetSelIndex ( tmpHnd,tmpAtoms, tmpNAtoms );
      while (tmpNAtoms > 0) {
        strncpy(aL,tmpAtoms[0]->altLoc,20);
        //cout << "tmpNAtoms " << tmpNAtoms << " *" << aL << "*" << endl;
        molHnds[iset]->Select(tmpHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,
		    "*","*","*","*",aL,SKEY_CLR);
        molHnds[iset]->GetSelIndex ( tmpHnd,tmpAtoms, tmpNAtoms );
        if ( strcmp(aL,use_altLoc) != 0 ) {
          //cout << "Removing " << aL << endl; 
          molHnds[iset]->Select(selHnds[iset],STYPE_ATOM,0,"*",ANY_RES,
             "*",ANY_RES, "*","*","*","*",aL,SKEY_CLR);
        }
      }
      molHnds[iset]->DeleteSelection ( tmpHnd );
    }
  }
  else {
    selHnds[iset] = selHndin;
  }
  //cout << "SetSelHandle selHnds " << selHnds[0] << " " << selHnds[1] << endl;
  return 0;
}

//------------------------------------------------------------------
int CMMANBase::SetSelHandle (const int selHndin ) {
//------------------------------------------------------------------
    return SetSelHandle (0, selHndin );
}


//--------------------------------------------------------------------
int  CMMANBase::GetOneModel(const int iset, const int selHin, 
          PPCAtom &atomTable , int &nAtoms, const int model) {
//--------------------------------------------------------------------
  // If model > 0 return the selected atoms in the appropriate NMR model
  int selH;
  if ( selHin > 0 ) 
    selH = selHin;
  else
    selH = selHnds[iset];
        
  //cout << "GetOneModel model " << model << " " << molHnds[iset]->GetNumberOfModels() << "selH" << selH << endl;

  if ( model > 0 && model <= molHnds[iset]->GetNumberOfModels() ) { 
    if (nmrSelHnds[iset] > 0) molHnds[iset]->DeleteSelection(nmrSelHnds[iset]); 
    nmrSelHnds[iset] = molHnds[iset]->NewSelection();
    molHnds[iset]->Select (nmrSelHnds[iset], STYPE_ATOM,selH,SKEY_OR);
    //molHnds[iset]->GetSelIndex ( nmrSelHnds[iset], atomTable, nAtoms );
    //cout << "GetOneModel initial model " << model << " nAtoms " << nAtoms << endl;
    molHnds[iset]->Select (nmrSelHnds[iset],STYPE_ATOM,model,"*",ANY_RES,"*",ANY_RES,
		    "*","*","*","*","*",SKEY_AND);
    molHnds[iset]->GetSelIndex ( nmrSelHnds[iset], atomTable, nAtoms );
    //cout << "GetOneModel model " << model << " nAtoms " << nAtoms << endl;
  }
  else {
    molHnds[iset]->GetSelIndex ( selH, atomTable, nAtoms );
  }
  return 0;
  
}
//--------------------------------------------------------------------
int  CMMANBase::GetOneModel(const int iset, const int selHin, 
     PPCResidue &resTable , int &nRes, const int model) {
//--------------------------------------------------------------------
  int nMod;
  int selH;

  //cout << "GetOneModel model " << model << " " << molHnds[iset]->GetNumberOfModels() << "selH" << selH << endl;

  if ( selHin > 0 ) 
    selH = selHin;
  else
    selH = selHnds[iset];

  nMod = molHnds[iset]->GetNumberOfModels();
  //cout << "nMod " << nMod << endl;
  if ( model > 0 && model <= molHnds[iset]->GetNumberOfModels() ) { 
    if (nmrResSelHnds[iset] >= 0) molHnds[iset]->DeleteSelection(nmrResSelHnds[iset]);
    nmrResSelHnds[iset] = molHnds[iset]->NewSelection();
    molHnds[iset]->Select (nmrResSelHnds[iset], STYPE_RESIDUE,selH,SKEY_OR); 
    molHnds[iset]->Select (nmrResSelHnds[iset],STYPE_RESIDUE,model,"*",ANY_RES,"*",
		    ANY_RES,"*","*","*","*","*",SKEY_AND);
    molHnds[iset]->GetSelIndex ( nmrResSelHnds[iset], resTable, nRes );
  }
  else {
    molHnds[iset]->GetSelIndex ( selH, resTable, nRes );
  }
  //cout << " first res" << resTable[0]->seqNum << " " << resTable[1]->seqNum << endl;

  return 0;
}

//------------------------------------------------------------------
int CMMANBase::GetSelection ( PPCAtom &atomTable , int &nAtoms,
                       const int model ) {
//------------------------------------------------------------------
  return GetSelection (0,atomTable,nAtoms, model);
}

//------------------------------------------------------------------
int CMMANBase::GetSelection ( const int iset ,PPCAtom &atomTable ,
                      int &nAtoms, const int model ) {
//------------------------------------------------------------------
  nAtoms = 0;
  if  (selHnds[iset] >= 0 ) {
     GetOneModel( iset ,selHnds[iset], atomTable, nAtoms, model );
     //cout << "GetSelection from GetOneModel " <<selHnds[iset] << " " << atomTable[1] << endl;
  }
  else if ( iset == 0 || selMode[iset] == SELECT_ALL ) {
    selHnds[iset] = molHnds[iset]->NewSelection();
    molHnds[iset]->Select (selHnds[iset],STYPE_ATOM,0,"*",ANY_RES,"*",
		ANY_RES,"*","*","*","*","*",SKEY_OR);
    GetOneModel(iset,selHnds[iset],atomTable,nAtoms,model);
    //cout << "GetSelection new selHnd  " << selHnds[iset] << atomTable[1] << endl;
  } else if ( selHnds[iset-1] <= 0 || molHnds[iset] != molHnds[iset-1] ) 
    return 1;
  else {
    if (selMode[iset] == SELECT_SAME )
      // Want the same as previous set
      GetOneModel(iset,selHnds[iset-1],atomTable,nAtoms,model);
    else if (selMode[iset] == SELECT_NOT ) {
      // Want what is excluded from previous set
      selHnds[iset] = molHnds[iset]->NewSelection();
      molHnds[iset]->Select(selHnds[iset],STYPE_ATOM,selHnds[iset-1],SKEY_XOR);
      GetOneModel( iset ,selHnds[iset], atomTable, nAtoms, model );
    } 
  }

  return 0;
}

//------------------------------------------------------------------
int CMMANBase::GetSelection ( PPCResidue &resTable, int &nRes,
                      const int model ) {
//------------------------------------------------------------------
  return GetSelection(0,resTable,nRes, model);
}

//------------------------------------------------------------------
int CMMANBase::GetSelection ( const int iset, PPCResidue &resTable,
                          int &nRes, const int model ) {
//------------------------------------------------------------------
  nRes = 0;

  if ( resSelHnds[iset] >= 0 ) {
    // The residue selection handle is already defined 
    GetOneModel( iset,resSelHnds[iset], resTable, nRes, model );
    //cout << "from GetOneModel resSelH " << resSelH << " nRes " << nRes << endl;
  }
  else if ( selHnds[iset] >= 0 ) {
    // The atom selection handle is already defined 
    resSelHnds[iset] = molHnds[iset]->NewSelection();
    molHnds[iset]->Select(resSelHnds[iset], STYPE_RESIDUE, selHnds[iset], SKEY_OR);
    GetOneModel( iset ,resSelHnds[iset], resTable, nRes, model );     
  }
  // Return all residues by default for first set or if previous
  // set is undefined
  else if ( iset == 0 || selMode[iset] == SELECT_ALL || selHnds[iset-1] < 0 ) {
    resSelHnds[iset] = molHnds[iset]->NewSelection();
    molHnds[iset]->Select(resSelHnds[iset],STYPE_RESIDUE, 0, "*",ANY_RES,
       "*", ANY_RES, "*","*","*","*","*",SKEY_OR);
    GetOneModel( iset ,resSelHnds[iset], resTable, nRes, model );
  }
  else {
    // Is dependent on previous selection 
    if ( molHnds[iset] != GetMolHnd(iset-1) )return 1; 
      // If do not have resSelHnd for the first selection then create it
    // from the atom selection
    if ( resSelHnds[iset-1] < 0 ) {
      resSelHnds[iset-1] = molHnds[iset-1]->NewSelection();
      molHnds[iset]->Select(resSelHnds[iset-1], STYPE_RESIDUE, selHnds[iset-1], SKEY_NEW);
    }
    // The second set is either the same as or an exclusive or of the first
    if (selMode[iset] == SELECT_SAME ) 
       GetOneModel( iset ,resSelHnds[iset-1], resTable, nRes, model );
    else if (selMode[iset] == SELECT_NOT ) {
      molHnds[iset]->Select(resSelHnds[iset], STYPE_RESIDUE, resSelHnds[iset-1], SKEY_XOR)
;
      GetOneModel( iset ,resSelHnds[iset], resTable, nRes, model );
    } 
  }
  return 0;
    
}

//-----------------------------------------------------------------
void CMMANBase::ClearSelection(const int iset, const int clear_selHnd) {
//-----------------------------------------------------------------
  int ifst=0;
  int ilst=MAXSETS ;
  if ( iset >= 0 ) {
    ifst = iset;
    ilst = iset+1;
  } 
  for (int i =ifst;i<ilst;i++) {
    if (resSelHnds[i]>0) molHnds[i]->DeleteSelection(resSelHnds[i]);
    if (nmrResSelHnds[i]>0) molHnds[i]->DeleteSelection(nmrResSelHnds[i]);
    if (nmrSelHnds[i]>0) molHnds[i]->DeleteSelection(nmrSelHnds[i]);
    resSelHnds[i] = -1;
    nmrResSelHnds[i] = -1;
    nmrSelHnds[i]=-1;
    // Default clear_selHnd=1, if class hasits own selHnds then
    // they are deleted
    if ( own_selHnds && clear_selHnd ) {
      if(selHnds[i]>0) molHnds[i]->DeleteSelection(selHnds[i]);
      selHnds[i]=-1;
    }
  }
}

