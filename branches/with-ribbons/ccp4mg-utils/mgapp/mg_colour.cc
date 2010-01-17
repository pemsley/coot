/*
     mgapp/mg_colour.cc: CCP4MG Molecular Graphics Program
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


#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <mman_manager.h>
#include <mg_colour.h>
#include <mgapp_base.h>
#include <string>
#include <rgbreps.h>
#include <vector>
#include "atom_util.h"


//*********************************************************************
// MOLECULE COLOUR CLASS
//************************************************************************

  /*! In context of molecular graphics one instance of CMolColour is
     instantiated for each molecule object.
  \param molHndin A pointer to the MMDB manager - i.e. the molecule
  \param selHndin A handle for the MMDB selection - i.e. the molecule object. See SetSelHandle.
  \param scheme A pointer to the CColourSchemes class which has the definition of the colour schemes.  There is only one instance of this class created by the MG. 
  */
//-----------------------------------------------------------------------
CMolColour::CMolColour(PCMMANManager molHndin, int selHndin,
   PCColourSchemes schemein ) : CMgAppBase ( molHndin , selHndin ) {
//-----------------------------------------------------------------------
  //data  - needs initiallising BEFORE SetSelHandle
  reapply = 1;
  if_continuous = 0;
  rules_colours = NULL;
  rules_selHnds = NULL;
  atom_colour_vector = NULL;
  //input data
  colourSchemes = schemein;

  //default params
  oneColour = 1;
  multColours = 0;
  mode = BYATOMTYPE;
  non_C_mode = BYATOMTYPE;

  resUDDHnd = -1;
  atomUDDHnd = -1;

}

//----------------------------------------------------------------------
CMolColour::~CMolColour(){
//----------------------------------------------------------------------
  if (rules_colours ) {
    FreeVectorMemory (rules_colours,0);
    rules_colours = NULL;
  }
  if (rules_selHnds ) {
    FreeVectorMemory (rules_selHnds,0);
    rules_selHnds = NULL;
  }
}

/*! Set the colour to be applied if the colouring mode is ONECOLOUR
    \param icol An integer value which sould be a valid colour code as defined in the CColour class 
    */ 

//----------------------------------------------------------------------
int CMolColour::SetOneColour (int icol) {
//----------------------------------------------------------------------
  oneColour = icol;
  if_continuous = 0;
  if (mode == ONECOLOUR) reapply = 1;
  return 0;
}

/*! Set the colouring mode to be applied. The valid values for the mode are
    - ONECOLOUR  Apply one colour as given by SetOneColour
    - BYATOMTYPE Colour by the atom element type
    - BYRESTYPE  Colour by residue type
    - BYCHAIN    Colour each chain a different colour
    - THRUCHAIN  Colour by position in chain
    - BLENDTHRU      Blend colours thru chains by user customised criteria
    - SECSTR     Colour by secondary structure
    - BVALUE     Colour by B value
    - OCCUPANCY  Colour by atom occupancy
    - CHARGE     Colour by atomic charge
    - SASRES     Colour by residue SAS area
    - SASATOM    Colour by atom SAS area
    - CONTACTRES     Colour by residue contact area
    - CONTACTATOM    Colour by atom contact area
   The detailed interpretation of each of these modes is defined within the CMolColourScheme class
    \param mod Colouring mode   
*/


//---------------------------------------------------------------------
void CMolColour::SetSelHandle ( const int selHndin) {
//---------------------------------------------------------------------
  reapply = 1;
  CMgAppBase::SetSelHandle(selHndin);
}

//----------------------------------------------------------------------
int CMolColour::SetMode ( int force, int mod , int non_C_mod, int icol,
                  int ncols, int icols[]) {
//----------------------------------------------------------------------
  //std::cout << "SetMode " << mod << std::endl;
  if (force ) reapply = 1;
  if ( mod>=0 && mod <=NCOLOURMODES && mod != mode ) {
    mode = mod;
    reapply = 1;
  }
    
  if (icol>=0) {
    oneColour = icol;
    if_continuous = 0;
    if (mode == ONECOLOUR) reapply = 1;
  }
  if (icols) {
    multColours = icols;
    nMultColours = ncols;
    if (mode == BYCHAIN || mode == BYMODEL) reapply = 1;
  }
  //cout << "nMultColours " << nMultColours << " multColours ";
  //for ( int i=0; i<ncols; i++ ) cout << multColours[i] << " "; 
  //cout << endl;

  if ( mode == BVALUE || mode == OCCUPANCY || mode == CHARGE ||
       mode == SASRES || mode == SASATOM ||
       mode == CONTACTRES || mode == CONTACTATOM ||
       mode == THRU_CHAIN || mode == BLENDTHRU ) 
    if_continuous = 1;
  else 
    if_continuous = 0;

  if (non_C_mod != non_C_mode) {
    non_C_mode = non_C_mod;
    if (if_continuous!=1) reapply = 1;
  }
    
  //std::cout << "Colour mode " << mode << std::endl;
  return 0;
}

/*! Return the per atom colour. The values of the vector icol correspond to a colour code defined in CColour. The number of elements in the vector is the number of atoms in the selection. Do NOT delete this vector.
  \param nat Number of atoms in the selection and therefore the number of elements in the vector icol.
  \param icol A pointer to an integer vector.
*/


//-----------------------------------------------------------------------
AtomColourVector *CMolColour::GetAtomColourVector() {
//----------------------------------------------------------------------
  int RC;
  //std::cout << "GetAtomColourVector reapply " << if_continuous << " " << reapply << std::endl;
  if ( reapply || !atom_colour_vector ) RC = ReColour();
  if (!RC) {
    return atom_colour_vector;
  } else {
    return NULL;
  }
}
 


 /*! Write a list of the atom id and colour to standard out.  
    This is primarilly for diagnostic purposes
 */


//-----------------------------------------------------------------------
int CMolColour::Print( ) {
//-----------------------------------------------------------------------
  int i,nat;
  PPCAtom selAtoms;
  int *colour;
  int RC;

  if ( reapply || !atom_colour_vector ) {
    RC =  ReColour();
    if(RC) return RC;
  }
  GetSelection( selAtoms, nat );
  printf("Print nat %i\n",nat);
  colour = new int[nat];
  for ( i = 0; i < nat; i++ ) {
       printf ( "colour %i %s %s %i\n",i,selAtoms[i]->residue->name,
         selAtoms[i]->name, colour[i]);
  }
  delete [] colour;
  return 0;
}


//-------------------------------------------------------------------------
int CMolColour::SetScheme ( int nrules, int def_colour, int rule_selHnd[],
                             int rule_colour[] , int mode_in ) {
//-------------------------------------------------------------------------
  // Set the list of selHnd and colours for user defined rules
  if (rules_colours ) {
    FreeVectorMemory (rules_colours,0);
    rules_colours = NULL;
  }
  if (rules_selHnds ) {
    FreeVectorMemory (rules_selHnds,0);
    rules_selHnds = NULL;
  }

  GetVectorMemory ( rules_colours, nrules, 0);
  GetVectorMemory ( rules_selHnds, nrules, 0);

  n_rules = nrules;
  rules_default_colour = def_colour;
  rules_mode = mode_in;

  //cout << "SetScheme n_rules " << n_rules << endl;
  for ( int  i = 0; i < nrules; i++ ) {
    rules_colours[i] = rule_colour[i];
    rules_selHnds[i] = rule_selHnd[i];
  }
  if_continuous = 0;
  mode = USER_SCHEME;
  reapply = 1;
  return 0;

}


//--------------------------------------------------------------------
// PRIVATE METHODS
//--------------------------------------------------------------------


/* Recalculate the atom colour vector
 */

//----------------------------------------------------------------------
int CMolColour::ReColour() {
//----------------------------------------------------------------------
  int RC;
  //std::cout << "ReColour mode " << mode << std::endl;
  if ( atom_colour_vector ) delete atom_colour_vector;
  atom_colour_vector = new AtomColourVector();
  switch ( mode ) {
  case ONECOLOUR:
    RC = OneColour();
    break;
  case BYRESTYPE:
    RC = ByResidueType();
    break;
  case BYCHAIN:
    RC = ByChain();
    break;
  case BYMODEL:
    RC = ByModel();
    break;
  case SECSTR:
    RC = BySecondaryStructure();
    break;
  case MAINSIDE:
    RC = ByMainSide();
    break;
  case BYATOMTYPE:
    RC = ByAtomType();
    break;
  case SASRES:
    RC = BySASArea(1);
    break;
  case SASATOM:
    RC = BySASArea(0);
    break;
  case CONTACTRES:
    RC = BySASArea(1,1);
    break;
  case CONTACTATOM:
    RC = BySASArea(0,1);
    break;
  case BVALUE:
    RC = ByBValue();
    break;
  case OCCUPANCY:
    RC = ByOccupancy();
    break;
  case CHARGE:
    RC = ByCharge();
    break;
  case USER_SCHEME:
    RC = UserScheme();
    break;
  case THRU_CHAIN:
    RC = ByThruChain(-1);
    break;
  case BLENDTHRU:
    RC = BlendThruChain();
    break;
  default:
    RC = 1;
  }
  //std::cout << "non_C_mode " << non_C_mode << " " << BYATOMTYPE << std::endl; 
  //if (non_C_mode == BYATOMTYPE && mode != BYATOMTYPE) 
  //  RC = ByAtomType();
  //cout << "ReColour mode" << mode << RC << endl;
  reapply = 0;
  return RC;
}


//-----------------------------------------------------------------------
int CMolColour::OneColour( ) {
//-----------------------------------------------------------------------
  int nat,i;
  PPCAtom selAtoms;
  int *colour;
  
  GetSelection (selAtoms, nat);
  colour = new int[nat];
  
  for ( i = 0; i < nat; i++ ) { colour[i] = oneColour; }

  if(non_C_mode == BYATOMTYPE) ByAtomType ( 0, colour, nat, selAtoms);

  atom_colour_vector->SetAtomColours(nat,colour);
  delete [] colour;

  return 0;
}

//----------------------------------------------------------------------
int CMolColour::EditColour (int icol, const int selHndin) {
//----------------------------------------------------------------------
  return 0;
}

//-----------------------------------------------------------------------
int CMolColour::ByAtomType( ) {
//-----------------------------------------------------------------------
  int *colour;
  PPCAtom selAtoms;
  int nat;

  GetSelection (selAtoms, nat);
  colour = new int[nat];

  int rv = ByAtomType ( 1,colour, nat, selAtoms);
  if (rv == 0) atom_colour_vector->SetAtomColours(nat,colour);
  delete [] colour;
  return rv;

}
  
//-----------------------------------------------------------------------
int CMolColour::ByAtomType( int colour_carbon, int *colour , 
                       int nat , PPCAtom selAtoms) {
//-----------------------------------------------------------------------
  int i,j,jC=-1;

   // Find the C colour 
  if (colour_carbon) {
    for ( j = 0; j < colourSchemes->AtomType->nTypes; j++ ) {
      if (colourSchemes->AtomType->strtypes[j].compare(" C")==0) jC = j;
    }
  }

  for ( i = 0; i < nat; i++ ) {
    if ( strcmp(selAtoms[i]->element," C")== 0) {
      if (jC>= 0 ) colour[i] = colourSchemes->AtomType->codes[jC];
      // Otherwise we are skipping colouring the carbon
    } else {
      for ( j = 1; j < colourSchemes->AtomType->nTypes; j++ ) {
        if ( colourSchemes->AtomType->strtypes[j].compare(selAtoms[i]->element)==0) {
          colour[i] =colourSchemes->AtomType->codes[j];
          break;
        }
        colour[i] = colourSchemes->AtomType->defColour;
      }
    }
  }

  return 0;
}


//-----------------------------------------------------------------------
int CMolColour::ByResidueType( ) {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int nat,i,j;
  int *colour;
 
  GetSelection (selAtoms, nat);
  colour = new int[nat];
  
  for ( i = 0; i < nat; i++ ) {
    for ( j = 0; j < colourSchemes->ResType->nTypes; j++ ) {
      if ( colourSchemes->ResType->strtypes[j].compare(
                  selAtoms[i]->residue->name)==0) {
        colour[i] = colourSchemes->ResType->codes[j];
        //cout << selAtoms[i]->residue->name << " " <<  colour[i] << endl;
        break;
      }
      //std::cout << "*" << selAtoms[i]->residue->name << "*" <<  std::endl;
      colour[i] = colourSchemes->ResType->defColour;
    }
  }

  if(non_C_mode == BYATOMTYPE) ByAtomType ( 0, colour, nat, selAtoms);
  atom_colour_vector->SetAtomColours(nat,colour);
  delete [] colour;
  return 0;
}

//-----------------------------------------------------------------------
int CMolColour::ByChain( ) {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int nat;
  int i,j;
  int nChainTotal;
  PPCChain chainTable;
  int *colour;
  int *cols;

  GetSelection (selAtoms, nat);
  colour = new int[nat];

  molHnd->GetChainTable(1,chainTable,nChainTotal);
  cols = new int[nChainTotal];

  //cout << " cols ";
  for ( j = 0; j < nChainTotal ; j++ ) {
    cols[j] = j % nMultColours;
    //cout << cols[j] << " " <<  multColours[cols[j]] << "";
  }
  //cout << endl;

  for ( i = 0; i < nat; i++ ) {
    for ( j = 0; j < nChainTotal ; j++ ) {
      if  ( selAtoms[i]->GetChain() == chainTable[j] )  {
	colour[i] = multColours[cols[j]];
        break;
      }
      colour[i] = 0;
    } 
  }

  if(non_C_mode == BYATOMTYPE) ByAtomType ( 0, colour, nat, selAtoms);
  return atom_colour_vector->SetAtomColours(nat,colour);
}

//-----------------------------------------------------------------------
int CMolColour::ByModel( ) {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int nat;
  int i,j;
  int nModels;
  int *colour;
  int *cols;

  GetSelection (selAtoms, nat);
  colour = new int[nat];

  nModels = molHnd->GetNumberOfModels();
  //cout << "ByModel nModels " << nModels << endl;
  cols = new int[nModels];
  for ( j = 0; j < nModels ; j++ ) {
    cols[j] = j % nMultColours;
  }

  for ( i = 0; i < nat; i++ ) {
    colour[i] = multColours[cols[selAtoms[i]->GetModelNum()-1]];
  }
  if(non_C_mode == BYATOMTYPE) ByAtomType ( 0, colour, nat, selAtoms);
  int rv = atom_colour_vector->SetAtomColours(nat,colour);
  delete [] colour;
  delete [] cols;
  return rv;
}

//-----------------------------------------------------------------------
int CMolColour::BySecondaryStructure( ) {
//-----------------------------------------------------------------------
  int nat,is,i,j;
  PPCAtom selAtoms;
  int *colour;
   
  GetSelection (selAtoms, nat);
  colour = new int[nat];
 
  for ( i = 0; i < nat ; i++ ) {
    is = selAtoms[i]->residue->SSE;
    colour[i] =  colourSchemes->SecStr->defColour;
    //printf ( "is %i\n",is);
    for ( j = 0; j < colourSchemes->SecStr->nTypes; j++) {
      if ( is == colourSchemes->SecStr->itypes[j]) {
        colour[i] = colourSchemes->SecStr->codes[j];
        break;
      }
    }
  }
  if(non_C_mode == BYATOMTYPE) ByAtomType ( 0, colour, nat, selAtoms);
  int rv = atom_colour_vector->SetAtomColours(nat,colour);
  delete [] colour;
  return rv;
}

//-----------------------------------------------------------------------
int CMolColour::BySASArea(int res_mode, int ifcontact ) {
//-----------------------------------------------------------------------
  int nat,i,RC;
  PPCAtom selAtoms;
  realtype sas;
  CColourScheme *scheme;
  double *fcolour;

  GetSelection (selAtoms, nat);
  if (nat <= 0 ) return 1;
  fcolour = new double[nat];
  //std::cout << "BySASArea nat" << nat << " " <<atomUDDHnd << " " <<resUDDHnd  << std::endl;

  if ( atomUDDHnd <= 0 || resUDDHnd <= 0 ) return 1;
 
  if (res_mode ) {
    if (ifcontact) 
      scheme = colourSchemes->ResContact;
    else
      scheme = colourSchemes->ResSAS;
  } else {
    if (ifcontact) 
      scheme = colourSchemes->AtomContact;
    else
      scheme = colourSchemes->AtomSAS;
  }
  //std::cout << "scheme " << scheme << std::endl;
  
  for ( i = 0; i < nat ; i++ ) {
    if ( res_mode )
      RC = selAtoms[i]->residue->GetUDData(resUDDHnd,sas);
    else
      RC = selAtoms[i]->GetUDData(atomUDDHnd,sas);
    fcolour[i] = scheme->GetFColour(sas);
    //std::cout << "sas,fcolour " << sas << " " << fcolour[i] << std::endl;
  }

  int rv = atom_colour_vector->SetAtomColours(
     scheme->GetInterpolationMode(),nat,fcolour,scheme->GetSchemeCodes(),
     scheme->GetColourWheelDirection() );
  delete [] fcolour;
  return rv;

}



//-----------------------------------------------------------------------
int CMolColour::ByBValue () {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int nat,i;
  double *fcolour;
  
  GetSelection (selAtoms, nat);
  if (nat <= 0 ) return 0;
  fcolour = new double[nat];

  for ( i = 0; i < nat ; i++ ) {
    fcolour[i] =  colourSchemes->BValue->GetFColour(selAtoms[i]->tempFactor);
  }
  int rv = atom_colour_vector->SetAtomColours( 
               colourSchemes->BValue->GetInterpolationMode(), 
               nat,fcolour,colourSchemes->BValue->GetSchemeCodes(),
               colourSchemes->BValue->GetColourWheelDirection());
  delete fcolour;
  return rv;

}

//----------------------------------------------------------------------------
int CMolColour::ByOccupancy() {
//----------------------------------------------------------------------------
  PPCAtom selAtoms;
  int nat,i;
  double *fcolour;

  GetSelection (selAtoms, nat);
  if (nat <= 0 ) return 0;
  fcolour = new double[nat];

  
  for ( i = 0; i < nat ; i++ ) {
    fcolour[i] =  colourSchemes->Occupancy->GetFColour(selAtoms[i]->occupancy);
  }
  int rv = atom_colour_vector->SetAtomColours(
            colourSchemes->Occupancy->GetInterpolationMode(),
            nat,fcolour,colourSchemes->Occupancy->GetSchemeCodes(),
            colourSchemes->Occupancy->GetColourWheelDirection());
  delete [] fcolour;
  return rv;
}

//----------------------------------------------------------------------------
int CMolColour::ByCharge() {
//----------------------------------------------------------------------------
  PPCAtom selAtoms;
  int nat,i,udd;
  double *fcolour,charge;
  //char AtomID[50];
  
  GetSelection (selAtoms, nat);
  if (nat <= 0 ) return 0;
  fcolour = new double[nat];

  // Is there charge in UDData - if not then load it
  molHnd->LoadCharge("electrostatic_potential_charge");
  //molHnd->LoadCharge("partial_charge")
  udd = molHnd->GetUDDHandle( UDR_ATOM,"atom_charge");
  if (udd<0) return -1;

  for ( i = 0; i < nat ; i++ ) {
    selAtoms[i]->GetUDData(udd,charge);
    fcolour[i] =  colourSchemes->Charge->GetFColour(charge);
    /*
    if (i<50) {
      selAtoms[i]->GetAtomID(AtomID);
      std::cout << AtomID << " " << charge << " " << fcolour[i] << std::endl;
    }
    */
  }
  
  int rv = atom_colour_vector->SetAtomColours(
               colourSchemes->Charge->GetInterpolationMode(),
               nat,fcolour,colourSchemes->Charge->GetSchemeCodes(),
               colourSchemes->Charge->GetColourWheelDirection());
  delete [] fcolour;
  return rv;
}

//----------------------------------------------------------------------------
int CMolColour::ByThruChain (int chain_no , int ncolours, int *colour) {
//----------------------------------------------------------------------------
  PPCAtom selAtoms,atomTable;
  PPCResidue resTable;
  int nranges,nat,tot_nat,i,ir,first_chain,nchains,ic,nres,isel,ir_amino;
  realtype rr,incr;
  double *fcolour;

  GetSelection (selAtoms, nat);
  if (nat <= 0 ) return 0;
  fcolour = new double[nat];
  
  if ( chain_no >= 0 ) { 
    nchains = 1;
    first_chain = chain_no;
  }
  else {
    first_chain = 0;
    nchains =  molHnd->GetNumberOfChains(1);

  }

  nranges =  colourSchemes->ThruChain->nTypes - 3;
  isel = -1;
  for ( ic = 0; ic < molHnd->GetNumberOfChains(1) ; ic++ ) {
    incr = fmod (float(ic - first_chain),float(nranges));
    //cout << "chain,incr " << ic << " " << incr << endl;
    molHnd->GetResidueTable(1,ic,resTable,nres);
    ir_amino = 0;
    for (ir = 0; ir< nres; ir++ ) {
      if ( resTable[ir]->isAminoacid() || resTable[ir]->isNucleotide() ) {
        rr = incr + float(ir_amino++)/float(nres);
      }
      else {
        rr = -1;
      }
      //std::cout << "ByThruChain " << ir << " " << rr << std::endl; 
      resTable[ir]->GetAtomTable ( atomTable,tot_nat );
      for ( i = 0; i < tot_nat; i++ ) {
        if ( atomTable[i]->isInSelection(selHnd) ) 
	  fcolour[++isel]= rr;  
      }
    }
  }
  int rv = atom_colour_vector->SetAtomColours(HSVCOLOURREP,nat,fcolour,
                     colourSchemes->ThruChain->GetSchemeCodes());
  delete [] fcolour;
  return rv;
}

void CMolColour::SetColourWheelDirection(std::vector<int> direction) {
  colour_wheel_direction.clear();
  int ndir = direction.size();
  for (int n=0;n<ndir;n++) {
    colour_wheel_direction.push_back(direction[n]); }
}

//----------------------------------------------------------------------------
int CMolColour::BlendThruChain() {
//----------------------------------------------------------------------------
  PPCAtom selAtoms;
  PPCResidue selRes=NULL,resTable=NULL;
  int RC,nat,nres,i,j,k;
  int udd,resSelHnd;
  float fmin,fincr;
  realtype ff;
  byte SSE_mapping[9],last_SSE;
  int nofSSE;
  double *fcolour;


  RC = GetSelection (selAtoms, nat);
  fcolour = new double[nat];
  for ( i = 0; i < nat; i++ ) fcolour[i] = -1.0;

  //std::cout << "BlendThruChain rules_mode " << rules_mode << std::endl; 

  // Initialise a temporary residue UDD for colour code
  udd = molHnd->GetUDDHandle ( UDR_RESIDUE,"tmp_res_real" );
  if (udd <= 0 ) {
    udd = -1;
    udd = molHnd->RegisterUDReal ( UDR_RESIDUE,"tmp_res_real" );
    if (udd <= 0 ) return udd;
  }
  molHnd->GetResidueTable (resTable, nres);
  for (i=0;i<nres;i++) resTable[i]->PutUDData(udd,-1.0);
  delete [] resTable;

  resSelHnd = molHnd->NewSelection();

  switch (rules_mode) {

  case BLEND_RESIDUE: 
  // Assign appropriate colour code to each residue - save in temporary UDD
    i = 1;
    fmin = 0.0;
    
    while (i < n_rules-1) {
      j = i+1;
      while ( j+1 < n_rules-1 &&  rules_selHnds[j+1]<0) j++;
      //std::cout << "i,j " << i << " " << j << " " << rules_selHnds[i]<< std::endl;
      molHnd->Select(resSelHnd,STYPE_RESIDUE,rules_selHnds[i],SKEY_NEW);
      molHnd->GetSelIndex(resSelHnd,selRes, nres);
      //std::cout << "nres " << nres << std::endl;
      if ( nres >= 2) {
        fincr = float(j-i)/float(nres-1);
        //std::cout << " nres,fmin,fincr " << nres << " " << fmin << " " << fincr << std::endl;       
      } else {
        fincr = 1.0;
      }
      for (k=0;k<nres-1;k++) {
        selRes[k]->PutUDData(udd,fmin+float(k)*fincr);
        //if ( k == 0 || k == nres-1 ) std::cout << " " << k << " " << fmin+float(k)*fincr << std::endl;
      }
      selRes[nres-1]->PutUDData(udd,fmin+(float(nres-1)*fincr)-0.0001);
      fmin = fmin + float(j-i) + 1.0;
      i = j+1;
    }
    break;

  case BLEND_SECSTR:
  // Assign appropriate colour code to each SSE - save in temporary UDD
    i = 1;
    fmin = 0.0;
    SSE_mapping[SSE_None] = SSE_None;
    SSE_mapping[SSE_Strand] = SSE_Strand;
    SSE_mapping[SSE_Bulge] = SSE_Strand;
    SSE_mapping[SSE_3Turn] = SSE_None;
    SSE_mapping[SSE_4Turn] = SSE_None;
    SSE_mapping[SSE_5Turn] = SSE_None;
    SSE_mapping[SSE_Helix] = SSE_Helix;
    SSE_mapping[SSE_Bridge] = SSE_None;
    SSE_mapping[SSE_Bend] = SSE_None;

    while (i < n_rules-1) {
      j = i+1;
      while ( j+1 < n_rules-1 &&  rules_selHnds[j+1]<0) j++;
      //cout << "i,j " << i << " " << j << endl;
      molHnd->Select(resSelHnd,STYPE_RESIDUE,rules_selHnds[i],SKEY_NEW);
      molHnd->GetSelIndex(resSelHnd,selRes, nres);

      // How many SSEs in this selection?
      nofSSE = 0;
      last_SSE = SSE_mapping[selRes[0]->SSE];
      if (last_SSE != SSE_None) nofSSE++;
      for (k=1;k<nres;k++) {
        if (SSE_mapping[selRes[k]->SSE]!=last_SSE) {
          last_SSE = SSE_mapping[selRes[k]->SSE];
          if (last_SSE != SSE_None) nofSSE++;
        }
      }

      // Calc the fincr
      //cout << "nofSSE " << nofSSE << endl;
      if ( nofSSE >= 2) {
        fincr = float(j-i)/float(nofSSE);
        //cout << " nofSSE,fmin,fincr " << nofSSE << " " << fmin << " " << fincr << endl;       
      } else {
        fincr = 1.0;
      }
      
      nofSSE = 0;
      last_SSE = SSE_None;
      for (k=0;k<nres;k++) {
        if (SSE_mapping[selRes[k]->SSE]!=last_SSE) {
          last_SSE = SSE_mapping[selRes[k]->SSE];
          if (last_SSE != SSE_None) nofSSE++;
        }
        if (SSE_mapping[selRes[k]->SSE]!=SSE_None) {
          selRes[k]->PutUDData(udd,fmin+float(nofSSE-1)*fincr);
        }
      }
      fmin = fmin + float(j-i) + 1.0;    
      i = j+1;
    }
    break;
  }
  
  // Loop over atoms to take colour code from their residue
  for (i=0;i<nat;i++) {
    selAtoms[i]->residue->GetUDData(udd,ff);
    fcolour[i]= ff;
    //cout << " " << i << " " << fcolour[i] ;
  }

  molHnd->DeleteSelection(resSelHnd);

  if (colour_wheel_direction.size() < n_rules ) {
    for (i=colour_wheel_direction.size();i<n_rules;i++) {
      colour_wheel_direction.push_back(COLOUR_WHEEL_CLOCK);
    }
  }
  int rv = atom_colour_vector->SetAtomColours(HSVCOLOURREP,nat,fcolour,GetSchemeCodes(),colour_wheel_direction);
  delete [] fcolour;
  return rv;

}

//------------------------------------------------------------------------
  std::vector <int> CMolColour::GetSchemeCodes() { 
//----------------------------------------------------------------------
    // This method is a kludge to return the colour codes for the BLEND mode
    // (ie from rules_colours) in the same form as 
    // CColourScheme::GetSchemeCodes does for all other colour modes  


    std::vector <int> cols(n_rules);
    //cols[0] = rules_default_colour;
    for (int i=0;i<n_rules;i++) cols[i] = rules_colours[i];
    return cols;
}


//-----------------------------------------------------------------------------
int CMolColour::ByMainSide() {
//-----------------------------------------------------------------------------
  PPCAtom selAtoms;
  int nat,i;
  int *colour;
  
  GetSelection (selAtoms, nat);
  colour = new int[nat];

  for ( i = 0; i < nat ; i++ ) {
    //printf ("Isam %s %i\n",selAtoms[i]->residue->name,selAtoms[i]->residue->isAminoacid());
    if (selAtoms[i]->residue->isAminoacid()) {
      if (molHnd->isMainChain(selAtoms[i])) {
        colour[i] = colourSchemes->MainSide->codes[0];
      }
      else {
        colour[i] = colourSchemes->MainSide->codes[1];
      } 
    }
    else {
      colour[i] = colourSchemes->MainSide->codes[2];
    }
  }
  if(non_C_mode == BYATOMTYPE) ByAtomType ( 0, colour, nat, selAtoms);
  int rv =  atom_colour_vector->SetAtomColours(nat,colour);
  delete [] colour;
  return rv;
}
//-----------------------------------------------------------------------------
int CMolColour::UserScheme() {
//-----------------------------------------------------------------------------
  PPCAtom selAtoms;
  int RC,nat,i,j;
  int *colour;

  RC = GetSelection (selAtoms, nat);
  colour = new int[nat];
  //cout << "CMolColour::UserScheme nat " << nat << endl;
  
  for ( i = 0; i < nat; i++ ) {
    colour[i] = rules_default_colour;
    for (j = n_rules-1 ; j>= 0; j--) {
      if (selAtoms[i]->isInSelection(rules_selHnds[j])) {
        colour[i] = rules_colours[j];
        //cout << "UserScheme " << i << " " << j << endl;
        break;
      }
    }
  }
  int rv= atom_colour_vector->SetAtomColours(nat,colour);
  delete [] colour;
  return rv;
}



//************************************************************************
//  COLOUR SCHEMES
//************************************************************************


//----------------------------------------------------------------------
CColourSchemes::CColourSchemes () {
//----------------------------------------------------------------------
  int RC;

  char *atmtyps[7] = {"*"," C"," O", " N", " S"," H"," P" };
  char *atmcols[7] = { "grey","green", "red", "blue" , "yellow", "grey","magenta" };
  char *restyps[35] = { "*","PHE", "TRP", "TYR", "PRO", "VAL",
		     "ALA", "ILE", "LEU", "SER", "THR",
		     "ASN", "GLN", "ARG", "LYS", "ASP",
		     "GLU", "CYS", "MET", "GLY", "HIS",
		     "A",   "T"  , "G"  , "C"  , "U",
		     "DA",   "DT"  , "DG"  , "DC"  ,
		     "ADE", "THY", "GUA", "CYT", "URA" } ;
  char *rescols[35] = {  "grey","magenta", "magenta", "magenta", "coral", "coral",
		      "coral", "coral", "coral", "cyan", "cyan",
		      "cyan", "cyan",  "blue", "blue", "red",
		      "red", "yellow", "yellow", "white", "light blue",
		      "red", "yellow", "green", "blue", "magenta",
		      "red", "yellow", "green", "blue",
		      "red", "yellow", "green", "blue", "magenta" } ;

  //Secondary structure colouring
  int secstrcods [9] = { SSE_None, SSE_Strand, SSE_Bulge, SSE_3Turn,
			 SSE_4Turn, SSE_5Turn, SSE_Helix, SSE_Bridge, SSE_Bend };
  char *secstrcol[9] = { "grey", "yellow", "green", "pink", "tan", "coral","purple", "green", "light green" }; 

  //Default Bvalue blue->red in range 0->40
  // and red->white in range 40->80
  // and colour below and above that range green and white
  realtype bvalrngs [5] = { 0.0, 0.0, 50.0, 100.0 ,0.0} ;
  char *bvalcol[5] = { "green", "blue","red", "white", "yellow" };

  //Default Res SAS red->blue in range 0->200
  // and colour below and above that range white and yellow
  realtype rsasrngs [5] = { 0.0, 0.0, 50.0, 200.0 ,0.0} ;
  char *rsascol[5] = { "green", "blue","white", "red", "yellow" };

  //Default Atom SAS red->blue in range 0->60
  // and colour below and above that range white and yellow
  realtype asasrngs [5] = { 0.0, 0.0, 10.0, 20.0 ,0.0} ;
  char *asascol[5] = { "green", "blue","white", "red", "yellow" };

  //Default Res Contact red->blue in range 0->100
  // and colour below and above that range white and yellow
  realtype rconrngs [5] = { 0.0, 0.0, 25.0, 100.0 ,0.0} ;
  char *rconcol[5] = { "green", "blue","white", "red", "yellow" };

  //Default Atom Contact red->blue in range 0->50
  // and colour below and above that range white and yellow
  realtype aconrngs [5] = { 0.0, 0.0, 10.0, 50.0 ,0.0} ;
  char *aconcol[5] = { "green", "blue","white", "red", "yellow" };

  //Default Occupancy red->white in range 0->1
  // and colour below and above that range blue and yellow
  realtype occrngs [4] = { 0.0, -0.0001, 1.0001 , 0.0} ;
  char *occcol[4] = { "green", "red","white", "yellow" };

  //Default Charge blue->white in range -1->0
  //white-> red in range 0->1
  // and colour below and above that range cyan and yellow
  realtype chgrngs [5] = { 0.0, -1.0 , 0.0, 1.0 , 0.0} ;
  char *chgcol[5] = { "green","red", "white","blue","yellow" };

  //Default ThroughChain
  realtype thrurngs [7] = { 0.0, 0.0, 0.25, 0.5, 0.75, 1.0 , 0.0} ;
  char *thrucol[7] = { "white", "red", "yellow", "green","blue","purple", "white" };

  // Main-side chain colouring
  char *grptyps[3] = { "AMINO_MAIN", "AMINO_SIDE", "OTHER" };
  //char *grptxts[3] = { "main chain", "side chain", "other" };
  char *grpcols[3] = { "red", "blue", "white" };

  AtomType = new CColourScheme();
  RC = AtomType->SetSchemeString ( 7, atmtyps, atmcols);
  ResType =  new CColourScheme();
  RC = ResType->SetSchemeString ( 31, restyps, rescols );
  BValue = new  CColourScheme();
  RC = BValue->SetSchemeFloat ( 5, bvalrngs, bvalcol );
  Occupancy = new  CColourScheme();
  RC = Occupancy->SetSchemeFloat ( 4, occrngs, occcol );
  Charge = new  CColourScheme();
  RC = Charge->SetSchemeFloat ( 5, chgrngs, chgcol );
  SecStr = new  CColourScheme();
  RC = SecStr->SetSchemeInt ( 9, secstrcods,secstrcol );
  MainSide = new  CColourScheme();
  RC = MainSide->SetSchemeString ( 3, grptyps, grpcols );
  //RC = MainSide->SetSchemeString ( 3, grptyps, grpcols, grptxts );
  ResSAS = new  CColourScheme();
  RC = ResSAS->SetSchemeFloat(5,rsasrngs,rsascol);
  AtomSAS = new  CColourScheme();
  RC = AtomSAS->SetSchemeFloat(5,asasrngs,asascol);
  ResContact = new  CColourScheme();
  RC = ResContact->SetSchemeFloat(5,rconrngs,rconcol);
  AtomContact = new  CColourScheme();
  RC = AtomContact->SetSchemeFloat(5,aconrngs,aconcol);
  ThruChain = new  CColourScheme();
  RC = ThruChain->SetSchemeFloat(7,thrurngs,thrucol);
}


CColourScheme *CColourSchemes::GetScheme ( std::string mode ) {
  if (mode.compare("atomtype")==0) {
    return  AtomType; 
  } else if (mode.compare("restype")==0) {
    return  ResType;
  } else if (mode.compare("secstr")==0) {
    return  SecStr;
  } else if (mode.compare("bvalue")==0) {
    return  BValue;
  } else if (mode.compare("occupancy")==0) {
    return  Occupancy;
  } else if (mode.compare("charge")==0) {
    return  Charge;
  } else if (mode.compare("atom_sas")==0) {
    return  AtomSAS;
  } else if (mode.compare("res_sas")==0) {
    return  ResSAS;
  } else if (mode.compare("atom_contact")==0) {
    return  AtomContact;
  } else if (mode.compare("res_contact")==0) {
    return  ResContact;
  } else if (mode.compare("mainside")==0) {
    return  MainSide;
  } else if (mode.compare("thru_chain")==0) {
    return  ThruChain;
  } else {
    return NULL;
  }
}
std::string CColourSchemes::GetMode ( std::string mode ) {
  if (mode.compare("atomtype")==0) {
    return  "str"; 
  } else if (mode.compare("restype")==0) {
    return  "str";
  } else if (mode.compare("secstr")==0) {
    return  "int";
  } else if (mode.compare("bvalue")==0) {
    return  "float";
  } else if (mode.compare("occupancy")==0) {
    return  "float";
  } else if (mode.compare("charge")==0) {
    return  "float";
  } else if (mode.compare("atom_sas")==0) {
    return  "float";
  } else if (mode.compare("res_sas")==0) {
    return  "float";
  } else if (mode.compare("atom_contact")==0) {
    return  "float";
  } else if (mode.compare("res_contact")==0) {
    return  "float";
  } else if (mode.compare("mainside")==0) {
    return  "str";
  } else if (mode.compare("thru_chain")==0) {
    return  "float";
  } else {
    return "";
  }
}






CColourScheme::CColourScheme() {
  defColour = 0;
  nTypes = 0;
  interpolation_mode = HSVCOLOURREP;
  colour_wheel_direction = COLOUR_WHEEL_CLOCK;
}

CColourScheme::~CColourScheme() {
}

int CColourScheme::SetSchemeInt ( const std::vector<int>& typ , const std::vector<std::string>& cols ) {
  if (typ.size() != cols.size() ) return 1;
  codes = RGBReps::GetColourNumbers (cols);
  nTypes =  int(cols.size());
  itypes = typ;
  colours = cols;
  mode = "int";
  return 0;
}

int CColourScheme::SetSchemeInt ( int n, int typ[] , char *cols[] ) {
  nTypes = n;
  itypes = std::vector<int> (nTypes);
  colours = std::vector<std::string> (nTypes);
  for ( int i = 0; i < n; i++ ) {
    itypes[i] = typ[i];
    colours[i] = cols[i];
  }
  codes = RGBReps::GetColourNumbers (colours);
  mode = "int";
  return 0;
}

int CColourScheme::SetSchemeFloat ( const std::vector<float>& typ , const std::vector<std::string>& cols ) {
  if (typ.size() != cols.size() ) return 1;
  codes = RGBReps::GetColourNumbers (cols);
  //std::cout << "SetSchemeFloat " << codes[0] << " " << codes[1] << " " << codes[2]<< std::endl;
  nTypes =  cols.size();
  ranges = typ;
  colours = cols;
  mode = "float";
  return 0;
}
int CColourScheme::SetSchemeFloat ( int n, realtype typ[] , char *cols[] ) {
  nTypes = n;
  ranges = std::vector<float> (nTypes);
  colours = std::vector<std::string> (nTypes);
  for ( int i = 0; i < n; i++ ) {
    ranges[i] = typ[i];
    colours[i] = cols[i];
  }
  codes = RGBReps::GetColourNumbers (colours);
  //std::cout << "SetSchemeFloat " << codes[0] << " " << codes[1] << " " << codes[2]<< std::endl;

  mode = "float";
  return 0;
}


int CColourScheme::SetSchemeString ( const std::vector<std::string>& typ , const std::vector<std::string>& cols ) {
  if (typ.size() != cols.size() ) return 1; 
  codes = RGBReps::GetColourNumbers (cols);
  nTypes =  cols.size();
  strtypes = typ;
  colours = cols;
  mode = "string";
  return 0;
}

int CColourScheme::SetSchemeString ( int n, char *typ[] , char *cols[] ) {
  nTypes = n;
  strtypes = std::vector<std::string> (nTypes);
  colours = std::vector<std::string> (nTypes);
  for ( int i = 0; i < n; i++ ) {
    strtypes[i] = typ[i];
    colours[i] = cols[i];
  }
  codes = RGBReps::GetColourNumbers (colours);
  mode = "string";
  return 0;
}

 std::vector<int> CColourScheme::GetSchemeInt ( ) {
  return itypes;
}
 std::vector<float> CColourScheme::GetSchemeFloat ( ) {
  return ranges;
}
 std::vector<std::string> CColourScheme::GetSchemeString ( ) {
  return strtypes;
}

 std::vector<std::string> CColourScheme::GetSchemeColours ( ) {
  return colours;
 }

 std::vector<int> CColourScheme::GetSchemeCodes ( ) {
   //std::cout << "GetSchemeCodes " << codes[0] << " " << codes[1] << " " << codes[3] <<std::endl;

   std::vector <int> new_codes(codes.size());
   for (unsigned i=0;i< codes.size();i++) { new_codes[i]=codes[i]; }
   return new_codes;
}

 void CColourScheme::Print() {
   int n = colours.size();
   std::cout << "Print " << n << std::endl;
   for ( int i=0; i< n;i++) { 
     std::cout << ranges[i] << " " << colours[i] << std::endl;
   }
 }

//-----------------------------------------------------------------------
double CColourScheme::GetFColour(double value) {
//-----------------------------------------------------------------------
  int n;
  double fcolour;

  if ( value < ranges[1] ) { 
    fcolour = -1.0;
  } else if ( value > ranges[nTypes-2] ) {
    fcolour = float(nTypes-2);
  } else {
    for ( n = 1; n<nTypes-2;n++) {
      if ( value < ranges[n+1] ) {
        fcolour = float(n-1) + 
	  ((value - ranges[n])/ (ranges[n+1] - ranges[n]));
          break;
      } 
    }
  }
  return fcolour;
}

//-----------------------------------------------------------------------
std::vector<double> CColourScheme::GetRGB(double value) {
//-----------------------------------------------------------------------
  int n;
  double frac;
  std::vector<double> col(4);

  //cout << "GetRGB " << blend_mode << " " << HSVCOLOURREP << endl;
 
  if ( value < ranges[1] ) { 
    return RGBReps::GetColour(codes[0]);
  } else if ( value > ranges[nTypes-2] ) {
    return RGBReps::GetColour(codes[nTypes-1]);
  } else {
    for ( n = 1; n<nTypes-2;n++) {
      if ( value < ranges[n+1] ) {
        frac = (value - ranges[n])/ (ranges[n+1] - ranges[n]);
        if ( blend_mode == HSVCOLOURREP ) {
          std::vector<double> col1 = RGBReps::GetColourHSV(codes[n]);        
          std::vector<double> col2 = RGBReps::GetColourHSV(codes[n+1]);
          for (int i=0;i<3;i++) 
            col1[i] = (frac * col2[i]) + ((1.0-frac)*col1[i]);
          col1[3] = 255.0;
          col = RGBReps::hsvtorgb(col1);
        } else {
          std::vector<double> col1 = RGBReps::GetColour(codes[n]);        
          std::vector<double> col2 = RGBReps::GetColour(codes[n+1]);
          //printf("col1 %f %f %f\n",col1[0],col1[1],col1[2]);
          //printf("col2 %f %f %f\n",col2[0],col2[1],col2[2]);
          for (int i=0;i<3;i++) 
            col[i] = (frac * col2[i]) + ((1.0-frac)*col1[i]);
          col[3] = 255.0;
        }
      break;
      } 
    }
  }
  return col;
}

