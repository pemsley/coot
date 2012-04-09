/*
     mmut/mmut_colour.cc: CCP4MG Molecular Graphics Program
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
#include <iostream>
#include <mmdb_manager.h>
#include <mmut_colour.h>

//************************************************************************
// MOLECULE COLOUR CLASS
//************************************************************************

  /*! In context of molecular graphics one instance of CMolColour is
     instantiated for each molecule object.
  \param molHndin A pointer to the MMDB manager - i.e. the molecule
  \param selHndin A handle for the MMDB selection - i.e. the molecule object. See SetSelHandle.
  \param scheme A pointer to the CColourSchemes class which has the definition of the colour schemes.  There is only one instance of this class created by the MG. 
  */
//-----------------------------------------------------------------------
CMolColour::CMolColour(PCMMDBManager molHndin, int selHndin,
		        PCColourSchemes schemein ) {
//-----------------------------------------------------------------------
  //data  - needs initiallising BEFORE SetSelHandle
  natoms = 0;
  colour = NULL;

  //input data
  molHnd = molHndin;
  colourSchemes = schemein;
  SetSelHandle(selHndin);

  //default params
  oneColour = 1;
  firstChainColour = 1;
 
}

//----------------------------------------------------------------------
CMolColour::~CMolColour(){
//----------------------------------------------------------------------
  Clear();
}

  /*! Update the selection handle which specifies a subset of atoms within an MMDBManager object for which the colour will be derived. If the value of selHndin is less than or equal to zero (i.e. an invalid selection handle) then a selection of all atoms within the MMDBManager object will be created and used
  \param selHndin An MMDB selection handle.  
*/

//----------------------------------------------------------------------
int CMolColour::SetSelHandle(int selHndin) {
//----------------------------------------------------------------------
  int RC;
  selHnd = selHndin;
  printf ("SelSelHandle %i\n",selHnd);
  // No selection specified so select all
  if ( selHnd <= 0 ) {
    selHnd = molHnd->NewSelection();
    RC = molHnd->Select ( selHnd,STYPE_ATOM,"*/*/*/*",SKEY_NEW );
  }
 
  Clear();
  return 0;
}

/*! Set the colour to be applied if the colouring mode is ONECOLOUR
    \param icol An integer value which sould be a valid colour code as defined in the CColour class 
    */ 

//----------------------------------------------------------------------
int CMolColour::SetOneColour (int icol) {
//----------------------------------------------------------------------
  oneColour = icol;
  if (mode == ONECOLOUR) Clear();
  return 0;
}

/*! Set the colouring mode to be applied. The valid values for the mode are
    - ONECOLOUR  Apply one colour as given by SetOneColour
    - BYATOMTYPE Colour by the atom element type
    - BYRESTYPE  Colour by residue type
    - BYCHAIN    Colour each chain a different colour
    - SECSTR     Colour by secondary structure
    - BVALUE     Colour by B value
   The detailed interpretation of each of these modes is defined within the CMolColourScheme class
    \param mod Colouring mode   
*/

//----------------------------------------------------------------------
int CMolColour::SetMode ( int mod ) {
//----------------------------------------------------------------------
  if ( mod < 0 || mod >= NCOLOURMODES ) return 1;
  if ( mod != mode ) {
    mode = mod;
    Clear();
  }
  printf("Colour mode %i\n",mode);
  return 0;
}

/*! Return the per atom colour. The values of the vector icol correspond to a colour code defined in CColour. The number of elements in the vector is the number of atoms in the selection. Do NOT delete this vector.
  \param nat Number of atoms in the selection and therefore the number of elements in the vector icol.
  \param icol A pointer to an integer vector.
*/


//----------------------------------------------------------------------
int CMolColour::GetAtomColourVector ( int &nat, int *icol ) {
//----------------------------------------------------------------------
  if ( !colour ) {
    if ( ReColour()) return 1;
  }
  nat = natoms;
  printf("CMolColour::GetAtomColourVector nat %i\n",nat);
  printf("colours %i %i %i\n",colour[0],colour[1],colour[2]);
  icol = colour;
  // printf ("colour %i icol %i\n",colour,icol); colour is and int?
  // No it isn't. comment it.
  return 0;
}

/*! Write a list of the atom id and colour to standard out.  This is primarilly for diagnostic purposes
 */

//-----------------------------------------------------------------------
int CMolColour::Print( ) {
//-----------------------------------------------------------------------
  int i;
  PPCAtom selAtoms;

  if ( !colour ) ReColour();

  molHnd->GetSelIndex (selHnd, selAtoms, natoms );

  for ( i = 0; i < natoms; i++ ) {
    printf ( "colour %i %s %s %i\n",i,selAtoms[i]->residue->name,
	    selAtoms[i]->name, colour[i]);
  }
  return 0;
}

//--------------------------------------------------------------------
// PRIVATE METHODS
//--------------------------------------------------------------------

/* Delete the atom colour vector. 
 */

//----------------------------------------------------------------------
int CMolColour::Clear() {
  //---------------------------------------------------------------------
  if (colour ) {
    //printf ("colour %i\n",colour);
    FreeVectorMemory (colour,0);
    colour = NULL;
    natoms = 0;
  }
  return 0;
}

/* Recalculate the atom colour vector
 */

//----------------------------------------------------------------------
int CMolColour::ReColour() {
//----------------------------------------------------------------------
  int RC;
  switch ( mode ) {
  case ONECOLOUR:
    RC = OneColour();
    break;
  case BYATOMTYPE:
    RC = ByAtomType();
    break;
  case BYRESTYPE:
    RC = ByResidueType();
    break;
  case BYCHAIN:
    RC = ByChain();
    break;
  case SECSTR:
    RC = SecondaryStructure();
    break;
  case BVALUE:
    RC = BValue();
    break;
  default:
    RC = 1;
  }
  printf ("Recolour %i %i\n",mode,RC);
  return RC;
}


//-----------------------------------------------------------------------
int CMolColour::OneColour( ) {
//-----------------------------------------------------------------------
  int i;
  PPCAtom selAtoms;

  molHnd->GetSelIndex (selHnd, selAtoms, natoms );
  FreeVectorMemory( colour,0);
  GetVectorMemory ( colour, natoms, 0);

  for ( i = 0; i < natoms; i++ ) {
    colour[i] = oneColour;
  }

  return 0;
}

//-----------------------------------------------------------------------
int CMolColour::ByAtomType( ) {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int i,j;

  molHnd->GetSelIndex (selHnd, selAtoms, natoms );
  FreeVectorMemory( colour,0);
  GetVectorMemory ( colour, natoms, 0);
  printf("ByAtomType %i\n",natoms);
  for ( i = 0; i < natoms; i++ ) {
    for ( j = 0; j < colourSchemes->AtomType.nTypes; j++ ) {
      if ( strcmp(selAtoms[i]->element,colourSchemes->AtomType.types[j] )==0) {
        colour[i] =colourSchemes->AtomType.codes[j];
        break;
      }
      colour[i] = colourSchemes->AtomType.defColour;
    }
  }
  return 0;
}


//-----------------------------------------------------------------------
int CMolColour::ByResidueType( ) {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int i,j;
 
  molHnd->GetSelIndex (selHnd, selAtoms, natoms );
  FreeVectorMemory( colour,0);
  GetVectorMemory ( colour, natoms, 0);
 
  for ( i = 0; i < natoms; i++ ) {
    for ( j = 0; j < colourSchemes->ResType.nTypes; j++ ) {
      if ( strcmp(selAtoms[i]->residue->name,
		  colourSchemes->ResType.types[j] )==0) {
        colour[i] = colourSchemes->ResType.codes[j];
        break;
      }
      colour[i] = colourSchemes->ResType.defColour;
    }
  }
  return 0;
}

//-----------------------------------------------------------------------
int CMolColour::ByChain( ) {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int natoms;
  int i,j;
  int chnHnd, nChains;
  PPCChain selChains;
  ivector cols;
  int ic = firstChainColour;

  molHnd->GetSelIndex (selHnd, selAtoms, natoms );
  chnHnd = molHnd->NewSelection(); 
  molHnd->Select(chnHnd,STYPE_CHAIN,selHnd,SKEY_NEW);
  molHnd->GetSelIndex (chnHnd, selChains, nChains );
  if ( nChains <= 0 ) return 1;

  GetVectorMemory ( cols, nChains, 0);
  for ( j = 0; j < nChains; j++ ) {
    if ( ic >= CColours::nColours ) ic = 1;
    cols[j] = ic;
    ic = ic++;
  }

  FreeVectorMemory( colour,0);
  GetVectorMemory ( colour, natoms, 0);

  for ( i = 0; i < natoms; i++ ) {
    for ( j = 0; j < nChains; j++ ) {
       if  ( strcmp(selAtoms[i]->GetChain()->GetChainID(),
		   selChains[j]->GetChainID() ) == 0 )  {
	colour[i] = cols[j];
      }
    } 
  }
  FreeVectorMemory(cols,0);
  molHnd->DeleteSelection(chnHnd);
  return 0;
}

//-----------------------------------------------------------------------
int CMolColour::SecondaryStructure( ) {
//-----------------------------------------------------------------------

  return 0;
}

//-----------------------------------------------------------------------
int CMolColour::BValue () {
//-----------------------------------------------------------------------
  PPCAtom selAtoms;
  int i,j;
  float frac;
  ivector colourbin;
  molHnd->GetSelIndex (selHnd, selAtoms, natoms );
  FreeVectorMemory( colour,0);
  GetVectorMemory ( colour, natoms, 0);
  GetVectorMemory ( colourbin, natoms, 0);

  for ( i = 0; i < natoms ; i++ ) {
    for ( j = 1; j < colourSchemes->BValue.nTypes-1; j++) { 
      if ( selAtoms[i]->tempFactor < colourSchemes->BValue.ranges[j] ) {
        if ( colourSchemes->BValue.bins[j-1] < 1 ) 
          colour[i] = colourSchemes->BValue.defColour ;
        else {
          colour[i] =  colourSchemes->BValue.codes[j];
          frac=(selAtoms[i]->tempFactor - colourSchemes->BValue.ranges[j-1])
	    /(colourSchemes->BValue.ranges[j] - colourSchemes->BValue.ranges[j-1]);
          colourbin[i] = floor (frac *
             static_cast<float>(colourSchemes->BValue.bins[j-1]));
        }
        break;
      }
    }
    // atom value is greater than max of the ranges specified
    if  ( colourSchemes->BValue.bins[colourSchemes->BValue.nTypes] >= 1 ) 
      colour[i] = colourSchemes->BValue.codes[colourSchemes->BValue.nTypes];
    else
      colour[i] = colourSchemes->BValue.defColour;
    printf ( "Atom %i  bvalue  %f10.2 , colour %i, colourbin %i\n",
	       i,selAtoms[i]->tempFactor, colour[i], colourbin[i] );
  }
  return 0;
}

//----------------------------------------------------------------------
//                                              COLOURS CLASS
//----------------------------------------------------------------------

/*! This class is instatiated once for MG.  It contains the current colour definitions and methods to access and edit those definitions.  For each colour there is a name and RGB colour definition.
*/

int CColours::nColours = 1;
psvector CColours::names;

//-----------------------------------------------------------------------
CColours::CColours () {
//-----------------------------------------------------------------------
   const char *cols[] = { "white", "blue", "red", "green", "grey", "yellow",
			  "magenta", "royal blue", "cyan", "coral", "pale green",
			  "pink", "lemon", "purple", "tan", "black" };
  GetVectorMemory(names, 1, 0);
  SetColours ( 16, cols );
  
}

//-----------------------------------------------------------------------
CColours::~CColours () {
//-----------------------------------------------------------------------
  FreeVectorMemory(names, 0);
}

//-----------------------------------------------------------------------
int CColours::SetColours (int n, const char *cols[]) {
//-----------------------------------------------------------------------
  int i;
  nColours = n;
  FreeVectorMemory(names,0);
  GetVectorMemory(names,nColours,0);
  for ( i = 0; i < nColours; i++) {
     names[i] = (char *) cols[i];
    //  printf ("SetColours names %i %s\n",i,names[i]);
  }
  return 0;
}

//-----------------------------------------------------------------------
int CColours::GetCode ( int nAppCol, const psvector appCol, ivector code ) {
//-----------------------------------------------------------------------
  //*****BEWARE code is in range 1 to nColours
  int i,j;
  for ( j = 0; j < nAppCol; j++ ) {
    code[j] = 0;
    if ( names[j] ) {                   // check the name is not NULL
      for ( i = 0; i < nColours; i++ ) {
        if ( strcmp(names[i],appCol[j] )== 0 ) {
          code[j] = i+1;
          break;
        }
      }
    }
  }
  return 0;
}

//************************************************************************
//  COLOUR SCHEMES
//************************************************************************


//----------------------------------------------------------------------
CColourSchemes::CColourSchemes () :
   AtomType(),
   ResType(),
   BValue(),
   SecStr() {
//----------------------------------------------------------------------
  int RC;

  const char *atmtyps[6] = {" C"," O", " N", " S"," H"," P" };
  const char *atmcols[6] = { "blue", "red", "green" , "magenta", "yellow", "white" };
  const char *restyps[25] = { "PHE", "TRP", "TYR", "PRO", "VAL",
		     "ALA", "ILE", "LEU", "SER", "THR",
		     "ASN", "GLN", "ARG", "LYS", "ASP",
		     "GLU", "CYS", "MET", "GLY", "HIS",
		     "A",   "T"  , "G"  , "C"  , "U" } ;
  const char *rescols[25] = {  "magenta", "magenta", "coral", "coral", "coral",
		      "coral", "coral", "cyan", "cyan", "cyan",
		      "cyan", "blue",  "blue", "blue", "red",
		      "red", "yellow", "yellow", "white", "royal blue",
		      "red", "yellow", "green", "blue", "magenta" } ;

  //Secondary structure colouring
  int secstrtyps [7] = { 0, 1, 2, 3, 4, 5, 6 };
  const char *secstrcol[7] = { "white", "yellow", "green", "pink",
			       "tan", "coral","purple" }; 

  //Default Bvalue blue->red in 10 bins between 0->100
  // and colour below and above that range white and yellow
  realtype bvalrngs [4] = { NULL, 0.0,50.0 , NULL} ;
  int bvalbns[4] = { 1, 10, 1, NULL };
  const char *bvalcol[4] = { "white", "blue","red", "yellow" };
  
  RC = AtomType.SetScheme ( 3, atmtyps, atmcols);
  RC = ResType.SetScheme ( 25, restyps, rescols );
  RC = BValue.SetScheme ( 4, bvalrngs, bvalbns, bvalcol );
  RC = SecStr.SetScheme ( 7, secstrtyps, secstrcol );
}

//-----------------------------------------------------------------------
CColourScheme::CColourScheme() {
//-----------------------------------------------------------------------
  defColour = 0;
  nTypes = 0;
  colours = NULL;
  types = NULL;
  itypes = NULL;
  codes = NULL;
  bins = NULL;
  ranges = NULL;
  iranges = NULL;
}

//-----------------------------------------------------------------------
CColourScheme::~CColourScheme() {
//-----------------------------------------------------------------------
  FreeMemory();
}

//-----------------------------------------------------------------------
int CColourScheme::SetScheme (int n, const char *typs[], const char *cols[] ) {
//-----------------------------------------------------------------------
  // Set colour scheme which is dependent value of a string attribute
  int RC,i;
  //Free any existing vector memory and reassign
  FreeMemory();
  GetVectorMemory(colours,n,0);
  GetVectorMemory(types,n,0);
  GetVectorMemory(codes,n,0);

  nTypes = n;
  for ( i = 0; i < n; i++ ) {
    types[i] = (char *) typs[i];
    colours[i] = (char *) cols[i];
  } 

  //Convert the test colour string to integer code
  //printf ("SetScheme types %s %s\n",types[0],types[1]);
  RC = CColours::GetCode (nTypes, colours, codes);
  //printf ("codes %i %i \n",codes[0],codes[1]);
  return 0;
}

//-----------------------------------------------------------------------
int CColourScheme::SetScheme (int n, int typs[], const char *cols[] ) {
//-----------------------------------------------------------------------
  // Set colour scheme which is dependent value of a integer attribute
  int RC,i;
  //Free any existing vector memory and reassign
  FreeMemory();
  GetVectorMemory(colours,n,0);
  GetVectorMemory(itypes,n,0);
  GetVectorMemory(codes,n,0);

  nTypes = n;
  for ( i = 0; i < n; i++ ) {
    itypes[i] = typs[i];
    colours[i] = (char *) cols[i];
  } 

  //Convert the test colour string to integer code
  RC = CColours::GetCode (nTypes, colours, codes);
  //printf ("codes %i %i \n",codes[0],codes[1]);
  return 0;
}

//-----------------------------------------------------------------------
int CColourScheme::SetScheme (int n, realtype rngs[], int bns[], const char *cols[]){
//-----------------------------------------------------------------------
  // Set colour scheme which is dependent value of a string attribute
  int RC,i;
  //Free any existing vector memory and reassign
  FreeMemory();
  GetVectorMemory(colours,n,0);
  GetVectorMemory(ranges,n,0);
  GetVectorMemory(bins,n,0);
  GetVectorMemory(codes,n,0);

  nTypes = n;
  for ( i = 0; i < n; i++ ) {
    ranges[i] = rngs[i];
    colours[i] = (char *) cols[i];
    bins[i] = bns[i];
  } 
  colours[i] = (char *) cols[i];
  

  //Convert the test colour string to integer code
  //printf ("SetScheme types %s %s\n",types[0],types[1]);
  RC = CColours::GetCode (nTypes, colours, codes);
  //printf ("codes %i %i \n",codes[0],codes[1]);
  return 0;
}

//-------------------------------------------------------------------------
void CColourScheme::FreeMemory () {
//-------------------------------------------------------------------------
  if (colours) FreeVectorMemory(colours,0);
  colours = NULL;
  if (types) FreeVectorMemory(types,0);
  types = NULL;
  if (itypes) FreeVectorMemory(itypes,0);
  itypes = NULL;
  if (codes) FreeVectorMemory(codes,0);
  ranges = NULL;
  if (ranges) FreeVectorMemory(ranges,0);
  ranges = NULL;
  if (iranges) FreeVectorMemory(iranges,0);
  iranges = NULL;
  if (bins) FreeVectorMemory(bins,0);
  bins = NULL;
}

