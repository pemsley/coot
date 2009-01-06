/*
     mmut/mmut_colour.h: CCP4MG Molecular Graphics Program
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


#ifndef __CCP4MolColour__ 
#define __CCP4MolColour__

#include <mmdb_manager.h>

#define NCOLOURMODES 6
enum MolColourMode { ONECOLOUR, BYATOMTYPE, BYRESTYPE, BYCHAIN, SECSTR, BVALUE };

DefineClass(CColours);

class CColours {
  friend class CMolColour;
public:
  CColours();
  ~CColours();
    
  static int SetColours(int n, const char *cols[] );
  static int GetCode ( int nAppCol, psvector appCol, ivector code );

  protected:
  static int nColours;
  static psvector names;

};

DefineClass(CColourScheme);




class CColourScheme {
  friend class CMolColour;
public:
  CColourScheme();
  ~CColourScheme();
  int SetScheme ( int n, const char *typs[], const char *cols[]);
  int SetScheme ( int n,  int ityps[], const char *cols[]);
  int SetScheme ( int n, realtype rngs[], int bns[], const char *cols[]);	

 protected:
  int defColour;
  int nTypes;
  psvector types;
  ivector itypes; 
  rvector ranges;
  ivector bins;
  ivector iranges;
  psvector colours;
  ivector codes;
  void FreeMemory();
};

DefineClass(CColourSchemes);

class CColourSchemes {

 public:
  CColourSchemes();
  CColourScheme AtomType;
  CColourScheme ResType;
  CColourScheme BValue;
  CColourScheme SecStr;

};

DefineClass(CMolColour);

//! Define atomic colours for an MMDB selection of atoms
/*! This class is an extention of MMDB functionality to derive per atom
colour code.  The class supports a range of different colouring modes.
*/
 
class CMolColour {
public :
 
  // Constructors
  /*! In context of molecular graphics one instance of CMolColour is
     initialised for each molecule object.  This means that if there
     are multiple molecule objects for one molecule then each can be
     coloured independently.  
  \param molHndin A pointer to the MMDB manager - i.e. the molecule
  \param selHndin A handle for the MMDB selection - i.e. the molecule object
  \param scheme A pointer to the CColourSchemes class which has the definition of the colour schemes.  There is only one instance of this class created by the MG. 
  */
  CMolColour( PCMMDBManager molHndin , int selHndin ,
    PCColourSchemes scheme);
 
 // Destructor
  ~CMolColour();

  int SetSelHandle(int selHndin );
  int GetAtomColourVector ( int &natoms, int *colour ); 
  int SetMode ( int mod );
  int SetOneColour (int icol);
  int Print();



private:

  // the input data
  PCMMDBManager molHnd;
  CColourSchemes *colourSchemes;
  int selHnd;

  // algorithm parameters
  int mode;
  int oneColour;
  int firstChainColour;

  // the derived data
  int natoms;
  ivector colour;

  int Clear();
  int ReColour ();
  int OneColour ();
  int ByAtomType ();
  int ByResidueType ();
  int ByChain();
  int SecondaryStructure ();
  int BValue ();
};
#endif
