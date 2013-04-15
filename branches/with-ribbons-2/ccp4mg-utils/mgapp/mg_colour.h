/*
     mgapp/mg_colour.h: CCP4MG Molecular Graphics Program
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


#ifndef __MgAppColour__ 
#define __MgAppColour__

#include <mman_manager.h>
#include "mgapp_base.h"
#include <mmut_sasarea.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "atom_util.h"

#define NCOLOURMODES 18
enum MolColourMode { ONECOLOUR, BYATOMTYPE, BYRESTYPE, BYCHAIN, SECSTR, 
    SASRES, SASATOM, CONTACTRES, CONTACTATOM, BVALUE, OCCUPANCY, ALTLOC, CHARGE, 
    MAINSIDE, USER_SCHEME, THRU_CHAIN, BYMODEL, BLENDTHRU };

enum MolColourBlend { BLEND_RESIDUE, BLEND_SECSTR, BLEND_CHAIN };


class CColourScheme {
  friend class CMolColour;
 public:
  CColourScheme();
  ~CColourScheme();
  int SetSchemeInt (  const std::vector<int>& ityp , const std::vector<std::string>& cols);
  int SetSchemeInt ( int n, int typ[] , char *cols[] );
  int SetSchemeFloat (  const std::vector<float>& rtyp , const std::vector<std::string>& cols);
  int SetSchemeFloat ( int n, realtype typ[] , char *cols[] );
  int SetSchemeString (  const std::vector<std::string>& chtyp , const std::vector<std::string>& cols);
  int SetSchemeString ( int n, char *typ[] , char *cols[] );
  int SetInterpolationMode(int n, int m=COLOUR_WHEEL_CLOCK){ interpolation_mode = n;colour_wheel_direction=m; return 0; }
  int GetInterpolationMode() { return interpolation_mode ;}
  int GetColourWheelDirection() { return colour_wheel_direction ;}

  std::vector<int> GetSchemeInt ( );
  std::vector<float> GetSchemeFloat ( );
  std::vector<std::string>GetSchemeString ( );
  std::vector<std::string> GetSchemeColours ( );
  std::vector<int> GetSchemeCodes ( );
  void Print(void);
  double GetFColour(double value);
  std::vector<double> GetRGB(double value);
  std::string GetMode() { return mode; }

 protected:
  int defColour;
  int blend_mode;
  int nTypes;
  std::vector<int> itypes;
  std::vector<float> ranges;
  std::vector<std::string> strtypes;
  std::vector<std::string> colours;
  std::vector<int> codes;
  std::string mode;
  int interpolation_mode;
  int colour_wheel_direction;

};


DefineClass(CColourSchemes);

class CColourSchemes {

 public:
  CColourSchemes();
  CColourScheme *AtomType;
  CColourScheme *ResType;
  CColourScheme *BValue;
  CColourScheme *ResSAS;
  CColourScheme *AtomSAS;
  CColourScheme *ResContact;
  CColourScheme *AtomContact;
  CColourScheme *Occupancy;
  CColourScheme *Charge;
  CColourScheme *SecStr;
  CColourScheme *MainSide;
  CColourScheme *ThruChain;

  CColourScheme *GetScheme ( std::string mode );
  std::string GetMode ( std::string mode );

};

DefineClass(CMolColour);

//! Define atomic colours for an MMDB selection of atoms
/*! This class is an extention of MMDB functionality to derive per atom
colour code.  The class supports a range of different colouring modes.
*/
 
class CMolColour: public CMgAppBase {
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
  CMolColour( PCMMANManager molHndin , int selHndin , PCColourSchemes scheme);
 
 // Destructor
  ~CMolColour();

  AtomColourVector GetAtomColourVector ();
  void SetSelHandle ( const int selHndin);
  int SetMode ( int force, int mod=-1 , int non_C_mod=-2, int icol = -1,
            int ncols=0, int icols[] = 0 );
  int SetOneColour (int icol);
  void SetUDD ( int AtomUDDin,int ResUDDin ) { resUDDHnd = ResUDDin; atomUDDHnd =AtomUDDin; }
  int SetScheme ( int nrules, int def_colour, int selHnd[],
                         int rule_colour[], int mode=0);
  void SetColourWheelDirection (std::vector<int> direction  );
   
  
  int GetSchemeN (){ return n_rules; }
  std::vector <int> GetSchemeCodes();
  int EditColour (int icol, const int selHndin);
  int Print();
  int ReColour ();

  private:

  // the input data
  CColourSchemes *colourSchemes;


  // algorithm parameters
  int reapply;
  int if_continuous;
  int mode;
  int non_C_mode;
  int oneColour;
  int nMultColours;
  int *multColours;
  int rules_default_colour;
  int n_rules;
  int rules_mode;
  int resUDDHnd, atomUDDHnd; 
  ivector rules_colours;
  ivector rules_selHnds;
  AtomColourVector atom_colour_vector;
  std::vector<int> colour_wheel_direction;

  int OneColour ();
  int ByAtomType ( );
  int ByAtomType (int colour_carbon, int *colour, int nat, PPCAtom selAtoms);
  int ByResidueType ();
  int ByChain();
  int BySecondaryStructure ();
  int BySASArea(int res_mode,int ifcontact=0);
  int ByBValue ();
  int ByOccupancy ();
  int ByCharge ();
  int ByMainSide();
  int ByThruChain(int chain_no,int ncolours=0,int *colours =0);
  int UserScheme();
  int ByModel();
  int BlendThruChain();

};
#endif
