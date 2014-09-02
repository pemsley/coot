/*
     pygl/atom_util.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_ATOM_UTIL_
#define _CCP4MG_ATOM_UTIL_

#include <vector>
#include "cartesian.h"
#include "mgtree.h"
#include <mman_manager.h>

enum enum_AtomColourVector { ACV_COLOUR_MODE_UNSET, ACV_COLOUR_MODE_ICOLOUR, ACV_COLOUR_MODE_PROPERTY,  ACV_COLOUR_MODE_BOTH };

enum enum_ColourWheelDirection { COLOUR_WHEEL_CLOCK, COLOUR_WHEEL_ANTICLOCK };


std::vector <Cartesian> CartesiansFromAtoms(PPCAtom clip_atoms, int nclip_atoms);
int TreeCartesiansToAtoms(Tree *tree, PPCAtom clip_atoms, int nclip_atoms);
PPCAtom GetAtomPair(const std::vector<std::pair<PCAtom,PCAtom> > &pair, int i);
std::vector<std::vector <Cartesian> > GetExternalCartesians(CMMDBManager *molhnd, const std::vector<std::vector<int> > &conn_lists, int side_to_ribbon=0, int side_to_worm=0);

class AtomColourVector {
 private:
  std::vector <int> icolour;
  std::vector <double> property;
  std:: vector <int> res_to_atom;
  std:: vector <int> colour_codes;
  std:: vector <int> colour_wheel_direction;
  std::vector< std:: vector <double> > colours;
  int interpolation_mode;
  int colour_mode;
  int nColours;
  double max_cutoff;
  int last_i;
  
 public:
  AtomColourVector ();
  ~AtomColourVector () {}
  int SetAtomColours(int nat, int *icolour );
  int SetAtomColours(int mode, int nat, double *propertyin,std::vector<int> cols , int direction=COLOUR_WHEEL_CLOCK);
  int SetAtomColours(int mode, int nat, double *propertyin,std::vector<int> cols,std::vector<int> direction  );
  void UpdateColours();
  double* GetRGB(int i);
  double* GetResRGB(int i);
  std::vector<double*> GetRGBVector();
  int SetupResidueColourVector( PCMMDBManager molhnd, int selHnd );
  void UnSetResidueColourVector();
};
  
#endif
