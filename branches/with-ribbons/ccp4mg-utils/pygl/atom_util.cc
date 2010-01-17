/*
     pygl/atom_util.cc: CCP4MG Molecular Graphics Program
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
#ifdef _WIN32
#include <windows.h>
#endif

#include <vector>
#include <iostream>
#include <stdlib.h>
#include "mgtree.h"
#include "cartesian.h"
#include <mmdb_manager.h>
#include <catmull.h>
#include "cbuild.h"
#include "atom_util.h"
#include "rgbreps.h"

PPCAtom GetAtomPair(const std::vector<std::pair<PCAtom,PCAtom> > &pair, int i){
  PPCAtom atom_pair = new PCAtom[2];
  atom_pair[0] = pair[i].first;
  atom_pair[1] = pair[i].second;
  return atom_pair;
}

std::vector <Cartesian> CartesiansFromAtoms(PPCAtom clip_atoms, int nclip_atoms){
  std::vector <Cartesian> clip_pos;

  for(int i=0;i<nclip_atoms;i++){
    clip_pos.push_back(Cartesian(clip_atoms[i]->x,clip_atoms[i]->y,clip_atoms[i]->z));
  }
  
  return clip_pos;
}

int TreeCartesiansToAtoms(Tree *tree, PPCAtom clip_atoms, int nclip_atoms){

  std::vector<Cartesian> coords = tree->GetAllCartesians();

  
  for(int i=0;i<nclip_atoms;i++){
    clip_atoms[i]->x = coords[i].get_x();
    clip_atoms[i]->y = coords[i].get_y();
    clip_atoms[i]->z = coords[i].get_z();
  }
  return 0;
}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// AtomColourVector
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

AtomColourVector::AtomColourVector () { 
  colour_mode = ACV_COLOUR_MODE_UNSET; 
}

//--------------------------------------------------------------------------
int AtomColourVector::SetAtomColours (  int nat, int *icolin) { 
//--------------------------------------------------------------------------
  for ( int n=0;n<nat;n++) { icolour.push_back(icolin[n]); }
  
  switch (colour_mode) { 
    case ACV_COLOUR_MODE_PROPERTY:
      colour_mode = ACV_COLOUR_MODE_BOTH;
      break;
    default:
      colour_mode = ACV_COLOUR_MODE_ICOLOUR;
      break;   
  }
  
  interpolation_mode = NOCOLOURREP;
  return 0;
}

//--------------------------------------------------------------------------
int AtomColourVector::SetAtomColours ( int mode, int nat,
        double *propertyin, std::vector<int> cols , int direction ) {
//--------------------------------------------------------------------------
  int n;
  interpolation_mode = mode;
  //std::cout << "SetAtomColours interpolation " << interpolation_mode << std::endl;
  property.clear();
  for ( n=0;n<nat;n++) { property.push_back(propertyin[n]); }
  nColours = cols.size();
  colour_codes.clear();
  for ( n=0;n<nColours;n++) { colour_codes.push_back(cols[n]); }
  colour_wheel_direction.clear();

  // Set colour wheel direction the same for all ranges
  for ( int n=0;n<nColours;n++) { colour_wheel_direction.push_back(direction); }
  max_cutoff = (float)(nColours-3);
  //std::cout << "max_cutoff " << max_cutoff << std::endl;
  
  switch (colour_mode) {
    case ACV_COLOUR_MODE_ICOLOUR:
      colour_mode = ACV_COLOUR_MODE_BOTH;
      break;
    default:
      colour_mode = ACV_COLOUR_MODE_PROPERTY;
      break;   
  }
  
  
  return 0;
}
//--------------------------------------------------------------------------
int AtomColourVector::SetAtomColours ( int mode, int nat,double *propertyin, 
      std::vector<int> cols,  std::vector<int> direction ) {
//--------------------------------------------------------------------------
  int n;
  interpolation_mode = mode;
  //std::cout << "interpolation_mode " << interpolation_mode << " RGB " <<  RGBCOLOURREP << std::endl;
  property.clear();
  for ( n=0;n<nat;n++) { property.push_back(propertyin[n]); }
  //std::cout << "SetAtomColours nat " << property.size() <<  std::endl;
  nColours = cols.size();
  colour_codes.clear();
  for ( n=0;n<nColours;n++) { colour_codes.push_back(cols[n]); }
  colour_wheel_direction.clear();
  for ( int n=0;n<nColours;n++) { colour_wheel_direction.push_back(direction[n]); }
 
  max_cutoff = (float)(nColours-3);
  switch (colour_mode) {
    case ACV_COLOUR_MODE_ICOLOUR:
      colour_mode = ACV_COLOUR_MODE_BOTH;
      break;
    default:
      colour_mode = ACV_COLOUR_MODE_PROPERTY;
      break;   
  }
  
  
  return 0;
}

//-------------------------------------------------------------------------
double *AtomColourVector::GetResRGB ( int i ) {
//-------------------------------------------------------------------------
  //std::cout << "GetResRGB " << i << " " << res_to_atom[i] << std::endl;

  if ( i >= 0 && i < int(res_to_atom.size()) )
    return GetRGB(res_to_atom[i]);
  else {
    double *col = new double[4];
    col[0] = 1.0;
    col[1] = 1.0;
    col[2] = 1.0;
    col[3] = 1.0;
    return col;
  }
}

//----------------------------------------------------------------------
void AtomColourVector::UpdateColours() {
//----------------------------------------------------------------------
  //std::cout << "UpdateColours" << std::endl;
  colours.clear();
  colours = RGBReps::GetColourDefinitions( colour_codes,interpolation_mode);
  //std::cout << "UpdateColours " << colours.size() << std::endl;
}

//-------------------------------------------------------------------------
double *AtomColourVector::GetRGB ( int i ) {
//-------------------------------------------------------------------------
  int j,ir;
  float fColour;
 
  //std::cout << "GetRGB interpolation_mode " << interpolation_mode << std::endl;

  if ( colour_mode == ACV_COLOUR_MODE_ICOLOUR) {
    if(i>=icolour.size())
      return  RGBReps::GetColourP(0);
    else
      return  RGBReps::GetColourP(icolour[i]);

  } else {
  
    if (colours.size()<=0) UpdateColours();

    
    double prop = 0.0;
    if(i<property.size())
       prop = property[i];
    double *col = new double[4];
    col[0] = 0.0;
    col[1] = 0.0;
    col[2] = 0.0;
    col[3] = 1.0;
    if (prop < 0.0&& colours.size()>0) {
      for (j=0;j<3;j++) col[j]= colours[0][j];
    } else if ( interpolation_mode == NOCOLOURREP && (nColours-2)<colours.size()) {
      if ( prop >max_cutoff  ) {
        for (j=0;j<3;j++) {col[j] =  colours[nColours-2][j]; }
      } else if((ir+1)<colours.size()){    
        ir = int(prop);
        //std::cout << prop<< " " << ir << std::endl;
        for (j=0;j<3;j++) {
          col[j] =  colours[ir+1][j];
        }
      }
    } else if ( prop >max_cutoff &&(nColours-1)<colours.size()) {
      for (j=0;j<3;j++) {
        col[j]= colours[nColours-1][j];
      }
    } else if ( interpolation_mode == RGBCOLOURREP &&(ir+2)<colours.size()) {
      ir = int(prop);
      fColour=prop-float(ir);
      for (j=0;j<3;j++) {
        col[j] =  fColour * colours[ir+2][j] +
	    (1-fColour) * colours[ir+1][j];
      }
    } else {
      ir = int(prop);
      fColour=prop-float(ir);
      //std::cout << "GetRGB " << i << " " << prop << " " << ir << " " << fColour << " " << colour_wheel_direction[ir+1] <<  std::endl;
      //std::cout << fColour << " " << ir << " " << colour_wheel_direction[ir+1];
      // dir is the direction round the wheel that we will get if
      // we do the simplest interpolation 
      int dir = COLOUR_WHEEL_CLOCK;
      if (colours[ir+2][0] < colours[ir+1][0]) dir = COLOUR_WHEEL_ANTICLOCK;

      if ( colour_wheel_direction[ir+1] == dir &&(ir+2)<colours.size()) {
          for (j=0;j<3;j++) {
          col[j] =  fColour * colours[ir+2][j] +
	    (1-fColour) * colours[ir+1][j]; }
      } else if((ir+2)<colours.size()){
        if ( colours[ir+1][0] > colours[ir+2][0] ) {
          col[0] =  fColour * colours[ir+2][0] +
            (1-fColour) * ( colours[ir+1][0] - 360.0);
          if (col[0] < 0.0) col[0] = col[0] + 360.0;
        } else {
          col[0] =  fColour * colours[ir+2][0] +
            (1-fColour) * ( colours[ir+1][0] + 360.0);
          if (col[0] > 360.0) col[0] = col[0] - 360.0;
        }
        for (j=1;j<3;j++) {
          col[j] =  fColour * colours[ir+2][j] +
                  (1-fColour) * colours[ir+1][j];
        }
      }
    }
    //std::cout << " " << col[0] << " " << col[1] << " " << col[2];
    //double *col1 = new double[4];
    //col1 = RGBReps::hsvtorgb(col);
    //std::cout << " : " <<  col1[0] << " " << col1[1] << " " << col1[2]<< " " << col1[3] << std::endl;

    if ( interpolation_mode == HSVCOLOURREP ) {
      //return RGBReps::hsvtorgb(col);
      double *col2 = RGBReps::hsvtorgb(col);
      delete col;
      return col2;

    } else {
      return col;
    }
  }
    /*
  default :
    std::cout << "ERROR in GetRGB no colour_mode defined" << std::endl;
    double *col = new double[4];
    col[0] = 1.0;
    col[1] = 1.0;
    col[2] = 1.0;
    col[3] = 1.0;
    break;
    }
    return col;
    */
}


//-----------------------------------------------------------------------
int AtomColourVector::SetupResidueColourVector( PCMMDBManager molHnd, 
                                   int selHnd ) {
//-----------------------------------------------------------------------
  /*
  Invoked from SplineInfo to set up res_to_atom vector which will
  give an index into the atom-based icolour and property vectors
  for any given residue index
  */
 
  PPCAtom pAtoms;
  int CAselHnd,nAtoms,idx;
  PCResidue pRes;

  
  molHnd->GetSelIndex ( selHnd,pAtoms,nAtoms);
  if (nAtoms<=0) {
    std::cout << "No atoms in SetupResidueColourVector" << std::endl;
    return -1;
  }

  // Diagnostic - how mant residues?
  //int res_selHnd = molHnd->NewSelection();
  //molHnd->Select(res_selHnd,STYPE_RESIDUE,selHnd,SKEY_NEW);
  //PPCResidue resTable;
  //int nRes;
  //molHnd->GetSelIndex ( res_selHnd, resTable, nRes );


  CAselHnd = molHnd->NewSelection();
  molHnd->Select(CAselHnd,STYPE_ATOM,0,"*",
		 ANY_RES,"*",ANY_RES,"*","*","CA","C","*",SKEY_NEW);

  if (interpolation_mode == NOCOLOURREP ) {
    if ( nAtoms != int(icolour.size())) {
      std::cout << "Number of atoms mismatch in SetupResidueColourVector, expect "<< icolour.size()<< " but got " << nAtoms << std::endl;
      molHnd->DeleteSelection(CAselHnd);
      return -2;
    }
  } else {
    if ( nAtoms != int(property.size())) {
      std::cout << "Number of atoms mismatch in SetupResidueColourVector, expect "<< property.size()<< " but got " << nAtoms << std::endl;
      molHnd->DeleteSelection(CAselHnd);
      return -2;
    }
  }

  res_to_atom.clear();
  //std::cout << "Initial Length of res_to_atom " << res_to_atom.size() << " nAtoms " << nAtoms << std::endl;
  pRes =  pAtoms[0]->GetResidue();
  idx = 0;
  
  for (int i=1;i<nAtoms;i++) {
    //std::cout << " i " << i <<  std::endl; std::cout.flush();
    // Override existing colour if this atom is a CA
    if ( pAtoms[i]->isInSelection(CAselHnd)) idx = i;
    if (  pAtoms[i]->GetResidue() != pRes ) {
      // Save the atom index of the previous residue
      res_to_atom.push_back(idx);
      // Initialise pRes and idx for this residue
      pRes =  pAtoms[i]->GetResidue();
      idx = i;
    }
  }
  // Save the atom index of the last residue
  res_to_atom.push_back(idx);
  idx = res_to_atom.size();
  //std::cout << "Length of res_to_atom " << idx << " " << nRes << std::endl; std::cout.flush();


  molHnd->DeleteSelection(CAselHnd);
  return nAtoms;
}

void AtomColourVector::UnSetResidueColourVector() {
  res_to_atom.clear();
}

std::vector<double*> AtomColourVector::GetRGBVector() {
  std::vector<double*> vcol;  
  double *rgb;
  
  for(int i=0;i<nColours;i++){
    rgb = GetRGB(i);
    vcol.push_back(rgb);
  }
  return vcol;
}
