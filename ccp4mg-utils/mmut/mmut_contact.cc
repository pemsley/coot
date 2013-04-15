/*
     mmut/mmut_contact.cc: CCP4MG Molecular Graphics Program
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


/*

Find close contacts
Liz Potterton Aug03

*/

#include <iostream>
#include <string.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include "mman_base.h"
#include <mmdb/mmdb_manager.h>
#include "mmut_manager.h"
#include "mman_manager.h"
#include "mmut_contact.h"
#include "mmut_sbase.h"
#include <connect.h>
#include <vector>
#include <mgutil.h>

using namespace std;

//---------------------------------------------------------------------
CContact::CContact( PCMMUTManager molHndin, int selHndin,
               PCMMUTManager molHndin2, int selHndin2) :
        CMMANBase ( molHndin, selHndin, molHndin2, selHndin2 ) {
//-------------------------------------------------------------------
  InitParams();
}

//---------------------------------------------------------------------
CContact::CContact( PCMMUTManager molHndin) :  CMMANBase ( molHndin) {
//-------------------------------------------------------------------
  InitParams();
}

//-------------------------------------------------------------------
CContact::CContact( PCMMUTManager molHndin, int selHndin) : 
                  CMMANBase ( molHndin, selHndin) {
//-------------------------------------------------------------------
  InitParams();
}

//-------------------------------------------------------------------
CContact::~CContact( ) {
//-------------------------------------------------------------------
}

//-------------------------------------------------------------------
void CContact::InitParams() {
//-------------------------------------------------------------------
  //params
  test_VDW_radius = 1;
  label_VDW_radius = 0;
  exclude_hbondable = 0;
  simple_max_cutoff = 4.0;
  simple_min_cutoff = 0.0;
  VDW_fraction_max = 0.9;
  VDW_fraction_min = 0.0;

}
//-------------------------------------------------------------------
void CContact::SetParams(int nv, double *value, int niv, int *ivalue) {
//-------------------------------------------------------------------
  //params

  test_VDW_radius = ivalue[0];
  label_VDW_radius = ivalue[1];
  exclude_hbondable = ivalue[2];
  simple_min_cutoff = value[0];
  simple_max_cutoff = value[1];
  VDW_fraction_min = value[2];
  VDW_fraction_max = value[3];

}

//-------------------------------------------------------------------
int CContact::Calculate(bool separate_models )  {
//-------------------------------------------------------------------
  //Calculate for all NMR models - each needs separate analysis
  int rv = 0;
  //cout << "separate_models " << separate_models << endl;
  if (separate_models && molHnds[0]->GetNumberOfModels() > 1 ) {
    for (int i = 1; i <= molHnds[0]->GetNumberOfModels(); i++) {
      rv = rv + Calculate0(i);
    }
  }
  else 
    rv = rv + Calculate0();
  return rv;
}

//-------------------------------------------------------------------
int CContact::Calculate0(int model)  {
//-------------------------------------------------------------------
  PPCAtom selAtoms;
  PPCAtom selAtoms2;
  int nSelAtoms,nSelAtoms2; 

  realtype min_cutoff,max_cutoff,frac;
  realtype vdw1,vdw2;
  PSContact contacts = NULL;
  int ncontacts;

  int n;
  int bb;

  //cout << " vdw frac min/max " << VDW_fraction_min << " " << VDW_fraction_max << endl;

  if ( model != 0 && molHnds[1]->GetModel(model) == NULL) return 0;

  GetSelection(0,selAtoms,nSelAtoms,model);

  if (molHnds[1]) 
    GetSelection(1,selAtoms2,nSelAtoms2,model);
  else 
    GetSelection(0,selAtoms2,nSelAtoms2,model);
  
  //cout << "nSelAtoms " << model << " " << nSelAtoms << " " << nSelAtoms2 << endl;
  if ( nSelAtoms <= 0 || nSelAtoms2 <= 0 ) return 1;
  // Find the close contacts between two sets of atoms
  mat44 * TMatrix=0;
  if(test_VDW_radius==1) {
    min_cutoff = 5.0 * VDW_fraction_min;
    max_cutoff = 5.0 * VDW_fraction_max;
  }
  else {
    min_cutoff = simple_min_cutoff;
    max_cutoff = simple_max_cutoff;
  }  
  molHnds[0]->SeekContacts ( selAtoms,nSelAtoms,selAtoms2,nSelAtoms2,
        min_cutoff,max_cutoff,0,contacts,ncontacts,0,TMatrix, 0 , 0);

  //cout << "ncontacts " <<min_cutoff << " " << max_cutoff << " " << ncontacts << endl;
  if ( !contacts ||  ncontacts <= 0 ) {
    //std::cout << "Returning due to dodgy ncontacts\n";
    return 2;
  }  

  // Loop over the contacts testing VDW distance
  if(test_VDW_radius==1 || label_VDW_radius) {
    for (n=0;n<ncontacts;n++) {
      vdw1=dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius(selAtoms[contacts[n].id1]);
      if (molHnds[1])
        vdw2=dynamic_cast<PCMMANManager>(molHnds[1])->GetAtomVDWRadius(selAtoms2[contacts[n].id2]);
      else
        vdw2=dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius(selAtoms2[contacts[n].id2]);
      frac = contacts[n].dist/(vdw1+vdw2);
      //cout << contacts[n].dist << " " << vdw1 << " " << vdw2 << " " << frac << endl;
      if ( !test_VDW_radius ||
         ( frac < VDW_fraction_max && frac > VDW_fraction_min &&
         molHnds[0]->doAltLocMatch(selAtoms[contacts[n].id1],selAtoms2[contacts[n].id2]))) {
         bb = dynamic_cast<PCMMANManager>(molHnds[0])->TestBonding( 
          selAtoms[contacts[n].id1],selAtoms2[contacts[n].id2]);
        if (dynamic_cast<PCMMANManager>(molHnds[0])->TestBonding( 
	  selAtoms[contacts[n].id1],selAtoms2[contacts[n].id2]) == 0) {
	  if (label_VDW_radius) {
            close_contacts.AddUniqueConnection(selAtoms[contacts[n].id1],
              selAtoms2[contacts[n].id2],FloatToString(frac,"%.2f"));
          } else {
            close_contacts.AddUniqueConnection(selAtoms[contacts[n].id1],
              selAtoms2[contacts[n].id2],FloatToString(contacts[n].dist,"%.1f"));
          }
	}
      }
    }
  }
  else {
    //int exclude_hbondable = 1;
    for (n=0;n<ncontacts;n++) {
      bb = dynamic_cast<PCMMANManager>(molHnds[0]) ->TestBonding( 
	     selAtoms[contacts[n].id1],selAtoms2[contacts[n].id2]);

      if ( (bb == 0 || (exclude_hbondable==0 && bb==5)) &&
         molHnds[0]->doAltLocMatch(selAtoms[contacts[n].id1],selAtoms2[contacts[n].id2]) )
        close_contacts.AddUniqueConnection(selAtoms[contacts[n].id1],
        selAtoms2[contacts[n].id2],FloatToString(contacts[n].dist,"%.1f"));
    }
  }
    	
 
  return 0;

}

//-----------------------------------------------------------------------
std::string CContact::Print(bool geometry) {
//-----------------------------------------------------------------------

  std::ostringstream output;
  std::vector<PCAtom>::iterator i = close_contacts.pAtom1.begin();
  std::vector<PCAtom>::iterator j = close_contacts.pAtom2.begin();
  std::string first,second;
  realtype dist,vdw1,vdw2,sum_vdw,frac;
  output.setf(ios::fixed);
  output.setf(ios::showpoint);
  output.setf(ios::left,ios::adjustfield);
  output.precision(2);

  output << "Atom in set 2        " << 
            "Atom in set 1        " << 
            "Distance        " << 
            "Sum of VDWs     " << 
            "Fraction of VDWs" <<  
         endl;
 
  while(i!=close_contacts.pAtom1.end()){

    first = molHnds[0]->AtomLabel_atom(*i);
    second =  molHnds[0]->AtomLabel_atom(*j);
    dist = molHnds[0]->BondLength(*i,*j);
    vdw1=dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius(*i);
      if (molHnds[1])
        vdw2=dynamic_cast<PCMMANManager>(molHnds[1])->GetAtomVDWRadius(*j);
      else
        vdw2=dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius(*j);
      sum_vdw = vdw1+vdw2;
      frac = dist/sum_vdw;
      output << first << std::setw(21-first.length()) << "  "; 
      output << second << std::setw(21-second.length()) << " "; 
      output << std::setw(15) << dist << " " ;
      output << std::setw(15) << sum_vdw << " ";
      output << std::setw(15) << frac << endl;
      i++;j++;
  }
  
  return output.str();
}
