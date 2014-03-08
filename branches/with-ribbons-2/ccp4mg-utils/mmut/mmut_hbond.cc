/*
     mmut/mmut_hbond.cc: CCP4MG Molecular Graphics Program
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

Find hydrogen bonds
Liz Potterton Oct02

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
#include "mmut_hbond.h"
#include "mmut_sbase.h"
#include <connect.h>
#include <vector>
#include <mgutil.h>
using namespace std;



//---------------------------------------------------------------------
CHBond::CHBond( PCMMUTManager molHndin, int selHndin,
               PCMMUTManager molHndin2, int selHndin2) :
        CMMANBase ( molHndin, selHndin, molHndin2, selHndin2 ) {
//-------------------------------------------------------------------
  InitParams();
}

//---------------------------------------------------------------------
CHBond::CHBond( PCMMUTManager molHndin) :  CMMANBase ( molHndin) {
//-------------------------------------------------------------------
  InitParams();
}

//-------------------------------------------------------------------
CHBond::CHBond( PCMMUTManager molHndin, int selHndin) : 
                  CMMANBase ( molHndin, selHndin) {
//-------------------------------------------------------------------
  InitParams();
}

//-------------------------------------------------------------------
CHBond::~CHBond( ) {
//-------------------------------------------------------------------
}

//-------------------------------------------------------------------
void CHBond::InitParams() {
//-------------------------------------------------------------------
  //params
  min_D_A = 2.0;
  max_D_A = 3.6;
  max_D_A_S = 3.9;
  max_H_A = 2.5;
  min_DD_D_A = 90 * PI/180.0;
  min_D_A_AA = 90 * PI/180.0;
  min_H_A_AA = 90 * PI/180.0;
  min_D_H_A = 90 * PI/180.0;

}
//-------------------------------------------------------------------
void CHBond::SetParams(int nv, double *value) {
//-------------------------------------------------------------------
  //params

  min_D_A = value[0];
  max_D_A = value[1];
  max_D_A_S = value[2];
  max_H_A = value[3];
  min_DD_D_A = value[4] * PI/180.0;
  min_D_A_AA = value[5] * PI/180.0;
  min_H_A_AA = value[6] * PI/180.0;
  min_D_H_A = value[7] * PI/180.0;

}

//-------------------------------------------------------------------
int CHBond::Calculate(bool separate_models)  {
//-------------------------------------------------------------------
  //Calculate for all NMR models - each needs separate analysis
  if (separate_models ) {
    for (int i = 1; i <= GetMolHnd(0)->GetNumberOfModels(); i++) {
      Calculate0(i);
    }
  }
  else 
    Calculate0();
  return 0;
}

//-------------------------------------------------------------------
int CHBond::Calculate0(int model)  {
//-------------------------------------------------------------------
  PPCAtom donorAtoms = NULL;
  PPCAtom acceptorAtoms = NULL;
  PPCAtom overlapAtoms = NULL;
  int nDonors,nAcceptors;
  int udd,udd2;
  int loop,nLoops = 1;
  int uddD,uddA,selHndD,selHndA;
  PCMMUTManager molHndD,molHndA;
  int selHndDonors,selHndAcceptors,selHndHydrogens,selHndOverlap;
  PSContact contacts = NULL;
  int ncontacts;
  int nOverlap = 0;

  int RC,hbtype1,hbtype2,accept;
  int i,n;
  SAtomBond *donorBond = 0;
  SAtomBond *acceptorBond = 0;
  int nDBonds,nABonds,nH;
  realtype H_A, D_H_A;
  int donor_angles_OK;

  //std::ostringstream label;
  //label.precision(3);

  udd = LoadUDDHBType(molHnds[0]);
  if ( udd <= 0 ) return 2;
  //cout << "udd " << udd;
  
  if (molHnds[1]) {
    udd2 = LoadUDDHBType(molHnds[1]);
    if ( udd2 <= 0 ) return 2;
  } 
  else
    udd2 = udd;

  if ( selHnds[1] >= 0 || molHnds[1] ) 
    nLoops = 2;


  // For two sets of atoms - do the sets overlap?
  if (nLoops == 2 && molHnds[0] == molHnds[1]) {
    selHndOverlap = molHnds[0]->NewSelection();
    molHnds[0]->Select ( selHndOverlap,STYPE_ATOM,selHnds[0],SKEY_NEW );    
    molHnds[0]->Select ( selHndOverlap,STYPE_ATOM,selHnds[1],SKEY_AND );    
    molHnds[0]->GetSelIndex(selHndOverlap,overlapAtoms,nOverlap);
    //cout << "nOverlap " << nOverlap << endl;
  }

  //cout << "nLoops " << nLoops << " udd,udd2 " << udd << " " << udd2 << endl;
  for (loop=0;loop<nLoops;loop++) {
    
    if (nLoops == 1) {
      molHndD = molHnds[0];
      molHndA = molHnds[0];
      selHndA = selHnds[0];
      selHndD = selHnds[0];
      uddD = udd;
      uddA = udd;
    }
    else if ( loop == 0 ) {
      molHndD = molHnds[0];
      selHndD = selHnds[0];
      uddD = udd;
      molHndA = molHnds[1];
      selHndA = selHnds[1];
      uddA = udd2;
    }
    else if ( loop == 1 ) {
      molHndD = molHnds[1];
      selHndD = selHnds[1];
      uddD = udd2;
      molHndA = molHnds[0];
      selHndA = selHnds[0];
      uddA = udd;
    }

    //cout << "molHndD " << molHndD << " selHndD " << selHndD << endl;
    //cout << "molHndA " << molHndA << " selHndA " << selHndA << endl;

  selHndDonors = molHndD->NewSelection();
  molHndD->Select ( selHndDonors,STYPE_ATOM,selHndD,SKEY_NEW );
  if (model > 0 && model <= molHndD->GetNumberOfModels() )
    molHndD->Select (selHndDonors,STYPE_ATOM,model,"*",ANY_RES,
         "*",ANY_RES, "*","*","*","*","*",SKEY_AND);
  molHndD->SelectUDD (selHndDonors,STYPE_ATOM,uddD,
            HBTYPE_DONOR,HBTYPE_BOTH,SKEY_AND);
  molHndD->GetSelIndex(selHndDonors,donorAtoms,nDonors);

  selHndAcceptors = molHndA->NewSelection();
  molHndA->Select ( selHndAcceptors,STYPE_ATOM,selHndA,SKEY_NEW );
  if (model > 0 && model <= molHndA->GetNumberOfModels() )
    molHndA->Select (selHndAcceptors,STYPE_ATOM,model,"*",ANY_RES,
         "*",ANY_RES, "*","*","*","*","*",SKEY_AND);
  molHndA->SelectUDD (selHndAcceptors,STYPE_ATOM,uddA,
           HBTYPE_BOTH,HBTYPE_ACCEPTOR,SKEY_AND);
  molHndA->GetSelIndex(selHndAcceptors,acceptorAtoms,nAcceptors);


  // And are there any actual real hydrogen atoms?
  selHndHydrogens = molHndD->NewSelection();
  molHndD->Select ( selHndHydrogens,STYPE_ATOM,selHndD,SKEY_NEW );    
  molHndD->SelectUDD (selHndHydrogens,STYPE_ATOM,uddD,
            HBTYPE_HYDROGEN,HBTYPE_HYDROGEN,SKEY_AND);

  /*
  cout << "nDonors " << nDonors << " ";
  for (int ii=0; ii<nDonors; ii++) {
    cout << donorAtoms[ii]->name << " "; }
  cout << endl;
    
  cout << "nAcceptors " << nAcceptors << " ";
  for (int ii=0; ii<nAcceptors; ii++) {
    cout << acceptorAtoms[ii]->name << " "; }
  cout << endl;
  */

  //cout << " nDonors " << nDonors << " nAcceptors " << nAcceptors << endl;
  if ( nDonors > 0 && nAcceptors > 0 ) {
    // Find close contacts beteen donors and acceptors
    mat44 * TMatrix=0;

    realtype max_D = max_D_A;
    if(max_D_A_S>max_D)
      max_D = max_D_A_S;
    molHnds[0]->SeekContacts ( donorAtoms,nDonors,acceptorAtoms,nAcceptors,
        min_D_A,max_D,0,contacts,ncontacts, 0,TMatrix, 0 , 0);

    //cout << "ncontacts " << ncontacts << endl;
    if ( !contacts ||  ncontacts <= 0 ) {
      if(selHndAcceptors>=0)molHndA->DeleteSelection(selHndAcceptors);
      if(selHndDonors>=0)molHndD->DeleteSelection(selHndDonors);
      //std::cout << "Returning due to dodgy ncontacts\n";
    }  

    // Loop over the contacts testing the geometry
  
    for ( n = 0; n<ncontacts;n++) {

    //  cout <<  donorAtoms[contacts[n].id1]->name << " " << acceptorAtoms[contacts[n].id2]->name << endl;

    if( molHnds[0]->doAltLocMatch( donorAtoms[contacts[n].id1],
                         acceptorAtoms[contacts[n].id2]) ) {
      donorBond = 0;
      donorAtoms[contacts[n].id1]->GetBonds(donorBond,nDBonds);
      acceptorBond = 0;
      acceptorAtoms[contacts[n].id2]->GetBonds(acceptorBond,nABonds);
      nH = 0;
      // Is donor bonded to any H atoms
      for ( i = 0; i < nDBonds; i++ ) {
        // If there is a hydrogen then test geometry
        if (donorBond[i].atom->isInSelection(selHndHydrogens) ) {
          nH=nH+1;
          H_A = molHnds[0]->BondLength(donorBond[i].atom,
				 acceptorAtoms[contacts[n].id2]);
          D_H_A = molHnds[0]->BondAngle(donorAtoms[contacts[n].id1],
				  donorBond[i].atom,
				  acceptorAtoms[contacts[n].id2]);
          if ( H_A < max_H_A && D_H_A > min_D_H_A ) {
  	    if ( acceptorBond <= 0 ||
            molHnds[0]->BondAngle( donorBond[i].atom,
            acceptorAtoms[contacts[n].id2],acceptorBond[0].atom ) 
             > min_H_A_AA ) {
	      // If the hydrogen atom is part of the original atom selection
              // then record the bond as being to the hydrogen
	        if (donorBond[i].atom->isInSelection(selHndD) )
                  if (loop==1)
                    hbonds.AddConnection(acceptorAtoms[contacts[n].id2],
                        donorBond[i].atom,FloatToString(H_A,"%.1f"));
                  else
                     hbonds.AddConnection(donorBond[i].atom,
                        acceptorAtoms[contacts[n].id2],FloatToString(H_A,"%.1f"));
                else
                   if (loop==1)
                     hbonds.AddConnection(acceptorAtoms[contacts[n].id2],
		       donorAtoms[contacts[n].id1],FloatToString(contacts[n].dist,"%.1f"));
                   else
                     hbonds.AddConnection( donorAtoms[contacts[n].id1],
		       acceptorAtoms[contacts[n].id2],FloatToString(contacts[n].dist,"%.1f"));

            }
	  }
        }
      }	
      if ( nH <= 0 ) {
        // If there were no hydrogens check that bonds between atoms that
        // have HBTYPE_BOTH are only recorded once
        // and check donor-acceptor geometry
	if ( acceptorAtoms[contacts[n].id2]->serNum >  donorAtoms[contacts[n].id1]->serNum )
          accept = 1;
        else {
          RC = acceptorAtoms[contacts[n].id2]->GetUDData(uddA,hbtype1);
          RC = donorAtoms[contacts[n].id1]->GetUDData(uddD,hbtype2);
          accept = hbtype1 -hbtype2 ;
        }
        if ( accept && nOverlap>0 && loop==1 &&
             acceptorAtoms[contacts[n].id2]->isInSelection(selHndOverlap) &&
             donorAtoms[contacts[n].id1]->isInSelection(selHndOverlap) ) {
          //cout << "overlapped " << acceptorAtoms[contacts[n].id2]->name << " " << donorAtoms[contacts[n].id1]->name << endl;
	  accept = 0;
        }

        //cout <<  acceptorAtoms[contacts[n].id2]->serNum << 
	//  acceptorAtoms[contacts[n].id2]->name << " " <<  donorAtoms[contacts[n].id1]->serNum <<  donorAtoms[contacts[n].id1]->name << " " << accept << endl;
        if (   accept != 0 &&
             ( nABonds <= 0 ||
             molHnds[0]->BondAngle (donorAtoms[contacts[n].id1] ,
                        acceptorAtoms[contacts[n].id2],
			acceptorBond[0].atom )  > min_D_A_AA ) ) {
      
          donor_angles_OK = 1;
          for ( i = 0; i < nDBonds; i++ ) {  
            //cout << " " << molHnds[0]->AtomLabel_atom(donorBond[i].atom);
	    //cout << " " << molHnds[0]->AtomLabel_atom(donorAtoms[contacts[n].id1]);
            //cout << " " << molHnds[0]->AtomLabel_atom(acceptorAtoms[contacts[n].id2]);
            //cout  << " " <<  molHnds[0]->BondAngle (donorBond[i].atom,
	    //                 donorAtoms[contacts[n].id1],
	    //                 acceptorAtoms[contacts[n].id2]) << endl;

	    if ( strcmp(donorBond[i].atom->element," H")!= 0 && 
                 molHnds[0]->BondAngle (donorBond[i].atom,
                        donorAtoms[contacts[n].id1] ,
			acceptorAtoms[contacts[n].id2] ) < min_DD_D_A ) {
              donor_angles_OK = 0;
              break;
            }
	  }
  	  if  ( donor_angles_OK > 0) {
            if (loop==1){
              if(
                 (((strcmp(acceptorAtoms[contacts[n].id2]->element," S"))==0||(strcmp(donorAtoms[contacts[n].id1]->element," S")==0))&&contacts[n].dist<=max_D_A_S)
               ||(((strcmp(acceptorAtoms[contacts[n].id2]->element," S"))!=0&&(strcmp(donorAtoms[contacts[n].id1]->element," S")!=0))&&contacts[n].dist<=max_D_A)
               )
              hbonds.AddConnection(acceptorAtoms[contacts[n].id2],donorAtoms[contacts[n].id1],FloatToString(contacts[n].dist,"%.1f"));
            }else{
              if(
                 (((strcmp(acceptorAtoms[contacts[n].id2]->element," S"))==0||(strcmp(donorAtoms[contacts[n].id1]->element," S")==0))&&contacts[n].dist<=max_D_A_S)
               ||(((strcmp(acceptorAtoms[contacts[n].id2]->element," S"))!=0&&(strcmp(donorAtoms[contacts[n].id1]->element," S")!=0))&&contacts[n].dist<=max_D_A)
               )
              hbonds.AddConnection(donorAtoms[contacts[n].id1],acceptorAtoms[contacts[n].id2],FloatToString(contacts[n].dist,"%.1f"));
            }
          }
        }
      }
    }

  } }
  
  //cout << "selHndAcceptors " << selHndAcceptors << selHndDonors << endl;

  if (contacts) {
    delete [] contacts;
    contacts = 0;
  }   
  if(selHndAcceptors>=0)molHndA->DeleteSelection(selHndAcceptors);
  if(selHndDonors>=0)molHndD->DeleteSelection(selHndDonors);
  if(selHndHydrogens>=0)molHndD->DeleteSelection(selHndHydrogens);
  

  }  // End of the big loop


  //cout << "End of Calculate0" << endl;
  return 0;

}

//-----------------------------------------------------------------------
int CHBond::LoadUDDHBType(PCMMUTManager molH ) { 
//-----------------------------------------------------------------------
  int nchains = molH->GetNumberOfChains(1);
  int nr0,nat0;
  int udd;
  //char *foo;

  udd = molH->GetUDDHandle ( UDR_ATOM,"hbtype" );
  if (udd <= 0 ) {
    udd = -1;
    udd = molH->RegisterUDInteger ( UDR_ATOM,"hbtype" );
    if (udd <= 0 ) return udd;
  }
  //cout << "HBond udd " << udd << endl;
  PPCAtom atomTable = NULL;
  for (int im=1;im<=molH->GetNumberOfModels();im++) {
    for (int ic=0;ic<nchains;ic++) {
      nr0 = molH->GetNumberOfResidues(im,ic);
      for (int i=0;i<nr0;i++) {
        molH->GetAtomTable(im,ic,i,atomTable,nat0);
        for (int j=0;j<nat0;j++) {
          atomTable[j]->PutUDData(udd,dynamic_cast<PCMMANManager>(molH)->
            GetAtomHBondType1(atomTable[j]) );
	  //cout << atomTable[j]->GetAtomID(foo) << " " <<
	  //    dynamic_cast<PCMMANManager>(molH)->GetAtomHBondType1(
	  //  atomTable[j]) << endl;
        }
      }
    }
  }
  return udd; 
}

//-----------------------------------------------------------------------
std::string CHBond::Print(bool geometry) {
//-----------------------------------------------------------------------

  SAtomBond *donorBond = 0;
  SAtomBond *acceptorBond = 0;
  int nDBonds,nABonds;

  std::ostringstream output;
  std::string first,second,first_bonded,second_bonded;
  float first_angle, second_angle, distance;
  std::vector<PCAtom>::iterator i = hbonds.pAtom1.begin();
  std::vector<PCAtom>::iterator j = hbonds.pAtom2.begin();
  PCMMUTManager molHnd0,molHnd1;
  output.setf(ios::fixed);
  output.setf(ios::showpoint);
  output.setf(ios::left,ios::adjustfield);
  output.precision(2);

  output << "Neighb set1 (N1)    " 
         << "Hbond atom set1(B1) "
         << "Hbond atom set2 (B2)"
         << "Neighb set2 (N2)    " << endl << "         "
         << "N1-B1-B2 angle      " 
         << "B1-B2 distance      " 
         << "B1-B2-N2 angle" << endl << endl;
  molHnd0 = molHnds[0];
  if (molHnds[1])
    molHnd1 = molHnds[1];
  else
    molHnd1 = molHnds[0];

  while(i!=hbonds.pAtom1.end()){

    first = molHnd0->AtomLabel_atom(*i);
    second =  molHnd1->AtomLabel_atom(*j);

    if ( geometry ) {
      // Ignore the naming of donor/acceptor - the lists indexed
      // by i and j are the atoms from the two separate selections
      donorBond = 0;
      acceptorBond = 0;
      first_bonded = "";
      first_angle = 0.0;
      second_bonded = "";
      second_angle = 0.0;
      (*i)->GetBonds(donorBond,nDBonds);
      if ( nDBonds > 0 ) {
        first_bonded =  molHnd0->AtomLabel_atom(donorBond[0].atom);
        first_angle  = molHnd0->BondAngle (donorBond[0].atom ,
					 *i,*j)*180.0/PI;
      }

      (*j)->GetBonds(acceptorBond,nABonds);
      if ( nABonds > 0 ) {
        second_bonded =  molHnd1->AtomLabel_atom(acceptorBond[0].atom);
        second_angle = molHnd1->BondAngle (*i,*j,
                                 acceptorBond[0].atom )*180.0/PI;
      } 

      distance =  molHnd0->BondLength(*i,*j);

      output << first_bonded << std::setw(20-first_bonded.length()) << " " 
             <<  first << std::setw(20-first.length()) << "  "  
             << second  << std::setw(20-second.length()) <<  " " 
             << second_bonded << endl << "         " 
             << std::setw(20) << first_angle 
             << std::setw(20) << distance 
             << std::setw(20) << second_angle << endl;

    } 
    else {  
      output << second << std::setw(20-second.length()) << "  "  << first  << endl;
    }
    i++;
j++;
  }
  //cout << output.str();
  
  return output.str();
}
