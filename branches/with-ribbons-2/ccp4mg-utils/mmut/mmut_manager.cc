/*
     mmut/mmut_manager.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2011 University of York
     Copyright (C) 2012 STFC

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
#include <vector>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <sstream>
#include "mmut_manager.h"
#include <mmdb/mmdb_utils.h>
#include <mmdb/mmdb_tables.h>
#include <mmdb/mmdb_cifdefs.h>

#include "mmut_util.h"

using namespace std;

//  =====================   CMMUTManager   =======================

 
  CMMUTManager::CMMUTManager() : CMMDBManager()  {

    iatom_types=NULL;
    iatom_type_lookup=NULL;

  }


  CMMUTManager::~CMMUTManager()  { 
   
  }


int const nhydro[nResNames] = { // No of hydrogen atoms in the 27 residues.
5,13,6,4,5,5,8,6,3,8,11,11,13,9,9,7,5,7,10,9,9,17,2,0,0,2
};



void CMMUTManager::ListAtomInfo (int selHnd) {

  PCAtom atom;
  int i;

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  realtype atomweight, molweight=0.0;


  for(i=0;i<natoms_sel;i++){
    atom = selected_atoms[i];
    cout << "Atom(" << atom->GetModelNum() << "," << atom->GetChainID() << "," << atom->GetResName() <<  "," <<  atom->GetSeqNum() << "): ";
    cout << " (Name " << atom->name;
    cout << ") (Symbol " << atom->element;
    cout << ") (Charge " << atom->charge;
    cout << ") (Atomic Number " << getElementNo(atom->element)+1;
    atomweight = MolecWeight[getElementNo(atom->element)];
    cout << ") (Atomic Weight " << atomweight;
    cout << ") (Co-ords " << atom->x << ", " << atom->y << ", " << atom->z;
    cout << ")" << endl;
    molweight += atomweight;
  }

  cout << "Total Molecular Weight " << molweight << endl;

}

realtype CMMUTManager::MolWeight(int selHnd){

  int l;

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  realtype molweight=0.0;

  for(l=0;l<natoms_sel;l++){
    molweight += MolecWeight[getElementNo((selected_atoms[l])->element)];
  }
  return molweight;
}

realtype CMMUTManager::MolWeightWithH(int selHnd){

  return MolWeight(selHnd)+
    NumberOfHydrogens(selHnd)*MolecWeight[0];

}

int CMMUTManager::ResNoLookup(pstr resname){

  int k;

  for(k=0;k<nResNames;k++){
    if(!strcmp(resname,ResidueName[k])) return k;
  }
  
  return 0;
}

int CMMUTManager::NumberOfHydrogens(int selHnd){

  int nh=0;
  int i;

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  nh = nhydro[ResNoLookup(selected_atoms[0]->GetResName())];

  for(i=1;i<natoms_sel;i++){
    if(selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum()){
      nh += nhydro[ResNoLookup(selected_atoms[i]->GetResName())];
    }
  }

  return nh;

}

//----------------------------------------------------------------------------
int CMMUTManager::NumberOfAtoms(int selHnd) {
//----------------------------------------------------------------------------

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  return natoms_sel;
}

//----------------------------------------------------------------------------
realtype CMMUTManager::Mass(int selHnd) {
//----------------------------------------------------------------------------

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int l;
  realtype weight = 0.0;
  for(l=0;l<natoms_sel;l++){
    weight += MolecWeight[getElementNo((selected_atoms[l])->element)];
  }

  return weight;
}


Cartesian CMMUTManager::CentreOfCoordinatesAsCartesian(int selHnd) {
  realtype* com = CentreOfCoordinates(selHnd);
  Cartesian cart;
  cart.set_x(com[0]);
  cart.set_y(com[1]);
  cart.set_z(com[2]);
  delete [] com;
  return cart;
}

Cartesian CMMUTManager::CentreOfMassAsCartesian(int selHnd) {
  realtype* com = CentreOfMass(selHnd);
  Cartesian cart;
  cart.set_x(com[0]);
  cart.set_y(com[1]);
  cart.set_z(com[2]);
  delete [] com;
  return cart;
}

//----------------------------------------------------------------------------
realtype* CMMUTManager::CentreOfMass(int selHnd) {
//----------------------------------------------------------------------------

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int l;
  realtype *com;
  realtype atweight;
  realtype weight;

  weight = 0.0;
  com = new realtype[3];

  com[0] = 0.0;
  com[1] = 0.0;
  com[2] = 0.0;

  for(l=0;l<natoms_sel;l++){
    atweight = MolecWeight[getElementNo((selected_atoms[l])->element)];
    com[0] +=  atweight * (selected_atoms[l])->x;
    com[1] +=  atweight * (selected_atoms[l])->y;
    com[2] +=  atweight * (selected_atoms[l])->z;
    weight += atweight;
  }

  com[0] = com[0]/weight;
  com[1] = com[1]/weight;
  com[2] = com[2]/weight;

  return com;

}

//----------------------------------------------------------------------------
realtype* CMMUTManager::CentreOfCoordinates(int selHnd) {
//----------------------------------------------------------------------------

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int l;
  realtype *com;
  realtype atweight;
  realtype weight;

  weight = 0.0;
  com = new realtype[3];

  com[0] = 0.0;
  com[1] = 0.0;
  com[2] = 0.0;

  for(l=0;l<natoms_sel;l++){
    com[0] +=  (selected_atoms[l])->x;
    com[1] +=  (selected_atoms[l])->y;
    com[2] +=  (selected_atoms[l])->z;
    weight += 1.0;
  }

  com[0] = com[0]/weight;
  com[1] = com[1]/weight;
  com[2] = com[2]/weight;

  return com;

}


//-------------------------------------------------------------------------
realtype CMMUTManager::ExtentSize(int selHnd) {
//-------------------------------------------------------------------------
  realtype* rtde = Extent(selHnd);

  realtype mine[3] = {rtde[0],rtde[1],rtde[2]};
  realtype maxe[3] = {rtde[3],rtde[4],rtde[5]};
  realtype theSize = fabs(mine[0]-maxe[0]);
  if(fabs(mine[1]-maxe[1])>theSize)
     theSize = fabs(mine[1]-maxe[1]);
  if(fabs(mine[2]-maxe[2])>theSize)
     theSize = fabs(mine[2]-maxe[2]);
  return theSize;
}

//-------------------------------------------------------------------------
realtype* CMMUTManager::Extent(int selHnd) {
//-------------------------------------------------------------------------

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int l;
  realtype *com;
  com = new realtype[6];
  realtype xmin = 9999999.9;
  realtype ymin = 9999999.9;
  realtype zmin = 9999999.9;
  realtype xmax = -9999999.9;
  realtype ymax = -9999999.9;
  realtype zmax = -9999999.9;


  for(l=0;l<natoms_sel;l++){
    if (selected_atoms[l]->x < xmin ) xmin = selected_atoms[l]->x;
    if (selected_atoms[l]->y < ymin ) ymin = (selected_atoms[l])->y;
    if (selected_atoms[l]->z < zmin ) zmin = (selected_atoms[l])->z;
    if (selected_atoms[l]->x > xmax ) xmax = (selected_atoms[l])->x;
    if (selected_atoms[l]->y > ymax ) ymax = (selected_atoms[l])->y;
    if (selected_atoms[l]->z > zmax ) zmax = (selected_atoms[l])->z;
  }
  
  com[0] = xmin;
  com[1] = ymin;
  com[2] = zmin;
  com[3] = xmax;
  com[4] = ymax;
  com[5] = zmax;

  return com;

}

std::vector<Cartesian> CMMUTManager::GetPrincipalComponents(int selHnd){
  std::vector<Cartesian> pca;
  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);
  std::vector<Cartesian> carts = PPCAtomsToCartesians(natoms_sel,selected_atoms);
  pca = Cartesian::PrincipalComponentAnalysis(carts);
  return pca;
}

int* CMMUTManager::AtomicComposition(int selHnd){
	
  int l;

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int *atomic_comp;
  atomic_comp = new int[nElementNames];

  for(l=0;l<nElementNames;l++) atomic_comp[l] = 0;

  for(l=0;l<natoms_sel;l++){
    atomic_comp[getElementNo(selected_atoms[l]->element)]++;
  }

  return atomic_comp;

}

void CMMUTManager::PrintAtomicComposition(int selHnd){

  int m;
  int *atomic_comp;

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  atomic_comp = AtomicComposition(selHnd);

  cout << "\n          Atomic Composition\n";
  cout << "          ------------------\n\n";
  for(m=0;m<nElementNames;m++){
    if(atomic_comp[m]>0) cout << ElementName[m] << ": " << atomic_comp[m] << endl;
  }
}

int* CMMUTManager::ResidueComposition(int selHnd) {
	
  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i,l;
  pstr resname;
  int *residue_comp;

  residue_comp = new int[nResNames];
  for(l=0;l<nResNames;l++)
    residue_comp[l]=0;

  resname = selected_atoms[0]->GetResName();
  for(l=0;l<nResNames;l++){
    if(!strcmp(resname,ResidueName[l])){
      residue_comp[l]++;
      //	  sequence[i] = ResidueName1[l];
    }
  }
  if(!strcmp(resname,"HOH")){
    residue_comp[22]++;        // HARDWIRE!!
    //	sequence[i] = 'O'; // Hackery since we have WAT in RNames
  }

  for(i=1;i<natoms_sel;i++){
    if(selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum()){
      resname = selected_atoms[i]->GetResName();
      for(l=0;l<nResNames;l++){
	if(!strcmp(resname,ResidueName[l])){
	  residue_comp[l]++;
	  //	  sequence[i] = ResidueName1[l];
	}
      }
      if(!strcmp(resname,"HOH")){
	residue_comp[22]++;        // HARDWIRE!!
	//	sequence[i] = 'O'; // Hackery since we have WAT in RNames
      }
    }
  }
  //  sequence[nrestot] = '\0';

  return residue_comp;

}

/*
pstr CMMUTManager::GetSequence0(int selHnd){
	
  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i,l;
  int numres;
  pstr resname,sequence;

  int maxresnum = selected_atoms[natoms_sel-1]->GetResidueNo();
  std::cout << maxresnum << "\n";

  numres = TotalNumRes(selHnd);
  //std::cout << "total residues : " << numres << "\n";
  sequence = new char[maxresnum+1];

  int ires = 0;
  for(i=0;i<natoms_sel;i++){
    if(selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum()){
      //std::cout << i << "\n";
      resname = selected_atoms[i]->GetResName();
      //std::cout << resname << "\n";
      bool haveThis = false;
      for(l=0;l<nResNames;l++){
	if(!strcmp(resname,ResidueName[l])){
	  //	  residue_comp[l]++;
          //std::cout << ResidueName1[l] << "(" << ires << ")\n";
	  if(ires<numres) sequence[ires] = ResidueName1[l];
          haveThis = true;
	}
      }
      if(!strcmp(resname,"HOH")){
	//	residue_comp[22]++;        // HARDWIRE!!
        //std::cout << 'O' << "(" << ires << ")\n";
	if(ires<numres) sequence[ires] = 'O'; // Hackery since we have WAT in RNames
        haveThis = true;
      }
      if(!haveThis&&ires<numres) sequence[ires] = '?'; 
      ires++;
    }
  }

  //std::cout << sequence[numres-1] << "\n";
  sequence[ires] = '\0';

  return sequence;

}
*/

void CMMUTManager::SelectOneAtomForNonAminoNucleotide(int selHnd, int atomSelHnd){
/* Given a selection of residues selHnd, a single atom is added to atom selection atomSelHnd,
 * if the residue is neither amino acid nor nucleotide. */
	
  PPCResidue selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);
  char AtomID[30];

  int i,l;
  int numres;
  pstr resname;

  numres = natoms_sel;
  if(numres==0||natoms_sel==0){
     if(getenv("CCP4MG_DEBUG_SEQUENCE")) std::cout << "No residues, returning\n";
     return;
  } 

  resname = selected_atoms[0]->GetResName();
  bool haveThis = false;
  for(l=0;l<nResNames;l++){
    if(!strcmp(resname,ResidueName[l])){
      haveThis = true;
    }
  }
  for(l=0;l<nNucleotideNames;l++){
    if(!strcmp(resname,NucleotideName[l])){
      if(strlen(NucleotideName[l])==1){
        haveThis = true;
      }else if(strlen(NucleotideName[l])==2){
        haveThis = true;
      }
    }
  }
  if(!strcmp(resname,"HOH")){
    haveThis = true;
  }
  if(!haveThis){
   if(!selected_atoms[0]->isNucleotide()&&!selected_atoms[0]->isAminoacid()&&selected_atoms[0]->GetNumberOfAtoms()>0){
     if(getenv("CCP4MG_DEBUG_SEQUENCE")) {
       std::cout <<  resname << " is not nucleotide nor amino acid\n";
     }
   }
   if(getenv("CCP4MG_DEBUG_SEQUENCE")) std::cout << selected_atoms[0]->GetAtom(0)->GetAtomID(AtomID) << "\n";
   if(!selected_atoms[0]->isInSelection(selHnd))
     Select(atomSelHnd,STYPE_ATOM,selected_atoms[0]->GetAtom(0)->GetAtomID(AtomID),SKEY_OR);
  }

  int ires = 1;
  for(i=1;i<natoms_sel;i++){
    if(selected_atoms[i]->seqNum!=selected_atoms[i-1]->seqNum
      ||(
      (selected_atoms[i]->seqNum==selected_atoms[i-1]->seqNum)&&((strlen(selected_atoms[i]->insCode)!=strlen(selected_atoms[i-1]->insCode))||(strncmp(selected_atoms[i]->insCode,selected_atoms[i-1]->insCode,strlen(selected_atoms[i]->insCode))))
      )
      ){
      resname = selected_atoms[i]->GetResName();
      bool haveThis = false;
      for(l=0;l<nResNames;l++){
	if(!strcmp(resname,ResidueName[l])){
          haveThis = true;
	}
      }
      for(l=0;l<nNucleotideNames;l++){
	if(!strcmp(resname,NucleotideName[l])){
	  if(ires<numres){
            if(strlen(NucleotideName[l])==1){
              haveThis = true;
            }else if(strlen(NucleotideName[l])==2){
              haveThis = true;
            }
          }
	}
      }
      if(!strcmp(resname,"HOH")){
        haveThis = true;
      }
      if(!haveThis&&ires<numres){
        if(!selected_atoms[ires]->isNucleotide()&&!selected_atoms[ires]->isAminoacid()&&selected_atoms[ires]->GetNumberOfAtoms()>0){
          if(getenv("CCP4MG_DEBUG_SEQUENCE")) {
            std::cout <<  resname << " is not nucleotide nor amino acid\n";
          }
       }
       if(getenv("CCP4MG_DEBUG_SEQUENCE")) std::cout << selected_atoms[ires]->GetAtom(0)->GetAtomID(AtomID) << "\n";
       if(!selected_atoms[ires]->isInSelection(selHnd))
         Select(atomSelHnd,STYPE_ATOM,selected_atoms[ires]->GetAtom(0)->GetAtomID(AtomID),SKEY_OR);
      }
      ires++;
    }
  }

}

pstr CMMUTManager::GetSequenceFromResidues(int selHnd){
	
  PPCResidue selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i,l;
  int numres;
  pstr resname,sequence;

  numres = natoms_sel;
  //std::cout << "total residues : " << natoms_sel << "\n";
  sequence = new char[numres+1];
  if(numres==0||natoms_sel==0){
    sequence[0] = '\0';
    return sequence;
  } 

  resname = selected_atoms[0]->GetResName();
  bool haveThis = false;
  for(l=0;l<nResNames;l++){
    if(!strcmp(resname,ResidueName[l])){
      //std::cout << "First one: " << ResidueName1[l] << "\n";
      //      residue_comp[l]++;
      sequence[0] = ResidueName1[l];
      haveThis = true;
    }
  }
  for(l=0;l<nNucleotideNames;l++){
    if(!strcmp(resname,NucleotideName[l])){
      //      residue_comp[l]++;
      if(strlen(NucleotideName[l])==1){
        sequence[0] = NucleotideName[l][0];
        haveThis = true;
      }else if(strlen(NucleotideName[l])==2){
        sequence[0] = NucleotideName[l][1];
        haveThis = true;
      }
    }
  }
  if(!strcmp(resname,"HOH")){
    //    residue_comp[22]++;        // HARDWIRE!!
    sequence[0] = 'O'; // Hackery since we have WAT in RNames
    haveThis = true;
  }
  if(!haveThis){
   sequence[0] = 'X'; 
  }
  sequence[1] = '\0';
  //std::cout << "Current sequence:" << sequence << "\n";

  int ires = 1;
  for(i=1;i<natoms_sel;i++){
    if(selected_atoms[i]->seqNum!=selected_atoms[i-1]->seqNum
      ||(
      (selected_atoms[i]->seqNum==selected_atoms[i-1]->seqNum)&&((strlen(selected_atoms[i]->insCode)!=strlen(selected_atoms[i-1]->insCode))||(strncmp(selected_atoms[i]->insCode,selected_atoms[i-1]->insCode,strlen(selected_atoms[i]->insCode))))
      )
      ){
      //std::cout << i << "\n";
      resname = selected_atoms[i]->GetResName();
      //std::cout << resname << "\n";
      bool haveThis = false;
      for(l=0;l<nResNames;l++){
	if(!strcmp(resname,ResidueName[l])/*&&selected_atoms[i]->isAminoacid()*/){
	  //	  residue_comp[l]++;
          //std::cout << ResidueName1[l] << "(" << ires << ")\n";
	  if(ires<numres){
            sequence[ires] = ResidueName1[l];
            sequence[ires+1] = '\0';
          }
          //std::cout << "Current sequence:" << sequence << "\n";
          haveThis = true;
	}
      }
      for(l=0;l<nNucleotideNames;l++){
	if(!strcmp(resname,NucleotideName[l])){
	  //	  residue_comp[l]++;
          //std::cout << ResidueName1[l] << "(" << ires << ")\n";
	  if(ires<numres){
            if(strlen(NucleotideName[l])==1/*&&selected_atoms[i]->isNucleotide()*/){
              sequence[ires] = NucleotideName[l][0];
              haveThis = true;
            }else if(strlen(NucleotideName[l])==2){
              sequence[ires] = NucleotideName[l][1];
              haveThis = true;
            }
          }
	}
      }
      if(!strcmp(resname,"HOH")){
	//	residue_comp[22]++;        // HARDWIRE!!
        //std::cout << 'O' << "(" << ires << ")\n";
	if(ires<numres) sequence[ires] = 'O'; // Hackery since we have WAT in RNames
        haveThis = true;
      }
      if(!haveThis&&ires<numres){
         sequence[ires] = 'X'; 
      }
      ires++;
    }
  }

  //std::cout << sequence[numres-1] << "\n";
  sequence[ires] = '\0';

  if(getenv("CCP4MG_DEBUG_SEQUENCE")) std::cout <<  "Sequence length: " << strlen(sequence) << ", sequence--" << sequence << "\n";
  return sequence;

}

pstr CMMUTManager::GetSequence(int selHnd){
	
  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i,l;
  int numres;
  pstr resname,sequence;

  numres = TotalNumRes(selHnd);
  std::cout << "total residues : " << numres << "\n";
  sequence = new char[numres+1];
  if(numres==0||natoms_sel==0){
    sequence[0] = '\0';
    return sequence;
  } 

  resname = selected_atoms[0]->GetResName();
  bool haveThis = false;
  for(l=0;l<nResNames;l++){
    if(!strcmp(resname,ResidueName[l])){
      //      residue_comp[l]++;
      sequence[0] = ResidueName1[l];
      haveThis = true;
    }
  }
  for(l=0;l<nNucleotideNames;l++){
    if(!strcmp(resname,NucleotideName[l])){
      //      residue_comp[l]++;
      if(strlen(NucleotideName[l])==1){
        sequence[0] = NucleotideName[l][0];
        haveThis = true;
      }else if(strlen(NucleotideName[l])==2){
        sequence[0] = NucleotideName[l][1];
        haveThis = true;
      }
    }
  }
  if(!strcmp(resname,"HOH")){
    //    residue_comp[22]++;        // HARDWIRE!!
    sequence[0] = 'O'; // Hackery since we have WAT in RNames
    haveThis = true;
  }
  if(!haveThis) sequence[0] = 'X'; 

  int ires = 1;
  for(i=1;i<natoms_sel;i++){
    if((selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum())
      ||(
      (selected_atoms[i]->GetSeqNum()==selected_atoms[i-1]->GetSeqNum())&&((strlen(selected_atoms[i]->residue->insCode)!=strlen(selected_atoms[i-1]->residue->insCode))||(strncmp(selected_atoms[i]->residue->insCode,selected_atoms[i-1]->residue->insCode,strlen(selected_atoms[i]->residue->insCode))))
      )
      ){
      //std::cout << i << "\n";
      resname = selected_atoms[i]->GetResName();
      //std::cout << resname << "\n";
      bool haveThis = false;
      for(l=0;l<nResNames;l++){
	if(!strcmp(resname,ResidueName[l])){
	  //	  residue_comp[l]++;
          //std::cout << ResidueName1[l] << "(" << ires << ")\n";
	  if(ires<numres) sequence[ires] = ResidueName1[l];
          haveThis = true;
	}
      }
      for(l=0;l<nNucleotideNames;l++){
	if(!strcmp(resname,NucleotideName[l])){
	  //	  residue_comp[l]++;
          //std::cout << ResidueName1[l] << "(" << ires << ")\n";
	  if(ires<numres){
            if(strlen(NucleotideName[l])==1){
              sequence[ires] = NucleotideName[l][0];
              haveThis = true;
            }else if(strlen(NucleotideName[l])==2){
              sequence[ires] = NucleotideName[l][1];
              haveThis = true;
            }
          }
	}
      }
      if(!strcmp(resname,"HOH")){
	//	residue_comp[22]++;        // HARDWIRE!!
        //std::cout << 'O' << "(" << ires << ")\n";
	if(ires<numres) sequence[ires] = 'O'; // Hackery since we have WAT in RNames
        haveThis = true;
      }
      if(!haveThis&&ires<numres) sequence[ires] = 'X'; 
      ires++;
    }
  }

  //std::cout << sequence[numres-1] << "\n";
  sequence[ires] = '\0';

  std::cout <<  "Sequence length: " << strlen(sequence) << ", sequence--" << sequence << "\n";
  return sequence;

}

void CMMUTManager::PrintResidueComposition(int selHnd){

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int k;
  int *residue_comp;

  residue_comp = ResidueComposition(selHnd);

  cout << "\n          Residue Composition\n";
  cout << "          -------------------\n\n";

  for(k=0;k<nResNames;k++){
    if(residue_comp[k]>0){
      cout << ResidueName[k] << ": " << residue_comp[k] << endl;
    }
  }
}

void CMMUTManager::PrintSequence(int selHnd){

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  pstr sequence;
  int numres;

  numres = TotalNumRes(selHnd);

  sequence = GetSequence(selHnd);

  cout << "\n        Amino Acid Sequence\n";
  cout << "        -------------------\n";

  for(int i=0;i<numres;i++)
    cout << sequence[i];
  cout << endl;

}


std::vector<double> CMMUTManager::GetBValuesDoubleVector(int selHnd){
  realtype* bvals = GetBValues(selHnd);
  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);
  std::vector<double> bvals_vec;
  for(int i=0;i<natoms_sel;i++)
     bvals_vec.push_back(bvals[i]);
  return bvals_vec;
}

realtype* CMMUTManager::GetBValues(int selHnd){

  // Handling anisotropic bvalues disabled to get compile working
  // Cryst class is protested
  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  realtype b;
  PCAtom atom;
  int i,k;
  int natomsres;

  realtype *bvalues;

  //if(!isCrystInfo())
    //return NULL;

  //Cryst.CalcOrthMatrices();

  bvalues = new realtype[TotalNumRes(selHnd)];

  k = 0;
  b = 0.0;
  natomsres = 1;
  atom = selected_atoms[0];
  //if (atom->WhatIsSet & ASET_Anis_tFac){  // Calculate from anisotropic data if available.
  //  b += 8 / 3 * PI*PI *(atom->u11+atom->u22+atom->u33);
  //}else{
    b += atom->tempFactor;
    //}

  for(i=1;i<natoms_sel;i++){
    if((selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum())
      ||(
      (selected_atoms[i]->GetSeqNum()==selected_atoms[i-1]->GetSeqNum())&&((strlen(selected_atoms[i]->residue->insCode)!=strlen(selected_atoms[i-1]->residue->insCode))||(strncmp(selected_atoms[i]->residue->insCode,selected_atoms[i-1]->residue->insCode,strlen(selected_atoms[i]->residue->insCode))))
      )
      ){
      bvalues[k] = b / natomsres;
      k++;
      b = 0.0;
      natomsres = 0;
    }
    natomsres++;
    atom = selected_atoms[i];
    //if (atom->WhatIsSet & ASET_Anis_tFac){  // Calculate from anisotropic data if available.
    //b += 8 / 3 * PI*PI *(atom->u11+atom->u22+atom->u33);
    //}else{
      b += atom->tempFactor;
      //}
  }
  bvalues[k] = b / natomsres;

  return bvalues;

}


void CMMUTManager::PrintBValues(int selHnd){

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i;
  int nrestot = TotalNumRes(selHnd);
  realtype *bvalues;
  //  ofstream plot("bplot");

  if(!isCrystInfo()){
    cout << "Cannot calculate B-values without crystallographic data\n";
    return;
}

  bvalues = GetBValues(selHnd);

  cout << "\nAverage B-Value per Residue\n";
  cout << "-----------------------------\n\n";

  cout << "nrestot: " << nrestot << endl;

  /*
  plot.flags(ios::fixed|ios::right);
  plot.precision(6);
  */

  for(i=0;i<nrestot;i++){
    cout << "Residue Number: " << i << " Average B-Value: " << bvalues[i] << endl;
    //    plot << i << setw(12) << bvalues[i] << endl;
  }
}

int CMMUTManager::TotalNumRes(int selHnd){

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i,nrestot;

  nrestot = 1; 

  for(i=1;i<natoms_sel;i++){
    if(selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum()){
      nrestot++; 
    }else if((selected_atoms[i]->GetSeqNum()==selected_atoms[i-1]->GetSeqNum())&&((strlen(selected_atoms[i]->residue->insCode)!=strlen(selected_atoms[i-1]->residue->insCode))||(strncmp(selected_atoms[i]->residue->insCode,selected_atoms[i-1]->residue->insCode,strlen(selected_atoms[i]->residue->insCode))))){
      //std::cout << selected_atoms[i]->GetSeqNum() << " " << selected_atoms[i]->residue->insCode << " does not match " << selected_atoms[i-1]->residue->insCode << "\n";
      nrestot++; 
    }
  }
  
  return nrestot;

}


realtype CMMUTManager::BondLength(PCAtom A, PCAtom B){

  realtype length;

  realtype x = (A->x-B->x);
  realtype y = (A->y-B->y);
  realtype z = (A->z-B->z);

  length= sqrt(x*x + y*y + z*z);
  
  return length;
}

realtype CMMUTManager::BondAngle(PCAtom A, PCAtom B, PCAtom C){

  realtype abx = (B->x-A->x);
  realtype aby = (B->y-A->y);
  realtype abz = (B->z-A->z);

  realtype acx = (A->x-C->x);
  realtype acy = (A->y-C->y);
  realtype acz = (A->z-C->z);

  realtype bcx = (C->x-B->x);
  realtype bcy = (C->y-B->y);
  realtype bcz = (C->z-B->z);

  realtype absq = abx*abx + aby*aby + abz*abz;
  realtype acsq = acx*acx + acy*acy + acz*acz;
  realtype bcsq = bcx*bcx + bcy*bcy + bcz*bcz;

  realtype ab = sqrt(absq);
  realtype bc = sqrt(bcsq);

  return  acos((bcsq + absq - acsq)/(2*bc*ab));
}

realtype CMMUTManager::TorsionAngle(PCAtom A, PCAtom B, PCAtom C, PCAtom D){

  /* Should really move all this vector manipulation stuff 
   * somewhere more general. In fact should use vectors ... */

  realtype abx = (B->x-A->x);      // AB
  realtype aby = (B->y-A->y);
  realtype abz = (B->z-A->z);

  realtype bcx = (C->x-B->x);      // BC 
  realtype bcy = (C->y-B->y);
  realtype bcz = (C->z-B->z);

  realtype cdx = (D->x-C->x);      // CD
  realtype cdy = (D->y-C->y);
  realtype cdz = (D->z-C->z);

  realtype qx = aby*bcz - bcy*abz;  // Q = AB X BC
  realtype qy = abz*bcx - bcz*abx;
  realtype qz = abx*bcy - bcx*aby;

  realtype tx = bcy*cdz - cdy*bcz;  // T = BC X CD
  realtype ty = bcz*cdx - cdz*bcx;
  realtype tz = bcx*cdy - cdx*bcy;

  realtype sx = qy*tz - ty*qz;      // S = Q X T
  realtype sy = qz*tx - tz*qx;
  realtype sz = qx*ty - tx*qy;

  realtype ss = sx*sx + sy*sy + sz*sz; // SS = S . S
  realtype qt = qx*tx + qy*ty + qz*tz; // QT = Q . T

  realtype angle = atan2(sqrt(ss),qt);

  realtype sbc = sx*bcx + sy*bcy + sz*bcz;

  if(sbc<0.0) {
    return -angle;
  } else {
    return angle;
  }

}

//----------------------------------------------------------------------------
Boolean CMMUTManager::isMainChain(PCAtom p_atom) {
//----------------------------------------------------------------------------
  char *mainchAtoms[5] = { "CA", "N", "C", "O", "HA" };
  if ( NameComparison(p_atom->name,5,mainchAtoms) >= 0 ) {
    return true;
  }
  else {
    return false;
  } 
}
//-----------------------------------------------------------------
Boolean CMMUTManager::doAltLocMatch ( PCAtom pa1, PCAtom pa2 ) {
//-----------------------------------------------------------------
 if ( strlen (pa1->altLoc) == 0 ||
        strlen (pa2->altLoc) == 0 ||
      strcmp ( pa1->altLoc, pa2->altLoc ) == 0 )
   return 1;
 else
    return 0;
}


//---------------------------------------------------------------------------------
int CMMUTManager::NameComparison ( const char *name , int ntypes , char *types[] ) {
//---------------------------------------------------------------------------------
  //Compare an fixed length char atom/residue/whatever name to a list of
  // variable length strings.  Return the position in the list of any match
  // Otherwise return -1
  // Allow for an initial blank character in the name (as in PDB atom name) 
  int n,m,il1,if1,strlength;

  if1 = 0;
  if (name[0] == ' ') if1 = 1;

  strlength = (int)strlen(name);
  il1 = strlength - if1;
  for(n=1;n<strlength;n++){
    if(name[n]==' ') {
      il1 = n-if1;
      break;
    }
  }
 
  //printf("strlength *%s* %i %i\n",name,if1,il1);

  for(m=0;m<ntypes;m++){
    int lenm = int(strlen(types[m]));
    if(lenm==il1){
      if(!strncmp(types[m],&name[if1],il1)){
        return m;
      }
    }
  }
  return -1;
}

//-------------------------------------------------------------------
std::string  CMMUTManager::TrimString(pstr input) {
//-------------------------------------------------------------------
  std::string inp = input;
  int f = inp.find(" ");
  while ( f >= 0 ) {
    inp.replace(f,1,"");
    f = inp.find(" ");
  }
  return inp;
}

//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_atom1(PCAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,1,1,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_atom(PCAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,1,1,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_residue(PCAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,0,0,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}
//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_chain(PCAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,0,0,0,0,0,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
Boolean CMMUTManager::ChainIDisDigit(PCChain p_ch) {
//--------------------------------------------------------------------
  char *digits[10] = { "0", "1", "2", "3", "4","5", "6", "7", "8", "9"  };
  if ( NameComparison(p_ch->GetChainID(),10,digits) >= 0 ) {
    return true;
  }
  else {
    return false;
  } 

}

//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_residue1(PCResidue p_res) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,0,0,0,0 };
  if (GetNumberOfModels()>1 || ChainIDisDigit(p_res->chain) ) mask[1]=1;
  PCAtom p_atom=p_res->GetAtom(0); 
  const char *al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
 const char* CMMUTManager::AtomLabel_mask(PCAtom p_atom, int mask[]) {
//--------------------------------------------------------------------
  const char *al = AtomLabel(p_atom,mask).c_str();
  return al;
}
 
//--------------------------------------------------------------------
std::string CMMUTManager::AtomLabel(PCAtom p_atom, int mask[]) {
//--------------------------------------------------------------------
  //Mask parameters
  //  0  Molecule
  //  1  Model
  //  2  Chain Id
  //  3  Sequence number
  //  4  Insertion code 
  //  5  Residue name
  //  6  Atom id
  //  7  Alternate position
  //  8  Element

  int masksum = 0;
  std::ostringstream label;

  //cout << "mask" << mask[0] << mask[1] << mask[2] << endl;
  //if (mask[0]){
  //}
  if (mask[1]) {
    label << "/";
    label << p_atom->GetModelNum();   
  }
  masksum = masksum + mask[1];
  if (mask[2]) {
    if (masksum || strlen(p_atom->GetChainID())==0 ) label << "/";
    //if (strlen(p_atom->GetChainID())==0 ) {
    //  label << " ";
    //} else {
      label << p_atom->GetChainID();
      //}
  }
  masksum = masksum + mask[2];
  if (masksum) label << "/";

  if (mask[3]>0) label << p_atom->GetSeqNum();
  if (mask[4] && strlen(p_atom->residue->insCode ) != 0 )   
          label << "." << p_atom->residue->insCode ;

  if (mask[5]>0) label << "(" << p_atom->GetResName() << ")";


  masksum = masksum + mask[3]+mask[5];
  if (mask[6] ) {
    if (masksum) label << "/";
    label << TrimString(p_atom->name);
    if ( mask[7] && strlen(p_atom->altLoc) != 0  )
      label << ":" << p_atom->altLoc; 
  }
  if ( mask[8] ) {
    label << "[" << p_atom->element << "]";
  }

  
  return label.str();

  // Return as pstr - may not be necessary
  //int len = label.length();
  //pstr chlabel;
  //label.copy(chlabel,len,0);
  //chlabel[len] = 0;
  //return label;

}

/*
int CMMUTManager::ApplyTransform(int selHnd,double rotmat[],
                                              double transv[]) {
  mat33 rot;
  vect3 tran;
  int i,j,ii;
  PPCAtom selAtoms;
  int nAtoms;

  if ( selHnd >= 0 ) 
    GetSelIndex ( selHnd, selAtoms, nAtoms);
  else
    GetAtomTable(selAtoms,nAtoms);

  //cout << "ApplyTransform " << nAtoms << endl;
  //cout << "rotmat " << rotmat[0] << " " << rotmat[1] << endl;
  
  for ( i=0; i<3; i++) tran[i]=transv[i];
  ii = 0;
  for ( i=0; i<3; i++) {  
    for (j=0; j<3; j++) rot[i][j] = rotmat[ii++];
  }
   
  for (int i = 0; i < nAtoms; i++) {
    selAtoms[i]->Transform(rot,tran);
  }  
  return 0;
}
*/


//-------------------------------------------------------------------------
int CMMUTManager::WriteSelection (int selHnd,char *file, char *format) {
//-------------------------------------------------------------------------
  PCMMUTManager mmdb2;
  PPCAtom selAtom;
  int nSelAtoms,i,RC;

  mmdb2 = new CMMUTManager();

    //         It looks like a good idea to copy the header and
    //         crystallographic information into output file
    //         as well. However this is not needed by MMDB to
    //         function. Check the copy function and its parameters
    //         in file  mmdb_manager.h .

  mmdb2->Copy (this ,MMDBFCM_Title | MMDBFCM_Cryst );

    //         Stuffing an empty MMDB with atoms. Check the
    //         corresponding function and its parameters in file
    //         mmdb_file.h .
  
  if (GetNumberOfModels()>1) {
    CopySelection(selHnd,mmdb2);
  } else {
  
    GetSelIndex ( selHnd, selAtom, nSelAtoms);
    for (i=0;i<nSelAtoms;i++)
      mmdb2->PutAtom ( i+1,selAtom[i] ); 
    }

  //         Purge the newly composed MMDB into the file
  //printf ( " ...writing PDB file %s\n",file );
  if (strcmp(format,"PDB")==0)
    RC = mmdb2->WritePDBASCII( file );
  else
    RC = mmdb2->WriteCIFASCII( file );
  delete mmdb2;
  return RC;
}
//-------------------------------------------------------------------------
int CMMUTManager::PutSelectedAtoms (int selHnd , const PCMMDBManager mmdb2) {
//-------------------------------------------------------------------------
  PPCAtom selAtom;
  int nSelAtoms,i;
  GetSelIndex ( selHnd, selAtom, nSelAtoms);
  //cout << "nSelAtoms " << nSelAtoms << endl;
  for (i=0;i<nSelAtoms;i++)
    mmdb2->PutAtom ( i+1,selAtom[i] );
  
  return 0;
}


//-------------------------------------------------------------------------
int CMMUTManager::CopySelection (int selHnd , const PCMMDBManager mmdb2 ) {
//-------------------------------------------------------------------------
  int          RC,i,j,k;
  int          nChains,nResidues,nAtoms;
  PPCChain     chain;
  PPCResidue   res;
  PPCAtom      atom;

  PCResidue newRes;
  PCChain newChain;
  PCModel newModel;

  int nAcopied=0, nRcopied=0, nCcopied=0, nMcopied=0;

  for (int model=1;model<GetNumberOfModels();model++) {
    GetChainTable ( model,chain,nChains );


    for (i=0;i<nChains;i++) {
      if (chain[i])  {
        chain[i]->GetResidueTable ( res,nResidues );

        for (j=0;j<nResidues;j++) {
          if (res[j])  {          
            res[j]->GetAtomTable ( atom,nAtoms );
            
            for (k=0;k<nAtoms;k++) {
              if (atom[k] && atom[k]->isInSelection(selHnd))  {
                if(nAcopied==0) {
                   newRes = new CResidue();
                   newRes->SetResID(res[j]->name,res[j]->seqNum, \
                            res[j]->GetInsCode());
                }
                nAcopied++;
                RC = newRes->AddAtom ( atom[k] );
              }
            }
            if (nAcopied>0) {
              if (nRcopied == 0) {
                newChain = new CChain();
                newChain->SetChainID(chain[i]->GetChainID());
              }
              newChain->AddResidue(newRes);
              nAcopied=0;
              nRcopied++;
            }
          }
        }
        if (nRcopied>0) {
          if (nCcopied==0) newModel = new CModel();
          newModel->AddChain(newChain);
          nRcopied = 0;
          nCcopied++;
        }
      }
    }
    if (nCcopied>0) {
      mmdb2->AddModel ( newModel );
      nCcopied=0;
      nMcopied++;
    }
  }

  return 0;
}

//-----------------------------------------------------------------------
PCResidue CMMUTManager::NextResidue( PCResidue pRes ,int increment) {
//-----------------------------------------------------------------------
  PCChain pc = pRes->GetChain();
  int resno = pc->GetResidueNo(pRes->GetSeqNum(),pRes->GetInsCode());
  if (resno < pc->GetNumberOfResidues()) {
    return pc->GetResidue(resno+increment);
  } else {
    return NULL;
  }
}


//-----------------------------------------------------------------------
int CMMUTManager::FindCloseAtomPairs ( int selHnd, double min_distance, 
      double max_distance) {
//-----------------------------------------------------------------------
// Find all atoms in set selHnd  within min_distance to max_distance of
// another member of the set and return a seelction handle that lists them
// Originally implemented to help fing SSbonds
  PPCAtom selAtoms,selAtoms2;
  int nSelAtoms,nSelAtoms2;
  PSContact contact = NULL; 
  int ncontacts; 

  int closeHnd = NewSelection();
  GetSelIndex ( selHnd, selAtoms, nSelAtoms);
  GetSelIndex ( selHnd, selAtoms2, nSelAtoms2);

  if(nSelAtoms<1||nSelAtoms2<1){
    return closeHnd;
  }

  SeekContacts ( selAtoms, nSelAtoms,  selAtoms2, nSelAtoms2, min_distance,
        max_distance , 0, contact,ncontacts);

  for (int n=0;n<ncontacts;n++) {
    if (doAltLocMatch( selAtoms[contact[n].id1],selAtoms2[contact[n].id2]) ) {
      SelectAtoms(closeHnd,selAtoms[contact[n].id1]->serNum,
                   selAtoms[contact[n].id1]->serNum,SKEY_OR);
      SelectAtoms(closeHnd,selAtoms2[contact[n].id2]->serNum,
                   selAtoms2[contact[n].id2]->serNum,SKEY_OR);
    }
  }

  return closeHnd;
  
}

const std::string str_trim(const std::string& pString, const std::string& pWhitespace = " \t")
{
    const size_t beginStr = pString.find_first_not_of(pWhitespace);
    if (beginStr == std::string::npos)
    {
        return "";
    }

    const size_t endStr = pString.find_last_not_of(pWhitespace);
    const size_t range = endStr - beginStr + 1;

    return pString.substr(beginStr, range);
}

int CMMUTManager::FixElementNames(){
  //clock_t t1 = clock();
  /*
 * This attempts to fix atoms names such as "  2HB "
 * which should be " 2HB  ".
  */
  int selHnd = NewSelection();
  SelectAtoms(selHnd,0,"*", ANY_RES,"*", ANY_RES,"*","*","*","*","*",SKEY_NEW);
  int nSelAtoms;
  PPCAtom SelAtoms=NULL;
  GetSelIndex(selHnd,SelAtoms,nSelAtoms);
  std::string numbers = "1234567890";
  int fixed_some = 0;
 
  for(int i=0;i<nSelAtoms;i++){
    if(getElementNo(SelAtoms[i]->element)==ELEMENT_UNKNOWN){
      ++fixed_some;
      // We do not know what this element is ...
      std::string trimmed = str_trim(std::string(SelAtoms[i]->element));
      if(strlen(SelAtoms[i]->name)>1&&trimmed.size()>0&&trimmed.substr(0,1).find_first_of(numbers)!=std::string::npos){
        // It begins with whitespace and some numbers so we try removing first character ...
        SelAtoms[i]->SetAtomName(SelAtoms[i]->name+1);
        strcpy(SelAtoms[i]->element,"  ");
        SelAtoms[i]->MakePDBAtomName();
        if(getElementNo(SelAtoms[i]->element)==ELEMENT_UNKNOWN){
          // That failed, so we try to find first non-whitespace, non-numeric character
          std::string name_trimmed = str_trim(std::string(SelAtoms[i]->name));
          if(name_trimmed.size()>0&&name_trimmed.find_first_not_of(numbers)!=std::string::npos){
             std::string firstLetter = name_trimmed.substr(name_trimmed.find_first_not_of(numbers),1);
             std::string tryElement = " "+ firstLetter;
             strcpy(SelAtoms[i]->element,tryElement.c_str());
          }
          if(getElementNo(SelAtoms[i]->element)==ELEMENT_UNKNOWN){
            // Even that failed, call it a C...
            strcpy(SelAtoms[i]->element," C");
          }
        }
      } else {
        std::string orig_name = std::string(SelAtoms[i]->name);
        SelAtoms[i]->SetAtomName((std::string(" ") + std::string(SelAtoms[i]->name).substr(0,10)).c_str());
        strcpy(SelAtoms[i]->element,"  ");
        SelAtoms[i]->MakePDBAtomName();
        // Then put original name back, potentially hazardous, I guess.
        SelAtoms[i]->SetAtomName(orig_name.c_str());
        if(getElementNo(SelAtoms[i]->element)==ELEMENT_UNKNOWN){
          strcpy(SelAtoms[i]->element," C");
        }
      }
    } else if(SelAtoms[i]->GetResidue()->isAminoacid()&&SelAtoms[i]->element[0]=='H'){
      std::cout << "Improbable He, Hf, Hg or Ho atom in amino acid (" << SelAtoms[i]->element << "), changing to H\n";
      strcpy(SelAtoms[i]->element," H");
    }
  }

  DeleteSelection(selHnd);
  //clock_t t2 = clock();
  //std::cout << "Time to fix element names " << ((t2-t1)*1000.0/CLOCKS_PER_SEC)/1000.0<< "\n";
  return fixed_some;
}

double CMMUTManager::Resolution(){
  double resolution = GetResolution();
  return resolution;
}

std::string CMMUTManager::StructureTitle(){
  
  pstr title_s=NULL;
  GetStructureTitle(title_s);
  if(!S) return std::string("");
  std::string title(title_s);
  return title;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

struct is_remark_800 : public std::unary_function<std::string, bool>
{
    bool operator() (const string & inValue) const { return inValue.compare(0,10,std::string("REMARK 800")); }
};

struct is_site : public std::unary_function<std::string, bool>
{
    bool operator() (const string & inValue) const { return inValue.compare(0,6,std::string("SITE  ")); }
};

std::string CMMUTManager::SiteInfo(){
  std::string s;

  
  std::string remarks = GetRemarksString();
  std::vector<std::string> remarks800 = split(remarks,'\n');
  std::vector<std::string>::iterator remarks800_iter = std::remove_if(remarks800.begin(),remarks800.end(),is_remark_800());
  remarks800.erase(remarks800_iter,remarks800.end());
  for(unsigned i=0;i<remarks800.size();i++){
    s += remarks800[i]+std::string("\n");
  }

  std::string unknowns = Unknowns();
  std::vector<std::string> site = split(unknowns,'\n');
  std::vector<std::string>::iterator site_iter = std::remove_if(site.begin(),site.end(),is_site());
  site.erase(site_iter,site.end());
  for(unsigned i=0;i<site.size();i++){
    s += site[i]+std::string("\n");
  }


  return s;
}

std::string CMMUTManager::GetRemarksString(){
  std::string s;
  char *tempfilename = tmpnam(0);
  CFile f;
  f.assign(tempfilename,True,False,False);
  f.rewrite();
  GetRemarks()->PDBASCIIDump(f);
  f.shut();

  std::ifstream ifstr(tempfilename);
  if(ifstr){
    std::stringstream buffer;
    buffer << ifstr.rdbuf();
    ifstr.close();
    return buffer.str();
  }
  
  return s;
}

std::string CMMUTManager::Unknowns(){
  std::string s;
  char *tempfilename = tmpnam(0);
  CFile f;
  f.assign(tempfilename,True,False,False);
  f.rewrite();
  SA.PDBASCIIDump(f); 
  f.shut();

  std::ifstream ifstr(tempfilename);
  if(ifstr){
    std::stringstream buffer;
    buffer << ifstr.rdbuf();
    ifstr.close();
    return buffer.str();
  }

  return s;
}

std::string CMMUTManager::Source(){

  pstr S=NULL;
  PCMMCIFData CIF =  new CMMCIFData();
  Title.MakeCIF(CIF);
  CIF->GetString(S,CIFCAT_STRUCT,CIFTAG_SOURCE,False);
  if(!S) return std::string("");
  std::string source(S);
  if(CIF) delete CIF;

  return source;
}

std::vector<double> CMMUTManager::GetCellInfo() {
  std::vector<double> cell_info;

  realtype a,b,c,alpha,beta,gamma,vol;
  int orth_code;

  GetCell(a,b,c,alpha,beta,gamma,vol,orth_code);

  cell_info.push_back(a);
  cell_info.push_back(b);
  cell_info.push_back(c);
  cell_info.push_back(alpha);
  cell_info.push_back(beta);
  cell_info.push_back(gamma);

  return cell_info;
}

std::string CMMUTManager::MMUTGetSpaceGroup(){
  PCMMDBCryst my_cryst_p = (CMMDBCryst *) &(get_cell());
  return std::string(my_cryst_p->spaceGroup);
}

int CMMUTManager::ApplyPDBSecStructure(int model){
  // Need to consider all the fields in
  // http://www.wwpdb.org/documentation/format33/sect5.html
  // as long as MMDB supports them.
  int nhelix = GetModel(model)->GetNumberOfHelices();
  int nsheets = GetModel(model)->GetNumberOfSheets();
  //std::cout << nhelix << " " << nsheets << "\n";
  for(int i=1;i<=nhelix;i++){
    PCHelix helix =GetModel(model)->GetHelix(i);
    if(helix){
      //std::cout << helix->initChainID << " " << helix->initSeqNum << " " << helix->endChainID << " " << helix->endSeqNum << "\n";
      if(strncmp(helix->initChainID,helix->endChainID,1)==0){
        for(int k=helix->initSeqNum;k<=helix->endSeqNum;k++){
          PCResidue res = GetModel(model)->GetResidue(helix->initChainID,k,helix->initICode);
          //std::cout << res->GetSeqNum() << "\n";
          //std::cout << helix->helixClass << "\n";
	  if(res&&(helix->helixClass==1||helix->helixClass==6||helix->helixClass==5)){
            res->SSE = SSE_Helix;
          }
        }
      }
    }
  }
  for(int i=1;i<=nsheets;i++){
    PCSheet sheet = GetModel(model)->GetSheet(i);
    if(sheet){
      //std::cout << sheet->nStrands << "\n"; std::cout.flush();
      for(int j=0;j<sheet->nStrands;j++){
        if(sheet->Strand[j]){
          //std::cout << sheet->Strand[j]->initChainID << " " << sheet->Strand[j]->initSeqNum << " " << sheet->Strand[j]->endChainID << " " << sheet->Strand[j]->endSeqNum << "\n"; std::cout.flush();
          if(strncmp(sheet->Strand[j]->initChainID,sheet->Strand[j]->endChainID,1)==0){
            for(int k=sheet->Strand[j]->initSeqNum;k<=sheet->Strand[j]->endSeqNum;k++){
              PCResidue res = GetModel(model)->GetResidue(sheet->Strand[j]->initChainID,k,sheet->Strand[j]->initICode);
              if(res)
                res->SSE = SSE_Strand;
            }
          }
        }
      }
    }
  }
  return nhelix+nsheets;
}


