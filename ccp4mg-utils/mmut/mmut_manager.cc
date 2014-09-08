/*
     mmut/mmut_manager.cc: CCP4MG Molecular Graphics Program
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

#if defined(__sgi) || defined(sgi) || defined(__OSF1__) || defined(__osf__)
#include <string.h>
#else
#include <cstring>
#endif
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <mmut_manager.h>
#include <mmdb_utils.h>
#include <mmdb_tables.h>

using namespace std;

//  =====================   CMMUTManager   =======================

 
  CMMUTManager::CMMUTManager() : mmdb::Manager()  {

    iatom_types=NULL;
    iatom_type_lookup=NULL;

  }


  CMMUTManager::~CMMUTManager()  { 
   
  }


int const nhydro[nResNames] = { // No of hydrogen atoms in the 27 residues.
5,13,6,4,5,5,8,6,3,8,11,11,13,9,9,7,5,7,10,9,9,17,2,0,0,2
};



void CMMUTManager::ListAtomInfo (int selHnd) {

  mmdb::PAtom atom;
  int i;

  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  mmdb::realtype atomweight, molweight=0.0;


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

mmdb::realtype CMMUTManager::MolWeight(int selHnd){

  int l;

  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  mmdb::realtype molweight=0.0;

  for(l=0;l<natoms_sel;l++){
    molweight += MolecWeight[getElementNo((selected_atoms[l])->element)];
  }
  return molweight;
}

mmdb::realtype CMMUTManager::MolWeightWithH(int selHnd){

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

  mmdb::PPAtom selected_atoms;
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
mmdb::realtype* CMMUTManager::CentreOfMass(int selHnd) {
//----------------------------------------------------------------------------

  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int l;
  mmdb::realtype *com;
  mmdb::realtype atweight;
  mmdb::realtype weight;

  //weight=MolWeight(selHnd);
  weight = 0.0;
  com = new mmdb::realtype[3];

  com[0] = 0.0;
  com[1] = 0.0;
  com[2] = 0.0;

  atweight = 1.0;
  for(l=0;l<natoms_sel;l++){
    //atweight = MolecWeight[getElementNo((selected_atoms[l])->element)];
    com[0] +=  atweight * (selected_atoms[l])->x;
    com[1] +=  atweight * (selected_atoms[l])->y;
    com[2] +=  atweight * (selected_atoms[l])->z;
    weight += 1.0;
  }

  com[0] = com[0]/weight;
  com[1] = com[1]/weight;
  com[2] = com[2]/weight;

  return com;

}

//-------------------------------------------------------------------------
mmdb::realtype* CMMUTManager::Extent(int selHnd) {
//-------------------------------------------------------------------------

  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  if (natoms_sel < 1) return NULL;

  int l;
  mmdb::realtype *com;
  com = new mmdb::realtype[6];
  mmdb::realtype xmin = 9999999.9;
  mmdb::realtype ymin = 9999999.9;
  mmdb::realtype zmin = 9999999.9;
  mmdb::realtype xmax = -9999999.9;
  mmdb::realtype ymax = -9999999.9;
  mmdb::realtype zmax = -9999999.9;


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

int* CMMUTManager::AtomicComposition(int selHnd){
	
  int l;

  mmdb::PPAtom selected_atoms;
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

  mmdb::PPAtom selected_atoms;
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
	
  mmdb::PPAtom selected_atoms;
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

pstr CMMUTManager::GetSequence(int selHnd){
	
  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i,l;
  int numres;
  cpstr resname;
  pstr sequence;

  numres = TotalNumRes(selHnd);
  sequence = new char[numres];

  resname = selected_atoms[0]->GetResName();
  for(l=0;l<nResNames;l++){
    if(!strcmp(resname,ResidueName[l])){
      //      residue_comp[l]++;
      sequence[0] = ResidueName1[l];
    }
  }
  if(!strcmp(resname,"HOH")){
    //    residue_comp[22]++;        // HARDWIRE!!
    sequence[0] = 'O'; // Hackery since we have WAT in RNames
  }

  for(i=1;i<natoms_sel;i++){
    if(selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum()){
      resname = selected_atoms[i]->GetResName();
      for(l=0;l<nResNames;l++){
	if(!strcmp(resname,ResidueName[l])){
	  //	  residue_comp[l]++;
	  sequence[selected_atoms[i]->GetSeqNum()-1] = ResidueName1[l];
	}
      }
      if(!strcmp(resname,"HOH")){
	//	residue_comp[22]++;        // HARDWIRE!!
	sequence[selected_atoms[i]->GetSeqNum()-1] = 'O'; // Hackery since we have WAT in RNames
      }
    }
  }

  sequence[numres] = '\0';

  return sequence;

}

void CMMUTManager::PrintResidueComposition(int selHnd){

  mmdb::PPAtom selected_atoms;
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

  mmdb::PPAtom selected_atoms;
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


mmdb::realtype* CMMUTManager::GetBValues(int selHnd){

  // Handling anisotropic bvalues disabled to get compile working
  // Cryst class is protested
  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  mmdb::realtype b;
  mmdb::PAtom atom;
  int i,k;
  int natomsres;

  mmdb::realtype *bvalues;

  if(!isCrystInfo())
    return NULL;

  //Cryst.CalcOrthMatrices();

  bvalues = new mmdb::realtype[TotalNumRes(selHnd)];

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
    if(selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum()) {
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

  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i;
  int nrestot = TotalNumRes(selHnd);
  mmdb::realtype *bvalues;
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

  mmdb::PPAtom selected_atoms;
  int natoms_sel;
  GetSelIndex ( selHnd, selected_atoms, natoms_sel);

  int i,nrestot;

  nrestot = 1; 

  for(i=1;i<natoms_sel;i++)
    if(selected_atoms[i]->GetSeqNum()!=selected_atoms[i-1]->GetSeqNum()) nrestot++; 
  
  return nrestot;

}


mmdb::realtype CMMUTManager::BondLength(mmdb::PAtom A, mmdb::PAtom B){

  mmdb::realtype length;

  mmdb::realtype x = (A->x-B->x);
  mmdb::realtype y = (A->y-B->y);
  mmdb::realtype z = (A->z-B->z);

  length= sqrt(x*x + y*y + z*z);
  
  return length;
}

mmdb::realtype CMMUTManager::BondAngle(mmdb::PAtom A, mmdb::PAtom B, mmdb::PAtom C){

  mmdb::realtype abx = (B->x-A->x);
  mmdb::realtype aby = (B->y-A->y);
  mmdb::realtype abz = (B->z-A->z);

  mmdb::realtype acx = (A->x-C->x);
  mmdb::realtype acy = (A->y-C->y);
  mmdb::realtype acz = (A->z-C->z);

  mmdb::realtype bcx = (C->x-B->x);
  mmdb::realtype bcy = (C->y-B->y);
  mmdb::realtype bcz = (C->z-B->z);

  mmdb::realtype absq = abx*abx + aby*aby + abz*abz;
  mmdb::realtype acsq = acx*acx + acy*acy + acz*acz;
  mmdb::realtype bcsq = bcx*bcx + bcy*bcy + bcz*bcz;

  mmdb::realtype ab = sqrt(absq);
  mmdb::realtype bc = sqrt(bcsq);

  return  acos((bcsq + absq - acsq)/(2*bc*ab));
}

mmdb::realtype CMMUTManager::TorsionAngle(mmdb::PAtom A, mmdb::PAtom B, mmdb::PAtom C, mmdb::PAtom D){

  /* Should really move all this vector manipulation stuff 
   * somewhere more general. In fact should use vectors ... */

  mmdb::realtype abx = (B->x-A->x);      // AB
  mmdb::realtype aby = (B->y-A->y);
  mmdb::realtype abz = (B->z-A->z);

  mmdb::realtype bcx = (C->x-B->x);      // BC 
  mmdb::realtype bcy = (C->y-B->y);
  mmdb::realtype bcz = (C->z-B->z);

  mmdb::realtype cdx = (D->x-C->x);      // CD
  mmdb::realtype cdy = (D->y-C->y);
  mmdb::realtype cdz = (D->z-C->z);

  mmdb::realtype qx = aby*bcz - bcy*abz;  // Q = AB X BC
  mmdb::realtype qy = abz*bcx - bcz*abx;
  mmdb::realtype qz = abx*bcy - bcx*aby;

  mmdb::realtype tx = bcy*cdz - cdy*bcz;  // T = BC X CD
  mmdb::realtype ty = bcz*cdx - cdz*bcx;
  mmdb::realtype tz = bcx*cdy - cdx*bcy;

  mmdb::realtype sx = qy*tz - ty*qz;      // S = Q X T
  mmdb::realtype sy = qz*tx - tz*qx;
  mmdb::realtype sz = qx*ty - tx*qy;

  mmdb::realtype ss = sx*sx + sy*sy + sz*sz; // SS = S . S
  mmdb::realtype qt = qx*tx + qy*ty + qz*tz; // QT = Q . T

  mmdb::realtype angle = atan2(sqrt(ss),qt);

  mmdb::realtype sbc = sx*bcx + sy*bcy + sz*bcz;

  if(sbc<0.0) {
    return -angle;
  } else {
    return angle;
  }

}

//----------------------------------------------------------------------------
bool CMMUTManager::isMainChain(mmdb::PAtom p_atom) {
//----------------------------------------------------------------------------
   const char *mainchAtoms[5] = { "CA", "N", "C", "O", "HA" };
   if ( NameComparison(p_atom->name,5,mainchAtoms) >= 0 ) {
    return true;
  }
  else {
    return false;
  } 
}
//-----------------------------------------------------------------
bool CMMUTManager::doAltLocMatch ( mmdb::PAtom pa1, mmdb::PAtom pa2 ) {
//-----------------------------------------------------------------
 if ( strlen (pa1->altLoc) == 0 ||
        strlen (pa2->altLoc) == 0 ||
      strcmp ( pa1->altLoc, pa2->altLoc ) == 0 )
   return 1;
 else
    return 0;
}


//---------------------------------------------------------------------------------
int CMMUTManager::NameComparison ( const char *name , int ntypes , const char *types[] ) {
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
const char* CMMUTManager::AtomLabel_atom1(mmdb::PAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,1,1,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_atom(mmdb::PAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,1,1,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_residue(mmdb::PAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,0,0,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}
//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_chain(mmdb::PAtom p_atom) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,0,0,0,0,0,0,0 };
  if (GetNumberOfModels()>1) mask[1]=1;
  const char* al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
bool CMMUTManager::ChainIDisDigit(mmdb::PChain p_ch) {
//--------------------------------------------------------------------
   const char *digits[10] = { "0", "1", "2", "3", "4","5", "6", "7", "8", "9"  };
  if ( NameComparison(p_ch->GetChainID(),10,digits) >= 0 ) {
    return true;
  }
  else {
    return false;
  } 

}

//--------------------------------------------------------------------
const char* CMMUTManager::AtomLabel_residue1(mmdb::PResidue p_res) {
//--------------------------------------------------------------------
  int mask[10] = { 0,0,1,1,1,1,0,0,0,0 };
  if (GetNumberOfModels()>1 || ChainIDisDigit(p_res->chain) ) mask[1]=1;
  mmdb::PAtom p_atom=p_res->GetAtom(0); 
  const char *al = AtomLabel(p_atom,mask).c_str();
  return al;
}

//--------------------------------------------------------------------
 const char* CMMUTManager::AtomLabel_mask(mmdb::PAtom p_atom, int mask[]) {
//--------------------------------------------------------------------
  const char *al = AtomLabel(p_atom,mask).c_str();
  return al;
}
 
//--------------------------------------------------------------------
std::string CMMUTManager::AtomLabel(mmdb::PAtom p_atom, int mask[]) {
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
  mmdb::PPAtom selAtoms;
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
  mmdb::PPAtom selAtom;
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
int CMMUTManager::PutSelectedAtoms (int selHnd , const mmdb::PManager mmdb2) {
//-------------------------------------------------------------------------
  mmdb::PPAtom selAtom;
  int nSelAtoms,i;
  GetSelIndex ( selHnd, selAtom, nSelAtoms);
  //cout << "nSelAtoms " << nSelAtoms << endl;
  for (i=0;i<nSelAtoms;i++)
    mmdb2->PutAtom ( i+1,selAtom[i] );
  
  return 0;
}


//-------------------------------------------------------------------------
int CMMUTManager::CopySelection (int selHnd , const mmdb::PManager mmdb2 ) {
//-------------------------------------------------------------------------
  int          RC,i,j,k;
  int          nChains,nResidues,nAtoms;
  Pmmdb::PChain     chain;
  mmdb::PPResidue   res;
  mmdb::PPAtom      atom;

  mmdb::PResidue newRes = 0;
  mmdb::PChain newChain = 0;
  mmdb::PModel newModel = 0;

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
                   newRes = new mmdb::Residue();
                   newRes->SetResID(res[j]->name,res[j]->seqNum, \
                            res[j]->GetInsCode());
                }
                nAcopied++;
                RC = newRes->AddAtom ( atom[k] );
              }
            }
            if (nAcopied>0) {
              if (nRcopied == 0) {
                newChain = new mmdb::Chain();
                newChain->SetChainID(chain[i]->GetChainID());
              }
              newChain->AddResidue(newRes);
              nAcopied=0;
              nRcopied++;
            }
          }
        }
        if (nRcopied>0) {
          if (nCcopied==0) newModel = new mmdb::Model();
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
mmdb::PResidue CMMUTManager::NextResidue( mmdb::PResidue pRes ,int increment) {
//-----------------------------------------------------------------------
  mmdb::PChain pc = pRes->GetChain();
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
  mmdb::PPAtom selAtoms,selAtoms2;
  int nSelAtoms,nSelAtoms2;
  mmdb::Contact *contact = NULL; 
  int ncontacts; 

  GetSelIndex ( selHnd, selAtoms, nSelAtoms);
  GetSelIndex ( selHnd, selAtoms2, nSelAtoms2);
  SeekContacts ( selAtoms, nSelAtoms,  selAtoms2, nSelAtoms2, min_distance,
        max_distance , 0, contact,ncontacts);

  int closeHnd = NewSelection();
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
