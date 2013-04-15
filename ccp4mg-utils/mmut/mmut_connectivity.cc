/*
     mmut/mmut_connectivity.cc: CCP4MG Molecular Graphics Program
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


#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include <string.h>

#include "connect.h"
#include "cartesian.h"

#include <mmdb/mmdb_manager.h>
#include <mmdb/mmdb_tables.h>
#include "mman_manager.h"
#include "mmut_connectivity.h"
#include "mman_base.h"
#include "mgutil.h"

//#include "mgtree.h"
//#include <atom_util.h>

using namespace std;

void Connection::AddConnection(int index){
  connected_atoms.push_back(index);
}

void Connection::AddExternalConnection(int index){
  //std::cout << "Pushing back external " << index << "\n";
  external_connected_atoms.push_back(index);
}

void Connection::AddExternalSplineConnection(int index){
  //std::cout << "Pushing back external " << index << "\n";
  external_connected_spline_atoms.push_back(index);
}

Connection::Connection(){
  isInTotalSelection = false;
}

Connection::~Connection(){
}

Connectivity::Connectivity(){
  nSelAtomsTot = 0;
  externalBonds = -1;
}

Connectivity::~Connectivity(){
  Clear();
  //cout << "Connectivity destructor" << endl;
}

int Connectivity::GetNumberOfAtoms(void) const{
  int i,j=0;
  //cout << "GetNumberOfAtoms nSelAtomsTot " << nSelAtomsTot << endl;
  for(i=0;i<nSelAtomsTot;i++){
    //cout << i << " " << flush;
    if(IsInTotalSelection(SelAtomsTot[i]->serNum-1)){
      j++;
    }
  }
  //cout << " j = " << j << endl;

  return j;
}

/* Returns a 1 if atom is C, 0 otherwise */
const std::vector<int> Connectivity::GetCAtomIndex(void) const {
  std::vector<int> indices;
  int i;

  const char* C = " C";
  for(i=0;i<nSelAtomsTot;i++){
    //if(IsInTotalSelection(SelAtomsTot[i]->serNum-1)){
    if(TotalSelection_new[SelAtomsTot[i]->serNum-1].isInTotalSelection!=0){
      //if(std::string(SelAtomsTot[i]->element)==std::string(" C")){
      if(strncmp(SelAtomsTot[i]->element,C,2)==0){
         indices.push_back(1);
      }else{
         indices.push_back(0);
      }
    }
  }

  return indices;
}

PPCAtom Connectivity::GetAtoms(void) const {
  PPCAtom atoms = new PCAtom[GetNumberOfAtoms()];
  int i,j=0;

  for(i=0;i<nSelAtomsTot;i++){
    if(SelAtomsTot[i]&&IsInTotalSelection(SelAtomsTot[i]->serNum-1)){
      //printf("Getting atom %d\n",j); fflush(stdout);
      atoms[j++] = SelAtomsTot[i];
    }
  }

  return atoms;
}

void Connectivity::Clear(void){
  nSelAtomsTot = 0;
  SelAtomsTot = 0;
  externalBonds = -1;
  TotalSelection_new.clear();
}

void Connectivity::AddBonds(PCMMANManager molhnd, int selhnd,
               const PPCAtom SelAtoms_in, const int nSelAtoms_in,
               int all_selHnd, int all_selHnd_spline){

  int i,j;
  SAtomBond *AtomBond;
  int nAtomBonds;

  std::vector<std::vector<int> > atom_connectivity;
  SelAtomsTot=0;
  externalBonds = 0;
  molhnd->GetAtomTable (SelAtomsTot,nSelAtomsTot);
  //cout << "AddBonds nSelAtomsTot" << nSelAtomsTot << endl;
  if(int(TotalSelection_new.size())!=nSelAtomsTot)
     TotalSelection_new = std::vector<Connection>(nSelAtomsTot);
	  
  /* Add the neccessary connections to the TotalSelection */
  for(i=0;i<nSelAtoms_in;i++){
    AtomBond=0;
    SelAtoms_in[i]->GetBonds ( AtomBond, nAtomBonds);
    Connection &conn = TotalSelection_new[SelAtoms_in[i]->serNum-1];
    if(!conn.isInTotalSelection){
      conn.isInTotalSelection = true;
    }
    //std::cout << "Connecting " << SelAtoms_in[i]->serNum-1 << " " << molhnd->AtomLabel_atom(SelAtoms_in[i]); std::cout.flush();
    for (j=0;j<nAtomBonds;j++) {
      if(AtomBond[j].atom->isInSelection(selhnd)){
        //std::cout << " to " << AtomBond[j].atom->serNum-1 << " " << molhnd->AtomLabel_atom(AtomBond[j].atom)<<  std::endl; std::cout.flush();
	conn.AddConnection(AtomBond[j].atom->serNum-1);
      } else {
         externalBonds++;
         if ( all_selHnd > 0 &&
		   AtomBond[j].atom->isInSelection(all_selHnd) ) {
            conn.AddExternalConnection(AtomBond[j].atom->serNum-1);
         }
         if ( all_selHnd_spline > 0 &&
		   AtomBond[j].atom->isInSelection(all_selHnd_spline) ) {
            conn.AddExternalSplineConnection(AtomBond[j].atom->serNum-1);
         }
      }
       
    }
  }

  //std::cout << "externalBonds " << externalBonds << std::endl;
  //printf("TotalSelection.size():%d\n",TotalSelection_new.size());

}

void Connectivity::AddContacts(PCMMANManager molhnd, int selhnd,
     const PPCAtom SelAtoms_in, const int nSelAtoms_in,
     const PSContact contacts_in, const int ncontacts_in){

  std::vector<std::vector<int> > atom_connectivity;

  SelAtomsTot=0;
  molhnd->GetAtomTable (SelAtomsTot,nSelAtomsTot);

  if(int(TotalSelection_new.size())!=nSelAtomsTot)
     TotalSelection_new = std::vector<Connection>(nSelAtomsTot);
	  
  for(int j=0;j<ncontacts_in;j++){
    Connection &conn = TotalSelection_new[SelAtoms_in[contacts_in[j].id1]->serNum-1];
    if(!conn.isInTotalSelection){
      conn.isInTotalSelection = true;
    }
    conn.AddConnection(SelAtoms_in[contacts_in[j].id2]->serNum-1);
  }

  //printf("Connectivity::AddContacts TotalSelection.size():%d\n",TotalSelection_new.size());

}

void Connectivity::AddTrace(PCMMANManager molhnd, const PPCAtom selAtoms,
              const int nSelAtoms, realtype cutoff ) {
  realtype cutoff2,dist;

  // Create a trace through a list of atoms - only allow connection between 
  // consecutive atoms which are within cutoff distance of each other
  SelAtomsTot=0;
  molhnd->GetAtomTable (SelAtomsTot,nSelAtomsTot);
  if(int(TotalSelection_new.size())!=nSelAtomsTot)
     TotalSelection_new = std::vector<Connection>(nSelAtomsTot);

  cutoff2 = cutoff*cutoff;

  // Loop over all atoms and tedt if within cutof distance of previous atom
  for (int j =1;j<nSelAtoms;j++) {
    dist = pow (( selAtoms[j-1]->x - selAtoms[j]->x ),2) +
           pow (( selAtoms[j-1]->y - selAtoms[j]->y ),2) +
           pow (( selAtoms[j-1]->z - selAtoms[j]->z ),2) ;
    
      // Add connection from j-1th to jth atom
    Connection &conn = TotalSelection_new[selAtoms[j-1]->serNum-1];
    if(!conn.isInTotalSelection){
      conn.isInTotalSelection = true;
    }
    if ( dist > 0.0001 && dist < cutoff2 ) {
      conn.AddConnection(selAtoms[j]->serNum-1);
      // Add connection from jth to j-1th atom
      Connection &conn2 = TotalSelection_new[selAtoms[j]->serNum-1];
      if(!conn2.isInTotalSelection){
        conn2.isInTotalSelection = true;
        //TotalSelection_new[selAtoms[j]->serNum-1] = conn;
      }
      //printf("Adding connection between %d and %d\n",j,j-1);
      conn2.AddConnection(selAtoms[j-1]->serNum-1);
    }
  }
  //printf("Connectivity::AddTrace TotalSelection.size():%d\n",TotalSelection_new.size());
}

std::vector<std::vector<int> > Connectivity::GetConnectivityLists(void) const{
  /* In here we have to map serNum to Index in our selection */

  std::vector<std::vector<int> > conn_lists;

  int i=0;
  std::vector<int> selection_map;

  selection_map.reserve(TotalSelection_new.size());
  conn_lists.reserve(TotalSelection_new.size());

  std::vector<Connection>::const_iterator k=TotalSelection_new.begin();
  while(k!=TotalSelection_new.end()){
    if(k->isInTotalSelection){
       //printf("Pushing back:%d\n",i);
       selection_map.push_back(i++);
    }else{
       //printf("Pushing back:%d\n",-1);
       selection_map.push_back(-1);
    }
    k++;
  }

  k=TotalSelection_new.begin();
  std::vector<int> connected_atoms;
  std::vector<int>::iterator mapiter=selection_map.begin();
  std::vector<int>::const_iterator atomiter;

  //std::cout << "TotalSelection_new[0].size " << TotalSelection_new[0].connected_atoms.size() << "\n";
  //std::cout << "TotalSelection_new[1].size " << TotalSelection_new[1].connected_atoms.size() << "\n";
  //std::cout << "TotalSelection_new[2].size " << TotalSelection_new[2].connected_atoms.size() << "\n";

  while(mapiter!=selection_map.end()){
    if(*mapiter!=-1){
      //printf("Got connection:%d\n",*mapiter);
      connected_atoms.clear();
      atomiter = k->connected_atoms.begin();
      while(atomiter!=k->connected_atoms.end()){
	//printf("Pushing atom connection:%d %d\n",*atomiter,selection_map[*atomiter]);      
        connected_atoms.push_back(selection_map[*atomiter]);
	atomiter++;
      }
      conn_lists.push_back(connected_atoms);
    }
    mapiter++; k++;
  }

  //printf("TotalSelection.size():%d\n",TotalSelection_new.size());
  //printf("conn_lists.size():%d\n",conn_lists.size());
  //printf("selection_map.size():%d\n",selection_map.size());

  //for(unsigned ii=0;ii<conn_lists.size();ii++){
     //std::cout << ii << ": ";
     //for(unsigned jj=0;jj<conn_lists[ii].size();jj++){
        //std::cout << conn_lists[ii][jj] << " ";
     //}
     //std::cout << "\n"; std::cout.flush();
  //}

  return conn_lists;
}

std::vector<std::vector<int> > Connectivity::GetExternalSplineConnectivityLists(void) const {
  /* In here we return  serNum-1 lists, ie index in total selection */

  std::vector<std::vector<int> > conn_lists;

  for(int i=0;i<nSelAtomsTot;i++){
    if(IsInTotalSelection(SelAtomsTot[i]->serNum-1)){
      conn_lists.push_back(TotalSelection_new[i].external_connected_spline_atoms);
    }
  }

  return conn_lists;
}

std::vector<std::vector<int> > Connectivity::GetExternalConnectivityLists(void) const {
  /* In here we return  serNum-1 lists, ie index in total selection */

  std::vector<std::vector<int> > conn_lists;

  for(int i=0;i<nSelAtomsTot;i++){
    if(IsInTotalSelection(SelAtomsTot[i]->serNum-1)){
      conn_lists.push_back(TotalSelection_new[i].external_connected_atoms);
    }
  }

  return conn_lists;
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
// Connectivity2 class to hold two vectors of atom pointers to store
// pairs of connected atoms

Connectivity2::Connectivity2(int data_mode_in,int crystal_axes_in,
                    int tagged_in) {

  data_mode = data_mode_in;
  if (data_mode == CONN_ATOM_POINT ||
      data_mode == CONN_POINT_ATOM || data_mode == CONN_ATOM_DISP )
    nAtomSets = 1;
  else if  (data_mode == CONN_ATOM_ATOM )
    nAtomSets = 2;
  else
    nAtomSets = 0;

  if (data_mode == CONN_ATOM_DISP || data_mode == CONN_POINT_DISP )
    nxyzSets = 1;
  else
    nxyzSets = 0;
  selection_on = false;
  crystal_axes = crystal_axes_in;
  tagged = tagged_in;
}

Connectivity2::~Connectivity2() {
  //cout << "Connectivity2 destructor" << endl;
}


void Connectivity2::Clear(void){

  connected.clear();
  if (nAtomSets>=1 ) {
    pAtom1.clear(); 
    if (nAtomSets>=2 ) pAtom2.clear(); 
  }
  if ( nxyzSets>= 1) XYZ1.clear();

  if (tagged>0) tags.clear();
 
}

void Connectivity2::AddConnection( PCAtom p_atom1, PCAtom p_atom2,
                    const std::string &label_in, int tag) {

  if (data_mode != CONN_ATOM_ATOM ) return; 
  
  pAtom1.push_back(p_atom1);
  pAtom2.push_back(p_atom2);
  Cartesian p1 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
  Cartesian p2 = Cartesian(p_atom2->x, p_atom2->y, p_atom2->z);
  connected.push_back(SimpleConnection(p1,p2,label_in));
  if (tagged>0 ) tags.push_back(tag);
  //cout << "Add connection " << p_atom1->GetModelNum() << " "
    //   << p_atom2->GetModelNum();
}
void Connectivity2::AddConnection( double xyz1[3] ,  double xyz2[3] ,
    const std::string &label_in, int tag ){

  if (data_mode == CONN_POINT_POINT) {
    Cartesian p1 = Cartesian(xyz1[0], xyz1[1], xyz1[2]);
    Cartesian p2 = Cartesian(xyz2[0], xyz2[1], xyz2[2]);
    connected.push_back(SimpleConnection(p1,p2,label_in));
    if (tagged>0 ) tags.push_back(tag);
  }
  else if  (data_mode == CONN_POINT_DISP) {
    Cartesian p1 = Cartesian(xyz1[0], xyz1[1], xyz1[2]);
    Cartesian p2 = Cartesian(xyz1[0]+xyz2[0],
                      xyz1[0]+xyz2[1], xyz1[2]+xyz1[2]);
    Cartesian pp =   Cartesian(xyz2[0], xyz2[1], xyz2[2]);
    connected.push_back(SimpleConnection(p1,p2,label_in));
    if (tagged>0 ) tags.push_back(tag);
    XYZ1.push_back(pp);
  }

}

void Connectivity2::AddConnection(  PCAtom p_atom1, double xyz[3] ,
    const std::string &label_in, int tag ){

  if (data_mode != CONN_ATOM_POINT && 
      data_mode != CONN_POINT_ATOM &&
      data_mode != CONN_ATOM_DISP  ) return; 
  Cartesian p1,p2,pp;
  
  switch (data_mode) {
  case CONN_ATOM_DISP:
    p1 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
    p2 = Cartesian(p_atom1->x+xyz[0],
                      p_atom1->y+xyz[1], p_atom1->z+xyz[2]);
    pp =   Cartesian(xyz[0], xyz[1], xyz[2]);
    XYZ1.push_back(pp);
    break;
  case CONN_POINT_ATOM:
     p1 = Cartesian(xyz[0], xyz[1], xyz[2]);
     p2 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
     break;
  default:
     p1 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
     p2 = Cartesian(xyz[0], xyz[1], xyz[2]);
     break;
  }
  pAtom1.push_back(p_atom1);
  connected.push_back(SimpleConnection(p1,p2,label_in));
  if (tagged>0 ) tags.push_back(tag);
  
}

void Connectivity2::InsertConnection( bool replace, unsigned int position,
                    PCAtom p_atom1, PCAtom p_atom2, 
                    const std::string &label_in, int tag) {

  if (data_mode != CONN_ATOM_ATOM ) return; 

  std::vector<PCAtom> new_atoms1;
  std::vector<PCAtom> new_atoms2;
  std::vector<SimpleConnection> new_conns;
  std::vector<int> new_tags;
  std::vector<bool> new_selected;

  for(unsigned int i=0;i<connected.size();i++){
    if(i == position){
      new_atoms1.push_back(p_atom1);
      new_atoms2.push_back(p_atom2);
      Cartesian p1 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
      Cartesian p2 = Cartesian(p_atom2->x, p_atom2->y, p_atom2->z);
      new_conns.push_back(SimpleConnection(p1,p2,label_in));
      if (tagged>0 ) new_tags.push_back(tag);
      if (selection_on) new_selected.push_back(selected[i]);
    }
    if ( i != position || !replace ) {
      new_atoms1.push_back(pAtom1[i]);
      new_atoms2.push_back(pAtom2[i]);
      new_conns.push_back(connected[i]);
      if (tagged>0) new_tags.push_back(tags[i]);
      if (selection_on) new_selected.push_back(selected[i]);
    }
  }
  pAtom1 = new_atoms1;
  pAtom2 = new_atoms2;
  connected = new_conns;
  if (tagged>0) tags = new_tags;
  if (selection_on) selected = new_selected; 
  
}


void Connectivity2::InsertConnection( bool replace, unsigned int position,
    double xyz1[3] ,  double xyz2[3] ,
    const std::string &label_in, int tag ){
  if (data_mode != CONN_POINT_POINT &&
      data_mode != CONN_POINT_DISP  ) return; 

  std::vector<Cartesian> new_xyz;
  std::vector<SimpleConnection> new_conns;
  std::vector<int> new_tags;
  std::vector<bool> new_selected;
  Cartesian p1,p2,pp;

  for(unsigned int i=0;i<connected.size();i++){
    if(i==position){
      if (data_mode == CONN_POINT_POINT) {
        p1 = Cartesian(xyz1[0], xyz1[1], xyz1[2]);
        p2 = Cartesian(xyz2[0], xyz2[1], xyz2[2]);
        new_conns.push_back(SimpleConnection(p1,p2,label_in));
        if (tagged>0 ) new_tags.push_back(tag);
      }
      else if  (data_mode == CONN_POINT_DISP) {
        p1 = Cartesian(xyz1[0], xyz1[1], xyz1[2]);
        p2 = Cartesian(xyz1[0]+xyz2[0],
                      xyz1[1]+xyz2[1], xyz1[2]+xyz2[2]);
        pp =   Cartesian(xyz2[0], xyz2[1], xyz2[2]);
        new_conns.push_back(SimpleConnection(p1,p2,label_in));
        if (tagged>0 ) new_tags.push_back(tag);
        new_xyz.push_back(pp);

      }
    }
    if ( i != position || !replace ) {
       if (nxyzSets>=1) new_xyz.push_back(XYZ1[i]);
       new_conns.push_back(connected[i]);
       if (tagged>0) new_tags.push_back(tags[i]);
       if (selection_on) new_selected.push_back(selected[i]);
    }
  }
  connected = new_conns;
  if (tagged>0) tags = new_tags;
  if (selection_on) selected = new_selected; 
  if (nxyzSets>=1) XYZ1 = new_xyz;
}

void Connectivity2::InsertConnection(  bool replace, unsigned int position,
     PCAtom p_atom1, double xyz[3] ,
    const std::string &label_in, int tag ){

  if (data_mode != CONN_ATOM_POINT && 
      data_mode != CONN_POINT_ATOM &&
      data_mode != CONN_ATOM_DISP  ) return; 
  Cartesian p1,p2,pp;

  std::vector<PCAtom> new_atoms1;
  std::vector<Cartesian> new_xyz;
  std::vector<SimpleConnection> new_conns;
  std::vector<int> new_tags;
  std::vector<bool> new_selected;

  for(unsigned int i=0;i<connected.size();i++){

    if ( i == position) {
      switch (data_mode) {
      case CONN_ATOM_DISP:
        p1 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
        p2 = Cartesian(p_atom1->x+xyz[0],
                      p_atom1->y+xyz[1], p_atom1->z+xyz[2]);
        pp =   Cartesian(xyz[0], xyz[1], xyz[2]);
        new_xyz.push_back(pp);
        break;
      case CONN_POINT_ATOM:
        p1 = Cartesian(xyz[0], xyz[1], xyz[2]);
        p2 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
        break;
      default:
        p1 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
        p2 = Cartesian(xyz[0], xyz[1], xyz[2]);
        break;
      }
      new_atoms1.push_back(p_atom1);
      connected.push_back(SimpleConnection(p1,p2,label_in));
      if (tagged>0 ) tags.push_back(tag);
    }
    if ( i != position || !replace ) {
      new_atoms1.push_back(pAtom1[i]);
       if (nxyzSets>=1) new_xyz.push_back(XYZ1[i]);
       new_conns.push_back(connected[i]);
       if (tagged>0) new_tags.push_back(tags[i]);
       if (selection_on) new_selected.push_back(selected[i]);
    }
  }
  pAtom1 = new_atoms1;
  connected = new_conns;
  if (nxyzSets>=1) XYZ1=new_xyz;
  if (tagged>0) tags = new_tags;
  if (selection_on) selected = new_selected; 
  
}


void Connectivity2::AddUniqueConnection( PCAtom p_atom1, PCAtom p_atom2,
                    const std::string &label_in,int tag){
  std::vector<unsigned int> hits;
  if (data_mode != CONN_ATOM_ATOM ) return; 

  hits = FindConnections( p_atom1,p_atom2, true );
  //cout << "AddUniqueConnection " << hits.size() << endl;
  if (hits.size() > 0) return;  

  Cartesian p1 = Cartesian(p_atom1->x, p_atom1->y, p_atom1->z);
  Cartesian p2 = Cartesian(p_atom2->x, p_atom2->y, p_atom2->z);
  
  pAtom1.push_back(p_atom1);
  pAtom2.push_back(p_atom2);
  connected.push_back(SimpleConnection(p1,p2,label_in));
  //cout << "Add connection " << p_atom1->GetModelNum() << " "
    //   << p_atom2->GetModelNum();
  if (tagged>0 ) tags.push_back(tag);

}

void Connectivity2::UpdateCoordinates(bool label_dist) {
  /*
  Update the connected list after atoms have moved
  */
  if (data_mode == CONN_POINT_POINT ||
      data_mode == CONN_POINT_DISP ) return;

  Cartesian p1,p2,p3;
  std::vector<SimpleConnection> new_conns;

  std::vector<SimpleConnection>::iterator old_conn=connected.begin();
  std::vector<PCAtom>::iterator at1=pAtom1.begin();
  std::vector<PCAtom>::iterator at2;
  std::vector<Cartesian>::iterator  xyz;

  switch (data_mode) {
  case CONN_ATOM_ATOM:
    at2=pAtom2.begin();
    while(at1!=pAtom1.end()){
      p1 = Cartesian((*at1)->x, (*at1)->y, (*at1)->z);
      p2 = Cartesian((*at2)->x, (*at2)->y, (*at2)->z);
      if (label_dist) {
        p3 = p1-p2;
        new_conns.push_back(SimpleConnection(p1,p2,FloatToString(p3.length(),"%.1f")));
      } else {
        new_conns.push_back(SimpleConnection(p1,p2,old_conn->label));
      }
      at1++;
      at2++;
      old_conn++;
    }
    break;
  case CONN_ATOM_POINT:
    while(at1!=pAtom1.end()){
      p1 = Cartesian((*at1)->x, (*at1)->y, (*at1)->z);
      new_conns.push_back(SimpleConnection(p1, old_conn->second,
                              old_conn->label));
      at1++;
      old_conn++;
    }
    break;
  case CONN_POINT_ATOM:
    while(at1!=pAtom1.end()){
      p2 = Cartesian((*at1)->x, (*at1)->y, (*at1)->z);
      new_conns.push_back( SimpleConnection(old_conn->first, 
                           p2,old_conn->label));
      at1++;
      old_conn++;
    }
    break;
  case CONN_ATOM_DISP:
    xyz = XYZ1.begin();
    while(at1!=pAtom1.end()){
      p1 = Cartesian((*at1)->x, (*at1)->y, (*at1)->z);
      p2 = Cartesian((*at1)->x+xyz->get_x(),(*at1)->y+xyz->get_y(), 
                           (*at1)->z+xyz->get_z());
      new_conns.push_back( SimpleConnection(p1,p2,old_conn->label));
      at1++;
      xyz++;
      old_conn++;
    }
  }
  connected = new_conns;
}
 
void Connectivity2::DeleteConnections(std::vector<unsigned int> indices, bool all){

  std::vector<PCAtom> new_atoms1;
  std::vector<PCAtom> new_atoms2;
  std::vector<Cartesian> new_xyz;
  std::vector<SimpleConnection> new_conns;
  std::vector<int> new_tags;
  std::vector<bool> new_selected;

  if (!all) {
    sort(indices.begin(),indices.end());

    for(unsigned int i=0;i<connected.size();i++){
      if(!binary_search(indices.begin(),indices.end(),i)){
        if (nAtomSets>=1) {
          new_atoms1.push_back(pAtom1[i]);
          if (nAtomSets>=2) new_atoms2.push_back(pAtom2[i]);
        }
        if (nxyzSets>=1) new_xyz.push_back(XYZ1[i]);
        new_conns.push_back(connected[i]);
        if (tagged>0) new_tags.push_back(tags[i]);
        if (selection_on) new_selected.push_back(selected[i]);
      }
    }
  }
  if (nAtomSets>=1) {
    pAtom1 = new_atoms1;
    if (nAtomSets>=2) pAtom2 = new_atoms2;
  }
  connected = new_conns;
  if (tagged>0) tags = new_tags;
  if (selection_on) selected = new_selected; 
  if (nxyzSets>=1) XYZ1=new_xyz;

}


void Connectivity2::DeleteConnection( int index ) {
  std::vector<unsigned int> indices;
  indices.push_back((unsigned int)index);
  DeleteConnections(indices);
}

void Connectivity2::DeleteConnections( ) {
  std::vector<unsigned int> indices;
  DeleteConnections(indices,true);
}

int Connectivity2::DeleteTaggedConnections(int first, int last ) {
  std::vector<unsigned int> indices;
  if (tagged==0) return 0;
  std::vector<int>::const_iterator i = tags.begin();
  while(i!=tags.end()){
    if (*i >= first && (last<0 || *i<=last)) indices.push_back((unsigned int)*i);
    i++;
  }
  int ndel = indices.size();
  if (ndel>0) DeleteConnections(indices,false);
  return ndel;
}

void Connectivity2::RemoveConnection( PCAtom p_atom1, int position ) {
  std::vector<unsigned int> indices = FindConnections(p_atom1,position);
  DeleteConnections(indices);
}

void Connectivity2::RemoveConnection( PCAtom p_atom1, PCAtom p_atom2 ) {
  std::vector<unsigned int> indices = FindConnections(p_atom1,p_atom2);
  //cout << "RemoveConnection indices " << indices << endl;
  DeleteConnections(indices);
}
  
std::vector<unsigned int> Connectivity2::FindConnections( PCAtom p_atom1,
                  PCAtom p_atom2, bool switchpos ) {
  // return index of all connections between the two specified atoms
  // If switch true then atoms could be inverted order
  std::vector <unsigned int> indices;
  if(nAtomSets<2) return indices;
  std::vector<PCAtom>::iterator i=pAtom1.begin();
  std::vector<PCAtom>::iterator j=pAtom2.begin();

  int n = 0;
  while(i!=pAtom1.end()){
    if (*i == p_atom1 && *j ==  p_atom2) indices.push_back(n);
    i++;
    j++;
    n++;
  }
  n = 0;
  if (switchpos ) {
    i=pAtom1.begin();
    j=pAtom2.begin();
    while(i!=pAtom1.end()){
      if (*i== p_atom2 && *j == p_atom1) indices.push_back(n);
      i++;
      j++;
      n++;
    }
  }
  return indices;
}

int Connectivity2::FindNofConnections( PCAtom p_atom1 ,int position) {
  std::vector <unsigned int> indices = FindConnections(p_atom1 ,position);
  return indices.size();
}

std::vector <unsigned int> Connectivity2::FindConnections( PCAtom p_atom1 ,
         int position) {
  
  //cout << "FindConnections " << p_atom1 << " " << position << endl;
  int n=0;
  std::vector <unsigned int> indices;
  if (nAtomSets<1) return indices;
  if ( position == 1 || nAtomSets<=1) {
    std::vector<PCAtom>::iterator i=pAtom1.begin();
    while(i!=pAtom1.end()){
      if (*i == p_atom1 ) indices.push_back(n);
      i++;
      n++;
    }
  } 
  else if ( position == 2 ) {
    std::vector<PCAtom>::iterator i=pAtom2.begin();
    while(i!=pAtom2.end()){
      if ( *i == p_atom1) indices.push_back(n);
      i++;
      n++;
    }
  } 
  else {
    std::vector<PCAtom>::iterator i=pAtom1.begin();
    std::vector<PCAtom>::iterator j=pAtom2.begin();
    while(i!=pAtom1.end()){
      if (*i == p_atom1 || *j ==  p_atom1) indices.push_back(n);
      i++;
      j++;
      n++;
    }
  }
  //cout << "indices " << indices.size() << endl;
  return indices;
}

double* Connectivity2::Extent() {
  /*
  Find the max/min of x/y/z of all vector ends
  */
  double coord;
  double *com;
  com = new double[6];
  com[0] = 9999999.9;
  com[1] = 9999999.9;
  com[2] = 9999999.9;
  com[3] = -9999999.9;
  com[4] = -9999999.9;
  com[5] = -9999999.9;
  std::vector<SimpleConnection>::const_iterator i = connected.begin();
  while(i!=connected.end()){
    coord = i->first.get_x();
    if (coord < com[0] ) com[0] = coord ;
    if (coord > com[3] ) com[3] = coord ;
    coord = i->second.get_x();
    if (coord < com[0] ) com[0] = coord ;
    if (coord > com[3] ) com[3] = coord ;
    coord = i->first.get_y();
    if (coord < com[1] ) com[1] = coord ;
    if (coord > com[4] ) com[4] = coord ;
    coord = i->second.get_y();
    if (coord < com[1] ) com[1] = coord ;
    if (coord > com[4] ) com[4] = coord ;
    coord = i->first.get_z();
    if (coord < com[2] ) com[2] = coord ;
    if (coord > com[5] ) com[5] = coord ;
    coord = i->second.get_z();
    if (coord < com[2] ) com[2] = coord ;
    if (coord > com[5] ) com[5] = coord ;
    i++;
  }
  return com;
}

int Connectivity2::SelectVectors (int nSelTags, int *selTags ) {
  /*
  Reset selected list.  For each vector seelcted is true
  if its tag value is one of those in the list selTags
  */
  int j;
  bool sel;
  if ( tagged==0) return 1;
  if ( selection_on ) selected.clear();
  if (nSelTags < 0 ) {
    selection_on = false;
     return 0;
  }
  std::vector<int>::const_iterator i = tags.begin();
  while(i!=tags.end()){
    sel = false;
    for (j=0;j< nSelTags;j++) {
      if (*i == selTags[j]) {
        sel = true;
        break;
      }
    }
    selected.push_back(sel);
    i++;
  }
  return 0;
  
}

Cartesian Connectivity2::GetCoordinate ( int i, int j) {
  if ( i < 0 || i>=int(connected.size()))
    return Cartesian(0,0,0);

  if (j==1)
    return connected[i].first;
  else if ( j == 2)
    return connected[i].second;
  else if (j == 3 && i < int(XYZ1.size()) )
    return XYZ1[i];

  return Cartesian(0,0,0); // Don't get here, just keeps compiler happy.

}


PCAtom Connectivity2::GetAtom ( int i, int j) {
  if ( i < 0 ) return NULL;
  if (j==1 || data_mode == CONN_POINT_ATOM) {
   if (i < int(pAtom1.size())) return pAtom1[i];
  }  else if (j==2) {
   if ( i < int(pAtom2.size())) return pAtom2[i];
  }
  return NULL;
}
std::string Connectivity2::GetAtomID ( int i, int j) {
  char tmp[200];
  if ( i < 0 ) return NULL;
  if (j==1 || data_mode == CONN_POINT_ATOM) {
   if (i < int(pAtom1.size())) return pAtom1[i]->GetAtomID(tmp);
  }  else if (j==2) {
   if ( i < int(pAtom2.size())) return pAtom2[i]->GetAtomID(tmp);
  }
  return "";
}

int Connectivity2::GetTag ( int iV ) {
  if ( iV >= 0 && iV < int(tags.size()) ) return tags[iV];
  return -1;
}

int Connectivity2::SetTag ( int iV, int i ) {
  if ( iV >= 0 && iV < int(tags.size()) ) {
    tags[iV]=i;
    return 0;
  } else {
    return 1;
  }
}

std::string Connectivity2::Print ( const std::string &tag1, const std::string &tag2 ) {
  PCMMUTManager M1,M2;

  std::ostringstream output;
  output.setf(ios::fixed);
  output.setf(ios::showpoint);
  output.setf(ios::left,ios::adjustfield);
  output.precision(1);
  std::vector<PCAtom>::iterator i;
  std::vector<PCAtom>::iterator j;
  
  if ( nAtomSets>=1 && pAtom1.size() < 1)  return output.str();

  std::vector<SimpleConnection>::const_iterator cc = connected.begin();
  
  if ( nAtomSets>=1) {
    i=pAtom1.begin();
    M1 = (PCMMUTManager)(PCMMDBFile)(*i)->GetCoordHierarchy();
    if ( nAtomSets>=2) {
      j=pAtom2.begin();
      M2 = (PCMMUTManager)(PCMMDBFile)(*j)->GetCoordHierarchy();
    }
  }
  if ( tag1.length()>=1 &&  (nAtomSets<=1 || tag1.length()>=1 )) {
    while(i<pAtom1.end()){
      //output << "<dataobj=" << tag1 << ">";
      output << tag1 << ": ";
      output << M1->AtomLabel_atom(*i) << " ";
      //output << "</dataobj>";
      if  ( nAtomSets>=2) {
        //output << "<dataobj=" << tag2 << ">";
        output << tag2 << ": ";
        output << M2->AtomLabel_atom(*j) << " ";
	// output << "</dataobj>";
        j++;
      }
      output << cc->label << endl;
      i++;
      cc++;
    }
  } else {
    while(i<pAtom1.end()){
      output <<  std::setw(20) << M1->AtomLabel_atom(*i) << " ";
      if  ( nAtomSets>=2) { 
        output << std::setw(20) <<  M2->AtomLabel_atom(*j) << " ";
        j++;
      }
      output << std::setw(10) << cc->label << endl;
      i++;
      cc++;
    }
  }
  return output.str();
 
}

int Connectivity2::AddContacts(PCMMANManager molHnd1,int selHnd1, PCMMANManager molHnd2,int selHnd2_in,realtype  dist1, realtype  dist2, int  seqDist , int inter_model, int closest_bonding, int handle_hbond) {
  // closest_bonding -- min 'bonding' number allowed from TestBonding 
  // -- can be 0, 2=bonded, 3=1-3bonded 1-4 = 1-4 bonded
  // handle_hbond - flag treatment of HBond -- not implemented yet
  PPCAtom selAtoms,selAtoms2;
  int selHnd2,nat1,nat2;
  PSContact  contacts = NULL;
  int ncontacts=0, nb_contacts=0;
  int bond,modno;

  molHnd1->GetSelIndex(selHnd1,selAtoms,nat1);
  if (nat1<=0) return 0;
  if (selHnd2_in<=0) {
    if (inter_model>0) {
      modno = 0;
    } else {
      modno = selAtoms[0]->GetModel()->GetSerNum();
    }
    selHnd2=molHnd2->NewSelection();
    molHnd2->SelectAtoms (selHnd2,modno,"*", ANY_RES,"*", ANY_RES,"*","*","*","*","*",SKEY_NEW);
  } else {
    selHnd2 = selHnd2_in;
  }
  molHnd2->GetSelIndex(selHnd2,selAtoms2,nat2);
  //std::cout << "nat1,nat2 " << nat1 << " " << nat2 <<endl;
  if (nat1>0 && nat2>0) {
    mat44 * TMatrix=0;
    molHnd1->SeekContacts( selAtoms,nat1,selAtoms2,nat2,dist1,dist2,seqDist,
			 contacts,ncontacts,0,TMatrix,0,0);
    
    //std::cout << "ncontacts " << ncontacts <<endl;
    for (int n=0;n<ncontacts;n++) {
      bond =  molHnd1->TestBonding(selAtoms[contacts[n].id1],
				   selAtoms2[contacts[n].id2]);
      // Exclude bonded, 1-3 contacts and 1-4 contacts
      if ( (bond == 0 || bond >= closest_bonding) &&
           molHnd1->doAltLocMatch(selAtoms[contacts[n].id1],
                                  selAtoms2[contacts[n].id2]) ) {
        AddUniqueConnection(selAtoms[contacts[n].id1],
                            selAtoms2[contacts[n].id2],
                           FloatToString(contacts[n].dist,"%.1f"));
        // Sum the non-bonded contacts
        nb_contacts++;
      }

    }
  }

  
  if (selHnd2_in<=0) molHnd2->DeleteSelection(selHnd2);
  return nb_contacts;
}

//----------------------------------------------------------------------------
int Connectivity2::AddRangeConnections( PCResidue res1,PCResidue res2, 
			     PCResidue mres1, PCResidue mres2,
			     const std::vector<std::string>&  mainchain_name,
                                       int tag ) {
//----------------------------------------------------------------------------
  // Set 1 or 2 - the set to which res1 and res2 belong
  // Match the range of residue res1->res2 (not necessarilly input
  // in correct sequence order) to mres1 (equivalent of res1)
  // selection_mode: CA/main/all etc
  PCResidue r1,r2;
  PCAtom pat1,pat2;
  int nconn = 0;
  int mainchain_max= mainchain_name.size();

  if (res1->GetChain() != res2->GetChain()) return -1;
  int inc  = 1;
  int irdif = res2->GetResidueNo() - res1->GetResidueNo();
  if (irdif<0) inc = -1;
  irdif =irdif + inc;


  int ir = 0;
  while (ir!=irdif) {
    r1 = res1->GetChain()->GetResidue(res1->GetResidueNo()+ir);
    r2 = mres1->GetChain()->GetResidue(mres1->GetResidueNo()+ir);
    if (r1 !=NULL && r2 != NULL) {
    //cout << "r1,r2" << r1 << " " << r2 <<endl;
      if (mainchain_max>0) {
        for (int ia=0;ia<mainchain_max;ia++) {         
          pat1 = r1->GetAtom(mainchain_name[ia].c_str(),"*","");
          pat2 = r2->GetAtom(mainchain_name[ia].c_str(),"*","");
          if (pat1!=NULL && pat2!=NULL) {
            nconn++;
            AddConnection(pat1,pat2,"",tag);
          }
        }
      } else {
        for (int ia=0;ia<r1->GetNumberOfAtoms();ia++) {
          pat2 = r2->GetAtom(r1->GetAtom(ia)->name,r1->GetAtom(ia)->element,r1->GetAtom(ia)->altLoc);
          //cout << "pat" << pat << endl;
          if (pat2!=NULL) {
            nconn++;
            AddConnection(r1->GetAtom(ia),pat2,"",tag);
          }
	}
      }
    }
    ir = ir + inc;
  }
  //cout << "nconn " << nconn <<endl;
  return nconn;
}

//----------------------------------------------------------------------------
int Connectivity2::AddRangeWithSameIdConnections( PCMMANManager molHnd1,int selHnd1, 
                                                  PCMMANManager molHnd2,int selHnd2, 
			     const std::vector<std::string>&  mainchain_name, int tag ) {
//----------------------------------------------------------------------------
  // Match the range of residue res1->res2 (not necessarilly input
  // in correct sequence order) to mres1 (equivalent of res1)
  // for the residues with the same sequence id 
  // selection_mode: CA/main/all etc
  int nconn = 0;
  int mainchain_max= mainchain_name.size();
  PPCResidue selRes1, selRes2;
  PCAtom pat1, pat2;
  int nr1,nr2;

  molHnd1->GetSelIndex(selHnd1,selRes1,nr1);
  molHnd2->GetSelIndex(selHnd2,selRes2,nr2);
  //cout << "AddRangeWithSameIdConnections " << nr1 << " " << nr2 << " " << mainchain_max<<endl;

  int targetHnd = molHnd2->NewSelection();
 
  for (int ir=0; ir< nr1; ir++) {
    
    molHnd2->Select ( targetHnd,STYPE_RESIDUE,selHnd2,SKEY_NEW);
    molHnd2->Select ( targetHnd,STYPE_RESIDUE,0,"*",selRes1[ir]->seqNum,selRes1[ir]->insCode,
                      selRes1[ir]->seqNum,selRes1[ir]->insCode,"*","*","*","*",SKEY_AND);
    molHnd2->GetSelIndex(targetHnd,selRes2,nr2);
    //cout << "AddRangeWithSameIdConnections " << ir << " " << selRes1[ir]->seqNum << selRes1[ir]->insCode << " " << nr2;
    //if (nr2>0) cout << " " << selRes2[0]->seqNum;
    //cout << endl;
    if (nr2>0) {
      if (mainchain_max>0) {
        for (int ia=0;ia<mainchain_max;ia++) {         
          pat1 = selRes1[ir]->GetAtom(mainchain_name[ia].c_str(),"*","");
          pat2 = selRes2[0]->GetAtom(mainchain_name[ia].c_str(),"*","");
          if (pat1!=NULL && pat2!=NULL) {
            nconn++;
            AddConnection(pat1,pat2,"",tag);
          }
        }
      } else {
        for (int ia=0;ia<selRes1[ir]->GetNumberOfAtoms();ia++) {
          pat1 = selRes1[ir]->GetAtom(ia);
          pat2 = selRes2[0]->GetAtom(pat1->name, pat1->element,pat1->altLoc);
          //cout << "pat" << pat << endl;
          if (pat2!=NULL) {
            nconn++;
            AddConnection(pat1,pat2,"",tag);
          }
	}
      }
    }
  }
  //cout << "nconn " << nconn <<endl;
  return nconn;
}

double Connectivity2::GetRMSD() {

  double rmsd,dist,distsum=0.0;
  int nat = 0;

  std::vector<PCAtom>::iterator at1=pAtom1.begin();
  std::vector<PCAtom>::iterator at2=pAtom2.begin();
  std::vector<Cartesian>::iterator  xyz;

  while(at1!=pAtom1.end()){
    //if (strcmp((*at1)->name,selected_atom)==0) {
      dist =   ((*at1)->x-(*at2)->x) *  ((*at1)->x-(*at2)->x) +
	     ((*at1)->y-(*at2)->y) *  ((*at1)->y-(*at2)->y) +
         	((*at1)->z-(*at2)->z) *  ((*at1)->z-(*at2)->z) ;
      //cout << "GetRMSD dist " << dist <<endl;
      at1++;
      at2++;
      distsum += dist;
      nat++;
      // }
  
  }
  rmsd = pow((distsum/nat),0.5);
  //cout << "Matching " << nat << " atoms with RMSD " << rmsd <<endl;
  return rmsd;
  
}

int Connectivity2::AddContactFromSelHandle( PCMMANManager molHnd1,int selHnd1, 
                                            PCMMANManager molHnd2,int selHnd2, int tag) {

  int na1,na2;
  PPCAtom pa1 = NULL;
  PPCAtom pa2 = NULL;
  double dist;
  molHnd1->GetSelIndex ( selHnd1,pa1,na1 );
  molHnd2->GetSelIndex ( selHnd2,pa2,na2 );
  //cout << "AddFromSelHandle na1,na2 " << na1 << " " << na2 << endl;
  if ( na1 > na2 ) na1 = na2;
  for (int i=0;i<na1;i++ ) {
    dist = sqrt ( ( pa1[i]->x-pa2[i]->x)* ( pa1[i]->x-pa2[i]->x) +
                  ( pa1[i]->y-pa2[i]->y)* ( pa1[i]->y-pa2[i]->y) +
                  ( pa1[i]->z-pa2[i]->z)* ( pa1[i]->z-pa2[i]->z)  ) ;
    AddConnection( pa1[i],pa2[i],FloatToString(dist,"%.1f"),tag);
  }
  return na1;
}

//------------------------------------------------------------------------
int Connectivity2::AddCloseAtoms( PCMMANManager molHnd1,int selHnd1_in, 
                                  PCMMANManager molHnd2,int selHnd2 , 
				  realtype central_cutoff, realtype cutoff , char *central, int tag ) {
//------------------------------------------------------------------------
  int selHnd1;
  PPCAtom SelAtoms1 = NULL;
  int NAtoms1;
  int ResHnd1;
  PPCResidue SelRes1 = NULL;
  int NRes1;
  int nConn=0;
  double dist;

  // Max possible number of superposed atoms
  selHnd1 =  molHnd1->NewSelection();
  if (selHnd1_in>0)
    molHnd1->Select(selHnd1,STYPE_ATOM,selHnd1_in,SKEY_NEW);
  else
    molHnd1->Select(selHnd1,STYPE_ATOM,0,"*",ANY_RES,
	       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);

  molHnd1->GetSelIndex(selHnd1,SelAtoms1,NAtoms1);
  if (NAtoms1<=0) return -1;

  
  //Get list of target residues
  ResHnd1 = molHnd1->NewSelection();
  molHnd1->Select(ResHnd1,STYPE_RESIDUE,selHnd1,SKEY_NEW);
  molHnd1->GetSelIndex(ResHnd1,SelRes1,NRes1);
  //cout << "fxNAtoms,fxNRes " << fxNAtoms << " " << fxNRes << endl;
  if (NRes1<=0) {
    molHnd1->DeleteSelection(ResHnd1);
    return -1;
  }
 

  //Get the 'central' atoms in the second model
  int centralHnd = molHnd2->NewSelection();
  PPCAtom centralAtoms;
  int centralNAtoms;
  if (selHnd2>0) {
    molHnd2->Select(centralHnd,STYPE_ATOM,selHnd2,SKEY_NEW);
    molHnd2->Select(centralHnd,STYPE_ATOM,0,"*",ANY_RES,
       "*",ANY_RES, "*","*","CA","C","*",SKEY_AND);
  } else
    molHnd2->Select(centralHnd,STYPE_ATOM,0,"*",ANY_RES,
       "*",ANY_RES, "*","*","CA","C","*",SKEY_NEW);
  molHnd2->GetSelIndex(centralHnd,centralAtoms,centralNAtoms);
  //cout << "centralNAtoms " << centralNAtoms << endl;

  // Loop over residues in fixed model

  PCAtom pCentralAtom;
  PCResidue closeRes;
  PCAtom fxAtom,closeAtom;
  int ncontacts,best_match;
  double score,best_score;

  for (int ir=0; ir<NRes1; ir++) {
    // Find the 'central' atoms in moving object that
    // are close to the central atom of fixed residue
    // -- there should just be one!

    pCentralAtom = SelRes1[ir]->GetAtom("CA","*","*");
    if (pCentralAtom) {
      PSContact contact= NULL;
      ncontacts = 0;
      molHnd2->SeekContacts (pCentralAtom,centralAtoms,centralNAtoms,
		    0.0, central_cutoff,0,contact,ncontacts,-1);
      best_match=-1;
      best_score=central_cutoff;
      if (ncontacts ==1 ) {
        best_match=0;
      } else if (ncontacts>1) {
        //cout <<  pCentralAtom->residue->seqNum << "ncontacts " << ncontacts << endl;
      
        for (int nc=0;nc<ncontacts;nc++) {
          score=contact[nc].dist -
           molHnd2->DeltaResidueOrientation(centralAtoms[contact[nc].id2]->residue,SelRes1[ir]);
          if(score<best_score && score <central_cutoff-1.5 ) {
            best_score=score;
            best_match=nc;
          }
        }
        //cout << "best " << best_match << " " << best_score << endl;
      }
      

      // -- there should be just one close residue!
      if (best_match >=0 ) {
        closeRes = centralAtoms[contact[best_match].id2]->residue;
        // Loop over atoms in the fixed residue
        // Is the atom in the original selection 
        for (int ia=0;ia<SelRes1[ir]->GetNumberOfAtoms();ia++) {
          fxAtom = SelRes1[ir]->GetAtom(ia);
          if ( fxAtom && fxAtom->isInSelection(selHnd1) ) {
	    // Is there an atom with same name in the close moving residue
	     closeAtom= closeRes->GetAtom(fxAtom->name,fxAtom->element,"*");
             if ( closeAtom && molHnd2->BondLength(fxAtom,closeAtom)<=cutoff) {
               dist = sqrt ( ( fxAtom->x-closeAtom->x)* ( fxAtom->x-closeAtom->x) +
                  ( fxAtom->y-closeAtom->y)* ( fxAtom->y-closeAtom->y) +
                  ( fxAtom->z-closeAtom->z)* ( fxAtom->z-closeAtom->z)  ) ;
               AddConnection(fxAtom,closeAtom,FloatToString(dist,"%.1f"),tag);
               nConn++;
             }
          }
        }
      }
      if (contact) delete contact;
    }
  }

  // Cleanup
  molHnd1->DeleteSelection(ResHnd1);
  molHnd2->DeleteSelection(centralHnd);

  return nConn;

}


//----------------------------------------------------------------------------
int Connectivity2::AddCloseRangeConnections(int set,PCResidue res1,
                    PCResidue res2,PCMMANManager M2, 
                    double central_cutoff, double cutoff, 
		    const std::vector<std::string>& mainchain_name,
                    const std::string &centralAtom, int tag ) {
//----------------------------------------------------------------------------
  // Set 1 or 2 - the set to which res1 and res2 belong
  // Match the range of residue res1->res2 (not necessarilly input
  // in correct sequence order) to 
  // M2 - the closest residues in this model
  if (res1->GetChain() != res2->GetChain()) return -1;
  int irdif = res2->GetResidueNo() - res1->GetResidueNo();
  //cout << "AddRangeConnections irdif " << irdif << endl;
  int nconn = 0;
  //Get the 'central' atoms in the 'unknown' model (without range selected)
  int centralHnd = M2->NewSelection();
  PPCAtom centralAtoms;
  PCAtom pat1,pat2;
  int centralNAtoms;
  M2->Select(centralHnd,STYPE_ATOM,0,"*",ANY_RES,
       "*",ANY_RES, "*","*",(char*)(centralAtom.c_str()),"*","*",SKEY_NEW);
  M2->GetSelIndex(centralHnd,centralAtoms,centralNAtoms);
  //cout << "AddCloseRangeConnections centralNAtoms " << centralNAtoms << endl;
  int mainchain_max= mainchain_name.size();


  // Loop over residues in known model
  PCChain pChn = res1->GetChain();
  PCAtom pCentralAtom;
  PCResidue closeRes;
  PCAtom fxAtom,closeAtom;
  int ifres,ilres,ncontacts,best_match;
  double score,best_score;

  if (irdif>0) {
    ifres = res1->GetResidueNo();
    ilres = res2->GetResidueNo();
  } else {
    ifres = res2->GetResidueNo();
    ilres = res1->GetResidueNo();
  }  
  for (int ir=ifres; ir<=ilres; ir++) {
    // Find the 'central' atoms in moving object that
    // are close to the central atom of fixed residue
    // -- there should just be one!

    pCentralAtom = pChn->GetResidue(ir)->GetAtom(centralAtom.c_str(),"*","*");
    if (pCentralAtom) {
      PSContact contact= NULL;
      ncontacts = 0;
      M2->SeekContacts (pCentralAtom,centralAtoms,centralNAtoms,
		    0.0, central_cutoff,0,contact,ncontacts,-1);
      best_match=-1;
      best_score=central_cutoff;
      if (ncontacts ==1 ) {
        best_match=0;
      } else if (ncontacts>1) {
        cout <<  pCentralAtom->residue->seqNum << "ncontacts " << ncontacts << endl;
      
        for (int nc=0;nc<ncontacts;nc++) {
          score=contact[nc].dist -
           M2->DeltaResidueOrientation(centralAtoms[contact[nc].id2]->residue,pChn->GetResidue(ir));
          if(score<best_score && score <central_cutoff-1.5 ) {
            best_score=score;
            best_match=nc;
          }
        }
          //cout << "best " << best_match << " " << best_score << endl;
      }
      

      // -- there should be just one close residue!
      if (best_match >=0 ) {
        closeRes = centralAtoms[contact[best_match].id2]->residue;
        if (mainchain_max>0) {
          for (int ia=0;ia<mainchain_max;ia++) {         
            pat1 = pChn->GetResidue(ir)->GetAtom(mainchain_name[ia].c_str(),"*","");
            pat2 = closeRes->GetAtom(mainchain_name[ia].c_str(),"*","");
            if (pat1!=NULL && pat2!=NULL && M2->BondLength(pat1,pat2)<=cutoff) {
              nconn++;
              if ( set == 1) {
                AddConnection(pat1,pat2,"",tag);
              } else {
                AddConnection(pat2,pat1,"",tag);
              }
            }
          }
        } else {
          // Loop over atoms in the 'known' residue
          // Is the atom in the original selection 
          for (int ia=0;ia< pChn->GetResidue(ir)->GetNumberOfAtoms();ia++) {
            fxAtom = pChn->GetResidue(ir)->GetAtom(ia);
	    // Is there an atom with same name in the close moving residue
	     closeAtom= closeRes->GetAtom(fxAtom->name,fxAtom->element,"*");
             if ( closeAtom && M2->BondLength(fxAtom,closeAtom)<=cutoff) {
               nconn++;
               if ( set == 1) {
                 AddConnection(fxAtom,closeAtom,"",tag);
               } else {
                 AddConnection(closeAtom,fxAtom,"",tag); 
               }
	     }
          }
        }
      }
      if (contact) delete contact;
    }
  }
  // Cleanup
  M2->DeleteSelection(centralHnd);
  //cout << "nconn " << nconn <<endl;
  return nconn;
}

int Connectivity2::GetSelection ( int mode) {
  PCMMUTManager M1;
  std::vector<PCAtom>::iterator i;

  if ( pAtom1.size() < 1)  return -1;
 
  if ( mode==1) {
    i=pAtom1.begin();
    M1 = (PCMMUTManager)(PCMMDBFile)(*i)->GetCoordHierarchy();
  } else if ( mode==2) {
    i=pAtom2.begin();
    M1 = (PCMMUTManager)(PCMMDBFile)(*i)->GetCoordHierarchy();
  } else {
    return -1;
  }
 
  
  int selHnd =  M1->NewSelection();

  while((mode==1&&i<pAtom1.end())||(mode==2&&i<pAtom2.end())) { 
    //cout << " Connectivity2::GetSelection " << (*i)->GetSeqNum()<< " "<< (*i)->GetInsCode() << " " << (*i)->name << " " << (*i)->altLoc << endl;
    M1->Select(selHnd,STYPE_ATOM,(*i)->GetModelNum(),(*i)->GetChainID(),(*i)->GetSeqNum(),(*i)->GetInsCode(),(*i)->GetSeqNum(),(*i)->GetInsCode(),"*",(*i)->name,"*",(*i)->altLoc,SKEY_OR);
    i++;
  }

  return selHnd;
}
//-------------------------------------------------------------------------
int Connectivity2::Superpose(int fixed) {
//-------------------------------------------------------------------------
  // Apply lsq fit of 'moving' set of atoms on 'fixed' atoms
  if ( pAtom1.size() < 3)  return -1;
  std::vector<PCAtom>::iterator i;
  PPCAtom A1 = new PCAtom[pAtom1.size()];
  PPCAtom A2 = new PCAtom[pAtom2.size()];
  i=pAtom1.begin(); int j = 0;
  while(i!=pAtom1.end()){  A1[j++] = *i; i++; }
  i=pAtom2.begin(); j = 0;
  while(i!=pAtom2.end()){  A2[j++] = *i; i++; }

  //cout << "Connectivity2::Superpose fixed " << fixed <<" " <<pAtom1.size() << endl;

  if (fixed == 1) {
    // Move the second set of atoms
    i = pAtom2.begin();
    PCMMANManager M2 = (PCMMANManager)(PCMMDBFile)(*i)->GetCoordHierarchy();
    M2->TransformToSuperposeAtoms( A2,pAtom2.size() ,A1 );

  } else {
    i = pAtom1.begin();
    PCMMANManager M1 = (PCMMANManager)(PCMMDBFile)(*i)->GetCoordHierarchy();
    M1->TransformToSuperposeAtoms( A1,pAtom1.size() ,A2 );

  }
  //cout << "Connectivity2::Superpose done " << endl;
  return 0;
 
}


int Connectivity2::MatchGraphs(PCResidue pRes1,const pstr altLoc1,
                               PCResidue pRes2,const pstr altLoc2,
                               int Hflag,int tag,float fracMinMatch,
                               bool keepmatch )  {
  
  PCGraph G2,G1;
  PCGraphMatch U = NULL;
  ivector      F1,F2;
  realtype     p1,p2;
  int     htype;
  int nInResidue1,nInResidue2,minMatch;
  PPCAtom atoms1,atoms2;
  int ii,nat1,nat2,nMatched,natMatch;
  bool keepmatchOK;
  std::vector<unsigned int> conns;
  std::vector <int> best_matches;
  std::vector <int>:: iterator i;
   
  G1 = new CGraph ( pRes1,altLoc1 );
  if (Hflag>=1) {
    htype = getElementNo(pstr("H"));
    if (Hflag==2)  G1->HideType    ( htype );
       else  G1->ExcludeType ( htype );
  }                                                                     
  G1->Build ( False );
  nInResidue1 = G1->GetNofVertices();
  if (nInResidue1<=0)  {
    delete G1;
    return -1;
  }
  G2 = new CGraph ( pRes2,altLoc2 );
  if (Hflag>=1) {
    htype = getElementNo(pstr("H"));
    if (Hflag==2)  G2->HideType    ( htype );
       else  G2->ExcludeType ( htype );
  }                                                                     
  G2->Build ( False );
  nInResidue2 = G2->GetNofVertices();
  if (nInResidue2<=0)  {
    delete G1;
    delete G2;
    return -2;
  }

  // Start with high minMatch (minimum number of matched atoms)
  // and work down
  for (int ifrac=9;ifrac>=1;ifrac--) {
    minMatch = (ifrac*min(nInResidue1,nInResidue2))/10;
    //cout << "minMatch " << minMatch <<endl;
    U = new CGraphMatch();
    U->MatchGraphs ( G1,G2,minMatch ); 
    nMatched =  U->GetNofMatches();
    if (nMatched>0) break;
    delete U;
    U = NULL;
  }


  //cout << "nMatched " << nMatched << endl;
  pRes1->GetAtomTable (atoms1, nat1);
  pRes2->GetAtomTable (atoms2, nat2);
  
    
  if (nMatched>0) {
    int maxM=0; int maxI=-1; int maxI2=-1;
    for (ii=0;ii<nMatched;ii++) {
      U->GetMatch ( ii,F1,F2,natMatch,p1,p2 );
      //cout << "ii,natMatch " << ii << " " << natMatch<< endl;
      if ( natMatch>=maxM) {
        keepmatchOK = true;
        if (keepmatch) {
          for (int j=0;j<natMatch;j++) {
            if (FindNofConnections(atoms1[F1[j]-1],1)>0 ||
                FindNofConnections(atoms2[F2[j]-1],2)>0 ) {
	        conns = FindConnections( atoms1[F1[j]-1],atoms2[F2[j]-1],false);
	      if  (conns.size()<1) keepmatchOK = false;
              //cout << "conns.size " << j << " " <<conns.size()<< " " << keepmatchOK << endl;
            }
	  }
        }
        if (keepmatchOK) {
          if (natMatch>maxM) {
            maxM=natMatch;
            maxI = ii;
            best_matches.clear();
            best_matches.push_back(ii);
          } else if (natMatch==maxM) {
            best_matches.push_back(ii);
          }
        }
      }
    }
    //cout << "maxM " << best_matches.size() << " " << maxM << endl;
    // Multiple best graph matches - decide by superposing and
    // chosing match with best lsq fit
    if (best_matches.size()>1 && maxM>3 ) {
      PPCAtom A1 = new PCAtom[maxM];
      PPCAtom A2 = new PCAtom[maxM];
      for (int j=0;j<natMatch;j++) {
        A1[j] = new CAtom();
        A2[j] = new CAtom();
      }
      mat44 TMatrix;
      double distSqMin=9999999.9,distSq;
      int rv;
      i=best_matches.begin();
      while (i!=best_matches.end()) {
        U->GetMatch (*i,F1,F2,natMatch,p1,p2 );
        //cout << "*i,natMatch " << *i << " " << natMatch << endl;
        for (int j=1;j<=natMatch;j++) {
          //cout << "j " << j << " " << F1[j] << " " << F2[j] << endl;
          A1[j-1]->Copy(atoms1[F1[j]-1]);
          A2[j-1]->Copy(atoms2[F2[j]-1]);
        }
        rv = SuperposeAtoms ( TMatrix, A1,maxM , A2 );
        //cout << "rv " << rv << endl;
        if (rv == SPOSEAT_Ok) {
          distSq = 0.0;
          for (int j=0;j<natMatch;j++) { 
            A1[j]->Transform(TMatrix);
            distSq = distSq +pow(A1[j]->x-A2[j]->x,2) +pow(A1[j]->y-A2[j]->y,2) +  pow(A1[j]->z-A2[j]->z,2);
          }
          //cout << "distSq " << distSq << endl;
	  if (distSq<distSqMin) {
            distSqMin = distSq;
            maxI2 = *i;
            //cout << "distSqMin,maxI2 " << distSqMin << " " << maxI2 << endl;
          }
        }
        i++;
      }
    }
    if (maxI2>=0)  maxI = maxI2;
    if (maxI>=0) {    
      U->GetMatch ( maxI,F1,F2,natMatch,p1,p2 );
      //cout << "maxI,natMatch " << maxI << " " << natMatch << endl;
      for (int j=1;j<=natMatch;j++) {
        //cout << "j,F1[j],F2[j]" << j << " " << F1[j] << " " << F2[j] << endl;
        AddConnection(atoms1[F1[j]-1],atoms2[F2[j]-1],"",tag);      
      }
    } else {
      natMatch = 0;
    } 
  }

  if (U) delete U;
  delete G2;
  delete G1;
  return natMatch;
}
