//   CCP4 Molecular Graphics Program
//
//   Copyright (C) Stuart McNicholas and Liz Potterton 2004.
//
//   This program is free software and is distributed under the terms
//   and conditions of the CCP4 licence agreement as `Part 0' (Annex 2)
//   software, which is version 2.1 of the GNU Lesser General Public
//   Licence (LGPL) with the following additional clause:
//
//      `You may also combine or link a "work that uses the Library"
//      to produce a work containing portions of the Library, and
//      distribute that work under terms of your choice, provided that
//      you give prominent notice w/ith each copy of the work that the
//      specified version of the Library is used in it, and that you
//      include or provide public access to the complete corresponding
//      machine-readable source code for the Library including whatever
//      changes were used in the work. (i.e. If you make changes to the
//      Library you must distribute those, but you do not need to
//      distribute source or object code to those portions of the work
//      not covered by this licence.)'
//
//   Note that this clause grants an additional right and does not
//   impose any additional restriction, and so does not affect
//   compatibility with the GNU General Public Licence (GPL). If you
//   wish to negotiate other terms, please contact the maintainer.
//   You can redistribute it and/or modify the program under the terms
//   of the GNU Lesser General Public License as published by the Free
//   Software Foundation; either version 2.1 of the License, or (at
//   your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
//
//   You should have received a copy of the CCP4 licence and/or GNU
//   Lesser General Public License along with this program; if not,
//   write to the CCP4 Secretary, Daresbury Laboratory, Warrington
//   WA4 4AD, UK. The GNU Lesser General Public can also be obtained
//   by writing to the Free Software Foundation, Inc., 59 Temple Place,
//   Suite 330, Boston, MA 02111-1307 USA


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <list>
#include <time.h>
#include "mgtree.h"
#include "quat.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

using namespace std;

int TreeVertex::FindDepth(void) const {
  TreeVertex* par;
  int depth=0;

  if(GetParentID()==-1)
    return depth;

  par = GetParent();
  depth++;

  while(par->GetParentID()!=-1){
    par = par->GetParent();
    depth++;
  }

  return depth;
}

std::ostream& operator<<(std::ostream& c, TreeVertex a){

   double rtod = 180.0 / M_PI;
   int i;

   for(i=0;i<a.GetNumberOfChildren();i++)
      c << a.GetChildren()[i];
   c << "\n";
   for(i=0;i<a.FindDepth();i++){c << " ";}
   c << "[Parent] Length:" << a.GetParentDistance() << " Angle:" << a.GetParentBondAngle()*rtod << " Torsion:" << a.GetParentDihedralAngle()*rtod << "\n";
   for(i=0;i<a.GetNumberOfChildren();i++)
      c << *(a.GetChildren()[i]);

   return c;

}

TreeVertex::TreeVertex(){
  parent_id = -1;
  parent_dist = 0.0;
  parent_bond_angle = 0.0;
  parent_dihedral_angle = 0.0;
  parent = 0;
}

int Tree::FindMaxDepth(){
  int depth, maxdepth = 0;

  for(unsigned int i=0;i<coords.size();i++){
    depth = coords[i]->FindDepth();
    if(depth>maxdepth)
      maxdepth = depth;
  }
  return maxdepth;
}
void TreeVertex::SetAngles(void){

  int depth = FindDepth();
  if(depth>0)
     SetParentDistance(LineLength(
     coord,
     parent->coord
     ));
  if(depth>1)
     SetParentBondAngle(Angle(
     coord,
     parent->coord,
     parent->parent->coord
     ));
  if(depth>2)
     SetParentDihedralAngle(DihedralAngle(
     coord,
     parent->coord,
     parent->parent->coord,
     parent->parent->parent->coord
     ));
  if(depth==1&&GetParent()->GetNumberOfChildren()>2)
     SetParentDihedralAngle(DihedralAngle(
     coord,
     parent->coord,
     parent->children[0]->coord,
     parent->children[1]->coord
     ));
  if(depth==2&&GetParent()->GetNumberOfChildren()>1)
     SetParentDihedralAngle(DihedralAngle(
     coord,
     parent->coord,
     parent->children[0]->coord,
     parent->parent->coord
     ));
}

TreeVertex::~TreeVertex(){
}

Tree::Tree(){
}

Tree::~Tree(){
}

Tree::Tree(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists, const  std::vector<std::vector<Cartesian> > &ext_cartesians){
  if(ext_cartesians.size()>0)
    SetCoords(SelAtoms_in, ioff, conn_lists, ext_cartesians);
  else
    SetCoords(SelAtoms_in, ioff, conn_lists);

}

void Tree::SetCoords(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists){
  std::vector<std::vector<Cartesian> > ext_carts(SelAtoms_in.size());
  SetCoords(SelAtoms_in,ioff,conn_lists,ext_carts);
}

void Tree::SetCoords(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists, const std::vector<std::vector<Cartesian> >&ext_carts_in){

  clock_t tv1, tv2;
  clock_t tvs;
  bool print_timing = false;

  tvs = clock();
  tv1 = clock();
  std::vector<std::vector<Cartesian> >ext_carts;

  //First thing to do is apply permutation to swap ioff(start pos) with 0.
  std::vector<Cartesian> SelAtoms = SelAtoms_in;
  connectivity = conn_lists;


  for(unsigned ii=0;ii<connectivity.size();ii++){
     for(unsigned jj=0;jj<connectivity[ii].size();jj++){
	
	if(std::find(connectivity[connectivity[ii][jj]].begin(),connectivity[connectivity[ii][jj]].end(),ii)!=connectivity[connectivity[ii][jj]].end()){
	   // PE std::cout << "OK\n";
	}else
	   if(connectivity[ii][jj]!=ii) { // That would be a disaster.
	      // PE std::cout << "Hmm: " << ii << " not in " << connectivity[ii][jj] << "\n";
 
	      connectivity[connectivity[ii][jj]].push_back(ii);
	   }
     }
  }
  
//   for(unsigned ii=0;ii<connectivity.size();ii++){
//      std::cout << ii << ": ";
//      for(unsigned jj=0;jj<connectivity[ii].size();jj++){
//         std::cout << connectivity[ii][jj] << " ";
//      }
//      std::cout << "\n";
//   }

  
  ext_carts = ext_carts_in;

  start = ioff;
  if(ioff!=0){
    std::vector<std::vector<int> > conn_perm(SelAtoms.size());
    std::vector<std::vector<Cartesian> > ext_perm(SelAtoms.size());
    std::vector<Cartesian> SelAtoms_perm(SelAtoms.size());
    permutation = std::vector<int>(SelAtoms.size());

    SelAtoms_perm[0] = SelAtoms[ioff];
    conn_perm[0] = connectivity[ioff];
    ext_perm[0] = ext_carts[ioff];
    permutation[ioff] = 0;
    for(unsigned i=0;i<SelAtoms.size();i++){
      if((int)i>ioff){
        SelAtoms_perm[i] = SelAtoms[i];
        conn_perm[i] = connectivity[i];
        ext_perm[i] = ext_carts[i];
        permutation[i] = i;
      }
      if((int)i<ioff){
        SelAtoms_perm[i+1] = SelAtoms[i];
        conn_perm[i+1] = connectivity[i];
        ext_perm[i+1] = ext_carts[i];
        permutation[i] = i+1;
      }
    }
    for(unsigned i=0;i<SelAtoms.size();i++){
      for(unsigned j=0;j<conn_perm[i].size();j++){
        conn_perm[i][j] = permutation[conn_perm[i][j]];
      }
    }
    ext_carts = ext_perm;
    connectivity = conn_perm;
    SelAtoms = SelAtoms_perm;
  }

  coords.clear();
  extra_bonded_pairs.clear();

  tv2 = clock();

  if(SelAtoms.size()==0)
    return;

  if ( print_timing) printf("Time to do setup: %f\n",(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0) / (double) 1000.0);

  std::vector<std::vector<Cartesian> >::iterator ext_iter = ext_carts.begin();
  tv1 = clock();
  for(unsigned i=0;i<SelAtoms.size();i++){
    coords.push_back(new TreeVertex());
    coords[i]->SetCoord(SelAtoms[i]);  // The above doesn't do much at moment.
    std::vector<Cartesian>::const_iterator this_ext_iter=ext_iter->begin();
    while(this_ext_iter!=ext_iter->end()){
      coords[i]->AddExternalChild(*this_ext_iter);
      this_ext_iter++;
    }
    ext_iter++;
  }
  tv2 = clock();
  if (print_timing) printf("Time to do copy atoms: %f\n",(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0) / (double) 1000.0);

  tv1 = clock();
  CalculateTree();
  tv2 = clock();
  if (print_timing) printf("Time to do calculate tree: %f\n",(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0) / (double) 1000.0);

  tv1 = clock();
  std::vector<TreeVertex*>::iterator k = coords.begin();
  while(k!=coords.end()){
    if((*k)->GetParentID()==-1) SetDummy(*k);
    k++;
  }
  tv2 = clock();
  if (print_timing) printf("Time to do set dummies: %f\n",(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0) / (double) 1000.0);
  if (print_timing) printf("Time to do all tree stuff: %f\n",(tv2 - tvs)/(CLOCKS_PER_SEC / (double) 1000.0) / (double) 1000.0);

}

void Tree::RecurseCalculateTree(TreeVertex* coord){

  //std::cout << "RecurseCalculateTree " << coord->GetID() << "\n"; std::cout.flush();
;
  TreeVertex* coordk;

  if(scanned[coord->GetID()]==1){
    //std::cout << "Already done this, returning\n"; std::cout.flush();
    return;
  }

  scanned[coord->GetID()] = 1;

  if(connectivity[coord->GetID()].size()==0){
    return;
  }

  std::vector<int>::iterator k=connectivity[coord->GetID()].begin();

  std::vector<int>::iterator kend = connectivity[coord->GetID()].end();
  while(k!=kend){
    if(coord->GetParentID()!=*k){
       coordk = coords[*k];
       if(coordk->GetParentID()==-1&&scanned[*k]!=1){
         coordk->SetParentID(coord->GetID()); 
         coordk->SetParent(coord);
         coord->AddChild(coordk);
         coordk->SetAngles();
//          std::cout << "Adding " << coordk->GetID() << " as a child of " << coord->GetID() << "\n"; std::cout.flush();
         RecurseCalculateTree(coordk);
       } else {  //added by Joel Bard
         //but coord may not have a parent so add it as
         //a child to coordk
         if( coord->GetParentID()==-1&&scanned[*k]!=1 ) {
         std::cout << coordk->GetID() << "\n"; std::cout.flush();
         std::cout << coord->GetID() << "\n"; std::cout.flush();
           coord->SetParentID(coordk->GetID());
           coord->SetParent(coordk);
           coordk->AddChild(coord);
           coord->SetAngles();
         }
       }
    }
    k++;
  } 
}

void Tree::CalculateTree(){

  unsigned int i;
  std::vector<int>::iterator k;
  std::vector<std::pair<int,int> > bp;
  TreeVertex* coord;
  TreeVertex* coordk;

  if(getenv("USE_OLD_TREE")){
  for(i=0;i<coords.size();i++){
    coord = coords[i];
    coord->SetID(i);
    k=connectivity[i].begin();
    std::vector<int>::iterator kend = connectivity[i].end();
    while(k!=kend){
      if(coord->GetParentID()!=*k){
	coordk = coords[*k];
        if(coordk->GetParentID()==-1){
	  //printf("Adding tree connection:%d -> %d\n",(int)i,*k);
           coordk->SetParentID(i); 
           coordk->SetParent(coord);
           coord->AddChild(coordk);
           coordk->SetAngles();
	}
      }
      k++;
    }
  }

  }else{

  for(i=0;i<coords.size();i++){
    coord = coords[i];
    coord->SetID(i);
  }

  scanned.clear();
  scanned.resize(coords.size());
  coord = coords[0];
  RecurseCalculateTree(coord);
  for(i=0;i<coords.size();i++)
    if(!scanned[i]){
      coord = coords[i];
      RecurseCalculateTree(coord);
    }
  }

  int unbonded = 0;
  for(i=0;i<coords.size();i++){
    if(coords[i]->FindDepth()!=1){
      coords[i]->SetAngles(); // Shouldn't need the if ...
    }else{
      unbonded++;
    }
  }

  for(i=0;i<coords.size();i++){
    coords[i]->SetID(i);
    k=connectivity[i].begin();
    while(k!=connectivity[i].end()){
      if(coords[i]->GetParentID()!=*k&&coords[*k]->GetParentID()!=int(i)){
	  if(*k>int(i)){
	    bp.push_back(std::pair<int,int>((int)i,*k));
	  }
	  else{
	    bp.push_back(std::pair<int,int>(*k,(int)i));
	  }
      }
      k++;
    }
  }

  std::sort(bp.begin(),bp.end(),bond_pair_cmp());
  unique_copy(bp.begin(),bp.end(),back_inserter(extra_bonded_pairs));

  //cout << "There are " << bp.size() << " extra bonded pairs.\n";
  //cout << "There are " << unbonded << " unbonded vertices.\n";

}

std::ostream& operator<<(std::ostream& c, Tree a){


  unsigned int i;

  for(i=0;i<a.coords.size();i++){
    if(a.coords[i]->GetParentID()>-1) {
      c << "Atom[" << i << "]: " << a.coords[i]->GetCoord() << " has parent atom[" << a.coords[i]->GetParent()->GetID() << "]: " << a.coords[i]->GetParent()->GetCoord() << "\n";
    }else{
      c << "Atom[" << i << "]: " << a.coords[i]->GetCoord() << " has no parent\n";
    }
  }

  c << "\nTrees ....\n";
  for(i=0;i<a.coords.size();i++){
    if(a.coords[i]->GetParentID()==-1)
      c << *(a.coords[i]);
  }


  return c;
}

int TreeVertex::GetNumberOfDescendants(void) const {
  int ndescendants = 0;
  ndescendants += GetNumberOfChildren();
  TreeVertex *child;
  for(int i=0;i<GetNumberOfChildren();i++){
    child = GetChild(i);
    ndescendants += child->GetNumberOfDescendants();
  }
  return ndescendants;
}
 
void Tree::SetDummy(TreeVertex *coord){

  int ndescendants = coord->GetNumberOfDescendants();

  /* 
  * This bit may be stuffed if first three atoms are collinear.
  */

  Cartesian p1;
  Cartesian p2;
  Cartesian p3;
  Cartesian result;
  Cartesian Dummy;
  Cartesian Dummy2;

  if(ndescendants>1){
    if(coord->GetChild(0)->GetNumberOfChildren()>0){
      p3 = Cartesian(coord->GetCoord().get_x(), coord->GetCoord().get_y(), coord->GetCoord().get_z());
      p2 = Cartesian(coord->GetChild(0)->GetCoord().get_x(), coord->GetChild(0)->GetCoord().get_y(), coord->GetChild(0)->GetCoord().get_z());
      p1 = Cartesian(coord->GetChild(0)->GetChild(0)->GetCoord().get_x(), coord->GetChild(0)->GetChild(0)->GetCoord().get_y(), coord->GetChild(0)->GetChild(0)->GetCoord().get_z());
      result = GetCartFrom3Carts(p3,1.0,p2,M_PI/2,p1,-M_PI/2);
      Dummy = result;
      p1 = p2;
      p2 = p3;
      p3 = result;
      result = GetCartFrom3Carts(p3,1.0,p2,M_PI/2,p1,M_PI/2);
      Dummy2 = result;
    }else{
      if(coord->GetNumberOfChildren()==2){
        p3 = Cartesian(coord->GetCoord().get_x(), coord->GetCoord().get_y(), coord->GetCoord().get_z());
        p2 = Cartesian(coord->GetChild(0)->GetCoord().get_x(), coord->GetChild(0)->GetCoord().get_y(), coord->GetChild(0)->GetCoord().get_z());
        p1 = Cartesian(coord->GetChild(1)->GetCoord().get_x(), coord->GetChild(1)->GetCoord().get_y(), coord->GetChild(1)->GetCoord().get_z());
	double angle = M_PI-Angle(coord->GetChild(1)->GetCoord(),coord->GetCoord(),coord->GetChild(0)->GetCoord())/2;
	result = GetCartFrom3Carts(p3,1.0,p2,angle,p1,M_PI);
	Dummy  = result;
	p1 = p2;
	p2 = p3;
	p3 = result;
	result = GetCartFrom3Carts(p3,1.0,p2,M_PI/2,p1,M_PI);
	Dummy2 = result;
      }else{
        //printf("Dodgy stuff....\n");
        p3 = Cartesian(coord->GetChild(0)->GetCoord().get_x(), coord->GetChild(0)->GetCoord().get_y(), coord->GetChild(0)->GetCoord().get_z());
        p2 = Cartesian(coord->GetChild(1)->GetCoord().get_x(), coord->GetChild(1)->GetCoord().get_y(), coord->GetChild(1)->GetCoord().get_z());
        p1 = Cartesian(coord->GetChild(2)->GetCoord().get_x(), coord->GetChild(2)->GetCoord().get_y(), coord->GetChild(2)->GetCoord().get_z());
	double x1 = p1.get_x(), x2 = p2.get_x(), x3 = p3.get_x(); 
	double y1 = p1.get_y(), y2 = p2.get_y(), y3 = p3.get_y(); 
	double z1 = p1.get_z(), z2 = p2.get_z(), z3 = p3.get_z(); 
	double A = y1 * (z2 - z3) + y2 * (z3 - z1) + y3 * (z1 - z2);
	double B = z1 * (x2 - x3) + z2 * (x3 - x1) + z3 * (x1 - x2);
	double C = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
	Dummy.set_x(coord->GetCoord().get_x() - A);
	Dummy.set_y(coord->GetCoord().get_y() - B);
	Dummy.set_z(coord->GetCoord().get_z() - C);
	p1 = p3;
	p3 = Cartesian(Dummy.get_x(), Dummy.get_y(), Dummy.get_z());
	p2 = Cartesian(coord->GetCoord().get_x(), coord->GetCoord().get_y(), coord->GetCoord().get_z());
	result = GetCartFrom3Carts(p3,1.0,p2,M_PI/2,p1,0);
	Dummy2 = result;
      }
    }
  }else{
    if(ndescendants==0){
      //char S[50];
      //printf("%s has zero descendants\n",coord->GetCoord()->GetCoordID(S));
      Dummy.set_x(coord->GetCoord().get_x()+1.0);
      Dummy.set_y(coord->GetCoord().get_y());
      Dummy.set_z(coord->GetCoord().get_z());
      Dummy2.set_x(coord->GetCoord().get_x()+1.0);
      Dummy2.set_y(coord->GetCoord().get_y()+1.0);
      Dummy2.set_z(coord->GetCoord().get_z());
    }else{
      //char S[50];
      //printf("%s has only one descendant\n",coord->GetCoord()->GetCoordID(S));
    // Need a better way !!!!!
    p3 = Cartesian(coord->GetCoord().get_x(), coord->GetCoord().get_y(), coord->GetCoord().get_z());
    p2 = Cartesian(coord->GetChild(0)->GetCoord().get_x(), coord->GetChild(0)->GetCoord().get_y(), coord->GetChild(0)->GetCoord().get_z());
    p1 = p3 - p2;
    p1.normalize();

    Cartesian xaxis(1,0,0);
    Cartesian yaxis(0,1,0);
    Cartesian zaxis(0,0,1);

    p3 = p1.CrossProduct(p1,xaxis);

    int axis = 0;
    double size = p3.length();
    double minsize = size;
    p3 = p1.CrossProduct(p1,yaxis);
    size = p3.length();
    if(size<minsize) { minsize = size; axis = 1; }
    p3 = p1.CrossProduct(p1,zaxis);
    size = p3.length();
    if(size<minsize) { minsize = size; axis = 2; }

    if(axis==0){
      Dummy.set_x(coord->GetCoord().get_x()+1.0);
      Dummy.set_y(coord->GetCoord().get_y());
      Dummy.set_z(coord->GetCoord().get_z());
    }
    if(axis==1){
      Dummy.set_x(coord->GetCoord().get_x());
      Dummy.set_y(coord->GetCoord().get_y()+1.0);
      Dummy.set_z(coord->GetCoord().get_z());
    }
    if(axis==2){
      Dummy.set_x(coord->GetCoord().get_x());
      Dummy.set_y(coord->GetCoord().get_y());
      Dummy.set_z(coord->GetCoord().get_z()+1.0);
    }            
    p1 = Cartesian(coord->GetChild(0)->GetCoord().get_x(), coord->GetChild(0)->GetCoord().get_y(), coord->GetChild(0)->GetCoord().get_z());
    p2 = Cartesian(coord->GetCoord().get_x(), coord->GetCoord().get_y(), coord->GetCoord().get_z());
    p3 = Cartesian(Dummy.get_x(), Dummy.get_y(), Dummy.get_z());
    result = GetCartFrom3Carts(p3,1.0,p2,M_PI/2,p1,M_PI/2);
    Dummy2 = result;
    }
  }
  coord->SetDummy(Dummy,Dummy2);
  coord->SetParentDistance(LineLength(Dummy,coord->GetCoord()));
  coord->SetParentBondAngle(Angle(Dummy2,Dummy,coord->GetCoord()));
  
  TreeVertex *child;
  TreeVertex *grandchild;

  for(int i=0;i<coord->GetNumberOfChildren();i++){
     child = coord->GetChildren()[i];
     child->SetParentBondAngle(Angle(
     child->GetCoord(),
     child->GetParent()->GetCoord(),
     Dummy));
     child->SetParentDihedralAngle(DihedralAngle(
     child->GetCoord(),
     coord->GetCoord(),
     Dummy,
     Dummy2));
     for(int j=0;j<child->GetNumberOfChildren();j++){
       grandchild = child->GetChildren()[j];
       grandchild->SetParentDihedralAngle(DihedralAngle(
       grandchild->GetCoord(),
       grandchild->GetParent()->GetCoord(),
       grandchild->GetParent()->GetParent()->GetCoord(),
       Dummy));
     }
  }

}

void Tree::ExtendBranchCartesians(const Cartesian &p1_in, const Cartesian &p2_in, const Cartesian &p3_in, TreeVertex* coord, std::vector<Cartesian> &cartesians)const {
  Cartesian result;
  Cartesian p1 = p1_in;
  Cartesian p2 = p2_in;
  Cartesian p3 = p3_in;

  result = GetCartFrom3Carts(p3,coord->GetParentDistance(),p2,coord->GetParentBondAngle(),p1,coord->GetParentDihedralAngle());
  cartesians[coord->GetID()] = result;
  p1 = p2;
  p2 = p3;
  p3 = result;
  for(int i=0;i<coord->GetNumberOfChildren();i++){
     ExtendBranchCartesians(p1,p2,p3,coord->GetChild(i),cartesians);
  }
}


Cartesian Tree::GetCartesian(TreeVertex *unknown) const {

  Cartesian result;
  Cartesian p1;
  Cartesian p2;
  Cartesian p3;

  std::vector <TreeVertex*> branch = unknown->GetBranch();
  std::vector<TreeVertex*>::iterator k = branch.begin();

  p1 =  Cartesian((*k)->Dummy2.get_x(), (*k)->Dummy2.get_y(), (*k)->Dummy2.get_z());
  p2 =  Cartesian((*k)->Dummy.get_x(), (*k)->Dummy.get_y(), (*k)->Dummy.get_z());
  p3 =  Cartesian((*k)->GetCoord().get_x(), (*k)->GetCoord().get_y(), (*k)->GetCoord().get_z());

  if(unknown->GetParentID()==-1){
    return result;
  }

  k++;

  while(k!=branch.end()){
    result = GetCartFrom3Carts(p3,(*k)->GetParentDistance(),p2,(*k)->GetParentBondAngle(),p1,(*k)->GetParentDihedralAngle());
    p1 = p2;
    p2 = p3;
    p3 = result;
    k++;
  }

  //printf("branch size:%d\n",branch.size());

  return result;
}

std::vector <TreeVertex*> TreeVertex::GetBranch(){


  std::vector<TreeVertex*> branch(FindDepth()+1);
  TreeVertex *par = this;

  for(int i=par->FindDepth();i>=0;i--){
    branch[i] = par;
    par = par->GetParent();
  }

  return branch;
}

void Tree::RecurseZMatrix(std::ostream &c,const TreeVertex *vertex, const std::vector<std::string> &labels){
  static int depth;
  depth++;
  int id = vertex->GetID();
  if(depth==1){
    //printf("%s,Q1,%f,Q2,%f\n",labels[id].c_str(),vertex->GetParentDistance(),-vertex->GetParentBondAngle()*180/M_PI);
    c << labels[id] << ",Q1," << vertex->GetParentDistance() << ",Q2," << -vertex->GetParentBondAngle()*180/M_PI << "\n";
  }
  if(depth==2){
    TreeVertex *parent = vertex->GetParent();
    int pid = parent->GetID();
    if(0){
      c << labels[id] << "," << labels[pid] << "," << vertex->GetParentDistance() << ",Q1," << vertex->GetParentBondAngle()*180/M_PI << "\n";
      /*
	printf("%s,%s,%f,Q1,%f\n",
	labels[id].c_str(),labels[pid].c_str(),
	vertex->GetParentDistance(),
	vertex->GetParentBondAngle()*180/M_PI);
      */
    }else{
      c << labels[id] << "," << labels[pid] << "," << vertex->GetParentDistance() << ",Q1," << vertex->GetParentBondAngle()*180/M_PI << ",Q2," <<  vertex->GetParentDihedralAngle()*180/M_PI << "\n";
      /*
	printf("%s,%s,%f,Q1,%f,Q2,%f\n",
	labels[id].c_str(),labels[pid].c_str(),
	vertex->GetParentDistance(),
	vertex->GetParentBondAngle()*180/M_PI,
	vertex->GetParentDihedralAngle()*180/M_PI);
      */
    }
  }
  if(depth==3){
    TreeVertex *parent = vertex->GetParent();
    int pid = parent->GetID();
    TreeVertex *grandparent = parent->GetParent();
    int gpid = grandparent->GetID();
    c << labels[id] << "," << labels[pid] << "," << vertex->GetParentDistance() << "," << labels[gpid] << "," << vertex->GetParentBondAngle()*180/M_PI << ",Q1," << vertex->GetParentDihedralAngle()*180/M_PI << "\n";
    /*
      printf("%s,%s,%f,%s,%f,Q1,%f\n",
      labels[id].c_str(),labels[pid].c_str(),
      vertex->GetParentDistance(),labels[gpid].c_str(),
      vertex->GetParentBondAngle()*180/M_PI,
      vertex->GetParentDihedralAngle()*180/M_PI);
    */
  }
  if(depth>3){
    TreeVertex *parent = vertex->GetParent();
    int pid = parent->GetID();
    TreeVertex *grandparent = parent->GetParent();
    int gpid = grandparent->GetID();
    TreeVertex *greatgrandparent = grandparent->GetParent();
    int ggpid = greatgrandparent->GetID();
    c << labels[id] << "," << labels[pid] << "," << vertex->GetParentDistance() << "," << labels[gpid] << "," << vertex->GetParentBondAngle()*180/M_PI << "," << labels[ggpid] << "," << vertex->GetParentDihedralAngle()*180/M_PI << "\n";
    /*
      printf("%s,%s,%f,%s,%f,%s,%f\n",
      labels[id].c_str(),labels[pid].c_str(),
      vertex->GetParentDistance(),labels[gpid].c_str(),
      vertex->GetParentBondAngle()*180/M_PI,labels[ggpid].c_str(),
      vertex->GetParentDihedralAngle()*180/M_PI);
    */
  }
   for(int i=0;i<vertex->GetNumberOfChildren();i++)
    RecurseZMatrix(c,vertex->GetChild(i),labels);
  depth--;
}

void Tree::PrintZMatrix(const std::vector<std::string> &labels){
  PrintZMatrix(std::cout,labels);
}

void Tree::PrintZMatrix(std::ostream &c, const std::vector<std::string> &labels){

  for(int i=0;i<GetNumberOfVertices();i++){
    TreeVertex *vertex = GetCoord(i);
    if(vertex->GetParentID()==-1){
      Cartesian Q2 = GetCoord(0)->Dummy2;
      Cartesian Q1 = GetCoord(0)->Dummy;
      c << "Q2\n" << "Q1,Q2," << LineLength(Q2,Q1) << "\n";
      RecurseZMatrix(c,vertex,labels);
    }
  }
}

TreeVertex* Tree::GetCoord(int i, bool permuted) const {
  if(start>0&!permuted)
    return coords[permutation[i]];
  return coords[i];
}

std::vector <TreeVertex*> Tree::GetCoords(bool permuted) const {

  if(start>0&!permuted){
    std::vector<TreeVertex*> perm_coords;
    for(int i=0;i<GetNumberOfVertices();i++){
       perm_coords.push_back(coords[permutation[i]]);
    }
    return perm_coords;
  }
    
  return coords;
}

Cartesian Tree::GetCartesian(int i, bool permuted) const {
  return GetCartesian(GetCoord(i,permuted));
}

std::vector<Cartesian> Tree::GetAllCartesians(bool permuted) const {

  std::vector<Cartesian> cartesians = std::vector<Cartesian>(coords.size());

  Cartesian p1;
  Cartesian p2;
  Cartesian p3;

  std::vector<TreeVertex*>::const_iterator k = coords.begin();

  while(k!=coords.end()){
    if((*k)->GetParentID()==-1) {
      p1 =  Cartesian((*k)->Dummy2.get_x(), (*k)->Dummy2.get_y(), (*k)->Dummy2.get_z());
      p2 =  Cartesian((*k)->Dummy.get_x(), (*k)->Dummy.get_y(), (*k)->Dummy.get_z());
      p3 =  Cartesian((*k)->GetCoord().get_x(), (*k)->GetCoord().get_y(), (*k)->GetCoord().get_z());
      cartesians[(*k)->GetID()]   = p3;
      for(int i=0;i<(*k)->GetNumberOfChildren();i++){
         ExtendBranchCartesians(p1,p2,p3,(*k)->GetChild(i),cartesians);
      }
    }
    k++;
  }

  if(start>0&!permuted){
    std::vector<Cartesian> perm_carts;
    for(int i=0;i<GetNumberOfVertices();i++){
       perm_carts.push_back(cartesians[permutation[i]]);
    }
    return perm_carts;
  }
  return cartesians;

}

void Tree::RotateAboutBond(int atom_in, int child_in, double TorsionAngle, bool permuted){
  int atom = atom_in;
  int child = child_in;
  
  if(start>0&!permuted){
    atom = permutation[atom];
    child = permutation[child];
  }

  TreeVertex *C = coords[child];
  TreeVertex *B = coords[atom];

  if(C->GetParentID()==atom){
  } else if(B->GetParentID()==child){
    C = coords[atom];
    B = coords[child];
  }else {
    std::cout << "These are not related\n";
    return;
  }
  if(C->GetNumberOfChildren()==0)
    return;

  std::vector<TreeVertex*> children = C->GetChildren();
  std::vector<TreeVertex*>::const_iterator child_iter=children.begin();

  while(child_iter!=children.end()){
    //std::cout << "Original Dihedral Angle: " << (*child_iter)->GetParentDihedralAngle() << "\n";
    (*child_iter)->SetParentDihedralAngle((*child_iter)->GetParentDihedralAngle()+TorsionAngle);
    //std::cout << "New Dihedral Angle: " << (*child_iter)->GetParentDihedralAngle() << "\n";
    child_iter++;
  }
  //std::cout << "\n"; Rotated:

}

void Tree::SetDihedralAngle(int atom_in, int child_in, double TorsionAngle, bool permuted){
  //std::cout << "SetDihedralAngle\n\n";

  int atom = atom_in;
  int child = child_in;
  
  if(start>0&!permuted){
    atom = permutation[atom];
    child = permutation[child];
  }

  TreeVertex *C = coords[child];
  TreeVertex *B = coords[atom];

  if(C->GetParentID()==atom){
  } else if(B->GetParentID()==child){
    C = coords[atom];
    B = coords[child];
  }else {
    std::cout << "These are not related\n";
    return;
  }

  if(C->GetNumberOfChildren()==0)
    return;

  double angle_orig;
  double diff;

  std::vector<TreeVertex*> children = C->GetChildren();
  std::vector<TreeVertex*>::const_iterator child_iter=children.begin();

  angle_orig = (*child_iter)->GetParentDihedralAngle();
  diff = TorsionAngle - angle_orig;
  (*child_iter)->SetParentDihedralAngle(TorsionAngle);
  child_iter++;

  while(child_iter!=children.end()){
    //std::cout << "Original Dihedral Angle: " << (*child_iter)->GetParentDihedralAngle() << "\n";
    (*child_iter)->SetParentDihedralAngle((*child_iter)->GetParentDihedralAngle()+diff);
    //std::cout << "New Dihedral Angle: " << (*child_iter)->GetParentDihedralAngle() << "\n";
    child_iter++;
  }
  //std::cout << "\n"; Rotated:

}
