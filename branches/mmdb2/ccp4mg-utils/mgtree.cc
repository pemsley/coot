/*
     util/mgtree.cc: CCP4MG Molecular Graphics Program
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


#include <iostream>
#include <iomanip>
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

  c << std::setprecision(6);
  c << std::fixed;
   double rtod = 180.0 / M_PI;
   int i;

   for(i=0;i<a.GetNumberOfChildren();i++)
      c << a.GetChildren()[i];
   c << "\n";
   for(i=0;i<a.FindDepth();i++){c << " ";}
   c << std::setw(3) << a.GetID() << ": [Parent " << a.GetParentID() << "] Length:" << a.GetParentDistance() << " Angle:" << a.GetParentBondAngle()*rtod << " Torsion:" << a.GetParentDihedralAngle()*rtod << "\n";
   for(i=0;i<a.GetNumberOfChildren();i++)
      c << *(a.GetChildren()[i]);

  c.unsetf(std::ios::fixed | std::ios::scientific);
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
  ClearCoords();
}

double FindAngle(const TreeVertex *c1, const TreeVertex *c2, const TreeVertex *c3, const std::vector<std::vector<int> > &angles, const std::vector<double> &bond_angles){
  int v1 = c1->GetID();
  int v2 = c2->GetID();
  int v3 = c3->GetID();
  int atom1;
  int atom2;
  int atom3;
  double angle = 0.0;
  double rtod = 180.0 / M_PI;
  bool have_angle = false;
  for(unsigned ii=0;ii<angles.size();ii++){
    atom1 = angles[ii][0];
    atom2 = angles[ii][1];
    atom3 = angles[ii][2];
    if((atom1==v1&&atom2==v2&&atom3==v3)||(atom1==v1&&atom2==v3&&atom3==v2)||(atom1==v2&&atom2==v1&&atom3==v3)||(atom1==v2&&atom2==v3&&atom3==v1)||(atom1==v3&&atom2==v1&&atom3==v2)||(atom1==v3&&atom2==v2&&atom3==v1)){
      angle = bond_angles[ii]/rtod;
      //std::cout << "Angle: " << atom1 << " " << atom2 << " " << atom3 << ": " << angle << "\n"; std::cout.flush();
      have_angle = true;
      break;
    }
  }
  //if(!have_angle) std::cout << "Angle: " << v1 << " " << v2 << " " << v2 << " not found\n"; std::cout.flush();
  return angle;
}

void SetKnownTorsionPositions(std::vector<TreeVertex*> &coords, const std::vector<std::vector<int> > &torsions){
  for(unsigned ii=0;ii<torsions.size();ii++){
    int atom1 = torsions[ii][0];
    int atom2 = torsions[ii][1];
    int atom3 = torsions[ii][2];
    int atom4 = torsions[ii][3];
    if(coords[atom1]->GetParentID()==atom2&&coords[atom1]->GetParent()->GetParentID()==atom3&&coords[atom1]->GetParent()->GetParent()->GetParentID()==atom4){
      //coords[atom1]->SetParentDihedralAngle(torsion_angles[ii]/rtod);
      if(coords[atom1]->GetCoord().length()<1e-5&&coords[atom1]->GetParent()->GetCoord().length()>1e-5&&coords[atom1]->GetParent()->GetParent()->GetCoord().length()>1e-5){ // ie, it is coincident with first atom which must* be wrong (*chemically, not topologically).
         //std::cout << "SetCoord(SKTP1) of " << atom1 << "\n";
         coords[atom1]->SetCoord(GetCartFrom3Carts(coords[atom1]->GetParent()->GetCoord(),
           coords[atom1]->GetParentDistance(),
           coords[atom1]->GetParent()->GetParent()->GetCoord(),
           coords[atom1]->GetParentBondAngle(),
           coords[atom1]->GetParent()->GetParent()->GetParent()->GetCoord(),
           coords[atom1]->GetParentDihedralAngle()
         ));
      }
    }else if(coords[atom4]->GetParentID()==atom3&&coords[atom4]->GetParent()->GetParentID()==atom2&&coords[atom4]->GetParent()->GetParent()->GetParentID()==atom1){
      //coords[atom4]->SetParentDihedralAngle(torsion_angles[ii]/rtod);
      if(coords[atom4]->GetCoord().length()<1e-5&&coords[atom4]->GetParent()->GetCoord().length()>1e-5&&coords[atom4]->GetParent()->GetParent()->GetCoord().length()>1e-5){ // ie, it is coincident with first atom which must* be wrong (*chemically, not topologically).
         //std::cout << "SetCoord(SKTP2) of " << atom4 << "\n";
         coords[atom4]->SetCoord(GetCartFrom3Carts(coords[atom4]->GetParent()->GetCoord(),
           coords[atom4]->GetParentDistance(),
           coords[atom4]->GetParent()->GetParent()->GetCoord(),
           coords[atom4]->GetParentBondAngle(),
           coords[atom4]->GetParent()->GetParent()->GetParent()->GetCoord(),
           coords[atom4]->GetParentDihedralAngle()
         ));
        //std::cout << coords[atom4]->GetCoord() << "\n";
      }
    }
  }

}

void Tree::SetBondsAnglesTorsions(const int &nvertices,
                                  const std::vector<std::vector<int> > &bonds, const std::vector<double> &bond_lengths,
                                  const std::vector<std::vector<int> > &angles, const std::vector<double> &bond_angles,
                                  const std::vector<std::vector<int> > &torsions, const std::vector<double> &torsion_angles,
				  const std::vector<std::vector<int> > &chirals_in){

  std::vector<std::vector<int> > chirals = chirals_in;

  start = 0;
  if(nvertices<1) return;
  ClearCoords();
  extra_bonded_pairs.clear();

  std::vector<std::vector<int> > conn_lists(nvertices);
  for(unsigned ii=0;ii<bonds.size();ii++){
    conn_lists[bonds[ii][0]].push_back(bonds[ii][1]);
    conn_lists[bonds[ii][1]].push_back(bonds[ii][0]);
  }

  for(unsigned ii=0;ii<torsions.size();ii++){
    int atom3 = torsions[ii][2];
    int atom4 = torsions[ii][3];
    if(conn_lists[atom3].size()>1){
      std::vector<int>new_conn_lists;
      std::vector<int>rest_conn_lists;
      for(unsigned jj=0;jj<conn_lists[atom3].size();jj++){
        if(conn_lists[atom3][jj]==atom4)
           new_conn_lists.push_back(atom4);
        else
           rest_conn_lists.push_back(conn_lists[atom3][jj]);
      }
      new_conn_lists.insert(new_conn_lists.end(),rest_conn_lists.begin(),rest_conn_lists.end());
      conn_lists[atom3]=new_conn_lists;
    }
  }

  connectivity = conn_lists;

  /*
  for(unsigned ii=0;ii<connectivity.size();ii++){
    std::cout << ii << ": ";
    for(unsigned jj=0;jj<connectivity[ii].size();jj++){
      std::cout << connectivity[ii][jj] << " ";
    } std::cout << "\n";
  }
  */

  for(unsigned ii=0;ii<connectivity.size();ii++){
    for(unsigned jj=0;jj<connectivity[ii].size();jj++){
      if(std::find(connectivity[connectivity[ii][jj]].begin(),connectivity[connectivity[ii][jj]].end(),(int)ii)!=connectivity[connectivity[ii][jj]].end()){
      }else if(connectivity[ii][jj]!=(int)ii){
        connectivity[connectivity[ii][jj]].push_back(ii);
      }
    }
  }

  /*
  for(unsigned ii=0;ii<connectivity.size();ii++){
    std::cout << ii << ": ";
    for(unsigned jj=0;jj<connectivity[ii].size();jj++){
      std::cout << connectivity[ii][jj] << " ";
    } std::cout << "\n";
  }
  */

  std::vector<int>::iterator k;
  for(int i=0;i<nvertices;i++){
    coords.push_back(new TreeVertex());
    coords[i]->SetID(i);
    coords[i]->SetCoord(Cartesian(0,0,0));
  }
  
  scanned.clear();
  scanned.resize(coords.size());
  TreeVertex *coord = coords[0];
  RecurseCalculateTree(coord);
  for(unsigned i=0;i<coords.size();i++)
    if(!scanned[i]){
      coord = coords[i];
      RecurseCalculateTreeWithLengthsAndAngles(coord);
    }

  // Now we do lengths, distances and angles ....
  for(unsigned ii=0;ii<bonds.size();ii++){
    int atom1 = bonds[ii][0];
    int atom2 = bonds[ii][1];
    if(coords[atom1]->GetParentID()==atom2)
      coords[atom1]->SetParentDistance(bond_lengths[ii]);
    else
      coords[atom2]->SetParentDistance(bond_lengths[ii]);
  }

  // Well, bonds were easy and so are bond angles.
  double rtod = 180.0 / M_PI;
  for(unsigned ii=0;ii<angles.size();ii++){
    int atom1 = angles[ii][0];
    int atom2 = angles[ii][1];
    int atom3 = angles[ii][2];
    //std::cout << "Angle (" << atom1 << ")-(" << atom2 << ")-(" << atom3 << ") " << bond_angles[ii] << "\n";
    if(coords[atom3]->GetParentID()==atom2&&coords[atom3]->GetParent()->GetParentID()==atom1){
        coords[atom3]->SetParentBondAngle(bond_angles[ii]/rtod); // 1->2->3
    } else if(coords[atom1]->GetParentID()==atom2&&coords[atom1]->GetParent()->GetParentID()==atom3){
        coords[atom1]->SetParentBondAngle(bond_angles[ii]/rtod); // 3->2->1
    } else if(coords[atom2]->GetParentID()==atom1&&coords[atom2]->GetParent()->GetParentID()==atom3){
        coords[atom2]->SetParentBondAngle(bond_angles[ii]/rtod); // 3->1->2
    } else if(coords[atom3]->GetParentID()==atom1&&coords[atom3]->GetParent()->GetParentID()==atom2){
        coords[atom3]->SetParentBondAngle(bond_angles[ii]/rtod); // 2->1->3
    } else if(coords[atom1]->GetParentID()==atom3&&coords[atom1]->GetParent()->GetParentID()==atom2){
        coords[atom1]->SetParentBondAngle(bond_angles[ii]/rtod); // 2->3->1
    } else if(coords[atom2]->GetParentID()==atom3&&coords[atom2]->GetParent()->GetParentID()==atom1){
        coords[atom2]->SetParentBondAngle(bond_angles[ii]/rtod); // 1->3->2
    } else {
      //std::cout << "Ignore (" << atom1 << ")-(" << atom2 << ")-(" << atom3 << ")\n";
    }
  }

  for(unsigned i=0;i<coords.size();i++)
    coords[i]->SetParentDihedralAngle(-99999); // We can check later if a torsion is set using this absurd value.

  // Set up dummy positions (in the case we have top level node with 2 children. 
  if(coords[0]->GetNumberOfChildren()>0){
    Cartesian Dummy(0,-1,0);
    Cartesian Dummy2(-1,-2,0);
    coords[0]->SetDummy(Dummy,Dummy2);
    coords[0]->SetParentDistance(1);
    coords[0]->SetParentBondAngle(135/rtod);
    coords[0]->GetChild(0)->SetCoord(Cartesian(-coords[0]->GetChild(0)->GetParentDistance(),0,0)); // Move first child along x-axis.
    coords[0]->GetChild(0)->SetParentBondAngle(M_PI/2);
    coords[0]->GetChild(0)->SetParentDihedralAngle(0);
    // Replace next 17 lines with FindAngle(...)
    if(coords[0]->GetNumberOfChildren()>1){
    //std::cout << "FindAngle 1\n";
      /*
      double angle = FindAngle(coords[0],coords[0]->GetChild(0),coords[0]->GetChild(1),angles,bond_angles);
      coords[0]->GetChild(1)->SetCoord(Cartesian(coords[0]->GetChild(1)->GetParentDistance()*cos(M_PI-angle),coords[0]->GetChild(1)->GetParentDistance()*sin(M_PI-angle),0)); // Now position second child
      coords[0]->GetChild(1)->SetParentBondAngle(Angle(coords[0]->GetChild(1)->GetCoord(),coords[0]->GetCoord(),Dummy));
      coords[0]->GetChild(1)->SetParentDihedralAngle(DihedralAngle(coords[0]->GetChild(1)->GetCoord(),coords[0]->GetCoord(),Dummy,Dummy2));
     */
      for(int j=1;j<2;j++){
        //std::cout << "Now use black magic to determine torsion for " << coords[0]->GetChild(j)->GetID() << "\n";
        double angle1 = FindAngle(coords[0]->GetChild(j),coords[0],coords[0]->GetChild(j-1),angles,bond_angles);
        double angle2 = M_PI/2;
        //std::cout << angle1*rtod << "\n";
        coords[0]->GetChild(j)->SetCoord(GetCartFrom3Carts(coords[0]->GetCoord(),
        coords[0]->GetChild(j)->GetParentDistance(),
        coords[0]->GetChild(j-1)->GetCoord(),
        angle1,
        coords[0]->Dummy,
        angle2,
        1));
        coords[0]->GetChild(j)->SetParentBondAngle(Angle(coords[0]->GetChild(j)->GetCoord(),coords[0]->GetCoord(),Dummy));
        coords[0]->GetChild(j)->SetParentDihedralAngle(DihedralAngle(coords[0]->GetChild(j)->GetCoord(),coords[0]->GetCoord(),coords[0]->Dummy,coords[0]->Dummy2));
      }
      for(int j=2;j<coords[0]->GetNumberOfChildren();j++){
        //std::cout << "Now use black magic to determine torsion for " << coords[0]->GetChild(j)->GetID() << "\n";
        double angle1 = FindAngle(coords[0]->GetChild(j),coords[0],coords[0]->GetChild(j-1),angles,bond_angles);
        double angle2 = FindAngle(coords[0]->GetChild(j),coords[0],coords[0]->GetChild(0),angles,bond_angles);
        //std::cout << angle1*rtod << " " << angle2*rtod << "\n";
        coords[0]->GetChild(j)->SetCoord(GetCartFrom3Carts(coords[0]->GetCoord(),
        coords[0]->GetChild(j)->GetParentDistance(),
        coords[0]->GetChild(j-1)->GetCoord(),
        angle1,
        coords[0]->GetChild(0)->GetCoord(),
        angle2,
        1));
        coords[0]->GetChild(j)->SetParentBondAngle(Angle(coords[0]->GetChild(j)->GetCoord(),coords[0]->GetCoord(),Dummy));
        coords[0]->GetChild(j)->SetParentDihedralAngle(DihedralAngle(coords[0]->GetChild(j)->GetCoord(),coords[0]->GetCoord(),coords[0]->Dummy,coords[0]->Dummy2));
      }
    }
  }

  // All bond angles should be now OK. 
  // Next the torsions - tricky as we aren't given most of them, the info is in 
  // the bond angles which we haven't used yet. We have to work out torsions somehow.
  std::vector<int> scanned_torsions(torsions.size());
  for(unsigned ii=0;ii<torsions.size();ii++){
    int atom1 = torsions[ii][0];
    int atom2 = torsions[ii][1];
    int atom3 = torsions[ii][2];
    int atom4 = torsions[ii][3];
    //std::cout << "Torsion: (" << atom1 << ")-(" << atom2 << ")-(" << atom3 << ")-(" << atom4 << "), " << torsion_angles[ii] << "\n";
    if(coords[atom1]->GetParentID()==atom2&&coords[atom1]->GetParent()->GetParentID()==atom3&&coords[atom1]->GetParent()->GetParent()->GetParentID()==atom4){
      coords[atom1]->SetParentDihedralAngle(torsion_angles[ii]/rtod);
      scanned_torsions[ii] = 1;
    }else if(coords[atom4]->GetParentID()==atom3&&coords[atom4]->GetParent()->GetParentID()==atom2&&coords[atom4]->GetParent()->GetParent()->GetParentID()==atom1){
      coords[atom4]->SetParentDihedralAngle(torsion_angles[ii]/rtod);
      scanned_torsions[ii] = 1;
    }else{
      //std::cout << "Torsion not found\n";
    }
  }

  //std::cout << "Dummy1: " << coords[0]->Dummy << "\n";
  //std::cout << "Dummy2: " << coords[0]->Dummy2 << "\n";
  //std::cout << "coords[0]: " << coords[0]->GetCoord() << "\n";
  //std::cout << "coords[0] first child: " << coords[0]->GetChild(0)->GetCoord() << "\n";
  //std::cout << "coords[0] second child: " << coords[0]->GetChild(1)->GetCoord() << "\n";

  coords[0]->SetParentDihedralAngle(0); // By definition ??
  //for(int i=0;i<coords.size();i++){
    //if(coords[i]->FindDepth()>0&&fabs(coords[i]->GetParentDihedralAngle()+99999)<1e-5){
       //std::cout << "Set torsion for " << i << "\n";
    //}
  //}

  bool have_any_grandchild = false;
  // We have to set torsions of first lot of grandchildren, by some magic.
  for(int i=0;i<coords[0]->GetNumberOfChildren();i++){
    if(coords[0]->GetChild(i)->GetNumberOfChildren()>0){
      bool have_this_a_grandchild = false;
      for(int j=0;j<coords[0]->GetChild(i)->GetNumberOfChildren();j++){
        if(fabs(coords[0]->GetChild(i)->GetChild(j)->GetParentDihedralAngle()+99999)>1e-5){
          have_this_a_grandchild = true;
          have_any_grandchild = true;
        }
      }
      if(!have_this_a_grandchild){
        if(!have_any_grandchild){
          //std::cout << "We have to set to 180 (eg) for first child of child " << i << "\n";
          coords[0]->GetChild(i)->GetChild(0)->SetParentDihedralAngle(M_PI);
          // And now set its Cartesian based on this assumed torsion.
          //std::cout << "SetCoord(GC1) of " << coords[0]->GetChild(i)->GetChild(0)->GetID() << "\n";
          coords[0]->GetChild(i)->GetChild(0)->SetCoord(GetCartFrom3Carts(
          coords[0]->GetChild(i)->GetCoord(), // Atom 1
          coords[0]->GetChild(i)->GetChild(0)->GetParentDistance(), // blength
          coords[0]->GetCoord(), // Atom 2
          coords[0]->GetChild(i)->GetChild(0)->GetParentBondAngle(), // angle
          coords[0]->Dummy, // Atom 3
          M_PI // angle2 (torsion)
          ));
          //std::cout << coords[0]->GetChild(i)->GetID() << " has " << coords[0]->GetChild(i)->GetNumberOfChildren() << " children\n";

          //for(int j=1;j<2;j++){
          for(int j=1;j<coords[0]->GetChild(i)->GetNumberOfChildren();j++){
            //std::cout << "Now use black magic to determine torsion for " << i << ", " << j << "\n";
            double angle1 = FindAngle(coords[0]->GetChild(i)->GetChild(j),coords[0]->GetChild(i),coords[0]->GetChild(i)->GetChild(j-1),angles,bond_angles);
            double angle2 = FindAngle(coords[0]->GetChild(i)->GetChild(j),coords[0]->GetChild(i),coords[0],angles,bond_angles);
            //double angle3 = FindAngle(coords[0]->GetChild(i)->GetChild(j),coords[0]->GetChild(i),coords[0]->GetChild(i)->GetChild(0),angles,bond_angles);
            //std::cout << "angle3 " << angle3 << "\n";
            //std::cout << "SetCoord(GC2) of " << coords[0]->GetChild(i)->GetChild(j)->GetID() << "\n";
            //std::cout << "angle1 (" << coords[0]->GetChild(i)->GetChild(j)->GetID() << " " << coords[0]->GetChild(i)->GetID() << " " << coords[0]->GetChild(i)->GetChild(j-1)->GetID() << "), " << angle1*rtod << "\n";
            //std::cout << "angle2 (" << coords[0]->GetChild(i)->GetChild(j)->GetID() << " " << coords[0]->GetChild(i)->GetID() << " " << coords[0]->GetID() << "), " << angle2*rtod << "\n";
            coords[0]->GetChild(i)->GetChild(j)->SetCoord(GetCartFrom3Carts(coords[0]->GetChild(i)->GetCoord(),
            coords[0]->GetChild(i)->GetChild(j)->GetParentDistance(),
            coords[0]->GetChild(i)->GetChild(j-1)->GetCoord(),
            angle1,
            coords[0]->GetCoord(),
            angle2,
            1));
            //std::cout << "Distance: " << LineLength(coords[0]->GetChild(i)->GetChild(j)->GetCoord(),coords[0]->GetChild(i)->GetCoord()) << "\n";
            //std::cout << "Should be: " << coords[0]->GetChild(i)->GetChild(j)->GetParentDistance() << "\n";
            //std::cout << "Torsion: " << DihedralAngle(coords[0]->GetChild(i)->GetChild(j)->GetCoord(),coords[0]->GetChild(i)->GetCoord(),coords[0]->GetCoord(),coords[0]->Dummy)*rtod << "\n";
            coords[0]->GetChild(i)->GetChild(j)->SetParentDihedralAngle(DihedralAngle(coords[0]->GetChild(i)->GetChild(j)->GetCoord(),coords[0]->GetChild(i)->GetCoord(),coords[0]->GetCoord(),coords[0]->Dummy));
          }
          /* Why doesn't this work?
          for(int j=2;j<coords[0]->GetChild(i)->GetNumberOfChildren();j++){
            //std::cout << "Now use black magic to determine torsion for " << i << ", " << j << "\n";
            double angle1 = FindAngle(coords[0]->GetChild(i)->GetChild(j),coords[0]->GetChild(i),coords[0]->GetChild(i)->GetChild(j-1),angles,bond_angles);
            double angle2 = FindAngle(coords[0]->GetChild(i)->GetChild(j),coords[0]->GetChild(i),coords[0]->GetChild(i)->GetChild(0),angles,bond_angles);
            coords[0]->GetChild(i)->GetChild(j)->SetCoord(GetCartFrom3Carts(coords[0]->GetChild(i)->GetCoord(),
            coords[0]->GetChild(i)->GetChild(j)->GetParentDistance(),
            coords[0]->GetChild(i)->GetChild(j-1)->GetCoord(),
            angle1,
            coords[0]->GetChild(i)->GetChild(0)->GetCoord(),
            angle2,
            1));
            coords[0]->GetChild(i)->GetChild(j)->SetParentDihedralAngle(DihedralAngle(coords[0]->GetChild(i)->GetChild(j)->GetCoord(),coords[0]->GetChild(i)->GetCoord(),coords[0]->GetCoord(),coords[0]->Dummy));
          }
          */
          SetKnownTorsionPositions(coords,torsions);
        } else {
          // We'll leave this *for now* we generall have N as first atom, with perhaps H, then CA.
          // So this will be important but not for our fisrt Alanine case.
          //std::cout << "Hmm, second case, more than one child has children " << i << "\n";
        }
      } // else { We already have a grandchild dihedral set, we must use it somehow ... }
    }
  }

  int max_depth = FindMaxDepth();
  // And now we have to up the tree, working out torsions of greatgrandchildren etc.
  for(int depth=3;depth<=max_depth;depth++){
  for(unsigned i=0;i<coords.size();i++){
    if(coords[i]->FindDepth()==depth&&fabs(coords[i]->GetParentDihedralAngle()+99999)<1e-5){
    //std::cout << i << "\n";
      //std::cout << coords[i]->GetParentDihedralAngle() << "\n";
      int prev_set_tor = -1;
      int prev_set_tor_num = 0;
      for(int j=0;j<coords[i]->GetParent()->GetNumberOfChildren();j++){
        if(fabs(coords[i]->GetParent()->GetChild(j)->GetParentDihedralAngle()+99999)>1e-5){
           prev_set_tor = coords[i]->GetParent()->GetChild(j)->GetID();
           prev_set_tor_num = j;
           break;
        }
      }
      if(prev_set_tor>-1){
        //std::cout << prev_set_tor << " has a torsion we can use\n"; std::cout.flush(); // But we only want to do this once. Well, we will when done.
        //std::cout << "SetCoord(1) of " << coords[prev_set_tor]->GetID() << "\n";
        coords[prev_set_tor]->SetCoord(GetCartFrom3Carts(
          coords[prev_set_tor]->GetParent()->GetCoord(), // Atom 1
          coords[prev_set_tor]->GetParentDistance(), // blength
          coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(), // Atom 2
          coords[prev_set_tor]->GetParentBondAngle(), // angle
          coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord(), // Atom 3
          coords[prev_set_tor]->GetParentDihedralAngle()
          ));
          SetKnownTorsionPositions(coords,torsions);
          //std::cout << "Torsion: " << DihedralAngle(coords[prev_set_tor]->GetCoord(),coords[prev_set_tor]->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord())*rtod << "\n";
          //std::cout << "Should be: " << coords[prev_set_tor]->GetParentDihedralAngle()*rtod << "\n";
      } else {
        prev_set_tor = coords[i]->GetParent()->GetChild(0)->GetID();
        int prev_set_tor_parent = coords[i]->GetParent()->GetID();
        //std::cout << "SetCoord(2a) of " << coords[prev_set_tor]->GetID() << " with parent " << prev_set_tor_parent << "\n";

        bool ok=false;
        int this_angle = -1;
        int this_atom = -1;
        unsigned ii;
        for(ii=0;ii<angles.size();ii++){
          int atom1 = angles[ii][0];
          int atom2 = angles[ii][1];
          int atom3 = angles[ii][2];
          if((atom3==prev_set_tor&&atom2==prev_set_tor_parent)){
            bool is_a_child = false;
            for(int kk=0;kk<coords[prev_set_tor]->GetNumberOfChildren();kk++){
              if(coords[prev_set_tor]->GetChild(kk)->GetID()==atom1) {
                is_a_child = true;
              }
            }
            bool is_grandparent = false;
            if((coords[prev_set_tor]->GetParent()->GetParentID()==atom1)) {
              is_grandparent = true;
            }
            bool isset = true;
            if((coords[0]->GetCoord()-coords[atom1]->GetCoord()).length()<1e-5&&atom1!=0) isset = false;
            if(isset&&!is_a_child&&!is_grandparent){
              //std::cout << "Consider non-tree a(" << atom1 << ")-(" << atom2 << ")-(" << atom3 << ")\n";
              this_atom = atom1; ok = true;
            }
          }else if((atom1==prev_set_tor&&atom2==prev_set_tor_parent)){
            bool is_a_child = false;
            for(int kk=0;kk<coords[prev_set_tor]->GetNumberOfChildren();kk++)
              if(coords[prev_set_tor]->GetChild(kk)->GetID()==atom3){
                is_a_child = true;
              }
            bool is_grandparent = false;
            if((coords[prev_set_tor]->GetParent()->GetParentID()==atom3)){
              is_grandparent = true;
            }
            bool isset = true;
            if((coords[0]->GetCoord()-coords[atom3]->GetCoord()).length()<1e-5&&atom3!=0) isset = false;
            if(isset&&!is_a_child&&!is_grandparent){
              //std::cout << "Consider non-tree a(" << atom1 << ")-(" << atom2 << ")-(" << atom3 << ")\n";
              this_atom = atom3; ok = true;
            }
          }
          if(ok) {this_angle = ii;/* std::cout << "this_atom:" << this_atom << "\n"*/; break; }
        }
        if(this_atom!=0&&coords[prev_set_tor_parent]->GetParentID()!=-1&&coords[prev_set_tor_parent]->GetParent()->GetParentID()!=-1){
          bool have_off_tree_torsion = false;
          //std::cout << "We have hope\n";
          for(unsigned jj=0;jj<torsions.size();jj++){
            int atom1 = torsions[jj][0];
            int atom2 = torsions[jj][1];
            int atom3 = torsions[jj][2];
            int atom4 = torsions[jj][3];
            if((!scanned_torsions[jj])&&((atom1==this_atom&&atom2==prev_set_tor_parent&&atom3==coords[prev_set_tor_parent]->GetParentID()&&atom4==coords[prev_set_tor_parent]->GetParent()->GetParentID())||(atom4==this_atom&&atom3==prev_set_tor_parent&&atom2==coords[prev_set_tor_parent]->GetParentID()&&atom1==coords[prev_set_tor_parent]->GetParent()->GetParentID()))){
              double thisparentdist=0.0;
              for(unsigned kk=0;kk<bonds.size();kk++){
                int atom1 = bonds[kk][0];
                int atom2 = bonds[kk][1];
                if((atom1==this_atom&&atom2==prev_set_tor_parent)||(atom2==this_atom&&atom1==prev_set_tor_parent)){
                  thisparentdist = bond_lengths[kk];
                  //std::cout << "Use distance: " << thisparentdist<< "\n";
                }
              }
              if((coords[this_atom]->GetCoord()-coords[0]->GetCoord()).length()<1e-5){
              double thisparentparentangle = FindAngle(coords[this_atom],coords[prev_set_tor]->GetParent(),coords[prev_set_tor]->GetParent()->GetParent(),angles,bond_angles);
              coords[this_atom]->SetCoord(GetCartFrom3Carts(coords[prev_set_tor]->GetParent()->GetCoord(),
              thisparentdist,
              coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),
              thisparentparentangle,
              coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord(),
              torsion_angles[jj]/rtod));
              }
              scanned_torsions[jj] = 1;
              have_off_tree_torsion = true;
              break;
            }
          }
          if(!have_off_tree_torsion){
            ok = false;
          }
        }
        /*
        for(unsigned jj=0;jj<torsions.size();jj++){
          if(!scanned_torsions[jj])std::cout << "Torsion remaining: (" << torsions[jj][0] << ")-(" << torsions[jj][1] << ")-(" << torsions[jj][2] << ")-(" << torsions[jj][3] << "), " << torsion_angles[jj] << "\n";
        }
        */
        if(ok){ // Now the real work
          //std::cout << "Using horrible stuff to set coord of " << prev_set_tor << "\n";
          //std::cout << "this_angle " << this_angle<< "\n";
          double angle1 = coords[prev_set_tor]->GetParentBondAngle();
          double angle2 = bond_angles[this_angle]/rtod;
          //std::cout << angle1*rtod << " " << angle2*rtod << "\n";
            coords[prev_set_tor]->SetCoord(GetCartFrom3Carts(coords[prev_set_tor]->GetParent()->GetCoord(),
            coords[prev_set_tor]->GetParentDistance(),
            coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),
            angle1,
            coords[this_atom]->GetCoord(),
            angle2,
            -1)); // This should depend on chirality. This however gets things better for Proline.
          coords[prev_set_tor]->SetParentDihedralAngle(DihedralAngle(coords[prev_set_tor]->GetCoord(),coords[prev_set_tor]->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord()));
          int rem_tor = 0;
          for(unsigned jj=0;jj<torsions.size();jj++){
            if(!scanned_torsions[jj]) rem_tor++;
          }
          
        }else{
          // We'll see if we are in a ring and there is a torsion that goes "other way".
          bool set_this_prev_by_another_torsion=false;
          for(unsigned jj=0;jj<torsions.size();jj++){
            if(!scanned_torsions[jj]) {
             int atom1 = torsions[jj][0];
             int atom2 = torsions[jj][1];
             int atom3 = torsions[jj][2];
             int atom4 = torsions[jj][3];
             if(atom1==coords[prev_set_tor]->GetID()){
               if(coords[atom4]->GetCoord().length()>1e-5&&coords[atom3]->GetCoord().length()>1e-5){ // Position is set
                 //if(coords[atom2]->GetParentID()==prev_set_tor)
                   //std::cout << "Possible Torsion(1): (" << atom1 << ")-(" << atom2 << ")-(" << atom3 << ")-(" << atom4 << "), " << torsion_angles[jj] << "\n";
                   // Hope this is never true.
               }
             }else if(atom4==coords[prev_set_tor]->GetID()){
               if(coords[atom1]->GetCoord().length()>1e-5&&coords[atom2]->GetCoord().length()>1e-5){ // Position is set
                 if(coords[atom3]->GetParentID()==prev_set_tor){
                   //std::cout << "Possible Torsion(2): (" << atom1 << ")-(" << atom2 << ")-(" << atom3 << ")-(" << atom4 << "), " << torsion_angles[jj] << "\n";
                   if(coords[atom3]->GetCoord().length()<1e-5){
                     //std::cout << "First we have to set atom 3\n";
                     if(coords[atom2]->GetNumberOfChildren()>0){
                       //std::cout << "Good news: 2 has a child\n";
                       // Is parent chiral?
                       unsigned ichiral;
                       for(ichiral=0;ichiral<chirals.size();ichiral++){
                         //std::cout << chirals[ichiral][0] <<  " " << atom2 << "\n";
                           if(chirals[ichiral][0]==atom2)
                             break;
                       }
                       if(ichiral<chirals.size()){
                         //std::cout << "Centre " << chirals[ichiral][0] << " " << atom2 << " is chiral, but we ignore that.\n";
		       }
                       double angle1 = FindAngle(coords[atom3],coords[atom2],coords[atom2]->GetParent(),angles,bond_angles);
                       double angle2 = FindAngle(coords[atom3],coords[atom2],coords[atom2]->GetChild(0),angles,bond_angles);
                       double dist = 0.0;
		       bool positive_chiral = true;
                       for(unsigned kk=0;kk<bonds.size();kk++){
                         int bond_atom1 = bonds[kk][0];
                         int bond_atom2 = bonds[kk][1];
                         if((bond_atom1==atom3&&bond_atom2==atom2)||(bond_atom2==atom3&&bond_atom1==atom2)){
                           dist = bond_lengths[kk];
                         }
                       }
                       //std::cout << angle1*rtod << " " << angle2*rtod << ", " << dist << "\n";
                       coords[atom3]->SetCoord(GetCartFrom3Carts(coords[atom2]->GetCoord(),
                       dist,
                       coords[atom2]->GetParent()->GetCoord(),
                       angle1,
                       coords[atom2]->GetChild(0)->GetCoord(),
                       angle2,
                       1));
                       //std::cout << "Set coord to " << coords[atom3]->GetCoord() << "\n";
                       double another_dist = 0.0;
                       for(unsigned kk=0;kk<bonds.size();kk++){
                         int bond_atom1 = bonds[kk][0];
                         int bond_atom2 = bonds[kk][1];
                         if((bond_atom1==atom4&&bond_atom2==atom3)||(bond_atom2==atom4&&bond_atom1==atom3)){
                           another_dist = bond_lengths[kk];
                         }
                       }
                       double angle3 = FindAngle(coords[atom4],coords[atom3],coords[atom2],angles,bond_angles);
                       //std::cout << "SetCoord(2b) of " << coords[prev_set_tor]->GetID() << "\n";
                       coords[prev_set_tor]->SetCoord(GetCartFrom3Carts(
                         coords[atom3]->GetCoord(), // Atom 1
                         another_dist,
                         coords[atom2]->GetCoord(), // Atom 2
                         angle3,
                         coords[atom1]->GetCoord(), // Atom 3
                         -torsion_angles[jj]/rtod
                         ));
		       scanned_torsions[jj] = 1;
                       //std::cout << "Unknown coord set to " << coords[atom4]->GetCoord() << "\n";
		       //std::cout << "Check torsion:\n";
		       //std::cout << DihedralAngle(coords[prev_set_tor]->GetCoord(),coords[atom3]->GetCoord(),coords[atom2]->GetCoord(),coords[atom1]->GetCoord())*rtod << "\n";
                       //std::cout << prev_set_tor << " " << atom4 << " " << coords[prev_set_tor]->GetID() << "\n";
                       //std::cout << "Ok, the parent: " << coords[prev_set_tor]->GetParent()->GetCoord() << " " << GetCartesian(coords[prev_set_tor]->GetParent()) << "\n";
                       coords[prev_set_tor]->SetParentDihedralAngle(DihedralAngle(coords[prev_set_tor]->GetCoord(),coords[prev_set_tor]->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord()));
		       if((GetCartesian(coords[atom4])-coords[atom4]->GetCoord()).length()>1e-1){

                         //std::cout << "this is a bugger " << (GetCartesian(coords[atom4])-coords[atom4]->GetCoord()).length() << ", try the other way:\n"; std::cout.flush();
                         coords[atom3]->SetCoord(GetCartFrom3Carts(coords[atom2]->GetCoord(),
                         dist,
                         coords[atom2]->GetParent()->GetCoord(),
                         angle1,
                         coords[atom2]->GetChild(0)->GetCoord(),
                         angle2,
                         -1));
                         //std::cout << "Set coord to " << coords[atom3]->GetCoord() << "\n";
                         //std::cout << "SetCoord(2b2) of " << coords[prev_set_tor]->GetID() << "\n";
                         coords[prev_set_tor]->SetCoord(GetCartFrom3Carts(
                           coords[atom3]->GetCoord(), // Atom 1
                           another_dist,
                           coords[atom2]->GetCoord(), // Atom 2
                           angle3,
                           coords[atom1]->GetCoord(), // Atom 3
                           -torsion_angles[jj]/rtod
                           ));
                         //std::cout << "Unknown coord set to " << coords[atom4]->GetCoord() << "\n";
                         coords[prev_set_tor]->SetParentDihedralAngle(DihedralAngle(coords[prev_set_tor]->GetCoord(),coords[prev_set_tor]->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord()));
		         positive_chiral = false;

		       }
		       if((GetCartesian(coords[atom4])-coords[atom4]->GetCoord()).length()>1e-1){
			  //std::cout << "This is really unfortunate!\n";
		       }
		       for(int ichild=0;ichild<coords[atom2]->GetNumberOfChildren();ichild++){
			  if(coords[atom2]->GetChild(ichild)->GetID()!=atom2){
			     double angle_got = Angle(coords[atom3]->GetCoord(),coords[atom2]->GetCoord(),coords[atom2]->GetChild(ichild)->GetCoord());
			     double angle_wanted = FindAngle(coords[atom3],coords[atom2],coords[atom2]->GetChild(ichild),angles,bond_angles);
			     if(fabs(angle_got-angle_wanted)>1e-3){
			       int chiral_direction = -1;
			       if(positive_chiral) chiral_direction = 1;
			       //std::cout << "Other child " << coords[atom2]->GetChild(ichild)->GetID() << " has dodgy angle\n"; 
                               coords[atom2]->GetChild(ichild)->SetCoord(GetCartFrom3Carts(coords[atom2]->GetCoord(),
                               coords[atom2]->GetChild(ichild)->GetParentDistance(),
                               coords[atom2]->GetParent()->GetCoord(),
                               coords[atom2]->GetChild(ichild)->GetParentBondAngle(),
                               coords[atom3]->GetCoord(),
                               angle_wanted,
                               -1));
			       int id = coords[atom2]->GetChild(ichild)->GetID();
                               coords[id]->SetParentDihedralAngle(DihedralAngle(coords[id]->GetCoord(),coords[id]->GetParent()->GetCoord(),coords[id]->GetParent()->GetParent()->GetCoord(),coords[id]->GetParent()->GetParent()->GetParent()->GetCoord()));
			     }
			  }
		       }

			 /* Is this ever necessary ?
			 // Now go backwards around the ring until we are one before the first child of "atom2".
                         int atom1_new = atom2;
                         int atom2_new = atom3;
                         int atom3_new = atom4;
                         int atom4_new = coords[atom4]->GetParent()->GetID();
			 bool have_torsions = true;
			 std::vector<int> set_positions;
			 while(atom4_new!=atom2&&have_torsions){
			   std::cout << "this old pos: " << coords[atom4_new]->GetCoord() << " " << GetCartesian(coords[atom4_new]) << "\n";
			   std::cout << atom4_new << "\n"; std::cout.flush();
			   unsigned kk;
                           for(kk=0;kk<torsions.size();kk++){
                             if(!scanned_torsions[kk]) {
                               if((atom1_new == torsions[kk][0]&& atom2_new == torsions[kk][1]&&
                                  atom3_new == torsions[kk][2]&& atom4_new == torsions[kk][3])||
			         (atom1_new == torsions[kk][3]&& atom2_new == torsions[kk][2]&&
                                  atom3_new == torsions[kk][1]&& atom4_new == torsions[kk][0])){
                                  scanned_torsions[kk] = 1;
                                  std::cout << "Possible Torsion(3): (" << atom1_new << ")-(" << atom2_new << ")-(" << atom3_new << ")-(" << atom4_new << "), " << torsion_angles[kk] << "\n"; std::cout.flush();
                                 coords[atom4_new]->SetCoord(GetCartFrom3Carts(
                                    coords[atom3_new]->GetCoord(),
                                    coords[atom3_new]->GetParentDistance(),
                                    coords[atom2_new]->GetCoord(),
                                    coords[atom2_new]->GetParentBondAngle(),
                                    coords[atom1_new]->GetCoord(),
                                    -torsion_angles[kk]/rtod));
				 set_positions.push_back(atom4_new);
				 atom1_new = atom2_new;
				 atom2_new = atom3_new;
				 atom3_new = atom4_new;
				 atom4_new = coords[atom4_new]->GetParent()->GetID();
			         break;
			       }
			     }
			   }
			   std::cout << kk << " " << torsions.size() << "\n";
			   if(kk==torsions.size()) have_torsions = false;
		         }
			 std::cout << "Now set dihedrals\n"; std::cout.flush();
			 while(set_positions.size()>0){
			   coords[set_positions.back()]->SetParentDihedralAngle(
			   DihedralAngle(coords[set_positions.back()]->GetCoord(),coords[set_positions.back()]->GetParent()->GetCoord(),coords[set_positions.back()]->GetParent()->GetParent()->GetCoord(),coords[set_positions.back()]->GetParent()->GetParent()->GetParent()->GetCoord()));
                           set_positions.pop_back();
			   std::cout << "this new pos: " << coords[atom4_new]->GetCoord() << " " << GetCartesian(coords[atom4_new]) << "\n";
			 }
                         coords[prev_set_tor]->SetParentDihedralAngle(DihedralAngle(coords[prev_set_tor]->GetCoord(),coords[prev_set_tor]->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord()));
			 */

                       //std::cout << "Unknown coord gotten as " << GetCartesian(coords[atom4]) << "\n";
                       coords[atom3]->SetParentDihedralAngle(DihedralAngle(coords[atom3]->GetCoord(),coords[atom3]->GetParent()->GetCoord(),coords[atom3]->GetParent()->GetParent()->GetCoord(),coords[atom3]->GetParent()->GetParent()->GetParent()->GetCoord()));
                       //std::cout << "coord gotten as " << GetCartesian(coords[atom3]) << "\n";
                       //std::cout << (GetCartesian(coords[atom3])-coords[atom3]->GetCoord()).length() << "\n";
                       //std::cout << (GetCartesian(coords[atom4])-coords[atom4]->GetCoord()).length() << "\n";
                       //std::cout << "Attempt to set torsion to " << coords[atom3]->GetParentDihedralAngle()*rtod << "\n";
                       set_this_prev_by_another_torsion=true; break;
                     }
                   }
                 }
               }
             }
            }
          } 
          for(unsigned jj=0;jj<torsions.size();jj++){
            if(!scanned_torsions[jj]) {
              int atom1 = torsions[jj][0];
              int atom2 = torsions[jj][1];
              int atom3 = torsions[jj][2];
              int atom4 = torsions[jj][3];
	      if(atom4==prev_set_tor&&atom3==coords[prev_set_tor]->GetParent()->GetID()&&atom2==coords[prev_set_tor]->GetParent()->GetParent()->GetID()&&coords[atom1]->GetParentID()==atom2){
	        //std::cout << "Possible remaining torsion: " << atom1 << " " << atom2 << " " << atom3 << " " << atom4 << "\n";
                coords[prev_set_tor]->SetCoord(GetCartFrom3Carts(
                  coords[prev_set_tor]->GetParent()->GetCoord(), // Atom 1
                  coords[prev_set_tor]->GetParentDistance(), // blength
                  coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(), // Atom 2
                  coords[prev_set_tor]->GetParentBondAngle(), // angle
                  coords[atom1]->GetCoord(), // Atom 3
                  torsion_angles[jj]/rtod
                  ));
                coords[prev_set_tor]->SetParentDihedralAngle(DihedralAngle(coords[prev_set_tor]->GetCoord(),coords[prev_set_tor]->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(),coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord()));
                set_this_prev_by_another_torsion=true; break;
	      }
            } 
          } 
        if(!set_this_prev_by_another_torsion){ // New magic
          //std::cout << "SetCoord(2) of " << coords[prev_set_tor]->GetID() << "\n";
          //So this is guesswork. We don't really want to be here, this is only for amino acid Os at the moment.
          coords[prev_set_tor]->SetParentDihedralAngle(M_PI);
          coords[prev_set_tor]->SetCoord(GetCartFrom3Carts(
            coords[prev_set_tor]->GetParent()->GetCoord(), // Atom 1
            coords[prev_set_tor]->GetParentDistance(), // blength
            coords[prev_set_tor]->GetParent()->GetParent()->GetCoord(), // Atom 2
            coords[prev_set_tor]->GetParentBondAngle(), // angle
            coords[prev_set_tor]->GetParent()->GetParent()->GetParent()->GetCoord(), // Atom 3
            coords[prev_set_tor]->GetParentDihedralAngle()
            ));
          }
        }
          //SetKnownTorsionPositions(coords,torsions);
      }
      // Is parent chiral?
      unsigned ichiral;
      for(ichiral=0;ichiral<chirals.size();ichiral++){
        //std::cout << chirals[ichiral][0] <<  " " << coords[i]->GetParent()->GetID() << "\n";
          if(chirals[ichiral][0]==coords[i]->GetParent()->GetID())
            break;
      }
      if(ichiral<chirals.size()){
        for(int j=0;j<coords[i]->GetParent()->GetNumberOfChildren();j++){
            if(j!=prev_set_tor_num){
              //std::cout << "SetCoord(3c) of " << coords[i]->GetParent()->GetChild(j)->GetID() << "\n";
	      int chiral_direction = 1;
              //std::cout << "Centre " << chirals[ichiral][0] << " " << coords[i]->GetParent()->GetID() << " is chiral\n";
	      //std::cout << chirals[ichiral][1] << " " << chirals[ichiral][2] << " " << chirals[ichiral][3] << "\n";
	      //std::cout << coords[i]->GetParent()->GetParent()->GetID() << " " << coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetID() << " " << coords[i]->GetParent()->GetChild(j)->GetID() << "\n";
	      if(coords[i]->GetParent()->GetParent()->GetID()==chirals[ichiral][1]&&coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetID()==chirals[ichiral][2]&&chirals[ichiral][3]==coords[i]->GetParent()->GetChild(j)->GetID()){
		//std::cout << "case1\n";
       	        chiral_direction *=-1;
	      }
	      if(coords[i]->GetParent()->GetParent()->GetID()==chirals[ichiral][2]&&coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetID()==chirals[ichiral][3]&&chirals[ichiral][1]==coords[i]->GetParent()->GetChild(j)->GetID()){
		//std::cout << "case2\n";
       	        chiral_direction *=-1;
	      }
	      if(coords[i]->GetParent()->GetParent()->GetID()==chirals[ichiral][3]&&coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetID()==chirals[ichiral][1]&&chirals[ichiral][2]==coords[i]->GetParent()->GetChild(j)->GetID()){
		//std::cout << "case3\n";
       	        chiral_direction *=-1;
	      }
	      if(coords[i]->GetParent()->GetParent()->GetID()==chirals[ichiral][2]&&coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetID()==chirals[ichiral][1]&&chirals[ichiral][3]!=coords[i]->GetParent()->GetChild(j)->GetID()){
		//std::cout << "case4\n";
       	        chiral_direction *=-1;
	      }
	      if(coords[i]->GetParent()->GetParent()->GetID()==chirals[ichiral][3]&&coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetID()==chirals[ichiral][2]&&chirals[ichiral][1]!=coords[i]->GetParent()->GetChild(j)->GetID()){
		//std::cout << "case5\n";
       	        chiral_direction *=-1;
	      }
	      if(coords[i]->GetParent()->GetParent()->GetID()==chirals[ichiral][1]&&coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetID()==chirals[ichiral][3]&&chirals[ichiral][2]!=coords[i]->GetParent()->GetChild(j)->GetID()){
		//std::cout << "case6\n";
       	        chiral_direction *=-1;
	      }
	      if(chirals[ichiral][4]==-1)
       	        chiral_direction *=-1;
	      //std::cout << "Chiral direction: " << chiral_direction << "\n";
              double angle1 = FindAngle(coords[i]->GetParent()->GetChild(j),coords[i]->GetParent(),coords[i]->GetParent()->GetChild(prev_set_tor_num),angles,bond_angles);
              double angle2 = FindAngle(coords[i]->GetParent()->GetChild(j),coords[i]->GetParent(),coords[i]->GetParent()->GetParent(),angles,bond_angles);
              coords[i]->GetParent()->GetChild(j)->SetCoord(GetCartFrom3Carts(coords[i]->GetParent()->GetCoord(),
              coords[i]->GetParent()->GetChild(j)->GetParentDistance(),
              coords[i]->GetParent()->GetChild(prev_set_tor_num)->GetCoord(),
              angle1,
              coords[i]->GetParent()->GetParent()->GetCoord(),
              angle2,
              chiral_direction)); // Depends on chirality ??
              coords[i]->GetParent()->GetChild(j)->SetParentDihedralAngle(DihedralAngle(coords[i]->GetParent()->GetChild(j)->GetCoord(),coords[i]->GetParent()->GetCoord(),coords[i]->GetParent()->GetParent()->GetCoord(),coords[i]->GetParent()->GetParent()->GetParent()->GetCoord()));
            SetKnownTorsionPositions(coords,torsions);
	  }
	}
      } else {
        for(int j=prev_set_tor_num+1;j<coords[i]->GetParent()->GetNumberOfChildren();j++){
            //std::cout << "SetCoord(3) of " << coords[i]->GetParent()->GetChild(j)->GetID() << "\n";
              double angle1 = FindAngle(coords[i]->GetParent()->GetChild(j),coords[i]->GetParent(),coords[i]->GetParent()->GetChild(j-1),angles,bond_angles);
              double angle2 = FindAngle(coords[i]->GetParent()->GetChild(j),coords[i]->GetParent(),coords[i]->GetParent()->GetParent(),angles,bond_angles);
              coords[i]->GetParent()->GetChild(j)->SetCoord(GetCartFrom3Carts(coords[i]->GetParent()->GetCoord(),
              coords[i]->GetParent()->GetChild(j)->GetParentDistance(),
              coords[i]->GetParent()->GetChild(j-1)->GetCoord(),
              angle1,
              coords[i]->GetParent()->GetParent()->GetCoord(),
              angle2,
              1)); // Depends on chirality ??
              coords[i]->GetParent()->GetChild(j)->SetParentDihedralAngle(DihedralAngle(coords[i]->GetParent()->GetChild(j)->GetCoord(),coords[i]->GetParent()->GetCoord(),coords[i]->GetParent()->GetParent()->GetCoord(),coords[i]->GetParent()->GetParent()->GetParent()->GetCoord()));
	    }
            //std::cout << "Distance: " << LineLength(coords[i]->GetParent()->GetChild(j)->GetCoord(),coords[i]->GetParent()->GetCoord()) << "\n";
            //std::cout << "Should be: " << coords[i]->GetParent()->GetChild(j)->GetParentDistance() << "\n" << "\n";
        //}
        }
        SetKnownTorsionPositions(coords,torsions);
        for(int j=prev_set_tor_num-1;j>=0;j--){
        //if(coords[i]->GetParent()->GetChild(j)->GetID()!=prev_set_tor){
            //std::cout << "set torsion (-) for  " << coords[i]->GetParent()->GetChild(j)->GetID() << "\n";

            double angle1 = FindAngle(coords[i]->GetParent()->GetChild(j),coords[i]->GetParent(),coords[i]->GetParent()->GetChild(j+1),angles,bond_angles);
            double angle2 = FindAngle(coords[i]->GetParent()->GetChild(j),coords[i]->GetParent(),coords[i]->GetParent()->GetParent(),angles,bond_angles);
	    int chiral_direction = -1;
            //std::cout << "SetCoord(4) of " << coords[i]->GetParent()->GetChild(j)->GetID() << "\n";
            coords[i]->GetParent()->GetChild(j)->SetCoord(GetCartFrom3Carts(coords[i]->GetParent()->GetCoord(),
            coords[i]->GetParent()->GetChild(j)->GetParentDistance(),
            coords[i]->GetParent()->GetChild(j+1)->GetCoord(),
            angle1,
            coords[i]->GetParent()->GetParent()->GetCoord(),
            angle2,
            chiral_direction)); // Depends on chirality ??
            SetKnownTorsionPositions(coords,torsions);

            coords[i]->GetParent()->GetChild(j)->SetParentDihedralAngle(DihedralAngle(coords[i]->GetParent()->GetChild(j)->GetCoord(),coords[i]->GetParent()->GetCoord(),coords[i]->GetParent()->GetParent()->GetCoord(),coords[i]->GetParent()->GetParent()->GetParent()->GetCoord()));
            //std::cout << "Distance: " << LineLength(coords[i]->GetParent()->GetChild(j)->GetCoord(),coords[i]->GetParent()->GetCoord()) << "\n";
            //std::cout << "Should be: " << coords[i]->GetParent()->GetChild(j)->GetParentDistance() << "\n" << "\n";
        //}
        }
        SetKnownTorsionPositions(coords,torsions);
      }
    }
  }
  
  // A check of stored cordinates and ones we get from traversing tree.
  for(unsigned i=0;i<coords.size();i++){
     Cartesian calced = GetCartesian(coords[i]);
     Cartesian stored = coords[i]->GetCoord();
     if((calced-stored).length()>1e-2)
     std::cout << "Potential problem in Tree calculation[" << i << "]: " << calced << ", " << stored << " (" << (calced-stored).length() << ")\n";
  }

  // Check the bond lengths ...
  for(unsigned ii=0;ii<bonds.size();ii++){
    int atom1 = bonds[ii][0];
    int atom2 = bonds[ii][1];
    if(fabs(bond_lengths[ii]-(coords[atom1]->GetCoord()-coords[atom2]->GetCoord()).length())>1e-2){
      std::cout << "Bad bond length(" << atom1 << " " << atom2 << "): " << bond_lengths[ii] << " " << (coords[atom1]->GetCoord()-coords[atom2]->GetCoord()).length() << "\n";
    }
  }

  // ... and angles ...
  for(unsigned ii=0;ii<angles.size();ii++){
    int atom1 = angles[ii][0];
    int atom2 = angles[ii][1];
    int atom3 = angles[ii][2];
    if(fabs(bond_angles[ii]-Angle(coords[atom1]->GetCoord(),coords[atom2]->GetCoord(),coords[atom3]->GetCoord())*rtod)>1){
        std::cout << "Bad bond angle(" << atom1 << " " << atom2 << " " << atom3 << "): " << bond_angles[ii] << " " << Angle(coords[atom1]->GetCoord(),coords[atom2]->GetCoord(),coords[atom3]->GetCoord())*rtod << "\n";
    }
  }

  // ... and torsions ...
  for(unsigned ii=0;ii<torsions.size();ii++){
    int atom1 = torsions[ii][0];
    int atom2 = torsions[ii][1];
    int atom3 = torsions[ii][2];
    int atom4 = torsions[ii][3];
    if(fabs(torsion_angles[ii]- DihedralAngle(coords[atom1]->GetCoord(),coords[atom2]->GetCoord(),coords[atom3]->GetCoord(),coords[atom4]->GetCoord())*rtod)>1e-5&&fabs(torsion_angles[ii]- 180.0)>1){
      std::cout << "Bad torsion angle(" << atom1 << " " << atom2 << " " << atom3 << " " << atom4 << "): " << torsion_angles[ii] << " " << DihedralAngle(coords[atom1]->GetCoord(),coords[atom2]->GetCoord(),coords[atom3]->GetCoord(),coords[atom4]->GetCoord())*rtod << "\n";
    }
  }

  // Then extra bonded pair stuff. Probably either useless or broken.
  std::vector<std::pair<int,int> > bp;
  for(unsigned i=0;i<coords.size();i++){
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

}

void Tree::RecurseCalculateTreeWithLengthsAndAngles(TreeVertex* coord){

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
         RecurseCalculateTreeWithLengthsAndAngles(coordk);
       } else {  //added by Joel Bard
         //but coord may not have a parent so add it as
         //a child to coordk
         if( coord->GetParentID()==-1&&scanned[*k]!=1 ) {
           coord->SetParentID(coordk->GetID());
           coord->SetParent(coordk);
           coordk->AddChild(coord);
         }
       }
    }
    k++;
  } 
}


Tree::Tree(const int &nvertices,
           const std::vector<std::vector<int> > &bonds, const std::vector<double> &bond_lengths,
           const std::vector<std::vector<int> > &angles, const std::vector<double> &bond_angles,
           const std::vector<std::vector<int> > &torsions, const std::vector<double> &torsion_angles,
	   const std::vector<std::vector<int> > &chirals){
  SetBondsAnglesTorsions(nvertices,bonds,bond_lengths,angles,bond_angles,torsions,torsion_angles,chirals);
}

Tree::Tree(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists, const  std::vector<std::vector<Cartesian> > &ext_cartesians){
  std::vector<std::vector<int> > forced_connections;
  if(ext_cartesians.size()>0)
    SetCoords(SelAtoms_in, ioff, conn_lists, ext_cartesians,  forced_connections);
  else
    SetCoords(SelAtoms_in, ioff, conn_lists, forced_connections);
}

Tree::Tree(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists, const  std::vector<std::vector<Cartesian> > &ext_cartesians, const std::vector<std::vector<int> > &forced_connections){
  if(ext_cartesians.size()>0)
    SetCoords(SelAtoms_in, ioff, conn_lists, ext_cartesians, forced_connections);
  else
    SetCoords(SelAtoms_in, ioff, conn_lists, forced_connections);
}

void Tree::SetCoords(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists){
  std::vector<std::vector<Cartesian> > ext_carts(SelAtoms_in.size());
  std::vector<std::vector<int> > forced_connections;
  SetCoords(SelAtoms_in,ioff,conn_lists,ext_carts,forced_connections);
}

void Tree::SetCoords(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists, const std::vector<std::vector<int> > &forced_connections){
  std::vector<std::vector<Cartesian> > ext_carts(SelAtoms_in.size());
  SetCoords(SelAtoms_in,ioff,conn_lists,ext_carts,forced_connections);
}

void Tree::ForceEarlyConnection(int parent,int child){
	std::vector<int> parent_conns = connectivity[parent];
	std::vector<int> child_conns = connectivity[child];

	bool isGood = false;
	for(unsigned i=0;i<parent_conns.size();i++){
		if(parent_conns[i]==child){
			isGood=true;
			break;
		}
	}
	for(unsigned i=0;i<child_conns.size();i++){
		if(child_conns[i]==parent){
			isGood=true;
			break;
		}
	}
	if(!isGood) return;

	std::vector<int> new_parent_conns;
	new_parent_conns.push_back(child);
	for(unsigned i=0;i<parent_conns.size();i++){
		if(parent_conns[i]!=child) new_parent_conns.push_back(parent_conns[i]);
	}

	std::vector<int> new_child_conns;
	new_child_conns.push_back(parent);
	for(unsigned i=0;i<child_conns.size();i++){
		if(child_conns[i]!=parent) new_child_conns.push_back(child_conns[i]);
	}
	connectivity[parent] = new_parent_conns;
	connectivity[child] = new_child_conns;
}

void Tree::SetCoords(const std::vector<Cartesian> &SelAtoms_in, int ioff, const std::vector<std::vector<int> > &conn_lists, const std::vector<std::vector<Cartesian> >&ext_carts_in, const std::vector<std::vector<int> > &forced_connections){

  clock_t tv1, tv2;
  clock_t tvs;
  bool print_timing = false;

  tvs = clock();
  tv1 = clock();
  std::vector<std::vector<Cartesian> >ext_carts;

  //First thing to do is apply permutation to swap ioff(start pos) with 0.
  std::vector<Cartesian> SelAtoms = SelAtoms_in;
  connectivity = conn_lists;
  ext_carts = ext_carts_in;

  for(unsigned ii=0;ii<connectivity.size();ii++){
    for(unsigned jj=0;jj<connectivity[ii].size();jj++){
      if(std::find(connectivity[connectivity[ii][jj]].begin(),connectivity[connectivity[ii][jj]].end(),(int)ii)!=connectivity[connectivity[ii][jj]].end()){
      }else if(connectivity[ii][jj]!=(int)ii){
        connectivity[connectivity[ii][jj]].push_back(ii);
      }
    }
  }

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

  ClearCoords();
  extra_bonded_pairs.clear();

  tv2 = clock();

  if(SelAtoms.size()==0)
    return;

  if ( print_timing) printf("Time to do setup: %f\n",(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0) / (double) 1000.0);

  std::vector<std::vector<Cartesian> >::iterator ext_iter = ext_carts.begin();
  tv1 = clock();
  for(unsigned i=0;i<SelAtoms.size();i++){
    coords.push_back(new TreeVertex());
  }

  std::vector<std::vector<int> >::const_iterator forced_iter=forced_connections.begin();
  while(forced_iter!=forced_connections.end()){
	if(forced_iter->size()>1){
		bool isGood=true;
		for(unsigned iforce=0;iforce<forced_iter->size();iforce++){
			std::cout << (*forced_iter)[iforce] << "\n";
			if(((*forced_iter)[iforce]>0)&&(unsigned)((*forced_iter)[iforce])>(SelAtoms.size())) isGood=false;
		}
		if(isGood){
		   // std::cout << "We have a possibly good forced conn\n";
		   std::vector<int>::const_iterator this_iter = forced_iter->begin();
		   int parent = *this_iter;
		   // std::cout << "----" << parent << "----\n";
		   while(this_iter!=forced_iter->end()-1){
		      this_iter++;
		      int child = *this_iter;
		      // std::cout << parent << " " << child << "\n";
		      ForceEarlyConnection(parent,child);
		      parent = *this_iter;
		   }
		}
	}
	forced_iter++;
  }
	
  for(unsigned i=0;i<SelAtoms.size();i++){
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
         RecurseCalculateTree(coordk);
       } else {  //added by Joel Bard
         //but coord may not have a parent so add it as
         //a child to coordk
         if( coord->GetParentID()==-1&&scanned[*k]!=1 ) {
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
	//std::cout << i << " unbonded\n";
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

  c << std::setprecision(6);
  c << std::fixed;

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


  c.unsetf(std::ios::fixed | std::ios::scientific);
  return c;
}

void TreeVertex::GetNonDescendants(std::vector<TreeVertex *> &desc_vertices, const std::vector<TreeVertex*> stop_nodes_in) const {
  TreeVertex *top = (TreeVertex *)(this);
  std::vector<TreeVertex* > stop_nodes = stop_nodes_in;
  //std::cout << "This id: " << GetID() << "\n";
  while(top->GetParent()){
    top=top->GetParent();
    for(unsigned ii=0;ii<stop_nodes.size();ii++){
      if(top->GetID()==stop_nodes[ii]->GetID())
        return;
    }
    //std::cout << "parent id: " << top->GetID() << "\n";
  }
  stop_nodes.push_back((TreeVertex *)(this));
  if(top==(TreeVertex *)(this)) return;
  desc_vertices.push_back(top);
  top->GetDescendants(desc_vertices,stop_nodes);
}

void TreeVertex::GetDescendants(std::vector<TreeVertex *> &desc_vertices, const std::vector<TreeVertex*> stop_nodes) const {
  for(unsigned ii=0;ii<stop_nodes.size();ii++){
    if((TreeVertex *)(this)==stop_nodes[ii])
      return;
  }
  TreeVertex *child;
  for(int i=0;i<GetNumberOfChildren();i++){
    child = GetChild(i);
    bool this_child_ok = true;
    for(unsigned ii=0;ii<stop_nodes.size();ii++){
      if(child==stop_nodes[ii]){
        this_child_ok = false;
      }
    }
    if(this_child_ok){
      desc_vertices.push_back(child);
      child->GetDescendants(desc_vertices,stop_nodes);
    }
  }
}

void TreeVertex::GetNonDescendants(std::vector<TreeVertex *> &desc_vertices) const {
  TreeVertex *top = (TreeVertex *)(this);
  while(top->GetParent()){
    top=top->GetParent();
  }
  top->GetDescendants(desc_vertices,this);
}

void TreeVertex::GetDescendants(std::vector<TreeVertex *> &desc_vertices, const TreeVertex *stop_node) const {
  if((this)==stop_node){
    return;
  }
  TreeVertex *child;
  for(int i=0;i<GetNumberOfChildren();i++){
    child = GetChild(i);
    if(child!=stop_node){
      desc_vertices.push_back(child);
      child->GetDescendants(desc_vertices,stop_node);
    }
  }
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

void Tree::RecurseZMatrix(std::ostream &c,const TreeVertex *vertex, const std::vector<std::string> &labels, const std::string &separater){
  static int depth;
  depth++;
  int id = vertex->GetID();
  if(depth==1){
    c << labels[id] << separater << "X1" << separater << vertex->GetParentDistance() << separater << "X2" << separater << -vertex->GetParentBondAngle()*180/M_PI << "\n";
  }
  if(depth==2){
    TreeVertex *parent = vertex->GetParent();
    int pid = parent->GetID();
    if(0){
      c << labels[id] << separater  << labels[pid] << separater  << vertex->GetParentDistance() << separater << "X1" << separater << vertex->GetParentBondAngle()*180/M_PI << "\n";
    }else{
      c << labels[id] << separater << labels[pid] << separater << vertex->GetParentDistance() << separater << "X1" << separater << vertex->GetParentBondAngle()*180/M_PI << separater << "X2" << separater <<  vertex->GetParentDihedralAngle()*180/M_PI << "\n";
    }
  }
  if(depth==3){
    TreeVertex *parent = vertex->GetParent();
    int pid = parent->GetID();
    TreeVertex *grandparent = parent->GetParent();
    int gpid = grandparent->GetID();
    c << labels[id] << separater << labels[pid] << separater << vertex->GetParentDistance() << separater << labels[gpid] << separater << vertex->GetParentBondAngle()*180/M_PI << separater << "X1" << separater << vertex->GetParentDihedralAngle()*180/M_PI << "\n";
  }
  if(depth>3){
    TreeVertex *parent = vertex->GetParent();
    int pid = parent->GetID();
    TreeVertex *grandparent = parent->GetParent();
    int gpid = grandparent->GetID();
    TreeVertex *greatgrandparent = grandparent->GetParent();
    int ggpid = greatgrandparent->GetID();
    c << labels[id] << separater << labels[pid] << separater << vertex->GetParentDistance() << separater << labels[gpid] << separater << vertex->GetParentBondAngle()*180/M_PI << separater << labels[ggpid] << separater << vertex->GetParentDihedralAngle()*180/M_PI << "\n";
  }
   for(int i=0;i<vertex->GetNumberOfChildren();i++)
    RecurseZMatrix(c,vertex->GetChild(i),labels,separater);
  depth--;
}

void Tree::PrintZMatrix(const std::vector<std::string> &labels, const std::string &separater){
  PrintZMatrix(std::cout,labels,separater);
}

void Tree::PrintZMatrix(std::ostream &c, const std::vector<std::string> &labels, const std::string &separater){
  c << std::setprecision(6);
  c << std::fixed;
  for(int i=0;i<GetNumberOfVertices();i++){
    TreeVertex *vertex = GetCoord(i);
    if(vertex->GetParentID()==-1){
      Cartesian Q2 = GetCoord(0)->Dummy2;
      Cartesian Q1 = GetCoord(0)->Dummy;
      c << "X2\n" << "X1" << separater << "X2" << separater << LineLength(Q2,Q1) << "\n";
      RecurseZMatrix(c,vertex,labels,separater);
    }
  }
  c.unsetf(std::ios::fixed | std::ios::scientific);
}

TreeVertex* Tree::GetCoord(int i, bool permuted) const {
  if(start>0&&permuted)
    return coords[permutation[i]];
  return coords[i];
}

Tree& Tree::operator=(const Tree &t){

  start = t.start;
  permutation = t.permutation;
  scanned = t.scanned;
  connectivity = t.connectivity;
  std::vector<TreeVertex*> tcoords = t.GetCoords();
  ClearCoords();
  for(int i=0;i<t.GetNumberOfVertices();i++){
    TreeVertex *new_coord = new TreeVertex();
    new_coord->SetParentID(tcoords[i]->GetParentID());
    new_coord->SetID(tcoords[i]->GetID());
    new_coord->SetCoord(tcoords[i]->GetCoord());
    new_coord->SetDummy(tcoords[i]->Dummy,tcoords[i]->Dummy2);
    new_coord->SetParentDistance(tcoords[i]->GetParentDistance());
    new_coord->SetParentBondAngle(tcoords[i]->GetParentBondAngle());
    new_coord->SetParentDihedralAngle(tcoords[i]->GetParentDihedralAngle());
    coords.push_back(new_coord);
  }
  for(int i=0;i<GetNumberOfVertices();i++){
    coords[i]->SetParent(coords[coords[i]->GetParentID()]);
  }
  for(int i=0;i<GetNumberOfVertices();i++){
    int nchild = tcoords[i]->GetNumberOfChildren();
    for(int j=0;j<nchild;j++){
      coords[i]->AddChild(coords[tcoords[i]->GetChild(j)->GetID()]);
    }
  }
  for(int i=0;i<GetNumberOfVertices();i++){
    int nchild = tcoords[i]->GetNumberOfExternalChildren();
    for(int j=0;j<nchild;j++){
      coords[i]->AddExternalChild(tcoords[i]->GetExternalChild(j));
    }
  }
  return *this;
}

Tree::Tree(const Tree &t){

  permutation = t.permutation;
  start = t.start;
  scanned = t.scanned;
  connectivity = t.connectivity;
  std::vector<TreeVertex*> tcoords = t.GetCoords();
  ClearCoords();
  for(int i=0;i<t.GetNumberOfVertices();i++){
    TreeVertex *new_coord = new TreeVertex();
    new_coord->SetParentID(tcoords[i]->GetParentID());
    new_coord->SetID(tcoords[i]->GetID());
    new_coord->SetCoord(tcoords[i]->GetCoord());
    new_coord->SetDummy(tcoords[i]->Dummy,tcoords[i]->Dummy2);
    new_coord->SetParentDistance(tcoords[i]->GetParentDistance());
    new_coord->SetParentBondAngle(tcoords[i]->GetParentBondAngle());
    new_coord->SetParentDihedralAngle(tcoords[i]->GetParentDihedralAngle());
    coords.push_back(new_coord);
  }
  for(int i=0;i<GetNumberOfVertices();i++){
    coords[i]->SetParent(coords[coords[i]->GetParentID()]);
  }
  for(int i=0;i<GetNumberOfVertices();i++){
    int nchild = tcoords[i]->GetNumberOfChildren();
    for(int j=0;j<nchild;j++){
      coords[i]->AddChild(coords[tcoords[i]->GetChild(j)->GetID()]);
    }
  }
  for(int i=0;i<GetNumberOfVertices();i++){
    int nchild = tcoords[i]->GetNumberOfExternalChildren();
    for(int j=0;j<nchild;j++){
      coords[i]->AddExternalChild(tcoords[i]->GetExternalChild(j));
    }
  }
    
}

std::vector <TreeVertex*> Tree::GetCoords(bool permuted) const {

  if(start>0&&permuted){
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

  if(start>0&&permuted){
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
  
  if(start>0&&permuted){
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

void Tree::SetDihedralAngle(int baseAtom, int atom, int child, int movingAtom, double TorsionAngle){
   SetDihedralAngle(atom,child,TorsionAngle,false,movingAtom,baseAtom);
}

void Tree::SetDihedralAngle(int atom_in, int child_in, double TorsionAngle, bool permuted, int movingAtom, int baseAtom){
//    std::cout << "SetDihedralAngle 2\n";

//    std::cout << "   atom indices: atom " << atom_in << " child: " << child_in
// 	     << " movingAtom " << movingAtom << " baseAtom " << baseAtom
// 	     << " permuted: " << permuted << " to angle: " << TorsionAngle << std::endl; 

  int atom = atom_in;
  int child = child_in;
  
  if(start>0&&permuted){
    atom = permutation[atom];
    child = permutation[child];
  }

  TreeVertex *E = NULL;
  if(movingAtom>-1)
	  E = coords[movingAtom];

  TreeVertex *O = NULL;
  if(baseAtom>-1)
	  O = coords[baseAtom];

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

  if(C->GetNumberOfChildren()==0||!E)
    return;

  double angle_orig;
  double diff;

  if(E){
	  if(O){
  		angle_orig = DihedralAngle(E->GetCoord(),C->GetCoord(),B->GetCoord(),O->GetCoord());
	  }else{
  		angle_orig = E->GetParentDihedralAngle();
	  }
  }else{
	  if(O){
  		angle_orig = DihedralAngle(C->GetChildren()[0]->GetCoord(),C->GetCoord(),B->GetCoord(),O->GetCoord());
	  }else{
  		angle_orig = C->GetChildren()[0]->GetParentDihedralAngle();
	  }
  }
  diff = TorsionAngle - angle_orig;

  RotateAboutBond(atom_in,child_in,diff,permuted);

  /*
  if(!E){
  	std::vector<TreeVertex*> children = C->GetChildren();
  	std::vector<TreeVertex*>::const_iterator child_iter=children.begin();

  	(*child_iter)->SetParentDihedralAngle(TorsionAngle);
  	child_iter++;

  	while(child_iter!=children.end()){
    		(*child_iter)->SetParentDihedralAngle((*child_iter)->GetParentDihedralAngle()+diff);
    		child_iter++;
  	}
  } else {
  	E->SetParentDihedralAngle(TorsionAngle);
  	std::vector<TreeVertex*> children = C->GetChildren();
  	std::vector<TreeVertex*>::const_iterator child_iter=children.begin();
  	while(child_iter!=children.end()){
    		if(E!=(*child_iter))(*child_iter)->SetParentDihedralAngle((*child_iter)->GetParentDihedralAngle()+diff);
    		child_iter++;
  	}
  }
  */

}

std::vector<std::vector<int> > Tree::FindLongBranches(int req_depth){
  std::vector<std::vector<int> > branches;

  //int max_depth = FindMaxDepth();
  //std::cout << "Max tree depth " << max_depth << "\n";
  //std::cout << "Searching for branches >= " << req_depth << "\n";

  std::vector<int> nodes;

  for(unsigned ii=0;ii<coords.size();ii++){
    int this_depth = coords[ii]->FindDepth();
    if(this_depth>= req_depth) {
      //std::cout << ii << " fits the bill " << this_depth << "\n";
      if(coords[ii]->GetNumberOfChildren()==0) { nodes.push_back(ii); /*std::cout << " branch end\n"; */}
    }
  }
  //std::cout << nodes.size() << " branches\n";
  for(unsigned ii=0;ii<nodes.size();ii++){
    TreeVertex *coord = coords[nodes[ii]];
    branches.push_back(std::vector<int>(0));
    branches.back().push_back(nodes[ii]);
    //std::cout << "--" << nodes[ii] << "\n";
    while((coord=coord->GetParent())){
      branches.back().push_back(coord->GetID());
      //std::cout << coord->GetID() << "\n";
    }
    std::reverse(branches.back().begin(),branches.back().end());
  }

  return branches;
}

void Tree::AddVertex(int pid, double dist, int gpid, double angle, int ggpid, double dihedral, int chiral){
  double rtod = 180.0 / M_PI;
  TreeVertex *new_coord = new TreeVertex();
  new_coord->SetParentID(pid);
  new_coord->SetParent(coords[pid]);
  coords.push_back(new_coord);
  coords[pid]->AddChild(new_coord);
  new_coord->SetParentDistance(dist);
  new_coord->SetID(coords.size()-1);
  if(gpid==-1){
    new_coord->SetParentBondAngle(angle/rtod);
    new_coord->SetParentDihedralAngle(dihedral/rtod);
    new_coord->SetCoord(GetCartFrom3Carts(new_coord->GetParent()->GetCoord(),
       dist,
       new_coord->GetParent()->GetParent()->GetCoord(),
       angle/rtod,
       new_coord->GetParent()->GetParent()->GetParent()->GetCoord(),
       dihedral/rtod));
  } else if(chiral!=0) {
     new_coord->SetCoord(GetCartFrom3Carts(new_coord->GetParent()->GetCoord(),
        dist,
        coords[gpid]->GetCoord(),
        angle/rtod,
        coords[ggpid]->GetCoord(),
        dihedral/rtod,
        chiral));
      new_coord->SetParentBondAngle(Angle(new_coord->GetCoord(),new_coord->GetParent()->GetCoord(),new_coord->GetParent()->Dummy));
      new_coord->SetParentDihedralAngle(DihedralAngle(new_coord->GetCoord(),new_coord->GetParent()->GetCoord(),new_coord->GetParent()->Dummy,new_coord->GetParent()->Dummy2));
  } else {
    // .... ?
  }
  
}

void Tree::ClearCoords(){
  for(unsigned i=0;i<coords.size();i++)
    delete coords[i];
  coords.clear();
}
