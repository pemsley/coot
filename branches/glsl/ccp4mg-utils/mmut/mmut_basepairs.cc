/*
     mmut/mmut_basepairs.cc: CCP4MG Molecular Graphics Program
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

#include <string.h>
#include <vector>
#include <utility>
#include <math.h>
#include <mman_manager.h>
#include "mmut_basepairs.h"
#include "plane.h"

CNABasePairs::CNABasePairs(CMMANManager *molHnd, int selHnd, PPCAtom selAtoms, int nSelAtoms, AtomColourVector *atom_colour_vector){
  Calculate(molHnd,selHnd,selAtoms,nSelAtoms,atom_colour_vector);
}

void CNABasePairs::Calculate(CMMANManager *molHnd, int selHnd, PPCAtom selAtoms, int nSelAtoms, AtomColourVector *atom_colour_vector){
  
  int C5sel = molHnd->NewSelection();
  
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*","C","*",SKEY_NEW);
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5'","C","*",SKEY_OR);

  std::vector<std::pair<PCResidue,Plane> > plane_res_pairs;
  std::vector<double*> plane_res_colours;
  
  for(int i=0;i<nSelAtoms;i++){
    if(selAtoms[i]->isInSelection(C5sel)){
      PCResidue res = selAtoms[i]->GetResidue();
      if(res){
        int restype = molHnd->GetRestypeCode(res);
        if(restype==RESTYPE_NUCL||restype==RESTYPE_DNA||restype==RESTYPE_RNA){
          //std::cout << "Considering atom " << selAtoms[i]->name << " in base " << selAtoms[i]->GetChainID() << "/" << selAtoms[i]->GetResidueNo() << "\n";
          PCAtom c4 = res->GetAtom("C4");
          PCAtom c5 = res->GetAtom("C5");
          PCAtom c6 = res->GetAtom("C6");
          if(c4&&c5&&c6){
            Cartesian cart_c4(c4->x,c4->y,c4->z);
            Cartesian cart_c5(c5->x,c5->y,c5->z);
            Cartesian cart_c6(c6->x,c6->y,c6->z);
            Plane plane(cart_c4,cart_c5,cart_c6);
            // We give all planes a +ve D. This means that planes
            // with similar values of D and large positive
            // overlaps of plane normals are paired. Those with
            // large -ve overlaps are not paired.
            if(plane.get_D()<0.0){
              plane.set_A(-plane.get_A());
              plane.set_B(-plane.get_B());
              plane.set_C(-plane.get_C());
              plane.set_D(-plane.get_D());
            }
            //std::cout << "Gives Plane " << plane << "\n";
            plane_res_pairs.push_back(std::pair<PCResidue,Plane>(res,plane));
            if(atom_colour_vector){
              double *col_tmp = atom_colour_vector->GetRGB(i);
              plane_res_colours.push_back(col_tmp);
            }else{
              double *dum;
              plane_res_colours.push_back(dum);
            }
          }
        }
      }
    }
  }
  std::vector<int> paired(plane_res_pairs.size());
  for(unsigned ii=0;ii<plane_res_pairs.size();ii++){
    if(paired[ii]==0){
    //std::cout << ii << " " << plane_res_pairs[ii].second << "\n";
    double min_dot = 1.0;
    bool have_pair = false;
    //for(unsigned jj=0;jj<plane_res_pairs.size()&&jj!=ii;jj++){
    int final_jj=-1;
    for(unsigned jj=0;jj<plane_res_pairs.size();jj++){
	    if(jj!=ii){
      Plane p1 = plane_res_pairs[ii].second;
      Plane p2 = plane_res_pairs[jj].second;
      Cartesian n1 = p1.get_normal();
      Cartesian n2 = p2.get_normal();
      n1.normalize();
      n2.normalize();
      double overlap = Cartesian::DotProduct(n1,n2);
      if(overlap>0.7){
        PCResidue res1 = plane_res_pairs[ii].first;
        PCResidue res2 = plane_res_pairs[jj].first;
        PCAtom c41 = res1->GetAtom("C4");
        PCAtom c42 = res2->GetAtom("C4");
        Cartesian cart1(c41->x,c41->y,c41->z);
        Cartesian cart2(c42->x,c42->y,c42->z);
	Cartesian c1c2 = cart2 - cart1;
	c1c2.normalize();
	double dot = fabs(Cartesian::DotProduct(c1c2,(n1+n2)*.5));
	//std::cout << c41->GetChainID() << "/" << c41->GetResidueNo() << ", " << c42->GetChainID() << "/" << c42->GetResidueNo()  << ": " << overlap << " (" << dot << ")\n";
        if(dot<min_dot){
	  Cartesian avg_n = (n1+n2)*0.5;
	  Cartesian p2p1 = cart1-cart2;
	  p2p1.normalize();
	  //std::cout << Cartesian::DotProduct(p2p1,avg_n) << "\n";
          if((cart1-cart2).length()<7.9){
		  //std::cout << "short enough\n";
            if(have_pair){ 
              base_pairs.pop_back();
              colours.pop_back();
            }
            double *col1 = plane_res_colours[ii];
            double *col2 = plane_res_colours[jj];
	    //std::cout << res1->GetResidueNo() << ", "  << res2->GetResidueNo() << "\n";
            base_pairs.push_back(std::pair<PCResidue,PCResidue>(res1,res2));
            colours.push_back(std::pair<double*,double*>(col1,col2));
            have_pair = true;
            min_dot = dot;
	    final_jj = jj;
          }
	  //else std::cout << "not short enough" << (cart1-cart2).length() << "\n";
        }
      }
    }
    }
    if(final_jj>=0)
      paired[final_jj] = 1;
    }
  }
}

PCResidue CNABasePairs::GetPairedResidue(const PCResidue res_in) const {
  return res_in;
}

int CNABasePairs::GetPairedResidueIndex(const int res_in) const {
  return res_in;
}
