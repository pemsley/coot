/*
     mmut/mmut_morph.cc: CCP4MG Molecular Graphics Program
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
#include <vector>
#include <iostream>
#include <cartesian.h>
#include <mgtree.h>
#include "mmut_morph.h"

#include <math.h>
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#define PIBY2 (M_PI * 2)
#endif

std::vector<std::vector<Cartesian> > CMMUTMorph::MorphCartesians(const std::vector<Cartesian> &a, const std::vector<Cartesian> &b, int nsteps){
  std::vector<std::vector<Cartesian> > morphs;
  return morphs;
}

std::vector<Tree> CMMUTMorph::MorphTrees(const Tree &a, const Tree &b, int nsteps){
  std::vector<Tree> morphs;
  std::vector <TreeVertex*> coords_a = a.GetCoords();
  std::vector <TreeVertex*> coords_b = b.GetCoords();

  unsigned len;
  if(coords_a.size()<coords_b.size()) len = coords_a.size(); else len = coords_b.size();

  std::vector<double> diff_lengths;
  std::vector<double> diff_angles;
  std::vector<double> diff_torsions;
  std::vector<Cartesian> diff_posns;
  std::vector<Cartesian> diff_dummys;
  std::vector<Cartesian> diff_dummy2s;

  for(unsigned i=0;i<len;i++){
    if(i<3){
        diff_posns.push_back((coords_b[i]->GetCoord()-coords_a[i]->GetCoord())/nsteps);
      if(i<1){
        diff_dummys.push_back((coords_b[i]->Dummy-coords_a[i]->Dummy)/nsteps);
        diff_dummy2s.push_back((coords_b[i]->Dummy2-coords_a[i]->Dummy2)/nsteps);
      }
    }
    diff_lengths.push_back((coords_b[i]->GetParentDistance()-coords_a[i]->GetParentDistance())/nsteps);
    diff_angles.push_back((coords_b[i]->GetParentBondAngle()-coords_a[i]->GetParentBondAngle())/nsteps);
    if(fabs(coords_b[i]->GetParentDihedralAngle()-coords_a[i]->GetParentDihedralAngle())>M_PI){
      if(coords_b[i]->GetParentDihedralAngle()<0.0){
        diff_torsions.push_back((2*M_PI+coords_b[i]->GetParentDihedralAngle()-coords_a[i]->GetParentDihedralAngle())/nsteps);
      } else {
        diff_torsions.push_back((-2*M_PI+coords_b[i]->GetParentDihedralAngle()-coords_a[i]->GetParentDihedralAngle())/nsteps);
      }
    } else {
      diff_torsions.push_back((coords_b[i]->GetParentDihedralAngle()-coords_a[i]->GetParentDihedralAngle())/nsteps);
    }
  }
  
  Tree c=a;
  morphs.reserve(nsteps+1);
  for(int j=0;j<=nsteps;j++){
    for(unsigned i=0;i<len;i++){
      TreeVertex* coord = c.GetCoord(i);
      coord->SetParentDistance(coords_a[i]->GetParentDistance()+j*diff_lengths[i]);
      coord->SetParentBondAngle(coords_a[i]->GetParentBondAngle()+j*diff_angles[i]);
      coord->SetParentDihedralAngle(coords_a[i]->GetParentDihedralAngle()+j*diff_torsions[i]);
      if(i<3){
        coord->SetCoord(coords_a[i]->GetCoord()+j*diff_posns[i]);
        if(i<1){
          Cartesian d1 = coords_a[i]->Dummy+j*diff_dummys[i];
          Cartesian d2 = coords_a[i]->Dummy2+j*diff_dummy2s[i];
          coord->SetDummy(d1,d2);
        }
      }
    }
    morphs.push_back(c);
  }

  return morphs;
}
