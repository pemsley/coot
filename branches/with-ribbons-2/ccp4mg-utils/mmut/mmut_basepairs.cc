/*
     mmut/mmut_basepairs.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

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
#include <stdlib.h>
#include <math.h>
#include "mman_manager.h"
#include "mmut_basepairs.h"
#include "mmut_hbond.h"
#include "plane.h"

CNABasePairs::CNABasePairs(CMMANManager *molHnd, int selHnd, PPCAtom selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector){
  Calculate(molHnd,selHnd,selAtoms,nSelAtoms,atom_colour_vector);
}

void CNABasePairs::Calculate(CMMANManager *molHnd, int selHnd, PPCAtom selAtoms, int nSelAtoms, const AtomColourVector &atom_colour_vector){
  
  int C5sel = molHnd->NewSelection();
  
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*","C","*",SKEY_NEW);
  molHnd->Select(C5sel,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5'","C","*",SKEY_OR);

  CHBond hbonds(molHnd, selHnd);
  hbonds.Calculate();

  std::vector<std::pair<PCResidue,Plane> > plane_res_pairs;
  std::vector<double*> plane_res_colours;

  //std::cout << "In base pairs calc " << nSelAtoms << " atoms\n"; std::cout.flush();
  //std::cout << "In base pairs calc " << hbonds.hbonds.GetNofConnections() << " hbonds\n"; std::cout.flush();
  
  char ID[50];
  char ID1[50];
  char ID2[50];
  char AtomID1[50];
  char AtomID2[50];
  std::map<std::string,std::map<std::string,int> > alreadyDone;
  for(int i=0;i<nSelAtoms;i++){
    if(selAtoms[i]->isInSelection(C5sel)){
      PCResidue res = selAtoms[i]->GetResidue();
      if(res){
        int restype = molHnd->GetRestypeCode(res);
        if(restype==RESTYPE_NUCL||restype==RESTYPE_DNA||restype==RESTYPE_RNA){
          //std::cout << "Considering atom " << selAtoms[i]->name << " in base " << selAtoms[i]->GetChainID() << "/" << selAtoms[i]->GetResidueNo() << "\n";
	  for(int j=0;j<hbonds.hbonds.GetNofConnections();j++){
	     const char* thisID = res->GetResidueID(ID);
             PCAtom atom1 = hbonds.hbonds.pAtom1[j];
             PCAtom atom2 = hbonds.hbonds.pAtom2[j];
             PCResidue res1 = atom1->GetResidue();
             PCResidue res2 = atom2->GetResidue();
	     const char* id1 = res1->GetResidueID(ID1);
	     const char* id2 = res2->GetResidueID(ID2);
	     //const char* atomid1 = atom1->GetAtomID(AtomID1);
	     //const char* atomid2 = atom2->GetAtomID(AtomID2);
             int restype1 = molHnd->GetRestypeCode(res1);
             int restype2 = molHnd->GetRestypeCode(res2);
	     if((strncmp(thisID,id2,50)==0&&strncmp(id1,id2,50)!=0&&abs(res1->GetSeqNum()-res2->GetSeqNum())!=1&&(restype1==RESTYPE_NUCL||restype1==RESTYPE_DNA||restype1==RESTYPE_RNA))||(strncmp(thisID,id1,50)==0&&strncmp(id1,id2,50)!=0&&abs(res1->GetSeqNum()-res2->GetSeqNum())!=1&&(restype2==RESTYPE_NUCL||restype2==RESTYPE_DNA||restype2==RESTYPE_RNA))){
		     if(
			   ((
			   std::string(res2->name)==std::string("C")||
			   std::string(res2->name)==std::string("DC")||
			   std::string(res2->name)==std::string("Cd")||
			   std::string(res2->name)==std::string("CYT")
			   )&&
			   (
			   std::string(res1->name)==std::string("T")||
			   std::string(res1->name)==std::string("DT")||
			   std::string(res1->name)==std::string("Td")||
			   std::string(res1->name)==std::string("THY")
			   ))||
			   ((
			   std::string(res2->name)==std::string("T")||
			   std::string(res2->name)==std::string("DT")||
			   std::string(res2->name)==std::string("Td")||
			   std::string(res2->name)==std::string("THY")
			   )&&
			   (
			   std::string(res1->name)==std::string("C")||
			   std::string(res1->name)==std::string("DC")||
			   std::string(res1->name)==std::string("Cd")||
			   std::string(res1->name)==std::string("CYT")
			   ))||
			   ((
			   std::string(res2->name)==std::string("A")||
			   std::string(res2->name)==std::string("DA")||
			   std::string(res2->name)==std::string("Ad")||
			   std::string(res2->name)==std::string("ADE")
			   )&&
			   (
			   std::string(res1->name)==std::string("G")||
			   std::string(res1->name)==std::string("DG")||
			   std::string(res1->name)==std::string("Gd")||
			   std::string(res1->name)==std::string("GUA")
			   ))||
			   ((
			   std::string(res2->name)==std::string("G")||
			   std::string(res2->name)==std::string("DG")||
			   std::string(res2->name)==std::string("Gd")||
			   std::string(res2->name)==std::string("GUA")
			   )&&
			   (
			   std::string(res1->name)==std::string("A")||
			   std::string(res1->name)==std::string("DA")||
			   std::string(res1->name)==std::string("Ad")||
			   std::string(res1->name)==std::string("ADE")
			   ))||
			   (
			    molHnd->BondLength(atom1,atom2)>3.3
			   )||
			   (
			    std::string(atom1->name).find('\'')!=std::string::npos
			   )||
			   (
			    std::string(atom1->name).find('*')!=std::string::npos
			   )||
			   (
			    std::string(atom2->name).find('\'')!=std::string::npos
			   )||
			   (
			    std::string(atom2->name).find('*')!=std::string::npos
			   )){
                              continue;
                 }
                 /*
                 if(alreadyDone[atomid1].find(std::string(atomid2))==alreadyDone[atomid1].end()){
                    alreadyDone[atomid1][atomid2] = 1;
                 } else {
                    std::cout << "Rejecting duplicate " << atomid1 << " " << atomid2 << ", with length: " << molHnd->BondLength(atom1,atom2) << "\n"; std::cout.flush();
                    continue;
                 }
                 */
		 const double *col = atom_colour_vector.GetRGB(i); 
		 colours.push_back(std::pair<const double*,const double*>(col,col));
                 //std::cout << "Accepting base pair " << atomid1 << " " << atomid2 << ", with length: " << molHnd->BondLength(atom1,atom2) << "\n"; std::cout.flush();
	         if(strncmp(thisID,id2,50)==0&&strncmp(id1,id2,50)!=0&&abs(res1->GetSeqNum()-res2->GetSeqNum())!=1&&(restype1==RESTYPE_NUCL||restype1==RESTYPE_DNA||restype1==RESTYPE_RNA)){
		     base_pairs.push_back(std::pair<PCResidue,PCResidue>(res,res1));
	         }
	         if(strncmp(thisID,id1,50)==0&&strncmp(id1,id2,50)!=0&&abs(res1->GetSeqNum()-res2->GetSeqNum())!=1&&(restype2==RESTYPE_NUCL||restype2==RESTYPE_DNA||restype2==RESTYPE_RNA)){
		     base_pairs.push_back(std::pair<PCResidue,PCResidue>(res,res2));
	         }
	     }
	  }
	  //std::cout.flush();
        }
      }
    }
  }

  //std::cout << base_pairs.size() << " " << colours.size() << "\n"; std::cout.flush();
  //  Now let's deal with the mess
  std::vector<std::map<int,int> > GCmap(nSelAtoms);
  for(unsigned i=0;i<base_pairs.size();i++){
     PCResidue res1 = base_pairs[i].first;
     PCResidue res2 = base_pairs[i].second;
     int seq1 = res1->GetSeqNum();
     int seq2 = res2->GetSeqNum();
     if(GCmap[seq2].find(seq1)==GCmap[seq2].end()){
       GCmap[seq2][seq1]=0;
     }
     GCmap[seq2][seq1]++;
     if(GCmap[seq1].find(seq2)==GCmap[seq1].end()){
       GCmap[seq1][seq2]=0;
     }
     GCmap[seq1][seq2]++;
  }
  
  std::map<int,int> popOff;
  for(unsigned i=0;i<GCmap.size();i++){
    if(GCmap[i].size()!=0) {
      //std::cout << i << "\n";
      std::map<int,int>::iterator map_iter = GCmap[i].begin();
      int max_conn = -1;
      int max_key = -1;
      while(map_iter!=GCmap[i].end()){
        //std::cout << "   " << map_iter->first << " " << map_iter->second << "\n";
        if(map_iter->second>max_conn){
           max_key = map_iter->first;
           max_conn = map_iter->second;
        }
        ++map_iter;
      }
      if(max_key>-1){
        map_iter = GCmap[i].begin();
        while(map_iter!=GCmap[i].end()){
          if(max_key!=map_iter->first){
             //std::cout << "Pop off " << i << " " << map_iter->first << "\n";
             popOff[i]=map_iter->first;
          }
          ++map_iter;
        }
      }
    }
  }

  std::vector<std::pair<PCResidue,PCResidue> > base_pairs_new;
  std::vector<std::pair<const double*,const double*> > colours_new;

  for(unsigned i=0;i<base_pairs.size();i++){
     PCResidue res1 = base_pairs[i].first;
     PCResidue res2 = base_pairs[i].second;
     int seq1 = res1->GetSeqNum();
     int seq2 = res2->GetSeqNum();
     if(popOff.find(seq1)!=popOff.end()&&popOff[seq1]==seq2&&popOff.find(seq2)!=popOff.end()&&popOff[seq2]==seq1){
        //std::cout << "Really pop off " << seq1 << " " << seq2 << "\n";
        continue;
     } else {
        base_pairs_new.push_back(base_pairs[i]);
        colours_new.push_back(colours[i]);
     }
  }
  base_pairs=base_pairs_new;
  colours=colours_new;

}

PCResidue CNABasePairs::GetPairedResidue(const PCResidue res_in) const {
  return res_in;
}

int CNABasePairs::GetPairedResidueIndex(const int res_in) const {
  return res_in;
}
