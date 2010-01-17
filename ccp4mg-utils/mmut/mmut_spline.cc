/*
     mmut/mmut_spline.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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

#include <math.h>
#include "mgtree.h"
#include "cartesian.h"
#include <mman_manager.h>
#include <catmull.h>
#include "splineinfo.h"
#include "atom_util.h"
#include "mmut_basepairs.h"
#include <list>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#ifndef HALF_PI
#define HALF_PI (M_PI * 0.5)
#endif

Cartesian AtToCart(PCAtom at){
  return Cartesian(at->x,at->y,at->z);
}

bool IsInBetaIndex(const int &i, const std::vector<std::vector<int> > &indices){
  return true;
}

void ReplaceSmoothly(std::vector<Cartesian> &old_vec, const std::vector<Cartesian> &new_vec, const int start, const int end){

  int ii;
  if(end-start>20){
    for(ii=start;ii<start+11&&ii<end&&ii<(int)old_vec.size();ii++){
      double frac = sin(double (ii-start)/10.0 * HALF_PI);
      old_vec[ii] = frac * new_vec[ii-start] + (1.0-frac) * old_vec[ii];
    }
    for(;ii<end-11&&ii<(int)old_vec.size();ii++){
      old_vec[ii] = new_vec[ii-start];
    }
    for(;ii<end&&ii<(int)old_vec.size();ii++){
      double frac = sin(double (end-ii-1)/10.0 * HALF_PI);
      old_vec[ii] = frac * new_vec[ii-start] + (1.0-frac) * old_vec[ii];
    }
  }

}

void Replace(std::vector<Cartesian> &old_vec, const std::vector<Cartesian> &new_vec, const int start, const int end){
   old_vec.erase(old_vec.begin()+start,old_vec.begin()+end);
   old_vec.insert(old_vec.begin()+start,new_vec.begin(),new_vec.begin()+end-start);
}

std::vector<Cartesian> GetBasePairEnds(PCResidue res1, PCResidue res2);

int GetCAFromSelection(CMMANManager *molH, int atom_selHnd_in){

  //std::cout << "GetSplineInfo spline_accu " << spline_accu << "\n";
  //std::cout << "GetSplineInfo udd_CA " << udd_CA  << "\n";
  //std::cout << "GetSplineInfo udd_chain " << udd_chain  << "\n";


  int CAselHnd;

  int atom_selHnd = molH->NewSelection();
  molH->Select(atom_selHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_NEW);
  molH->Select(atom_selHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);

  // Find all CA - to use as quick check if atom is in this set

  int CALocAselHnd = molH->NewSelection();
  molH->Select(CALocAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","A",SKEY_NEW);
  molH->Select(CALocAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  int CALocBselHnd = molH->NewSelection();
  molH->Select(CALocBselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","B",SKEY_NEW);
  molH->Select(CALocBselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  int CALocCselHnd = molH->NewSelection();
  molH->Select(CALocCselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","C",SKEY_NEW);
  molH->Select(CALocCselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  
  PPCAtom CALocAAtomTable;
  int nCALocAAtoms;
  PPCAtom CALocBAtomTable;
  int nCALocBAtoms;
  PPCAtom CALocCAtomTable;
  int nCALocCAtoms;

  molH->GetSelIndex ( CALocAselHnd, CALocAAtomTable, nCALocAAtoms );
  molH->GetSelIndex ( CALocBselHnd, CALocBAtomTable, nCALocBAtoms );
  molH->GetSelIndex ( CALocCselHnd, CALocCAtomTable, nCALocCAtoms );

  CAselHnd = molH->NewSelection();
  molH->Select(CAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","*",SKEY_NEW);
  molH->Select(CAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);

  for(int ia=0;ia<nCALocAAtoms;ia++){
	  for(int ib=0;ib<nCALocBAtoms;ib++){
		  if(CALocAAtomTable[ia]->isInSelection(CAselHnd)&&CALocBAtomTable[ib]->isInSelection(CAselHnd)&&CALocAAtomTable[ia]->residue==CALocBAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(CAselHnd,CALocBAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }
  for(int ia=0;ia<nCALocAAtoms;ia++){
	  for(int ib=0;ib<nCALocCAtoms;ib++){
		  if(CALocAAtomTable[ia]->isInSelection(CAselHnd)&&CALocCAtomTable[ib]->isInSelection(CAselHnd)&&CALocAAtomTable[ia]->residue==CALocCAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(CAselHnd,CALocCAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }
  for(int ia=0;ia<nCALocCAtoms;ia++){
	  for(int ib=0;ib<nCALocBAtoms;ib++){
		  if(CALocCAtomTable[ia]->isInSelection(CAselHnd)&&CALocBAtomTable[ib]->isInSelection(CAselHnd)&&CALocCAtomTable[ia]->residue==CALocBAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(CAselHnd,CALocBAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }

  molH->ExcludeOverlappedAtoms(CAselHnd,0.8); // Final catchall

  int NAselHnd; // nucleic selection

  int NALocAselHnd = molH->NewSelection();
  //molH->SelectProperty(NALocAselHnd,SELPROP_Nucleotide,STYPE_RESIDUE,SKEY_NEW);
  //molH->Select(NALocAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","A",SKEY_AND);
  molH->Select(NALocAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*'","C","*",SKEY_NEW);
  molH->Select(NALocAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5\'","C","*",SKEY_OR);
  molH->Select(NALocAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3*'","C","*",SKEY_OR);
  molH->Select(NALocAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3\'","C","*",SKEY_OR);
  molH->Select(NALocAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  int NALocBselHnd = molH->NewSelection();
  //molH->SelectProperty(NALocBselHnd,SELPROP_Nucleotide,STYPE_RESIDUE,SKEY_NEW);
  //molH->Select(NALocBselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","B",SKEY_AND);
  molH->Select(NALocBselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*'","C","*",SKEY_NEW);
  molH->Select(NALocBselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5\'","C","*",SKEY_OR);
  molH->Select(NALocBselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3*'","C","*",SKEY_OR);
  molH->Select(NALocBselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3\'","C","*",SKEY_OR);
  molH->Select(NALocBselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  int NALocCselHnd = molH->NewSelection();
  //molH->SelectProperty(NALocCselHnd,SELPROP_Nucleotide,STYPE_RESIDUE,SKEY_NEW);
  //molH->Select(NALocCselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","C",SKEY_AND);
  molH->Select(NALocCselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*'","C","*",SKEY_NEW);
  molH->Select(NALocCselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5\'","C","*",SKEY_OR);
  molH->Select(NALocCselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3*'","C","*",SKEY_OR);
  molH->Select(NALocCselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3\'","C","*",SKEY_OR);
  molH->Select(NALocCselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);

  PPCAtom NALocAAtomTable;
  int nNALocAAtoms;
  PPCAtom NALocBAtomTable;
  int nNALocBAtoms;
  PPCAtom NALocCAtomTable;
  int nNALocCAtoms;

  molH->GetSelIndex ( NALocAselHnd, NALocAAtomTable, nNALocAAtoms );
  molH->GetSelIndex ( NALocBselHnd, NALocBAtomTable, nNALocBAtoms );
  molH->GetSelIndex ( NALocCselHnd, NALocCAtomTable, nNALocCAtoms );

  NAselHnd = molH->NewSelection();
  //molH->SelectProperty(NAselHnd,SELPROP_Nucleotide,STYPE_RESIDUE,SKEY_NEW);
  //molH->Select(NAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_AND);
  molH->Select(NAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5*'","C","*",SKEY_NEW);
  molH->Select(NAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C5\'","C","*",SKEY_OR);
  molH->Select(NAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3*'","C","*",SKEY_OR);
  molH->Select(NAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","C3\'","C","*",SKEY_OR);
  molH->Select(NAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);

  for(int ia=0;ia<nNALocAAtoms;ia++){
	  for(int ib=0;ib<nNALocBAtoms;ib++){
		  if(NALocAAtomTable[ia]->isInSelection(NAselHnd)&&NALocBAtomTable[ib]->isInSelection(NAselHnd)&&NALocAAtomTable[ia]->residue==NALocBAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(NAselHnd,NALocBAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }
  for(int ia=0;ia<nNALocAAtoms;ia++){
	  for(int ib=0;ib<nNALocCAtoms;ib++){
		  if(NALocAAtomTable[ia]->isInSelection(NAselHnd)&&NALocCAtomTable[ib]->isInSelection(NAselHnd)&&NALocAAtomTable[ia]->residue==NALocCAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(NAselHnd,NALocCAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }
  for(int ia=0;ia<nNALocCAtoms;ia++){
	  for(int ib=0;ib<nNALocBAtoms;ib++){
		  if(NALocCAtomTable[ia]->isInSelection(NAselHnd)&&NALocBAtomTable[ib]->isInSelection(NAselHnd)&&NALocCAtomTable[ia]->residue==NALocBAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(NAselHnd,NALocBAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }

  /*
  PPCAtom NAAtomTable; 
  int nNAAtoms;
  molH->GetSelIndex ( NAselHnd, NAAtomTable, nNAAtoms );
  std::cout << nNAAtoms << "\n"; std::cout.flush();
  */

  //molH->MakeSelIndex(NAselHnd);
  molH->Select(CAselHnd,STYPE_ATOM,NAselHnd,SKEY_OR);
  molH->MakeSelIndex(CAselHnd);
  
  molH->DeleteSelection(CALocAselHnd);  
  molH->DeleteSelection(CALocBselHnd);  
  molH->DeleteSelection(CALocCselHnd);  
  molH->DeleteSelection(NALocAselHnd);  
  molH->DeleteSelection(NALocBselHnd);  
  molH->DeleteSelection(NALocCselHnd);
  
  return CAselHnd;
}

SplineInfo GetSplineInfo (CMMANManager *molH, int atom_selHnd_in ,AtomColourVector *atom_colour_vector,int spline_accu, int udd_chain, int udd_CA, int flatten_beta_sheet, int flatten_loop,int smooth_helix){

  //std::cout << "GetSplineInfo spline_accu " << spline_accu << "\n";
  //std::cout << "GetSplineInfo udd_CA " << udd_CA  << "\n";
  //std::cout << "GetSplineInfo udd_chain " << udd_chain  << "\n";

  int CAselHnd;
  PPCAtom atomTable;
  PCResidue pRes;
  int nAtoms;
  PCAtom pCAprev,pCA;
  float dist2;
  int newchain;

  molH->GetSelIndex ( atom_selHnd_in, atomTable, nAtoms );

  SplineInfo splineinfo;
  std::vector<std::vector<PCAtom> > cavertices;
  std::vector<std::vector<PCAtom> > nac5vertices;

  Cartesian pos;

  double colour_array[4] = { 255.0,255.0,255.0,255.0 };
  colour_array[0] = 0.0;
  colour_array[1] = 0.0;
  colour_array[2] = 0.0;
  colour_array[3] = 255.0;


  int atom_selHnd = molH->NewSelection();
  if(nAtoms>0){
  molH->Select(atom_selHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_NEW);
  molH->Select(atom_selHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  }

  molH->GetSelIndex ( atom_selHnd, atomTable, nAtoms );

  // Find all CA - to use as quick check if atom is in this set

  int CALocAselHnd = molH->NewSelection();
  molH->Select(CALocAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","A",SKEY_NEW);
  molH->Select(CALocAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  //molH->ExcludeOverlappedAtoms(CALocAselHnd,0.8);
  int CALocBselHnd = molH->NewSelection();
  molH->Select(CALocBselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","B",SKEY_NEW);
  molH->Select(CALocBselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  //molH->ExcludeOverlappedAtoms(CALocBselHnd,0.8);
  int CALocCselHnd = molH->NewSelection();
  molH->Select(CALocCselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","C",SKEY_NEW);
  molH->Select(CALocCselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);
  //molH->ExcludeOverlappedAtoms(CALocCselHnd,0.8);
  
  PPCAtom CALocAAtomTable;
  int nCALocAAtoms;
  PPCAtom CALocBAtomTable;
  int nCALocBAtoms;
  PPCAtom CALocCAtomTable;
  int nCALocCAtoms;

  molH->GetSelIndex ( CALocAselHnd, CALocAAtomTable, nCALocAAtoms );
  molH->GetSelIndex ( CALocBselHnd, CALocBAtomTable, nCALocBAtoms );
  molH->GetSelIndex ( CALocCselHnd, CALocCAtomTable, nCALocCAtoms );

  CAselHnd = molH->NewSelection();
  molH->Select(CAselHnd,STYPE_ATOM,0,"*",ANY_RES,"*",ANY_RES,"*","*","CA","C","*",SKEY_NEW);
  molH->Select(CAselHnd,STYPE_ATOM,atom_selHnd_in,SKEY_AND);

  for(int ia=0;ia<nCALocAAtoms;ia++){
	  for(int ib=0;ib<nCALocBAtoms;ib++){
		  if(CALocAAtomTable[ia]->isInSelection(CAselHnd)&&CALocBAtomTable[ib]->isInSelection(CAselHnd)&&CALocAAtomTable[ia]->residue==CALocBAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(CAselHnd,CALocBAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }
  for(int ia=0;ia<nCALocAAtoms;ia++){
	  for(int ib=0;ib<nCALocCAtoms;ib++){
		  if(CALocAAtomTable[ia]->isInSelection(CAselHnd)&&CALocCAtomTable[ib]->isInSelection(CAselHnd)&&CALocAAtomTable[ia]->residue==CALocCAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(CAselHnd,CALocCAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }
  for(int ia=0;ia<nCALocCAtoms;ia++){
	  for(int ib=0;ib<nCALocBAtoms;ib++){
		  if(CALocCAtomTable[ia]->isInSelection(CAselHnd)&&CALocBAtomTable[ib]->isInSelection(CAselHnd)&&CALocCAtomTable[ia]->residue==CALocBAtomTable[ib]->residue){
			  // Exclude one of them!
			  molH->SelectAtom(CAselHnd,CALocBAtomTable[ib],SKEY_XOR,False);
		  }
	  }
  }
  molH->MakeSelIndex(CAselHnd);


  molH->ExcludeOverlappedAtoms(CAselHnd,0.8); // Final catchall

  //if (atom_colour_vector)
  //  int rv = atom_colour_vector->SetupResidueColourVector( molH,atom_selHnd );

  molH->GetSelIndex ( atom_selHnd, atomTable, nAtoms );
   
  // Start the vector for the first 'chain'
  splineinfo.colours.push_back(std::vector<Cartesian>(0));
  cavertices.push_back(std::vector<PCAtom>(0));
  nac5vertices.push_back(std::vector<PCAtom>(0));
  std::vector<int> prev_atoms;
  std::vector<int> next_atoms;

  // Loop over all residues creating vector of CA atom pointers
  // and the colour & secstr_indices vector

  newchain = 1;
  int nPept = -1;
  int nChain = 0;
  int nCA = -1;
  for(int j=0;j<nAtoms;j++){
    // Is it an amino acid with CA and O atoms
    //if  (atomTable[j]->isInSelection(CAselHnd) && atomTable[j]->residue->isAminoacid()) {
      if  (atomTable[j]->isInSelection(CAselHnd) && molH->isAminoacid(atomTable[j]->residue) ) {
        pCA = atomTable[j];
        pRes = pCA->residue;
        nPept++;
        nCA++;
        // If the next residue is in a different chain then start a
        // new chain vector
        if(cavertices.back().size()>0){
          pCAprev = cavertices.back().back();
          dist2 = (pCA->x - pCAprev->x)*(pCA->x - pCAprev->x) +
                  (pCA->y - pCAprev->y)*(pCA->y - pCAprev->y)+
	          (pCA->z - pCAprev->z)*(pCA->z - pCAprev->z);
          if(dist2 > 20.0 || 
                pCA->GetChainID()!= pCAprev->GetChainID()){
	    cavertices.push_back(std::vector<PCAtom>(0));
	    splineinfo.colours.push_back(std::vector<Cartesian>(0));
            newchain = 1;
            nChain++;
            nCA = 0;
          }
        }
        if (udd_CA>0) {
          //std::cout << "PutUDData " << nChain << " " << nCA << std::endl;
          pCA->PutUDData(udd_chain,nChain);
          pCA->PutUDData(udd_CA,nCA);
        }
       
        // At start of chain or when there is a change of SSE pushback
        // the residue index and the SSE
        //  NBNB Should treat bulge same as sheet
        // On first pass through (newchain=1) this should not
        // test the secstr_indices, right?
        if ( atom_colour_vector ){
	  if(j==0||newchain){
            double * atcol = atom_colour_vector->GetRGB(j);
            splineinfo.colours.back().push_back(Cartesian(atcol));
	    delete [] atcol;
	  } else {
            if((splineinfo.secstr_indices.back().back().back()==(int)SSE_Strand
              &&((int)pRes->SSE!=(int)SSE_Strand)&&(int)pRes->SSE!=(int)SSE_Bulge)
              ||(splineinfo.secstr_indices.back().back().back()==(int)SSE_Helix
              &&((int)pRes->SSE!=(int)SSE_Helix))
               ){
	      if((int)pRes->SSE==(int)SSE_Helix||(int)pRes->SSE==(int)SSE_Strand||(int)pRes->SSE==(int)SSE_Bulge){
		// This is fine except for case of colour by secondary structure.
                double * atcol1 = atom_colour_vector->GetRGB(j-1);
                splineinfo.colours.back().push_back(Cartesian(atcol1));
	        //delete [] atcol;
	        delete [] atcol1;
	      } else {
                double * atcol = atom_colour_vector->GetRGB(j);
                splineinfo.colours.back().pop_back();
                splineinfo.colours.back().push_back(Cartesian(atcol));
                splineinfo.colours.back().push_back(Cartesian(atcol));
	        delete [] atcol;
	      }
	    } else {
              double * atcol = atom_colour_vector->GetRGB(j);
              splineinfo.colours.back().push_back(Cartesian(atcol));
	      delete [] atcol;
	    }
          }
	}
        if(newchain) {
          splineinfo.secstr_indices.push_back(std::vector<std::vector<int> >(0));
        }
        if(newchain|| ( int(pRes->SSE)!=splineinfo.secstr_indices.back().back()[1] &&
          !( ( int(pRes->SSE)== SSE_Strand || int(pRes->SSE)== SSE_Bulge) &&
              (splineinfo.secstr_indices.back().back()[1] == SSE_Strand 
                 || splineinfo.secstr_indices.back().back()[1] == SSE_Bulge )) )  ) {
          splineinfo.secstr_indices.back().push_back(std::vector<int>(0));
          splineinfo.secstr_indices.back().back().push_back(nPept);
          if (pRes->SSE == SSE_Bulge) {
	    splineinfo.secstr_indices.back().back().push_back(int(SSE_Strand));
	  } else {
            splineinfo.secstr_indices.back().back().push_back(int(pRes->SSE));
          }
          newchain = 0;
        }
        // Save the CA pointer and aternate O pointers
	cavertices.back().push_back(pCA);
      }
  }


  int CA_nAtoms;
  PPCAtom CA_atomTable;
  // Find all CA - to use as quick check if atom is in this set
  for(unsigned ii=0; ii<cavertices.size(); ii++){
   prev_atoms.push_back(0);
   std::list<PCAtom> prev_atom_p;
   if(cavertices[ii].size()>1){
    
    int CA_this_chain_selHnd = molH->NewSelection();
    molH->Select(CA_this_chain_selHnd,STYPE_ATOM,0,cavertices[ii][0]->GetChainID(),
		 ANY_RES,"*",ANY_RES,"*","*","CA","C","*",SKEY_NEW);
    molH->ExcludeOverlappedAtoms(CA_this_chain_selHnd,0.8);
    //molH->Select(CA_this_chain_selHnd,STYPE_ATOM,atom_selHnd,SKEY_AND);

    molH->GetSelIndex ( CA_this_chain_selHnd, CA_atomTable, CA_nAtoms );
    //std::cout << "Selected " << CA_nAtoms << " CA atoms\n"; std::cout.flush();
    //std::cout << "Sequence number of first atom in table " << CA_atomTable[0]->GetResidueNo() << "\n";

    //std::cout << "Sequence number of first atom "  << cavertices[ii][0]->GetResidueNo() << "\n";
    PCAtom first_at = cavertices[ii][0];
    int resno = first_at->GetResidueNo();

    PCAtom curr_at = first_at;
    int curr_resno = resno;

    bool in_this_chain = true;
    //if(CA_atomTable[0]->GetResidueNo()>cavertices[ii][0]->GetResidueNo()){
    while(in_this_chain&&curr_resno>0&&curr_at->isInSelection(CA_this_chain_selHnd)){
      if(CA_atomTable[curr_resno-1]&&CA_atomTable[curr_resno]){
         PCAtom pCA1 = CA_atomTable[curr_resno];
         PCAtom pCA2 = CA_atomTable[curr_resno-1];
         dist2 = (pCA1->x - pCA2->x)*(pCA1->x - pCA2->x) +
                 (pCA1->y - pCA2->y)*(pCA1->y - pCA2->y)+
	         (pCA1->z - pCA2->z)*(pCA1->z - pCA2->z);
         if(dist2 > 20.0 || pCA1->GetChainID()!= pCA2->GetChainID()){
           in_this_chain = false;
         } else {
           curr_at = CA_atomTable[curr_resno-1];
           //std::cout << "Add CA of res " << curr_at->GetResidueNo() << " to beginning\n";
           curr_resno--; if(udd_CA<0) prev_atoms[ii]++;
           prev_atom_p.push_front(curr_at);
         }
      }
    }
   }
   if(prev_atoms[ii]>0){
    std::vector<PCAtom> new_ca0;
    std::list<PCAtom>::iterator patom_iter = prev_atom_p.begin();
    while(patom_iter!=prev_atom_p.end()){
      new_ca0.push_back(*patom_iter);
      patom_iter++;
    }
    new_ca0.insert(new_ca0.end(),cavertices[ii].begin(),cavertices[ii].end());
    cavertices[ii] = new_ca0;
   }
  }
  //}

  for(unsigned ii=0; ii<cavertices.size(); ii++){
   next_atoms.push_back(0);
   if(cavertices[ii].size()>1){
    int CA_this_chain_selHnd = molH->NewSelection();
    molH->Select(CA_this_chain_selHnd,STYPE_ATOM,0,cavertices[ii].back()->GetChainID(),
		 ANY_RES,"*",ANY_RES,"*","*","CA","C","*",SKEY_NEW);
    molH->ExcludeOverlappedAtoms(CA_this_chain_selHnd,0.8);
    //molH->Select(CA_this_chain_selHnd,STYPE_ATOM,atom_selHnd,SKEY_AND);

    molH->GetSelIndex ( CA_this_chain_selHnd, CA_atomTable, CA_nAtoms );
    //std::cout << "Sequence number of last atom " << cavertices[ii].back()->GetResidueNo() << "\n";
    PCAtom first_at = cavertices[ii].back();
    int resno = first_at->GetResidueNo();

    PCAtom curr_at = first_at;
    int curr_resno = resno;

    bool in_this_chain = true;
    while(in_this_chain&&curr_resno+1<CA_nAtoms&&curr_at->isInSelection(CA_this_chain_selHnd)){
      if(CA_atomTable[curr_resno+1]){
         PCAtom pCA1 = CA_atomTable[curr_resno];
         PCAtom pCA2 = CA_atomTable[curr_resno+1];
         dist2 = (pCA1->x - pCA2->x)*(pCA1->x - pCA2->x) +
                 (pCA1->y - pCA2->y)*(pCA1->y - pCA2->y)+
	         (pCA1->z - pCA2->z)*(pCA1->z - pCA2->z);
         if(dist2 > 20.0 || pCA1->GetChainID()!= pCA2->GetChainID()){
           in_this_chain = false;
         } else {
           curr_at = CA_atomTable[curr_resno+1];
           //std::cout << "Add CA of res " << curr_at->GetResidueNo() << " to end\n";
           curr_resno++; if(udd_CA<0) next_atoms[ii]++;
	   cavertices[ii].push_back(curr_at);
         }
      }
    }
   }
  }

  //std::cout << "prev.size(): " << prev_atoms.size() << std::endl;
  //std::cout << "next.size(): " << next_atoms.size() << std::endl;
  //std::cout << "prev[0]: " << prev_atoms[0] << std::endl;
  //std::cout << "next[0]: " << next_atoms[0] << std::endl;

  //std::cout << "splineinfo.colours.back().size(): " << splineinfo.colours.back().size() << "\n";

  //CNABasePairs bp(molH,atom_selHnd,atomTable,nAtoms,0);
  //std::vector<std::pair<PCResidue,PCResidue> > base_pairs = bp.GetPairs();

  splineinfo.nacolours.push_back(std::vector<Cartesian>(0));
  PCAtom c5prev = 0;
  for(int j=0;j<nAtoms;j++){
    PCResidue res = atomTable[j]->GetResidue();
    if(res){
      int restype = molH->GetRestypeCode(res);
      if(restype==RESTYPE_NUCL||restype==RESTYPE_DNA||restype==RESTYPE_RNA){
        PCAtom c5 = res->GetAtom("C5\'");
        if(!c5) c5 = res->GetAtom("C5*");
        //PCAtom c5 = res->GetAtom("P");
        if(nac5vertices.back().size()>0&&c5){
          PCAtom nac5prev = nac5vertices.back().back();
          double distnca = (AtToCart(c5) - AtToCart(nac5prev)).length();
          if(distnca> 9.5 ||(c5prev&&c5->GetChainID()!=c5prev->GetChainID())){
            PCResidue resp = c5prev->GetResidue();
            PCAtom c3p = resp->GetAtom("C3\'");
            if(!c3p) c3p = resp->GetAtom("C3*");
            if(!c3p) c3p = resp->GetAtom("C5\'");
            if(!c3p) c3p = resp->GetAtom("C5*");
	    nac5vertices.back().pop_back(); 
	    nac5vertices.back().push_back(c3p); 
            nac5vertices.push_back(std::vector<PCAtom>(0));
	    splineinfo.nacolours.push_back(std::vector<Cartesian>(0));
          }
        }
        if(c5){
          Cartesian cartc5 = AtToCart(c5);
          if(nac5vertices.back().size()>0&&(AtToCart(nac5vertices.back().back())-cartc5).length()>1e-6){
            c5prev = c5;
            nac5vertices.back().push_back(c5);
            if (atom_colour_vector){
              double * atcol = atom_colour_vector->GetRGB(j);
              splineinfo.nacolours.back().push_back(Cartesian(atcol));
	      delete [] atcol;
	    }
          }else if(nac5vertices.back().size()==0&&c5){
            c5prev = c5;
            nac5vertices.back().push_back(c5);
            if(atom_colour_vector){
              double * atcol = atom_colour_vector->GetRGB(j);
              splineinfo.nacolours.back().push_back(Cartesian(atcol));
	      delete [] atcol;
	    }
          }
        }
      }
    }
  }
  
  if(c5prev&&nac5vertices.size()>0&&nac5vertices.back().size()>0){
    PCAtom c3 = c5prev->GetResidue()->GetAtom("C3\'");
    if(!c3) c3 = c5prev->GetResidue()->GetAtom("C5*");
    if(c3){
      nac5vertices.back().pop_back(); 
      nac5vertices.back().push_back(c3); 
    }
  }

  //std::cout << nac5vertices.size() << "\n";
  //for(unsigned ii=0;ii<nac5vertices.size();ii++)
    //std::cout << "  " << ii << " " << nac5vertices[ii].size() << "\n";

  std::vector<std::vector<PCAtom> >::iterator mna=nac5vertices.begin();  
  while(mna!=nac5vertices.end()){
    std::vector<PCAtom>::iterator l=(*mna).begin();
    splineinfo.nasplines.push_back(std::vector<Cartesian>(0));
    splineinfo.n1_nasplines.push_back(std::vector<Cartesian>(0));
    splineinfo.n2_nasplines.push_back(std::vector<Cartesian>(0));
    if(mna->size()>2){
      std::vector <Cartesian> ctrlpoints;
      std::vector <Cartesian> n1_ctrlpoints;
      std::vector <Cartesian> spline;
      std::vector <Cartesian> n1_spline;
      std::vector <Cartesian> n2_spline;
      Cartesian posm1;
      Cartesian posm2;
      Cartesian posp1;
      int spaccu = spline_accu;
      int ii=0;
      while(l!=(*mna).end()){
        if(l<(*mna).end()-1)
          posp1 = AtToCart(*(l+1));
        if(l>=(*mna).begin()+1)
          posm1 = AtToCart(*(l-1));
        if(l>=(*mna).begin()+2)
          posm2 = AtToCart(*(l-2));
        pos= AtToCart(*l);
        PCAtom n1 = (*l)->GetResidue()->GetAtom("N1");
        PCAtom c2 = (*l)->GetResidue()->GetAtom("C2");
        PCAtom n3 = (*l)->GetResidue()->GetAtom("N3");
        PCAtom c4 = (*l)->GetResidue()->GetAtom("C4");
        PCAtom c5 = (*l)->GetResidue()->GetAtom("C5");
        PCAtom c6 = (*l)->GetResidue()->GetAtom("C6");
        PCAtom n9 = (*l)->GetResidue()->GetAtom("N9");
        bool reverse=false;
        if(n9&&strncmp((*l)->GetResidue()->name,"EDA",3))
          reverse = true;
	bool have_near = false;
	Cartesian near_atom;
	Cartesian at_to_near;
	if(n1&&c2&&n3&&c4&&c5&&c6){
	   have_near = true;
	   Cartesian n1_cart = AtToCart(n1);
	   Cartesian c2_cart = AtToCart(c2);
	   Cartesian n3_cart = AtToCart(n3);
	   Cartesian c4_cart = AtToCart(c4);
	   Cartesian c5_cart = AtToCart(c5);
	   Cartesian c6_cart = AtToCart(c6);
           Cartesian v1 = n1_cart - c2_cart;
           Cartesian v2 = n1_cart - c4_cart;
           v1.normalize();
           v2.normalize();
           Cartesian norm1 = Cartesian::CrossProduct(v1,v2); // Normal 1 not N1, Oh dear.
           if(strncmp((*l)->GetResidue()->name,"PSU",3)==0){
             PCAtom c1p = (*l)->GetResidue()->GetAtom("C1\'");
             if(!c1p) c1p = (*l)->GetResidue()->GetAtom("C1*");
             PCAtom c2p = (*l)->GetResidue()->GetAtom("C2\'");
             if(!c2p) c2p = (*l)->GetResidue()->GetAtom("C2*");
             Cartesian c1c2 = AtToCart(c1p) - AtToCart(c2p);
             if(Cartesian::DotProduct(c1c2,norm1)>0.0) reverse = true;
           }
           if(reverse)
             norm1 = -norm1;
           if(l<(*mna).end()-1&&l>=(*mna).begin()+1){
             Cartesian pmm = posp1 - posm1;
             pmm.normalize();
             Cartesian n2_temp = Cartesian::CrossProduct(pmm,norm1);
             norm1 = Cartesian::CrossProduct(pmm,n2_temp);
           } else if(l<(*mna).end()-1&&l==(*mna).begin()) {
             Cartesian pmm = posp1 - pos;
             pmm.normalize();
             Cartesian n2_temp = Cartesian::CrossProduct(pmm,norm1);
             norm1 = Cartesian::CrossProduct(pmm,n2_temp);
           } else if(l==(*mna).end()-1&&l>=(*mna).begin()+1) {
             Cartesian pmm = pos - posm1;
             pmm.normalize();
             Cartesian n2_temp = Cartesian::CrossProduct(pmm,norm1);
             norm1 = Cartesian::CrossProduct(pmm,n2_temp);
           }
           //if(n1_ctrlpoints.size()>0&&Cartesian::DotProduct(norm1,n1_ctrlpoints.back())<0.0)
           n1_ctrlpoints.push_back(norm1);
           ctrlpoints.push_back(pos);
	   //near_atom = (n1_cart + c2_cart + n3_cart + c4_cart + c5_cart + c6_cart)/6.0;
	   //at_to_near = near_atom-pos;
	   //at_to_near.normalize();
	}

        /*

        if(l<(*mna).begin()+1){
          ctrlpoints.push_back(pos);
        }else if(l>=(*mna).begin()+1&&l<(*mna).end()-1){

          Cartesian pmm = posm1-posp1;
          pmm.normalize();
          n1_ctrlpoints.push_back(Cartesian::CrossProduct(pmm,at_to_near));
          n1_ctrlpoints.back().normalize();
          n1_ctrlpoints.push_back(n1_ctrlpoints.back());

          ctrlpoints.push_back(pos);

        }else if(l==(*mna).end()-1){
          n1_ctrlpoints.push_back(n1_ctrlpoints.back());

          if(Cartesian::DotProduct(norm1,n1_ctrlpoints.back())<0.0)
            norm1 = -norm1;
          n1_ctrlpoints.push_back(norm1);

          n1_ctrlpoints.back().normalize();

          ctrlpoints.push_back(pos);
        }
        */
        l++; ii++;
      }
      spline = SplineCurve(ctrlpoints,(ctrlpoints.size()-1)*spaccu,3,1);
      n1_spline = SplineCurve(n1_ctrlpoints,(n1_ctrlpoints.size()-1)*spaccu,3,1);
      std::vector<Cartesian>::iterator spline_iter = spline.begin();
      std::vector<Cartesian>::iterator n1_iter = n1_spline.begin();
      while(spline_iter!=spline.end()-1&&n1_iter!=n1_spline.end()-1){
        Cartesian n2 = Cartesian::CrossProduct(*(spline_iter+1)-(*spline_iter),*n1_iter);
        n2.normalize();
        n2_spline.push_back(n2);
        spline_iter++; n1_iter++;
      }
      n2_spline.push_back(n2_spline.back());
      splineinfo.nasplines.back() = spline;
      splineinfo.n1_nasplines.back() = n1_spline;
      splineinfo.n2_nasplines.back() = n2_spline;
    }
    mna++;
  }
  if(splineinfo.nasplines.back().size()==0){
    splineinfo.nasplines.pop_back();
    splineinfo.n1_nasplines.pop_back();
    splineinfo.n2_nasplines.pop_back();
  }
  //std::cout << splineinfo.nasplines.size() << "\n";
  //for(unsigned jj=0;jj<splineinfo.nasplines.size();jj++){
    //std::cout << splineinfo.nasplines[jj].size() << "\n";
  //}

  // Loop over the chains of CA creating a vector of Cartesians of 
  // the control points which are the CA position.
  // Then create the spline for the chain.
  //std::cout << "Starting to build splines\n"; std::cout.flush();
  std::vector<std::vector<PCAtom> >::iterator m=cavertices.begin();  
  int ichain = 0;
  while(m!=cavertices.end()){
    //std::cout << "Starting to build a spline\n"; std::cout.flush();
    std::vector<PCAtom>::iterator l=(*m).begin();
    splineinfo.splines.push_back(std::vector<Cartesian>(0));
    splineinfo.n1_splines.push_back(std::vector<Cartesian>(0));
    splineinfo.n2_splines.push_back(std::vector<Cartesian>(0));
    if(m->size()>2){
      std::vector <Cartesian> ctrlpoints;
      std::vector <Cartesian> n1_ctrlpoints;
      std::vector <Cartesian> spline;
      std::vector <Cartesian> n1_spline;
      std::vector <Cartesian> n2_spline;

      bool is_beta=false;
      bool is_loop=false;
      bool is_helix=false;
      std::vector<std::vector<int> > beta_indices;
      std::vector<std::vector<int> > loop_indices;
      std::vector<std::vector<int> > helix_indices;
      int im = 0;
      while(l!=(*m).end()){
        pos= AtToCart(*l);
        /* With this we build a vector of start and end points of beta sheets in the current spline.*/
        bool curr_beta = (int((*l)->GetResidue()->SSE)==SSE_Strand || int((*l)->GetResidue()->SSE)==SSE_Bulge);
        if(l+1<(*m).end()) curr_beta = (int((*(l+1))->GetResidue()->SSE)==SSE_Strand || int((*(l+1))->GetResidue()->SSE)==SSE_Bulge)&&(int((*l)->GetResidue()->SSE)==SSE_Strand || int((*l)->GetResidue()->SSE)==SSE_Bulge);
        if(curr_beta&&!is_beta){
          std::vector<int> first_index(1);
          first_index[0] = im;
          beta_indices.push_back(first_index);
          is_beta = true;
        }else if(!curr_beta){
          if(is_beta)
            beta_indices.back().push_back(im);
          is_beta = false;
        }
        if(curr_beta&&l==(*m).end()-1){
          beta_indices.back().push_back(im);
        }
        /* With this we build a vector of start and end points of helices in the current spline.*/
        if(smooth_helix){
          bool curr_helix = (int((*l)->GetResidue()->SSE)==SSE_Helix );
          if(curr_helix&&!is_helix){
            std::vector<int> first_helix_index;
            first_helix_index.push_back(im);
            helix_indices.push_back(first_helix_index);
            is_helix = true;
          }else if(!curr_helix){
            if(is_helix){
              helix_indices.back().push_back(im);
            }
            is_helix = false;
          }
          if(curr_helix&&l==(*m).end()-1){
            helix_indices.back().push_back(im);
          }
        }
        if(flatten_loop){
          bool curr_loop = (int((*l)->GetResidue()->SSE)==SSE_None || (int((*l)->GetResidue()->SSE)!=SSE_Strand && int((*l)->GetResidue()->SSE)!=SSE_Bulge && int((*l)->GetResidue()->SSE)!=SSE_Helix) );
          if(curr_loop&&!is_loop){
            std::vector<int> first_loop_index;
            first_loop_index.push_back(im);
            loop_indices.push_back(first_loop_index);
            is_loop = true;
          }else if(!curr_loop){
            if(is_loop){
              loop_indices.back().push_back(im);
            }
            is_loop = false;
          }
          if(curr_loop&&l==(*m).end()-1){
            loop_indices.back().push_back(im);
          }
        }
        im++;
        /* 
         * Now we calculate the normals(n1) at a C alpha position. The normal of CA(n) is (CA(n)-CA(n-1)) X (CA(n+1)-CA(n)).
         * Additionally we calculate the normals(n2) of the normals and the bonds. 
         */
        if(l<(*m).begin()+2){
          ctrlpoints.push_back(pos);
        }else if(l==(*m).begin()+2){
          n1_ctrlpoints.push_back(Cartesian::CrossProduct(AtToCart(*(l-1))-AtToCart(*(l-2)),AtToCart(*(l))-AtToCart(*(l-1))));
          n1_ctrlpoints.back().normalize();
          n1_ctrlpoints.push_back(n1_ctrlpoints.back());
          ctrlpoints.push_back(pos);
        }else if((l>(*m).begin()+2)){
          Cartesian n1 = Cartesian::CrossProduct(AtToCart(*(l-1))-AtToCart(*(l-2)),AtToCart(*(l))-AtToCart(*(l-1)));
          n1.normalize();
          if(Cartesian::DotProduct(n1,n1_ctrlpoints.back())<0.0){
            if((l!=(*m).begin()&&(*l)->GetResidue()->SSE==SSE_Helix&&(*(l-1))->GetResidue()->SSE==SSE_Helix)){
              n1_ctrlpoints.push_back(n1);
            }
            else
              n1_ctrlpoints.push_back(-n1);
          }else
            n1_ctrlpoints.push_back(n1);
          n1_ctrlpoints.back().normalize();
          ctrlpoints.push_back(pos);
        }
        if(l==(*m).end()-1){
          n1_ctrlpoints.push_back(n1_ctrlpoints.back());
          n1_ctrlpoints.back().normalize();
        }
        l++;
      }

      spline = SplineCurve(ctrlpoints,(ctrlpoints.size()-1)*spline_accu,3,1);
      //std::cout << ctrlpoints.size() << " control points gives " << spline.size() << " spline points\n";
      n1_spline = SplineCurve(n1_ctrlpoints,(n1_ctrlpoints.size()-1)*spline_accu,3,1);

      /* Now we fit the CA positions in beta-sheets and the normals to them to smooth curves. */
      //std::cout << "Before flattening beta sheets " << ctrlpoints.size() << " " << n1_ctrlpoints.size() << "\n";
      std::vector<Cartesian> spline_smoothed = spline;
      for(unsigned ii=0; ii<beta_indices.size(); ii++){
        int start_beta = beta_indices[ii][0];
        int end_beta   = beta_indices[ii][1];
        if(end_beta<int(ctrlpoints.size())-1);
           end_beta++;

        if(end_beta-start_beta>1){
         /* Take the section of the whole CA trace which is this beta sheet. */
         std::vector<Cartesian> beta_ctrlpoints(ctrlpoints.begin()+start_beta,ctrlpoints.begin()+end_beta);
         std::vector<Cartesian> beta_n1_ctrlpoints(n1_ctrlpoints.begin()+start_beta,n1_ctrlpoints.begin()+end_beta);
         /* Now calculate a Bezier curve through these points. */
         std::vector<Cartesian> beta_spline = BezierCurve(beta_ctrlpoints,spline_accu);
         std::vector<Cartesian> beta_n1_spline = BezierCurve(beta_n1_ctrlpoints,spline_accu);
         //std::cout << end_beta-start_beta << ": " << beta_spline.size() << "\n";
         beta_spline.pop_back();
         beta_n1_spline.pop_back();
         /* Then splice the new points back in.*/
         if(flatten_beta_sheet){
           Replace(spline,beta_spline,start_beta*spline_accu,(end_beta-1)*spline_accu);
           Replace(n1_spline,beta_n1_spline,start_beta*spline_accu,(end_beta-1)*spline_accu);
         }

	 if(!flatten_beta_sheet&&beta_spline.size()>20){
	   for(unsigned ii=0;ii<11&&ii;ii++){
             double frac = double(ii)/10.0;
             spline_smoothed[start_beta*spline_accu+ii] = (1.0-frac) * beta_spline[ii] + frac * spline_smoothed[start_beta*spline_accu+ii];
           }
	   for(unsigned ii=0;ii<11&&ii;ii++){
             double frac = double(ii)/10.0;
             spline_smoothed[end_beta*spline_accu-1-ii] = frac * beta_spline[beta_spline.size()-1-ii] + (1.0-frac) * spline_smoothed[end_beta*spline_accu-1-ii];
           }
         } else {
           Replace(spline_smoothed,beta_spline,start_beta*spline_accu,(end_beta-1)*spline_accu);
         }

         if(next_atoms[ichain]>0&&(end_beta-1)*spline_accu<(int)spline.size())spline[(end_beta-1)*spline_accu]=ctrlpoints[end_beta-1];
         if(next_atoms[ichain]>0&&(end_beta-1)*spline_accu<(int)spline.size())n1_spline[(end_beta-1)*spline_accu]=n1_ctrlpoints[end_beta-1];
         if(prev_atoms[ichain]>0)spline[(start_beta)*spline_accu]=ctrlpoints[start_beta];
         if(prev_atoms[ichain]>0)n1_spline[(start_beta)*spline_accu]=n1_ctrlpoints[start_beta];
        }
      }

      //std::cout << "Before flattening loops " << ctrlpoints.size() << " " << n1_ctrlpoints.size() << "\n";
      if(flatten_loop){
        for(unsigned ii=0; ii<loop_indices.size(); ii++){
          int start_loop = loop_indices[ii][0];
          int end_loop   = loop_indices[ii][1];
          if(end_loop<int(ctrlpoints.size())-1);
             end_loop++;
  
          if(end_loop-start_loop>2){
           std::vector<Cartesian> loop_ctrlpoints(ctrlpoints.begin()+start_loop,ctrlpoints.begin()+end_loop);
           std::vector<Cartesian> loop_n1_ctrlpoints(n1_ctrlpoints.begin()+start_loop,n1_ctrlpoints.begin()+end_loop);

           std::vector<Cartesian> loop_spline = BezierCurve(loop_ctrlpoints,spline_accu);
	   if(start_loop==0||end_loop==(int)ctrlpoints.size())
             Replace(spline,loop_spline,start_loop*spline_accu,(end_loop-1)*spline_accu);
	   else
             ReplaceSmoothly(spline,loop_spline,start_loop*spline_accu,(end_loop-1)*spline_accu);

           for(unsigned jj=1;jj<loop_n1_ctrlpoints.size()&&jj+start_loop<ctrlpoints.size()-1;jj++){
             Cartesian pp1 = spline[(start_loop+jj+1)*spline_accu-1];
             int low_index = (start_loop+jj-1)*spline_accu-1;
	     if(low_index<0) low_index = 0;
             Cartesian pm1 = spline[low_index];
             pp1.normalize();
	     /*
	     if(pp1.length()<1e-5) {
	       std::cout << "pp1: " << jj << " " << start_loop << " " << end_loop << "\n";
	       std::cout << "pp1: " << (start_loop+jj+1)*spline_accu << " "  << spline.size()<< "\n";
	     }
	     */
             pm1.normalize();
	     /*
	     if(pm1.length()<1e-5) {
	       std::cout << "pm1: " << jj << " " << start_loop << " " << end_loop << "\n";
	       std::cout << "pm1: " << int((start_loop+jj-1)*spline_accu) << " "  << spline.size()<< "\n";
	     }
	     */
             loop_n1_ctrlpoints[jj] = Cartesian::CrossProduct(pm1,pp1);
	     //if(loop_n1_ctrlpoints[jj].length()<1e-5) std::cout << "n1: " << jj << " " << start_loop << " " << end_loop << "\n";
             loop_n1_ctrlpoints[jj].normalize();
             if(Cartesian::DotProduct(loop_n1_ctrlpoints[jj],loop_n1_ctrlpoints[jj-1])<0.0)
               loop_n1_ctrlpoints[jj] = -loop_n1_ctrlpoints[jj];
           }
           
           std::vector<Cartesian> loop_n1_spline = BezierCurve(loop_n1_ctrlpoints,spline_accu);
           Replace(n1_spline,loop_n1_spline,start_loop*spline_accu,(end_loop-1)*spline_accu);

           Replace(spline_smoothed,loop_spline,start_loop*spline_accu,(end_loop-1)*spline_accu);

           if(next_atoms[ichain]>0&&(end_loop-1)*spline_accu<(int)spline.size())spline[(end_loop-1)*spline_accu]=ctrlpoints[end_loop-1];
           if(next_atoms[ichain]>0&&(end_loop-1)*spline_accu<(int)spline.size())n1_spline[(end_loop-1)*spline_accu]=n1_ctrlpoints[end_loop-1];
           if(prev_atoms[ichain]>0)spline[(start_loop)*spline_accu]=ctrlpoints[start_loop];
           if(prev_atoms[ichain]>0)n1_spline[(start_loop)*spline_accu]=n1_ctrlpoints[start_loop];
          }
        }
      }
      //std::cout << "After flattening " << ctrlpoints.size() << " " << n1_ctrlpoints.size() << "\n";
      //std::cout << "Before smoothing helices " << ctrlpoints.size() << " " << n1_ctrlpoints.size() << "\n"; std::cout.flush();
      //std::cout << "Before smoothing helices(spline) " << spline.size() << " " << n1_spline.size() << "\n"; std::cout.flush();
      if(smooth_helix){
        for(unsigned ii=0; ii<helix_indices.size(); ii++){
          int start_helix = helix_indices[ii][0];
          int end_helix   = helix_indices[ii][1];
          //if(next_atoms[ichain]>0&&end_helix<int(ctrlpoints.size())-1) {
            //end_helix++;
          //}
          //if(prev_atoms[ichain]>0&&start_helix>0) start_helix--;
          //if(end_helix<int(ctrlpoints.size())-1);
             //end_helix++;
  
          if(end_helix-start_helix>2){
           std::vector<Cartesian> helix_ctrlpoints(ctrlpoints.begin()+start_helix,ctrlpoints.begin()+end_helix);
           std::vector<Cartesian> helix_n1_ctrlpoints(n1_ctrlpoints.begin()+start_helix,n1_ctrlpoints.begin()+end_helix);

           std::vector<Cartesian> helix_spline;
           helix_spline = BezierCurve(helix_ctrlpoints,spline_accu);
           bool is_straight=false;
           if(helix_ctrlpoints.size()<16){
             is_straight = true;
           }else{
             int half = (end_helix-start_helix)/2;
             std::vector<Cartesian> pca1 = Cartesian::PrincipalComponentAnalysis(ctrlpoints.begin()+start_helix,ctrlpoints.begin()+start_helix+half);
             std::vector<Cartesian> pca2 = Cartesian::PrincipalComponentAnalysis(ctrlpoints.begin()+start_helix+half,ctrlpoints.begin()+end_helix);
             double olap = Cartesian::DotProduct(pca1[0],pca2[0]);
             if(fabs(olap)>0.7) is_straight = true;
           }
           if(is_straight){ 
             std::vector<Cartesian> old_helix_spline = helix_spline;
             helix_spline.clear();
             std::vector<Cartesian> pca = Cartesian::PrincipalComponentAnalysis(ctrlpoints.begin()+start_helix,ctrlpoints.begin()+end_helix);
             Cartesian mean = Mean(helix_ctrlpoints);
             double length = (helix_ctrlpoints[0] - helix_ctrlpoints[helix_ctrlpoints.size()-1]).length();
             Cartesian orig_dir = helix_ctrlpoints[0] - helix_ctrlpoints[helix_ctrlpoints.size()-1];
             orig_dir.normalize();
             Cartesian new_dir = pca[0];
             new_dir.normalize();
             Cartesian start;
             Cartesian end;
             if(Cartesian::DotProduct(orig_dir,new_dir)>0.0){
               start =  mean + 0.5*length*pca[0];
               end =  mean - 0.5*length*pca[0];
             }else{
               start =  mean - 0.5*length*pca[0];
               end =  mean + 0.5*length*pca[0];
               new_dir = -new_dir;
             }
             
             for(unsigned ifrac=0;ifrac<(helix_ctrlpoints.size()-1)*spline_accu+2;ifrac++){
               double frac = double(ifrac)/((helix_ctrlpoints.size()-1)*spline_accu+1);
               helix_spline.push_back((1.0-frac)*start + frac*end);
             }

             std::vector<Cartesian> helix_n1_spline = BezierCurve(helix_n1_ctrlpoints,spline_accu);
             std::vector<Cartesian> old_helix_n1_spline = helix_n1_spline;
             helix_spline = std::vector<Cartesian>(helix_spline.begin()+1,helix_spline.end());
             int ismooth;
             for(ismooth=0;ismooth<(int)spline_accu-2;ismooth++){
                double frac = sin(double(ismooth)/(spline_accu-1)*HALF_PI);
                helix_spline[ismooth] = frac*helix_spline[ismooth] + (1.0-frac) * old_helix_spline[ismooth];
             }
             for(ismooth=helix_spline.size()-spline_accu+2;ismooth<(int)helix_spline.size();ismooth++){
                double frac = sin(double(helix_spline.size()-ismooth-1)/(spline_accu-1)*HALF_PI);
                helix_spline[ismooth] = frac*  helix_spline[ismooth] +  (1.0-frac)* old_helix_spline[ismooth];
             }

             Replace(spline,helix_spline,start_helix*spline_accu,(end_helix-1)*spline_accu);

             Cartesian ntry = helix_n1_spline[helix_n1_spline.size()/2];
             ntry.normalize();
             Cartesian ntryp = Cartesian::CrossProduct(new_dir,ntry);
             ntryp.normalize();
             ntry = Cartesian::CrossProduct(ntryp,new_dir);
             for(unsigned in1=0;in1<helix_n1_spline.size();in1++){
               helix_n1_spline[in1] = ntry;
             }
             for(ismooth=0;ismooth<spline_accu-2;ismooth++){
                double frac = sin(double(ismooth)/(spline_accu-1)*HALF_PI);
                helix_n1_spline[ismooth] = frac*helix_n1_spline[ismooth] + (1.0-frac) * old_helix_n1_spline[ismooth];
                if(Cartesian::DotProduct(helix_n1_spline[ismooth],helix_n1_spline[ismooth+1])<0.0)
                  helix_n1_spline[ismooth+1] = -helix_n1_spline[ismooth+1];
             }
             for(ismooth=helix_n1_spline.size()-spline_accu+2;ismooth<(int)helix_n1_spline.size();ismooth++){
                double frac = sin(double(helix_n1_spline.size()-ismooth-1)/(spline_accu-1)*HALF_PI);
                helix_n1_spline[ismooth] = frac*  helix_n1_spline[ismooth] +  (1.0-frac)* old_helix_n1_spline[ismooth];
                if(Cartesian::DotProduct(helix_n1_spline[ismooth],helix_n1_spline[ismooth-1])<0.0)
                  helix_n1_spline[ismooth] = -helix_n1_spline[ismooth];
             }

             if(end_helix<int(ctrlpoints.size())){
               helix_n1_spline.push_back(helix_n1_spline.back()); // Pointless??
               helix_n1_spline.push_back(helix_n1_spline.back());
             }
             Replace(n1_spline,helix_n1_spline,start_helix*spline_accu,(end_helix-1)*spline_accu);
             Replace(spline_smoothed,helix_spline,start_helix*spline_accu,(end_helix-1)*spline_accu);
           } else {
             helix_spline = BezierCurve(helix_ctrlpoints,spline_accu);
             std::vector<Cartesian> old_helix_spline = helix_spline;
             int ismooth;
             for(ismooth=0;ismooth<(int)spline_accu-2;ismooth++){
                double frac = sin(double(ismooth)/(spline_accu-1)*HALF_PI);
                helix_spline[ismooth] = frac*helix_spline[ismooth] + (1.0-frac) * old_helix_spline[ismooth];
             }
             for(ismooth=helix_spline.size()-spline_accu+2;ismooth<(int)helix_spline.size();ismooth++){
                double frac = sin(double(helix_spline.size()-ismooth-1)/(spline_accu-1)*HALF_PI);
                helix_spline[ismooth] = frac*  helix_spline[ismooth] +  (1.0-frac)* old_helix_spline[ismooth];
             }
             for(unsigned jj=1;jj<helix_n1_ctrlpoints.size()-1;jj++){
               Cartesian p = helix_spline[jj*spline_accu];
               Cartesian pp1 = helix_spline[(jj+1)*spline_accu];
               Cartesian pm1 = helix_spline[(jj-1)*spline_accu];
               p.normalize();
               pp1.normalize();
               pm1.normalize();
               helix_n1_ctrlpoints[jj] = Cartesian::CrossProduct(p-pm1,pp1-p);
               helix_n1_ctrlpoints[jj].normalize();
               if(Cartesian::DotProduct(helix_n1_ctrlpoints[jj],helix_n1_ctrlpoints[jj-1])<0.0)
                 helix_n1_ctrlpoints[jj] = -helix_n1_ctrlpoints[jj];
             }
             helix_n1_ctrlpoints[0] = helix_n1_ctrlpoints[1];
             std::vector<Cartesian> helix_n1_spline = BezierCurve(helix_n1_ctrlpoints,spline_accu);
             std::cout << helix_n1_spline.size() << "\n";
             helix_n1_spline[0] = n1_spline[spline_accu];
             for(unsigned jj=1;jj<helix_n1_spline.size();jj++){
               Cartesian p = helix_spline[jj];
               Cartesian pm1 = helix_spline[jj-1];
               Cartesian pp1 = helix_spline[jj+1];
               helix_n1_spline[jj] = Cartesian::CrossProduct(p-pm1,pp1-p);
               if(Cartesian::DotProduct(helix_n1_spline[jj],helix_n1_spline[jj-1])<0.0)
                 helix_n1_spline[jj] = -helix_n1_spline[jj];
             }
             
             Replace(spline,helix_spline,start_helix*spline_accu,(end_helix-1)*spline_accu);
             Replace(n1_spline,helix_n1_spline,start_helix*spline_accu,(end_helix-1)*spline_accu);
             Replace(spline_smoothed,helix_spline,start_helix*spline_accu,(end_helix-1)*spline_accu);
           }
          }
        }
      }
      //std::cout << "After flattening " << ctrlpoints.size() << " " << n1_ctrlpoints.size() << "\n"; std::cout.flush();
      //std::cout << "After flattening(spline) " << spline.size() << " " << n1_spline.size() << "\n"; std::cout.flush();

      //std::cout << ctrlpoints.size() << " " << n1_ctrlpoints.size() << "\n";

      std::vector<Cartesian>::iterator spline_smoothed_iter = spline_smoothed.begin();
      std::vector<Cartesian>::iterator spline_iter = spline.begin();
      std::vector<Cartesian>::iterator n1_iter = n1_spline.begin();
      //std::cout << "Calculate n2\n"; std::cout.flush();
      if(spline.size()>0){
        int jj=0;
        while(spline_smoothed_iter!=spline_smoothed.end()-1&&spline_iter!=spline.end()-1&&n1_iter!=n1_spline.end()-1){
          n1_iter->normalize();
           // Shortening loops may cause there to be a 180 deg change in normals at start of next
	   // secondary structure element. So we fix that here.
          if(Cartesian::DotProduct(*(n1_iter),*(n1_iter+1))<0.0) *(n1_iter+1) = - *(n1_iter+1);
          Cartesian n2;
          n2 = Cartesian::CrossProduct(*(spline_smoothed_iter+1)-(*spline_smoothed_iter),*n1_iter);
          n2.normalize();
          /*
          if(n2_spline.size()>0&&smooth_helix){ //Helix smooth is the most destructive....
            double olap = Cartesian::DotProduct(n2,n2_spline.back());
            if(fabs(olap)<0.3){
              Cartesian spmmp = Cartesian::MidPoint(*(spline_iter-1),*(spline_iter+1));
              *(spline_iter) = Cartesian::MidPoint(*(spline_iter),spmmp);
              Cartesian pmmp = Cartesian::MidPoint(*(n1_iter-1),*(n1_iter+1));
              *(n1_iter) = Cartesian::MidPoint(*(n1_iter),pmmp);
              Cartesian n2mp=Cartesian::MidPoint(n2,n2_spline.back());
              n2 = n2_spline.back() = n2mp;
            }
          }
          */
          n2_spline.push_back(n2);
          spline_smoothed_iter++; spline_iter++; n1_iter++; jj++;
        }
        n2_spline.push_back(n2_spline.back());
      }

      //std::cout << "Spline points: " << spline.size() << "\n";
      //std::cout << "First spline point: " << spline[0] << "\n";
      if(next_atoms[ichain]>0&&int(ctrlpoints.size())>next_atoms[ichain]+1){
         spline.resize((ctrlpoints.size()-1-next_atoms[ichain])*spline_accu+1);
         //std::cout << ctrlpoints[ctrlpoints.size()-1-next_atoms[ichain]] << "\n";
         n1_spline.resize((ctrlpoints.size()-1-next_atoms[ichain])*spline_accu+1);
         n2_spline.resize((ctrlpoints.size()-1-next_atoms[ichain])*spline_accu+1);
         //std::cout << "This(End): " << ctrlpoints[ctrlpoints.size()-1-next_atoms[ichain]] << "\n";
         //std::cout << "Previous(End): " << ctrlpoints[ctrlpoints.size()-2-next_atoms[ichain]] << "\n";
         //std::cout << "Next(Begin): " << ctrlpoints[ctrlpoints.size()-next_atoms[ichain]] << "\n";
         //std::cout << "This(N): " << n1_ctrlpoints[ctrlpoints.size()-1-next_atoms[ichain]] << "\n";
         //std::cout << "Next(N): " << n1_ctrlpoints[ctrlpoints.size()-next_atoms[ichain]] << "\n";
         if((!flatten_loop)&&ctrlpoints.size()>3){
           spline[spline.size()-1] = ctrlpoints[ctrlpoints.size()-1-next_atoms[ichain]];
           n1_spline[spline.size()-1] = n1_ctrlpoints[ctrlpoints.size()-1-next_atoms[ichain]];
           if(Cartesian::DotProduct(n1_spline[spline.size()-2],n1_spline[spline.size()-1])<0.0)
             n1_spline[spline.size()-1] = -n1_spline[spline.size()-1];
	   Cartesian diff = ctrlpoints[ctrlpoints.size()-1-next_atoms[ichain]+1]-ctrlpoints[ctrlpoints.size()-2-next_atoms[ichain]];
           diff.normalize();
           n2_spline[spline.size()-1] = Cartesian::CrossProduct(diff,n1_spline[spline.size()-1]);
           if(Cartesian::DotProduct(n2_spline[spline.size()-2],n2_spline[spline.size()-1])<0.0)
             n2_spline[spline.size()-1] = -n2_spline[spline.size()-1];
         }
         //std::cout << "End: " << n1_spline[spline.size()-1] << "\n";

      }
      //std::cout << "Spline(2) points: " << spline.size() << "\n";
      //std::cout << "First spline point: " << spline[0] << "\n";
      if(prev_atoms[ichain]>0&&(int)spline.size()>=prev_atoms[ichain]*spline_accu&&(int)n1_spline.size()>=prev_atoms[ichain]*spline_accu){
         if(flatten_loop){
           int extraCut = 1;
           int prev_idx = prev_atoms[ichain];
           // Are we in a flattenned loop?
           for(unsigned iloop=0; iloop<loop_indices.size(); iloop++){
             int start_loop = loop_indices[iloop][0];
             int end_loop   = loop_indices[iloop][1];
             if(prev_idx>=start_loop&&prev_idx<=end_loop){
               //std::cout << prev_idx << " in loop " << start_loop << ", " << end_loop << "\n";
               extraCut = 0;
               break;
             }
           }
           spline.erase(spline.begin(),spline.begin()+prev_atoms[ichain]*spline_accu+extraCut);
           n1_spline.erase(n1_spline.begin(),n1_spline.begin()+prev_atoms[ichain]*spline_accu+extraCut);
           n2_spline.erase(n2_spline.begin(),n2_spline.begin()+prev_atoms[ichain]*spline_accu+extraCut);
         }else{
           spline.erase(spline.begin(),spline.begin()+prev_atoms[ichain]*spline_accu+1);
           n1_spline.erase(n1_spline.begin(),n1_spline.begin()+prev_atoms[ichain]*spline_accu+1);
           n2_spline.erase(n2_spline.begin(),n2_spline.begin()+prev_atoms[ichain]*spline_accu+1);
         }
         //std::cout << "This(Begin): " << ctrlpoints[prev_atoms[ichain]] << "\n";
         //std::cout << "Previous(Begin): " << ctrlpoints[prev_atoms[ichain]-1] << "\n";
         //std::cout << "This(N): " << n1_ctrlpoints[prev_atoms[ichain]] << "\n";
         if((!flatten_loop)&&ctrlpoints.size()>3&&prev_atoms[ichain]+1<ctrlpoints.size()){
           Cartesian diff = ctrlpoints[prev_atoms[ichain]+1] - ctrlpoints[prev_atoms[ichain]-1];
           n1_spline[0] = n1_ctrlpoints[prev_atoms[ichain]];
           if(Cartesian::DotProduct(n1_spline[0],n1_spline[1])<0.0)
             n1_spline[0] = -n1_spline[0];
           diff.normalize();
           n2_spline[0] = Cartesian::CrossProduct(diff,n1_spline[0]);
           if(Cartesian::DotProduct(n2_spline[0],n2_spline[1])<0.0)
             n2_spline[0] = -n2_spline[0];
         }
         //std::cout << "Previous(N): " << n1_ctrlpoints[prev_atoms[ichain]-1] << "\n";
         //std::cout << "Begin: " << n1_spline[0] << "\n";
      }
      //std::cout << "Spline(3) points: " << spline.size() << "\n";
      //std::cout << "First spline point: " << spline[0] << "\n";
      //std::cout << "Last spline point: " << spline[spline.size()-1] << "\n";

      //std::cout << spline.size() << " " << n1_spline.size() << " " << n2_spline.size() << "\n";
      //std::cout << "first spline point: " <<  spline[0] << "\n";
      //std::cout << "Calculated n2\n"; std::cout.flush();
 
      // LP ?
      // Stick the last of cavertices coords on the end of each spline 
      // so there is something there for GetExtCartesians to join the 
      // side chains to.  
      // NO! -- this causes problems somewhere else 
      // probably because not all vectors are the same length
      // SJM 19/08/05
      // YES? - If we extend n1 and n2 as well? 
      /*
      if (spline.size()>0 && (*m).size()>1 && next_atoms[ichain]<1){
        spline.push_back(AtToCart((*m)[(*m).size()-1]));
        n1_spline.push_back(n1_spline.back());
        n2_spline.push_back(n2_spline.back());
      }*/
      //std::cout << "Last spline point: " << spline[spline.size()-1] << "\n";
      //std::cout << "First spline point: " << spline[0] << "\n";
      //std::cout << "First control point: " << AtToCart((*m)[prev_atoms[ichain]]) << "\n";
      if(spline.size()>0&&n1_spline.size()>0&&n2_spline.size()>0&&(spline[0]-AtToCart((*m)[prev_atoms[ichain]])).length()>1e-3){
        //std::cout << "Inserting at beginning\n";
        if(!flatten_loop){
        spline.insert(spline.begin(),AtToCart((*m)[prev_atoms[ichain]]));
        n1_spline.insert(n1_spline.begin(),n1_spline[0]);
        n2_spline.insert(n2_spline.begin(),n2_spline[0]);
        }
      }
      //std::cout << "Extended\n"; std::cout.flush();

      splineinfo.splines.back() = spline;
      splineinfo.n1_splines.back() = n1_spline;
      splineinfo.n2_splines.back() = n2_spline;
      if(ctrlpoints.size()>0){
      for(unsigned int ii=0;ii<ctrlpoints.size()-1;ii++){
	double length = (ctrlpoints[ii+1]-ctrlpoints[ii]).length();
	if(length>15.0){
          std::cout << "Potentially bad CA-CA distance" << length << "\n";
          std::cout << "Most likely due to beta-sheet flattening" << length << "\n";
        }
      }
      }
      //std::cout << "Done\n"; std::cout.flush();
    }
    m++; ichain++;
  }


  while(splineinfo.splines.size()>0&&splineinfo.splines.back().size()==0){
    splineinfo.splines.pop_back();
    splineinfo.n1_splines.pop_back();
    splineinfo.n2_splines.pop_back();
    splineinfo.colours.pop_back();
  }
  while(splineinfo.nasplines.size()>0&&splineinfo.nasplines.back().size()==0){
    splineinfo.nasplines.pop_back();
    splineinfo.n1_nasplines.pop_back();
    splineinfo.n2_nasplines.pop_back();
    splineinfo.nacolours.pop_back();
  }

  // Clean up the selection handle 
  molH->DeleteSelection(CALocAselHnd);
  molH->DeleteSelection(CALocBselHnd);
  molH->DeleteSelection(CALocCselHnd);
  molH->DeleteSelection(CAselHnd);
  molH->DeleteSelection(atom_selHnd);

  std::vector<std::vector<std::vector<int> > > secstr_indices = splineinfo.secstr_indices;
  std::vector<std::vector<std::vector<int> > > secstr_indices_new;

  // Now we subtract one from secondary structure "end" of beta strands and helices.
  // If a strand/helix runs straight into another we insert a bit of "no-structure".
  for(unsigned ichain=0;ichain<secstr_indices.size();ichain++){
    secstr_indices_new.push_back(std::vector<std::vector<int> >(0));
    std::vector<std::vector<int> >::iterator sec_ind_iter=secstr_indices[ichain].begin();
    while(sec_ind_iter!=secstr_indices[ichain].end()){
      if(sec_ind_iter!=secstr_indices[ichain].end()-1){
	secstr_indices_new.back().push_back(*sec_ind_iter);
        if((*sec_ind_iter)[1]==SSE_Strand||(*sec_ind_iter)[1]==SSE_Bulge||(*sec_ind_iter)[1]==SSE_Helix){
          if((*(sec_ind_iter+1))[1]==SSE_Strand||(*(sec_ind_iter+1))[1]==SSE_Bulge||(*(sec_ind_iter+1))[1]==SSE_Helix){
            std::vector<int> empty;
            empty.push_back((*(sec_ind_iter+1))[0]-1);
            empty.push_back(0);
            secstr_indices_new.back().push_back(empty);
          } else {
            (*(sec_ind_iter+1))[0] -= 1;
	  }
        }
      }else {
	secstr_indices_new.back().push_back(*sec_ind_iter);
      }
      sec_ind_iter++;
    }
  }

  splineinfo.secstr_indices = secstr_indices_new;
 
  return splineinfo;

}

std::vector<std::vector<Cartesian> > GetExternalCartesians(PCMMDBManager molhnd, const std::vector<std::vector<int> > &ext_conn_lists, int side_to_ribbon, int side_to_worm ){

  PPCAtom atomTable=0;
  int nAtoms;
  molhnd->GetAtomTable ( atomTable, nAtoms );  
  std::vector<std::vector<Cartesian> > carts;
  PCAtom atom;
  int udd_chain,udd_CA,idx,ids;
  int spline_accu = 4;
  
  // side_to_ribbon is a selection handle for an different object drawn 
  // as a spline.  We want any connection to CA atoms in that spline to
  // use the coordinates of the spline and not the 'true' CA coords
  if (side_to_ribbon>0||side_to_worm>0) {
    // Create a spline through all atoms 
    // Set up UDD for GetSplineInfo to load with index to spline data
    // for each selected residue
    udd_chain = molhnd->GetUDDHandle ( UDR_ATOM,"tmp_atom_int" );
    if (udd_chain<=0) udd_chain = molhnd->RegisterUDInteger(UDR_ATOM,"tmp_atom_int" );
    udd_CA = molhnd->GetUDDHandle ( UDR_ATOM,"tmp_atom_int1" );
    if (udd_CA<=0) udd_CA = molhnd->RegisterUDInteger(UDR_ATOM,"tmp_atom_int1" );
    for (int i=0;i<nAtoms;i++)atomTable[i]->PutUDData(udd_CA,-1);

    SplineInfo splineinfo;
    PPCAtom SelAtoms=0;
    PCAtom pCA; 
    int nSelAtoms=0;
    molhnd->GetSelIndex(side_to_ribbon,SelAtoms,nSelAtoms);
    bool isworm=false;
    if (side_to_ribbon>0&&nSelAtoms>0){
      splineinfo = GetSplineInfo ((PCMMANManager)molhnd, side_to_ribbon ,NULL,spline_accu, udd_chain, udd_CA, 1 );
    }else{
      splineinfo = GetSplineInfo ((PCMMANManager)molhnd, side_to_worm ,NULL,spline_accu, udd_chain, udd_CA, 0 );
      isworm=true;
    }
    //if(SelAtoms) delete [] SelAtoms;
    
    std::vector<std::vector<int> >::const_iterator ext_iter = ext_conn_lists.begin();
    while(ext_iter!=ext_conn_lists.end()){
      std::vector<int>::const_iterator this_ext_iter=ext_iter->begin();
      carts.push_back(std::vector<Cartesian>());
      while(this_ext_iter!=ext_iter->end()){
        atom = atomTable[*this_ext_iter];
        pCA =  atom->GetResidue()->GetAtom("CA");
        if((!strncmp((atom)->GetResidue()->name,"PRO",3)) && pCA != NULL &&
             pCA->GetUDData(udd_CA,idx) == UDDATA_Ok && 
             pCA->GetUDData(udd_chain,ids) == UDDATA_Ok && 
             idx >= 0 && 
             (!strncmp(atom->name," N",2)) && 
             unsigned(ids) < splineinfo.splines.size() &&
             unsigned(idx * spline_accu - spline_accu/2) < (splineinfo.splines[ids].size())){
          idx = idx * spline_accu - spline_accu/2; // This is of course an approximation.
          carts.back().push_back(Cartesian(splineinfo.splines[ids][idx].get_x(),
					   splineinfo.splines[ids][idx].get_y(),
                                           splineinfo.splines[ids][idx].get_z() ));
        }else if ((!isworm) && 
                   atom->GetUDData(udd_CA,idx) == UDDATA_Ok && 
                   idx >= 0 && 
                   atom->GetUDData(udd_chain,ids) == UDDATA_Ok &&
                   unsigned(ids) < splineinfo.splines.size() &&
                   unsigned(idx * spline_accu) < (splineinfo.splines[ids].size()) ) {
          // Get the coords out of the SplineInfo
          idx = idx * spline_accu;
          carts.back().push_back(Cartesian(splineinfo.splines[ids][idx].get_x(),
					   splineinfo.splines[ids][idx].get_y(),
                                           splineinfo.splines[ids][idx].get_z() ));
          if((carts.back().back()-Cartesian(atom->x,atom->y,atom->z)).length()<0.17)
            carts.back().back() = Cartesian(atom->x,atom->y,atom->z);
        } else {
          carts.back().push_back(Cartesian(atom->x,atom->y,atom->z));
        }
        this_ext_iter++;
      }
      ext_iter++;
    }

  } else {

    std::vector<std::vector<int> >::const_iterator ext_iter = ext_conn_lists.begin();
    while(ext_iter!=ext_conn_lists.end()){
      std::vector<int>::const_iterator this_ext_iter=ext_iter->begin();
      carts.push_back(std::vector<Cartesian>());
      while(this_ext_iter!=ext_iter->end()){
        atom = atomTable[*this_ext_iter];
        carts.back().push_back(Cartesian(atom->x,atom->y,atom->z));
        this_ext_iter++;
      }
      ext_iter++;
    }
  }
  return carts;
}

