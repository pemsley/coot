/*
     mmut/mmut_sasarea.cc: CCP4MG Molecular Graphics Program
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

Calculate Solvent Accessible Surface Area by methods of:
  Surface Points
  Lee and Richards
  Wodak and Janin

Converted from original CCP4 Fortran code by Stuart McNicholas 2001
Hacked to CSASArea class structure Liz Potterton 2002

*/

#include <iostream>
#include <string.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <mman_base.h>
#include <mmut_manager.h>
#include <mman_manager.h>
#include <mmut_sasarea.h>
#include "mginterrupt.h"

#define MAXPNT 1500      
#define ZSTEP 0.25
#define ICT 2000
#define MAXNNBR 1500      // specify array sizes. Ugh!
#define PI 3.141592653589793238462643

using namespace std;

//---------------------------------------------------------------------
CSASArea::CSASArea( PCMMUTManager molHndin, int selHndin) : 
                             CMMANBase ( molHndin) {
//-------------------------------------------------------------------
  InitParams();
  // Set CMMANBase to exclude hydrogen, water and alt confs
  int exclude_solvent = 1;
  int exclude_hydrogens = 1;
  int exclude_alternate = 1;
  SetExclusions(exclude_solvent,exclude_hydrogens,exclude_alternate);
  SetSelHandle(selHndin);
}

//---------------------------------------------------------------------
CSASArea::CSASArea( PCMMUTManager molHndin, int selHndin,
                    PCMMUTManager molHndin1, int selHndin1) : 
                       CMMANBase ( molHndin,-1,molHndin1,-1) {
//-------------------------------------------------------------------
  // INitialise with two sets of atoms - for calculation of
  // contact buried SAS
  InitParams();
  // Set CMMANBase to exclude hydrogen, water and alt confs
  int exclude_solvent = 1;
  int exclude_hydrogens = 1;
  int exclude_alternate = 1;
  SetExclusions(exclude_solvent,exclude_hydrogens,exclude_alternate);
  SetSelHandle(0,selHndin,molHndin);
  SetSelHandle(1,selHndin1,molHndin1);
}

//---------------------------------------------------------------------
CSASArea::CSASArea( PCMMUTManager molHndin) :
                             CMMANBase ( molHndin) {
//-------------------------------------------------------------------
  InitParams();
  // Set CMMANBase to exclude hydrogen, water and alt confs
  int exclude_solvent = 1;
  int exclude_hydrogens = 1;
  int exclude_alternate = 1;
  SetExclusions(exclude_solvent,exclude_hydrogens,exclude_alternate);  
  SetSelHandle(-1);
}


//-------------------------------------------------------------------
CSASArea::~CSASArea( ) {
//-------------------------------------------------------------------
}

//-------------------------------------------------------------------
void CSASArea::InitParams() {
//-------------------------------------------------------------------
  //params
  HOHrad = 1.4;
  brick_margin = 5.0;
  point_density = 1.0;
  method =   LEERICHARDS;
  exclude_solvent = 1;
  strcpy(keep_altLoc,"A"); 
  atomUDDHnd = -1;
  resUDDHnd = -1;
}
//-------------------------------------------------------------------
void CSASArea::SetParams(int nv, double *value, int niv, int *ivalue) {
//-------------------------------------------------------------------
  //params
  HOHrad = value[0];
  point_density = value[1];
  brick_margin = value[2];
  exclude_solvent = ivalue[0];
  //cout << " CSASArea::SetParams exclude_solvent " << exclude_solvent << endl;
}



//------------------------------------------------------------------
void  CSASArea::SetUDD(int atomUDD,int resUDD) {
//------------------------------------------------------------------
  atomUDDHnd = atomUDD;
  resUDDHnd = resUDD;
}



//--------------------------------------------------------------------
int CSASArea::SetMethod ( int meth ) {
//--------------------------------------------------------------------
  if ( meth < 0 || meth > 2 ) return 1;
  if ( meth != method ) {
    method = meth;
  }
  return 0;
}

//-------------------------------------------------------------------
int CSASArea::Calculate_Contact ( void )  {
//-------------------------------------------------------------------
  int rv,i,j;
  int tmp_selHnd,local_selHnd, nSelAtoms,nSelAtoms1;
  PPCAtom selAtoms,selAtoms1;
  rvector radwithhoh;

  if ( !selHnds[0] || ! selHnds[1]) return 1;


  molHnds[0]->GetSelIndex(selHnds[0],selAtoms,nSelAtoms);
  molHnds[1]->GetSelIndex(selHnds[1],selAtoms1,nSelAtoms1);
  //cout << " nSelAtoms,nSelAtoms1 " <<  nSelAtoms << " " << nSelAtoms1 << endl;
  if ( nSelAtoms <= 0 || nSelAtoms1 <= 0 ) return 1;

  // Make sure the two selections do not overlap
  // TO BE DONE

  // Initialise UDD
  int nchains = molHnds[0]->GetNumberOfChains(0);
  int nr0,nat0,nat_set;
  PPCResidue resTable;
  PPCAtom atomTable;
  for (int ic=0;ic<nchains;ic++) {
    resTable = NULL;
    molHnds[0]->GetResidueTable(0,ic,resTable,nr0);
    if(resUDDHnd>0) {
      for (i=0;i<nr0;i++) {
        resTable[i]->PutUDData(resUDDHnd,-1.0);
      }
    }
    for (i=0;i<nr0;i++) {
      molHnds[0]->GetAtomTable(0,ic,i,atomTable,nat0);
      for (j=0;j<nat0;j++) {
        atomTable[j]->PutUDData(atomUDDHnd,-1.0);
      }
      nat_set+= nat0;
    }
  }

  // Find the limited number of atoms in selHnd[0] which are close
  // to those in selHnd[1] - these are the only atoms which might
  // have some buried surface area
  tmp_selHnd=molHnds[0]->NewSelection();
  local_selHnd = molHnds[0]->NewSelection();
  // Initially select whole res 
  molHnds[0]->SelectNeighbours(tmp_selHnd,STYPE_RESIDUE,selAtoms1,
			      nSelAtoms1,0.0,10.0,SKEY_NEW);


  molHnds[0]->Select(local_selHnd,STYPE_ATOM,tmp_selHnd,SKEY_NEW);
  molHnds[0]->Select(local_selHnd,STYPE_ATOM,selHnds[0],SKEY_AND);

  // Diagnostic
  //molHnds[0]->GetSelIndex(local_selHnd,tmpAtoms,nTmp);

  // Get the atom radii 
  GetVectorMemory ( radwithhoh, nSelAtoms+nSelAtoms1, 0);
  for(i=0;i<nSelAtoms;i++){
    radwithhoh[i] =  dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius
                               (selAtoms[i]) + HOHrad;
    if (radwithhoh[i] < 1.0 || radwithhoh[i]>4.0)
        cout << "Suspect radius " << radwithhoh[i] << " ";
  }
  j =nSelAtoms; 
  for(i=0;i<nSelAtoms1;i++){
    radwithhoh[j] =  dynamic_cast<PCMMANManager>(molHnds[1])->GetAtomVDWRadius
                               (selAtoms1[i]) + HOHrad;
    if (radwithhoh[j] < 1.0 || radwithhoh[j]>4.0)
        cout << "Suspect radius " << radwithhoh[j] << " ";
    j++;
  }

  // Calculate the SAS of the first model  
  rv = LeeAndRichards(radwithhoh,0,0,local_selHnd);
  if ( rv>0) return rv;
  // Subtract off the SAS of the complex
  rv = LeeAndRichards(radwithhoh,0,1,local_selHnd);
  if ( rv>0) return rv;

  // Calculate the residue SAS area
  if (resUDDHnd>0) ResSASArea(0);


  // Cleanup 
  molHnds[0]->DeleteSelection(tmp_selHnd);
  FreeVectorMemory ( radwithhoh , 0);

  return rv;

}


//-------------------------------------------------------------------
int CSASArea::Calculate (int imodel,bool separate_models)  {
//-------------------------------------------------------------------
   int rv = 0;
  if ( molHnds[0]->GetNumberOfModels() == 1)
    rv = Calculate0(1);
  else if (imodel>0 && imodel <=  molHnds[0]->GetNumberOfModels()) 
    rv = Calculate0(imodel);
  else if (imodel > molHnds[0]->GetNumberOfModels())
    return 1;
  // There is more than 1 model in structure and no specific
  // model specified by imodel
  else if (separate_models ) {
    for (int i = 1; i <= molHnds[0]->GetNumberOfModels(); i++) {
      rv = Calculate0(i);
    }
  }
  else
      rv = Calculate0(0);
  return rv;
}

//-------------------------------------------------------------------
int CSASArea::Calculate0 (int imodel)  {
//-------------------------------------------------------------------

  int RC = 0, i,j,nat;
  PPCAtom selected_atoms;
  rvector radwithhoh;

  // Initiallise mmdb udd
  int nr0,nat0;
  PPCAtom atomTable;
  PPCResidue resTable;
  int im,fm,lm,nchains;
  int nat_set = 0;

  if (imodel<1) {
    fm=1;
    lm =  molHnds[0]->GetNumberOfModels();
  } else {
    fm = imodel;
    lm = imodel;
  }

  for (im=fm;im<=lm;im++) {
    nchains = molHnds[0]->GetNumberOfChains(im);

    for (int ic=0;ic<nchains;ic++) {
      resTable = NULL;
      molHnds[0]->GetResidueTable(im,ic,resTable,nr0);
      if(resUDDHnd>0) {
        for (i=0;i<nr0;i++) {
          resTable[i]->PutUDData(resUDDHnd,-1.0);
        }
      }
      for (i=0;i<nr0;i++) {
        molHnds[0]->GetAtomTable(im,ic,i,atomTable,nat0);
        for (j=0;j<nat0;j++) {
          atomTable[j]->PutUDData(atomUDDHnd,-1.0);
        }
        nat_set+= nat0;
      }
    }
  }
 
  if (imodel>0) {
    GetSelection (selected_atoms, nat,imodel);
  } else {
    molHnds[0]->GetSelIndex(selHnds[0],selected_atoms, nat);
  }
  //cout << "Calculate0 imodel " << imodel << " nat " << nat << endl; 

  GetVectorMemory ( radwithhoh, nat, 0);
  for(i=0;i<nat;i++){
    radwithhoh[i] =  dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius
                               (selected_atoms[i]) + HOHrad;
    //cout << " " << radwithhoh[i];
    if (radwithhoh[i] < 1.0 || radwithhoh[i]>4.0)
        cout << "Suspect radius " << radwithhoh[i] << " ";
  }

  total_area = 0.0;

  // Calculate by the required method
  switch ( method ) {
  case SURFACEPOINTS:
    RC = SurfacePoints(radwithhoh,imodel);
    break;
  case LEERICHARDS:
    RC = LeeAndRichards(radwithhoh,imodel);
    break;
  case WODAKJANIN:
    RC = WodakAndJanin(radwithhoh,imodel);
    break;
  }
  
  FreeVectorMemory ( radwithhoh , 0);

  // If method failed then bail out here
  if ( RC ) return RC;
    
  // Calculate the residue SAS area
  if (resUDDHnd>0) ResSASArea(imodel);
 
  return 0;
}


//-----------------------------------------------------------------------
int CSASArea::LeeAndRichards(realtype *radwithhoh,int imodel,int imode,
                     int local_selHnd ){
//-----------------------------------------------------------------------

  PCBrick brick;
  realtype total_carea;
  int i,j,k,l,nnbr;
  realtype aarc,xr,yr,zr,rr,rrx2,rrsq;
  int x,y,z,ix,iy,iz;
  realtype sumr, dx, dy, dz, dijsq,rj;
  realtype *rnnbr;
  realtype **xyznnbr;
  int nzp;
  realtype zgrid,rsec2r,rsecr,rsec2n,rsecn;
  realtype *arci, *arcf;
  realtype arcsum;
  int karc;
  realtype *dsq,  *d;
  realtype *dxn,  *dyn;
  realtype alpha, beta, ti, tf;
  realtype aarea,carea,aarea0;
  int skip_plane;
  PPCAtom selected_atoms, tmp_sel_atoms[2];
  int natoms_sel, tmp_natoms[2];

  if (imode == 0) {
    // Straightforward SAS for one set of atoms
    if ( imodel > 0 ) {
      GetSelection ( selected_atoms, natoms_sel,imodel);
    } else {
      molHnds[0]->GetSelIndex(selHnds[0],selected_atoms, natoms_sel);
    }
  } else {
    // Second part of calculating buried area - merge two sets of atoms
    molHnds[0]->GetSelIndex(selHnds[0],tmp_sel_atoms[0],tmp_natoms[0]);
    molHnds[1]->GetSelIndex(selHnds[1],tmp_sel_atoms[1],tmp_natoms[1]);
    //cout << "tmp_natoms " << tmp_natoms[0] << " " << tmp_natoms[1] << endl;
    natoms_sel = tmp_natoms[0]+tmp_natoms[1];
    selected_atoms= new PCAtom[natoms_sel];
    for (i=0;i<tmp_natoms[0];i++) selected_atoms[i] = tmp_sel_atoms[0][i];
    j = tmp_natoms[0];
    for (i=0;i<tmp_natoms[1];i++) {
      selected_atoms[j] = tmp_sel_atoms[1][i];
      j++;
    }
  }
  //cout << "LeeAndRichards natoms_sel" << natoms_sel << endl;

  arci = new realtype[ICT];
  arcf = new realtype[ICT];
  dxn = new realtype[MAXNNBR];
  dyn = new realtype[MAXNNBR];
  d = new realtype[MAXNNBR];
  dsq = new realtype[MAXNNBR];
  rnnbr = new realtype[MAXNNBR];
  xyznnbr = new realtype*[MAXNNBR];
  for(i=0;i<MAXNNBR;i++)
    xyznnbr[i] = new realtype[3];

  total_carea = 0.0;

  // cout << "Lee and Richards Calculation\n";

  molHnds[0]->MakeBricks(selected_atoms,natoms_sel,brick_margin,(HOHrad +2.0)*2.0);
  //cout << "natoms_sel " << natoms_sel << endl;
  for(i=0;i<natoms_sel;i++){
    if(!selected_atoms[i]->Ter && 
      (local_selHnd<=0 || selected_atoms[i]->isInSelection(local_selHnd) ) ) {

      if (mginterrupt::GetStatus()>0) {
        cout << "SAS interrupted" << endl;
        return 1;
      }

    aarc = 0;
    //cout << "i " << i << selected_atoms[i]->name << endl;
    xr = selected_atoms[i]->x;
    yr = selected_atoms[i]->y;
    zr = selected_atoms[i]->z;
    rr = radwithhoh[i];
    rrx2 = rr * 2;
    rrsq = rr * rr;

    molHnds[0]->GetBrickCoor(selected_atoms[i],x,y,z);
    
    nnbr = -1;
    for(iz=z-1;iz<z+2;iz++){
      for(iy=y-1;iy<y+2;iy++){
        for(ix=x-1;ix<x+2;ix++){
          brick = NULL;
          brick = molHnds[0]->GetBrick(ix,iy,iz);
          //cout << "brick " << iz << " " << iy << " " << ix << " " << brick;
          //cout << " " << brick->nAtoms << endl;
          if ( brick ) { 
	  for(k=0;k<brick->nAtoms;k++){ // Loop over atoms in bin.
            j = brick->id[k];
            if(j!=i){
              rj = radwithhoh[j];
              sumr = rr + rj;
              dx = selected_atoms[j]->x - xr;
              dy = selected_atoms[j]->y - yr;
              dz = selected_atoms[j]->z - zr;
	      //              if(1>0){
	      if(fabs(dx)<sumr && fabs(dy)<sumr && fabs(dz)<sumr){
		//  dijsq = dx*dx + dy*dy;
		dijsq = dx*dx + dy*dy + dz*dz;
                if(dijsq<(sumr*sumr)){
                  nnbr++;
                  if(nnbr>MAXNNBR){
                    cout << "nnbr > MAXNNBR in solvent accessible area calc" << endl;
                    cout << "Please recompile mmut_sasarea with larger MAXNNBR" << endl;
                    cout << "and report to CCP4" << endl;
		    delete [] xyznnbr;
                    delete [] rnnbr;
                    delete [] dxn;
                    delete [] dyn;
                    delete [] d;
                    delete [] dsq;
                    delete [] arci;
                    delete [] arcf;

                    return 1;
                  }
                  xyznnbr[nnbr][0] = selected_atoms[j]->x;
                  xyznnbr[nnbr][1] = selected_atoms[j]->y;
                  xyznnbr[nnbr][2] = selected_atoms[j]->z;
                  rnnbr[nnbr] = rj;
                  dsq[nnbr] = dx*dx + dy*dy;
                  d[nnbr] = sqrt(dsq[nnbr]);
		  dxn[nnbr] = dx;
		  dyn[nnbr] = dy;
                }
              }
            }
          } }
	}
      }
    }
    //cout << "nnbr" << nnbr << endl;
    if(nnbr>-1){
      nzp = (int)(rrx2/ZSTEP + 0.5);
      zgrid = zr - rr -ZSTEP/2.0;
      for(j=0;j<nzp;j++){
	skip_plane = 0;
	zgrid = zgrid + ZSTEP;
	rsec2r = rrsq - (zgrid-zr)*(zgrid-zr);
	if(rsec2r<0) rsec2r = 0.000001;
	rsecr = sqrt(rsec2r);
	for(k=0;k<ICT;k++){
	  arci[k] = 0.0;
	  arcf[k] = 0.0;
	}
	karc = -1;
	arcsum = 0.0;
	for(l=0;l<=nnbr;l++){
	  rsec2n = rnnbr[l] * rnnbr[l] - (zgrid - xyznnbr[l][2])*(zgrid - xyznnbr[l][2]);
	  if(rsec2n>0){
	    rsecn = sqrt(rsec2n);
	    if(d[l]<rsecr+rsecn){
	      if(d[l]>fabs(rsecr-rsecn)){
		karc++;
		if(karc>ICT){
		  cout << "karc > ICT in solvent accessible area calc\n";
		  cout << "Please recompile mmut_sasarea with larger ICT and\n";
		  cout << "report to CCP4.\n";
		  delete [] xyznnbr;
                  delete [] rnnbr;
                  delete [] dxn;
                  delete [] dyn;
                  delete [] d;
                  delete [] dsq;
                  delete [] arci;
                  delete [] arcf;
		  return 1;
		}else{
		  if((dsq[l]+rsec2r-rsec2n)/(2*d[l]*rsecr)>1){
		    
		    
		    alpha = 0.999;
		  }else{
		    alpha = acos((dsq[l]+rsec2r-rsec2n)/(2*d[l]*rsecr));
		  }
		  beta = atan2(dyn[l],dxn[l]);
		  if(dyn[l]<0) beta += 2*PI;
		  ti = beta - alpha;
		  tf = beta + alpha;
		  if(ti<0) ti += 2*PI;
		  if(tf>(2*PI)) tf -= 2*PI;
		  arci[karc] = ti;
		  if(tf<ti){
		    arcf[karc] = 2*PI;
		    karc++;
		  }
		  arcf[karc] = tf;
		}
	      }else{
		if((rsecr-rsecn)<0){
		  // GOTO 170
		  skip_plane = 1;
		  l = nnbr;
		}
	      } // Endif d[l] > fabs(rsecn-rsecr)
	    } // Endif d[l] < rsecn+rsecr
	  } // Endif rsec2n > 0
	} // End l loop
	if(!skip_plane){
	  if(karc>-1){
	    SortAB(arci,arcf,karc);
	    ArcLap(arci,arcf,karc);
	    arcsum  = arci[0];
	    if(karc>0){
	      for(l=1;l<=karc;l++){
		arcsum += arci[l] - arcf[l-1];
	      }
	    }
	    arcsum += 2*PI - arcf[karc];
	  }else{
	    arcsum = 2*PI;
	  }
	}
	aarc += arcsum; //170
      } // Loop over NZP
      aarea = aarc*rr*ZSTEP;
    }else{
      aarea = 4*PI*rr*rr;
    } // if nnbr > -1
    carea = (rr - HOHrad) * (rr - HOHrad) * aarea/(rr*rr);
    total_area += aarea;
    total_carea += carea;
    if ( imode == 0 ) {
      selected_atoms[i]->PutUDData(atomUDDHnd,aarea);
    } else {
      selected_atoms[i]->GetUDData(atomUDDHnd,aarea0);
      selected_atoms[i]->PutUDData(atomUDDHnd,aarea0-aarea);
    }
  
  } else {
    // Atom excluded from local_selHnd but should be given a value
    selected_atoms[i]->PutUDData(atomUDDHnd,0.0);
  } } // end of atom loop


  //cout << "Total solvent accessible area: " << total_area << endl;
 
  delete [] xyznnbr;

  delete [] rnnbr;
  delete [] dxn;
  delete [] dyn;
  delete [] d;
  delete [] dsq;
  delete [] arci;
  delete [] arcf;
  if (imode ==1 && selected_atoms) delete[] selected_atoms;

  return 0;
}

//-------------------------------------------------------------------------
int CSASArea::WodakAndJanin(realtype *radwithhoh,int imodel){
//-------------------------------------------------------------------------

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelection ( selected_atoms, natoms_sel,imodel);
  // This doesn't seem to want to work
  realtype b,bpr,ri,rj,dij;
  int i,j;
 
  realtype A;
  realtype B;
  realtype S;

  for(i=0;i<natoms_sel;i++){
    if(!selected_atoms[i]->Ter) {
    S = 4 * PI * radwithhoh[i]*radwithhoh[i];
    A = S;
    //    cout << "Original S: " << S << endl;
    B = 0.0;
    for(j=0;j<natoms_sel;j++){
    if(!selected_atoms[j]->Ter) {
      if(i!=j){
	ri =  dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius(selected_atoms[i]);
	rj =  dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomVDWRadius(selected_atoms[j]);

      dij = molHnds[0]->BondLength(selected_atoms[i],selected_atoms[j]);
      if(dij<(ri + rj + 2*HOHrad)){
	//	cout <<"dij, i,j: " << dij << " " << i << ", " << j << endl;
	//	cout << "i: " << selected_atoms[i]->x << " " << selected_atoms [i]->y << " "<< selected_atoms[i]->z << endl;
	//	cout << "j: " << selected_atoms[j]->x << " " << selected_atoms[j]->y << " "<< selected_atoms[j]->z << endl;
	// Maximum area cut out of i by j.
	b = PI * (ri + HOHrad) * (ri + rj + 2*HOHrad - dij)*(1 + (rj - ri)/dij);
	// Minimum area assuming hard spheres s = 2*HOHRAD
	bpr = PI * (ri + HOHrad) * (ri + rj - dij)*(1 + (rj - 2*HOHrad - ri)/dij);
      }else{
	b = 0.0;
	bpr = 0.0;
      }

      //      cout << "b, bpr: " << b << " " << bpr << endl;
      //      cout << "b - bpr: " <<  b-bpr << endl;
      //      cout << "(b - bpr)/S: " <<  (b-bpr)/S << endl;

      if(bpr<0) bpr = 0;
      A *= (1 - (b - bpr)/S); 
      B += bpr;
      //      cout << "Current A: " << A << ", Current B: " << B << endl;
      }
    }
    }
    //    cout << "A,B: " << A << B << endl;
    if(A<B) {
      //      cout << "Area on atom " << i << ": 0\n";
    }else{
            cout << "Area on atom " << i << ": " << A - B << endl;
      total_area = total_area + A - B;
      selected_atoms[i]->PutUDData(atomUDDHnd,A-B);
    }    
    }else{
      selected_atoms[i]->PutUDData(atomUDDHnd,0.0);
    }
  }

  return 0;

}

//-------------------------------------------------------------------------
int CSASArea::SurfacePoints(realtype *radwithhoh,int imodel) {
//-------------------------------------------------------------------------

  PPCAtom selected_atoms;
  int natoms_sel;
  GetSelection ( selected_atoms, natoms_sel,imodel);
  PCBrick brick;

  int i,j,k,x,y,z;
  int nAtomTypes,atom_type;
  realtype ri,rj,rk,sarea,points,equat,vert,phi,theta;
  realtype hor,trunc,rxy,rz;
  int nvert,ipt,nhor;
  realtype ***xyzsph;
  int nnbr, ix, iy, iz;
  realtype sumr, dx, dy, dz, dijsq;
  realtype **xyznnbr;
  realtype *rnnbr;
  int *nptsph;
  int issurface;
  int *npoint;

 
  npoint = new int[natoms_sel];

  nAtomTypes = dynamic_cast<PCMMANManager>(molHnds[0])->GetMGSBase()->GetNofLibAtoms();
  //cout << "nAtomTypes" << nAtomTypes << endl;
  
  xyzsph = new realtype**[nAtomTypes];
  for(i=0;i<nAtomTypes;i++){
    xyzsph[i] = new realtype*[MAXPNT];
    for(j=0;j<MAXPNT;j++)
      xyzsph[i][j] = new realtype[3];
  }

  rnnbr = new realtype[MAXNNBR];
  xyznnbr = new realtype*[MAXNNBR];
  for(i=0;i<MAXNNBR;i++)
    xyznnbr[i] = new realtype[3];

  nptsph = new int[nAtomTypes];


  //  cout << "Calculating Solvent accessible area by surface points\n";

  // Brick size from areaimol and margin 5.

  molHnds[0]->MakeBricks(selected_atoms,natoms_sel,brick_margin,(HOHrad+2.0)*2.0);

  for(i=0;i<nAtomTypes;i++){
      ri =  dynamic_cast<PCMMANManager>(molHnds[0])->GetMGSBase()->libAtom[i].vdwRadius;
      //cout << "Atom type rad " << ri << endl;
      sarea = 4*PI*ri*ri;
      points = sarea*point_density;
      equat = sqrt(points*PI);
      vert = equat/2.0;
      nvert = (int)vert;

      ipt = -1;
      trunc = 0.0;

      for(j=0;j<nvert;j++){
        phi = (realtype)j*PI/vert;
        hor = sin(phi)*equat - trunc;
        nhor = (int)hor;
        if(j==0) nhor =1;
        trunc = (realtype)nhor - hor;
        rxy = sin(phi)*ri;
        rz = cos(phi)*ri;

        for(k=0;k<nhor-1;k++){
          theta = 2.0*PI*(realtype)k/(realtype)nhor;
          ipt++;
          if(ipt>MAXPNT){
            cout << "npt > MAXPNT in solvent accessible area calc\n";
            cout << "Please recompile with larger MAXPNT and\n";
            cout << "report this as ugly code.\n";
            exit(1);
          }
          xyzsph[i][ipt][0] = cos(theta)*rxy;
          xyzsph[i][ipt][1] = sin(theta)*rxy;
          xyzsph[i][ipt][2] = rz;
        }
      }
      nptsph[i] = ipt;
  }

  for(i=0;i<nAtomTypes;i++)
    if(nptsph[i]>0)
      //      printf("Number of points on atom type %d[%s]: %d\n",i,AtomTypes[i],nptsph[i]);

  for(i=0;i<natoms_sel;i++){
    npoint[i]=0;
    atom_type = dynamic_cast<PCMMANManager>(molHnds[0])->GetAtomEnergyType(selected_atoms[i]);
    if(!selected_atoms[i]->Ter) {
    // Should be an if water, etc condition here.
    ri = radwithhoh[i];
    molHnds[0]->GetBrickCoor(selected_atoms[i],x,y,z);
    nnbr = -1;
    for(iz=z-1;iz<z+2;iz++){
      for(iy=y-1;iy<y+2;iy++){
        for(ix=x-1;ix<x+2;ix++){
          brick = molHnds[0]->GetBrick(ix,iy,iz);
          for(k=0;k<brick->nAtoms;k++){ // Loop over atoms in bin.
            j = brick->id[k];
            if(j!=i){
              rj = radwithhoh[j];
              sumr = ri + rj;
              dx = selected_atoms[i]->x - selected_atoms[j]->x;
              dy = selected_atoms[i]->y - selected_atoms[j]->y;
              dz = selected_atoms[i]->z - selected_atoms[j]->z;
	      if(fabs(dx)<sumr && fabs(dy)<sumr && fabs(dz)<sumr){
		dijsq = dx*dx + dy*dy +dz*dz;
                if(dijsq<(sumr*sumr)){
                  nnbr++;
                  if(nnbr>MAXNNBR){
                    cout << "nnbr > MAXNNBR in solvent accessible area calc\n";
                    cout << "Please recompile with larger MAXNNBR and\n";
                    cout << "report this as ugly code.\n";
                    exit(1);
                  }
                  xyznnbr[nnbr][0] = selected_atoms[j]->x;
                  xyznnbr[nnbr][1] = selected_atoms[j]->y;
                  xyznnbr[nnbr][2] = selected_atoms[j]->z;
                  rnnbr[nnbr] = rj;
                }
              }
            }
          }
        }
      }
    }
    for(j=0;j<nptsph[atom_type];j++){ // Loop over points on sphere.
      issurface = 1;
      for(k=0;k<=nnbr;k++){
        rk = rnnbr[k];
        dx=0.0;
        dy=0.0;
        dz=0.0;
        dx = selected_atoms[i]->x + xyzsph[atom_type][j][0] - xyznnbr[k][0];
        dy = selected_atoms[i]->y + xyzsph[atom_type][j][1] - xyznnbr[k][1];
        dz = selected_atoms[i]->z + xyzsph[atom_type][j][2] - xyznnbr[k][2];
	if(fabs(dx)<rk && fabs(dy)<rk && fabs(dz)<rk){
          dijsq = dx*dx + dy*dy +dz*dz;
          if(dijsq<rk*rk) {
            k = nnbr;
            issurface = 0;
          }
        }
      }
      if(issurface) npoint[i]++;
    }
    }
  }

  delete [] xyznnbr;

  for(i=0;i<nAtomTypes;i++){
    delete [] xyzsph[i];
  }
  delete [] xyzsph;

  delete [] rnnbr;  
  delete [] nptsph;
  delete [] xyzsph;
  delete [] xyznnbr;

  for(i=0;i<natoms_sel;i++){
    total_area += (realtype)npoint[i];
    if(!selected_atoms[i]->Ter) {
      selected_atoms[i]->PutUDData(atomUDDHnd,(realtype)npoint[i]/point_density);
    }else{
      selected_atoms[i]->PutUDData(atomUDDHnd,0.0);
    }
  }
  total_area = total_area/point_density;
  delete [] npoint;

  return 0;
  
}

//--------------------------------------------------------------------
int CSASArea::ResSASArea (int imodel) {
//--------------------------------------------------------------------
  // Sum over all atoms in a residue for residue SAS area
  realtype sas,tot;
  int nr0,nat0, ntot, RC;
  PPCAtom atomTable;
  PPCResidue resTable;
  int first_model=1;
  int last_model=1;
  if (imodel == 0 && molHnds[0]->GetNumberOfModels() > 1) {
    first_model = 1;
    last_model =  molHnds[0]->GetNumberOfModels();
  }
  else if (imodel>0) {
    first_model = imodel;
    last_model = imodel;
  }
 
  for (int im= first_model; im <=last_model;im++) {
    //cout << " ResSASArea imodel " << im << endl;
    int nchains = molHnds[0]->GetNumberOfChains(im);
    for (int ic=0;ic<nchains;ic++) {
      molHnds[0]->GetResidueTable(im,ic,resTable,nr0);
      for (int i=0;i<nr0;i++) {
        molHnds[0]->GetAtomTable(im,ic,i,atomTable,nat0);
        tot = 0.0;
        ntot = 0;
        for (int j=0;j<nat0;j++) {
          RC = atomTable[j]->GetUDData(atomUDDHnd,sas);
          if (RC == 0 && sas >= 0.0 ) {
             tot = tot + sas;
             ntot = ntot + 1;
          }
        }
        if ( ntot > 0 )
          resTable[i]->PutUDData(resUDDHnd,tot);
        else
          resTable[i]->PutUDData(resUDDHnd,-1.0);
        //cout << "ResSASArea " << i << " " << tot << endl;
      }
    }
  } 
  return 0;
}

//-----------------------------------------------------------------------
std::string CSASArea::Print(int imodel, string selection1, int set_selHnd1, string selection2, int set_selHnd2 ) {
//-----------------------------------------------------------------------
  PPCAtom selected_atoms;
  PPCResidue selected_res;
  int RC,i,natoms_sel,nres_sel;
  realtype sas,total, set_total1,set_total2 ;
  std::ostringstream output;

  output.setf(ios::fixed);
  output.setf(ios::showpoint);
  output.setf(ios::left,ios::adjustfield);
  output.precision(2);
  int first_model=1;
  int last_model=1;

 
  if (imodel == 0 && molHnds[0]->GetNumberOfModels() > 1) {
    first_model = 1;
    last_model =  molHnds[0]->GetNumberOfModels();
  }
  else if (imodel>0) {
    first_model = imodel;
    last_model = imodel;
  }

  for (int im=first_model;im<=last_model;im++) {

    if (im != 1) output << "For model number " << im << endl;


    if ( resUDDHnd > 0 ) {  
      GetSelection ( selected_res, nres_sel,im);
      output << "\nSolvent Accessible area by residue\n" << "----------------------------------\n";
    
      for (i=1;i<nres_sel;i++){
        RC = selected_res[i]->GetUDData(resUDDHnd,sas);
        if ( RC == 0 && sas > 0.001) output << std::setw(15) <<molHnds[0]->AtomLabel_residue1(selected_res[i]) << sas << endl;
      }
    }

    if ( atomUDDHnd > 0 ) {
      total = 0.0;
      set_total1 = 0.0;
      set_total2 = 0.0;
      GetSelection ( selected_atoms, natoms_sel,im);
      output << "\nSolvent Accessible area by atom\n"
             << "-------------------------------" << endl;
       
      for(i=1;i<natoms_sel;i++) {
        RC = selected_atoms[i]->GetUDData(atomUDDHnd,sas);
        if ( RC == 0 && sas > 0.001) { 
          output << std::setw(15) <<
            molHnds[0]->AtomLabel_atom(selected_atoms[i]) << sas << endl;
          total+= sas;
          if (set_selHnd1>0 && selected_atoms[i]->isInSelection(set_selHnd1)) set_total1 += sas;
          if (set_selHnd2>0 && selected_atoms[i]->isInSelection(set_selHnd2)) set_total2 += sas;
        }
      }
      output << "Total solvent accessible area: " << total << "Angstoem*2" << endl << endl;
      if (set_selHnd1>0)  output << "Solvent accessible area of " << selection1 << ": " << set_total1 << "Angstoem*2" << endl << endl;
      if (set_selHnd2>0)  output << "Solvent accessible area of " << selection2 << ": " << set_total2 << "Angstoem*2" << endl << endl;
    }
  }
  return output.str();
}


//----------------------------------------------------------------
void CSASArea::SortAB(realtype *arci, realtype *arcf, int &karc)  {
//----------------------------------------------------------------
  int i,j,it;
  realtype aa,bb,temp;

  for(j=0;j<=karc;j++){
    aa = arci[j];
    bb = arcf[j];
    temp = 10.0;
    it = j;
    for(i=j;i<=karc;i++){
      if(arci[i]<=temp){
	temp = arci[i];
	it = i;
      }
    }
    arci[j] = arci[it];
    arcf[j] = arcf[it];
    arci[it] = aa;
    arcf[it] = bb;
  }

}


//--------------------------------------------------------------
void CSASArea::ArcLap(realtype *arci, realtype *arcf, int &karc)  {
//--------------------------------------------------------------
  realtype t;
  int k,l,ks;//,m;

  t = arcf[0];
  ks = 1;

  // This might work .... See Fortran

  for(k=ks;k<=karc;){
    if(t>arci[k]){
      if(arcf[k]>arcf[k-1]) arcf[k-1] = arcf[k];
      t = arcf[k-1];
      for(l=k+1;l<=karc;l++){
	arci[l-1] = arci[l];
	arcf[l-1] = arcf[l];
      }
      karc--;
    }else{
      t = arcf[k];
      k++;
    }
  }

}
