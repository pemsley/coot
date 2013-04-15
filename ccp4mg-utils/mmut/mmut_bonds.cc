/*
     mmut/mmut_bonds.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York
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


#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <mmdb/mmdb_manager.h>
#include <mmdb/mmdb_sbase.h>
#include <mmdb/mmdb_graph.h>
#include <mmdb/mattype_.h>
#include "mmut_bonds.h"
#include "mmut_sbase.h"
#include "mmut_manager.h"
#include "mman_manager.h"
#include <mmdb/mmdb_tables.h>

using namespace std;

//-----------------------------------------------------------------------
CMolBondParams::CMolBondParams( CMGSBase *sbasein ) {
//-----------------------------------------------------------------------
   sbase = sbasein;
   //printf ("Sbase maxAtomInRes %i\n",sbase->maxAtomInRes);
   interResCut = 2.4; 
   intraResCut = 2.4;
   maxBondRad = 2.4;
   // Mx bond length is ( at1->maxBondRad + at2->maxBondRad)
   // Make some allowance upwards
   maxBondRadFactor = 1.2; 
}

//-----------------------------------------------------------------------
CMolBondParams::~CMolBondParams() {
//-----------------------------------------------------------------------
 
}


//-----------------------------------------------------------------------
CMolBonds::CMolBonds(const PCMMUTManager molHndin, CMolBondParams *paramsin) :
  CMMANBase ( molHndin ) {
//-----------------------------------------------------------------------

  params = paramsin;
  
  sbaseCompoundID = NULL;
  sbaseAtomIndex = NULL;

  nAtoms = 0;
  nRes = 0;

}
//----------------------------------------------------------------------
CMolBonds::~CMolBonds(){
//----------------------------------------------------------------------
  //cout << "CMolBonds destructor" << endl;
}
 

//----------------------------------------------------------------------
std::string CMolBonds::FindBonds ( int udd_sbaseCompoundID,
        int udd_sbaseAtomOrdinal, int udd_atomEnergyType ) {
//----------------------------------------------------------------------

  std::ostringstream output;
  char AtomID[30];

  PPCAtom selAtom = NULL;
  PPCAtom selAtom0 = NULL;
  PPCResidue selRes = NULL;
  PCChain pCh;
  PCModel pMdl;
  PCSBStructure pSbRes=NULL;
  ivector nMatchAtom = NULL;
  imatrix matchAtom = NULL;
  int nAtominRes,nAtoms,ir,i,j,na1,na2,ia1,ia2,ib,selH;
  int nModels,nm,nAtominModel; 

  int nr;
  PSContact contacts = NULL;
  int ic,nContacts;
  PCAtom pa1,pa2;
  PCResidue pr1;

  int nAlt;
  int nAltTot = 0 ;
  pstr sbaseCompoundID;
  bool doContacts;  
  bool modelBondsSame = true;
  int maxModels = 1;
  int firstModel = 1;
  int lastModel = 1;
  int maxMatch = 5;

  int restype1,restype2;

  // Pointer to atoms in residue for all models
  // Big bovver if >100 models
  /*
  PPCAtom modelSelAtom[100];
  int nSelAtom[100];
  for (nm=1;nm<=100;nm++) {  modelSelAtom[nm]=NULL; }
  */

  // For NMR structure with multiple models we need to 
  // find the bonds in one model but populate the bonds data
  // structure for all models - but first make sure all models
  // have the same composition 
  nModels = molHnds[0]->GetNumberOfModels();

  // Avoid big bovver. May use a lot of memory.
  PPCAtom *modelSelAtom = new PPCAtom[nModels+1];
  int* nSelAtom = new int[nModels+1];
  for (nm=0;nm<=nModels;nm++) {  modelSelAtom[nm]=NULL; }

  if (nModels>1) {
    lastModel = nModels;
    firstModel = 0;
    PCModel pMdl = NULL;
    while ( firstModel < lastModel && pMdl == NULL ) {
      firstModel++;
      pMdl = molHnds[0]->GetModel(firstModel);
    }
    //cout << "firstModel " << firstModel << endl;
  }
  nAtominModel = molHnds[0]->GetModel(firstModel)->GetNumberOfAtoms(0);
  if (lastModel>firstModel) {
    for (nm=firstModel+1;nm<=lastModel;nm++) {
      pMdl = molHnds[0]->GetModel(nm);
      if (pMdl != NULL && pMdl->GetNumberOfAtoms(0) != nAtominModel) {
        modelBondsSame = false;
        output << "Models do not have the same number of atoms" << endl;
      }
    }
  }
  if ( !modelBondsSame) {
    maxModels = nModels;
  } else {
    maxModels = firstModel;
  }
  //cout << endl << "modelBondsSame " << modelBondsSame << " " << maxModels << endl;

  // Get memory for   imatrix matchAtom; 
  // the pointers to atom in the Sbase structure
  GetMatrixMemory(matchAtom,maxMatch ,params->sbase->maxAtomInRes,0,0);
  GetVectorMemory(nMatchAtom,params->sbase->maxAtomInRes,0);

  // Loop over the unique models (normally just one model)
  for (int iMod=firstModel;iMod<=maxModels;iMod++) {
  if  (molHnds[0]->GetModel(iMod)!= NULL) {

    //cout << "iMod " << iMod << endl;
    for (int ich = 0; ich< molHnds[0]->GetModel(iMod)->GetNumberOfChains();ich++) {
    //Get the residue selections
    pCh = NULL;
    pCh = molHnds[0]->GetModel(iMod)->GetChain(ich);
    selRes = NULL;
    molHnds[0]->GetModel(iMod)->GetResidueTable (pCh->GetChainID(), selRes, nr);


    // To find the INTRA-residue bonds
    // Loop over all residues finding the required Sbase residue
    nAtominRes=0;
    for (ir = 0; ir < nr; ir++ ) {
      molHnds[0]->GetAtomTable1(iMod,pCh->GetChainID(),ir,selAtom,nAtominRes);
      //cout << endl << "ich,ir,nAtominRes " << ich << " " << ir << " " << selAtom[0]->serNum << " " << nAtominRes;
      if (modelBondsSame && nModels>1) {
        for (nm=firstModel+1;nm<=lastModel;nm++) {
          molHnds[0]->GetAtomTable1(nm,pCh->GetChainID(),ir,modelSelAtom[nm],nSelAtom[nm]);
        }
      }
      doContacts = false;
      sbaseCompoundID = 0;
      selRes[ir]->GetUDData(udd_sbaseCompoundID,sbaseCompoundID);
      if ( strlen(sbaseCompoundID) >= 1 ) {
        pSbRes = params->sbase->GetStructure(sbaseCompoundID, dynamic_cast<PCMMANManager>(molHnds[0])->monlib, dynamic_cast<PCMMANManager>(molHnds[0])->GetUnremediated());
        if(sbaseCompoundID) delete [] sbaseCompoundID;
        //cout << "FindBonds " << ir << " " << sbaseCompoundID << " " << pSbRes <<  endl;
        for ( j = 0; j < pSbRes->nAtoms; j++ ) nMatchAtom[j] = 0;
          for (ia1=0;ia1 < nAtominRes;ia1++) {

            selAtom[ia1]->GetUDData( udd_sbaseAtomOrdinal,j );
            //cout << selAtom[ia1]->name << " " << selAtom[ia1]->residue->seqNum << " " << selAtom[ia1]->residue->chain->GetChainID() << " " << j << endl;
            if ( j >= 0 ) {
              if ( j <  pSbRes->nAtoms && nMatchAtom[j] < maxMatch) {
                matchAtom[nMatchAtom[j]][j]=ia1;
	        nMatchAtom[j]++;
              } else {
                //cout << "Error in FindBonds for " << /*selAtom[ia1]->GetAtomID(atomid) << */" j=" << j << " nMatchAtom[j]=" << nMatchAtom[j] << endl; 
              }
            }
            else 
              doContacts = true; 
          }
          nAlt = 0;
          for ( j = 0; j < pSbRes->nAtoms; j++ ) { 
            if (nMatchAtom[j]>1) nAlt++;
          }
          nAltTot = nAltTot + nAlt;
    
          // The residue has alternate locations so there may be more
          // than one atom of a given atom name and more than one 
          // equivalent bond  
          if ( nAlt >= 1 ) {
            for ( ib = 0; ib < pSbRes->nBonds; ib++ ) {
              realtype dist = 0;
              na1 = nMatchAtom[pSbRes->Bond[ib]->atom1-1];
              na2 = nMatchAtom[pSbRes->Bond[ib]->atom2-1];
              if (na1>0 && na2>0) {
                for ( i = 0; i < na1; i++ ) {
                  ia1 = matchAtom[i][pSbRes->Bond[ib]->atom1-1];
                  for ( j = 0; j< na2; j++ ) {
                    ia2 = matchAtom[j][pSbRes->Bond[ib]->atom2-1];
                    realtype x1 = selAtom[ia1]->x;
                    realtype y1 = selAtom[ia1]->y;
                    realtype z1 = selAtom[ia1]->z;
                    realtype x2 = selAtom[ia2]->x;
                    realtype y2 = selAtom[ia2]->y;
                    realtype z2 = selAtom[ia2]->z;
                    dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
                    if (molHnds[0]->doAltLocMatch (selAtom[ia1],selAtom[ia2] )&&params->sbase->CheckCovalentDistance(selAtom[ia1]->element,selAtom[ia2]->element,dist)!=0){ 
                      AddConnection(ia1, ia2, selAtom);
                      if (modelBondsSame && nModels > 1 ) {
                        for (nm = firstModel+1; nm <= lastModel; nm++) {
                          if (molHnds[0]->GetModel(nm) != NULL) {
                            AddConnection(ia1, ia2, modelSelAtom[nm]);
                          }
                        }
                      }
                    }
                  }  
                } 
              }
            }
          }
           else {
            // No alternate locations - should only be one bond of each type
            for ( ib = 0; ib < pSbRes->nBonds; ib++ ) {
              realtype dist = 0;
              if ( nMatchAtom[pSbRes->Bond[ib]->atom1-1] > 0 &&
		     nMatchAtom[pSbRes->Bond[ib]->atom2-1] > 0){
                 realtype x1 = selAtom[matchAtom[0][pSbRes->Bond[ib]->atom1-1]]->x;
                 realtype y1 = selAtom[matchAtom[0][pSbRes->Bond[ib]->atom1-1]]->y;
                 realtype z1 = selAtom[matchAtom[0][pSbRes->Bond[ib]->atom1-1]]->z;
                 realtype x2 = selAtom[matchAtom[0][pSbRes->Bond[ib]->atom2-1]]->x;
                 realtype y2 = selAtom[matchAtom[0][pSbRes->Bond[ib]->atom2-1]]->y;
                 realtype z2 = selAtom[matchAtom[0][pSbRes->Bond[ib]->atom2-1]]->z;
                 dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
              }
              if ( nMatchAtom[pSbRes->Bond[ib]->atom1-1] > 0 &&
		     nMatchAtom[pSbRes->Bond[ib]->atom2-1] > 0 &&
                  params->sbase->CheckCovalentDistance(selAtom[matchAtom[0][pSbRes->Bond[ib]->atom1-1]]->element,selAtom[matchAtom[0][pSbRes->Bond[ib]->atom2-1]]->element,dist)!=0
                ) {
                AddConnection ( matchAtom[0][pSbRes->Bond[ib]->atom1-1],
			  matchAtom[0][pSbRes->Bond[ib]->atom2-1],selAtom) ;
                if (modelBondsSame && nModels > 1 ) {
                  for (nm = firstModel+1; nm <= lastModel; nm++) {
		  if (molHnds[0]->GetModel(nm) != NULL) {
                    if (nSelAtom[nm] == nAtominRes) {
                      //cout << "  AddConnection " << nm << " " << selRes[ir]->seqNum << " " <<  pSbRes->Bond[ib]->atom1-1 << " " <<  pSbRes->Bond[ib]->atom2-1 << endl;   cout.flush();       
                      AddConnection(matchAtom[0][pSbRes->Bond[ib]->atom1-1],
		      matchAtom[0][pSbRes->Bond[ib]->atom2-1],modelSelAtom[nm]);
                    } else {
                      pr1 = molHnds[0]->GetResidue(nm,pCh->GetChainID(),ir);
                      if (pr1 != NULL) {
                        pr1->GetResidueID(AtomID);
                        output << "Residue " << AtomID << " in model " << nm << " does not match equivalent in the first model" << endl;
                        // Resideu in model does not match same residue in
                        // first model - just do sledgehammer search for bonds
                        IntraResContacts (pr1 , 0 );
		      } else {
                        output << "Model " <<  nm << " does not match the first model" << endl;
                      }
                    }
                  } } 
                }
              }
	    }
          }
        }

        else
          doContacts = true;

      // There was no match in the Structure database so find
      // bonds on basis of close contacts
      if ( doContacts )  {
         if (modelBondsSame && nModels > 1 ) {
           IntraResContacts ( selRes[ir], 1, modelSelAtom, nSelAtom, firstModel, lastModel );
         } else {
           IntraResContacts ( selRes[ir], 1);
         }
      }
    }
  } }

  //INTER-residue bonds

  selH = molHnds[0]->NewSelection();
  molHnds[0]->SelectAtoms(selH,iMod,"*", ANY_RES,"*", ANY_RES,"*","*","*","*","*",SKEY_NEW);
  molHnds[0]->SelectAtoms(selH,iMod,"*", ANY_RES,"*", ANY_RES,"*","*","*","H","*",SKEY_CLR);
  selAtom0 = NULL;
  molHnds[0]->GetSelIndex ( selH,selAtom0,nAtoms );

  // Find inter-res close contacts
  if (contacts) delete [] contacts;
  contacts = NULL;
  nContacts = 0;
  molHnds[0]->SeekContacts(selAtom0,nAtoms,selAtom0,nAtoms,
		 0.0,params->interResCut,1,contacts,nContacts,0,NULL,0);
  //cout << "FindBonds INTER-residue bonds nContacts " << nAtoms << " " << nContacts << endl;

  // Check if each contact is between possible altLoc matches
  // and that this corresponds to recognised chemical link
  // if ( contacts[ic].id1 < contacts[ic].id2 ) {
  if ( contacts && nContacts > 0 ) { 
    for ( ic = 0; ic < nContacts; ic++) {
      if ( contacts[ic].id1 < contacts[ic].id2 ) {
        pa1 = selAtom0[contacts[ic].id1];
        pa2 = selAtom0[contacts[ic].id2];
        if (  nAltTot <= 0 || molHnds[0]->doAltLocMatch( pa1, pa2 )  ) {
          //cout << "testing contact " << pa1->residue->seqNum << " " << pa2->residue->seqNum << endl;
          restype1 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa1->residue);
          restype2 = dynamic_cast<PCMMANManager>( molHnds[0])->GetRestypeCode (pa2->residue);
          //cout << "testing interres " <<  pa1->residue->seqNum << pa1->residue->name << " " << restype1 << " "  << pa2->residue->seqNum <<  pa2->residue->name  << " " << restype2 << endl;
          if (restype1 ==  RESTYPE_PEPTIDE && restype2 == RESTYPE_PEPTIDE ) {
            if (isInterResBond(pa1,pa2)) 
               AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
          } else if ( (restype1==RESTYPE_NUCL || restype1==RESTYPE_RNA || restype1==RESTYPE_DNA) && 
                 (restype2 == RESTYPE_NUCL || restype2 == RESTYPE_RNA || restype2 == RESTYPE_DNA) ) {
 
            if (isInterResBond(pa1,pa2)) 
               AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
          } else if ( (restype1==RESTYPE_SACH || restype1==RESTYPE_DSACH || restype1==RESTYPE_LSACH) && 
                 (restype2 == RESTYPE_SACH || restype2 == RESTYPE_DSACH || restype2 == RESTYPE_LSACH) ) {
 
            if (isInterResBond(pa1,pa2)) 
               AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 

	  } else if ((restype1<RESTYPE_SOLVENT || restype1>RESTYPE_NONPOLY) &&
               (restype2<RESTYPE_SOLVENT || restype2>RESTYPE_NONPOLY ) ) {
            if (ltBondDistance(pa1,pa2,contacts[ic].dist) )  {
	      //cout << "non-peptide" << pa1->name << " " << pa2->name << endl;
              if((isMetal(pa1->name)&&!isMetal(pa2->name))||(isMetal(pa2->name)&&!isMetal(pa1->name))){
	        //cout << "metal coordination, ignoring" << endl;
              } else {
                AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0);
              }
            }
          } else {
            if (isInterResBond(pa1,pa2)) {
               AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
            } else {
               double dist = contacts[ic].dist;
               if(params->sbase->CheckCovalentDistance(pa1->element,pa2->element,dist)==1){
                 //std::cout << "Adding possible covalent contact\n";
                 //std::cout << pa1->element << " " <<  pa2->element << " " << dist << "\n";
                 AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0); 
               }
            }
          }
        }
      }
    }
    if (modelBondsSame && nModels > 1 ) {
      //int nb;
      for (nm = firstModel+1; nm <= lastModel; nm++) {
      if (molHnds[0]->GetModel(nm) != NULL) {
        //nb = 0;
        molHnds[0]->SelectAtoms(selH,nm,"*", ANY_RES,"*", ANY_RES,"*","*","*","*","*",SKEY_NEW);
        molHnds[0]->SelectAtoms(selH,nm,"*", ANY_RES,"*", ANY_RES,"*","*","*","H","*",SKEY_CLR);
        selAtom0 = NULL;
        molHnds[0]->GetSelIndex ( selH,selAtom0,nAtoms );
        //cout << "selAtom0 " << nm << " nAtoms " << nAtoms << " nContacts " << nContacts << endl;
        for ( ic = 0; ic < nContacts; ic++) {
          if ( contacts[ic].id1 < contacts[ic].id2 && contacts[ic].id2 < nAtoms && contacts[ic].id2 > 0 ) {
            pa1 = selAtom0[contacts[ic].id1];
            pa2 = selAtom0[contacts[ic].id2];
            if ( nAltTot <= 0 || molHnds[0]->doAltLocMatch( pa1, pa2 )) { 
              restype1 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa1->residue);
              restype2 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa2->residue);
	      if (restype1 ==  RESTYPE_PEPTIDE && 
                      restype2 == RESTYPE_PEPTIDE) { 
                if (isInterResBond(pa1,pa2)) {
                  //nb++;
                  AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0);
                }
              } else if (ltBondDistance(pa1,pa2,contacts[ic].dist)) {
		//nb++; 
                  AddConnection (contacts[ic].id1,contacts[ic].id2,selAtom0);
	      }
            }
	  }
        }
        //cout << "nm " << nm << " nb " << nb << endl; 
      } }
    }                     
  }

  }  // End of loop over models

  molHnds[0]->DeleteSelection(selH);
  if (nMatchAtom) FreeVectorMemory(nMatchAtom, 0);
  nMatchAtom = NULL;
  if (matchAtom) FreeMatrixMemory(matchAtom,maxMatch,0,0);
  matchAtom = NULL;
  if (contacts) delete [] contacts;
  if (selAtom) delete [] selAtom;
  /*
  for (nm=2;nm<=nModels;nm++) {
    if (modelSelAtom[nm]) delete [] modelSelAtom[nm];
  }
  */
  //cout << "done FindBonds" << endl;
  if(modelSelAtom) delete [] modelSelAtom;
  if(nSelAtom) delete [] nSelAtom;


  return output.str();
}


//-----------------------------------------------------------
void CMolBonds::AddConnection ( PCAtom pa1 , PCAtom pa2) {
//-----------------------------------------------------------
  
  //cout << "AddConnection " <<  GetMolHnd(0)->AtomLabel_atom(pa1) << " " 
  //	<< GetMolHnd(0)->AtomLabel_atom(pa2) << endl;

  if ( fabs(pa1->x - pa2->x) < 0.0001 &&
       fabs(pa1->y - pa2->y) < 0.0001 &&
       fabs(pa1->z - pa2->z)  < 0.0001 ) {
    cout << "Bonded atoms superposed " <<
        GetMolHnd(0)->AtomLabel_atom(pa1) << " " <<
        GetMolHnd(0)->AtomLabel_atom(pa2) << endl;
    return;
  }
      
   if ( pa1->GetNBonds() == 0 )
     pa1->AddBond(pa2,1,3);
   else 
     pa1->AddBond(pa2,1,1);
   if ( pa2->GetNBonds() == 0 )
     pa2->AddBond(pa1,1,3);
   else 
     pa2->AddBond(pa1,1,1);

}


//-----------------------------------------------------------
int CMolBonds::DeleteConnection ( PCAtom pa1 , PCAtom pa2) {
//-----------------------------------------------------------
// Delete the connection between two atom - need to zap
// from the bond list of both atoms

  SAtomBond *AtomBond;
  int nAtomBonds;
  PCAtom atoms[10];
  byte order[10];
  PCAtom pA,pB;

  //char AtomID1[30];
  //char AtomID2[30];
  
  int idel,j;

  for (int nn=0;nn<=1;nn++) {
    if (nn == 0) {
      pA = pa1;
      pB = pa2;
    } else {
      pA = pa2;
      pB = pa1;
    }
    idel = -1;
    j = 0;

    pA->GetBonds ( AtomBond, nAtomBonds);
    if(nAtomBonds>0) {
      //pA->GetAtomID ( AtomID1 );
      //pB->GetAtomID ( AtomID2 );
      //cout << AtomID1 << " " << AtomID2 << " nAtomBonds " << nAtomBonds;
      for (int i=0;i<nAtomBonds;i++) {
	//AtomBond[i].atom->GetAtomID(AtomID2);
        //cout << " " << AtomID2;
        if (AtomBond[i].atom == pB) {
          idel = i;
        } else {
          atoms[j] = AtomBond[i].atom;
          order[j] = AtomBond[i].order;
          j++;
        }
      }
      //cout << endl << "idel " << idel << endl;
      if (idel>=0) {
        pA->FreeBonds();
        for (int i=0;i<nAtomBonds-1;i++) {
          //cout << "adding " << i << " " <<  pA->GetNBonds() << " " <<  atoms[i]->name << endl;
          if ( pA->GetNBonds() == 0 )
            pA->AddBond( atoms[i],order[i],3);
          else 
            pA->AddBond( atoms[i],order[i],1);
        }
      }
    }
  }
  if (idel>=0)
    return 0;
  else
    return 1;
 
}

//-----------------------------------------------------------
void CMolBonds::AddConnection (
    int ia1 , int ia2, PPCAtom selAtom, int offset ) {
//-----------------------------------------------------------
  //cout << "ia1,ia2,offset " << ia1 << " " << ia2 << " " << offset << endl;
  //printf ( "Bond %s %s %i %s %s %i\n",selAtom[ia1]->residue->name,
  //selAtom[ia1]->name,selAtom[ia1]->serNum,selAtom[ia2]->residue->name,selAtom[ia2]->name),selAtom[ia2]->serNum;
  int ja1,ja2;
  ja1 = ia1 + offset;
  ja2 = ia2 + offset;
  //cout << "ja1,ja2 " << ja1 << " " << ja2 << endl;
  if ( fabs(selAtom[ja1]->x - selAtom[ja2]->x) < 0.0001 &&
       fabs(selAtom[ja1]->y - selAtom[ja2]->y) < 0.0001 &&
       fabs(selAtom[ja1]->z - selAtom[ja2]->z)  < 0.0001 ) {
    //cout << "Bonded atoms superposed " <<
    //  GetMolHnd(0)->AtomLabel_atom(selAtom[ja1]) <<
    // " " <<  GetMolHnd(0)->AtomLabel_atom(selAtom[ja2]) << endl;
    return;
  }
  int rv1,rv2;

   if ( selAtom[ja1]->GetNBonds() == 0 ) {
     rv1 = selAtom[ja1]->AddBond(selAtom[ja2],1,3);
   }
   else {
     rv1 = selAtom[ja1]->AddBond(selAtom[ja2],1,1);
   }
   if ( selAtom[ja2]->GetNBonds() == 0 ) {
     rv2 = selAtom[ja2]->AddBond(selAtom[ja1],1,3);
   }
   else {
     rv2 = selAtom[ja2]->AddBond(selAtom[ja1],1,1);
   }
   //cout << "rv " << rv1 << " " << rv2 << endl;
}
//-----------------------------------------------------------
void CMolBonds::AddConnection ( int ia1 , int ia2,
                    PPCAtom selAtom1, PPCAtom selAtom2,
                    int offset1 , int offset2 ) {
//-----------------------------------------------------------
  int ja1 = ia1 + offset1;
  int ja2 = ia2 + offset2;
  //printf ( "Bond %s %s %s %s\n",selAtom1[ja1]->residue->name,
  // selAtom1[ja1]->name,selAtom2[ja2]->residue->name,selAtom2[ja2]->name);

  if ( fabs(selAtom1[ja1]->x - selAtom2[ja2]->x) < 0.0001 &&
       fabs(selAtom1[ja1]->y - selAtom2[ja2]->y) < 0.0001 &&
       fabs(selAtom1[ja1]->z - selAtom2[ja2]->z)  < 0.0001 ) {
    //cout << "Bonded atoms superposed " <<
    //  GetMolHnd(0)->AtomLabel_atom(selAtom1[ja1]) <<
    // " " <<  GetMolHnd(0)->AtomLabel_atom(selAtom2[ja2]) << endl;
    return;
  }


   if ( selAtom1[ja1]->GetNBonds() == 0 ) {
     selAtom1[ja1]->AddBond(selAtom2[ja2],1,3);
   }
   else {
     selAtom1[ja1]->AddBond(selAtom2[ja2],1,1);
   }
   if ( selAtom2[ja2]->GetNBonds() == 0 ) {
     selAtom2[ja2]->AddBond(selAtom1[ja1],1,3);
   }
   else {
     selAtom2[ja2]->AddBond(selAtom1[ja1],1,1);
   }
}

//---------------------------------------------------------------------
int CMolBonds::IntraResContacts ( PCResidue pRes, int nAlt, 
     PPCAtom modelSelAtom[], int nSelAtom[], int firstModel, int lastModel ) {
//---------------------------------------------------------------------
  PPCAtom pAtom1=0;
  PPCAtom pAtom2=0;
  int nAtom1,nAtom2,nContacts,ic;
  PSContact contacts = NULL;
  PCAtom pa1,pa2;

  pAtom1 = NULL;
  pAtom2 = NULL;
  pRes->GetAtomTable1(pAtom1,nAtom1);
  pRes->GetAtomTable1(pAtom2,nAtom2);
  //cout << "IntraResContacts nAtom1 " << nAtom1 << " " << firstModel << " " << lastModel << endl;

  char AtomID1[30];
  char AtomID2[30];

  // One atom in residue - probably water
  if ( nAtom1 == 1 ){
    if(pAtom1) delete [] pAtom1;
    if(pAtom2) delete [] pAtom2;
    return 0;
  }

  GetMolHnd(0)->SeekContacts(pAtom1,nAtom1,pAtom2,nAtom2,
		 0.0,params->intraResCut,0,
		 contacts,nContacts,0,NULL,0);

  //printf ( "nContacts %i\n",nContacts);
  if ( contacts && nContacts > 0 ) { 
    for ( ic = 0; ic < nContacts; ic++) {
      if ( contacts[ic].id1 < contacts[ic].id2 ) {
        pa1 = pAtom1[contacts[ic].id1];
        pa2 = pAtom1[contacts[ic].id2];
        pa1->GetAtomID ( AtomID1 );
        pa2->GetAtomID ( AtomID2 );
        //cout << "IntraResContacts " << AtomID1 << " " << AtomID2 << endl;
        if ( ( nAlt <= 0 || (molHnds[0]->doAltLocMatch( pa1, pa2 )) ) 
            && (ltBondDistance(pa1,pa2,contacts[ic].dist)) ) {
          AddConnection (contacts[ic].id1,contacts[ic].id2,pAtom1,pAtom2);
          if (lastModel>0) {
            for (int nm=firstModel+1;nm<=lastModel;nm++) {
            if (GetMolHnd(0)->GetModel(nm) != NULL) {
              if (nSelAtom[nm] == nAtom1) AddConnection (contacts[ic].id1,contacts[ic].id2,modelSelAtom[nm],modelSelAtom[nm]);
            }}
          }
        } 
      }
    }
  }
  
  if(pAtom1) delete [] pAtom1;
  if(pAtom2) delete [] pAtom2;
  if ( contacts ) delete [] contacts;
  return 0;
}

//---------------------------------------------------------------------
bool CMolBonds::ltBondDistance ( PCAtom pa1, PCAtom pa2, realtype dist ) {
//---------------------------------------------------------------------
  PCLibElement la1,la2;
  realtype rmax;
  if ( (la1 = params->sbase->LibElement(pa1->element)) &&
       (la2 = params->sbase->LibElement(pa2->element)) ) {
    rmax = (la1->maxBondRad + la2->maxBondRad) * params->maxBondRadFactor;
    //printf( "%s %s dist %f rmax %f\n",pa1->name,pa2->name,dist,rmax);
    return dist < rmax;
  }
  else
    return dist < params->maxBondRad;
}

//----------------------------------------------------------------------
bool CMolBonds::isInterResBond ( PCAtom pa1, PCAtom pa2 ) {
//----------------------------------------------------------------------
  int type1,type2,n;
  AtomName a1,a2;

  /*
  std::string first,second;
  first =  std::string(GetMolHnd(0)->AtomLabel_atom(pa1));// Why is this dodgy (well valgrind thinks it is)?
  second = std::string(GetMolHnd(0)->AtomLabel_atom(pa2));
  */

  type1 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa1->residue);
  type2 =  dynamic_cast<PCMMANManager>(molHnds[0])->GetRestypeCode (pa2->residue);

  //cout << "CMolBonds::isInterResBond " << first << " " <<  second << " " << type1 << " " << type2 << " " << params->sbase->nLinks << endl;
  
  strcpy_css(a1,pa1->name);
  strcpy_css(a2,pa2->name);
  
  for ( n = 0; n <= params->sbase->nLinks; n++ ) {
    if (strcmp(params->sbase->link[n]->id,"symmetry") != 0 &&
        strcmp(params->sbase->link[n]->id,"gap") != 0 ) {
      //cout << "CMolBonds::isInterResBond " << n << " " << params->sbase->link[n]->id << endl;
      if ( params->sbase->link[n]->lg1.Match(type1, pa1->residue->name, a1)&&
	  params->sbase->link[n]->lg2.Match(type2, pa2->residue->name, a2)) {
	//printf ("Type of link: %i\n",n);
        return true;
      } else if ( params->sbase->link[n]->lg2.Match(type1, pa1->residue->name, a1)&&
	  params->sbase->link[n]->lg1.Match(type2, pa2->residue->name, a2)) {
	//printf ("Type of link: %i\n",n);
          return true;
      }
    }
  }
  return false;
}



