/*
     mmut/mmut_secstr.cc: CCP4MG Molecular Graphics Program
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


//
// Assign secondary structure based on Kabsch and Sanders
// Biopolymers 22 2577-2637 (1983)
// For all residues that are recognised as amino acids assign an integer
// which represents the secondary structure type (see enum SecStr below)
//

#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <mmdb/mmdb_manager.h>
#include "mmut_manager.h"
#include "mmut_secstr.h"

using namespace std;

//------------------------------------------------------------------
CSecStructure::CSecStructure (PCMMUTManager molHndin) : CMMANBase( molHndin) {
//------------------------------------------------------------------
  InitParams();

}

//------------------------------------------------------------------
CSecStructure::CSecStructure (const PCMMUTManager molHndin, const int selHndin) : CMMANBase( molHndin, selHndin) {
//------------------------------------------------------------------
  InitParams();

}

//------------------------------------------------------------------
CSecStructure::~CSecStructure () {
//------------------------------------------------------------------
  ClearMemory();
}

//------------------------------------------------------------------
void CSecStructure::InitParams () {
//------------------------------------------------------------------
  //Parameters
  NOmaxdist = 3.5;
  NOmaxdist2 = 12.25;
  NOCanglemin = PI/2;
  flagBulge = 1;

  //Output data 
  // NB nres used as a flag to indicate if data is currently set
  nRes = -1;
  secstr = NULL;
  hbonds = NULL;

}

//------------------------------------------------------------------
void CSecStructure::SetParams (int nv,double *value, int niv,int *ivalue) {
//------------------------------------------------------------------
  //Parameters
  NOmaxdist = value[0];
  NOmaxdist2 =NOmaxdist* NOmaxdist;
  NOCanglemin =value[1]* PI/180;
  flagBulge = ivalue[0];
  //cout <<  "CSecStructure::SetParams NOmaxdist 2NOCanglemin flagBulge "
    //   << NOmaxdist2 << " " << NOCanglemin << " " << flagBulge << endl;
  ClearMemory();
}

//------------------------------------------------------------------
int CSecStructure::GetSecondaryStructure ( int &nresout,
		  ivector &secstrout, imatrix &hbondsout, int imodel ) {
//------------------------------------------------------------------
  int RC;
  if (!secstr || !hbonds || nRes < 0 ) 
       RC = CalculateSecondaryStructure(imodel);
  nresout = nRes;
  secstrout = secstr;
  hbondsout = hbonds;
  return RC;
}

//--------------------------------------------------------------------
int **CSecStructure::GetHBonds (int imodel) {
//--------------------------------------------------------------------
  int RC;
  if (!secstr || !hbonds || nRes < 0 ) 
     RC = CalculateSecondaryStructure(imodel);
  return hbonds;
}
//--------------------------------------------------------------------
PPCAtom *CSecStructure::GetHBondAtoms (int imodel) {
//--------------------------------------------------------------------
  int RC;
  if (!secstr || !hbonds || nRes < 0 ) 
         RC = CalculateSecondaryStructure(imodel);
  return hbond_atoms;
}
//--------------------------------------------------------------------
int *CSecStructure::GetSecStr (int imodel) {
//--------------------------------------------------------------------
  int RC;
  if (!secstr || !hbonds || nRes < 0 ) 
            RC = CalculateSecondaryStructure(imodel);
  return secstr;
}



//------------------------------------------------------------------
int CSecStructure::SetFlagBulge ( int flag ) {
//------------------------------------------------------------------
  if ( flag < 0 || flag > 1 ) return 1;
  if ( flag == flagBulge ) return 0;
  flagBulge = flag;
  nRes = -1;
  return 0;
}

//------------------------------------------------------------------
int CSecStructure::InitMemory( int nres ) {
//------------------------------------------------------------------
  int i,j;
  // First clear memory
  ClearMemory();
  nRes = nres;
  // Reserve memory for hbonds and secstr 
  GetMatrixMemory ( hbonds, nres, 3, 0, 0);
  hbond_atoms = new PPCAtom[nres];
  hbondsN = nres;
  for ( i = 0; i < nres; i++) {
     hbond_atoms[i] = new PCAtom[6];
     for ( j = 0; j < 6; j++) hbond_atoms[i][j] = 0;
     for ( j = 0; j <= 2; j++) hbonds[i][j] = 0;
  }
  GetVectorMemory ( secstr, nres, 0);
  for ( i = 0; i < nres; i++) secstr[i] = 0;
  return 0;
  
}


//------------------------------------------------------------------
void CSecStructure::ClearMemory ( ) {
//------------------------------------------------------------------
  nRes = -1;
  if (secstr ) FreeVectorMemory ( secstr, 0);
  if (hbonds ) FreeMatrixMemory (hbonds,hbondsN,0,0);
  secstr = NULL;
  hbonds = NULL;
}    



//------------------------------------------------------------------
int CSecStructure::CalculateSecondaryStructure (int imodel ) {
//------------------------------------------------------------------
// Define a secondary structure type of each amino acid residue in the 
// structure.  
// Procedure:
// Find all amino acids
// Find all pairs of amino acids which have inter-Ca distance  < 10.0A
// Test for hydrogen bonds between the main chain N and O of the close residues
//   and store the information in the hbonds matrix
// Analyse the info in hbonds matrix to assign secondary structure to secstr vector

  PPCResidue    selRes;
  int		nres,nres2;
  PSContact     contact = NULL;
  int           ncontacts;
  int		ir1, ir2, irdif;
  int		i,j,k,l;
  int           namino = 0; 

  PCMMUTManager molH = GetMolHnd(0);

  // Bware the non-existant model
  if ( molH->GetModel(imodel) == NULL) return 1;

  // Get the selected residues
  GetSelection ( selRes, nres, imodel);
  if ( nres <= 0 ) return 1;
  //cout << "CalculateSecondaryStructure " << imodel << " " << nres << " " << selRes[0]->GetModel()->GetSerNum() << endl;

  // Create array of pointers to CA atoms in amino acids 
  PPCAtom selCa = new PCAtom[nres];
  for ( i = 0; i < nres; i++ ) {
    if ( selRes[i]->isAminoacid() ) {
      namino++; 
      selCa[i] = selRes[i]->GetAtom("CA", " C", "*");
    }
    else 
      selCa[i] = NULL;
  }
  if ( namino < 4 ) return 1;
        
  // Second copy of the same data
  nres2 = nres;
  PPCAtom selCa2 = new PCAtom[nres];
  for ( i = 0; i < nres; i++ ) selCa2[i] = selCa[i];

  // Find all close Ca's - i.e. find the contacts between the two
  // equivalent sets of Ca atoms
  molH->SeekContacts ( selCa,nres, selCa2,nres2,
               2.0,10.0,2,contact,ncontacts, 0 );

  //printf ("Number of contacts %d\n", ncontacts);
  if ( ncontacts <= 0 ) return 1;

  // Init memory for data
  InitMemory(nres);

  // Loop over all close (in space) residues - excluding those
  // that are close in sequence

  for ( i=0;i<ncontacts;i++) {
    if ( (irdif = (ir1 = contact[i].id2) - (ir2 = contact[i].id1) ) > 2  ) {
  // Test if there is donor Hbond from residue ir1
      if ( IsHBond ( selRes[ir1], selRes[ir2] ) ) {
        k=0;
        while ( hbonds[ir1][k] != 0  && k < 2 ) k++;
        hbonds[ir1][k] = -irdif;
	hbond_atoms[ir1][k] = selRes[ir1]->GetAtom( "N");
	hbond_atoms[ir1][k+3] = selRes[ir2]->GetAtom( "O");
      }
  // Test if there is donor Hbond from residue ir2
      if ( IsHBond ( selRes[ir2], selRes[ir1] ) ) {        
	k=0;
        while ( hbonds[ir2][k] != 0  && k < 2 ) k++;
        hbonds[ir2][k] = irdif;
	hbond_atoms[ir2][k] = selRes[ir2]->GetAtom( "N");
	hbond_atoms[ir2][k+3] = selRes[ir1]->GetAtom( "O");
      }
    }
  }

  //  Assign the turns - if there is bifurcated bond then the 4-turn takes
  // precedence - read the paper to make sense of this

  for ( i=0; i < nres; i++ ) { 
    k = 0;
    while ( k <= 2 &&  hbonds[i][k] != 0 ) {
      if (hbonds[i][k] == -5  ) {
	secstr[i-1] = TURN5;
	secstr[i-2] = TURN5;
	secstr[i-3] = TURN5;
	secstr[i-4] = TURN5;
      }
      if (hbonds[i][k] == -3  ) {
	secstr[i-1] = TURN3;
	secstr[i-2] = TURN3;
      }
      k++;
    }
  }
  for ( i=0; i < nres; i++ ) {
    k = 0;
    while ( k <= 2 &&  hbonds[i][k] != 0 ) {
      if (hbonds[i][k] == -4  ) {
        secstr[i-1] = TURN4;
        secstr[i-2] = TURN4;
        secstr[i-3] = TURN4;
      }
      k++;
    }
  }

  // Look for consecutive 4-turns which make alpha helix

  for ( i = 1; i < (nres - 3); i++ ) {
    if ( ( secstr[i] == ALPHA || secstr[i] == TURN4 ) && 
	 ( secstr[i+1] == ALPHA || secstr[i+1] == TURN4 ) &&
	 ( secstr[i+2] == ALPHA || secstr[i+2] == TURN4 ) &&
	 ( secstr[i+3] == ALPHA || secstr[i+3] == TURN4 ) ) {
        for ( j = i; j <= i+3; j++ ) secstr[j] = ALPHA;
    }
  }

  for (i=0; i < nres; i++ ) {
    k = 0;
    while (  k <= 2 && hbonds[i][k] != 0 ) {
      irdif = hbonds[i][k];
  // Test for 'close' hbond 
      j = i + irdif;
      l = 0;
      while (  l <= 2 && hbonds[j][l] != 0 ) {
      // Antiparallel strands
        if ( hbonds[j][l] == -irdif ) {
          secstr[i] = BETA;
          secstr[j] = BETA;
        }
      // Parallel strand
        if ( hbonds[j][l] == -irdif-2 ) {
          secstr[i-1] = BETA;
          secstr[j] = BETA;
        }
      // Parallel beta bulge
        if ( hbonds[j][l] == -irdif-3 ) {
          if (flagBulge ) {
            if (secstr[i-1] == NOSECSTR) secstr[i-1] = BULGE;
            if (secstr[i-2] == NOSECSTR) secstr[i-2] = BULGE;
            if (secstr[j] == NOSECSTR) secstr[j] = BULGE;
          }
          else {
            if (secstr[i-1] == NOSECSTR) secstr[i-1] = BETA;
            if (secstr[i-2] == NOSECSTR) secstr[i-2] = BETA;
            if (secstr[j] == NOSECSTR) secstr[j] = BETA;
          }
        }
        l++;
      }
   // Test for 'wide' hbond
      if ( (j = i + hbonds[i][k] +2) < nres ) {
        l = 0;
        while ( l <= 2  && hbonds[j][l] != 0 ) {
        // Antiaprallel strands
          if ( hbonds[j][l] == -irdif-4 ) {
             secstr[i-1] = BETA;
             secstr[j-1] = BETA;
          }
          // Parallel strands
          if ( hbonds[j][l] == -irdif-2 ) {
            secstr[i] = BETA;
	    secstr[j-1] = BETA;
          }
          l++;
        }
      }

   // test for anti-parallel B-bulge between 'close' hbonds
      if ( (j = i + hbonds[i][k] -1) >= 0 ) {
        l = 0;
        while ( l <= 2 &&  hbonds[j][l] != 0 ) {
          if ( hbonds[j][l] == -irdif +1 ) {
            if ( flagBulge ) {
	      if (secstr[i] == NOSECSTR) secstr[i] = BULGE;
	      if (secstr[j+1] == NOSECSTR) secstr[j+1] = BULGE;
	      if (secstr[j] == NOSECSTR) secstr[j] = BULGE;
            }
            else {
              if (secstr[i] == NOSECSTR) secstr[i] = BETA;
              if (secstr[j+1] == NOSECSTR) secstr[j+1] = BETA;
              if (secstr[j] == NOSECSTR) secstr[j] = BETA;
            } 
          }
          l++;
        }
      }

    // test for anti-parallel B-bulge between 'wide' hbonds
      if ( (j = i + hbonds[i][k] + 3 ) < nres ) {
        l = 0;
        while ( l <= 2 &&  hbonds[j][l] != 0 ) {
          if ( hbonds[j][l] == -irdif + 5) {
            if ( flagBulge ) {
              if (secstr[i-1] == NOSECSTR) secstr[i-1] = BULGE;
              if (secstr[j-1] == NOSECSTR) secstr[j-1] = BULGE;
              if (secstr[j-2] == NOSECSTR) secstr[j-2] = BULGE;
            }
            else {
              if (secstr[i-1] == NOSECSTR) secstr[i-1] = BETA;
              if (secstr[j-1] == NOSECSTR) secstr[j-1] = BETA;
              if (secstr[j-2] == NOSECSTR) secstr[j-2] = BETA;
            }
          }
         // and bulge in parallel strand
          else if ( hbonds[j][l] == -irdif -3 ) {
	    if ( flagBulge ) {
              if ( secstr[i] == NOSECSTR) secstr[i] = BULGE;
              if (secstr[j-1] == NOSECSTR) secstr[j-1] = BULGE;
              if (secstr[j-2] == NOSECSTR) secstr[j-2] = BULGE;
            }
            else {
              if ( secstr[i] == NOSECSTR) secstr[i] = BETA;
              if (secstr[j-1] == NOSECSTR) secstr[j-1] = BETA;
              if (secstr[j-2] == NOSECSTR) secstr[j-2] = BETA;
            }
          }
          l++;
        }
      }
      k++;
    } // Finish looping over Hbonds for residue (k loop)
  }  // Finish looping over residues ( i loop)

  //for ( i = 0; i <=20; i++ ) {
  //printf ("CalcSec Res %i secstr %i \n", i, secstr[i] );
  //}

  //ClearSelection();

  delete [] selCa;
  delete [] selCa2;

  return 0;
 
}

//-----------------------------------------------------------------------
std::string CSecStructure::Print (int imodel) { 
//-----------------------------------------------------------------------
  int		i,k;
  const char 	secstr_text [7][12] = { "           ",
                                        "beta strand",
                                        "beta bulge ",
			                "3-turn     ",
                                        "4-turn     ",
                                        "5-turn     ",
                                        "alpha helix"};
  std::ostringstream output;
  int           nr;
  PPCResidue    selRes;
  PCResidue     j;
  std::string resid;
    int first_model=1;
  int last_model=1;

  PCMMUTManager molH = GetMolHnd(0);

  if (imodel == 0 && molH->GetNumberOfModels() > 1) {
    first_model = 1;
    last_model =  molH->GetNumberOfModels();
  }
  else if (imodel>0) {
    first_model = imodel;
    last_model = imodel;
  }
  output << "Secondary Structure Assignment based on Kabsch & Sanders"
         << endl
         << "Residue        "
         << "SecStruct   "
         << "Hydrogen bonded to.."
         << endl;


  for (int im=first_model;im<=last_model;im++) {
  // Bware the non-existant model
  if ( molH->GetModel(im) != NULL) {


    if (im != 1) output << endl << "For model number " << im << endl;


    CalculateSecondaryStructure(im);
    //Get all residues with atoms in this atom selection
    GetSelection( selRes, nr, im );
  
 
    for ( i = 0; i < nr; i++) {
      if ( selRes[i]->isAminoacid() ) {
        resid = molH->AtomLabel_residue1(selRes[i]);
        output << endl 
             <<  resid <<  std::setw(15-resid.length()) <<" " 
             << secstr_text[secstr[i]] << " ";  
        if ( hbonds ) {
          k = 0;
          while (  k <= 2 && hbonds[i][k] != 0 ) {
            j = selRes[i + hbonds[i][k]];
            resid = molH->AtomLabel_residue1(j);
            output << resid << std::setw(15-resid.length()) << " ";
            k++;
          }
        }
      }
    }
    output << endl; 
  } }
  return output.str();
}
 
//-----------------------------------------------------------------------
Boolean CSecStructure::IsHBond ( PCResidue PCRes1, PCResidue PCRes2 ) {
//-----------------------------------------------------------------------
  PCAtom        NAtom, OAtom , CAtom;
  realtype	dx,dy,dz;
   
// This probably need the option of supporting alternative criteria

  // Test if there is main chain Hbond between PCRes1 (donor) and
  // PCRes2 (acceptor) 
  // As defined Kabsch & Sanders

  if ( (NAtom = PCRes1->GetAtom( "N")) != NULL &&
       (OAtom = PCRes2->GetAtom( "O")) != NULL &&
       (CAtom = PCRes2->GetAtom( "C")) != NULL &&
       ( fabs (dx = (NAtom->x - OAtom->x ))) < NOmaxdist &&
       ( fabs (dy = (NAtom->y - OAtom->y ))) < NOmaxdist &&
       ( fabs (dz = (NAtom->z - OAtom->z ))) < NOmaxdist &&
       	dx*dx+dy*dy+dz*dz <= NOmaxdist2 &&
            GetMolHnd(0)->BondAngle ( NAtom, OAtom, CAtom ) >= NOCanglemin ) {
    return 1;
  } else {
    return 0;
  }
}
