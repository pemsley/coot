/*
     mmut/mman_manager.cc: CCP4MG Molecular Graphics Program
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


#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#if defined(__sgi) || defined(sgi) || defined(__OSF1__) || defined(__osf__)
#include <string.h>
#else
#include <cstring>
#endif
#include <mmdb/mmdb_manager.h>
#include "mman_manager.h"
#include "mmut_manager.h"
#include "mmut_bonds.h"
#include "mmut_sbase.h"
#include "cartesian.h"
#include "matrix.h"
#include "mmdb/mmdb_tables.h"

using namespace std;


//  =====================   CMMANManager   =======================


CMMANManager::CMMANManager(CMMDBManager* molHnd){
  p_bonds = 0;
  p_sas = 0;
  loaded_charge = "None";
  p_sbase = 0;
  customResTypes.clear();
  Copy(molHnd,MMDBFCM_All);
}

CMMANManager::CMMANManager(){
  p_bonds = 0;
  p_sas = 0;
  loaded_charge = "None";
  p_sbase = 0;
  customResTypes.clear();
}

int CMMANManager::SetSBaseAndBondParams(PCMGSBase p_sbase_in,
           PCMolBondParams p_bond_params_in)  {
  p_bond_params = p_bond_params_in;
  p_sbase = p_sbase_in;
  /*
  if (!p_sbase) {
    p_sbase = new CMGSBase();
    p_sbase->InitSBase();
  }

  if (!p_bond_params) p_bond_params = new CMolBondParams(p_sbase);
  //printf ("MMAN Sbase maxAtomInRes %i\n",p_bond_params->sbase->maxAtomInRes);
  */

  udd_sbaseCompoundID = -1;
  udd_sbaseAtomOrdinal = -1;
  udd_atomEnergyType = -1;

  p_sas = NULL;
  udd_atomSAS = -1;
  udd_resSAS = -1;

  p_bonds = NULL;

  int mask[20] = {0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
  for (int i=0;i<20;i++) { label_mask[i] = mask[i]; }

  loaded_charge = "None";

  isTransformed = false;
  transform_com_set = false;
  unremediated = false;
  Mat4Init(current_transform);
  return 0;
}

//------------------------------------------------------------------- 
CMMANManager::CMMANManager(PCMGSBase p_sbase_in,
           PCMolBondParams p_bond_params_in) : CMMUTManager()  {
//-------------------------------------------------------------------
  SetSBaseAndBondParams(p_sbase_in,p_bond_params_in);
}

//-------------------------------------------------------------------------
CMMANManager::~CMMANManager()  { 
//-------------------------------------------------------------------------
  if (p_bonds) delete p_bonds;
  //cout << " CMMANManager destructor" << endl;
}


//-----------------------------------------------------------------------
PCSASArea CMMANManager::GetSASArea ( int selHndin ) {
//-----------------------------------------------------------------------
  if ( !p_sas ) 
    p_sas = new CSASArea( dynamic_cast<PCMMUTManager>(this),selHndin );

  udd_atomSAS = GetUDDHandle( UDR_ATOM,"atom_sas");
  if (udd_atomSAS <= 0 )
    udd_atomSAS = RegisterUDReal ( UDR_ATOM,"atom_sas" );
  //cout << "udd_atomSAS " << udd_atomSAS << endl;
  if ( udd_atomSAS <= 0 ) {
    printf ( "ERROR registering UDD atom SAS data.\n" );
    return NULL;
  }

  udd_resSAS = GetUDDHandle( UDR_RESIDUE,"residue_sas");
  if (udd_resSAS <= 0 )
    udd_resSAS = RegisterUDReal ( UDR_RESIDUE,"residue_sas" );
  if ( udd_resSAS <= 0 ) {
    printf ( "ERROR registering UDD residue SAS data.\n" );
    return NULL;
  }
  p_sas->SetUDD(udd_atomSAS,udd_resSAS);
  return p_sas;
}

//-----------------------------------------------------------------------
CMGSBase *CMMANManager::GetMGSBase() {
//------------------------------------------------------------------------
return p_sbase;
}

//-----------------------------------------------------------------------
std::string CMMANManager::GetMolBonds (std::string monlib_file) {
//-----------------------------------------------------------------------
  int nres,natoms,i;
  PPCAtom atomTable = NULL;
  PPCResidue selRes = NULL;
  std::string output;
  int selHnd = -1;

  //cout << "monlib_file " << monlib_file << endl;
  if (monlib_file.size()>0 ) {
    const char* lib1 = monlib_file.c_str();
    char lib[500];
    strcpy(lib,lib1);
    if (monlib.size()>0) {
      //out << "clearing monlib" << endl;
      //LoadedPCSBStructure_iter p = monlib.find(resn);
      //if (p!=monlib.end()) return p->second;
      monlib.clear();
    } 
    int nmon = p_sbase->LoadMonomerLibrary(lib,monlib);
    if (nmon>0) {
      output += "Using monomer library for this model: ";
      output += monlib_file;
      output +="\n\n";
    }
    //cout << "from LoadMonomerLibrary " << nmon << " " << monlib.size() <<endl;
  }

  if ( selHnd > 0 ) 
    GetSelIndex(selHnd,atomTable,natoms);
  else
    GetAtomTable1(atomTable,natoms);

  if ( udd_sbaseCompoundID < 0 ) {
    udd_sbaseCompoundID = RegisterUDString(UDR_RESIDUE,"sbaseCompoundID" );
    if ( udd_sbaseCompoundID < 0 ) {
      printf ( "ERROR registering sbase CompoundID data.\n" );
      return "";
    }
  }

  if (udd_sbaseAtomOrdinal < 0) {
    udd_sbaseAtomOrdinal = RegisterUDInteger(UDR_ATOM,"sbaseAtomOrdinal" );
    if ( udd_sbaseAtomOrdinal < 0 ) {
      printf ( "ERROR registering sbase AtomOrdinal data.\n" );
      return "";
    }
  }

  if (udd_atomEnergyType < 0) {
    udd_atomEnergyType = RegisterUDInteger(UDR_ATOM,"atomEnergyType" );
    if ( udd_atomEnergyType < 0 ) {
       printf ( "ERROR registering atom EnergyType data.\n" );
       return "";
    }
  }


  for (i=0;i<natoms;i++)
         atomTable[i]->PutUDData(udd_sbaseAtomOrdinal,-1);
  for (i=0;i<natoms;i++)
         atomTable[i]->PutUDData(udd_atomEnergyType,0);

  
  //int selH = NewSelection();
  //Select(selH,STYPE_ATOM,"/*/*/(ASP,ASN,LYS,TYR,LEU,GLU,GLN,HIS,SER,PHE,TRP,PRO,ARG,MET,CYS)/1HB,1HD1,1HE2,1HH1",SKEY_NEW);
  //Select(selH,STYPE_ATOM,"/*/*/(A,C,T,G,U)/C2*",SKEY_OR);
  //int nat;
  //PPCAtom selat;
  //GetSelIndex(selH, selat, nat );
  //cout << "unremediated nat " << nat << endl;
  //if (nat>0) { 
  //  unremediated = true;
  //  output += "Treating this as UNREMEDIATED coordinate file\n";
  //} else {
  //   output += "Treating this as REMEDIATED coordinate file\n";
  //}
  //DeleteSelection(selH);
  


  if (selHnd < 0 ) {
    GetResidueTable(selRes,nres);
  } else {
    int resSelHnd;
    resSelHnd = NewSelection();
    Select (resSelHnd,STYPE_RESIDUE,selHnd,SKEY_OR); 
    GetSelIndex ( resSelHnd, selRes, nres );
  }

  

  /*
  std::map<std::string,std::string>::iterator q = customResSynonym.begin();
  while (q !=  customResSynonym.end()) {
    //cout << "GetMolBonds " << q->first.c_str() << " " <<  q->second.c_str()<< endl;
    q++;
  }
  */



  for (i=0;i<nres;i++) {
    output += p_sbase->AssignAtomType(selRes[i], monlib,customResSynonym,
     udd_sbaseCompoundID,udd_sbaseAtomOrdinal, udd_atomEnergyType, unremediated);
      //cout << p_sbase->ListAtomType(dynamic_cast<PCMMANManager>(this),
	//         selRes[i], udd_sbaseCompoundID,
	//	udd_sbaseAtomOrdinal, udd_atomEnergyType);
  }

  if (selHnd<0) {
    if (selRes) delete [] selRes;
    if (atomTable) delete [] atomTable;
  }
  
  if (!p_bonds)
    p_bonds = new CMolBonds(dynamic_cast<PCMMUTManager>(this), p_bond_params );
  
  output += p_bonds->FindBonds( udd_sbaseCompoundID,
		udd_sbaseAtomOrdinal, udd_atomEnergyType);
  //printf ("init pbonds %i\n",p_bonds);
  
  return output;
  
}

//-----------------------------------------------------------------------
int CMMANManager::EditBonds (int mode, PCAtom p_atom1, PCAtom p_atom2) {
//-----------------------------------------------------------------------
  if (!p_bonds) return 1;
  if (mode < 0) 
    return p_bonds->DeleteConnection(p_atom1,p_atom2);
  else {
    // Before attempting to add a bond make sure that it is not
    // there already
    SAtomBond *AtomBond;
    int nAtomBonds;
    p_atom1->GetBonds ( AtomBond, nAtomBonds);
    for (int i=0;i<nAtomBonds;i++) {
      if (AtomBond[i].atom == p_atom2) return 2;
    }
    p_bonds->AddConnection(p_atom1,p_atom2);
  }
  return 0;
}

//---------------------------------------------------------------------
int CMMANManager::SetupAtomEnergyTypes () {
//---------------------------------------------------------------------
  PPCResidue resTable = NULL;
  PPCAtom atomTable = NULL;
  int i,j,nres,nat;
  PCResidue p_res;
  pstr compoundID = NULL;
  PCSBStructure p_sbase_struct;
  pstr atomType;
  int atomOrdinal;
  PCSBAtom p_sbase_atom;

  udd_atomEnergyType = RegisterUDInteger(UDR_ATOM,"atomEnergyType" );
  if ( udd_atomEnergyType < 0 ) {
    printf ( "ERROR registering atomEnergyType data.\n" );
    return 1;
  }


  // Loop over ever residue in the molecule 
  GetResidueTable(resTable,nres);
  
  for (i=0;i<nres;i++) {
    p_res = resTable[i];
    // Find the sbase structure class for this residue 
    p_res->GetUDData(udd_sbaseCompoundID,compoundID);
    atomTable = NULL;
    p_res->GetAtomTable1(atomTable,nat);
    if( (p_sbase_struct = p_sbase->GetStructure(compoundID,monlib))) {
      //printf("p_sbase_struct %i\n",p_sbase_struct);
      // loop over all atoms in the residue to find its energy type
      for (j=0;j<nat;j++) {
        atomTable[j]->GetUDData(udd_sbaseAtomOrdinal,atomOrdinal);
        //printf("compound %s j %i atom %s ordinal %i\n",compoundID,j,atomTable[j]->name,atomOrdinal);
        if( atomOrdinal >= 0 &&
            (p_sbase_atom = p_sbase_struct->Atom[atomOrdinal]) )
          atomType = p_sbase_atom->energyType;
        else
          atomType = "";
    
	//printf ( "resType %s index %i atomType %s %i\n"
	    //  ,compoundID,atomOrdinal,atomType,strlen(atomType));
        if ( strlen(atomType) <  1 ) 
          //atomTable[j]->PutUDData(udd_atomEnergyType,-1);
          atomTable[j]->PutUDData(udd_atomEnergyType,
		p_sbase->LibAtom("",atomTable[j]->element));
        else 
          atomTable[j]->PutUDData(udd_atomEnergyType,
                              p_sbase->LibAtom(atomType));
	
      }  // End of loop over atoms
    } 
    // Residue type not recognised so assign a default type for the
   // atom element type    
    else { 
      for (j=0;j<nat;j++) {
        atomTable[j]->PutUDData(udd_atomEnergyType,
		p_sbase->LibAtom("",atomTable[j]->element));
        //printf ( "element %s type %i\n",atomTable[j]->element,
          //      p_sbase->LibAtom(1,atomTable[j]->element));
      }
    }
  }        // End of loop over residues
  return 0; 
}

//--------------------------------------------------------------------------
int CMMANManager::GetAtomTypeData ( int selHnd, int type,
                                rvector &dataout, int &nat ) {
//--------------------------------------------------------------------------

  PPCAtom atomTable;
  int i;
  int atomOrdinal;

  if ( selHnd > 0 ) 
    GetSelIndex(selHnd,atomTable,nat);
  else 
    GetAtomTable1(atomTable,nat);

  if ( dataout ) FreeVectorMemory(dataout,0);
  dataout = NULL;
  GetVectorMemory(dataout,nat,0);

  switch (type) {
  case VDWRADIUS: 
    for (i=0;i<nat;i++) {
      atomTable[i]->GetUDData(udd_atomEnergyType, atomOrdinal);
      dataout[i] =p_sbase->libAtom[atomOrdinal].vdwRadius;
    }
    break;
  case VDWHRADIUS:
    for (i=0;i<nat;i++) {
      atomTable[i]->GetUDData(udd_atomEnergyType, atomOrdinal);
      dataout[i] =p_sbase->libAtom[atomOrdinal].vdwHRadius;
    }
    break;
  case IONRADIUS:
    for (i=0;i<nat;i++) {
      atomTable[i]->GetUDData(udd_atomEnergyType, atomOrdinal);
      dataout[i] =p_sbase->libAtom[atomOrdinal].ionRadius;
    }
    break;
  }
  return 0;
}
//--------------------------------------------------------------------------
std::vector<double> CMMANManager::GetAtomRadii ( int selHnd, int type, double scale ) {
//--------------------------------------------------------------------------
/*
Return an array of radii for selected atoms - expected to be used in 
drawing spheres etc.
*/

  PPCAtom atomTable=NULL;
  int nat,i;
  int atomOrdinal;
  std::vector<double> atomRadii;

  //cout << "scale " << scale << endl;
  if ( selHnd > 0 ) 
    GetSelIndex(selHnd,atomTable,nat);
  else 
    GetAtomTable1(atomTable,nat);

  atomRadii = std::vector<double> (nat);

  switch (type) {
  case VDWRADIUS: 
    for (i=0;i<nat;i++) {
      atomTable[i]->GetUDData(udd_atomEnergyType, atomOrdinal);
      atomRadii[i] =p_sbase->libAtom[atomOrdinal].vdwRadius * scale;
      //cout << "atomRadii " << atomRadii[i] << endl;
    }
    break;
  case VDWHRADIUS:
    for (i=0;i<nat;i++) {
      atomTable[i]->GetUDData(udd_atomEnergyType, atomOrdinal);
      atomRadii[i] =p_sbase->libAtom[atomOrdinal].vdwHRadius * scale;
    }
    break;
  case IONRADIUS:
    for (i=0;i<nat;i++) {
      atomTable[i]->GetUDData(udd_atomEnergyType, atomOrdinal);
      atomRadii[i] =p_sbase->libAtom[atomOrdinal].ionRadius * scale;
    }
    break;
  }
  /*
  for (i=0;i<nat;i++) {
    cout << "atomRadii " << atomRadii[i] << endl;
  }
  */
  return atomRadii;
}

//--------------------------------------------------------------------
int CMMANManager::GetAtomEnergyType(PCAtom p_atom) {
//--------------------------------------------------------------------
  int atomType;
  p_atom->GetUDData(udd_atomEnergyType, atomType);
  return atomType;
}

//---------------------------------------------------------------------
realtype CMMANManager::GetAtomVDWRadius(PCAtom p_atom) {
//---------------------------------------------------------------------
  int atomType;
  p_atom->GetUDData(udd_atomEnergyType, atomType);
  return p_sbase->libAtom[atomType].vdwRadius;
  //return 1.5;
}

//--------------------------------------------------------------------------
int CMMANManager::GetAtomTypeData ( int selHnd, int type,
                                ivector &dataout, int &nat ) {
//--------------------------------------------------------------------------
  PPCAtom atomTable;
  int i;
  int atomOrd;

  if ( selHnd > 0 ) 
    GetSelIndex(selHnd,atomTable,nat);
  else 
    GetAtomTable1(atomTable,nat);

  if ( dataout ) FreeVectorMemory(dataout,0);
  dataout = NULL;
  GetVectorMemory(dataout,nat,0);

  switch (type) {
  case HBTYPE: 
    for (i=0;i<nat;i++) {
      atomTable[i]->GetUDData(udd_atomEnergyType, atomOrd);
      dataout[i] =p_sbase->libAtom[atomOrd].hbType;
    }
    break;
  }
  return 0;
}

//---------------------------------------------------------------------
int CMMANManager::GetAtomHBondType1(PCAtom p_atom) {
//---------------------------------------------------------------------
  int atomOrd;
  std::string name;
  p_atom->GetUDData(udd_atomEnergyType, atomOrd);
  if(atomOrd<0) return HBTYPE_UNKNOWN;
  //name = AtomLabel(p_atom);
  return p_sbase->libAtom[atomOrd].hbType;
}
//---------------------------------------------------------------------
char* CMMANManager::GetAtomHBondType(PCAtom p_atom) {
//---------------------------------------------------------------------
  int atomOrd;
  p_atom->GetUDData(udd_atomEnergyType, atomOrd);
  if (atomOrd >= 0) {
    return p_sbase->libAtom[atomOrd].getHBType();
  }
  else {
    return "U";
  }
}
//---------------------------------------------------------------------
int CMMANManager::LoadCharge(std::string loadfrom ) {
//---------------------------------------------------------------------
  PPCAtom atomTable = NULL;
  int i,nat,udd;
  int atomOrd,atomEtype;
  PCSBStructure pSbaseRes=NULL;
  bool useUDD=0;


  if (useUDD) {
    udd = GetUDDHandle( UDR_ATOM,"atom_charge");
    if (udd <= 0 )
      udd = RegisterUDReal ( UDR_ATOM,"atom_charge" );
    //cout << "udd " << udd << endl;
    if ( udd <= 0 ) {
      printf ( "ERROR registering UDD atom charge data.\n" );
      return 1;
    }
  
    if ( loadfrom == "partial_charge" && loaded_charge != "partial_charge") {
      // Get charge info from the monomer library description of the residue
      GetAtomTable1(atomTable,nat);
      for (i=0;i<nat;i++) {
        atomTable[i]->GetUDData(udd_sbaseAtomOrdinal,atomOrd);
        if (atomOrd>=0) {
          pSbaseRes = p_sbase->GetStructure(atomTable[i]->residue->name,monlib);
          if (pSbaseRes) {
            atomTable[i]->PutUDData(udd,pSbaseRes->Atom[atomOrd]->ccp4_charge);
          } else {
            atomTable[i]->PutUDData(udd,0.0);
          }
        } else
          atomTable[i]->PutUDData(udd,0.0);
      }
      loaded_charge = "partial_charge";
    } else if ( loadfrom == "surface_potential_charge" && 
              loaded_charge != "surface_potential_charge") {
      // Get the charge implied by the atom type 
      // Read from ccp4mg/data/ener_lib.cif labelled surface_potential_charge
      GetAtomTable1(atomTable,nat);
      for (i=0;i<nat;i++) {
        atomTable[i]->GetUDData(udd_atomEnergyType,atomEtype);
        atomTable[i]->PutUDData(udd,p_sbase->libAtom[atomEtype].charge);
      }
      loaded_charge =  "surface_potential_charge";
    }
  } else {
    if ( loadfrom == "partial_charge" && loaded_charge != "partial_charge") {
      // Get charge info from the monomer library description of the residue
      GetAtomTable1(atomTable,nat);
      for (i=0;i<nat;i++) {
        atomTable[i]->GetUDData(udd_sbaseAtomOrdinal,atomOrd);
        if (atomOrd>=0) {
          pSbaseRes = p_sbase->GetStructure(atomTable[i]->residue->name,monlib);
          if (pSbaseRes) {
            atomTable[i]->charge=pSbaseRes->Atom[atomOrd]->ccp4_charge;
          } else {
            atomTable[i]->charge=0.0;
          }
        } else
          atomTable[i]->charge=0.0;
      }
      loaded_charge = "partial_charge";
    } else if ( loadfrom == "surface_potential_charge" && 
              loaded_charge != "surface_potential_charge") {
      // Get the charge implied by the atom type 
      // Read from ccp4mg/data/ener_lib.cif labelled surface_potential_charge
      GetAtomTable1(atomTable,nat);
      for (i=0;i<nat;i++) {
        if (atomTable[i]->charge < -0.001 || atomTable[i]->charge > 0.001 ) {
          if (atomTable[i]->charge > 98.0) atomTable[i]->charge=0.0;
        } else {
          atomTable[i]->GetUDData(udd_atomEnergyType,atomEtype);
          atomTable[i]->charge=p_sbase->libAtom[atomEtype].charge;
        }
      }
      loaded_charge =  "surface_potential_charge";
    }
  }
  return 0;

}

//-------------------------------------------------------------------- 
std::string CMMANManager::PrintCharges(void) {
//-------------------------------------------------------------------- 

  std::ostringstream output;
  //char AtomID[30];
  std::string AtomID;
  bool useUDD =0;
  realtype charge;
  PPCAtom atomTable = NULL;
  int i,nat,udd;

  if (useUDD) {
    udd = GetUDDHandle( UDR_ATOM,"atom_charge");
    if (udd <= 0 ) {
      output << "No charges loaded";
      return output.str();
    }
  }

  output << "Loaded charge for " << loaded_charge << endl;
  GetAtomTable1(atomTable,nat);

  for (i=0;i<nat;i++) {
    AtomID=AtomLabel(atomTable[i]);
    //atomTable[i]->GetAtomID(AtomID);
    if (useUDD) {
      atomTable[i]->GetUDData(udd,charge);
    } else {
      charge = atomTable[i]->charge;
    }
    if (charge < -0.01 || charge > 0.01)
      output <<AtomID << " "  << charge << endl;
  }
  output << std::endl; 
  return output.str();
  
}
 
//-------------------------------------------------------------------- 
pstr CMMANManager::GetSbaseCompoundID(PCAtom p_atom) {
//--------------------------------------------------------------------
  pstr sbaseCompoundID = NULL;
  p_atom->residue->GetUDData(udd_sbaseCompoundID,sbaseCompoundID);
  return sbaseCompoundID;
}

//---------------------------------------------------------------------
pstr CMMANManager::GetSbaseCompoundID(PCResidue p_res) {
//---------------------------------------------------------------------
  pstr sbaseCompoundID = NULL;
  p_res->GetUDData(udd_sbaseCompoundID,sbaseCompoundID);
  return sbaseCompoundID; 
}

//-------------------------------------------------------------------
int CMMANManager::GetSbaseAtomOrdinal(PCAtom p_atom) {
//-------------------------------------------------------------------
  int sbaseAtomOrdinal;
  p_atom->GetUDData(udd_sbaseAtomOrdinal,sbaseAtomOrdinal);
  return  sbaseAtomOrdinal;
}

//-------------------------------------------------------------------
void CMMANManager::SetLabelMask(int i,int value) {
//-------------------------------------------------------------------
  label_mask[i] = value;
}

//---------------------------------------------------------------------
std::string CMMANManager::AtomLabel(PCAtom p_atom) {
//--------------------------------------------------------------------
  return AtomLabel(p_atom, &label_mask[0]);
}

//--------------------------------------------------------------------
std::string CMMANManager::AtomLabel(PCAtom p_atom,int mask[]) {
//--------------------------------------------------------------------
/*
There is a simpler version of this method in CMMUTManager which
only outputs the atom identifier without option of additional
atom attributes 
Mask parameters
 0  Molecule
 1  Model
 2  Chain Id
 3  Sequence number
 4  Insertion code
 5  Residue name
 6  Atom id
 7  Alternate position
 8  Element
 9  Index
 10 Atom energy type
 11 xyz
 12 Sig x,y,z 
 13 Bval
 14 Sig bval
 15 occ
 16 Sig occ
 17 charge
 18 Secondary structure 
*/
  const char secstr_text [7][12] = { " ", "beta strand", "beta bulge",
			"3-turn", "4-turn", "5-turn", "alpha helix"};
  int i;
  int masksum = 0;
  std::ostringstream label;

//cout << "mask" ;
//for ( i = 0; i < 20; i++) { cout << " " << mask[i]; }
//cout << endl;

  //if (mask[0]){
  //}

  int nmask = 0;
  for (i=0;i<=9;i++) { 
    if (mask[i]) nmask++;
  }

  if (nmask > 0 ) {
    if (mask[1]>0) {
      label << "/";
      label << p_atom->GetModelNum();   
    }
    masksum = mask[1];

    if (mask[2]>0) {
      if (masksum || strlen(p_atom->GetChainID())==0 ) label << "/";
      //if (strlen(p_atom->GetChainID())==0 ) {
      //  label << " ";
      //} else {
        label << p_atom->GetChainID();
      //}
    }

    masksum = masksum + mask[2];
    if (masksum ) label << "/";
    if (mask[3]>0) label << p_atom->GetSeqNum();

    if (mask[4] && strlen(p_atom->residue->insCode ) != 0 )   
          label << "." << p_atom->residue->insCode ;
 
    if (mask[5]>0) label << "(" << p_atom->GetResName() << ")";


    masksum = masksum + mask[3]+mask[4]+mask[5];
    if (mask[6]>0 ) {
      if (masksum) label << "/"; 
      label << TrimString(p_atom->name);
      if ( mask[7]>0 && strlen(p_atom->altLoc ) != 0  )
        label << ":" << p_atom->altLoc; 
    }
    if ( mask[8]>0 ) {
      label << "[" << TrimString(p_atom->element) << "]";
    }
//if ( mask[9]>0 ) {
// Need GetIndex in mmdb
//     label << " " << p_atom->GetIndex();
//  }
  }

  if ( mask[10]>0 ) {
    int atomOrd;
    p_atom->GetUDData(udd_atomEnergyType, atomOrd);
    label << " " << p_sbase->libAtom[atomOrd].type;
  }

  if ( mask[11]>0 ) {
label << " " << p_atom->x << " " << p_atom->y << " " << p_atom->z;  
  }

  if ( mask[13]>0 ) {
    label << " " << p_atom->tempFactor;  
  }

  if ( mask[15]>0 ) {
    if ( p_atom->occupancy < 0.999 || p_atom->occupancy > 1.001)
      label << " " << p_atom->occupancy;  
    else
      label << " 1.0";
  }

  if ( mask[17]>0 ) {
    if (p_atom->charge < -0.001 || p_atom->charge > 0.001)
      label << " " << p_atom->charge;  
    else
      label << " 0.0";
  }

  if ( mask[18]>0 ) {
    int secStr;
    secStr = p_atom->residue->SSE;
    //cout << "secStr" << secStr << endl;
    label << " " << secstr_text[secStr];
  }

  //cout << "MMAN Label" << label.str() << endl;
  return label.str();

  // Return as pstr - may not be necessary
  //int len = label.length();
  //pstr chlabel;
  //label.copy(chlabel,len,0);
  //chlabel[len] = 0;
  //return label;

}

//-----------------------------------------------------------------------
std::string CMMANManager::ListSecStructure (int mask_in[] , PCAtom pAtom) {
//-----------------------------------------------------------------------
  // Return a string with one SSE defined per line -
  // for use  in GUI listing SSEs
  int           i;
  const char  secstr_text [7][2] = { "L", "S", "U", "3", "4", "5", "H" };
  int label_mask1[20] = {0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
  int label_mask2[20] = {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  std::ostringstream output;
  int           nr;
  PPCResidue  selRes=NULL;
  int mask[7];
  int last_sec = -1;
  int last_i = -1;
  int ilres=-1;
  int ifres=0;

  // Set a mask to filter only the secstr types we want to list
  for (i=0;i<7;i++) {
    mask[i] = -1;
    if (mask_in[i]>0) mask[i]=i;
  }
  // If bulge not required then list as strand
  if (mask[2]<0 && mask[1]>0) mask[2]=1;
  //for (i=0;i<7;i++) cout << mask[i] << " ";
  //cout << endl;

  if (pAtom != NULL) {
    // Return the SSE containing this atom
    if (!isAminoacid(pAtom->residue)) return "";
    last_sec = mask[pAtom->residue->SSE];
    if (last_sec < 0) return ""; 
    last_i =  pAtom->residue->index;
    //char S[100];
    //pAtom->GetAtomID(S);
    //cout << " ListSecStructure " << S << " " << last_sec << " " << last_i << endl;
    pAtom->residue->chain->GetResidueTable(selRes,nr);
    ilres = nr-1;
    for ( i = last_i; i < nr; i++) {
      if ( mask[selRes[i]->SSE] != last_sec ) {
        ilres = i-1;
        break;
      }
    }
    for ( i = last_i; i >=0; i--) {
      if ( mask[selRes[i]->SSE] != last_sec ) {
        ifres = i+1;
        break;
      }
    }
    //cout << "ifres,ilres " << ifres << " " << ilres << " " <<last_sec << endl;
    output << secstr_text[last_sec] << " ";
    output << AtomLabel(selRes[ifres]->GetAtom(0),label_mask1)<< "-";
    output << AtomLabel(selRes[ilres]->GetAtom(0),label_mask2) << endl;
  } else { 
    // Return list of all SSEs
    //Get all residues with atoms in this atom selection
    GetResidueTable(selRes,nr);
    last_i = 0;
    last_sec = mask[selRes[0]->SSE];
    for ( i = 1; i < nr; i++) {
      //cout << i << " " << selRes[i]->seqNum << " " << selRes[i]->SSE << endl;
      if (mask[selRes[i]->SSE] != last_sec ) {
        if (last_sec>=0) {
           output << secstr_text[last_sec] << " ";
           output << AtomLabel(selRes[last_i]->GetAtom(0),label_mask1)<< "-";
           output << AtomLabel(selRes[i-1]->GetAtom(0),label_mask2) << endl;
        }
        last_i = i;
        last_sec = mask[selRes[i]->SSE];
      }
    }
    if (last_sec>=0) {
      output << secstr_text[last_sec] << " ";
      output << AtomLabel(selRes[last_i]->GetAtom(0),label_mask1)<< "-";
      output << AtomLabel(selRes[i-1]->GetAtom(0),label_mask2) << endl;
    }
  }

  return output.str();
}


//----------------------------------------------------------------------
int CMMANManager::TestBonding ( PCAtom patom1, PCAtom patom2, int max ) {
//----------------------------------------------------------------------
// For pair of atoms return 2=bonded, 3=1,3 liked, 4=1-4 linked, 
// 5=compatible Hbond types
// Only test for up to return value of max (default=5)
  PSAtomBond Bonds1;
  PSAtomBond Bonds2;
  PSAtomBond Test;
  int nBonds1,nBonds2,nTest;
  int hbtype1, hbtype2;
  int ret = 0;
 
  patom1->GetBonds(Bonds1,nBonds1);
  for ( int i=0; i< nBonds1;i++) {
    if (Bonds1[i].atom->serNum == patom2->serNum) ret = 2;
  }

  if (ret <= 0 && max>2) {
  
    patom2->GetBonds(Bonds2,nBonds2);
    for ( int i=0; i< nBonds1;i++) {
      for ( int j=0; j< nBonds2;j++) {
        if (Bonds1[i].atom->serNum == Bonds2[j].atom->serNum ) ret = 3;
      }
    }

    if (ret <= 0 && max>3) { 
  
      for ( int n=0; n< nBonds1;n++) {
        Bonds1[n].atom->GetBonds(Test,nTest);
        for  ( int i=0; i< nTest;i++) {
          for ( int j=0; j< nBonds2;j++) {
            if (Test[i].atom->serNum == Bonds2[j].atom->serNum ) ret= 4;
          }
        }
      }

      // Are they potential Hbond? 
      if (ret <= 0 && max>4) { 
        hbtype1 = GetAtomHBondType1(patom1);
        hbtype2 = GetAtomHBondType1(patom2);
        //cout << "hbtypes " << hbtype1 << " " << hbtype2 << endl;
        if (hbtype1!=HBTYPE_UNKNOWN&&hbtype2!=HBTYPE_UNKNOWN){
        if ( hbtype1 >= HBTYPE_HYDROGEN && hbtype2 >= HBTYPE_HYDROGEN ) {
          if ( (hbtype1 <= HBTYPE_BOTH && hbtype2 >= HBTYPE_BOTH) ||
	        (hbtype2 <= HBTYPE_BOTH && hbtype1 >= HBTYPE_BOTH) ) 
            ret = 5;
        }
       }
      }
    }
  }

  //cout << "ret " << ret << endl;
  return ret;
}

//---------------------------------------------------------------------
void CMMANManager::ListBonds(int selHnd,int natoms,PPCAtom selAtom) {
//---------------------------------------------------------------------
  // List the bonds
  PSAtomBond AtomBond;
  int nAtomBonds;
  int mask[20] = {0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
  int i,j;

  std::ostringstream output;

  for (i=0;i<natoms;i++) {
    //sbaseCompoundID = molHnd->GetSbaseCompoundID(selAtom[i]);
    //sbaseAtomIndex = molHnd->GetSbaseAtomIndex(selAtom[i]);
    selAtom[i]->GetBonds ( AtomBond, nAtomBonds);
    output << AtomLabel(selAtom[i],mask) << " bonded to \n";
    for (j=0;j<nAtomBonds;j++) {
      if (AtomBond[j].atom->isInSelection(selHnd)) {
        output <<  "           " << AtomLabel(AtomBond[j].atom,mask) << "\n";
      }
    }
  }
  std::cout << output.str();
}

//--------------------------------------------------------------------
int CMMANManager::RestoreData (PCMMDBManager restore_molHnd, int mode) {
//--------------------------------------------------------------------
/*
Restore data from an already loaded instance of CMMANManager
This assumes that the molecule composition is unchanged 
mode will determine the properties that are restored 
mode MMAN_COORDINATES - restore coordinates only
*/
  int nAtoms,nRestoreAtoms,i;
  PPCAtom atoms = NULL;
  PPCAtom restoreAtoms = NULL;
  GetAtomTable1(atoms,nAtoms);
  restore_molHnd->GetAtomTable1(restoreAtoms,nRestoreAtoms);
  if ( nAtoms != nRestoreAtoms) return 1;

  if ( mode == MMAN_COORDINATES) {
    for ( i = 0; i < nAtoms; i++) {
      atoms[i]->x = restoreAtoms[i]->x;
      atoms[i]->y = restoreAtoms[i]->y;
      atoms[i]->z = restoreAtoms[i]->z;
    }
  } else {
    for ( i = 0; i < nAtoms; i++) {
      atoms[i]->x = restoreAtoms[i]->x;
      atoms[i]->y = restoreAtoms[i]->y;
      atoms[i]->z = restoreAtoms[i]->z;
      atoms[i]->occupancy = restoreAtoms[i]->occupancy;
      atoms[i]->tempFactor = restoreAtoms[i]->tempFactor;
      if (restoreAtoms[i]->WhatIsSet & ASET_CoordSigma){
         atoms[i]->sigX = restoreAtoms[i]->sigX;
         atoms[i]->sigY = restoreAtoms[i]->sigY;
         atoms[i]->sigZ = restoreAtoms[i]->sigZ;
      }    
      if (restoreAtoms[i]->WhatIsSet & ASET_OccSigma){
         atoms[i]->sigOcc = restoreAtoms[i]->sigOcc;
      }    
      if (restoreAtoms[i]->WhatIsSet & ASET_tFacSigma){
         atoms[i]->sigTemp = restoreAtoms[i]->sigTemp;
      }    
      if (restoreAtoms[i]->WhatIsSet & ASET_Anis_tFac){
         atoms[i]->u11 = restoreAtoms[i]->u11;
         atoms[i]->u22 = restoreAtoms[i]->u22;
         atoms[i]->u33 = restoreAtoms[i]->u33;
         atoms[i]->u12 = restoreAtoms[i]->u12;
         atoms[i]->u23 = restoreAtoms[i]->u23;
         atoms[i]->u13 = restoreAtoms[i]->u13;
         if (restoreAtoms[i]->WhatIsSet & ASET_Anis_tFSigma){
           atoms[i]->su11 = restoreAtoms[i]->su11;
           atoms[i]->su22 = restoreAtoms[i]->su22;
           atoms[i]->su33 = restoreAtoms[i]->su33;
           atoms[i]->su12 = restoreAtoms[i]->su12;
           atoms[i]->su23 = restoreAtoms[i]->su23;
           atoms[i]->su13 = restoreAtoms[i]->su13;
         }    
      }    
    }
  }
  return 0;
}

//---------------------------------------------------------------------------
int CMMANManager::MoveFragment(int nMove, PPCAtom moveAtoms, Cartesian dxyz){
//--------------------------------------------------------------------------- 
  /*
  Move a defined set of atoms
  */
  int i;
  for (i=0;i<nMove;i++) {
    moveAtoms[i]->x =  moveAtoms[i]->x + dxyz.get_x();
    moveAtoms[i]->y =  moveAtoms[i]->y + dxyz.get_y();
    moveAtoms[i]->z =  moveAtoms[i]->z + dxyz.get_z();
  }
  return 0;  
}


//-----------------------------------------------------------------------
int CMMANManager::LoadUDDData( const int property ) { 
//-----------------------------------------------------------------------
  int selHnd;
  PPCAtom atomTable = NULL;
  PPCResidue resTable = NULL;
  int nAtoms,nRes;
  int udd;

  //cout << "LoadUDDData " << property << endl;

  switch (property) {
  case PROPERTY_SEC:
    udd = GetUDDHandle ( UDR_RESIDUE,"tmp_res_int" );
    if (udd <= 0 ) {
      udd = -1;
      udd = RegisterUDInteger ( UDR_RESIDUE,"tmp_res_int" );
      if (udd <= 0 ) return udd;
    }
    selHnd = NewSelection();
    Select(selHnd,STYPE_RESIDUE,0,"*",ANY_RES,
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
    GetSelIndex(selHnd,resTable,nRes);
    for (int j=0;j<nRes;j++) {
       resTable[j]->PutUDData(udd,resTable[j]->SSE);
    }
    DeleteSelection(selHnd);
    //cout << "Property udd " << udd << endl;
    return udd;
    break;

  case PROPERTY_ATOM_SAS:
    udd = GetUDDHandle( UDR_ATOM,"atom_sas");
    return udd;
    break;
  case PROPERTY_RES_SAS:
    udd = GetUDDHandle( UDR_RESIDUE,"residue_sas");
    return udd;
    break;
  case PROPERTY_ATOM_CONTACT:
    udd = GetUDDHandle( UDR_ATOM,"atom_contact");
    return udd;
    break;
  case PROPERTY_RES_CONTACT:
    udd = GetUDDHandle( UDR_RESIDUE,"residue_contact");
    return udd;
    break; 
  }


  //Copy some data such as b or x into a UDD
  udd = GetUDDHandle ( UDR_ATOM,"tmp_atom_real" );
  if (udd <= 0 ) {
    udd = -1;
    udd = RegisterUDReal ( UDR_ATOM,"tmp_atom_real" );
    if (udd <= 0 ) return udd;
  }
  //cout << "Property udd " << udd << endl;

  // Loop over all atoms to do a copy 
  selHnd = NewSelection();
  Select(selHnd,STYPE_ATOM,0,"*",ANY_RES,
         "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  GetSelIndex(selHnd,atomTable,nAtoms);
 
  switch (property) {
  case PROPERTY_B:
    for (int j=0;j<nAtoms;j++) {
      atomTable[j]->PutUDData(udd,atomTable[j]->tempFactor);
    }
    break;
  case PROPERTY_OCC:
    for (int j=0;j<nAtoms;j++) {
      atomTable[j]->PutUDData(udd,atomTable[j]->occupancy);
    }
    break;
//case PROPERTY_CHARGE:
//  for (int j=0;j<nAtoms;j++) {
//    atomTable[j]->PutUDData(udd,atomTable[j]->charge);
//  }
//  break;
  case PROPERTY_X:
    for (int j=0;j<nAtoms;j++) {
      atomTable[j]->PutUDData(udd,atomTable[j]->x);
    }
    break;
  case PROPERTY_Y:
    for (int j=0;j<nAtoms;j++) {
      atomTable[j]->PutUDData(udd,atomTable[j]->y);
    }
    break;
  case PROPERTY_Z:
    for (int j=0;j<nAtoms;j++) {
      atomTable[j]->PutUDData(udd,atomTable[j]->z);
    }
    break;
  case PROPERTY_SERIAL:
    for (int j=0;j<nAtoms;j++) {
      atomTable[j]->PutUDData(udd,float(atomTable[j]->serNum));
    }
    break;
  }
  DeleteSelection(selHnd);
  return udd; 
}



CMMANManager* GetMMANManager(PCAtom pAtom) {
  return  (PCMMANManager)(PCMMDBFile)pAtom->GetCoordHierarchy();
}

/*
//----------------------------------------------------------------------
std::string CMMANManager::GetSymOpTitle(int nsym,int i,int j,int k) {
//---------------------------------------------------------------------
  std::ostringstream output;
  CSymOp SymOp;      
  mat44 TMatrix;
  char symOpTitle[100]; 

  if (nsym < 0 || nsym >= GetNumberOfSymOps()) return output.str();

  SymOp.SetSymOp ( GetSymOp(nsym) );

  TMatrix[0][3] += i;
  TMatrix[1][3] += j;
  TMatrix[2][3] += k;
  SymOp.CompileOpTitle ( symOpTitle,TMatrix,False ); 

  output << symOpTitle;
  return output.str();
}
*/

//----------------------------------------------------------------------
std::string CMMANManager::GetSymOpTitle(int nsym,int i,int j,int k) {
//---------------------------------------------------------------------
  std::ostringstream output;
  pstr symOpTitle;
  CSymOps SymOps;

  SymOps.SetGroup(GetSpaceGroup());
  if (nsym < 0 || nsym >= SymOps.GetNofSymOps()) return output.str();

  symOpTitle = SymOps.GetSymOp(nsym);

  output << symOpTitle;
  return output.str();
}


//------------------------------------------------------------------------
int CMMANManager::GenerateTransformedModel(int model,realtype *vmat) {
//------------------------------------------------------------------------
  // Apply transformation
  int new_model;
  PCModel p_model;
  mat44 TMatrix;
  mat44 TMatrixT;

  new_model =  CopyModel(model);
  p_model = GetModel(new_model);
  if (!p_model) {
    if(GetNumberOfModels()>0){
        p_model = GetModel(1);
        if (!p_model) {
	  return -1;
        }
      }
  }
  TMatrixT[0][0] = 1;
  TMatrixT[1][1] = 1;
  TMatrixT[2][2] = 1;
  TMatrixT[3][3] = 1;
  TMatrixT[0][1] = 0;
  TMatrixT[0][2] = 0;
  TMatrixT[1][0] = 0;
  TMatrixT[1][2] = 0;
  TMatrixT[2][0] = 0;
  TMatrixT[2][1] = 0;
  TMatrix[0][3] = 0;
  TMatrix[1][3] = 0;
  TMatrix[2][3] = 0;
  TMatrix[3][0] = 0;
  TMatrix[3][1] = 0;
  TMatrix[3][2] = 0;
  TMatrix[3][3] = 1;
  int kk = 0;
  for (int ii = 0;ii<4;ii++) {
    for (int jj = 0;jj<4;jj++) {
      if(ii<3&&jj<3)
      TMatrix[jj][ii] = vmat[kk];
      else
      TMatrixT[jj][ii] = -vmat[kk];
      kk++;
    }
  }
  /*
  for (int ii = 0;ii<4;ii++) {
    for (int jj = 0;jj<4;jj++) {
      std::cout << TMatrix[ii][jj] << " ";
    } std::cout << "\n";
  } std::cout << "\n";
  for (int ii = 0;ii<4;ii++) {
    for (int jj = 0;jj<4;jj++) {
      std::cout << TMatrixT[ii][jj] << " ";
    } std::cout << "\n";
  } std::cout << "\n";
  */
  p_model->ApplyTransform (TMatrixT);
  p_model->ApplyTransform (TMatrix);
  return new_model;
}

//------------------------------------------------------------------------
int CMMANManager::ApplyTransformtoModel(int model,realtype *vmat,Boolean undo) {
//------------------------------------------------------------------------
  mat44 TMatrix,invTMatrix;
  PCModel p_model = GetModel(model);
  //std::cout << "model num: " << model << "\n";
  //std::cout << "model: " << p_model << "\n";
  //std::cout << "number of models: " << GetNumberOfModels() << "\n";
  if (!p_model) {
    if(GetNumberOfModels()>0){
        p_model = GetModel(1);
        if (!p_model) {
	  return 1;
        }
      }
  }
  int kk = 0;
  for (int ii = 0;ii<4;ii++) {
    for (int jj = 0;jj<4;jj++) {
      TMatrix[jj][ii] = vmat[kk];
      kk++;
    }
  }
  
  if (undo) {
    Mat4Inverse (TMatrix,invTMatrix);
    p_model->ApplyTransform (invTMatrix);
  } else
    p_model->ApplyTransform (TMatrix);
  return 0;
}

//------------------------------------------------------------------------
int CMMANManager::GetLibTMatrix(mat44 &TMatrix,int nsym,int i,int j,int k) {
//------------------------------------------------------------------------
  CSymOps SymOps;
  //mat44 transmat,rotmat;
  int rv = 0;
  //cout << "GetLibTMatrix " << nsym << " " << i << " " << j << " " << k << endl;
  
  /*
  This uses symops from syminfo but does not seem to get the translation
  component correct
  SymOps.SetGroup(GetSpaceGroup());
  cout << "GetLibTMatrix SymOps " << SymOps.GetSymOp(nsym) <<endl;
  rv = SymOps.GetTMatrix(rotmat,nsym);

  GetTMatrix(transmat,0,i,j,k);  
  cout << "transmat" << endl;
  for (int ii=0;ii<4;ii++) 
    cout << transmat[0][ii] << " " <<  transmat[1][ii] << " "  << transmat[2][ii] << " " << transmat[3][ii] << endl;
  Mat4Mult(TMatrix,transmat,rotmat);
  */
  
  
  GetTMatrix(TMatrix,nsym,i,j,k);
  /* cout << "Tmatrix" << endl;
  for (int ii=0;ii<4;ii++) 
    cout << TMatrix[0][ii] << " " <<  TMatrix[1][ii] << " "  << TMatrix[2][ii] << " " << TMatrix[3][ii] << endl;
  */
  rv = 0;
  

  return rv;
}

//------------------------------------------------------------------------
int CMMANManager::GenerateSymmetryModel(int model,int nsym,int i,int j,int k) {
//------------------------------------------------------------------------
  // Apply transformation
  int new_model;
  PCModel p_model;
  mat44 TMatrix;

  new_model =  CopyModel(model);
  p_model = GetModel(new_model);
  if (!p_model) {
	  return -1;
  }
  GetLibTMatrix(TMatrix,nsym,i,j,k);
  p_model->ApplyTransform (TMatrix);
  return new_model;
}


//------------------------------------------------------------------------
int CMMANManager::ApplySymmetrytoModel(int model,int nsym,int i,int j,int k,Boolean undo) {
//------------------------------------------------------------------------
  // Apply transformation
  mat44 TMatrix, invTMatrix;
  PCModel p_model;

  p_model = GetModel(model);
  if (p_model == NULL) return 1;
  GetLibTMatrix(TMatrix,nsym,i,j,k);
  if (undo) {
    GetLibTMatrix(invTMatrix,nsym,i,j,k);
    Mat4Inverse (invTMatrix,TMatrix);
  } else 
     GetLibTMatrix(TMatrix,nsym,i,j,k);
     
  p_model->ApplyTransform (TMatrix);
  return 0;
}

//-----------------------------------------------------------------------
int CMMANManager::IfSymmetryNeighbours(int selHnd, int model, int nsym, 
                 int i, int j, int k, double dist ) {
//-----------------------------------------------------------------------
  mat44 TMatrix;
  PSContact contact = NULL; 
  int ncontacts; 
  int   maxlen = 0;
  PPCAtom selAtoms=NULL, modelAtoms=NULL;
  int modelSelHnd, nSelAtoms, nModelAtoms;

  if ( dist < 0.1 || dist > 100.0 ) dist = 10.0;

  if (selHnd>0) {
    GetSelIndex(selHnd, selAtoms,nSelAtoms );
    if (nSelAtoms <= 0) return -1;
  } 

  modelSelHnd = NewSelection();
  Select(modelSelHnd,STYPE_ATOM,model,"*",ANY_RES,
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  GetSelIndex(modelSelHnd,modelAtoms,nModelAtoms);
  if ( nModelAtoms <= 0 ) {
    DeleteSelection(modelSelHnd);
    return -2;
  }

  // If selection is not sepecified then take contacts to all atoms
  // in model
  if ( selHnd <= 0 )  GetSelIndex(modelSelHnd,selAtoms,nSelAtoms);  

  GetLibTMatrix(TMatrix,nsym,i,j,k);
  
  SeekContacts ( selAtoms, nSelAtoms, modelAtoms, nModelAtoms, 0.0, dist, 0,
			 contact,ncontacts,maxlen,&TMatrix);
  
  if (contact)  delete [] contact;
  if (modelSelHnd)  DeleteSelection(modelSelHnd);
  return ncontacts;
}

//------------------------------------------------------------------------
int CMMANManager::CopyModel(int model) {
//------------------------------------------------------------------------
  PCModel new_model;
  
  new_model = new CModel();
  new_model->Copy(GetModel(model));
  AddModel(new_model);
  PDBCleanup(PDBCLEAN_SERIAL);
  //cout << "CopyModel serNum " << new_model->GetSerNum() <<endl;

  // Copy atom types and bonds
  int selHnd1,selHnd2;
  PPCAtom p_atom1,p_atom2;
  PPCResidue p_res1,p_res2;
  int ia2,nat1,nat2,nb,eType,RC,udd;
  char S[100];
  PSAtomBond bonds;

  selHnd1 = NewSelection();
  selHnd2 = NewSelection();

  Select(selHnd1,STYPE_RESIDUE,model,"*",ANY_RES, 
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  Select(selHnd2,STYPE_RESIDUE,new_model->GetSerNum(),"*",ANY_RES, 
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  GetSelIndex(selHnd1, p_res1, nat1);
  GetSelIndex(selHnd2, p_res2, nat2);
  // Copy residue UDD
  for ( int ia=0;ia<nat1;ia++) {
    p_res1[ia]->GetUDData(udd_sbaseCompoundID, eType);
    RC = p_res2[ia]->PutUDData(udd_sbaseCompoundID, eType);
  }  

  // Copy atom UDD
  Select(selHnd1,STYPE_ATOM,model,"*",ANY_RES, 
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  Select(selHnd2,STYPE_ATOM,new_model->GetSerNum(),"*",ANY_RES, 
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  GetSelIndex(selHnd1, p_atom1, nat1);
  GetSelIndex(selHnd2, p_atom2, nat2);
  //cout << "nat1,nat2 " << nat1 << " " << nat2 << endl;
  
  for ( int ia=0;ia<nat1;ia++) {
    p_atom1[ia]->GetUDData(udd_atomEnergyType, eType);
    RC = p_atom2[ia]->PutUDData(udd_atomEnergyType, eType);
    p_atom1[ia]->GetUDData(udd_sbaseAtomOrdinal, eType);
    RC = p_atom2[ia]->PutUDData(udd_sbaseAtomOrdinal, eType);
  }

  udd = GetUDDHandle ( UDR_ATOM,"tmp_atom_int" );
  if (udd <= 0 ) {
    udd = -1;
    udd = RegisterUDInteger( UDR_ATOM,"tmp_atom_int" );
  }

  // Put incremental index in udd 
  for ( int ia=0;ia<nat1;ia++)
    RC = p_atom1[ia]->PutUDData(udd,ia);

  // Copy bonds 

  for ( int ia=0;ia<nat1;ia++) {
    //p_atom1[ia]->GetAtomID(S);
    p_atom1[ia]->GetBonds(bonds,nb);
    if (nb>0) {
      bonds[0].atom->GetUDData(udd,ia2);
      //cout << "ia2 " << ia2 << " " ;
      p_atom2[ia]->AddBond(p_atom2[ia2], bonds[0].order,nb);
      if (nb>1) {
        for ( int ib=1;ib<nb;ib++) {
          bonds[ib].atom->GetUDData(udd,ia2);
          //cout << ia2 << " " ;
          p_atom2[ia]->AddBond(p_atom2[ia2], bonds[ib].order,nb);
        }
      }
    }
    //cout << endl;
    if ( ia < 5 ) {
      p_atom2[ia]->GetAtomID(S);
      p_atom2[ia]->GetUDData(udd_atomEnergyType, eType);
      p_atom2[ia]->GetBonds(bonds,nb);
      //cout << S << " " << eType << " " << nb << endl;
      for (int ib =0; ib<nb;ib++) {
        bonds[ib].atom->GetAtomID(S);
        //cout << S << " ";
      }
      //cout << endl;
    }
  }
  //cout << "CopyModel before DeleteSelection" << endl;
  DeleteSelection(selHnd1);
  DeleteSelection(selHnd2);
  //cout << "CopyModel after DeleteSelection" << endl;

  return GetNumberOfModels();
  
}

//--------------------------------------------------------------------------------
int CMMANManager::SelectChainTermini( void ) {
//--------------------------------------------------------------------------------

  /*
  Return a selection handle for set of CA atoms of N/C chain termini residues
  This is used for labelling the termini
  */
  PPCChain chainTable = NULL;
  PPCResidue resTable = NULL;
  int nChains,nRes,iser1,iser2;

  int tmpHnd = NewSelection();

  for (int nm=1;nm<=GetNumberOfModels();nm++) {
    chainTable=NULL;
    GetChainTable(nm,chainTable,nChains);
    //cout << "nChains " << nChains << endl;
    for (int nc=0;nc<nChains;nc++) {
      resTable = NULL;
      chainTable[nc]->GetResidueTable(resTable,nRes);
      for (int nr=0;nr<nRes;nr++) {
        if ( isAminoacid(resTable[nr]) && resTable[nr]->GetAtom("CA","C","*") ) {
          iser1 = iser2 = resTable[nr]->GetAtom("CA","C","*")->serNum;
          SelectAtoms(tmpHnd,iser1,iser2,SKEY_OR);
          break;
        }
      }
      for (int nr=nRes-1;nr>=1;nr--) {
        if ( isAminoacid(resTable[nr]) && resTable[nr]->GetAtom("CA","C","*") ) {
          iser1 = iser2 = resTable[nr]->GetAtom("CA","C","*")->serNum;
          SelectAtoms(tmpHnd,iser1,iser2,SKEY_OR);
          break;
        }
      }
    }
  }
  return tmpHnd;
}

//--------------------------------------------------------------------------------
int CMMANManager::SelectSSETermini( int selHnd ) {
//--------------------------------------------------------------------------------

  /*
  Return a selection handle for set of CA atoms of SSE termini residues
  This is used for labelling the termini
  */
  PPCResidue resTable = NULL;
  PCAtom pCA,pCA0;
  int nRes,iser1,iser2;
  int nr, nr0, current_SSE;
  int resSelHnd;
  bool noSelHnd = false;
  int tmpHnd = NewSelection();
  
  if (selHnd>=0) {
    resSelHnd = NewSelection();
    Select(resSelHnd,STYPE_RESIDUE,selHnd,SKEY_NEW);
    GetSelIndex(resSelHnd,resTable,nRes);
  } else {
    noSelHnd = true;
    GetResidueTable(resTable,nRes);
  }
  nr = 0;
  while (nr<nRes-1) {
    if ( !isAminoacid(resTable[nr]) || !resTable[nr]->GetAtom("CA","C","*") ) {
      nr++;
    } else {
      pCA = resTable[nr]->GetAtom("CA","C","*");
      if ( (noSelHnd || pCA->isInSelection(selHnd)) &&    
             ( resTable[nr]->SSE ==  SSE_Helix || 
               resTable[nr]->SSE ==  SSE_Bulge ||
               resTable[nr]->SSE == SSE_Strand ) ) {
         // Select first residue of SSE
         iser1 = iser2 = pCA->serNum;       
         SelectAtoms(tmpHnd,iser1,iser2,SKEY_OR);
         if (resTable[nr]->SSE ==  SSE_Helix) {
           current_SSE = SSE_Helix;
         } else {
           current_SSE = SSE_Strand;
         }
         // Look on from the next residue for the end of the SSE
         nr0 = nr + 1;
         while (   nr0 < nRes &&
              (resTable[nr0]->SSE == current_SSE || 
                (resTable[nr0]->SSE ==SSE_Bulge && current_SSE == SSE_Strand ))) {
            nr0++;
         }
         // Select  nr0-1 which is last residue of SSE
         pCA0 =  resTable[nr0-1]->GetAtom("CA","C","*");
         if ( pCA0 && ( noSelHnd || pCA0->isInSelection(selHnd)) ) {
           iser1 = iser2 = resTable[nr0-1]->GetAtom("CA","C","*")->serNum;
           SelectAtoms(tmpHnd,iser1,iser2,SKEY_OR);
         }
         // Set nr to the number of the first residue after this SSE
         nr = nr0;
      } else {
         nr++;
      }
    }
  }
  return tmpHnd;
}

//-----------------------------------------------------------------
bool CMMANManager::isAminoacid (PCResidue pres) {
//-----------------------------------------------------------------
  pstr sbaseCompoundID=NULL;
  PCSBStructure  pSbRes;
  // Is it a user defined 'custom' restype?
  if(customResTypes.size()>0){
  std::map<std::string,int>::iterator pos = customResTypes.find(pres->name);
  if (pos!=customResTypes.end()) {
    if ( pos->second == RESTYPE_PEPTIDE || pos->second == RESTYPE_LPEPTIDE || pos->second == RESTYPE_DPEPTIDE )
      return true;
    else
      return false;
  }
  }

  if(!p_sbase) return pres->isAminoacid();
  if (strlen(p_sbase->monomers_dir)<1) {
    return pres->isAminoacid();
  } else {
    if ( pres->chain->GetModel()->GetSerNum()>1) {
      PCChain pc = NULL;
      if (GetModel(1) != NULL) pc = GetModel(1)->GetChain(pres->chain->GetChainID());
      if (!pc)  return  pres->isAminoacid();
      PCResidue pres1= pc->GetResidue(pres->seqNum,pres->insCode);
      if (!pres1) return  pres->isAminoacid();
      pres1->GetUDData(udd_sbaseCompoundID,sbaseCompoundID);
    } else {
      pres->GetUDData(udd_sbaseCompoundID,sbaseCompoundID);
    }


    if(!sbaseCompoundID) return false;

    if ( strlen(sbaseCompoundID) < 1 ) {
     delete [] sbaseCompoundID;
     return false;
    }
    pSbRes = p_sbase->GetStructure(sbaseCompoundID,monlib);
    delete [] sbaseCompoundID;
    if (!pSbRes) return false;
    std::string formula = pSbRes->Formula;
    //cout << "isAminoacid " << sbaseCompoundID << " " << formula << endl;
    if ( strcmp(pSbRes->Formula, "L-peptide")==0 ||
       strcmp(pSbRes->Formula, "D-peptide")==0 ||
       strcmp(pSbRes->Formula, "peptide")==0 ) {
      return true;
    } else {
      return false;
    }
  }
}

//------------------------------------------------------------------------
int CMMANManager::GetRestypeCode ( PCResidue pres ) {
//------------------------------------------------------------------------
  pstr sbaseCompoundID=NULL;
  PCSBStructure  pSbRes;
  // Is it a user defined 'custom' restype?
  if(customResTypes.size()>0){
  std::map<std::string,int>::iterator pos = customResTypes.find(pres->name);
  if (pos!=customResTypes.end()) return pos->second;
  }

  if (p_sbase&&(strlen(p_sbase->monomers_dir)>1) ) {
    pres->GetUDData(udd_sbaseCompoundID,sbaseCompoundID);

    //cout << "sbaseCompoundID " << sbaseCompoundID << endl;
    if ( sbaseCompoundID && strlen(sbaseCompoundID)>1 ) {
      pSbRes = p_sbase->GetStructure(sbaseCompoundID,monlib);
      delete [] sbaseCompoundID;
      if (!pSbRes) {
         map<std::string,std::string>::iterator p = p_sbase->synonyms.find(pres->name);
         if (p!=p_sbase->synonyms.end()) pSbRes = p_sbase->GetStructure(p->second.c_str(),monlib);
         //cout << "testing synonyms " << pres->name << " " << pSbRes << endl;
      }
      if (pSbRes) {
        //cout << "GetRestypeCode " << pres->name << " " << pSbRes->Formula << endl;
        if ( strcmp(pSbRes->Formula, "L-peptide")==0 ||
           strcmp(pSbRes->Formula, "D-peptide")==0 ||
           strcmp(pSbRes->Formula, "peptide")==0 ) 
         return RESTYPE_PEPTIDE;
        else if ( strcmp(pSbRes->Formula, "DNA")==0 ||
              strcmp(pSbRes->Formula, "RNA")==0 ||
              strcmp(pSbRes->Formula, "DNA/RNA")==0 )
          return RESTYPE_NUCL;
        else if (strcmp(pSbRes->Formula, "solvent")==0)
          return RESTYPE_SOLVENT;
        else if ( strcmp(pSbRes->Formula, "saccharide")==0 ||
           strcmp(pSbRes->Formula, "pyranose")==0 ||
           strcmp(pSbRes->Formula, "D-saccharid")==0 ||
           strcmp(pSbRes->Formula, "L-saccharid")==0 ) 
          return RESTYPE_SACH;
        else {
          for(int iElementMetals=0;iElementMetals<nElementMetals;iElementMetals++){
            if(pres->GetNumberOfAtoms()==1&&pres->GetAtom(ElementMetal[iElementMetals])&&(!strncmp(ElementMetal[iElementMetals],pres->GetAtom(ElementMetal[iElementMetals])->name,2)||!strncmp(ElementMetal[iElementMetals],pres->GetAtom(ElementMetal[iElementMetals])->element,2))){
              //std::cout << "Metal: " << ElementMetal[iElementMetals] << "\n";
              return RESTYPE_METAL;
            }
          }
          return RESTYPE_NONPOLY;
        }
      }
    }
  }
  //cout << "GetRestypeCode defaulting" << endl;
  
  if ( pres->isAminoacid())
    return RESTYPE_PEPTIDE;
  else if (pres->isNucleotide())
    return RESTYPE_NUCL;
  else if (pres->isSolvent() )
    return RESTYPE_SOLVENT;
  else if (pres->isSugar())
    return RESTYPE_SACH;

  if (pres->GetAtom("CA","*","*") &&  
      pres->GetAtom("N","*","*") &&  
      pres->GetAtom("C","*","*")) {
     return RESTYPE_PEPTIDE;
  } else if ((pres->GetAtom("C3'","*","*") || pres->GetAtom("C3*","*","*")) &&
             (pres->GetAtom("C5'","*","*") || pres->GetAtom("C5*","*","*")) &&
	      pres->GetAtom("P","*","*") ) {
    return RESTYPE_NUCL;
  }
  
  for(int iElementMetals=0;iElementMetals<nElementMetals;iElementMetals++){
    if(pres->GetNumberOfAtoms()==1&&pres->GetAtom(ElementMetal[iElementMetals])&&(!strncmp(ElementMetal[iElementMetals],pres->GetAtom(ElementMetal[iElementMetals])->name,2)||!strncmp(ElementMetal[iElementMetals],pres->GetAtom(ElementMetal[iElementMetals])->element,2))){
      //std::cout << "Metal: " << ElementMetal[iElementMetals] << "\n";
      return RESTYPE_METAL;
    }
  }
  return RESTYPE_NONPOLY;

}

int CMMANManager::SetCustomRestype ( const std::string &resname , const int &restype , Boolean clear ) {
  //cout << "SetCustomRestype " << resname << " " << restype << endl;
  if (clear) customResTypes.clear();
  customResTypes[resname]=restype;
  return customResTypes.size();
}
int CMMANManager::SetCustomResSynonym ( const std::string &resname , const std::string &alias , Boolean clear ) {
  if (clear) customResSynonym.clear();
  customResSynonym[resname]=alias;
  return customResSynonym.size();
}

int ExcludeOverlappedAtoms ( CMMDBManager* molHnd, const int selHnd ,
                                const realtype cutoff, int theModel ) {
  // Remove atoms closer that cutoff distance from a selection
  //std::cout << "ExcludeOverlappedAtoms\n";
  int tmpHnd1,tmpHnd2;
  PPCAtom selAtoms1=NULL;
  PPCAtom selAtoms2=NULL;
  int nat1,nat2;
  PSContact contacts = NULL;
  int ncontacts;
  mat44 * TMatrix=0;
  int iser1,iser2;
  int ndel = 0;


  tmpHnd1 = molHnd->NewSelection();
  tmpHnd2 = molHnd->NewSelection();
  int nmodels = molHnd->GetNumberOfModels();
  for(int imodel=1;imodel<=nmodels;imodel++){
    if(theModel>0&&imodel!=theModel) continue;

    // Copy the selection
    molHnd->Select(tmpHnd1,STYPE_ATOM,imodel,"*", ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_NEW);
    molHnd->Select(tmpHnd1,STYPE_ATOM,selHnd,SKEY_AND);
    molHnd->Select(tmpHnd2,STYPE_ATOM,imodel,"*", ANY_RES,"*",ANY_RES,"*","*","*","*","*",SKEY_NEW);
    molHnd->Select(tmpHnd2,STYPE_ATOM,selHnd,SKEY_AND);
    selAtoms1 = NULL;
    selAtoms2 = NULL;
    
    molHnd->GetSelIndex(tmpHnd1,selAtoms1,nat1);  
    molHnd->GetSelIndex(tmpHnd2,selAtoms2,nat2);
    //std::cout << "model: " << imodel << ", nat1: " << nat1 << ", nat2: " << nat2 << "\n"; std::cout.flush();
    if(nat1>0&&nat2>0) {
   

      contacts = NULL;
      molHnd->SeekContacts(selAtoms1,nat1,selAtoms2,nat2,0.0,cutoff,0,
		           contacts,ncontacts,0,TMatrix, 0 , 0);

      //std::cout << "ncontacts: " << ncontacts << "\n"; std::cout.flush();
      if (contacts && ncontacts>0) {
      for (int i=0;i<ncontacts;i++) {
        iser1 = selAtoms1[contacts[i].id1]->serNum;
        iser2 = selAtoms2[contacts[i].id2]->serNum;
        //cout << i << "(" << ncontacts << ")" << ", iser1,iser2 " << iser1 << " " << iser2 << endl;
        // Arbitrarilly exclude the atom with higher serial number
        // beware each contact listed twice
        if ( iser1 > iser2 ) {
          iser2 = iser1;
          //SelectAtoms(selHnd,iser1,iser2,SKEY_CLR);
          molHnd->UnselectAtoms(selHnd,iser1,iser2);
          ndel++;
        }
      } }

    }

  }
  molHnd->DeleteSelection(tmpHnd1);
  molHnd->DeleteSelection(tmpHnd2);

  return ndel;
  
}

//------------------------------------------------------------------------
int CMMANManager::ExcludeOverlappedAtoms ( const int selHnd ,  \
                                const realtype cutoff, int theModel ) {
//------------------------------------------------------------------------
  return ::ExcludeOverlappedAtoms(this,selHnd,cutoff,theModel);
}

int CMMANManager::SetTransform ( const matrix tMat, const bool reset) {
   mat44 TMatrix;
   for (int i = 0;i<4;i++) {
     for (int j = 0;j<4;j++) {
       TMatrix[j][i]=tMat(j,i);
       //cout << TMatrix[j][i] << "  " ; 
     }
   }
   //cout << endl;
  return SetTransform ( TMatrix , reset );
}

int CMMANManager::SetTransform ( const std::vector<float>& transf , const std::vector<float>& transl , const bool reset) {

  mat44 TMatrix;

  if ( transf.size() != 9 || transl.size() != 3) {
    if (reset) {
      UnSetTransform();
    } else {
      return 1;
    }
  } else {

    int n=0;
    for (int i = 0;i<3;i++) {
     TMatrix[3][i]=0.0;
     for (int j = 0;j<3;j++) {
        TMatrix[i][j] = transf[n];      
       n++;
     }
   } 
    for (int i = 0;i<3;i++) {
      TMatrix[i][3] = transl[i];
   }
    TMatrix[3][3]=1.0;
  }

  if (transform_com_set ) {
    for (int i=0;i<3;i++) {
      transform_com[i]=  transform_com[i] + transl[i];
    }
  }

  return SetTransform ( TMatrix , reset );
}


void CMMANManager::UnSetTransform(bool apply_inverse) {
  //cout << "UnSetTransform isTransformed " << isTransformed << " apply inverse " << apply_inverse << endl;
  
  if (isTransformed && apply_inverse) {
    mat44 InvMatrix;
    Mat4Inverse(current_transform,InvMatrix);
    ApplyTransform(InvMatrix);
  }
  Mat4Init(current_transform);
  isTransformed = false;
  transform_com_set = false;
  //cout << "unset transform" << endl;
}


int CMMANManager::SetTransform ( const realtype rot ,const std::vector<float>& axis , const int selHndin) {
  mat44 Tmat;
  int selHnd;

  if ((!transform_com_set) || selHndin>0 ) {
    PPCAtom atomTable;
    int nAtoms;
    realtype xx=0.0,yy=0.0,zz=0.0;
    if (selHndin > 0 ) {
      selHnd = selHndin;
    } else {
      selHnd = NewSelection();
      Select(selHnd,STYPE_ATOM,0,"*",ANY_RES,
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
    }
    GetSelIndex(selHnd,atomTable,nAtoms);
    //cout <<" selHndin " << selHndin << " nAtoms " << nAtoms << endl;
    for (int i =0;i<nAtoms;i++) {
      xx = xx+atomTable[i]->x;
      yy = yy+atomTable[i]->y;
      zz = zz+atomTable[i]->z;
    }
    transform_com[0] = xx/nAtoms;
    transform_com[1] = yy/nAtoms;
    transform_com[2] = zz/nAtoms;
  
    if (selHnd<=0) DeleteSelection(selHnd);
    transform_com_set = true;
  }

  if (rot < 0.0001 && rot > -0.0001) return 0;

  //cout << "COM " << transform_com[0]<< " " << transform_com[1]<< " " <<transform_com[2] << endl;

  GetVecTMatrix (Tmat , rot*PI/180,  axis[0], axis[1], axis[2], 
               transform_com[0],transform_com[1],transform_com[2] ); 

  SetTransform ( Tmat , 0 );
  return 0;
}



int CMMANManager::SetTransform (const std::vector<float> &TMatrix , const bool reset) {
  mat44 newMatrix;
  int kk=0;
  for(int k=0;k<4;k++){
    for(int i=0;i<4;i++){
       newMatrix[k][i] = TMatrix[kk];
       kk++;
    }
  }
  SetTransform(newMatrix,reset);
}

int CMMANManager::SetTransform (mat44 &TMatrix , const bool reset) {

  if (reset) UnSetTransform();

  mat44 newMatrix;
  
  ApplyTransform(TMatrix);
  
  isTransformed = true;
  for(int k=0;k<4;k++){
    for(int s=0;s<4;s++){
      newMatrix[k][s]=0.0;
      for(int i=0;i<4;i++){
        newMatrix[k][s] += TMatrix[k][i] * current_transform[i][s];
      }
    }
  }
  for(int k=0;k<4;k++){
    for(int i=0;i<4;i++){
      current_transform[k][i]=newMatrix[k][i];
    } 
  } 
 
  /*
  cout << "\nSetTransform:\n\n";
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) cout << current_transform[i][j] << " ";
    cout << endl;
  }
  */
  
  return 0;
}

int CMMANManager::ReApplyTransform( const bool reset) {

  if (reset) UnSetTransform();  
  ApplyTransform(current_transform);

  return 0;
}

//------------------------------------------------------------------
std::vector<float> CMMANManager::GetTransform() {
//------------------------------------------------------------------
  // Return transform matrix as standard vector - for Python interface
  std::vector<float> mtrx;
  mtrx = std::vector<float> (16);
  int k=0;
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      mtrx[k]=current_transform[i][j];
      k++; 
    }
  }
  return mtrx;
}

//----------------------------------------------------------------------------
std::string CMMANManager::GetTransformString() {
//----------------------------------------------------------------------------
  std::ostringstream output;
  output.setf(ios::fixed);
  output.setf(ios::showpoint);
  output.setf(ios::left,ios::adjustfield);
  output.precision(3);

  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      output << std::setw(17) <<  current_transform[i][j];
    }
    output << endl;
  }
  
  //output << endl << "translation" << endl;
  //for (int j=0;j<3;j++) {
  //  output  << std::setw(20) <<  current_transform[j][3];
  //}
  output << endl;
  return output.str();
}

//----------------------------------------------------------------------------
int CMMANManager::TransformToSuperposeAtoms (  PPCAtom A1, int nA, PPCAtom A2) {
//-----------------------------------------------------------------------------
  /*
  Atom set A1 should be in this model and are the 'moving' atoms
  Atom set A2 can be in another model and are the 'fixed' atoms
  This method avoids need to pass a mat44 to Python
  */
  mat44 TMatrix;
  //cout << "TransformToSuperposeAtoms nA " << nA << endl;
  int rv = SuperposeAtoms ( TMatrix, A1, nA, A2 );
  //cout << "TransformToSuperposeAtoms rv " << rv <<  endl;
  if (rv == SPOSEAT_Ok) {
    SetTransform (TMatrix,0);
  }
  return rv;
}
//----------------------------------------------------------------------------
double CMMANManager::DeltaResidueOrientation (PCResidue pRes,PCResidue pResFx) {
//-----------------------------------------------------------------------------
  // Create selection of atoms that will be superposed
  PCAtom mvAtoms[5];
  PCAtom fxAtoms[5];
  int nat=0;
  char *names[] = { "CA" , "N", "C", "CB" };

  for (int n=0;n<4;n++ ){
    mvAtoms[nat]=pRes->GetAtom(names[n],"*","*");
    fxAtoms[nat]=pResFx->GetAtom(names[n],"*","*");
    if (mvAtoms[nat] && fxAtoms[nat] ) nat++;
  }
  if (nat<0) return -10.0;
  
  mat44 TMatrix;
  int rv = SuperposeAtoms ( TMatrix, mvAtoms, nat,fxAtoms );
  if (rv != SPOSEAT_Ok)return -10.0;
  double ld=0.0;
  for (int i=0;i<3;i++) {
    //cout <<"matx " << TMatrix[i][0]<< " " << TMatrix[i][1]<< " " << TMatrix[i][2]<< endl;
    ld=ld+fabs(TMatrix[i][i]);
  }
  //cout << "SuperposeResidue " << pResFx->seqNum << " " << pRes->seqNum << " " << nat << " " << ld <<endl;
  return ld;
}

//------------------------------------------------------------------------
double CMMANManager::AtomicRMSDistance( PPCAtom A1, int nA, PPCAtom A2) {
//------------------------------------------------------------------------
  double dd,ddtot = 0.0;

  for (int n=0;n<nA;n++) {
    dd =  (A1[n]->x-A2[n]->x)*(A1[n]->x-A2[n]->x) +  
          (A1[n]->y-A2[n]->y)*(A1[n]->y-A2[n]->y) + 
          (A1[n]->z-A2[n]->z)*(A1[n]->z-A2[n]->z);
    //cout << n << " " << dd << endl;
    ddtot = ddtot+dd;
    }
  ddtot = sqrt(ddtot/nA);
  return ddtot;
    }


int CMMANManager::CopyCoordinates(const PCMMDBManager fromMolHnd,int fromModel) {
  /*
  Copy coordinates - assume same number of atoms in both models
  */
  PPCAtom toAtoms=NULL, fromAtoms=NULL;
  int toSelHnd,fromSelHnd, toNat, fromNat;

  if (GetNumberOfModels() == 1) {
    GetAtomTable1(toAtoms,toNat);
    toSelHnd = -1;
  } else { 
    toSelHnd = NewSelection();
    Select(toSelHnd,STYPE_ATOM,1,"*",ANY_RES,
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
    GetSelIndex(toSelHnd,toAtoms,toNat);
  }
  fromSelHnd = fromMolHnd->NewSelection();
  fromMolHnd->Select(fromSelHnd,STYPE_ATOM,fromModel,"*",ANY_RES,
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  fromMolHnd->GetSelIndex(fromSelHnd,fromAtoms,fromNat);
  if (fromNat != toNat ) {
    cout << "ERROR copying coordinates - different number of atoms" << endl;
    cout << "Atoms in source: " << fromNat << " atoms in target " << toNat << endl;
    fromMolHnd->DeleteSelection(fromSelHnd);
    if(toAtoms) delete [] toAtoms;
    return 1;
  }

  for (int i=0;i<toNat;i++) {
    toAtoms[i]->x = fromAtoms[i]->x;
    toAtoms[i]->y = fromAtoms[i]->y;
    toAtoms[i]->z = fromAtoms[i]->z;
  }
  fromMolHnd->DeleteSelection(fromSelHnd);
  if (toSelHnd>=0) {
    DeleteSelection(toSelHnd);
  } else {
    if(toAtoms) delete [] toAtoms;
  }
  return 0;
}

//---------------------------------------------------------------------------
int CMMANManager::LoadSerial(const PCMMDBManager fromMolHnd ) {
//---------------------------------------------------------------------------

  /*
  Copy serial numbers from one MMDBManager to UDD of another MMDBManager
  */
  PPCAtom toAtoms=NULL, fromAtoms=NULL;
  int toSelHnd,fromSelHnd, toNat, fromNat;

  toSelHnd = NewSelection();
  Select(toSelHnd,STYPE_ATOM,0,"*",ANY_RES,
     "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  GetSelIndex(toSelHnd,toAtoms,toNat);

  fromSelHnd = fromMolHnd->NewSelection();
  fromMolHnd->Select(fromSelHnd,STYPE_ATOM,0,"*",ANY_RES,
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  fromMolHnd->GetSelIndex(fromSelHnd,fromAtoms,fromNat);
  if (fromNat != toNat ) {
    cout << "ERROR copying serial numbers - different number of atoms" << endl;
    cout << "Atoms in source: " << fromNat << " atoms in target " << toNat << endl;
    fromMolHnd->DeleteSelection(fromSelHnd);
    DeleteSelection(toSelHnd);
    return -1;
  }

  // Get a UDD identifier
  int udd_serial =  GetUDDHandle( UDR_ATOM,"atomSerial");
  if (udd_serial<=0)  udd_serial = RegisterUDInteger(UDR_ATOM,"atomSerial" );
  if ( udd_serial < 0 ) {
     printf ( "ERROR registering atom serial data.\n" );
     return -1;
  }
  
  for (int i=0;i<toNat;i++) {
    toAtoms[i]->PutUDData(udd_serial,fromAtoms[i]->serNum);
  }
  fromMolHnd->DeleteSelection(fromSelHnd);
  DeleteSelection(toSelHnd);
  return udd_serial;
}

//---------------------------------------------------------------------------
int CMMANManager::LoadSerialFromDifferentModel(const PCMMDBManager fromMolHnd ,int udd ) {
//---------------------------------------------------------------------------

  /*
  Copy serial numbers from one MMDBManager to UDD of another MMDBManager
  */
  PCAtom pAtom;
  PPCAtom  fromAtoms=NULL,atomTable=NULL;
  int fromSelHnd, fromNat;
  int natoms=0,nCopy = 0;

  GetAtomTable1(atomTable,natoms);
  for (int ia=0; ia<natoms;ia++) atomTable[ia]->PutUDData(udd,0);

  fromSelHnd = fromMolHnd->NewSelection();
  fromMolHnd->Select(fromSelHnd,STYPE_ATOM,0,"*",ANY_RES,
       "*",ANY_RES, "*","*","*","*","*",SKEY_NEW);
  fromMolHnd->GetSelIndex(fromSelHnd,fromAtoms,fromNat);

  for ( int ia=0; ia<fromNat; ia++ ) {
    pAtom = GetAtom ( fromAtoms[ia]->GetModelNum(),
		      fromAtoms[ia]->GetChainID(),
                      fromAtoms[ia]->residue->seqNum,
                      fromAtoms[ia]->residue->insCode,
                      fromAtoms[ia]->name,
                      fromAtoms[ia]->element,
                      fromAtoms[ia]->altLoc );

    if (pAtom) {
      //if (fromAtoms[ia]->serNum==3 ||fromAtoms[ia]->serNum==6 ) cout << "LoadSerialFromDifferentModel " << fromAtoms[ia]->residue->seqNum << " " <<fromAtoms[ia]->serNum << endl;
      pAtom->PutUDData(udd,fromAtoms[ia]->serNum);
      nCopy++;
    }
  }
  return nCopy;
  
}

int CMMANManager::GetNumberOfSecStructure (int type) { 

  int ntype=0;

  int           nr;
  PPCResidue    selRes=NULL;

  GetResidueTable(selRes,nr);
  for( int i = 0; i < nr; i++) {
    if( selRes[i]->isAminoacid() ) {
      if(selRes[i]->SSE==type){
        ntype++;
      }
    }
  }
  return ntype;
}

//-----------------------------------------------------------------------
std::string CMMANManager::PrintSecStructure (void) { 
//-----------------------------------------------------------------------
  int		i;
  const char 	secstr_text [9][12] = { "           ",
                                        "beta strand",
                                        "beta bulge ",
			                "3-turn     ",
                                        "4-turn     ",
                                        "5-turn     ",
                                        "alpha helix",
                                        "H-bond turn",
                                        "bend       " };

  std::ostringstream output;
  int           nr;
  PPCResidue    selRes=NULL;
  std::string resid;

  output << "Secondary Structure Assignment based on Kabsch & Sanders"
         << "\n"
         << "Residue        "
         << "SecStruct   \n";
     

   //   << "Hydrogen bonded to.."


    //Get all residues with atoms in this atom selection
    GetResidueTable(selRes,nr);
 
    for ( i = 0; i < nr; i++) {
      if ( selRes[i]->isAminoacid() ) {
        resid = AtomLabel_residue1(selRes[i]);
        output << std::endl 
             <<  resid <<  std::setw(15-resid.length()) <<" " 
             << secstr_text[selRes[i]->SSE] << " ";
        /*  
        if ( hbonds ) {
          int k = 0;
          while (  k <= 2 && hbonds[i][k] != 0 ) {
            j = selRes[i + hbonds[i][k]];
            resid = molH->AtomLabel_residue1(j);
            output << resid << std::setw(15-resid.length()) << " ";
            k++;
          }
        }
        */
      }
    }
    output << std::endl; 
    return output.str();
}


//------------------------------------------------------------------------
int CMMANManager::TransformToSuperposeCloseAtoms( PCMMANManager fxMolHnd, 
         int fxSelHnd , realtype central_cutoff, realtype cutoff ,
         int mvSuperposeHnd,int fxSuperposeHnd ) {
//------------------------------------------------------------------------
  PPCAtom fxSelAtoms = NULL;
  int fxNAtoms;
  int fxResHnd;
  PPCResidue fxRes = NULL;
  int fxNRes;

  // Max possible number of superposed atoms
  fxMolHnd->GetSelIndex(fxSelHnd,fxSelAtoms,fxNAtoms);
  if (fxNAtoms<=0) return -1;

  
  //Get list of target residues
  fxResHnd = fxMolHnd->NewSelection();
  fxMolHnd->Select(fxResHnd,STYPE_RESIDUE,fxSelHnd,SKEY_NEW);
  fxMolHnd->GetSelIndex(fxResHnd,fxRes,fxNRes);
  //cout << "fxNAtoms,fxNRes " << fxNAtoms << " " << fxNRes << endl;
  if (fxNRes<=0) {
    fxMolHnd->DeleteSelection(fxResHnd);
    return -1;
  }
 
  // Create selection of atoms that will be superposed
  //mvSuperposeHnd = NewSelection();
  //fxSuperposeHnd = fxMolHnd->NewSelection();
  PPCAtom fxSuperpose = NULL;
  PPCAtom mvSuperpose = NULL;
  int fxNSuperpose,mvNSuperpose;


  //Get the 'central' atoms in the moving model
  int centralHnd = NewSelection();
  PPCAtom centralAtoms;
  int centralNAtoms;
  Select(centralHnd,STYPE_ATOM,0,"*",ANY_RES,
       "*",ANY_RES, "*","*","CA","C","*",SKEY_NEW);
  GetSelIndex(centralHnd,centralAtoms,centralNAtoms);
  //cout << "centralNAtoms " << centralNAtoms << endl;

  // Loop over residues in fixed model

  PCAtom pCentralAtom;
  PCResidue closeRes;
  PCAtom fxAtom,closeAtom;
  int ncontacts,best_match;
  double score,best_score;

  for (int ir=0; ir<fxNRes; ir++) {
    // Find the 'central' atoms in moving object that
    // are close to the central atom of fixed residue
    // -- there should just be one!

    pCentralAtom = fxRes[ir]->GetAtom("CA","*","*");
    if (pCentralAtom) {
      PSContact contact= NULL;
      ncontacts = 0;
      SeekContacts (pCentralAtom,centralAtoms,centralNAtoms,
		    0.0, central_cutoff,0,contact,ncontacts,-1);
      best_match=-1;
      best_score=central_cutoff;
      if (ncontacts ==1 ) {
        best_match=0;
      } else if (ncontacts>1) {
        //cout <<  pCentralAtom->residue->seqNum << "ncontacts " << ncontacts << endl;
      
        for (int nc=0;nc<ncontacts;nc++) {
          score=contact[nc].dist -
           DeltaResidueOrientation(centralAtoms[contact[nc].id2]->residue,fxRes[ir]);
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
        for (int ia=0;ia<fxRes[ir]->GetNumberOfAtoms();ia++) {
          fxAtom = fxRes[ir]->GetAtom(ia);
          if ( fxAtom && fxAtom->isInSelection(fxSelHnd) ) {
	    // Is there an atom with same name in the close moving residue
	     closeAtom= closeRes->GetAtom(fxAtom->name,fxAtom->element,"*");
             if ( closeAtom && BondLength(fxAtom,closeAtom)<=cutoff) {
               SelectAtoms(mvSuperposeHnd,closeAtom->serNum,closeAtom->serNum,SKEY_OR);
               fxMolHnd->SelectAtoms(fxSuperposeHnd,fxAtom->serNum,fxAtom->serNum,SKEY_OR);
               
             }
          }
        }
      }
      if (contact) delete [] contact;
    }
  }

  // Cleanup
  fxMolHnd->DeleteSelection(fxResHnd);
  DeleteSelection(centralHnd);


  //Get the list of selected atoms to superpose and do the
  // superposition and transform
  GetSelIndex(mvSuperposeHnd,mvSuperpose,mvNSuperpose);
  fxMolHnd->GetSelIndex(fxSuperposeHnd,fxSuperpose,fxNSuperpose);
  //cout << "mvNSuperpose,fxNSuperpose " << mvNSuperpose << " " << fxNSuperpose << endl;

  if (mvNSuperpose != fxNSuperpose ) {
    return -2;
  } else if ( mvNSuperpose<3 ) {
    return -3;
  }
  int rv = TransformToSuperposeAtoms(mvSuperpose,mvNSuperpose,fxSuperpose);
  if (rv != SPOSEAT_Ok) {
   return -4;  
  } else {
    return mvNSuperpose;
  }
}

int CMMANManager::ApplyCartesiansDeltas(const std::vector<Cartesian> &dxyz, int selHnd, double scale){
  PPCAtom atoms = NULL;
  int natoms;

  GetSelIndex(selHnd,atoms,natoms);
  //std::cout << "natoms: " << natoms << ", displacements size: " << dxyz.size() << "\n";
  if(natoms>0&&((unsigned)natoms!=dxyz.size())){
     std::cout << "Error, wrong number of atom deltas: natoms: " << natoms << ", displacements size: " << dxyz.size() << "\n";
     DeleteSelection(selHnd);
     return 1;
  }
  for(unsigned i=0;i<dxyz.size();i++){
    atoms[i]->x += dxyz[i].get_x()*scale;
    atoms[i]->y += dxyz[i].get_y()*scale;
    atoms[i]->z += dxyz[i].get_z()*scale;
  }

  return 1;
}
