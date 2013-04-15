/*
     mmut/mmut_sbase.cc: CCP4MG Molecular Graphics Program
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
#include <string>
#include <map>
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <mmdb/mmdb_manager.h>
#include "mmut_manager.h"
#include "mmut_sbase.h"
#include <mmdb/mmdb_sbase.h>
#include <mmdb/mmdb_graph.h>
#include <mmdb/mmdb_tables.h>
#include "mgtree.h"
#include "mgutil.h"

using namespace std;

CSBase *CMGSBase::SBase=0;

// Definitions for the compound group class
bool CCompoundGroup::groupMatch[13][13] = { 
 {true ,true ,true ,false,false,false,false,false,false,false,false,false,true}, 
 {false,true ,false,false,false,false,false,false,false,false,false,false,true}, 
 {false,false,true ,false,false,false,false,false,false,false,false,false,true}, 
 {false,false,false,true ,true ,true ,false,false,false,false,false,false,true}, 
 {false,false,false,false,true ,false,false,false,false,false,false,false,true}, 
 {false,false,false,false,false,true ,false,false,false,false,false,false,true}, 
 {false,false,false,false,false,false,true ,true ,true ,false,false,false,true}, 
 {false,false,false,false,false,false,false,true ,true ,true ,true ,false,true}, 
 {false,false,false,false,false,false,false,true ,true ,true ,true ,false,true}, 
 {false,false,false,false,false,false,false,true ,true ,true ,true ,false,true}, 
 {false,false,false,false,false,false,false,true ,true ,true ,true ,false,true}, 
 {false,false,false,false,false,false,false,false,false,false,false,true ,true}, 
 {true ,true ,true ,true ,true ,true ,true ,true ,true ,true ,true ,true ,true} };


char *CCompoundGroup::cifGroupNames[12] = { "peptide", "D-peptide", "L-peptide", 
			     "DNA/RNA", "DNA", "RNA", 
		 "saccharide", "pyranose", "D-saccharid","L-saccharid", 
			     "solvent","non-polymer" };
int CCompoundGroup::groupCode[12] = { RESTYPE_PEPTIDE, RESTYPE_DPEPTIDE,  RESTYPE_LPEPTIDE, 
		      RESTYPE_NUCL,RESTYPE_DNA, RESTYPE_RNA, 
	     RESTYPE_SACH, RESTYPE_SACH, RESTYPE_DSACH, RESTYPE_LSACH, 
		      RESTYPE_SOLVENT, RESTYPE_NONPOLY };

// definition of the CLibAtom atom hydrogen bonding types
int CLibAtom::nHbCodes = 6;
char *CLibAtom::hbCharCode[6] = { "U", "N", "H", "D", "B", "A" };
int CLibAtom::hbCode[6] = { HBTYPE_UNKNOWN,HBTYPE_NEITHER, HBTYPE_HYDROGEN, 
                    HBTYPE_DONOR, HBTYPE_BOTH, HBTYPE_ACCEPTOR };

// definition of CLibBond bond types
 
int CLibBond::nBondCodes = 6;
char *CLibBond::bondCharCode[6] = { "single" , "double", "triple", 
				    "aromatic", "deloc", "metal" };
int CLibBond::bondCode[6] = {  BONDTYPE_SINGLE, BONDTYPE_DOUBLE, 
			       BONDTYPE_TRIPLE, BONDTYPE_AROMATIC, 
                               BONDTYPE_DELOC, BONDTYPE_METAL};

Tree CMGSBase::GetMonomerLibraryTree(const char *monomer_name){
  Tree tree;
  if(SBase) return tree;
  LoadedPCSBStructure monlib;
  CSBStructure *pStruct = GetStructure(monomer_name,monlib);
  if(pStruct){
    //printf("Got %s structure\n",monomer_name);
    std::map<std::string,int> used_elements;
    std::vector<Cartesian> carts;
    std::vector<std::string> labels;
    std::vector<std::vector<int> > bond_lists(pStruct->nBonds);
    std::vector<std::vector<int> > angle_lists(pStruct->nAngles);
    std::vector<std::vector<int> > torsion_lists(pStruct->nTorsions);
    std::vector<double> bond_lengths(pStruct->nBonds);
    std::vector<double> bond_angles(pStruct->nAngles);
    std::vector<double> torsion_angles(pStruct->nTorsions);
    //printf("nAtoms:%d\n",pStruct->nAtoms);
    //printf("nBonds:%d\n",pStruct->nBonds);
    //printf("nAngles:%d\n",pStruct->nAngles);
    //printf("nTorsions:%d\n",pStruct->nTorsions);
    //std::ofstream map_output((std::string(monomer_name)+std::string(".atom_map")).c_str());
    for(int i=0;i<pStruct->nAtoms;i++){
      carts.push_back(Cartesian(pStruct->Atom[i]->x,pStruct->Atom[i]->y,pStruct->Atom[i]->z));
      std::string label = std::string(pStruct->Atom[i]->element);
      used_elements[label]++;
      labels.push_back(label+IntToString(used_elements[label]));
      //map_output << label+IntToString(used_elements[label]) << " --> " << pStruct->Atom[i]->pdb_name << "\n";
      //map_output.flush();
    }
    for(int i=0;i<pStruct->nBonds;i++){
      bond_lists[i].push_back((pStruct->Bond[i]->atom1)-1);
      bond_lists[i].push_back((pStruct->Bond[i]->atom2)-1);
      bond_lengths[i] = pStruct->Bond[i]->length;
      //printf("%d %d: %f\n",pStruct->Bond[i]->atom1-1,pStruct->Bond[i]->atom2-1,bond_lengths[i]);
    }
    for(int i=0;i<pStruct->nAngles;i++){
      angle_lists[i].push_back((pStruct->Angle[i]->atom1)-1);
      angle_lists[i].push_back((pStruct->Angle[i]->atom2)-1);
      angle_lists[i].push_back((pStruct->Angle[i]->atom3)-1);
      bond_angles[i] = pStruct->Angle[i]->angle;
      //printf("%d %d %d: %f\n",pStruct->Angle[i]->atom1-1,pStruct->Angle[i]->atom2-1,pStruct->Angle[i]->atom3-1,bond_angles[i]);
    }
    for(int i=0;i<pStruct->nTorsions;i++){
      torsion_lists[i].push_back((pStruct->Torsion[i]->atom1)-1);
      torsion_lists[i].push_back((pStruct->Torsion[i]->atom2)-1);
      torsion_lists[i].push_back((pStruct->Torsion[i]->atom3)-1);
      torsion_lists[i].push_back((pStruct->Torsion[i]->atom4)-1);
      torsion_angles[i] = pStruct->Torsion[i]->torsion;
      //printf("%d %d %d %d: %f\n",pStruct->Torsion[i]->atom1-1,pStruct->Torsion[i]->atom2-1,pStruct->Torsion[i]->atom3-1,pStruct->Torsion[i]->atom4-1,torsion_angles[i]);
    }
    tree.SetBondsAnglesTorsions(pStruct->nAtoms,bond_lists,bond_lengths,angle_lists,bond_angles,torsion_lists,torsion_angles,chirals);
   }
  return tree;
}

PPCAtom CMGSBase::GetMonomerLibraryStructure(const char *monomer_name){
  PPCAtom mon_atoms=0;
  if(SBase) return 0;
  LoadedPCSBStructure monlib;
  CSBStructure *pStruct = GetStructure(monomer_name,monlib);
  // Residue, chain, model, seqNum to do elsewhere.
  // serNum
  if(pStruct){
    Tree tree = GetMonomerLibraryTree(monomer_name);
    // Then make some atoms, populate mon_atoms, etc. After lunch.
    mon_atoms = new PCAtom[tree.GetNumberOfVertices()];
    std::vector<Cartesian> carts = tree.GetAllCartesians();
    //std::cout << carts.size() << " " << tree.GetNumberOfVertices() <<  " " << pStruct->nAtoms << "\n";
    for(unsigned i=0;i<carts.size();i++){
      //std::cout << carts[i] << ", "<< pStruct->Atom[i]->pdb_name << " " << pStruct->Atom[i]->element << " " << pStruct->Atom[i]->energyType << "\n";
      mon_atoms[i] = new CAtom();
      //mon_atoms[i]->x = carts[i].get_x();
      //mon_atoms[i]->y = carts[i].get_y();
      //mon_atoms[i]->z = carts[i].get_z();
      mon_atoms[i]->SetAtomName(i,i+1,pStruct->Atom[i]->pdb_name,"","",pStruct->Atom[i]->element);
      mon_atoms[i]->SetCoordinates(carts[i].get_x(),carts[i].get_y(),carts[i].get_z(),1.0,0.0); // We should no better than 0 for tFac
      //strncpy(mon_atoms[i]->name,pStruct->Atom[i]->pdb_name,20);
      //strncpy(mon_atoms[i]->element,pStruct->Atom[i]->element,10);
      //strncpy(mon_atoms[i]->energyType,pStruct->Atom[i]->energyType,10);
    }
  }
  return mon_atoms;
}


//-----------------------------------------------------------------------
CMGSBase::CMGSBase(char *monomers_dir_in, char *user_monomers_dir_in, char *sb, 
                    char *ener_lib,   char *mon_lib, char *ele_lib) {
//-----------------------------------------------------------------------
  //char *ener_lib = "ener_lib_m.lib";
  //char *mon_lib = "mon_lib_com.cif";
  
  // Parameters
  nLinks = 0;
  maxAtomInRes = 2000;
  graphSearchMode = 0;
  //Load Structure database
  strcpy(monomers_dir,monomers_dir_in);
  strcpy(user_monomers_dir,user_monomers_dir_in);
  //cout << "monomers_dir " << monomers_dir << endl;
  InitSBase(sb);
  //Load reference file for links/mods etc.
  LoadMonLib(mon_lib);
  LoadEnerLib(ener_lib);
  LoadEleLib(ele_lib);
  
  LoadSynonyms(mon_lib);


}

//-------------------------------------------------------------------
int CMGSBase::InitSBase(char *sb) {
//------------------------------------------------------------------
  if(SBase) return 0;
  int RC,i;
  int nload = 21;
  LoadedPCSBStructure dummy;
  char *load[] = { "ALA","GLY","SER","THR","ASP","GLU","ASN","GLN",
                   "LEU","ILE","PHE","TYR","HIS","CYS","MET","TRP",
		   "ARG","LYS","PRO","VAL","HOH" };
  if (strlen(sb)>0) {
    SBase  = new CSBase ();
    RC = SBase->LoadIndex1("CCP4_SBASE");
    RC = SBase->LoadIndex(sb);
    if (RC!=SBASE_Ok)  {
      printf ( " *****Structure database files not found.\n" );
      delete SBase;
     }

    for (i = 0; i < nload; i++ ) {
      SBase->LoadStructure ( load[i]);
    }

    for (i = 0; i < nload; i++ ) {
      GetStructure(load[i],dummy);
    }
  }  else if ( strlen(monomers_dir)>0 ) {
    
    // Unremediated?
    //for (i = 0; i < nload-1; i++ ) {
    //  GetStructure(load[i],dummy,true);
    //}
    for (i = 0; i < nload; i++ ) {
      GetStructure(load[i],dummy);
    }
    //cout << "InitSBase unrem " << unremedPStruct.size() << " " << loadedPStruct.size() << endl;

  } else {
    printf ( " *****Structure database files not found.\n" );
    return 1;
  }
  fracMatch = 0.9;

  return 0;

}

//------------------------------------------------------------------
PCSBStructure CMGSBase::GetStructure( const ResName resNam , 
                LoadedPCSBStructure monlib, const bool unremediated) {
//------------------------------------------------------------------
  PCSBStructure pStruct=NULL;
  char filename[500];
  int RC=1;


  //cout << "GetStructure " << resNam << " " << unremediated << user_monomers_dir << endl;
  std::string resn  (resNam);

  // Is it in the monomer library specific for the model?
  if ( monlib.size()>0) {
    LoadedPCSBStructure_iter p = monlib.find(resn);
    if (p!=monlib.end()) return p->second;
  }

  // Is it already loaded?
  LoadedPCSBStructure_iter p = loadedPStruct.find(resn);
  if (p!=loadedPStruct.end()) return p->second;

  PCMMCIFFile file = new CMMCIFFile;

  // Look in users monomer library first
  if ( strlen(user_monomers_dir)>1 ) {
    strcpy(filename,user_monomers_dir);
    if (monomers_dir[strlen(monomers_dir)-1] == '\\') 
      strcat(filename,"\\");
    else
      strcat(filename,"/");
    strcat(filename,resNam);
    strcat(filename,".cif\0");
    //cout << "Reading monomer filename " << filename << endl;
    RC = file->ReadMMCIFFile ( filename );
  }

  // Load new structure
  if ( RC &&  strlen(monomers_dir)>1 ) {
    strcpy(filename,monomers_dir);
    if (monomers_dir[strlen(monomers_dir)-1] == '\\') 
      strcat(filename,"\\");
    else
      strcat(filename,"/");
    strcat(filename,resNam);
    strcat(filename,".cif\0");
    //cout << "Reading monomer filename " << filename << endl;
    RC = file->ReadMMCIFFile ( filename );
    if (RC) {
      strcpy(filename,monomers_dir);
      char dir[2];
      *dir = tolower(resNam[0]);
      dir[1] = '\0';
      strcat(filename,dir);
      if (monomers_dir[strlen(monomers_dir)-1] == '\\') 
        strcat(filename,"\\");
      else
        strcat(filename,"/");
      strcat(filename,resNam);
      strcat(filename,".cif\0");
      //cout << "Reading monomer filename " << filename << endl;
      RC = file->ReadMMCIFFile ( filename );
    }
  }
  if (!RC) {
    pStruct = LoadCifMonomer(resNam,file);
  }
  if (file) delete file;

  if ( pStruct ) {
    loadedPStruct[resn]=pStruct;
    return pStruct;
  } else {
    return NULL;
  }


  /*
  // List diagnostic info
  char energyType[energy_type_len+1];
  char pdbName  [pdb_name_len+1];    // PDB atom name
  int nat = pStruct->nAtoms;
  for (int j = 0;j < nat;j++) {
    energyType = pStruct->Atom[j]->energyType;
    pdbName =  pStruct->Atom[j]->pdb_name;
    cout << resNam << " " << pdbName << " " << energyType << endl;
  }
  */
}
  
//-----------------------------------------------------------------------
CMGSBase::~CMGSBase() {
//-----------------------------------------------------------------------
  int i;
  LoadedPCSBStructure_iter p = loadedPStruct.begin();
  while (p != loadedPStruct.end()){
    delete p->second;
    p++;
  }
  if (SBase) delete SBase;
  for ( i = 0; i < nLinks; i++ ){
    if(link[i])
      delete link[i];
  }
}


//-----------------------------------------------------------------------------
PCSBStructure CMGSBase::LoadCifMonomer (const ResName mon, const PCMMCIFFile file ,
                        const bool unscramble ) {
//-----------------------------------------------------------------------------
  // Load the monomer info from CIF file
  //std::cout << "CMGSBase::LoadCifMonomer " << mon << "\n"; std::cout.flush();
  int RC,N;
  PCMMCIFLoop Loop;
  PCMMCIFData dataBlock;
  char dataName[200];
  char id[100];
  char comp_group[100];
  bool add_OXT = false;

  //cout << endl << "LoadCifMonomer " << mon << "*";

  // Get the compound type
  dataBlock = file->GetCIFData ( "comp_list" );
  if (dataBlock) {
    Loop = dataBlock->GetLoop ( CIFCAT_COMP );
     if ( !Loop ) {
      cerr << "ERROR reading CIF loop " << CIFCAT_COMP << " for " << mon << endl;
     } else {
       for (N=0;N<Loop->GetLoopLength();N++) {
         RC = CIFGetString ( id, Loop,CIFTAG_COMP_ID,N,sizeof(id),"");
         if (!RC ) {
           RC = CIFGetString ( comp_group, Loop,CIFTAG_COMP_GROUP,N,sizeof(comp_group),"");
           //cout << "LoadCifMonomer " << mon << " " << comp_group << endl;
           if (!RC) {
             if (strcmp(comp_group,"L-peptide")==0 || strcmp(comp_group,"D-peptide")==0) 
                                                           add_OXT=true;
           }
         }
       }
     }
     //cout << "comp_group " << mon << " " << comp_group << " " << add_OXT <<endl;
  } else 
    cerr << "ERROR reading data block comp_list for " << mon << endl;

  strcpy(dataName,"comp_");
  strcat(dataName,mon);
  //cout << "dataName " << dataName << endl;
  dataBlock = file->GetCIFData ( dataName );
  //cout << "dataBlock " << dataBlock << endl;
  if (!dataBlock) {
    cerr << "ERROR reading block data " << dataName << " for  " << mon << endl;
    return NULL;
  }

  CSBStructure *structure = new CSBStructure;
  strcpy(structure->compoundID,mon);
  structure->PutFormula (comp_group);

  // ----------------------------Get atom data
  Loop = dataBlock->GetLoop ( CIFCAT_COMP_ATOM );
  if ( !Loop ) {
    cerr << "ERROR reading CIF loop " << CIFCAT_COMP_ATOM << " for " << mon << endl;
    delete structure;
    return NULL;
  }
 
  // Create the CChemCompAtom objects and read definition from cif
  int nAtoms = Loop->GetLoopLength();
  Element el;

  for (int N=0;N<nAtoms;N++) {
    CSBAtom *new_atom  = new CSBAtom();
    RC = CIFGetString ( new_atom->pdb_name, Loop,CIFTAG_COMP_ATOM_ID,N,sizeof(new_atom->pdb_name),"");
    if ( RC || (!new_atom->pdb_name) ) {
      delete structure;
      delete new_atom;
      return NULL;
    }
    RC = CIFGetString ( new_atom->sca_name, Loop,CIFTAG_COMP_ALT_ATOM_ID,N,sizeof(new_atom->sca_name),"");
    RC =  CIFGetString (new_atom->energyType, Loop,CIFTAG_COMP_TYPE_ENERGY,N,sizeof(new_atom->energyType),"");
    RC = CIFGetString ( el,Loop,CIFTAG_COMP_TYPE_SYMBOL,N,sizeof(el),"");
    CLibElement::Justify(el,new_atom->element);
    RC = CIFGetReal1( new_atom->ccp4_charge,Loop,CIFTAG_COMP_CHARGE, N );
    structure->AddAtom(new_atom);
    if (add_OXT && strcmp(new_atom->pdb_name,"O")==0) {
      CSBAtom *oxt_atom = new CSBAtom();
      oxt_atom->Copy(new_atom);
      strcpy(oxt_atom->pdb_name,"OXT");
      structure->AddAtom(oxt_atom);
    }
  }

  // -------------------------------------------------Load bonds
  Loop = NULL;
  Loop = dataBlock->GetLoop ( CIFCAT_COMP_BOND );
  if ( !Loop && nAtoms>1) {
      cerr << "ERROR reading CIF loop " << CIFCAT_COMP_BOND << " for " << mon << endl;
      delete structure;
      return NULL;
  } else if ( Loop ) {
    int nBonds =  Loop->GetLoopLength();
    AtomName atom1,atom2;
    for (int N=0;N<nBonds;N++) {
      RC = CIFGetString(atom1,Loop,CIFTAG_COMP_BOND_ATOM_1,N,sizeof(atom1),"");
      if ( RC || (!atom1) ) {
        delete structure;
        return NULL;
      }
      RC = CIFGetString(atom2,Loop,CIFTAG_COMP_BOND_ATOM_2,N,sizeof(atom2),"");

      CSBBond *new_bond =  new CSBBond();
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_bond->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_bond->atom2 = j+1;
      }
      //cout << "bond " << structure->Atom[new_bond->atom1-1]->pdb_name << " " << structure->Atom[new_bond->atom2-1]->pdb_name<< endl;
     //RC = CIFGetString(bt,Loop,CIFTAG_COMP_BOND_TYPE,N,10,"");
     //bondType=encodeBondType(bt);
     realtype length;
     RC = CIFGetReal1( length,Loop, "value_dist", N );
      if (new_bond->atom1>0 && new_bond->atom2>0) {
        new_bond->length = length;
        structure->AddBond(new_bond);
      } else {
        delete new_bond;
      }
    }
    if (add_OXT) {
      strcpy(atom1,"C");
      strcpy(atom2,"OXT");
      CSBBond *new_bond =  new CSBBond();
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_bond->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_bond->atom2 = j+1;
      }
      if (new_bond->atom1>0 && new_bond->atom2>0) {
        new_bond->length = 1.231;
        structure->AddBond(new_bond);
      } else {
        delete new_bond;
      }
    }
  }

  Loop = dataBlock->GetLoop ( "_chem_comp_angle" );
  if ( !Loop && nAtoms>2) {
      //cerr << "ERROR reading CIF loop _chem_comp_angle for " << mon << endl;
      //delete structure;
      //return NULL;
  } else if ( Loop ) {
    int nAngles =  Loop->GetLoopLength();
    AtomName atom1,atom2,atom3;
    for (int N=0;N<nAngles;N++) {
      RC = CIFGetString(atom1,Loop,CIFTAG_COMP_BOND_ATOM_1,N,sizeof(atom1),"");
      if ( RC || (!atom1) ) {
        return structure;
      }
      RC = CIFGetString(atom2,Loop,CIFTAG_COMP_BOND_ATOM_2,N,sizeof(atom2),"");
      RC = CIFGetString(atom3,Loop,"atom_id_3",N,sizeof(atom3),"");
      CSBAngle *new_angle =  new CSBAngle();
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_angle->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_angle->atom2 = j+1;
        if ( strcmp (atom3, structure->Atom[j]->pdb_name) == 0 ) new_angle->atom3 = j+1;
      }
      realtype angle;
      RC = CIFGetReal1( angle,Loop, "value_angle", N );
      if (new_angle->atom1>0 && new_angle->atom2>0 && new_angle->atom3>0) {
        //std::cout << "Adding angle " << atom1 << " " << atom2 << " " << atom3 << "\n";
        new_angle->angle = angle;
        structure->AddAngle(new_angle);
      } else {
        delete new_angle;
      }
    }
    if (add_OXT) {
      strcpy(atom1,"CA");
      strcpy(atom2,"C");
      strcpy(atom3,"OXT");
      CSBAngle *new_angle =  new CSBAngle();
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_angle->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_angle->atom2 = j+1;
        if ( strcmp (atom3, structure->Atom[j]->pdb_name) == 0 ) new_angle->atom3 = j+1;
      }
      if (new_angle->atom1>0 && new_angle->atom2>0 && new_angle->atom3>0) {
        //std::cout << "Adding angle " << atom1 << " " << atom2 << " " << atom3 << "\n";
        new_angle->angle = 120.8;
        structure->AddAngle(new_angle);
      } else {
        delete new_angle;
      }
      strcpy(atom1,"OXT");
      strcpy(atom2,"C");
      strcpy(atom3,"O");
      CSBAngle *new_angle2 =  new CSBAngle();
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_angle2->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_angle2->atom2 = j+1;
        if ( strcmp (atom3, structure->Atom[j]->pdb_name) == 0 ) new_angle2->atom3 = j+1;
      }
      if (new_angle2->atom1>0 && new_angle2->atom2>0 && new_angle2->atom3>0) {
        //std::cout << "Adding angle " << atom1 << " " << atom2 << " " << atom3 << "\n";
        new_angle2->angle = 118.4;
        structure->AddAngle(new_angle2);
      } else {
        delete new_angle2;
      }
    }
  }

  Loop = dataBlock->GetLoop ( "_chem_comp_tor" );
  if ( Loop ) {
    int nTorsions =  Loop->GetLoopLength();
    AtomName atom1,atom2,atom3,atom4;
    for (int N=0;N<nTorsions;N++) {
      RC = CIFGetString(atom1,Loop,CIFTAG_COMP_BOND_ATOM_1,N,sizeof(atom1),"");
      if ( RC || (!atom1) ) {
        return structure;
      }
      RC = CIFGetString(atom2,Loop,CIFTAG_COMP_BOND_ATOM_2,N,sizeof(atom2),"");
      RC = CIFGetString(atom3,Loop,"atom_id_3",N,sizeof(atom3),"");
      RC = CIFGetString(atom4,Loop,"atom_id_4",N,sizeof(atom4),"");
      CSBTorsion *new_torsion =  new CSBTorsion();
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom2 = j+1;
        if ( strcmp (atom3, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom3 = j+1;
        if ( strcmp (atom4, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom4 = j+1;
      }
      realtype torsion;
      RC = CIFGetReal1( torsion,Loop, "value_angle", N );
      if (new_torsion->atom1>0 && new_torsion->atom2>0 && new_torsion->atom3>0 && new_torsion->atom4>0) {
        //std::cout << "Adding torsion " << atom1 << " " << atom2 << " " << atom3 << " " << atom4 << "\n";
        new_torsion->torsion = torsion;
        structure->AddTorsion(new_torsion);
      } else {
        delete new_torsion;
      }
    }
    if (add_OXT) {
      strcpy(atom1,"CA");
      strcpy(atom2,"C");
      strcpy(atom3,"O");
      strcpy(atom4,"OXT");
      //strcpy(atom1,"N");
      //strcpy(atom2,"CA");
      //strcpy(atom3,"C");
      //strcpy(atom4,"OXT");
      CSBTorsion *new_torsion =  new CSBTorsion();
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom2 = j+1;
        if ( strcmp (atom3, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom3 = j+1;
        if ( strcmp (atom4, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom4 = j+1;
      }
      if (new_torsion->atom1>0 && new_torsion->atom2>0 && new_torsion->atom3>0 && new_torsion->atom4>0) {
        //std::cout << "Adding torsion " << atom1 << " " << atom2 << " " << atom3 << " " << atom4 << "\n";
        new_torsion->torsion = 180.0;
        structure->AddTorsion(new_torsion);
      } else {
        delete new_torsion;
      }
    }
  }

  chirals.clear();
  Loop = dataBlock->GetLoop ( "_chem_comp_chir" );
  if ( Loop ) {
    int nChirals =  Loop->GetLoopLength();
    AtomName atom1,atom2,atom3,atom4;
    for (int N=0;N<nChirals;N++) {
      RC = CIFGetString(atom1,Loop,"atom_id_centre",N,sizeof(atom1),"");
      if ( RC || (!atom1) ) {
        return structure;
      }
      CSBTorsion *new_torsion =  new CSBTorsion();
      RC = CIFGetString(atom2,Loop,"atom_id_1",N,sizeof(atom2),"");
      RC = CIFGetString(atom3,Loop,"atom_id_2",N,sizeof(atom3),"");
      RC = CIFGetString(atom4,Loop,"atom_id_3",N,sizeof(atom4),"");
      for ( int j = 0; j < structure->nAtoms; j++ ) {
        if ( strcmp (atom1, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom1 = j+1;
        if ( strcmp (atom2, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom2 = j+1;
        if ( strcmp (atom3, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom3 = j+1;
        if ( strcmp (atom4, structure->Atom[j]->pdb_name) == 0 ) new_torsion->atom4 = j+1;
      }
      char chirality[8];
      RC = CIFGetString( chirality,Loop, "volume_sign",N, 8, "" );
      if(RC) std::cout << "error geting chirality sign\n";
      if (new_torsion->atom1>0 && new_torsion->atom2>0 && new_torsion->atom3>0 && new_torsion->atom4>0) {
        //std::cout << "Adding chirality " << atom1 << " " << atom2 << " " << atom3 << " " << atom4 << " " << chirality << "\n";
	std::vector <int> this_chiral;
	this_chiral.push_back(new_torsion->atom1-1);
	this_chiral.push_back(new_torsion->atom2-1);
	this_chiral.push_back(new_torsion->atom3-1);
	this_chiral.push_back(new_torsion->atom4-1);
	if(!strncmp("positiv",chirality,7))
	  this_chiral.push_back(1);
	else if(!strncmp("negativ",chirality,7))
	  this_chiral.push_back(-1);
	else
	  this_chiral.push_back(0);
        //std::cout << "                 " << new_torsion->atom1 << " " << new_torsion->atom2 << " " << new_torsion->atom3 << " " << new_torsion->atom4 << " " << this_chiral.back() << "\n";
	chirals.push_back(this_chiral);
        delete new_torsion;
      }
    }
  }


   // ------------------------ Scramble atom names into usual PDB order
 
  int l,k;
  AtomName tmp;
  char  names[500][21];

  bool scrambled;
  for (int j=0;j<nAtoms;j++) {
    strcpy(names[j],structure->Atom[j]->sca_name);
    strcat(names[j],"\0");
    //cout << "names *" << names[j] << "*\n";
  }


  for (int j=0;j<nAtoms;j++) {
    l = strlen(names[j]);

    if (unscramble && names[j][0] == 'H' ) {
      scrambled = true;
      for (int i=0;i<nAtoms;i++) {
        if (i!=j && names[i][0]!='H' && strlen( names[i])==l) {
          for (k = 1;k<l;k++) {
            if( names[i][k]!= names[j][k]) break;
          }
          if(k>=l) {
            scrambled = false;
            break;
          }
        }
      }
    } else {
      // not H - dont scramble
     scrambled = false;
    }

    if (scrambled) {
      tmp[0]=names[j][l-1];
      for (int i=1;i<l;i++) tmp[i]=names[j][i-1];
      for (int i=l;i<4;i++)tmp[i]=' ';
      tmp[4]='\0';
    } else if (strlen(names[j])<4 ) {
      strcpy(tmp," ");
      strcat(tmp,names[j]);
      for (int i=0;i<3-l;i++) strcat(tmp," ");
    } else {
      strcpy(tmp,names[j]);
    }
    strcpy(structure->Atom[j]->sca_name,tmp);
    
    //cout << structure->Atom[j]->pdb_name << "*";
  }
  
  //cout << "Load nAtoms " << structure->nAtoms << " "  << structure->nBonds << endl;
  
  return structure;
}

//------------------------------------------------------------------
  int  CMGSBase::LoadMonomerLibrary( char* filename, 
            LoadedPCSBStructure &monlist ) {
//------------------------------------------------------------------
  PCSBStructure pStruct;
  PCMMCIFFile file;
  PCMMCIFLoop Loop;
  PCMMCIFData dataBlock;
  int RC,N;
  char id[500];
  char code[500];
  std::map<std::string,std::string> codes;
  int nmon = 0;

  file = new CMMCIFFile;
  
  RC = file->ReadMMCIFFile ( filename );
  if (RC!=0) return nmon; 

  // Get list of monomers in library
  dataBlock = file->GetCIFData ( "comp_list" );
  if (dataBlock) {
    Loop = dataBlock->GetLoop ( CIFCAT_COMP );
     if ( !Loop ) {
      cerr << "ERROR reading CIF loop " << CIFCAT_COMP << " from " << filename << endl;
     } else {
       for (N=0;N<Loop->GetLoopLength();N++) {
         RC = CIFGetString ( id, Loop,CIFTAG_COMP_ID,N,sizeof(id),"");
         RC = CIFGetString ( code, Loop,CIFTAG_COMP_CODE,N,sizeof(code),"");
         //cout << "id " << id << " code " << code << endl;
         if (strlen(code)==0) strcpy(code,id);
         codes[id] = code;
       }
     }
  }

  delete file;
  file = new CMMCIFFile;  
  RC = file->ReadMMCIFFile ( filename );

  for ( std::map<std::string,std::string>::iterator p = codes.begin(); 
                                     p!=codes.end();p++) {
    const char* pid = p->first.c_str();
    //cout << "pid " << pid << endl;
    pStruct = LoadCifMonomer( pid,file);
    //cout << "pStruct " << pStruct << endl;
    if (pStruct)  {
      monlist[p->second] = pStruct;
      nmon++;
    }
  }

  delete file;
  return nmon;
}


//-------------------------------------------------------------------------
std::string  CMGSBase::ListMonomer(char *mon, bool unremediated) {
//-------------------------------------------------------------------------
  PCSBStructure SBS;
  PCSBAtom      Atom;
  PCSBBond      Bond;
  PCSBAngle     Angle;
  char          S[100];
  int           i;
  LoadedPCSBStructure dummy; 
  std::ostringstream output;
  output.setf(ios::fixed);
  output.setf(ios::showpoint);
  output.setf(ios::right,ios::adjustfield);
  output.precision(3);

  //cout << "ListMonomer mon " << mon << endl;
  SBS = GetStructure ( mon ,dummy,unremediated );
  if (!SBS) return output.str();

  if (SBS)  {
    output << endl << " Compound: " << SBS->compoundID << endl;
    output <<         " ~~~~~~~~~~~~~" << endl << endl ;
    if (SBS->Formula)
      output << " Formula: " << SBS->Formula  << endl;
    if (SBS->Name)
      output << " Name: " << SBS->Name << endl;
    if (SBS->Synonym)
      output << " Synonym: " << SBS->Synonym << endl;
    if (SBS->Charge)
      output << " Charge: " << SBS->Charge << endl;
 
    if (SBS->xyz_source=='N')
      output << endl << endl  << " XYZ coordinate data unavailable" << endl;
    else  {
      output << endl << endl  <<" XYZ coordinate data source: ";
      switch (SBS->xyz_source)  {
        case 'A' : output << "ACD Labs\n";   break;
        case 'R' : output << "RCSB\n" ;      break;
        case 'P' : output << "PDB\n";       break;
        default  : output << "unknown, report as a bug\n";
      }
    }
    output << endl <<
      "Atom|PDB  |Chem| Energy |Chi|    X      |    Y      |     Z\n" <<
      " No |Name |elem|  Type  |   |           |           |\n" <<
      "=================================================================\n" ;
    for (i=0;i<SBS->nAtoms;i++)  {
      Atom = SBS->Atom[i];
      //printf ( " %3i|%4s| %2s |%8s| %1c | %1c "
      output << " " <<
        std::setw(3) << i+1 << "  " << 
	std::setw(4) << Atom->pdb_name << "  " <<
	std::setw(2) << Atom->element << "  " << 
        std::setw(8) << Atom->energyType << "  " << 
	std::setw(1) << Atom->chirality << "  " ;   
	// std::setw(1) << Atom->leaving;
      if (Atom->x>-MaxReal)
	//  "|%10.3f|%10.3f|%10.3f\n"
        output << 
          std::setw(10) << Atom->x << "  " << 
          std::setw(10) << Atom->y << "  " << 
          std::setw(10) << Atom->z << endl ;
      else  output << "---------------------------------\n" ;
    }
 
    output << endl <<
      "Bond| Atom|Atom|Bond |  Length  | Length ESD\n" <<
      " No |  1  | 2  |order|          |\n" <<
      "===============|=====|==========|=============\n" ;
    for (i=0;i<SBS->nBonds;i++)  {
      Bond = SBS->Bond[i];
      switch (Bond->order)  {
        case BOND_SINGLE   : strcpy ( S,"sing" );  break;
        case BOND_DOUBLE   : strcpy ( S,"doub" );  break;
        case BOND_AROMATIC : strcpy ( S,"arom" );  break;
        case BOND_TRIPLE   : strcpy ( S,"trip" );  break;
        default            : strcpy ( S,"unkn" );
      }
      //printf ( "%4i| %3i  %3i|%5s|%10.3f|%10.3f\n",
      output << 
	std::setw(4) << i+1 << "  " << 
        std::setw(3) << Bond->atom1 <<  "  " <<
        std::setw(3) << Bond->atom2 << "  " <<
	std::setw(5) << S << " " << 
        std::setw(10) << Bond->length << " " << 
        std::setw(10) << Bond->length_esd << endl ;
    }
 
    if (SBS->nAngles>0)  {
      output << endl <<
	"Angle|Atom|Atom|Atom|  Angle   |  Angle ESD\n" <<
	"  No | 1  | 2  | 3  |          |\n" <<
	"=====|==============|==========|=============\n" ;
      for (i=0;i<SBS->nAngles;i++)  {
        Angle = SBS->Angle[i];
        //printf ( " %4i| %2i   %2i   %2i |%10.3f|%10.3f\n",
        output << " " <<  
          std::setw(4) << i+1 << "  " << 
          std::setw(3) << Angle->atom1 << "  " <<
          std::setw(3) << Angle->atom2 << "  " <<
	  std::setw(3) << Angle->atom3 << " " <<
          std::setw(10) << Angle->angle << " " << 
          std::setw(10) << Angle->angle_esd << endl;
      }
    }
 
    /*
    if (SBS->nLeavingAtoms>0)  {
      printf ( "\n Leaving atoms:\n" );
      for (i=0;i<SBS->nLeavingAtoms;i++)
        printf ( " %2i.  %3i -> %3i\n",i+1,SBS->leavingAtom[i],
                 SBS->bondedAtom[i] );
    }
    */
 
  }
 
  return  output.str();
}


//------------------------------------------------------------------------
std::string CMGSBase::AssignAtomType ( PCResidue pRes,
       LoadedPCSBStructure monlib,
       std::map<std::string,std::string> &customResSynonym,
       int udd_sbaseCompoundID,int udd_sbaseAtomOrdinal,
       int udd_atomEnergyType, const bool unremediated) {
//------------------------------------------------------------------------

  int RC;
  PPCAtom pAtom;
  int ia, j, k, nAtominRes,nmatch,alt_nmatch;
  bool got_alt_names,try_atom_match,use_alt_names = false;          
  PCSBStructure pSbaseRes = NULL;
  const char *nucleic_acid[] = { "A","G","C","T","U"};
  const char *nucleic_acid3[] = { "ADE","GUA","CYT","GUA","THY","URA" };
  const char *nucl_flags[] = { "+","r","d" }; 
  int n_nucleic_acid = 5;

  ivector imatch;

  int tolMatch;
  int Hflag;
  int itype;

  char AtomID[30];
  AtomName name1,name2,name3;
  std::ostringstream output;

  //PPCAtom atomTable;
  //int nAtomTable;

  pRes->GetResidueID ( AtomID );
  //cout << "AssignAtomType " << AtomID << " " << unremediated << endl;

  // Get table of atoms in residue
  pAtom = NULL;
  pRes->GetAtomTable1(pAtom,nAtominRes);
  pRes->PutUDData(udd_sbaseCompoundID,"" );


  // 1) Priorty is user entered aliases
  std::map<std::string,std::string>::iterator p = customResSynonym.find(pRes->name);
  if (p!=customResSynonym.end()) {
    pSbaseRes = GetStructure(p->second.c_str(),monlib, unremediated);
    pRes->GetResidueID ( AtomID );
    if (pSbaseRes) {
      output <<  AtomID << "using user specified synonymous monomer from database: " <<                                                            p->second << endl;
    } else {
      output <<  AtomID << "attempting to use synonymous monomer from database: " << p->second << ", but monomer not found in database" << endl;
    } 
  } 

  // 2) Is there structure of the same name in the database
  if (!pSbaseRes) {
    //Beware - for unremediated PDBs with nucleic acid don't use
    //single letter code monomers
    bool skip = false;
    if (unremediated && strlen(pRes->name)==1) {
      for ( int in = 0; in < n_nucleic_acid; in++ ) {
        if  (strncmp(pRes->name,nucleic_acid[in],1)==0 ) skip = true;
      }
    }
    if (!skip)
      pSbaseRes = GetStructure(pRes->name , monlib, unremediated);
      //  << " " << pSbaseRes->nAtoms << endl;
  }
  

    // 3) Try standard aliases
    if (!pSbaseRes) {
      map<std::string,std::string>::iterator p = synonyms.find(pRes->name);
      if (p!=synonyms.end()) {
        pSbaseRes = GetStructure(p->second.c_str(),monlib, unremediated);
        if (pSbaseRes) {
          pRes->GetResidueID ( AtomID );
          /*
          std::map<std::string,int>::iterator rep_err = reported_errors.find(resn);
          if (rep_err == reported_errors.end() or ) reported_errors[resn]=0;
	  reported_errors[resn] =  reported_errors[resn] + 1;
          if  (reported_errors[resn]<=3)
          */
            output <<  AtomID << "using synonymous monomer from database: " << 
                                                           p->second << endl;
        }
      }
    } 

    // 4) Is it DNA/RNA?
    // Assume its water if its one oxygen atom
    // Set residue and atom type irrespective of their names 
    if (!pSbaseRes) {
      ResName resn = "";
      switch (strlen(pRes->name)) {
      case 1:
        for ( int in = 0; in < n_nucleic_acid; in++ ) {
          if  (strncmp(pRes->name,nucleic_acid[in],1)==0 ) {
            strcpy(resn,nucleic_acid[in]) ;
            break;
          }
        }
        break; 
      case 2:
        if ( strncmp(pRes->name,nucl_flags[0],1) == 0) {
          for ( int in = 0; in < n_nucleic_acid; in++ ) {
            if  (strstr(pRes->name,nucleic_acid[in])!=NULL ) strcpy(resn,nucleic_acid[in]) ;
          }
        } else if ( strstr(pRes->name,nucl_flags[1]) != NULL) {
          for ( int in = 0; in < n_nucleic_acid; in++ ) {
            if  (strstr(pRes->name,nucleic_acid[in])!=NULL ) strcpy(resn,nucleic_acid[in]) ;
          }
        } else if ( strstr(pRes->name,nucl_flags[2]) != NULL) {
          for ( int in = 0; in < n_nucleic_acid; in++ ) {
            if  (strstr(pRes->name,nucleic_acid[in])!=NULL ) strcpy(resn,nucleic_acid[in]) ;
          }
        }
        break;
      case 3:
        for ( int in = 0; in < n_nucleic_acid; in++ ) {
          if  (strcmp(pRes->name,nucleic_acid3[in])!=0 ) strcpy(resn,nucleic_acid[in]) ;
        }
        break;
      }
    
      if ( strcmp(resn,"")!=0 ) {
        // DNA/RNA? is there an O2*
        char nucl_type = 'D';
        for ( ia = 0; ia < nAtominRes; ia++ ) {
          if (strcmp(pAtom[ia]->name," O2*")==0 || strcmp(pAtom[ia]->name," O2'")==0 ) {
            nucl_type = 'R';
            break;
          }
        }
        resn[1]=nucl_type;
        resn[2]='\0';
        //cout << "Nucl resn " << pRes->name << " *" << resn << "*" << nucl_type <<endl;
        pSbaseRes = GetStructure(resn,monlib, unremediated);
        //if (pSbaseRes) {
        //  pstr name = pSbaseRes->Name;
        //  cout << "pSbaseRes name " << name << endl;
        //} else { 
        //  cout << "no pSbaseRes" << endl;
        //}
      } else if (pRes->nAtoms == 1 || pRes->nAtoms == 3 ) {
        pSbaseRes = GetStructure("HOH",monlib, unremediated);
      }
    }
 
  
  if (!pSbaseRes) {
    pRes->GetResidueID ( AtomID );
    output << AtomID << " No monomer in database: " <<  pRes->name << endl;
  // There is matching residue name
  } else  {
    pRes->PutUDData(udd_sbaseCompoundID,pSbaseRes->compoundID );

    // Do atom names match?
    
    //nmatch = 0;
    //for ( ia = 0; ia < nAtominRes; ia++ ) {
    //  for ( j = 0; j < pSbaseRes->nAtoms; j++ ) {
    //    if ( strcmp (pAtom[ia]->name, pSbaseRes->Atom[j]->pdb_name) == 0 ||
    //         strcmp (pAtom[ia]->name, pSbaseRes->Atom[j]->sca_name) == 0) {
    //      pAtom[ia]->PutUDData(udd_sbaseAtomOrdinal,j);
    //      pAtom[ia]->PutUDData(udd_atomEnergyType, 
    //	    LibAtom(pSbaseRes->Atom[j]->energyType,pAtom[ia]->element));
    //      pAtom[ia]->GetAtomID(AtomID);
    //      nmatch++;
    //      break;
    //    }
    //  }
      //pAtom[ia]->GetUDData(udd_sbaseAtomOrdinal,itype); 
      //if ( itype < 0 ) cout << "unmatched " <<   pRes->name << pAtom[ia]->name << endl;
    //}

   // Can stop now if all matched
     //if ( nmatch == nAtominRes ) { 
     // delete [] pAtom; return output.str();
     //}
   // Try again ignoring spaces in names
    try_atom_match = true;
    got_alt_names = true;
    bool iflip = false;
    if (!strncmp(pSbaseRes->Atom[0]->sca_name,"   ",3))got_alt_names=false;
    while (try_atom_match) {
      nmatch = 0;
      alt_nmatch = 0;
      for ( ia = 0; ia < nAtominRes; ia++ ) {
        strcpy_css(name1,pAtom[ia]->name);
        j = 0; 
        while ( j < pSbaseRes->nAtoms) {
          if (use_alt_names && got_alt_names ) {
            strcpy_css(name3,pSbaseRes->Atom[j]->pdb_name);
            strcpy_css(name2,pSbaseRes->Atom[j]->sca_name);
          } else {
            strcpy_css(name2,pSbaseRes->Atom[j]->pdb_name);
            strcpy_css(name3,pSbaseRes->Atom[j]->sca_name);
          }
          if ( strcmp (name1,name2) == 0 ) {
            pAtom[ia]->PutUDData(udd_sbaseAtomOrdinal,j);
	    if(strcmp(pAtom[ia]->element,"  ")==0){
               strcpy(pAtom[ia]->element,pSbaseRes->Atom[j]->element);
            } else if(strncmp(pAtom[ia]->element,pSbaseRes->Atom[j]->element,2)!=0){
               output << "Warning, element name " << pAtom[ia]->element << " does not match dictionary name " << pSbaseRes->Atom[j]->element << " for " << AtomID <<  "\n";
               strcpy(pAtom[ia]->element,pSbaseRes->Atom[j]->element);
            }
            pAtom[ia]->PutUDData(udd_atomEnergyType, 
  	       LibAtom(pSbaseRes->Atom[j]->energyType,pAtom[ia]->element));
            nmatch++;
            j = 9999;
          } 
          if ( strcmp (name1,name3) == 0 ) {
	    alt_nmatch++;
          }
          j++;
        }
      }
      //cout <<  pRes->name << " " << use_alt_names << " " << got_alt_names<< " " << nAtominRes << " " << nmatch << " " << alt_nmatch << endl;
      if (alt_nmatch<=nmatch||iflip) {
        try_atom_match = false;
      } else {
        use_alt_names = !use_alt_names;
        iflip = true;
      }
    }



    // Can stop now if all matched
    if ( nmatch == nAtominRes ) {
      delete [] pAtom; 
      HandleTerminii(pRes,udd_atomEnergyType);
      return output.str();
    }

   // List the unmatched residues
    pRes->GetResidueID ( AtomID );

    output << AtomID << "   " << nmatch << " of " <<  nAtominRes << 
    " atom names match database " << pRes->name << " Unmatched atoms: ";
    for ( ia = 0; ia < nAtominRes; ia++ ) {
      pAtom[ia]->GetUDData(udd_sbaseAtomOrdinal,itype); 
      if ( itype < 0 ) output << pAtom[ia]->name << "  ";    
      //if ( itype < 0 ) cout << "unmatched " <<  pAtom[ia]->name << endl;
      
    }
    output << "\n";
    
    tolMatch = pRes->GetNumberOfAtoms(false) -
      (int)ceil(fracMatch * (float)pRes->GetNumberOfAtoms(false));
    if ( tolMatch < 2 )  tolMatch =2 ;

    // Clean up any matches set in the mmdb data structure
    /*
    pRes->PutUDData(udd_sbaseCompoundID,"" );
    for ( ia = 0; ia < nAtominRes; ia++ ) {
      pAtom[ia]->PutUDData(udd_sbaseAtomOrdinal,-1);
    }
    */
       
    // Atom names do not match - try matching graphs for the same
    // name monomer in the library
    // But only do this for limited number of atoms

    if ( nAtominRes> 1 && nAtominRes < 10 ) {
      
      GetVectorMemory ( imatch,2000 , 0);
      for ( j = 0; j < 2000 ; j++ ) imatch[j] = -1;
    
      // Are there H atoms present? - Set the Hflag for CheckResidue
      // And alternate conformers ???
      Hflag = 1;
      for  ( ia = 0; ia < nAtominRes; ia++ ) {
        if ( !strcmp(pAtom[ia]->element ," H") ) Hflag = 0;
      }

      rvector  occupancy;
      int nAlt,alflag;
      PAltLoc  aL=NULL;
      pRes->GetAltLocations ( nAlt,aL,occupancy,alflag );
      //cout << "Hflag " << Hflag << " nAlt " << nAlt << endl;
      // Get max number of residues mismatch 
      if (nAlt > 1) {
        pRes->GetResidueID ( AtomID );
        //cout <<AtomID << " nAlt " << nAlt << " " << alflag << " " << Hflag  << endl;
        for ( int iAlt = 0; iAlt < nAlt; iAlt++ ) {
          if (strlen( aL[iAlt])>0) { 
            RC = MatchGraphs( pRes, Hflag, false,aL[iAlt], 
                      pSbaseRes, nmatch,imatch,tolMatch );
            //cout << "Alt " <<  aL[iAlt] << " " << nSbase << " " << nmatch << endl;
          }
        }
      } else {
        RC = MatchGraphs( pRes, Hflag, false, "", pSbaseRes, nmatch,imatch,tolMatch );


        //printf("%s MatchGraphs Nmatch %i nAtominRes %i\n",AtomID,nmatch,nAtominRes);
      }
      delete [] aL;
      FreeVectorMemory (occupancy,0 );

      // The graphs match - so is good match for monomer but has
      // 'wrong' atom names 
      if ( nmatch > nAtominRes-tolMatch ) {
        output << "...graph matching matched " << nmatch << 
                    " atoms in  residue \n";
        pRes->PutUDData(udd_sbaseCompoundID,pSbaseRes->compoundID );
        for ( ia = 0; ia < nAtominRes; ia++ ) {
         k = imatch[ia];
         if ( k >= 0 ) {
            pAtom[ia]->PutUDData(udd_sbaseAtomOrdinal,k);
            pAtom[ia]->PutUDData(udd_atomEnergyType,
              LibAtom( pSbaseRes->Atom[k]->energyType,pAtom[ia]->element));
          } 
          else {
            // This atom is not part of the graph matching fragment
            pAtom[ia]->PutUDData(udd_atomEnergyType,LibAtom("",pAtom[ia]->element));
          }  //pAtom[ia]->GetUDData(udd_atomEnergyType,itype);
        }
        FreeVectorMemory(imatch,0);
	delete [] pAtom; 
        HandleTerminii(pRes,udd_atomEnergyType);
	return output.str();
      } 
      else {
        output <<  "...............graph matching failed to match residue\n";
        FreeVectorMemory(imatch,0);
      }

    }
  }

  // There was no name match - try graph matching

  // Graph matching switched off - so give up
  if ( graphSearchMode == 0) {
    // Residue type not recognised so assign a default type for the
    // atom element type
    for (j=0;j < nAtominRes;j++) {
      pAtom[j]->GetUDData(udd_atomEnergyType,itype);
      if (itype <= 0)
        pAtom[j]->PutUDData(udd_atomEnergyType,LibAtom("",pAtom[j]->element));
      //cout << "Assigned atom type " << LibAtom("",pAtom[j]->element) << endl;
    }
    //return 1;
  }


  //else {
    //RC = GraphSearch ( pRes, pSbaseRes, nAtominRes, nMatchAtom, matchAtom);
  //}

  HandleTerminii(pRes,udd_atomEnergyType);
  delete [] pAtom;
  return output.str();
}


int CMGSBase::HandleTerminii(PCResidue pRes, int udd_atomEnergyType ) {
  // Check for peptide chain terminii to fix the typing of N and O
  PCAtom pAtom;
  char *Nnames[] = { "N", "NT" };
  char *Onames[] = { "O", "OXT" , "OE" };
  if ( isAminoacid(pRes->name) ) {
    //cout << "isAminoacid " << pRes->name << pRes->seqNum << endl;
    if ( pRes->index == 0 ) {
      for (int ia=0; ia<2; ia++ ) {
        pAtom = pRes->GetAtom(Nnames[ia],"*","*");
        //cout << "first " << pRes->name << pRes->seqNum << pAtom << " " << LibAtom("NC2","N") <<endl;
        if (pAtom) pAtom->PutUDData(udd_atomEnergyType,LibAtom("NC2","N"));
      }
    }
    // Test if last in chain or the next residue is not amino acid
    else if ( (pRes->index == pRes->chain->GetNumberOfResidues()-1) ||
	      !( isAminoacid(pRes->chain->GetResidue(pRes->index+1)->name) ) ) {
      for (int ia=0; ia<3; ia++ ) {
        pAtom = pRes->GetAtom(Onames[ia],"*","*");
        //cout << "last " << pRes->name << pRes->seqNum << pAtom << " " << LibAtom("OC","O") << endl;
        if (pAtom) pAtom->PutUDData(udd_atomEnergyType,LibAtom("OC","O"));
      }
    }
  }
  return 0;
}


int CMGSBase::MatchGraphs(PCResidue pRes,int Hflag, Boolean Cflag, 
                          const pstr altLoc, 
                          PCSBStructure pSbaseRes, int &nMatched,
			  ivector match, int tolMatch ) {
  
  PCGraph G,G1;
  PCGraphMatch U;
  ivector      F1,F2;
  realtype     p1,p2;
  int     rc,htype;
  int nInResidue,minMatch,natMatch;

  /*
  if ( strlen(monomers_dir)<=0 ) {
    // Using the EBI monomers with pre-built graphs
    SBase->CheckResidue ( pRes, Hflag, Cflag,nInResidue, nInStructure,
		   nMatched, match, altLoc, minMatchSize);
    return 0;
  }
  */

  // Using LibCheck/Refmac libraries - need to do it manually
  // Create a graph of the model residue             
  //cout << " MatchGraphs " <<   pRes->name << " " << pRes->GetNumberOfAtoms() << " atoms in dbase frag " << pSbaseRes->nAtoms << endl;                                                               
  G = new CGraph ( pRes,altLoc );
  if (Hflag>=1) {
    htype = getElementNo(pstr("H"));
    if (Hflag==2)  G->HideType    ( htype );
       else  G->ExcludeType ( htype );
  }
                                                                                
  G->Build ( False );
  nInResidue = G->GetNofVertices();
  if (nInResidue<=0)  {
    rc = SBASE_NoAtomsFound;
    nMatched = 0;
    delete G;
    return rc;
  }

  // Create a graph of the SBase monomer
  G1 = new CGraph();
  std::vector<int> ia_index;
  int nn = 0;
  std::vector<int> ia_revert;
  for (int ia=0;ia<pSbaseRes->nAtoms;ia++) {
    if (Hflag<=0 || strcmp(pSbaseRes->Atom[ia]->element," H") != 0 ) {
      G1->AddVertex ( new CVertex(
                 getElementNo(pSbaseRes->Atom[ia]->element),
                 pSbaseRes->Atom[ia]->pdb_name) );
      ia_index.push_back(++nn);
      ia_revert.push_back(ia);
      //cout <<  ia << " " << getElementNo(pSbaseRes->Atom[ia]->element) << " " <<  pSbaseRes->Atom[ia]->pdb_name << " " << ia_index[ia] << endl;
    } else {
       ia_index.push_back(-1);
    }
  }
  int ia1,ia2;
  for (int ib=0; ib< pSbaseRes->nBonds; ib++) {
    ia1 = pSbaseRes->Bond[ib]->atom1 - 1;
    ia2 = pSbaseRes->Bond[ib]->atom2 - 1;
    if ( ia_index[ia1] >= 0 && ia_index[ia2] >= 0 ) {
      G1->AddEdge ( new CEdge(ia_index[ia1],ia_index[ia2], BOND_SINGLE) );
      //cout << "Edge " << ia_index[ia1] << " " << ia_index[ia2] << endl;
    }
  }

  G1->Build(false);

  //int nInMonomer =  G1->GetNofVertices();
  //cout << "Created graphs with " << nInResidue << " and " <<  nInMonomer << " vertices\n";

  if (nInResidue>10) 
    minMatch = nInResidue-tolMatch;
  else
    minMatch = nInResidue;

  // Try matching the two graphs
  U = new CGraphMatch();
  U->MatchGraphs ( G,G1,minMatch );

  
  nMatched =  U->GetNofMatches();

  
  for (int j=0;j<=nInResidue;j++) match[j]=-1;
  if (nMatched>0) {
    U->GetMatch ( 0,F1,F2,natMatch,p1,p2 );
    int maxM=natMatch; int maxI=0;
    for (int ii=1;ii<nMatched;ii++) {
      U->GetMatch ( 0,F1,F2,natMatch,p1,p2 );
      //cout << "MatchGraphs ii,natMatch " << ii << " " << natMatch<< endl;
      if (natMatch>maxM) {
        maxM=natMatch;
        maxI =ii;
      }
    }

    U->GetMatch (maxI ,F1,F2,natMatch,p1,p2 );
    for (int j=1;j<=natMatch;j++) {
      //cout << "j,F1[j],F2[j],iarevert" << j << " " << F1[j] << " " << F2[j] << " " <<ia_revert[F2[j]-1] << endl;
      match[F1[j]-1] =ia_revert[F2[j]-1];
    }
  }

  delete U;
  delete G;
  delete G1;
  return 0;
}

//------------------------------------------------------------------------
std::string CMGSBase::ListAtomType ( PCMMUTManager molHnd, PCResidue pRes,
        int udd_sbaseCompoundID,
        int udd_sbaseAtomOrdinal, int udd_atomEnergyType ) {
//------------------------------------------------------------------------
 
  PPCAtom pAtom;
  int nAtominRes,ia, j;
  pstr compoundID = NULL;
  int atomType;
  std::ostringstream output;

  pRes->GetUDData(udd_sbaseCompoundID,compoundID );
  //pstr res_label =  molHnd->AtomLabel_residue1(pRes);
  
  pAtom = NULL;
  pRes->GetAtomTable1(pAtom,nAtominRes);
  output << molHnd->AtomLabel_residue1(pRes) << "    "  << compoundID << endl;

  for (ia=0;ia<nAtominRes;ia++) {
    pAtom[ia]->GetUDData(udd_sbaseAtomOrdinal,j);
    pAtom[ia]->GetUDData(udd_atomEnergyType,atomType); 

    output << "    "  << pAtom[ia]->name;
    if ( strlen(pAtom[ia]->altLoc) == 0 ) 
      output << "   ";
    else
      output << "," << pAtom[ia]->altLoc;
    output <<  " " << pAtom[ia]->element << "  " << j << " " <<
           atomType << endl;
  }
  return output.str();
}

//-------------------------------------------------------------------
int CMGSBase::LoadSynonyms( pstr filename ) {
//-------------------------------------------------------------------
  int RC;
  CFile f;
  PCMMCIFLoop Loop=NULL;
  PCMMCIFFile file;
  PCMMCIFData dataBlock;
  ResName id,alt_id;
  std::string str_id, str_alt_id;

  file = new CMMCIFFile;
  RC = file->ReadMMCIFFile ( filename );
  //printf("file %s %i RC %i\n",filename,file,RC);
  dataBlock = file->GetCIFData ( "comp_synonym_list" );
  if (dataBlock) Loop = dataBlock->GetLoop (CIFCAT_COMP_SYNONYM  );
  if ( !Loop ) {
    printf ("ERROR reading CIF loop %s from %s\n",CIFCAT_COMP_SYNONYM,filename);
    if (file) delete file;
    return -1;
  }
  
  for (int n=0;n<Loop->GetLoopLength();n++) {
    RC = CIFGetString ( id, Loop,CIFTAG_COMP_SYNONYM_ID,n,sizeof(id),"");
    if (!RC)
      RC = CIFGetString ( alt_id, Loop,CIFTAG_COMP_SYNONYM_ALTERNATIVE,n,sizeof(alt_id),"");
    str_id = id;
    str_alt_id = alt_id;  
    if (!RC) synonyms[alt_id] = id;
    //cout << "synonym " << alt_id << " " << synonyms[alt_id] << endl;
  }
  if (file) delete file;
  return 0;

}

//------------------------------------------------------------------
int CMGSBase::LoadMonLib (pstr filename ) {
//------------------------------------------------------------------
  // Load the monomer link info from CIF file
  int RC,n;
  CFile f;
  PCMMCIFLoop Loop1=NULL;
  PCMMCIFFile file;
  PCMMCIFData dataBlock;
  char dataName[200];

 
  file = new CMMCIFFile;
  RC = file->ReadMMCIFFile ( filename );
  //printf("file %i RC %i\n",file,RC);
  dataBlock = file->GetCIFData ( "link_list" );
  if (dataBlock) Loop1 = dataBlock->GetLoop ( CIFCAT_LINK );
  if ( !Loop1 ) {
    printf ("ERROR reading CIF loop %s from %s\n",CIFCAT_LINK,filename);
    if (file) delete file;
    return -1;
  }
 
  // Create the link objects and read link definition from mon_lib_com.cif
  nLinks = -1;
  RC = 0;
  do {
    nLinks++;
    link[nLinks] = new MGCLink();
    RC = link[nLinks]->GetCif ( Loop1 , nLinks );
  }while (!RC);
  delete link[nLinks];
  nLinks = nLinks -1;
  // read details of each link from the appropriate datablock
  for ( n = 0; n < nLinks; n++ ) {
    strcpy(dataName,"link_");
    strcat(dataName,link[n]->id);
    dataBlock = file->GetCIFData ( dataName );
    //printf("dataBlock %s %i\n",dataName,dataBlock);
    if (dataBlock) {
      RC = link[n]->GetCifBond(dataBlock);
    }
  }

  //if (nLinks>=0) for ( n = 0; n <= nLinks; n++ ) link[n]->Print();

  if (file) delete file;
   

  return 0;
}

//------------------------------------------------------------------
int CMGSBase::AddLink ( pstr resName1 , pstr atmName1, pstr resName2 , pstr atmName2  ) {
//------------------------------------------------------------------
  nLinks++;
  //cout << "CMGSBase::AddLink " << nLinks << " " << resName1 << " " << atmName1 << " " <<  resName2 << " " << atmName2 << endl;
  link[nLinks] = new MGCLink();
  link[nLinks]->lg2.Set( resName1,"","*", atmName1);
  link[nLinks]->lg1.Set( resName2,"", "*", atmName2);
  return nLinks;
}

//-----------------------------------------------------------------
MGCLink::MGCLink() : lg1(),lg2(){
//-----------------------------------------------------------------
  strcpy(id,"");
  //name = "";
  
}

//-----------------------------------------------------------------
MGCLink::~MGCLink() {
//----------------------------------------------------------------

}

//----------------------------------------------------------------
int  MGCLink::GetCif(  PCMMCIFLoop Loop1, int N ) {
//----------------------------------------------------------------
  int RC;
  pstr cmp,modif,grp,atm;
 
  atm = "";
 
   if (!Loop1) return -1;
  if (N >= Loop1->GetLoopLength()) return -1;
   
  RC = CIFGetString (id,Loop1,CIFTAG_LINK_ID,N,sizeof(id),"");
  if ( RC || (!id) ) return -1;
  //name =  Loop1->GetString(CIFTAG_LINK_NAME,N,RC);
  cmp = Loop1->GetString (CIFTAG_LINK_COMP_ID_1,N,RC);
  modif = Loop1->GetString (CIFTAG_LINK_MOD_ID_1,N,RC);
  grp = Loop1->GetString (CIFTAG_LINK_GROUP_COMP_1,N,RC);
  lg1.Set( cmp, modif, grp, atm);
  cmp = Loop1->GetString (CIFTAG_LINK_COMP_ID_2,N,RC);
  modif = Loop1->GetString (CIFTAG_LINK_MOD_ID_2,N,RC);
  grp = Loop1->GetString (CIFTAG_LINK_GROUP_COMP_2,N,RC);
  lg2.Set( cmp, modif, grp, atm);
  return 0;
     
}

//---------------------------------------------------------------
int MGCLink::GetCifBond (PCMMCIFData dataBlock ){
//---------------------------------------------------------------
  int RC;
  PCMMCIFLoop Loop;
  pstr at;

  Loop = dataBlock->GetLoop ( CIFCAT_LINK_BOND );
  if(!Loop) {
    //printf("ERROR reading _link_bond\n");
    return -1;
  }
  at = Loop->GetString (CIFTAG_LINK_BOND_ATOM1,0,RC);
  if (!RC) strcpy(lg1.atom,at);
  at = Loop->GetString (CIFTAG_LINK_BOND_ATOM2,0,RC);
  if (!RC) strcpy(lg2.atom,at);
  return 0;
}

void MGCLink::Print() {

  printf ("Link %s *%s* %s %s   *%s* %s %s \n",id,lg1.compId,lg1.modId,
    lg1.atom, lg2.compId, lg2.modId,lg2.atom);
}


//------------------------------------------------------------------------
MGCLinkGroup::MGCLinkGroup ( ): group() {
//------------------------------------------------------------------------
  strcpy(compId,"");
  strcpy(modId,"");
  strcpy(atom,"");
  group.Set(-1);
}

//------------------------------------------------------------------------
void MGCLinkGroup::Set ( char *comp, char *modif, char *grp, char *atm ) {
//------------------------------------------------------------------------
  if ( comp) strcpy(compId,comp);
  if ( modif) strcpy(modId,modif);
  if (atm ) strcpy(atom,atm);
  group.Set( grp);
  
}

//-----------------------------------------------------------------------
MGCLinkGroup::~MGCLinkGroup () {
//-----------------------------------------------------------------------
 
}

void MGCLinkGroup::Print() {
//printf ("LinkGroup %i %s %s\n",compId,modId,atom);
}

//------------------------------------------------------------------------
bool MGCLinkGroup::Match ( int grp, pstr comp, pstr atm ) {
//------------------------------------------------------------------------
  //if (group.code < 0 || group.Match(grp)) {
  // printf ( "match group %i %i comp *%s* *%s* atom *%s* *%s* \n", group.code,grp,comp,compId,atom,atm);
  //  cout << group.Match(grp) << " " << strcmp(comp,compId) << " " << strcmp(atom,atm) << endl;
    //}
  if ( (group.code < 0 || group.Match(grp)) &&
       ((strlen(compId)<1 ) || (strcmp(comp,compId)==0) ) 
       && ((strlen(atom)<1) || (strcmp(atom,atm)==0) ) ) {
    //printf ("MATCH\n");
    return true;
  }
  else
    return false;
  
}

//------------------------------------------------------------------------
CCompoundGroup::CCompoundGroup ( ) {
//----------------------------------------------------------------------
}
//------------------------------------------------------------------------
void CCompoundGroup::Set ( int cd ) {
//------------------------------------------------------------------------
  code = cd;
}
//------------------------------------------------------------------------
void CCompoundGroup::Set ( pstr name ) {
//------------------------------------------------------------------------
  code = GetCifGroupCode ( name );
  //printf("Compound %s %i\n",name,code);
}

//------------------------------------------------------------------------
int CCompoundGroup::GetCifGroupCode ( pstr name ) {
//------------------------------------------------------------------------
  int i;
  int retCode = RESTYPE_UNKNOWN;

  if ( !name || strcmp(name,"*")==0) return -1;

  for ( i = 0; i < 11; i++ ) {
    if ( strcmp ( name , cifGroupNames[i]) == 0 ) {
      retCode = groupCode[i];
      break;
    }
  }
  return retCode;
}


//---------------------------------------------------------------
bool CCompoundGroup::Match ( int cd ) {
//---------------------------------------------------------------
  return groupMatch[code][cd];
}


//------------------------------------------------------------------
int CMGSBase::LoadEnerLib (pstr filename ) {
//------------------------------------------------------------------
  // Load the monomer link info from CIF file
  int RC;
  PCMMCIFLoop Loop1;
  CIF =  new CMMCIFData();
  RC = CIF->ReadMMCIFData ( filename );
  if ( RC ) {
    printf ("Error reading %s\n",filename);
    delete CIF;
    return -1;
  }
  Loop1 = CIF->GetLoop ( CIFCAT_LIBATOM );
  if ( !Loop1 ) {
    delete CIF;
    return -1;
  }
  
  nLibAtoms = 0;
  RC =0;
  do {
    CLibAtom new_libAtom = CLibAtom();
    RC =  new_libAtom.GetCif ( Loop1 , nLibAtoms );
    if ( !RC) {
      libAtom.push_back(new_libAtom);
      nLibAtoms++;
    }
  }while (!RC );
 
  /*
  printf ( "nLibAtoms %i\n",nLibAtoms);
  for ( int i = 0; i < nLibAtoms; i++) 
    printf ("Atom type %s %s \n",
       libAtom[i].type,libAtom[i].element  );
  */


  /*
  Loop2 = CIF->GetLoop ( CIFCAT_LIBBOND );
  if ( !Loop2 ) {
    printf("Error reading bond info from ener_lib.cif\n");
    delete CIF;
    return -1;
  }
  
  nLibBonds = 0;
  RC =0;
  do {
    libBond[nLibBonds] = new CLibBond();
    RC =  libBond[nLibBonds]->GetCif ( Loop2 , nLibBonds );
    if ( !RC) 
      nLibBonds++;
    else
      delete libBond[nLibBonds];
  }while (!RC && nLibBonds < MGSBASE_MAX_LIBBONDS );
  */

  /*
  for ( i = 0; i < nLibBonds; i++) 
    printf ("Bond type %i %s %s %f\n", i,
       libBond[i]->atomType1,libBond[i]->atomType2,
	    libBond[i]->length );
  */


  delete CIF;
  return 0;
}

//------------------------------------------------------------------
int CMGSBase:: CheckCovalentDistance (Element a1, Element a2,realtype d) const {
//------------------------------------------------------------------
  std::string s1 = std::string(a1);
  std::string s2 = std::string(a2);
  s1.erase(std::remove(s1.begin(), s1.end(),' '), s1.end() );
  s2.erase(std::remove(s2.begin(), s2.end(),' '), s2.end() );
  bool found = false;
  for(unsigned i=0;i<covalentDistances.size();i++){
     CMGCovalentDistance cd = covalentDistances[i];
     if((s1==cd.GetFirstAtom()&&s2==cd.GetSecondAtom())||(s2==cd.GetFirstAtom()&&s1==cd.GetSecondAtom())){
       found = true;
       if(d>=cd.GetMinLength()&&d<=cd.GetMaxLength()){
         return 1;
       }
     }
  }
  if(!found) return -1; // We have to trust dictionaries, or else we need complete covalent_distances file
  return 0;
}

//------------------------------------------------------------------
int CMGSBase::LoadCovalentDistances(pstr filename ) {
//------------------------------------------------------------------
  covalentDistances.clear();
  ifstream cdFile(filename);
  char line[1024];
  std::string a1,a2;
  realtype mind,maxd;
  if(cdFile){
    while(cdFile.getline(line,1024)){
      std::istringstream sl(std::string(line),std::istringstream::in);
      sl >> a1;
      sl >> a2;
      sl >> mind;
      sl >> maxd;
      CMGCovalentDistance cd(a1,a2,mind,maxd);
      covalentDistances.push_back(cd);
    }
    cdFile.close();
  }
  return 0;
}

//------------------------------------------------------------------
int CMGSBase::AddCovalentDistance(Element a1, Element a2,realtype mind,realtype maxd) {
//------------------------------------------------------------------
  CMGCovalentDistance cd(std::string(a1),std::string(a2),mind,maxd);
  covalentDistances.push_back(cd);
}

//------------------------------------------------------------------
int CMGSBase::LoadEleLib (pstr filename ) {
//------------------------------------------------------------------
  // Load the monomer link info from CIF file
  int RC;
  PCMMCIFLoop Loop2;
  CIF =  new CMMCIFData();
  RC = CIF->ReadMMCIFData ( filename );
  if ( RC ) {
    printf ("Error reading %s\n",filename);
    delete CIF;
    return -1;
  }

  Loop2 = CIF->GetLoop ( CIFCAT_LIBELEMENT );
  if ( !Loop2 ) {
    printf("Error reading element info from ener_lib.cif\n");
    delete CIF;
    return -1;
  }
  
  nLibElements = 0;
  RC =0;
  do {
    libElement[nLibElements] = new CLibElement();
    RC =  libElement[nLibElements]->GetCif ( Loop2 , nLibElements, this );
    if ( !RC) {
      //libElement[nLibElements]->SetDefaultLibAtom(this);
      nLibElements++;
    }
    else
      delete libElement[nLibElements];
  }while (!RC && nLibElements < MGSBASE_MAX_LIBELEMENTS );


  delete CIF;

  // Old code to derive element bonding from bonds in ener_lib
  // lead to overlong bond rad esp for hydrogen
  //CreateLibElements();
  //for ( i= 0; i < nLibElements; i++) {
  //printf ("Element %s %f\n",libElement[i]->name,libElement[i]->maxBondRad);
  //}
  return 0;
}

//-------------------------------------------------------------
PCLibElement CMGSBase::LibElement ( pstr el ) {
//-------------------------------------------------------------
  int i;
  PCLibElement ret = NULL;

  if (strlen(el)<1) return NULL;

  for (i = 0; i < nLibElements; i++ ) {
    if ( strcmp(libElement[i]->name,el) == 0 ) {
      ret =  libElement[i];
      return ret;
    }
  }
  for (i = 0; i < nLibElements; i++ ) {
    if ( strcmp(libElement[i]->bad_name,el) == 0 ) {
      ret =  libElement[i];
      return ret;
    }
  }
  return NULL;
}

//---------------------------------------------------------------
int  CMGSBase::GetNofLibAtoms() {
//---------------------------------------------------------------
  return nLibAtoms;
}

//--------------------------------------------------------------
//PCLibAtom CMGSBase::LibAtom ( int index ) {
//--------------------------------------------------------------
//  if ( index >= 0 && index < nLibAtoms )
//    return &libAtom[index];
//  else
//    return NULL;
//}




//--------------------------------------------------------------
int CMGSBase::LibAtom ( pstr atomType ) {
//--------------------------------------------------------------
  return LibAtom(atomType,"");
}

//--------------------------------------------------------------
int CMGSBase::LibAtom ( pstr atomType, pstr element ) {
//--------------------------------------------------------------
  int i;
  //cout << "element" << element << "*" << endl;
  // Input string is atomType
  if ( strlen(atomType) >= 1 ) {
    for (i = 0; i < nLibAtoms; i++ ) {
      if ( strcmp(libAtom[i].type,atomType)==0 ) return i;
    }
  }

  if ( strlen(element) >= 1 ) { 

    // Input string is an element name 
    for (i = 0; i < nLibElements; i++ ) {
      if ( strcmp(libElement[i]->name,element)==0 ) 
       return libElement[i]->defaultAtomIndex;
    }

    for (i = 0; i < nLibElements; i++ ) {
      if ( strcmp(libElement[i]->bad_name,element)==0 ) 
        return libElement[i]->defaultAtomIndex;
    }
  }
  
  //cout << "defaulting to C" << endl;
  element = " C";
  for (i = 0; i < nLibElements; i++ ) {
    if ( strcmp(libElement[i]->name,element)==0 ) 
        return libElement[i]->defaultAtomIndex;
  }
  return -1; 
}



//--------------------------------------------------------------
int CMGSBase::LibAtom ( pstr resType, int atomIndex ) {
//--------------------------------------------------------------
  int i;
  PCSBStructure p_sbase_struct;
  PCSBAtom p_sbase_atom;
  pstr atomType;
  //printf ("nLibAtoms %i\n",nLibAtoms);
  LoadedPCSBStructure dummy;

  atomType = "";
  if( (p_sbase_struct = GetStructure(resType,dummy))) {
    if( (p_sbase_atom = p_sbase_struct->Atom[atomIndex]) )
      atomType = p_sbase_atom->energyType;
    //printf ("Got here\n");
  }
  //printf ( "resType %s index %i atomType %s %i\n"
      //   ,resType,atomIndex,atomType,strlen(atomType));
  if ( strlen(atomType) < 1 ) return -1;
  for (i = 0; i < nLibAtoms; i++ ) {
    if ( strcmp(libAtom[i].type,atomType)==0 )
      return  i;
  }
  return -1;
}


//--------------------------------------------------------------
PCLibBond CMGSBase::LibBond ( pstr attype1, pstr attype2 ) {
//--------------------------------------------------------------
  int i;
  if ( strlen(attype1) < 1 || strlen(attype2) < 1) return NULL;
  for (i = 0; i < nLibBonds; i++ ) {
    if ( strcmp(libBond[i]->atomType1,attype1)==0 &&
         strcmp(libBond[i]->atomType2,attype2)==0 )
      return  libBond[i];
  }
  return NULL;
}



//--------------------------------------------------------------
void CMGSBase::CreateLibElements () {
//--------------------------------------------------------------
  bool match;
  int i,j;
  int at1,at2;
  
  nLibElements = 0;
  for (i=0; i<nLibAtoms; i++ ) {
    match = false;
    j = 0;
    while ( (!match) && j < nLibElements ) {
      if ( strcmp(libAtom[i].element,libElement[j]->name ) == 0) {
        match = true;
      }
      j++;
    }
    if ( !match && nLibElements < MGSBASE_MAX_LIBELEMENTS ) {
      libElement[nLibElements] = new CLibElement(libAtom[i].element,i);
      cout << "new libElement  " << libAtom[i].element << " " 
            << libAtom[i].type << endl;
      nLibElements++;
    }
  }

  for ( i =0; i < nLibBonds; i++ ) {
    at1 = LibAtom(libBond[i]->atomType1);
    at2 = LibAtom(libBond[i]->atomType2);
    for ( j = 0; j < nLibElements; j++ ) {
        //printf("j %i name %s\n",j,libElement[j]->name);
      if ( strcmp(libAtom[at1].element,libElement[j]->name) == 0) {
        if (libBond[i]->length > libElement[j]->maxBondRad ) {
            libElement[j]->maxBondRad = libBond[i]->length;
        }
        break;
      }
    }
  }
}


CLibAtom::CLibAtom() {
}
CLibAtom::~CLibAtom() {
}

//---------------------------------------------------------------
int CLibAtom::GetCif( PCMMCIFLoop Loop, int N ) {
//----------------------------------------------------------------
  int RC;
  char hbtype[10];
  Element el;
  if (!Loop) return -1;
  if (N >= Loop->GetLoopLength()) return -1;
   
  RC = CIFGetString ( type, Loop,CIFTAG_LIBATOM_TYPE,N,sizeof(type),"");
  if ( RC || (!type) ) return -1;
  RC =  CIFGetString (hbtype, Loop,CIFTAG_LIBATOM_HBTYPE,N,10,"");
  hbType = encodeHbType(hbtype);
  RC = CIFGetReal1(vdwRadius, Loop,CIFTAG_LIBATOM_VDWRAD, N);
  RC = CIFGetReal1(vdwHRadius, Loop,CIFTAG_LIBATOM_VDWHRAD, N);
  RC = CIFGetReal1( ionRadius,Loop,CIFTAG_LIBATOM_IONRAD, N );
  RC = CIFGetString ( el,Loop,CIFTAG_LIBATOM_ELEMENT,N,sizeof(el),"");
  CLibElement::Justify(el,element);
  RC = CIFGetReal1( charge,Loop,CIFTAG_LIBATOM_CHARGE, N );
  return 0;
}

//------------------------------------------------------------------
int CLibAtom::encodeHbType ( pstr hb ) {
//------------------------------------------------------------------
//Interpret one-letter atom hydrogen bonding type as an integer code
  int i;
  if (strlen(hb) < 1 ) return HBTYPE_UNKNOWN;
  for ( i = 0; i < nHbCodes; i++ ) {
    if ( strcmp ( hb, hbCharCode[i]) == 0 ) 
      return hbCode[i];
  }
  return HBTYPE_UNKNOWN;
}
//------------------------------------------------------------------
char* CLibAtom::getHBType () {
//------------------------------------------------------------------
//Interpret one-letter atom hydrogen bonding type as an integer code
  int i;
  for ( i = 0; i < nHbCodes; i++ ) {
    if ( hbType == hbCode[i] ) return hbCharCode[i];
  }
  return  hbCharCode[0];
}

   
CLibBond::CLibBond() {
}
CLibBond::~CLibBond() {
}

//---------------------------------------------------------------
int CLibBond::GetCif( PCMMCIFLoop Loop, int N ) {
//----------------------------------------------------------------
  int RC;
  char bt[10];
  
  if (!Loop) return -1;
  if (N >= Loop->GetLoopLength()) return -1;
   
  RC = CIFGetString(atomType1,Loop,CIFTAG_LIBBOND_ATOM1,N,sizeof(atomType1),"");
  if ( RC || (!atomType1) ) return -1;
  RC = CIFGetString(atomType2,Loop,CIFTAG_LIBBOND_ATOM2,N,sizeof(atomType2),"");
  //printf ("LibBond %s %s\n",atomType1,atomType2);
  RC = CIFGetString(bt,Loop,CIFTAG_LIBBOND_TYPE,N,10,"");
  bondType=encodeBondType(bt);
  RC = CIFGetReal1( length,Loop, CIFTAG_LIBBOND_LENGTH, N );
  return 0;
}

//-------------------------------------------------------------
int CLibBond::encodeBondType ( pstr ty ) {
//-------------------------------------------------------------
  int i;
  if ( !(ty))  return BONDTYPE_UNKNOWN;
  for ( i = 0; i < nBondCodes; i++) {
    if ( strcmp( ty, bondCharCode[i]) == 0 ) 
      return bondCode[i];      
  }
  return BONDTYPE_UNKNOWN;
}
 
//------------------------------------------------------------
CLibElement::CLibElement ( ) {
//------------------------------------------------------------
  defaultAtomIndex = 0;
}

//------------------------------------------------------------
CLibElement::CLibElement ( pstr el, int atomIndex ) {
//------------------------------------------------------------
  Justify(el,name);
  if ( strlen(el) == 1 ) {
    bad_name[0]=el[0];
    bad_name[1]=' ';
    bad_name[2]='\0';
  } else
  strcpy(bad_name,name);
  defaultAtomIndex = atomIndex;
}

//------------------------------------------------------------
void CLibElement::Justify( pstr el, pstr elo ) {
//--------------------------------------------------------------

  if ( strlen(el) == 1 ) {
    elo[0] =' ';
    elo[1] = el[0];
    elo[2] = '\0';
  }
  else
    strcpy(elo,el);
  
}
//---------------------------------------------------------------
int CLibElement::GetCif( PCMMCIFLoop Loop, int N,  CMGSBase *p_sbase ) {
//----------------------------------------------------------------
  int RC;
  char type [energy_type_len+1];
  
  if (!Loop) return -1;
  if (N >= Loop->GetLoopLength()) return -1;
   
  pstr el = Loop->GetString (CIFTAG_LIBELEMENT_NAME,N,RC);
  if (RC) return RC;
  Justify(el,name);  
  if ( strlen(el) == 1 ) {
    bad_name[0]=el[0];
    bad_name[1]=' ';
    bad_name[2]='\0';
  } else
    strcpy(bad_name,name);

  RC = CIFGetString ( type, Loop,CIFTAG_LIBELEMENT_LIBATOM_TYPE,N,sizeof(type),"");
  if (RC) return RC;
  for (int i=0;i<p_sbase->nLibAtoms;i++) {
    if (strcmp(p_sbase->libAtom[i].type,type)==0) {
      defaultAtomIndex=i;
      break;
    }
  }

  RC = CIFGetReal1(maxBondRad ,Loop, CIFTAG_LIBELEMENT_BONDRAD, N );
  //cout << " CLibElement::GetCif " << name << " " << maxBondRad << " " 
  // << type << " " << defaultAtomIndex << endl;
  return RC;
}

//----------------------------------------------------------------
void  CLibElement::SetDefaultLibAtom ( CMGSBase *p_sbase) {
//----------------------------------------------------------------
  //cout << "SetDefaultLibAtom " << p_sbase->nLibAtoms << " " << name<< endl;
  for ( int i = 0; i<p_sbase->nLibAtoms; i++ ) {
    if (strcmp(p_sbase->libAtom[i].element,name)==0) {
      defaultAtomIndex = i;
      return;
    }
  } 
}

