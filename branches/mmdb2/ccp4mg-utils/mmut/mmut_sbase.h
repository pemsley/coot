/*
     mmut/mmut_sbase.h: CCP4MG Molecular Graphics Program
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


#ifndef __CCP4Sbase__ 
#define __CCP4Sbase__
#include <mmdb_manager.h>
#include <mmut_manager.h>
#include <mmdb_sbase.h>
#include <string>
#include <vector>
#include <map>
#include <mgtree.h>

// The CIF tags for _chem_comp_atom
// and _chem_comp_bond as used in the LibCheck
// monomer library
#define CIFCAT_COMP "_chem_comp"
#define CIFTAG_COMP_ID "id"
#define CIFTAG_COMP_CODE "three_letter_code"
#define CIFTAG_COMP_GROUP "group"
#define CIFCAT_COMP_ATOM "_chem_comp_atom"
#define CIFTAG_COMP_ATOM_ID "atom_id"
#define CIFTAG_COMP_ALT_ATOM_ID "alt_atom_id"
#define CIFTAG_COMP_TYPE_SYMBOL "type_symbol"
#define CIFTAG_COMP_TYPE_ENERGY "type_energy"
#define CIFTAG_COMP_CHARGE "partial_charge"
#define CIFTAG_COMP_X "x" 
#define CIFTAG_COMP_Y "y" 
#define CIFTAG_COMP_Z "z" 
#define CIFCAT_COMP_BOND "_chem_comp_bond"
#define CIFTAG_COMP_BOND_ATOM_1 "atom_id_1"
#define CIFTAG_COMP_BOND_ATOM_2 "atom_id_2"

// Synonyms
#define CIFCAT_COMP_SYNONYM "_chem_comp_synonym"
#define CIFTAG_COMP_SYNONYM_ID "comp_id"
#define CIFTAG_COMP_SYNONYM_ALTERNATIVE "comp_alternative_id"

// The CIF tags from the mon_lib.cif file
// for the 'links' tables
#define CIFCAT_LINK "_chem_link"
#define CIFTAG_LINK_ID "id"
#define CIFTAG_LINK_NAME "name"
#define CIFTAG_LINK_COMP_ID_1 "comp_id_1"
#define CIFTAG_LINK_MOD_ID_1 "mod_id_1"
#define CIFTAG_LINK_GROUP_COMP_1 "group_comp_1"
#define CIFTAG_LINK_COMP_ID_2 "comp_id_2"
#define CIFTAG_LINK_MOD_ID_2 "mod_id_2"
#define CIFTAG_LINK_GROUP_COMP_2 "group_comp_2"

#define CIFCAT_LINK_BOND "_chem_link_bond"
#define CIFTAG_LINK_BOND_ATOM1 "atom_id_1"
#define CIFTAG_LINK_BOND_ATOM2 "atom_id_2"

// CIF tags for ener_lib.cif file the lib_atom table
#define CIFCAT_LIBATOM "_lib_atom"
#define CIFTAG_LIBATOM_TYPE "type"
#define CIFTAG_LIBATOM_WEIGHT "weight"
#define CIFTAG_LIBATOM_HBTYPE "hb_type"
#define CIFTAG_LIBATOM_VDWRAD "vdw_radius"
#define CIFTAG_LIBATOM_VDWHRAD "vdwh_radius"
#define CIFTAG_LIBATOM_IONRAD "ion_radius"
#define CIFTAG_LIBATOM_ELEMENT "element"
#define CIFTAG_LIBATOM_VALENCY "valency"
#define CIFTAG_LIBATOM_CHARGE "surface_potential_charge"

// CIF tags for  ener_lib.cif lib_element table
// (this is added by Lizp to carry bonding radii for
// bonding by contact search when all else fails)
#define CIFCAT_LIBELEMENT "_lib_element"
#define CIFTAG_LIBELEMENT_NAME "name"
#define CIFTAG_LIBELEMENT_BONDRAD "bond_radius"
#define CIFTAG_LIBELEMENT_LIBATOM_TYPE "lib_atom_type"

//CIF tags for the ener_lib.cif lib_bond table

#define CIFCAT_LIBBOND "_lib_bond"
#define CIFTAG_LIBBOND_ATOM1 "atom_type_1"
#define CIFTAG_LIBBOND_ATOM2 "atom_type_2"
#define CIFTAG_LIBBOND_TYPE "type"
#define CIFTAG_LIBBOND_CONST "const"
#define CIFTAG_LIBBOND_LENGTH "length"
#define CIFTAG_LIBBOND_LENGTH_ESD "value_esd"

#define MGSBASE_MAX_LIBELEMENTS 200
#define MGSBASE_MAX_LIBATOMS 300
#define MGSBASE_MAX_LIBBONDS 800


enum { RESTYPE_PEPTIDE, RESTYPE_DPEPTIDE, RESTYPE_LPEPTIDE,
       RESTYPE_NUCL,RESTYPE_DNA, RESTYPE_RNA,
       RESTYPE_SACH, RESTYPE_DSACH, RESTYPE_LSACH,
       RESTYPE_SOLVENT, RESTYPE_SOLUTE,RESTYPE_NONPOLY, RESTYPE_PSEUDO,
       RESTYPE_UNKNOWN };

// Atom hydrogen bonding type - following definitions in ener_lib.cif
enum { HBTYPE_UNKNOWN, HBTYPE_NEITHER, HBTYPE_HYDROGEN, HBTYPE_DONOR,
       HBTYPE_BOTH, HBTYPE_ACCEPTOR };

// Bond types
enum { BONDTYPE_UNKNOWN, BONDTYPE_SINGLE, BONDTYPE_DOUBLE, BONDTYPE_TRIPLE,
       BONDTYPE_AROMATIC, BONDTYPE_METAL, BONDTYPE_DELOC };

typedef std::map<std::string,PCSBStructure> LoadedPCSBStructure;
typedef std::map<std::string,PCSBStructure>::iterator LoadedPCSBStructure_iter;

DefineClass ( CCompundGroup);

class CCompoundGroup {
  friend class MGCLink;
  friend class MGCLinkGroup;
 public:
  CCompoundGroup();
  void Set ( int cd);
  void Set ( pstr name );
  static int GetCifGroupCode ( pstr name);
  bool Match ( int cd );
  static bool groupMatch[13][13];
  static const char *cifGroupNames[12];
  static int groupCode[12];

 protected:
  int code;

};


DefineClass ( MGCLinkGroup);

class MGCLinkGroup {
friend class MGCLink;

 public:  
  MGCLinkGroup ();
  void Set ( pstr comp, pstr modif, pstr grp, cpstr atm );
  ~MGCLinkGroup();
  void Print();
  bool Match( int grp, pstr comp, pstr atm );
 

 protected:
  CompoundID compId;
  CompoundID modId;
  CCompoundGroup group;
  AtomName atom;
  
};

DefineClass ( MGCLink );
DefineClass(CMGSBase);

class MGCLink {
  friend class CMGSBase;
  friend class CMolBonds;
 public:
  MGCLink();
  ~MGCLink();
  int GetCif(PCMMCIFLoop Loop, int N);
  void Print();
  
  //protected:
  CompoundID id;
  MGCLinkGroup lg1;
  MGCLinkGroup lg2;
  int GetCifBond ( PCMMCIFData dataBlock );
 
  //int cifGroupCode( pstr cifGroup);
 
};

DefineClass (CLibElement);

class CLibElement {
  friend class CMGSBase;
 public:
  CLibElement();
  CLibElement ( pstr el, int atomIndex );
  static void Justify(pstr el, pstr elo);
  int GetCif (PCMMCIFLoop Loop, int N, CMGSBase *p_sbase);
  void SetDefaultLibAtom( CMGSBase *p_sbase);
  Element name;
  Element bad_name;
  int defaultAtomIndex;  //default LibAtom index
  realtype maxBondRad;

};

DefineClass (CLibAtom);

class CLibAtom {

 public:
  CLibAtom();
  ~CLibAtom();
  int GetCif(PCMMCIFLoop Loop, int N);

  char type [energy_type_len+1];
  int hbType;
  realtype vdwRadius;
  realtype vdwHRadius;
  realtype ionRadius;
  Element element;
  realtype charge;

  static int nHbCodes;
  static const char *hbCharCode[6];
  static int hbCode[6];
  int encodeHbType( pstr hb);
  const char* getHBType();

};


DefineClass ( CLibBond );

class CLibBond {
  
 public:
  CLibBond();
  ~CLibBond();
  int GetCif (PCMMCIFLoop Loop, int N);
  
  char atomType1[energy_type_len+1];
  char atomType2[energy_type_len+1];
  int bondType;
  //realtype const;
  realtype length;
  //realtype length_esd;

  static int nBondCodes;
  static const char *bondCharCode[6];
  static int bondCode[6];
  static int encodeBondType( pstr ty );

};

class CMGSBase {
  friend class CMolBonds;
  std::vector<std::vector<int> > chirals;
public:
  CMGSBase(char *mon_dir ,char *user_mon_dir,char *sb, char *ener_lib, char *mon_lib, char *ele_lib );
  ~CMGSBase();
  int LoadMonLib ( pstr filename );
  int LoadEnerLib ( pstr filename );
  int LoadEleLib ( pstr filename );
  int LoadSynonyms (pstr filename );
  std::string  ListMonomer(char *mon, bool unremediated=false);
  void InitialiseErrorReporting () { reported_errors.clear(); }
  std::string AssignAtomType ( PCResidue pRes,
      LoadedPCSBStructure monlib,
      std::map<std::string,std::string> &customResSynonym,
      int udd_sbaseCompoundID,int udd_sbaseAtomOrdinal, 
      int udd_atomEnergyType, const bool unremediated= false );
  std::string ListAtomType ( PCMMUTManager molHnd, PCResidue pRes,
  int udd_sbaseCompoundID, int udd_sbaseAtomOrdinal, int udd_atomEnergyType );
  //int GraphSearch ( PCResidue pRes, PCSBStructure &pSbaseRes,
  //     int &nAtom, ivector &nMatchAtom, imatrix &matchAtom  );

  PCSBStructure GetStructure ( const ResName resNam , LoadedPCSBStructure monlib, const bool unremediated=false);
  int LoadMonomerLibrary( char* filename, LoadedPCSBStructure &monlib);
  PCSBStructure LoadCifMonomer ( const ResName resNam , const PCMMCIFFile file, const bool unscramble=true );
  int MatchGraphs(PCResidue pRes,int Hflag, Boolean Cflag, const pstr altLoc, 
		  PCSBStructure pSbaseRes, int &nMatched,
		  ivector match, int minMatchSize );
  //PCLibAtom LibAtom (char *);
  int LibAtom(char*);
  int LibAtom(char *, char *);
  int LibAtom ( pstr resType, int atomIndex );
  int GetNofLibAtoms();
  //PCLibAtom LibAtom ( int index );
  PCLibElement LibElement ( char *);
  PCLibBond LibBond ( char *, char *);
  int maxAtomInRes;
  Tree GetMonomerLibraryTree(const char *monomer_name);
  PPCAtom GetMonomerLibraryStructure(const char *monomer_name);

  //private:

  char monomers_dir[500];
  char user_monomers_dir[500];
  int InitSBase(char *sb);
  void CreateLibElements();
  //int MatchGraphs  (PCGraph G1, PCGraph G2);
    
  int graphSearchMode;
  static PCSBase SBase;
  //static int maxNLibElements;
  //static int maxNLibAtoms;
  //static int maxNLibBonds;
 
  LoadedPCSBStructure loadedPStruct;
  LoadedPCSBStructure unremedPStruct;
  std::map<std::string,std::string> synonyms;
  // Record residuesfor which error has been reported
  // This needs to be reinitiallised for each new loaded structure
  std::map<std::string,int> reported_errors;
  
  PCMMCIFData CIF;
  int nLinks;
  // Ooch again!
  MGCLink *link[100];

  int nLibElements;
  CLibElement *libElement[MGSBASE_MAX_LIBELEMENTS];

  int nLibAtoms;
  std::vector<CLibAtom> libAtom;

  int nLibBonds;
  CLibBond *libBond[MGSBASE_MAX_LIBBONDS];

  realtype fracMatch;
};


#endif
