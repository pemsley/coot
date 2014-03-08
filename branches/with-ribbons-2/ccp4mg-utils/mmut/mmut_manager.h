/*
     mmut/mmut_manager.h: CCP4MG Molecular Graphics Program
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


#ifndef __MMUT_Manager__
#define __MMUT_Manager__

#include <string>
#include <vector>
#include <mmdb/mmdb_manager.h>
#include "cartesian.h"

#define PI 3.141592653589793238462643

// =======================  CMMUTManager  ===========================


DefineClass(CMMUTManager);
DefineStreamFunctions(CMMUTManager) ;

class CMMUTManager : public CMMDBManager  {

    std::vector<matrix> biomatrices;
    bool doneBiomolecule;
  public :

    CMMUTManager();
    ~CMMUTManager();

    const CMMDBCryst& get_cell() { return Cryst; } 
    CMMDBCryst* get_cell_p() { return & Cryst; } 
  
    void ListAtomInfo(int selHnd);
    void PrintAtomicComposition(int selHnd);
    void PrintResidueComposition(int selHnd);
    void PrintSequence(int selHnd);
    void PrintBValues(int selHnd);
    //void PrintLengthsAndAngles(int selHnd);
    //void PrintLengthsAndAngles(realtype min, realtype max);
 
    int NumberOfHydrogens(int selHnd);
    int ResNoLookup(pstr resname);
    int TotalNumRes(int selHnd);
    pstr GetSequence(int selHnd);
    pstr GetSequenceFromResidues(int selHnd);
    void SelectOneAtomForNonAminoNucleotide(int selHnd, int atomSelHnd);
    int *AtomicComposition(int selHnd);
    int *ResidueComposition(int selHnd);
    PCResidue NextResidue( PCResidue pRes ,int increment=1);

    realtype  BondLength(PCAtom A, PCAtom B);
    realtype  BondAngle(PCAtom A, PCAtom B, PCAtom C);
    realtype  TorsionAngle(PCAtom A, PCAtom B, PCAtom C, PCAtom D);
    realtype  MolWeight(int selHnd);
    realtype  MolWeightWithH(int selHnd);
    realtype *GetBValues(int selHnd);
    std::vector<double> GetBValuesDoubleVector(int selHnd);
    realtype *CentreOfMass(int selHnd);
    Cartesian CentreOfMassAsCartesian(int selHnd);
    realtype  Mass(int selHnd);
    realtype *CentreOfCoordinates(int selHnd);
    Cartesian CentreOfCoordinatesAsCartesian(int selHnd);
    int  NumberOfAtoms(int selHnd);
    realtype *Extent(int selHnd);
    realtype ExtentSize(int selHnd);
    std::vector<Cartesian> GetPrincipalComponents(int selHnd);

    Boolean isMainChain(PCAtom p_atom);
    Boolean doAltLocMatch ( PCAtom pa1, PCAtom pa2 ); 
    int NameComparison ( const char *name , int ntypes , const char *types[] );
    std::string  TrimString(pstr inp);
    std::string AtomLabel(PCAtom p_atom, int mask[]);
    Boolean ChainIDisDigit(PCChain p_ch);
    const char* AtomLabel_atom1(PCAtom p_atom);
    const char* AtomLabel_atom(PCAtom p_atom);
    const char* AtomLabel_residue(PCAtom p_atom);
    const char* AtomLabel_chain(PCAtom p_atom);
    const char* AtomLabel_residue1(PCResidue p_res);
    const char* AtomLabel_mask(PCAtom p_atom, int mask[]);
    //int ApplyTransform(int selHnd,double rotmat[],double transv[]);

    //Editor
    int WriteSelection (int selHnd,char *file, const char *format="PDB");
    int PutSelectedAtoms (int selHnd , const PCMMDBManager mmdb2);
    int CopySelection (int selHnd,const PCMMDBManager mmdb2);
    int FindCloseAtomPairs ( int selHnd, double min_distance, 
			     double max_distance);
   int FixElementNames();

   std::string Source();
   std::string Unknowns();
   std::string GetRemarksString();
   std::string SiteInfo();
   double Resolution();
   std::string StructureTitle();

   int ApplyPDBSecStructure(int model);

   std::vector<double> GetCellInfo();
   std::string MMUTGetSpaceGroup();

   std::string SelectionToSCOP(int selHnd);

   std::vector<matrix> GetBiomoleculeAsMatrices(int nBiomol,int nModel=1);
   
   static CMMDBManager* GetCAModel(CMMDBManager *molHnd);
   int GenerateTransformedChain(const char *chainID, realtype *vmat, CMMDBManager *molHnd2);

   int RemoveSmallHelices(int model,int minHelices=4);

   void  SelectAminoNotHet (
             int   selHnd,   // must be obtained from NewSelection()
             int   selType,  // selection type STYPE_XXXXX
             int   iModel,   // model number; iModel=0 means
                             // 'any model'
             cpstr Chains,   // may be several chains "A,B,W"; "*"
                             // means 'any chain' (in selected
                             // model(s))
             int   ResNo1,   // starting residue sequence number
             cpstr Ins1,     // starting residue insertion code; "*"
                             // means 'any code'
             int   ResNo2,   // ending residue sequence number.
             cpstr Ins2,     // ending residue insertion code; "*"
                             // means 'any code'. Combination of
                             // ResNo1=ResNo2=ANY_RES and
                             // Ins1=Ins2="*" means 'any residue'
                             // (in selected chain(s))
             cpstr RNames,   // may be several residue names
                             // "ALA,GLU,CIS"; "*" means
                             // 'any residue name'
             cpstr ANames,   // may be several names "CA,CB"; "*"
                             // means 'any atom' (in selected
                             // residue(s))
             cpstr Elements, // may be several element types
                             // 'H,C,O,CU'; "*" means 'any element'
             cpstr altLocs,  // may be several alternative
                             // locations 'A,B'; "*" means
                             // 'any alternative location'
             int selKey=SKEY_OR // selection key
           );
   bool isNTerminusBound(PCResidue res);
   bool isPeptideBound(PCResidue res);

   private: 

     //Analysis
    int *iatom_types;
    int *iatom_type_lookup;

 
};

#endif
