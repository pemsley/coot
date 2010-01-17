/*
     mmut/mmut_manager.h: CCP4MG Molecular Graphics Program
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


#ifndef __MMUT_Manager__
#define __MMUT_Manager__

#include <string>
#include <mmdb_manager.h> 

#define PI 3.141592653589793238462643

// =======================  CMMUTManager  ===========================


DefineClass(CMMUTManager);
DefineStreamFunctions(CMMUTManager) ;

class CMMUTManager : public CMMDBManager  {

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
    int *AtomicComposition(int selHnd);
    int *ResidueComposition(int selHnd);
    PCResidue NextResidue( PCResidue pRes ,int increment=1);

    realtype  BondLength(PCAtom A, PCAtom B);
    realtype  BondAngle(PCAtom A, PCAtom B, PCAtom C);
    realtype  TorsionAngle(PCAtom A, PCAtom B, PCAtom C, PCAtom D);
    realtype  MolWeight(int selHnd);
    realtype  MolWeightWithH(int selHnd);
    realtype *GetBValues(int selHnd);
    realtype *CentreOfMass(int selHnd);
    realtype  Mass(int selHnd);
    realtype *CentreOfCoordinates(int selHnd);
    int  NumberOfAtoms(int selHnd);
    realtype *Extent(int selHnd);

    Boolean isMainChain(PCAtom p_atom);
    Boolean doAltLocMatch ( PCAtom pa1, PCAtom pa2 ); 
    int NameComparison ( char *name , int ntypes , char *types[] );
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
    int WriteSelection (int selHnd,char *file, char *format="PDB");
    int PutSelectedAtoms (int selHnd , const PCMMDBManager mmdb2);
    int CopySelection (int selHnd,const PCMMDBManager mmdb2);
    int FindCloseAtomPairs ( int selHnd, double min_distance, 
			     double max_distance);
   private: 

     //Analysis
    int *iatom_types;
    int *iatom_type_lookup;

 
};

#endif
