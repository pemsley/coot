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

class CMMUTManager : public mmdb::Manager  {

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
    //void PrintLengthsAndAngles(mmdb::realtype min, mmdb::realtype max);
 
    int NumberOfHydrogens(int selHnd);
    int ResNoLookup(pstr resname);
    int TotalNumRes(int selHnd);
    pstr GetSequence(int selHnd);
    int *AtomicComposition(int selHnd);
    int *ResidueComposition(int selHnd);
    mmdb::PResidue NextResidue( Pmmdb::Residue pRes ,int increment=1);

    mmdb::realtype  BondLength(mmdb::PAtom A, Pmmdb::Atom B);
    mmdb::realtype  BondAngle(mmdb::PAtom A, Pmmdb::Atom B, Pmmdb::Atom C);
    mmdb::realtype  TorsionAngle(mmdb::PAtom A, Pmmdb::Atom B, Pmmdb::Atom C, Pmmdb::Atom D);
    mmdb::realtype  MolWeight(int selHnd);
    mmdb::realtype  MolWeightWithH(int selHnd);
    mmdb::realtype *GetBValues(int selHnd);
    mmdb::realtype *CentreOfMass(int selHnd);
    mmdb::realtype *Extent(int selHnd);

    Boolean isMainChain(mmdb::PAtom p_atom);
    Boolean doAltLocMatch ( mmdb::PAtom pa1, Pmmdb::Atom pa2 ); 
    int NameComparison (const char *name , int ntypes , const char *types[] );
    std::string  TrimString(pstr inp);
    std::string AtomLabel(mmdb::PAtom p_atom, int mask[]);
    Boolean ChainIDisDigit(Pmmdb::Chain p_ch);
    const char* AtomLabel_atom1(mmdb::PAtom p_atom);
    const char* AtomLabel_atom(mmdb::PAtom p_atom);
    const char* AtomLabel_residue(mmdb::PAtom p_atom);
    const char* AtomLabel_chain(mmdb::PAtom p_atom);
    const char* AtomLabel_residue1(mmdb::PResidue p_res);
    const char* AtomLabel_mask(mmdb::PAtom p_atom, int mask[]);
    //int ApplyTransform(int selHnd,double rotmat[],double transv[]);

    //Editor
    int WriteSelection (int selHnd,char *file, char *format="PDB");
    int PutSelectedAtoms (int selHnd , const mmdb::PManager mmdb2);
    int CopySelection (int selHnd,const mmdb::PManager mmdb2);
    int FindCloseAtomPairs ( int selHnd, double min_distance, 
			     double max_distance);
   private: 

     //Analysis
    int *iatom_types;
    int *iatom_type_lookup;

 
};

#endif
