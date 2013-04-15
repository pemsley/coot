/*
     mmut/mmut_connectivity.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

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


#ifndef __MMUT_Connectivity__
#define __MMUT_Connectivity__


#include <vector>
#include <string>

#include "mgtree.h"
#include "cartesian.h"
#include "connect.h"

#include <mmdb/mmdb_manager.h>
#include "mman_manager.h"

//enum types { SINGLE, DOUBLE, TRIPLE, QUADRUPLE, HBOND };

// Modes of connectivity2 data a - vector between
// atoms and/or points or displacements from the atom or point
 enum connectivity2_types { CONN_POINT_POINT, CONN_ATOM_POINT,
			    CONN_POINT_ATOM,CONN_ATOM_ATOM, 
                            CONN_ATOM_DISP, CONN_POINT_DISP };

class Connection {
 public:
  bool isInTotalSelection;
  std::vector<int> external_connected_atoms;
  std::vector<int> external_connected_spline_atoms;
  std::vector<int> connected_atoms;
  std::vector<unsigned int> orders;
  Connection();
  ~Connection();
  void AddConnection(int serNum);
  void AddExternalConnection(int serNum);
  void AddExternalSplineConnection(int serNum);
};

class Connectivity {
  std::vector<Connection> TotalSelection_new;
  int IsInTotalSelection(int idx) const;
  int nSelAtomsTot;
  int externalBonds;
  PPCAtom SelAtomsTot;
  
 public:
  Connectivity();
  ~Connectivity();
  void AddBonds(CMMANManager *molhnd, int selhnd, const PPCAtom SelAtoms_in, const int nSelAtoms_in, int all_selHnd, int all_selHnd_spline);
  void AddContacts(CMMANManager *molhnd, int selhnd, const PPCAtom SelAtoms_in, const int nSelAtoms_in, const PSContact contacts_in, const int ncontacts_in);
  void AddTrace(PCMMANManager molhnd, const PPCAtom selAtoms, const int nSelAtoms,
		realtype cutoff );
  std::vector<std::vector<int> > GetConnectivityLists(void) const;
  std::vector<std::vector<int> > GetExternalConnectivityLists(void) const;
  std::vector<std::vector<int> > GetExternalSplineConnectivityLists(void) const;
  int GetNumberOfAtoms(void) const;
  int GetNumberOfExternalBonds(void) { return externalBonds; }
  PPCAtom GetAtoms(void) const;
  const std::vector<int> GetCAtomIndex(void) const;
  void Clear(void);
  //Tree GetTree(PCMMANManager molhnd) const;
};

// Connectivy2 class to handle 'connectivity' such as Hbonds which is 
// optionally between two separate molecules (ie. MMDBManagers)


class Connectivity2 {
  
  int data_mode;
  int crystal_axes;
  int tagged;
  bool selection_on;
  int nAtomSets;
  int nxyzSets;


 public:
  std::vector<PCAtom> pAtom1;
  std::vector<PCAtom> pAtom2;
  std::vector<Cartesian> XYZ1;
  std::vector<SimpleConnection> connected;
  std::vector<int> tags;
  std::vector<bool> selected;

  void DeleteConnections(std::vector<unsigned int> indices, bool all=false);

  Connectivity2(int data_mode=CONN_ATOM_ATOM,int crystal_axes=0,int tagged=0);
  ~Connectivity2();
  void AddConnection(PCAtom p_atom1,PCAtom p_atom2, const std::string &label_in="", int tag=-1);
  void AddConnection( double xyz1[3] ,  double xyz2[3] , const std::string &label_in,int tag=-1);
  void AddConnection(  PCAtom p_atom1, double xyz2[3] ,
    const std::string &label_in, int tag=-1);
  void InsertConnection(bool replace, unsigned int position,PCAtom p_atom1,PCAtom p_atom2, const std::string &label_in="", int tag=-1);
  void InsertConnection( bool replace, unsigned int position, double xyz1[3] ,  double xyz2[3] , const std::string &label_in,int tag=-1);
  void InsertConnection(  bool replace, unsigned int position, PCAtom p_atom1, double xyz2[3] ,
    const std::string &label_in, int tag=-1);
  void AddUniqueConnection(PCAtom p_atom1,PCAtom p_atom2, const std::string &label_in="",int tag = -1);
  void UpdateCoordinates(bool label_dist=false);
  void RemoveConnection(PCAtom p_atom1,PCAtom p_atom);
  void RemoveConnection(PCAtom p_atom1,int position);
  void DeleteConnection(int index);
  void DeleteConnections();
  int DeleteTaggedConnections(int first, int last=-1);
  int GetNofConnections() { return  connected.size(); }
  std::vector<unsigned int> FindConnections(PCAtom p_atom1,PCAtom p_atom,bool switchpos=1);
  std::vector<unsigned int> FindConnections(PCAtom p_atom1, int position=0);
  int FindNofConnections(PCAtom p_atom1,int position=0);
  void Clear(void);	
  int SelectVectors (int nSelTags , int *selTags );
  double* Extent();
  std::string Print( const std::string &tag1="", const std::string &tag2="" );
  int DataMode(void) { return data_mode; }
  std::string GetLabel(int i){ return connected[i].label ; }
  Cartesian GetCoordinate(int i, int j=1);
  PCAtom GetAtom(int i, int j=1);
  std::string GetAtomID ( int i, int j);
  int GetTag(int iV);
  int SetTag(int iV,int i);
  double GetRMSD();
  int AddContactFromSelHandle( PCMMANManager molHnd1,int selHnd1,PCMMANManager molHnd2,int selHnd2,int tag=-1);
  int AddContacts(PCMMANManager molHnd1,int selHnd1, PCMMANManager molHnd2,int selHnd2_in,realtype  dist1, realtype  dist2, int  seqDist, int inter_model=0, int closest_bonding=5, int handle_hbond=0 );
  int AddRangeConnections(PCResidue res1,PCResidue res2, 
	                  PCResidue mres1,  PCResidue mres2,
		          const std::vector<std::string>& mainchain_name,int tag );
  int AddCloseRangeConnections(int set,PCResidue res1,
                    PCResidue res2,PCMMANManager M2,
		      double central_cutoff, double cutoff,
		    const std::vector<std::string>& mainchain_name,
                    const std::string &centralAtom, int tag );

  int AddRangeWithSameIdConnections( PCMMANManager molHnd1,int selHnd1, PCMMANManager molHnd2,int selHnd2, 
				     const std::vector<std::string>&  mainchain_name, int tag );

  int AddCloseAtoms( PCMMANManager molHnd1,int selHnd1, 
                     PCMMANManager molHnd2,int selHnd2 , 
		     realtype central_cutoff, realtype cutoff , char *central, int tag );

  int GetSelection (int mode=1 );
  int Superpose(int fixed);
  int MatchGraphs(PCResidue pRes1,const pstr altLoc1,
                  PCResidue pRes2,const pstr altLoc2,
		  int Hflag,int tag, float fracMinMatch=-1.0,
                  bool keep_matches=true);
};

inline int Connectivity::IsInTotalSelection(int idx) const {
  if(!TotalSelection_new[idx].isInTotalSelection)
    return 0;
  else
    return 1;
}



#endif
