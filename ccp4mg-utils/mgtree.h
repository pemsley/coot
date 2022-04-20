/*
     util/mgtree.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_TREE_
#define _CCP4MG_TREE_
#include <vector>
#include "cartesian.h"
#include <utility>
#include "geomutil.h"

class bond_pair_cmp{
  public:
    int operator()(const std::pair<int,int> &p1, const std::pair<int,int> &p2) const
       { return (p1.first < p2.first)&&(p1.second < p2.second) ; }
};

class TreeVertex{
  private:
   int id;
   int parent_id;
   double parent_dist;
   double parent_bond_angle;
   double parent_dihedral_angle;
   TreeVertex* parent;
   std::vector<TreeVertex*> children;
   /* 
      Since this corresponds to something outside this tree, the actual child doesn't exist. 
      So this number corresponds to something beyond the scope of this object, ug. The application
      has to know what this means ....
   */
   std::vector<Cartesian> ext_children; 
   Cartesian coord;
  public:
   TreeVertex(); 
   ~TreeVertex(); 
   void SetParentID(const int &parent_id_in){parent_id = parent_id_in;};
   void SetID(const int &id_in){id = id_in;};
   void SetCoord(const Cartesian &coord_in){coord = coord_in;};
   void SetDummy(Cartesian &dummy_in, Cartesian &dummy2_in){Dummy = dummy_in;Dummy2 = dummy2_in;};
   void SetAngles(void);
   Cartesian Dummy;
   Cartesian Dummy2;
   Cartesian GetCoord(void) const {return coord;};

   void SetParentDistance(double parent_dist_in){parent_dist = parent_dist_in;};
   void SetParentBondAngle(double parent_bond_angle_in){parent_bond_angle = parent_bond_angle_in;};
   void SetParentDihedralAngle(double parent_dihedral_angle_in){parent_dihedral_angle = parent_dihedral_angle_in;};

   void AddParentDistance(double parent_dist_in){parent_dist += parent_dist_in;};
   void AddParentBondAngle(double parent_bond_angle_in){parent_bond_angle += parent_bond_angle_in;};
   void AddParentDihedralAngle(double parent_dihedral_angle_in){parent_dihedral_angle += parent_dihedral_angle_in;};

   void SetParent(TreeVertex* parent_in){parent = parent_in;};
   void AddChild(TreeVertex *child){children.push_back(child);};
   void AddExternalChild(const Cartesian &child){ext_children.push_back(child);}; // Add a child connection index which is outside this tree.
   int GetParentID() const {return parent_id;};
   int GetID() const {return id;};
   TreeVertex* GetParent() const {return parent;};
   std::vector<TreeVertex*> GetChildren() const {return children;};
   TreeVertex* GetChild(int i) const {return children[i];};
   Cartesian GetExternalChild(int i) const {return ext_children[i];};
   int GetNumberOfChildren() const {return children.size();};
   int GetNumberOfExternalChildren() const {return ext_children.size();};
   void PrintTree(void) const; // {std::cout << *this;};
   friend std::ostream& operator<<(std::ostream &c, TreeVertex a);
   int FindDepth(void) const ;
   double GetParentDistance() const {return parent_dist;};
   double GetParentBondAngle() const {return parent_bond_angle;};
   double GetParentDihedralAngle() const {return parent_dihedral_angle;};
   std::vector <TreeVertex*> GetBranch();
   int GetNumberOfDescendants() const ;
   void GetDescendants(std::vector<TreeVertex *> &desc_vertices, const TreeVertex *stop_node=0) const ;
   void GetDescendants(std::vector<TreeVertex *> &desc_vertices, const std::vector<TreeVertex*> stop_nodes) const ;
   void GetNonDescendants(std::vector<TreeVertex *> &desc_vertices) const ;
   void GetNonDescendants(std::vector<TreeVertex *> &desc_vertices, const std::vector<TreeVertex*> stop_nodes) const ;
};

class Tree{

  private:
   std::vector<TreeVertex*> coords;
   std::vector<int> scanned;
   std::vector<std::vector<int> > connectivity;
   void CalculateTree(void);
   void RecurseCalculateTree(TreeVertex* coord);
   void RecurseCalculateTreeWithLengthsAndAngles(TreeVertex* coord);
   void SetDummy(TreeVertex *coord);
   void ExtendBranchCartesians(const Cartesian &p1, const Cartesian &p2, const Cartesian &p3, TreeVertex* child, std::vector<Cartesian> &cartesians) const ;
   void RecurseZMatrix(std::ostream &c, const TreeVertex *vertex, const std::vector<std::string> &labels, const std::string &separater);
   int start;
   std::vector<int> permutation;
   void ClearCoords();
  public:
   Tree& operator=(const Tree &a);
   Tree(const Tree &t);
   Tree();
   // This method and following constructor build a tree from a set of existing bonds, angles and torsions (eg. RefMac monomer library).
   void SetBondsAnglesTorsions(const int &nvertices,
                               const std::vector<std::vector<int> > &bonds, const std::vector<double> &bond_lengths,
                               const std::vector<std::vector<int> > &angles, const std::vector<double> &bond_angles,
                               const std::vector<std::vector<int> > &torsions, const std::vector<double> &torsion_angles,
			       const std::vector<std::vector<int> > &chirals);
   Tree(const int &nvertices,
        const std::vector<std::vector<int> > &bonds, const std::vector<double> &bond_lengths,
        const std::vector<std::vector<int> > &angles, const std::vector<double> &bond_angles,
        const std::vector<std::vector<int> > &torsions, const std::vector<double> &torsion_angles,
        const std::vector<std::vector<int> > &chirals);
   Tree(const std::vector<Cartesian> &SelAtoms, int start, const std::vector<std::vector<int> > &conn_lists, const  std::vector<std::vector<Cartesian> > &ext_cartesians);
   Tree(const std::vector<Cartesian> &SelAtoms, int start, const std::vector<std::vector<int> > &conn_lists, const  std::vector<std::vector<Cartesian> > &ext_cartesians, const std::vector<std::vector<int> > &forced_connections);
   ~Tree();
   void SetCoords(const std::vector<Cartesian> &SelAtoms, int start, const std::vector<std::vector<int> > &conn_lists);
   void SetCoords(const std::vector<Cartesian> &SelAtoms, int start, const std::vector<std::vector<int> > &conn_lists, const std::vector<std::vector<int> > &forced_connections);
   void SetCoords(const std::vector<Cartesian> &SelAtoms, int start, const std::vector<std::vector<int> > &conn_lists, const  std::vector<std::vector<Cartesian> > &ext_cartesians, const std::vector<std::vector<int> > &forced_connections);
   friend std::ostream& operator<<(std::ostream &c, Tree a);
   void Print(void); // { std::cout << *this;};
   int FindMaxDepth(void);
   std::vector<std::pair<int,int> > extra_bonded_pairs;
   void PrintZMatrix(std::ostream &c, const std::vector<std::string> &labels, const std::string &separater=",");
   void PrintZMatrix(const std::vector<std::string> &labels, const std::string &separater=",");
   int GetNumberOfVertices() const {return coords.size();}
   Cartesian GetCartesian(TreeVertex *unknown) const ;

   /*
    * When start>0, there is a permutation applied to 
    * coords and connections. This allows different trees to be built
    * from same set of atoms. Essential for things like keeping different
    * parts of structure fixed when doing rotations.

    * These methods by default apply the permutation or its inverse where
    * appropriate so that you do not need to know that the nodes have been
    * permuted. So if you create from a set of coordinates and then request
    * rotation about bond 2,3, the rotation will be applied about 2,3 of the
    * original order. What you want usually. One can ask for the result to
    * be in the permuted space by passing optional argument. Will probably never 
    * be required.
    */

   TreeVertex* GetCoord(int i, bool permuted=false) const ;
   std::vector <TreeVertex*> GetCoords(bool permuted=false) const;
   Cartesian GetCartesian(int i, bool permuted=false) const;
   std::vector<Cartesian> GetAllCartesians(bool permuted=false) const ;
   void RotateAboutBond(int atom, int child, double TorsionAngleDiff, bool permuted=false);
   void SetDihedralAngle(int atom, int child, double TorsionAngle, bool permuted=false, int movingAtom=-1, int baseAtom=-1);
   std::vector<std::vector<int> > FindLongBranches(int req_depth);
   void AddVertex(int pid, double dist, int gpid, double angle, int ggpid, double dihedral, int chiral);
   void SetDihedralAngle(int baseAtom, int atom, int child, int movingAtom, double TorsionAngle);
   void ForceEarlyConnection(int parent,int child);
   
};

#endif
