/* high-res/test-atom-graph.cc
 * 
 * Copyright 2003, 2004  by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "coot-atom-graph.hh"

void simple_nodes();
void test_alignment(int argc, char**argv);

int
main(int argc, char **argv) {

   simple_nodes();

   test_alignment(argc, argv);

   return 0;
}


void simple_nodes() { 

   std::vector<std::vector<coot::node_info> > connection_indices;
   std::vector<clipper::Coord_orth> coords;

   int n_nodes = 13;
   connection_indices.resize(n_nodes);
   coords.resize(n_nodes, clipper::Coord_orth(0.0, 0.0, 0.0));

//    
//     1           10
//    / \           |
//   0---2         11
//   |    \         |
//   5     3-4     12
//   |
//   6-7
//  /   \    . shut up compiler
// 8     9
//
   
   std::vector<coot::node_info> a;
   std::vector<coot::node_info> b;
   std::vector<coot::node_info> c;
   std::vector<coot::node_info> d;
   std::vector<coot::node_info> e;
   std::vector<coot::node_info> f;
   std::vector<coot::node_info> g;
   std::vector<coot::node_info> h;
   std::vector<coot::node_info> i;
   std::vector<coot::node_info> j;
   std::vector<coot::node_info> k;
   std::vector<coot::node_info> l;
   std::vector<coot::node_info> m;
   
   // 0
   a.push_back(coot::node_info(1));
   a.push_back(coot::node_info(2));
   a.push_back(coot::node_info(5));

   // 1
   b.push_back(coot::node_info(0));
   b.push_back(coot::node_info(2));

   // 2
   c.push_back(coot::node_info(0));
   c.push_back(coot::node_info(1));
   c.push_back(coot::node_info(3));

   // 3
   d.push_back(coot::node_info(2));
   d.push_back(coot::node_info(4));

   // 4
   e.push_back(coot::node_info(3));

   // 5
   f.push_back(coot::node_info(0));
   f.push_back(coot::node_info(6));

   // 6
   g.push_back(coot::node_info(8));
   g.push_back(coot::node_info(7));
   g.push_back(coot::node_info(5));

   // 7
   h.push_back(coot::node_info(6));
   h.push_back(coot::node_info(9));

   // 8
   i.push_back(coot::node_info(6));

   // 9
   j.push_back(coot::node_info(7));

   // 10 
   k.push_back(coot::node_info(11));

   // 11
   l.push_back(coot::node_info(10));
   l.push_back(coot::node_info(12));

   // 12
   m.push_back(coot::node_info(11));

   connection_indices[0] = a;
   connection_indices[1] = b;
   connection_indices[2] = c;
   connection_indices[3] = d;
   connection_indices[4] = e;
   connection_indices[5] = f;
   connection_indices[6] = g;
   connection_indices[7] = h;
   connection_indices[8] = i;
   connection_indices[9] = j;
   connection_indices[10] = k;
   connection_indices[11] = l;
   connection_indices[12] = m;

   for(int i=0; i<n_nodes; i++) 
      coords[i] = clipper::Coord_orth(double(i), double(i), double(i));

   CMMDBManager *mol = new CMMDBManager;
   mol->SetCell(10.0, 10.0, 10.0, 90, 90, 90);
   mol->SetSpaceGroup("P 1");
   
   coot::atom_graph ag(mol,connection_indices, coords);

   ag.sort();

}


void test_alignment(int argc, char **argv) {

   // So the plan is the read a pdb file that is well-refined,
   // and build node graphs for each of the residues.
   // We then calculate distortion of each residue to each residue type.
   //
   // And we can try to slide the sequence too later.

   if (argc > 1) {
      std::string filename(argv[1]);
      coot::minimol::molecule t;
      t.read_file(filename);

      for (unsigned int ifrag=0; ifrag<t.fragments.size(); ifrag++) {
	 for (unsigned int ires=1; ires<t[ifrag].residues.size(); ires++) {

	    // Build a node graph for this residue
	    for (unsigned int iat=0; iat<t[ifrag][ires].atoms.size(); iat++) {
	    }
	 }
      }
   } else {
      std::cout << "Usage: " << argv[0] << " pdb-file-name" << std::endl;
   }
} 
