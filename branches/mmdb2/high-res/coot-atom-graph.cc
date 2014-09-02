/* high-res/coot-atom-graph.cc
 * 
 * Copyright 2004  The University of York
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


#include <iostream>
#include <string>
#include <algorithm> 

#include "coot-atom-graph.hh"

coot::atom_graph::atom_graph(CMMDBManager *mol,
			     const std::vector<std::vector<coot::node_info> > &connection_indices, 
			     const std::vector<clipper::Coord_orth> &coords) { 

   if (connection_indices.size() != coords.size()) { 
      std::cout << "ERROR:: atom_graph initialization error\n";
   } else { 

      // for every node there is an associated (set of atom_infos)
      int ic = connection_indices.size();
      atom_info.resize(ic);
      nodes = connection_indices;
      
      std::string dum_atom_name(" DUM");
      std::string dum_atom_ele (" C");
      std::string altloc("");
      coot::minimol::atom dummy_atom(dum_atom_name, dum_atom_ele, 
				     clipper::Coord_orth(0.0, 0.0, 0.0), altloc, 30.0);
      atoms.resize(ic, dummy_atom);
      for (int iat=0; iat<coords.size(); iat++)
      	 atoms[iat] = coot::minimol::atom(dum_atom_name, dum_atom_ele, coords[iat], altloc, 30.0);
      mol_internal_ptr_copy = mol;

      // assign cell and space group spgr
      cell.resize(6);
      realtype vol;
      int orthcode;
      mol->GetCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], vol, orthcode);
      spgr = mol->GetSpaceGroup();
   }
}


void
coot::atom_graph::sort() { 

   int ic = nodes.size();
   connectedness.resize(ic);

   for (int i=0; i<ic; i++)
      connectedness[i] = 0;
   std::cout << "There are " << nodes.size() << " nodes in atom_graph::sort\n";

//    std::cout << "##### DEBUG:: " << std::endl;
//    int n_connections;
//    for (int ii=0; ii<nodes.size(); ii++) {
//       // std::cout << " " << ii << "  ";
//       for (int jj=0; jj<nodes[ii].size(); jj++) {
// 	 // std::cout << nodes[ii][jj] << " ";
// 	 n_connections++;
//       }
//       // std::cout << "\n";
//    }
//    std::cout << "##### INFO:: there were " << n_connections 
//              << " connections for those nodes\n";

   int max_con;
   int i_max_con;
   std::vector <std::pair<int, int> > tips; // node number, connectedness

   for (int i=0; i<ic; i++) {
      if (nodes[i].size() > 0) {
	 if (connectedness[i] == 0) {
	    // std::cout << "starting from wrapper " << i << std::endl;
	    trace_along(i, 1);

	    max_con = 0;
	    i_max_con = -1;
	    for (int inode=0; inode<ic; inode++) { 
	       if (connectedness[inode] > max_con) { 
		  max_con = connectedness[inode];
		  i_max_con = inode;
	       }
	       if (connectedness[inode] > 0) { 
		  connectedness[inode] = -connectedness[inode];
	       }
	    }
	    // std::cout << "end of connection at node " << i_max_con << std::endl;
	    // connectedness set later.
	    tips.push_back(std::pair<int, int>(i_max_con, 0));
	 }
      }
   }

   // amusing debugging:
//    std::cout << "Node Connectedness\n";
//    for (int inode=0; inode<ic; inode++)
//       std::cout << "  " << inode << "       " << connectedness[inode] << "\n";

   for (int itip=0; itip<tips.size(); itip++) {
      for (int i=0; i<ic; i++)
	 connectedness[i] = 0;
      trace_along(tips[itip].first, 1);
      // now find the maximum connectedness
      max_con = 0;
      i_max_con = -1;
      for (int i=0; i<ic; i++) { 
	 if (connectedness[i] > max_con) {
	    max_con = connectedness[i];
	    i_max_con = i;
	 }
      }
      if (i_max_con >= 0) { 
	 std::cout << " Tip " << itip << " has max connectedness " 
		   <<  max_con << std::endl;
	 tips[itip].second = max_con;

	 // amusing debugging:
// 	 std::cout << "Node Connectedness\n";
// 	 for (int inode=0; inode<ic; inode++)
// 	    std::cout << "  " << inode << "       " << connectedness[inode] << "\n";

	 
	 std::vector<std::vector<coot::node_info> > t = get_trace(i_max_con);
	 std::cout << "Node Connections\n";
// 	 for (int it=0; it<t.size(); it++) { 
// 	    std::cout << "  " << it << "    ";
// 	    for (int in=0; in<t[it].size(); in++) { 
// 	       std::cout << " " << t[it][in];
// 	    }
// 	    std::cout << "\n";
// 	 }
// 	 std::cout << "-------\n";
	 coot::chain_helper_info chi;
	 chi = peptide_search(t, i_max_con, chi, 0);
	 chi.chain_number = itip;
	 chi = peptide_search(t, i_max_con, chi, 1);
      }
   }

   std::cout << " Atom info in sort:\n";
   for(int iat=0; iat<atom_info.size(); iat++) {
      for (int j=0; j<atom_info[iat].size(); j++) { 
	 std::cout << " atom: " << iat << " info node: " << j << " chain_number:"
		   << atom_info[iat][j].chain_number << " residue_number:" 
		   << atom_info[iat][j].residue_number << " " 
		   << atom_info[iat][j].atom.name << " " 
		   << "\n";
      }
   } 

   assign_c_betas();
   assign_waters();

   sort_tips(&tips);
   std::cout << "---- tips:----" << std::endl;
   for (int itip=0; itip<tips.size(); itip++)
      std::cout << itip << "   " <<  tips[itip].second << "\n";

   // amusing debugging:
//    std::cout << "Node Connectedness\n";
//    for (int inode=0; inode<ic; inode++)
//       std::cout << "  " << inode << "       " << connectedness[inode] << "\n";

   std::string filename("traced.pdb");
   write_molecule_from_atom_info(filename);
}


void
coot::atom_graph::trace_along(int ic, int connection_number) { 

   // std::cout << "trace along " << ic << " " << connection_number << std::endl;

   connectedness[ic] = connection_number;

   for (int i=0; i<nodes[ic].size(); i++) {
//       std::cout << "node " << nodes[ic][i] << " of parent node "
// 			<< ic << " has connectedness " 
// 			<< connectedness[nodes[ic][i]] << std::endl;
      if (connectedness[nodes[ic][i].index] == 0) {
// 	 if (nodes[ic][i].symm_trans_needed_flag)
// 	    std::cout << "Neighbour traced cross boundary at : " 
// 		      << atoms[ic].pos.format() << " with neighbour " 
// 		      << atoms[nodes[ic][i].index].pos.format() << "\n";
	 trace_along(nodes[ic][i].index, (connection_number+1));
      }
   }
}


void
coot::atom_graph::sort_tips(std::vector <std::pair<int, int> > *tips) const { 

   std::sort(tips->begin(), tips->end(), coot::atom_graph::tip_compare);
}

// static 
short int 
coot::atom_graph::tip_compare(const std::pair<int, int> &a, 
			      const std::pair<int, int> &b) { 
   return a.second > b.second ? 1 : 0; 
} 


// i_start_node is the starting node, it is the node with the maximum
// connectedness.  We want to find the connection to it's node that has
// connectedness 1.
// 
// This should be a directed graph.  Perhaps we should pass the
// direction too.  Another day.
//
std::vector<std::vector<coot::node_info> >
coot::atom_graph::get_trace(int i_start_node) const {

   std::vector<std::vector<coot::node_info> > trace(nodes.size());
   std::vector<int> con_local(nodes.size());

   // what's con_local?  Is it useful in digraph_trace_along?
   // 
   for (int i=0; i<con_local.size(); i++) 
      con_local[i] = 0;

//    std::cout << "INFO:: getting digraph_trace from node " << i_start_node 
// 	     << " which has connectedness " << connectedness[i_start_node] 
// 	     << std::endl;

   digraph_trace_along(i_start_node, 
		       connectedness[i_start_node],
		       &trace, &con_local);

   return trace;
}

void
coot::atom_graph::digraph_trace_along(int inode, 
				      int connection_number,
				      std::vector<std::vector<coot::node_info> > *trace,
				      std::vector<int> *con_local) const {

//    std::cout << "Node " << inode << " has connectedness " << connectedness[inode]
// 	     << " and connection_number is " << connection_number << "\n";

   for (int i=0; i<nodes[inode].size(); i++) {
      if (connectedness[nodes[inode][i].index] == (connection_number - 1)) {
	 (*trace)[inode].push_back(nodes[inode][i]);
	 digraph_trace_along(nodes[inode][i].index, 
			     (connection_number - 1), trace, con_local);
      }
   }
}


void
coot::atom_graph::assign_c_betas() {

   for (int index=0; index<atoms.size(); index++) {
      // first find a c-alpha:
      if (is_c_alpha_p(index)) {
	 if (nodes[index].size() == 3) { 
	    for (int j=0; j<nodes[index].size(); j++) {
	       if (atom_info[nodes[index][j].index].size() == 0) {
		  if (atom_info[index].size() > 0) { 
		     int resno =        atom_info[index][0].residue_number;
		     int chain_number = atom_info[index][0].chain_number;
		     std::cout << "DEBUG:: C beta assignment to residue: " << resno << std::endl;
		     // initially the cb tranformations are the ca transformations.
		     std::vector<clipper::RTop_orth> cb_tranformations =
			atom_info[index][0].rtops;
		     // if the link between the ca and the cb needed
		     // symmetry, add it to cb transformations (FYI:
		     // this extra transformation captures 2 more
		     // Cbetas in pad).
		     if (nodes[index][j].symm_trans_needed_flag)
			cb_tranformations.push_back(nodes[index][j].rtop);
		     clipper::Coord_orth tpos = get_transformed_atom(atoms[nodes[index][j].index].pos, cb_tranformations);
		     coot::minimol::atom at(" CB ", " C", tpos, "", 30.0);
		     atom_info[nodes[index][j].index].push_back(coot::graph_atom_info(resno, chain_number, at,
										      cb_tranformations, 0, 1));
		  } else {
		     std::cout << "DEBUG:: Oops! failure to find connections to C-alpha!? " << std::endl;
		  }
	       }
	    }
	 }
      }
   }
}

short int
coot::atom_graph::is_c_alpha_p(int index) const {

   short int istat = 0;
   if (index < atoms.size() && index >= 0) { 
      if (atom_info[index].size() > 0) {
	 if (atom_info[index][0].atom.name == " CA ")
	    istat = 1;
      }
   }
   return istat;
}



void
coot::atom_graph::assign_waters() {

   int chain_number = 22; // W
   int resno = 1;
   for (int i=0; i<nodes.size(); i++) {
      if (nodes[i].size() == 0) { 
	 // it has no neigbours.  A water then?
	 //
	 coot::minimol::atom at(" O  ", " O", atoms[i].pos, "", 30.0);
	 short int is_water_flag = 1;
	 std::vector<clipper::RTop_orth> water_transformations;
	 atom_info[i].push_back(coot::graph_atom_info(resno, chain_number,
						      at, water_transformations,
						      is_water_flag, 1));
	 resno++; // ready for next water.
      }
   }
} 

// Return the chain score.
//
// we call this function twice.  Once to get the direction, the other
// to make the assignments given that we know the direction and the
// length.
// 
coot::chain_helper_info
coot::atom_graph::peptide_search(const std::vector<std::vector<coot::node_info> > &t,
				 int i_node_start,
				 const coot::chain_helper_info &hi_in,
				 short int make_assignments_flag) {

   coot::chain_helper_info hi_out;
   int confirmed_direction = 0;
   int dir_offset = 0;
   int res_no = -1; // unset
   int chain_number = hi_in.chain_number;
   int made_an_assignment = 0;
   int skip = 0;
   

   if (make_assignments_flag) {
      confirmed_direction = hi_in.get_direction();
      if (confirmed_direction != 0) {
	 dir_offset = confirmed_direction;
	 if (confirmed_direction == 1){
	    res_no = 1;
	 } else {
	    res_no = hi_in.n_peptides + 2;
	 }
      }

      std::cout << "INFO:: peptide_search: chain " << (hi_in.chain_number)
		<< " start resno: " << res_no
		<< " had " << hi_in.n_peptides
		<< " residues " // starting at " << atoms[i_node_start].pos.format()
 		<< " has direction " << confirmed_direction << std::endl;
   }

   int i_node_peptide_start = i_node_start;
   // int t_index_1;
   int t_index_2, t_index_3, t_index_4;
   int n_tripep = 0;
   int direction_score = 0;
   int forward_score = 0;
   int backward_score = 0;
   int this_peptide_connectivity_score = 0;
   int peptide_atom_counter = 0;
   // dummy argments needed for construction... they are assigned
   // before they are used in anger.
   // coot::node_info t_node_1(-),
   coot::node_info t_node_2(-1), t_node_3(-1), t_node_4(-1);
   coot::node_info t_node_5(-1), t_node_6(-1), t_node_7(-1);
   std::vector<clipper::RTop_orth> running_ops;

   for (int i_node_peptide_start=i_node_start; 
	t[i_node_peptide_start].size()>0; 
	i_node_peptide_start = t[i_node_peptide_start][0].index) { // CHECKME

      peptide_atom_counter = 0;

      for (int ip1 = 0; ip1<t[i_node_peptide_start].size(); ip1++) {
	 t_index_2 = t[i_node_peptide_start][ip1].index;
	 t_node_2  = t[i_node_peptide_start][ip1];

	 for (int ip2 = 0; ip2<t[t_index_2].size(); ip2++) {
	    t_index_3 = t[t_index_2][ip2].index;
	    t_node_3  = t[t_index_2][ip2];

	    for (int ip3 = 0; ip3<t[t_index_3].size(); ip3++) {
	       t_index_4 = t[t_index_3][ip3].index;
	       t_node_4  = t[t_index_3][ip3];
			   
	       // 			      std::cout << "peptide indices: " 
	       // 					<< t_index_1 << " "
	       // 					<< t_index_2 << " "
	       // 					<< t_index_3 << " "
	       // 					<< t_index_4 << "\n";

	       if (is_possible_ca_c_n_c(i_node_peptide_start, 
					t_node_2, t_node_3, t_node_4)) {
		  
		  int t_index_5, t_index_6, t_index_7;

		  for (int ip4=0; ip4<t[t_index_4].size(); ip4++) { 
		     t_index_5 = t[t_index_4][ip4].index;
		     t_node_5 = t[t_index_4][ip4];
		     for (int ip5=0; ip5<t[t_index_5].size(); ip5++) { 
			t_index_6 = t[t_index_5][ip5].index;
			t_node_6 = t[t_index_5][ip5];
			for (int ip6=0; ip6<t[t_index_6].size(); ip6++) { 
			   t_index_7 = t[t_index_6][ip6].index;
			   t_node_7 = t[t_index_6][ip6];

			   short int direction = 1;

			   double d = peptide_distortion_score(direction,  
							       i_node_peptide_start,
							       t_node_2, 
							       t_node_3, 
							       t_node_4,
							       t_node_5, 
							       t_node_6, 
							       t_node_7);

			   std::cout << "peptide score: " << d << "\n";

			   n_tripep++;
			   this_peptide_connectivity_score = 0;


			   // --------------------------------------
			   //           forwards
			   // -------------------------------------- 
			   if (nodes[t_index_2].size() == 3) { 
			      // std::cout << "carbonyl backward\n";
			      direction_score++;
			      forward_score++;
			      this_peptide_connectivity_score++;
			      if (nodes[t_index_5].size() == 3) {
				 direction_score++;
				 forward_score++;
				 this_peptide_connectivity_score++;
				 // std::cout << "carbonyl backward backward\n";
			      }
			   }
			   if (nodes[t_index_3].size() == 2) {
			      // std::cout << "backward\n";
			      direction_score++;
			      forward_score++;
			      this_peptide_connectivity_score++;
			      if (nodes[t_index_6].size() == 2) {
				 direction_score++;
				 forward_score++;
				 // std::cout << "backward backward\n";
				 this_peptide_connectivity_score++;
			      }
			   }

			   // --------------------------------------
			   //           backwards
			   // -------------------------------------- 
			   if (nodes[t_index_2].size() == 2) { 
			      // std::cout << "forward\n";
			      direction_score--;
			      backward_score--;
			      this_peptide_connectivity_score--;
			      if (nodes[t_index_5].size() == 2) { 
				 // std::cout << "forward forward\n";
				 backward_score--;
				 direction_score--;
				 this_peptide_connectivity_score--;
			      }
			   }

			   if (nodes[t_index_3].size() == 3) { 
			      direction_score--;
			      backward_score--;
			      this_peptide_connectivity_score--;
			      if (nodes[t_index_6].size() == 3) { 
				 backward_score--;
				 direction_score--;
				 this_peptide_connectivity_score--;
			      }
			   }
			   
			   if (d < 0.1) {
			      if (this_peptide_connectivity_score >=  1 ||
				  this_peptide_connectivity_score <= -1) {

				 if (confirmed_direction != 0) { 
				    if (make_assignments_flag) {
				       running_ops = 
					  make_assignments(confirmed_direction,
							   d,
							   res_no, 
							   chain_number,
							   i_node_peptide_start,
							   running_ops,
							   t_node_2,
							   t_node_3,
							   t_node_4);
				       res_no += confirmed_direction;
				       if (peptide_atom_counter > 3)
					  res_no++;
				       peptide_atom_counter = 0;
				    }
				 }
				 made_an_assignment++; // need to count peptide the
				                       // first time round.
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   std::cout << "INFO:: storing chain with "
	     << made_an_assignment << " assignments "
	     << " make_assignments_flag " << make_assignments_flag << "\n";
   hi_out = coot::chain_helper_info(direction_score,
				    forward_score, backward_score,
				    made_an_assignment); // make_assignments is n_peptides

   if (n_tripep > 0) 
      std::cout << "chain of length " << n_tripep << " peptides: " << direction_score
		<< " forward score: " << forward_score << " backward_score: "
		<< backward_score << std::endl;
   return hi_out;
}

// Return +1 for positive, -1 for negative and 0 for undetermined.
// 
int
coot::chain_helper_info::get_direction() const {

   int direction = 0;
   if (direction_score > 0)
      if (direction_score > n_peptides)
	 direction = 1;
   if (direction_score < 0)
      if (direction_score < -n_peptides)
	 direction = -1;

   return direction;
}


// Fill the atom_info vector.
// 
// Get your thinking cap on, because this is complicated.
//
// We need to transform this peptide by the result of the
// transformations of all the other peptides in this chain - hence we
// pass running_ops (and the returned thing becomes the new running_ops).
//
// We need to keep a track of what transformations have been applied
// here, so we return a new running_ops which is transformations_4.
//
// We need to assign the carbonyl oxygens too.  Slightly tricky, in
// the forward direction the oxygens are connected to t_node_2 and in
// the backwards direction to t_node_3.  So the transformations to go
// from the t_node of the C to the O depends on the direction.
// 
std::vector<clipper::RTop_orth>
coot::atom_graph::make_assignments(int direction,
				   double weight,
				   int res_no,
				   int chain_number,
				   int i_node_peptide_start,
				   const std::vector<clipper::RTop_orth> &running_ops,
				   const coot::node_info &tnode_2,
				   const coot::node_info &tnode_3,
				   const coot::node_info &tnode_4) {

   std::cout << "make_assignments: with res_no: " << res_no << std::endl;

   int index_2 = tnode_2.index;
   int index_3 = tnode_3.index;
   int index_4 = tnode_4.index;

   std::vector<clipper::RTop_orth> transformations_1 = running_ops;
   std::vector<clipper::RTop_orth> transformations_2 = running_ops;
   std::vector<clipper::RTop_orth> transformations_3 = running_ops;
   std::vector<clipper::RTop_orth> transformations_4 = running_ops;
   std::vector<clipper::RTop_orth> transformations_o = running_ops;

   if (tnode_2.symm_trans_needed_flag) { 
      transformations_2.push_back(tnode_2.rtop);
      transformations_3.push_back(tnode_2.rtop);
      transformations_4.push_back(tnode_2.rtop);
   }

   if (tnode_3.symm_trans_needed_flag) { 
      transformations_3.push_back(tnode_3.rtop);
      transformations_4.push_back(tnode_3.rtop);
   }

   if (tnode_4.symm_trans_needed_flag) { 
      transformations_4.push_back(tnode_4.rtop);
   }


   int resno_ca_index_4;
   int resno_index_2;
   int resno_index_3;
   std::string index_2_at_name;
   std::string index_2_ele_name;
   std::string index_3_at_name;
   std::string index_3_ele_name;
   short int do_o_assign = 0;
   int oxygen_res_no = -1; // unset
   int oxygen_index = 0;   // unset
   short int is_water_flag = 0;
   clipper::Coord_orth tpos;
   
   if (direction > 0) {
      resno_index_2 = res_no;
      resno_index_3 = res_no + 1;
      resno_ca_index_4 = res_no + 1;
      index_2_at_name  = " C  ";
      index_2_ele_name = " C";
      index_3_at_name  = " N  ";
      index_3_ele_name = " N";

      // assign carboxyl oxygen if the C has 3 edges:

      if (nodes[index_2].size() == 3) {
	 for (int i=0; i<3; i++) {
	    if (nodes[index_2][i].index != i_node_peptide_start &&
		nodes[index_2][i].index != index_3) {
	       do_o_assign = 1;
	       oxygen_res_no = res_no;
	       oxygen_index = nodes[index_2][i].index;

	       transformations_o = running_ops;
	       if (tnode_2.symm_trans_needed_flag)
		  transformations_o.push_back(tnode_2.rtop);
	       
	       if (nodes[index_2][i].symm_trans_needed_flag) {
		  transformations_o.push_back(nodes[index_2][i].rtop);
	       }
	       
	       break;
	    }
	 }
      }
   } else {
      resno_index_2 = res_no;
      resno_index_3 = res_no - 1;
      resno_ca_index_4 = res_no - 1;
      index_2_at_name  = " N  ";
      index_2_ele_name = " N";
      index_3_at_name  = " C  ";
      index_3_ele_name = " C";
      // assign carboxyl oxygen if the C has 3 edges:

      if (nodes[index_3].size() == 3) {
	 for (int i=0; i<3; i++) {
	    if (nodes[index_3][i].index != index_4 &&
		nodes[index_3][i].index != index_2) {
	       do_o_assign = 1;
	       oxygen_res_no = res_no - 1;
	       oxygen_index = nodes[index_3][i].index;

	       transformations_o = running_ops;
	       if (tnode_2.symm_trans_needed_flag)
		  transformations_o.push_back(tnode_2.rtop);
	       if (tnode_3.symm_trans_needed_flag)
		  transformations_o.push_back(tnode_3.rtop);
	       
	       if (nodes[index_3][i].symm_trans_needed_flag) {
		  transformations_o.push_back(nodes[index_3][i].rtop);
	       }
	       
	       break;

	    }
	 }
      }
   }

   tpos = get_transformed_atom(atoms[i_node_peptide_start].pos, transformations_1);
   coot::minimol::atom at1(" CA ", " C", tpos, "", 30.0);
   atom_info[i_node_peptide_start].push_back(coot::graph_atom_info(res_no, chain_number, at1, transformations_1, is_water_flag, weight));

   tpos = get_transformed_atom(atoms[index_2].pos, transformations_2);
   coot::minimol::atom at2(index_2_at_name, index_2_ele_name, tpos, "", 30.0);
   atom_info[index_2].push_back(coot::graph_atom_info(resno_index_2, chain_number, at2, transformations_2, is_water_flag, weight));

   tpos = get_transformed_atom(atoms[index_3].pos, transformations_3);
   coot::minimol::atom at3(index_2_at_name, index_3_ele_name, tpos, "", 30.0);
   atom_info[index_3].push_back(coot::graph_atom_info(resno_index_3, chain_number, at3, transformations_3, is_water_flag, weight));
   
   tpos = get_transformed_atom(atoms[index_4].pos, transformations_4);
   coot::minimol::atom at4(" CA ", " C", tpos, "", 30.0);
   atom_info[index_4].push_back(coot::graph_atom_info(resno_ca_index_4, chain_number, at4, transformations_4, is_water_flag, weight));

   if (do_o_assign) {
      tpos = get_transformed_atom(atoms[oxygen_index].pos, transformations_o);
      coot::minimol::atom ox_at(" O  ", " O", tpos, "", 30.0);
      atom_info[oxygen_index].push_back(coot::graph_atom_info(oxygen_res_no, chain_number, ox_at, transformations_o, is_water_flag, weight));
   }

   return transformations_4;
}



short int 
coot::atom_graph::is_possible_ca_c_n_c(int index_start_peptide,
				       const coot::node_info &t_node_2, 
				       const coot::node_info &t_node_3, 
				       const coot::node_info &t_node_4) const {
   
   short int state;
   std::vector<clipper::RTop_orth> transformations;

//    std::cout << "indices: " 
// 	     << "                   " << t_node_2.index << " " 
// 	     << t_node_3.index << " " << t_node_4.index << "\n";

   clipper::Coord_orth a1 = atoms[index_start_peptide].pos;
   clipper::Coord_orth a2 = atoms[t_node_2.index].pos;
   if (t_node_2.symm_trans_needed_flag) { 
      transformations.push_back(t_node_2.rtop);
      a2 = get_transformed_atom(a2, transformations);
   }
      
   clipper::Coord_orth a3 = atoms[t_node_3.index].pos;
   if (t_node_3.symm_trans_needed_flag)
      transformations.push_back(t_node_3.rtop);
   if (transformations.size() > 0) 
      a3 = get_transformed_atom(a3, transformations);

   clipper::Coord_orth a4 = atoms[t_node_4.index].pos;
   if (t_node_4.symm_trans_needed_flag)
      transformations.push_back(t_node_4.rtop);
   if (transformations.size() > 0) 
      a4 = get_transformed_atom(a4, transformations);

//    std::cout << "checking omega on " 
// 	     << a1.format() << "\n" << "                  " 
// 	     << a2.format() << "\n" << "                  " 
// 	     << a3.format() << "\n" << "                  " 
// 	     << a4.format() << "\n";
      
// 	       double len1 = clipper::Coord_orth::length(a1, a2);
// 	       double len2 = clipper::Coord_orth::length(a2, a3);
// 	       double len3 = clipper::Coord_orth::length(a3, a4);
// 	       double len4 = clipper::Coord_orth::length(a1, a4);
// 	       double d_bond = 1.45;
// 	       double bit = 0.3;

// 		  if (len1 > (d_bond-bit) && len1 < (d_bond+bit)) { 
// 		     if (len2 > (d_bond-bit) && len2 < (d_bond+bit)) { 
// 			if (len3 > (d_bond-bit) && len3 < (d_bond+bit)) {
// 			   // forces trans peptide:
// 			   if (len4 > 3.5 && len4 < 4.0) {

   double torsion  = clipper::Coord_orth::torsion(a1, a2, a3, a4);
//    std::cout << "Note: used " << transformations.size() << " transformations in single peptide search\n";
   std::cout << "omega (torsion): " << clipper::Util::rad2d(torsion) << "\n";
   if (fabs(torsion) > 165.0*3.14/180.0) { 
      std::cout << "interesting Ca candidate at " << a1.format() 
		<< " " << clipper::Util::rad2d(torsion) << "\n";
      state = 1;
   } else { 
      state = 0;
   }
   return state;
}


clipper::Coord_orth
coot::atom_graph::get_transformed_atom(const clipper::Coord_orth &a,
				       const std::vector<clipper::RTop_orth> &transformations) const { 

   clipper::Coord_orth a_copy = a;

   for (int i=transformations.size()-1; i>=0; i--) 
      a_copy = a_copy.transform(transformations[i]);

//     std::cout << "Untransformed atom:     " << a.format() << "\n" 
// 	      << "  transformationed atom " << a_copy.format() << "\n"
// 	      << "  via " << transformations.size() << " rtops:\n "  << "\n";

//    for (int i=0; i<transformations.size(); i++) 
//       std::cout << transformations[i].format() << "\n";

   return a_copy;
}


double 
coot::atom_graph::peptide_distortion_score(short int direction,
					   int index_start,
					   const coot::node_info &t_node_2, 
					   const coot::node_info &t_node_3, 
					   const coot::node_info &t_node_4, 
					   const coot::node_info &t_node_5, 
					   const coot::node_info &t_node_6, 
					   const coot::node_info &t_node_7) const { 
   double distortion = -1;

   std::vector<clipper::RTop_orth> transformations;

   clipper::Coord_orth a1 = atoms[index_start].pos;
   clipper::Coord_orth a2 = atoms[t_node_2.index].pos;
   if (t_node_2.symm_trans_needed_flag) { 
      transformations.push_back(t_node_2.rtop);
      a2 = get_transformed_atom(a2, transformations);
   }
      
   clipper::Coord_orth a3 = atoms[t_node_3.index].pos;
   if (t_node_3.symm_trans_needed_flag)
      transformations.push_back(t_node_3.rtop);
   if (transformations.size() > 0) 
      a3 = get_transformed_atom(a3, transformations);

   clipper::Coord_orth a4 = atoms[t_node_4.index].pos;
   if (t_node_4.symm_trans_needed_flag)
      transformations.push_back(t_node_4.rtop);
   if (transformations.size() > 0) 
      a4 = get_transformed_atom(a4, transformations);

   clipper::Coord_orth a5 = atoms[t_node_5.index].pos;
   if (t_node_5.symm_trans_needed_flag)
      transformations.push_back(t_node_5.rtop);
   if (transformations.size() > 0) 
      a5 = get_transformed_atom(a5, transformations);

   clipper::Coord_orth a6 = atoms[t_node_6.index].pos;
   if (t_node_6.symm_trans_needed_flag)
      transformations.push_back(t_node_6.rtop);
   if (transformations.size() > 0) 
      a6 = get_transformed_atom(a6, transformations);

   clipper::Coord_orth a7 = atoms[t_node_7.index].pos;
   if (t_node_5.symm_trans_needed_flag)
      transformations.push_back(t_node_7.rtop);
   if (transformations.size() > 0) 
      a7 = get_transformed_atom(a7, transformations);
   
   distortion = peptide_distortion_score(direction, a1, a2, a3, a4, a5, a6, a7);

   std::cout << "Note: used " << transformations.size() << " transformations\n";
   return distortion;
}





// return distortion : 0 -> inf
// -1 means failure.
// 
double
coot::atom_graph::peptide_distortion_score(short int direction,
			      const clipper::Coord_orth &a1,
			      const clipper::Coord_orth &a2,
			      const clipper::Coord_orth &a3,
			      const clipper::Coord_orth &a4,
			      const clipper::Coord_orth &a5,
			      const clipper::Coord_orth &a6,
			      const clipper::Coord_orth &a7) const { 
   double distortion = -1;

//    std::cout << "DEBUG:: peptide atom positions: \n" 
// 	     << "      " << a1.format() << "\n"
// 	     << "      " << a2.format() << "\n"
// 	     << "      " << a3.format() << "\n"
// 	     << "      " << a4.format() << "\n"
// 	     << "      " << a5.format() << "\n"
// 	     << "      " << a6.format() << "\n"
// 	     << "      " << a7.format() << "\n";

// distances:  Ca-C   1-2 1.52
//             Ca-N   1-2 1.46
//             C-N    1-2 1.33
//             C-Ca   1-3 2.43
//             N-C    1-3 2.46
//             N-Ca   1-3 2.43

   if (direction == 1) { 

      // a1  a2 a3 a4 a5 a6 a7 
      // Ca   N  C Ca N  C  Ca
      // 
      double len1 = clipper::Coord_orth::length(a1, a2);
      double len2 = clipper::Coord_orth::length(a2, a3);
      double len3 = clipper::Coord_orth::length(a3, a4);
      double len4 = clipper::Coord_orth::length(a1, a4);

      double len5  = clipper::Coord_orth::length(a1, a3);
      double len6  = clipper::Coord_orth::length(a2, a4);
      double len7  = clipper::Coord_orth::length(a3, a5);
      double len8  = clipper::Coord_orth::length(a4, a6);
      double len9  = clipper::Coord_orth::length(a5, a7);

      distortion  = squared(len1 - 1.46);
      distortion += squared(len2 - 1.33);
      distortion += squared(len3 - 1.52);

      distortion += squared(len5  - 2.43);
      distortion += squared(len6  - 2.43);
      distortion += squared(len7  - 2.46);
      distortion += squared(len8  - 2.43);
      distortion += squared(len9  - 2.43);
   }

   if (direction == -1) {

   }

   return distortion;
}

void
coot::atom_graph::write_molecule_from_atom_info(const std::string &file_name) const {

   coot::minimol::molecule mol;
   std::string chain_id_local;
   int fragment_no;

   int i_ass_atoms = 0;
   for (int i=0; i<atom_info.size(); i++) {
      if (atom_info[i].size() > 0)
	 i_ass_atoms++;
   }
   std::cout << "INFO:: " << i_ass_atoms << " out of " << atom_info.size()
	     << " atoms have type assignments\n";
   
   for (int i=0; i<atom_info.size(); i++) {
      if (atom_info[i].size() > 0) {
	 chain_id_local = chain_id(atom_info[i][0].chain_number);
	 if (atom_info[i][0].residue_number > 0 ) {
	    std::cout << "Adding atom to residue " << atom_info[i][0].residue_number
		      << chain_id_local << " " << atom_info[i][0].atom.name << " "
		      << atom_info[i][0].atom.pos.format() << std::endl;
	    mol.addatom (chain_id_local, atom_info[i][0].residue_number, atom_info[i][0].atom, atom_info[i][0].is_water_flag);
	 } else {
	    std::cout << "Ooops! Can't add residue with residue number "
		      << atom_info[i][0].residue_number << std::endl;
	 } 
      }
      
      if (atom_info[i].size() > 1) {
	 std::cout << "We have " << atom_info[i].size() << " infos for atom at "
		   << atoms[i].pos.format() << std::endl;
	 for (int j=0; j<atom_info[i].size(); j++) {
	    std::cout << j << " " << chain_id(atom_info[i][j].chain_number) << " "
		      << atom_info[i][j].weight << " "
		      << atom_info[i][j].residue_number << " "
		      << atom_info[i][j].atom.name << " "
		      << atom_info[i][j].atom.pos.format()  << "\n";
	 }
      }
   }
   mol.set_cell(cell);
   mol.set_spacegroup(spgr);
   mol.write_file(file_name, 20.0);
}

std::string
coot::atom_graph::chain_id(int chain_number) const {

   std::string s("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

   std::string r;
   if (chain_number > 51) 
      r = "Z";
   else
      r = s[chain_number];

   return r;
}

					     
coot::minimol::molecule
coot::atom_graph::traced() const {

   coot::minimol::molecule mol;

   return mol;
}


// What's the plan?
//
// We need to fit every side chain type to every C-alpha.
// For each residue type, we return a score.
// 
// We store the score for that side chain type in a deeply indexed
// vector: side_chain_score[chain_number][residue_number][coot::side_chain_name_index]
//
// Later we will try to fit the scores to a sequence.  How do we do
// that?  See sequence-assignment.
// 
void
coot::atom_graph::sidechains_search() {

   // What's the plan?
   //
   // We adjust the side_chain_score vector/class here.

   double d;
   int i_node_ca; 
   std::vector<coot::node_info> side_chain_nodes;

   // The class of much consideration:
   coot::sequence_assignment::side_chain_score_t scs;

   // We need to get distortion scores like this:
   // d = distortion_score_side_chain(i_node_ca, std::string("SER"), side_chain_nodes);
   // That d needs to be also associated with a chain_number and a residue_number.
   // That's OK, they're part of the atom_info.  The tricky thing is to get the 
   // sidechain atoms.

   // Run through the graph_atom_info vector (atom_info) looking for C-alphas.
   // When we find one, we need to find the side chain atoms.
   // 
   for (int i=0; i<atom_info.size(); i++) { 
      for (int j=0; j<atom_info[i].size(); j++) { 
      
	 if (atom_info[i][j].atom.name == " CA ") { 
	    // we've got a C-alpha
	    side_chain_nodes = get_side_chain_nodes(i, coot::sequence_assignment::SER);
	    // add stuff to scs
	    score_all_side_chain_types(i, side_chain_nodes, &scs);
	 }
      }
   }

   scs.debug();
   scs.slider();

}

// So this is a C-alpha.
//
// Now we want the side chains ordered just so (like my notes)...
// 
std::vector<coot::node_info>
coot::atom_graph::get_side_chain_nodes(int i_node_ca, 
				       coot::sequence_assignment::side_chain_name_index scn) const { 

   // Actually we want the index into the nodes vector to make side
   // chain node vector...
   //

   // What is the relationship between nodes and atom_info?  atom_info
   // is the set of atom infos for each node (there can be more than
   // one atom assignment for each node).
   
   std::vector<coot::node_info> v;
   if (scn == coot::sequence_assignment::SER) {

      if (nodes[i_node_ca].size() > 2) {

	 // i.e. this Ca seem to have a side chain (not a GLY)
	 
	 for (int j=0; j<nodes[i_node_ca].size(); j++) {
	    
	    // OK, so what are the unassigned atoms?
	    // We start off at the CBeta, which (currently) has been assigned

	    // We want an atom that is not assigned as mainchain (C or N)
	    // 
	    for (int ati=0; ati<nodes[i_node_ca].size(); ati++) {
	       if (atom_info[i_node_ca][ati].atom.name != " C  " &&
		   atom_info[i_node_ca][ati].atom.name != " N  ") { 
		  std::cout << "AAAARRRRGGGHHH too complicated!\n";
	       }
	    }
	 }
      }
   }
   return v;
}

// Side chain independent atom ordering
// 
std::vector<coot::node_info>
coot::atom_graph::get_side_chain_nodes(int i_node_ca) const {
				       

   std::vector<coot::node_info> v;

   return v;
}


double
coot::atom_graph::abs_chiral_vol(const clipper::Coord_orth &a1,
				 const clipper::Coord_orth &a2,
				 const clipper::Coord_orth &a3,
				 const clipper::Coord_orth &a4) const {

   clipper::Coord_orth v1 = a2 - a1;
   clipper::Coord_orth v2 = a3 - a1;
   clipper::Coord_orth v3 = a4 - a1;

   return fabs(clipper::Coord_orth::dot(v1, clipper::Coord_orth::cross(v2, v3)));

}


