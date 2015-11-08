

#include "trace.hh"

void
coot::trace::next_vertex(const std::vector<scored_node_t> &path,
			 unsigned int depth, scored_node_t this_scored_vertex) {

   std::vector<scored_node_t> neighbour_vertices =
      get_neighbours_of_vertex_excluding_path(this_scored_vertex.atom_idx, path);

   if (false) {  // debug;
      std::cout << "next_vertex() depth: "
		<< depth << " this_scored_vertex: " << this_scored_vertex.atom_idx << " path (";
      for (unsigned int ip=0; ip<path.size(); ip++)
	 std::cout << path[ip].atom_idx << ",";
      std::cout << ")" << std::endl;
      std::cout << "neighbour_vertices of " << this_scored_vertex.atom_idx << " are : ";
      for (unsigned int in=0; in<neighbour_vertices.size(); in++) { 
	 std::cout << neighbour_vertices[in].atom_idx << "  ";
      }
      std::cout << std::endl;
   }

   for (unsigned int i=0; i<neighbour_vertices.size(); i++) {
	 std::vector<scored_node_t> new_path = path;
	 new_path.push_back(this_scored_vertex);
	 if (new_path.size() >= 4) { 
	    // print_tree(new_path); // for debugging
	    add_tree_maybe(new_path);
	 }

	 const scored_node_t &candidate_vertex = neighbour_vertices[i];
	 const unsigned int &candidate_vertex_idx = neighbour_vertices[i].atom_idx;

	 // if the path is greater than 1, then we can check the angle
	 // of addition of candidate_vertex. It can't be less than 90
	 // degrees (let's say 81 for wiggle room).
	 // 
	 bool angle_ok = true;

	 if (new_path.size() > 1) {
	    double angle = path_candidate_angle(new_path, candidate_vertex_idx);
	    if (angle < 0.9 * M_PI * 0.5)
	       angle_ok = false;
	 }

	 if (false) { // debug
	    if (angle_ok == false) {
	       std::cout << "vertex rejected by angle " << candidate_vertex_idx
			 << " given new_path ";
	       for (unsigned int ii=0; ii<new_path.size(); ii++)
		  std::cout << "  " << new_path[ii].atom_idx;
	       std::cout << std::endl;
	    }
	 }
	 
	 if (angle_ok) {
	    if (false) { // debug
	       std::cout << "calling next_vertex() for candidate-idx "
			 << candidate_vertex_idx << " given path ";
	       for (unsigned int ii=0; ii<new_path.size(); ii++)
		  std::cout << "  " << new_path[ii].atom_idx;
	       std::cout << std::endl;
	    }
	    next_vertex(new_path, depth+1, candidate_vertex);
	 }
   }
}


// not this_vertex or any of the vertices in path.
std::vector<coot::scored_node_t>
coot::trace::get_neighbours_of_vertex_excluding_path(unsigned int this_vertex,
						     const std::vector<scored_node_t> &path) {

   std::vector<scored_node_t> v;
   std::vector<scored_node_t> all_neighbs = tr[this_vertex];

   for (unsigned int ii=0; ii<all_neighbs.size(); ii++) {
      bool add = true;
      for (unsigned int jj=0; jj<path.size(); jj++) { 
	 if (path[jj] == all_neighbs[ii]) {
	    add = false;
	    break;
	 }
      }
      if (all_neighbs[ii].atom_idx == this_vertex)
	 add = false;
      if (add)
	 v.push_back(all_neighbs[ii]);
   }
   return v;
}

