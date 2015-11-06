

#include "trace.hh"

std::vector<unsigned int>
coot::trace::next_vertex(const std::vector<unsigned int> &path,
			 unsigned int depth, unsigned int this_vertex) {

   std::vector<unsigned int> v;
   std::vector<unsigned int> neighbour_vertices = get_neighbours_of_vertex_excluding_path(this_vertex, path);

   if (false) {  // debug;
      std::cout << "next_vertex() depth: "
		<< depth << " this_vertex: " << this_vertex << " path (";
      for (unsigned int ip=0; ip<path.size(); ip++)
	 std::cout << path[ip] << ",";
      std::cout << ")" << std::endl;
      std::cout << "neighbour_vertices of " << this_vertex << " are : ";
      for (unsigned int in=0; in<neighbour_vertices.size(); in++) { 
	 std::cout << neighbour_vertices[in] << "  ";
      }
      std::cout << std::endl;
   }

   for (unsigned int i=0; i<neighbour_vertices.size(); i++) {
	 std::vector<unsigned int> new_path = path;
	 new_path.push_back(this_vertex);
	 if (new_path.size() > 4) { 
	    print_tree(new_path); // for debugging
	    add_tree_maybe(new_path);
	 }

	 const unsigned int &candidate_vertex = neighbour_vertices[i];

	 // if the path is greater than 1, then we can check the angle
	 // of addition of candidate_vertex. It can't be less than 90
	 // degrees (let's say 81 for wiggle room).
	 // 
	 bool angle_ok = true;

	 if (path.size() > 1) {
	    double angle = path_candidate_angle(path, candidate_vertex);
	    if (angle < 0.9 * M_PI *0.5)
	       angle_ok = false;
	 } 
	 
	 if (angle_ok)
	    v = next_vertex(new_path, depth+1, candidate_vertex);
   }
   return v;
}

double
coot::trace::path_candidate_angle(const std::vector<unsigned int> &path,
				  unsigned int candidate_vertex) const {

   unsigned int l = path.size();
   const clipper::Coord_orth &pt_1 = sas[candidate_vertex]->pos;
   const clipper::Coord_orth &pt_2 = sas[path[l-1]]->pos;
   const clipper::Coord_orth &pt_3 = sas[path[l-2]]->pos;
   double angle = clipper::Coord_orth::angle(pt_1, pt_2, pt_3);

   return angle;
} 


// not this_vertex or any of the vertices in path.
   std::vector<unsigned int>
coot::trace::get_neighbours_of_vertex_excluding_path(unsigned int this_vertex,
						     const std::vector<unsigned int> &path) {

   std::vector<unsigned int> v;
   std::vector<unsigned int> all_neighbs = tr[this_vertex];

   for (unsigned int ii=0; ii<all_neighbs.size(); ii++) {
      bool add = true;
      for (unsigned int jj=0; jj<path.size(); jj++) { 
	 if (path[jj] == all_neighbs[ii]) {
	    add = false;
	    break;
	 }
      }
      if (all_neighbs[ii] == this_vertex)
	 add = false;
      if (add)
	 v.push_back(all_neighbs[ii]);
   }
   return v;
}


void
coot::trace::print_tree(const std::vector<unsigned int> &path) const {

   std::cout << "path: ";
   for (unsigned int i=0; i<path.size(); i++) {
      int res_no = nearbyint(sas[path[i]]->occupancy);
      std::cout << "  " << path[i] << " (" << sas[path[i]]->name << " " << res_no << ")";
   }
   std::cout << std::endl;

   if (false) 
      for (unsigned int i=0; i<path.size(); i++) {
	 const clipper::Coord_orth &pt = sas[path[i]]->pos;
	 std::cout << "long path " << i
		   << " " << pt.x()
		   << " " << pt.y()
		   << " " << pt.z()
		   << std::endl;
   }
}



void
coot::trace::add_tree_maybe(const std::vector<unsigned int> &path) {


   // add this path if there is not already a path that is longer that
   // contains at least 4 of the same path points.

   bool add_this = true;
   unsigned int n_match_crit = 4;
   unsigned int n_match_for_replace = 6; // at least this number of matches to replace an existing tree

   for (unsigned int itree=0; itree<interesting_trees.size(); itree++) {
      unsigned int n_match = 0;

      // do I already have something that's longer than path already
      // in interesting_trees? (if so, set add_this to false).
      
      if (path.size() < interesting_trees[itree].size()) {
	 for (unsigned int i=0; i<interesting_trees[itree].size(); i++) {
	    for (unsigned int j=0; j<path.size(); j++) {
	       if (path[j] == interesting_trees[itree][i]) {
		  n_match ++;
	       }

	       if (n_match > n_match_crit) {
		  add_this = false;
		  break;
	       }
	    }
	    if (add_this == false)
	       break;
	 }
      }
   }

   if (0) { 
      std::cout << "add status " << add_this << " for tree of length " << path.size()
		<< " because " << interesting_trees.size() << " trees of lengths ";
      for (unsigned int ii=0; ii<interesting_trees.size(); ii++) { 
	 std::cout << "  " << interesting_trees[ii].size();
      }
      std::cout << std::endl;
   }


   if (add_this) {

      // do we append or replace an existing one?
      int idx_replace = -1;
      unsigned int n_match_best = 0;

      // try to set replace index if we can
      // 
      for (unsigned int itree=0; itree<interesting_trees.size(); itree++) {
	 if (path.size() >= n_match_for_replace) { 
	    if (path.size() > interesting_trees[itree].size()) { 
	       unsigned int n_match = 0;
	       for (unsigned int i=0; i<interesting_trees[itree].size(); i++) {
		  for (unsigned int j=0; j<path.size(); j++) {
		     if (path[j] == interesting_trees[itree][i]) {
			n_match++;
		     }
		  }
	       }

	       std::cout << "n_match " << n_match << " n_match_best " << n_match_best
			 << " for interesting_trees[" << itree << "] of size "
			 << interesting_trees[itree].size() << " vs path size " << path.size() << std::endl;
	       if (n_match >= n_match_for_replace) { 
		  if (n_match > n_match_best) {
		     n_match_best = n_match;
		     idx_replace = itree;
		  }
	       }
	    }
	 }
      }
      

      if (idx_replace >= 0) {
	 std::cout << "replacing tree at idx " << idx_replace << " with tree of length " << path.size() << std::endl;
	 interesting_trees[idx_replace] = path;
	 // delete trees in interesting_trees that are "duplicates" of this one
      } else {
	 std::cout << "adding tree of length " << path.size() << std::endl;
	 interesting_trees.push_back(path);
      }
   } 
}

