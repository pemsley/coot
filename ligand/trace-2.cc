

#include "trace.hh"

std::vector<unsigned int>
coot::trace::next_vertex(unsigned int start_vertex, const std::vector<unsigned int> &path,
			 unsigned int depth, unsigned int this_vertex) {

   std::vector<unsigned int> v;
   std::vector<unsigned int> neighbour_vertices = get_neighbours_of_vertex_excluding_path(this_vertex, path);

   if (false) {  // debug;
      std::cout << "next_vertex() start_vertex: " << start_vertex << " depth: "
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
	 if (this_vertex != start_vertex)
	    new_path.push_back(this_vertex);
	 if (new_path.size() > 5)
	    print_tree(new_path);
	 v = next_vertex(start_vertex, new_path, depth+1, neighbour_vertices[i]);
   }
   return v;
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
      std::cout << "  " << path[i];
   }
   std::cout << std::endl;

   for (unsigned int i=0; i<path.size(); i++) {
      const clipper::Coord_orth &pt = sas[path[i]]->pos;
      std::cout << "long path " << i
		<< " " << pt.x()
		<< " " << pt.y()
		<< " " << pt.z()
		<< std::endl;
   }
}

