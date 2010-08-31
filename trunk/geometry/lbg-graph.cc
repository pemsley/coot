
#include <iostream>

#include "lbg-graph.hh"


coot::aromatic_graph_t::aromatic_graph_t(const std::vector<std::pair<std::string, std::string> > &bonds_in) {

   original_bonds = bonds_in;
   for (unsigned int i=0; i<bonds_in.size(); i++) {
      std::string name_1 = bonds_in[i].first;
      std::string name_2 = bonds_in[i].second;
      int index_1 = -1;
      for (unsigned int iv=0; iv<vertices.size(); iv++) {
	 if (vertices[iv].match_name(name_1)) {
	    index_1 = iv;
	    break;
	 }
      }
      if (index_1 == -1) {
	 coot::lbg_vertex vertex(name_1);
	 vertices.push_back(vertex);
	 index_1 = vertices.size() -1; // the index of vertices.back()
      }

      int index_2 = -1;
      for (unsigned int iv=0; iv<vertices.size(); iv++) {
	 if (vertices[iv].match_name(name_2)) {
	    index_2 = iv;
	    break;
	 }
      }
      if (index_2 == -1) {
	 coot::lbg_vertex vertex(name_2);
	 vertices.push_back(vertex);
	 index_2 = vertices.size() -1; // the index of vertices.back()
      }

      coot::lbg_edge edge(index_1, index_2);
      edges.push_back(edge);
      int edge_index = edges.size() - 1;
      vertices[index_1].add_edge_index(edge_index);
      vertices[index_2].add_edge_index(edge_index);
   }

   // print();
} 
      
std::vector<std::vector<std::string> >
coot::aromatic_graph_t::ring_list() {

   for (unsigned int iv=0; iv<vertices.size(); iv++) {
      std::vector<int> path;
      if (0)
	 std::cout << "============ ring check for vertex " << iv
		   << " ===================" << std::endl;
      next_vertex(iv, path, 0, iv);
   }

   // do any of the rings in the rings vector contain smaller rings?
   // If so, reject the bigger ring.
   // 
   std::vector<std::vector<int> > filtered_rings;

   for (unsigned int i=0; i<rings.size(); i++) {
      bool contains_sub_ring = 0;
      for (unsigned int j=0; j<rings.size(); j++) {
	 if (i!=j) {

	    // does ring i contain ring j?
	    // 
	    if (rings[i].size() > rings[j].size()) {
	       
	       // are 3 or more of the atoms in the smaller ring (j)
	       // contained in the bigger ring (i)?
	       // 
	       int nfound = 0;
	       for (unsigned int ii=0; ii<rings[i].size(); ii++) {
		  for (unsigned int jj=0; jj<rings[j].size(); jj++) {
		     if (rings[i][ii] == rings[j][jj]) {
			nfound++;
			break;
		     }
		  }
	       }
	       if (nfound > 2) {
		  // OK, so ring i contains ring j
		  contains_sub_ring = 1;
		  break;
	       }
	    }
	 }
      }
      if (contains_sub_ring == 0)
	 filtered_rings.push_back(rings[i]);
   }
   
   
   return indexes_to_names(filtered_rings); // converts rings to
					       // name from indices.
}


void
coot::aromatic_graph_t::print() const {

   for (unsigned int i=0; i<edges.size(); i++) { 
      std::cout << "edge " << i << ": "
		<< vertices[edges[i].get_vertex_index_1()].get_name() << " to "
		<< vertices[edges[i].get_vertex_index_2()].get_name()
		<< std::endl;
   }
   for (unsigned int i=0; i<vertices.size(); i++) {
      std::cout << "vertex " << i << ": " << vertices[i].get_name()
		<< " had edges ";
      for (unsigned int j=0; j<vertices[i].get_edges().size(); j++) { 
	 std::cout << vertices[i].get_edges()[j] << " ";
      }
      std::cout << std::endl;
   }

}

std::vector<int>
coot::aromatic_graph_t::next_vertex(int start_vertex, const std::vector<int> &path, int depth, int this_vertex) {

   std::vector<int> v;
   std::vector<int> neighbour_vertices = get_neighbours_of_vertex_excluding_path(this_vertex, path);

   if (0) {  // debug;
      std::cout << "next_vertex() start_vertex: " << start_vertex << " depth: " << depth << " this_vertex: "
		<< this_vertex << " path (";
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
      if (neighbour_vertices[i] == start_vertex) {
	 if (depth > 1) { 

	    // yay, a ring.
	    
	    std::vector<int> circular_path = path;
	    circular_path.push_back(this_vertex);
	    circular_path.push_back(start_vertex);
	    // print_path(circular_path);
	    add_path_maybe(circular_path);
	 }
      } else {

	 if (depth < 9 ) { 
	    std::vector<int> new_path = path;
	    if (this_vertex != start_vertex)
	       new_path.push_back(this_vertex);
	    next_vertex(start_vertex, new_path, depth+1, neighbour_vertices[i]);
	 }
      } 
   }
   return v;
}


// get neighbouring vertices of this_vertex
// 
std::vector<int>
coot::aromatic_graph_t::get_neighbours_of_vertex(int this_vertex) const {

   std::vector<int> v;
   std::vector<int> local_edges = vertices[this_vertex].get_edges();
   for (unsigned int i=0; i<local_edges.size(); i++) {
      int v_1 = edges[local_edges[i]].get_vertex_index_1();
      int v_2 = edges[local_edges[i]].get_vertex_index_2();
      if (v_1 != this_vertex)
	 v.push_back(v_1);
      if (v_2 != this_vertex)
	 v.push_back(v_2);
   }
   return v;
}

// not this_vertex or any of the vertices in path.
std::vector<int>
coot::aromatic_graph_t::get_neighbours_of_vertex_excluding_path(int this_vertex,
								      const std::vector<int> &path) const {

   std::vector<int> v;
   std::vector<int> local_edges = vertices[this_vertex].get_edges();

   if (0) { // debug
      std::cout << "local_edges (" << local_edges.size() << ") of " << this_vertex << " are : ";
      for (unsigned int in=0; in<local_edges.size(); in++) { 
	 std::cout << local_edges[in] << "  ";
      }
      std::cout << std::endl;
   }
   
   for (unsigned int i=0; i<local_edges.size(); i++) { 
      int v_1 = edges[local_edges[i]].get_vertex_index_1();
      int v_2 = edges[local_edges[i]].get_vertex_index_2();
      if (v_1 != this_vertex) {
	 bool ifound = 0;
	 for (unsigned int j=0; j<path.size(); j++) { 
	    if (path[j] == v_1) {
	       ifound = 1;
	       break; // no
	    }
	 }
	 if (ifound == 0) { // not in path (either)
	    v.push_back(v_1);
	 } 
      }
      if (v_2 != this_vertex) {
	 bool ifound = 0;
	 for (unsigned int j=0; j<path.size(); j++) { 
	    if (path[j] == v_2) {
	       ifound = 1;
	       break; // no
	    }
	 }
	 if (ifound == 0) { // not in path (either)
	    v.push_back(v_2);
	 }
      }
   }
   return v;
}

void
coot::aromatic_graph_t::print_path(std::vector<int> &path) const {

   std::cout << "================ path: =========== "; 
   for (unsigned int i=0; i<path.size(); i++) { 
      std::cout << vertices[path[i]].get_name() << " ";
   }
   std::cout << std::endl;

}

std::ostream &
coot::operator<<(std::ostream &s, lbg_edge e) {

   s << "edge{" << e.vertex_1 << "," << e.vertex_2 << "}";
   return s;
} 

// add circular_path only if it is not in rings already
void
coot::aromatic_graph_t::add_path_maybe(std::vector<int> circular_path) {

   std::sort(circular_path.begin(), circular_path.end());

   bool ifound = 0;
   for (unsigned int i=0; i<rings.size(); i++) { 
      std::vector<int> ring = rings[i];
      if (circular_path.size() == ring.size()) {
	 bool jfound = 1;
	 for (unsigned int j=0; j<circular_path.size(); j++) { 
	    if (circular_path[j] != ring[j]) {
	       jfound = 0;
	       break;
	    }
	 }
	 if (jfound) {
	    ifound = 1;
	    break;
	 } 
      }
   }
   if (! ifound) {
      rings.push_back(circular_path);
   } 
} 


// pass filtered_rings
std::vector<std::vector<std::string> >
coot::aromatic_graph_t::indexes_to_names(const std::vector<std::vector<int> > &filtered_rings) const {

   std::vector<std::vector<std::string> > sv;
   for (unsigned int i=0; i<filtered_rings.size(); i++) {
      std::vector<std::string> string_ring(filtered_rings[i].size());
      for (unsigned int j=0; j<filtered_rings[i].size(); j++) {
	 string_ring[j] = vertices[filtered_rings[i][j]].get_name();
      }
      sv.push_back(string_ring);
   }
   return sv;
}
