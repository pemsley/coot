
#include <string>
#include <vector>
#include <algorithm>

#include "bonded-quad.hh"

namespace coot {

   class lbg_edge {
      int vertex_1;  // indices
      int vertex_2;
   public:
      lbg_edge(int v1, int v2) {
	 vertex_1 = v1;
	 vertex_2 = v2;
      }
      int get_vertex_index_1() const {
	 return vertex_1;
      } 
      int get_vertex_index_2() const {
	 return vertex_2;
      }
      friend std::ostream &operator<<(std::ostream &s, lbg_edge e);
   };
   std::ostream &operator<<(std::ostream &s, lbg_edge e);
   
   class lbg_vertex {
      std::string name;
      std::vector<int> edge_indices;
   public:
      lbg_vertex(const std::string &n) {
	 name = n;
      }
      bool match_name(const std::string &test_name) const {
	 return (name == test_name);
      }
      std::string get_name() const { return name; }
      std::vector<int> get_edges() const { return edge_indices; }
      void add_edge_index(int ei) { edge_indices.push_back(ei); }
   };


   class aromatic_graph_t {
      std::vector<std::pair<std::string, std::string> > original_bonds;
      std::vector<lbg_edge> edges;
      std::vector<lbg_vertex> vertices; // unfiltered
      void print() const;
      void print_path(std::vector<int> &path) const;
      std::vector<int> next_vertex(int start_vertex,
				   const std::vector<int> &path,
				   int depth, int this_vertex);
      std::vector<int> get_neighbours_of_vertex(int this_vertex) const;
      std::vector<int> get_neighbours_of_vertex_excluding_path(int this_vertex,
							       const std::vector<int> &path) const;
      void add_path_maybe(std::vector<int> circular_path);
      std::vector<std::vector<int> > rings; // converted from ints to string in ring_list()
      std::vector<std::vector<std::string> > indexes_to_names(const std::vector<std::vector<int> > &filtered_rings) const;
      bool
      has_same_elements_p(const std::vector<int> &test_ring,
                          const std::vector<std::vector<int> > &filtered_rings) const;

   public:
      aromatic_graph_t(const std::vector<std::pair<std::string, std::string> > &bonds_in);
      std::vector<std::vector<std::string> > ring_list();
      std::vector<bonded_quad_atom_names> bonded_quad_ring_list();
   };

}
