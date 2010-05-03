
#include <iostream>
#include "lbg-graph.hh"

int main(int argc, char **argv) {

   std::vector<std::pair<std::string, std::string> > bonds;

    bonds.push_back(std::pair<std::string,std::string>(" C1 "," C2 "));
    bonds.push_back(std::pair<std::string,std::string>(" C2 "," C3 "));
    bonds.push_back(std::pair<std::string,std::string>(" C2 "," H2 "));
    bonds.push_back(std::pair<std::string,std::string>(" C3 "," H3 "));
    bonds.push_back(std::pair<std::string,std::string>(" C4 "," C3 "));
    bonds.push_back(std::pair<std::string,std::string>(" C4 "," C5 "));
    bonds.push_back(std::pair<std::string,std::string>(" C6 "," C5 "));
    bonds.push_back(std::pair<std::string,std::string>(" C6 "," C1 "));
    bonds.push_back(std::pair<std::string,std::string>(" C6 "," H6 "));
    bonds.push_back(std::pair<std::string,std::string>(" C5 "," C7 "));
    bonds.push_back(std::pair<std::string,std::string>(" C7 "," C8 "));
    bonds.push_back(std::pair<std::string,std::string>(" C8 "," C9 "));
    bonds.push_back(std::pair<std::string,std::string>(" C9 "," C6 "));

//    bonds.push_back(std::pair<std::string,std::string>(" C1 "," C2 "));
//    bonds.push_back(std::pair<std::string,std::string>(" C2 "," C3 "));
//    bonds.push_back(std::pair<std::string,std::string>(" C3 "," C4 "));
//    bonds.push_back(std::pair<std::string,std::string>(" C4 "," C1 "));
   
   coot::aromatic_graph_t graph(bonds);

   std::vector<std::vector<std::string> > ring_list = graph.ring_list();

   std::cout << "----------- " << ring_list.size() << " rings ---------- " << std::endl;
   for (unsigned int i=0; i<ring_list.size(); i++) {
      std::cout << "ring " << i << "\n   ";
      for (unsigned int j=0; j<ring_list[i].size(); j++) { 
	 std::cout << ring_list[i][j] << "  ";
      }
      std::cout << std::endl;
   }

   return 0;
}
