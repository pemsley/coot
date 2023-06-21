
#include <fstream>
#include "utils/coot-utils.hh"
#include "protein-geometry.hh"
#include "src/json.hpp"  // thanks Niels Lohmann.
using json = nlohmann::json;

std::map<std::string, std::vector<metal_ligand_t> > metals_store;

void
coot::protein_geometry::print_metal_store() const {

   for(const auto &kd : metals_store) {
      const std::string &key = kd.first;
      std::cout << "----- " << key << " -------" << std::endl;
      for (unsigned int i=0; i<kd.second.size(); i++) {
         const metal_ligand_t &ml = kd.second[i];
         std::cout << "   " << ml.element << std::endl;
         for (unsigned int j=0; j<ml.coordinated_atoms.size(); j++) {
            const coordinated_atom_t &ca = ml.coordinated_atoms[j];
            std::cout << "      " << ca.coordination_number << " " << ca.median << " " << ca.mad
                      << " " << ca.mean << " " << ca.std << " " << ca.count << std::endl;
         }
      }
   }

}

void
coot::protein_geometry::read_metal_distances(const std::string &file_name) {

   if (file_exists(file_name)) {

      // file to string
      std::fstream f(file_name);
      std::string s;
      f.seekg(0, std::ios::end);
      s.reserve(f.tellg());
      f.seekg(0, std::ios::beg);
      s.assign((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

      try {
         json j = json::parse(s);
         json ls = j["metal_coordination"];

         for (json::iterator it_1=ls.begin(); it_1!=ls.end(); ++it_1) {
            const std::string &metal_ele = it_1.key();
            std::cout << "metal_ele " << metal_ele << std::endl;
            json ls_2 = ls[metal_ele];
            for (json::iterator it_2=ls_2.begin(); it_2!=ls_2.end(); ++it_2) {
               const std::string &ligand_ele = it_2.key();
               metal_ligand_t metal_ligand(ligand_ele); // no, not ele here
               json &j_coord_atom_list = *it_2;
               unsigned int size = j_coord_atom_list.size();
               for (unsigned int i=0; i<size; i++) {
                  json &j_coord_atom = j_coord_atom_list[i];
                  int coord = -1; double median = -1; double mad = -1;
                  double mean = -1; double std = -1; int count = -1;
                  json::iterator it;
                  it = j_coord_atom.find(std::string("coord"));
                  if (it != j_coord_atom.end()) { coord = it.value(); }
                  it = j_coord_atom.find(std::string("median"));
                  if (it != j_coord_atom.end()) { median = it.value(); }
                  it = j_coord_atom.find(std::string("mad"));
                  if (it != j_coord_atom.end()) { mad = it.value(); }
                  it = j_coord_atom.find(std::string("mean"));
                  if (it != j_coord_atom.end()) { mean = it.value(); }
                  it = j_coord_atom.find(std::string("std"));
                  if (it != j_coord_atom.end()) { std = it.value(); }
                  it = j_coord_atom.find(std::string("count"));
                  if (it != j_coord_atom.end()) { count = it.value(); }

                  coordinated_atom_t ca(coord, median, mad, mean, std, count);
                  metal_ligand.add(ca);
               }
               metals_store[metal_ele].push_back(metal_ligand);
            }
         }
      }

      catch(const nlohmann::detail::type_error &e) {
         std::cout << "ERROR:: " << e.what() << std::endl;
      }
      catch(const nlohmann::detail::parse_error &e) {
         std::cout << "ERROR:: " << e.what() << std::endl;
      }
   }

   print_metal_store();
}
