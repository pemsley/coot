
#ifndef MAKE_ENHANCED_LIGAND_TOOLS
int main(int argc, char **argv) {return 0;}
#else 
#include "cod-types.hh"

#include <map>
#include <algorithm>

#include "geometry/protein-geometry.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "rdkit-interface.hh"

// ordered so that the table is written in decreasing order
bool string_int_pair_sorter(const std::pair<std::string, unsigned int> &p1,
			    const std::pair<std::string, unsigned int> &p2) {

   return (p2.second < p1.second);
} 

int find_string_in_vector(const std::vector<std::pair<std::string, unsigned int> > &v,
			  const std::string &t) {

   int r = -1;
   for (unsigned int i=0; i<v.size(); i++) {
      if (v[i].first == t) {
	 return i;
      } 
   }
   return r;
} 


int main(int argc, char **argv) {

   int r = 0;

   coot::protein_geometry geom;
   geom.set_verbose(false);

   int read_number = 0;
   for (unsigned int i=1; i<argc; i++) {
      std::string file_name = argv[i];
      int status = geom.init_refmac_mon_lib(file_name, read_number);
      read_number++;
   }


   if (1)
      std::cout << "Examining the atoms in the " << geom.size() << " entries."<< std::endl;

   // key: refmac enerygy type
   // data: vector of COD atom type pairs (type and n occurances)
   std::map<std::string, std::vector<std::pair<std::string, unsigned int> > > atom_map;
   std::map<std::string, std::vector<std::pair<std::string, unsigned int> > >::const_iterator it;

   for (unsigned int i=0; i<geom.size(); i++) {
      const coot::dictionary_residue_restraints_t &r = geom.get_monomer_restraints(i);
      // std::cout << "adding atoms of " << r.residue_info.comp_id << " to atoms map" << std::endl;

      const std::string &comp_id = r.residue_info.comp_id;
      bool idealised_flag = true;
      mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, idealised_flag);

      if (! mol) {
	 std::cout << "Null mol from mol_from_dictionary() for " <<  comp_id << std::endl;
      } else {
	 
	 mmdb::Residue *residue_p = coot::util::get_first_residue(mol);

	 if (! residue_p) {
	    // pretty strange
	    std::cout << "Null residue from mol from mol_from_dictionary() for "
		      << comp_id << std::endl;
	 } else { 

	    try { 
	       RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, geom);
	       std::vector<std::string> v = cod::get_cod_atom_types(rdkm);
	       if (v.size() == r.atom_info.size()) {
		  // std::cout << comp_id << " was good" << std::endl;
		  for (unsigned int iat=0; iat<r.atom_info.size(); iat++) {
		     const std::string &key = r.atom_info[iat].type_energy;
		     if (0) 
			std::cout << comp_id << "  "
				  << std::setw(4) << r.atom_info[iat].atom_id << "  "
				  << std::setw(6) << r.atom_info[iat].type_energy << "   " 
				  << v[iat] << "\n";

		     it = atom_map.find(key);
		     if (it == atom_map.end()) {
			std::pair<std::string, unsigned int> p(v[iat], 1);
			atom_map[key].push_back(p);
		     } else { 
			int idx = find_string_in_vector(atom_map[key], v[iat]);
			if (idx == -1) { // not found
			   std::pair<std::string, unsigned int> p(v[iat], 1);
			   atom_map[key].push_back(p);
			} else {
			   atom_map[key][idx].second++;
			}
		     }
		  }
	       } else {
		  std::cout << comp_id << " was not good" << v.size() << " " << r.atom_info.size()
			    << std::endl;
	       }
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "error:: rte " << rte.what() << " " << comp_id << std::endl;
	    } 
	    catch (const std::exception &e) {
	       std::cout << "error:: exception " << e.what() << " " << comp_id << std::endl;
	    } 
	 }
      }
   }

   std::map<std::string, std::vector<std::pair<std::string, unsigned int> > >::iterator itv;
   for (itv=atom_map.begin(); itv!=atom_map.end(); itv++) {
      std::sort(itv->second.begin(), itv->second.end(), string_int_pair_sorter);
      std::cout << itv->first << "       ";
      for (unsigned int i=0; i<itv->second.size(); i++) { 
	 std::cout << "   " << itv->second[i].first  << " " << itv->second[i].second;
      }
      std::cout << "\n";

   } 

   
   return r;
}

#endif // ENHANCED_LIGAND_TOOLS
