
#include <fstream>
#include "protein-geometry.hh"

#include "analysis/stats.hh"

void
write_bond_lengths(const std::string &file_name,
		   const std::vector<double> &d) {

   std::ofstream f(file_name.c_str());

   if (f) {
      for (unsigned int i=0; i<d.size(); i++) { 
	 f << d[i] << "\n";
      }
   }
   f.close();
}

int main(int argc, char **argv) {

   int r = 0;

   coot::protein_geometry geom;
   geom.set_verbose(false);

   int read_number = 0;
   for (unsigned int i=1; i<argc; i++) {
      std::string file_name = argv[i];
      // std::cout << "file_name " << file_name << std::endl;
      int status = geom.init_refmac_mon_lib(file_name, read_number);
      read_number++;
   }

   std::map<std::pair<std::string, std::string>, std::vector<double> > bonds;

   for (unsigned int i=0; i<geom.size(); i++) { 
      const coot::dictionary_residue_restraints_t r = geom.get_monomer_restraints(i);
      for (unsigned int ibond=0; ibond<r.bond_restraint.size(); ibond++) { 
	 const coot::dict_bond_restraint_t &br = r.bond_restraint[ibond];
	 if (!(r.is_hydrogen(br.atom_id_1_4c())) &&
	     !(r.is_hydrogen(br.atom_id_2_4c()))) { 
	    std::pair<std::string, std::string> p(r.type_energy(br.atom_id_1_4c()),
						  r.type_energy(br.atom_id_2_4c()));
	    bonds[p].push_back(br.value_dist());
	 }
      }
   }

   std::map<std::pair<std::string, std::string>, std::vector<double> >::const_iterator it;
   unsigned int min_counts = 250;
   
   for (it=bonds.begin(); it!=bonds.end(); it++) {
      if (it->second.size() > min_counts) {
	 
	 // std::cout << it->first.first << " " << it->first.second << " "
	 // << it->second.size() << std::endl;

	 std::string stub = it->first.first + "-" + it->first.second;
	 std::string file_name = stub + ".tab";
	 write_bond_lengths(file_name, it->second);

	 coot::stats::single stats;
	 for (unsigned int i=0; i<it->second.size(); i++)
	    stats.add(it->second[i]);

	 double var = stats.variance();
	 std::cout << stub << "   N: " << stats.size() << " mean: "
		   << stats.mean() << "   v: " << var << "   s: "
		   << sqrt(var) << "   k: "
		   << stats.kurtosis()  << std::endl;

	 
      }
   }

   return r;
} 
