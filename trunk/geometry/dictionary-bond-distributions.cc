
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "protein-geometry.hh"

#include "analysis/stats.hh"

class energy_type_key_t {

public:
   std::string energy_type_1; 
   std::string energy_type_2;
   std::string bond_type;
   energy_type_key_t(const std::string &e1,
		     const std::string &e2,
		     const std::string &bond_type_in) {
      bond_type = bond_type_in;
      energy_type_1 = e1;
      energy_type_2 = e2;
      if (energy_type_2 < energy_type_1) {
	 energy_type_1 = e2;
	 energy_type_2 = e1;
      }
   }
   energy_type_key_t() {} // needed?
   bool operator<(const energy_type_key_t &k) const {
      if (k.bond_type < bond_type) { 
	 return true;
      } else {
	 if (k.energy_type_2 < energy_type_2) {
	    return true;
	 } else {
	    return (k.energy_type_1 < energy_type_1);
	 }
      } 
   }
};

bool bonds_vec_sorter(const std::pair<energy_type_key_t, std::vector<double> > &k1,
		      const std::pair<energy_type_key_t, std::vector<double> > &k2) {
   return (k1.second.size() < k2.second.size());
} 

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
      int status = geom.init_refmac_mon_lib(file_name, read_number);
      // std::cout << "read " << file_name << " with status " << status << std::endl;
      read_number++;
   }

   // std::map<std::pair<std::string, std::string>, std::vector<double> > bonds;
   std::map<energy_type_key_t, std::vector<double> > bonds;

   std::cout << "Examining the bonds in the " << geom.size() << " entries."<< std::endl;
   for (unsigned int i=0; i<geom.size(); i++) {
      const coot::dictionary_residue_restraints_t r = geom.get_monomer_restraints(i);
      for (unsigned int ibond=0; ibond<r.bond_restraint.size(); ibond++) { 
	 const coot::dict_bond_restraint_t &br = r.bond_restraint[ibond];
	 if (!(r.is_hydrogen(br.atom_id_1_4c())) &&
	     !(r.is_hydrogen(br.atom_id_2_4c()))) { 
	    energy_type_key_t k(r.type_energy(br.atom_id_1_4c()),
				r.type_energy(br.atom_id_2_4c()),
				br.type());
	    bonds[k].push_back(br.value_dist());
	 }
      }
   }

   // convert to vector for sorting by N (counts)
   //
   std::map<energy_type_key_t, std::vector<double> >::const_iterator it;
   std::vector<std::pair<energy_type_key_t, std::vector<double> > > bonds_vec(bonds.size());
   std::vector<std::pair<energy_type_key_t, std::vector<double> > >::const_iterator itv;

   for (it=bonds.begin(); it!=bonds.end(); it++) {
      std::pair<energy_type_key_t, std::vector<double> > p(it->first, it->second);
      bonds_vec.push_back(p);
   }

   // now sort bonds_vec how you like:
   std::sort(bonds_vec.begin(), bonds_vec.end(), bonds_vec_sorter);
   

   unsigned int min_counts = 10;
   unsigned int n_found = 0;

   for (itv=bonds_vec.begin(); itv!=bonds_vec.end(); itv++) {
      
      if (itv->second.size() > min_counts) {

	 n_found++;
	 
	 std::string stub = itv->first.energy_type_1 + "-" + itv->first.energy_type_2;

	 std::string file_name = stub + ".tab";
	 write_bond_lengths(file_name, itv->second);

	 coot::stats::single stats;
	 for (unsigned int i=0; i<itv->second.size(); i++)
	    stats.add(itv->second[i]);
	 double var = stats.variance();


	 // Formatted output:
	 //
	 // For float/double numbers
	 // Using std::fixed the precision is the number of decimal places
	 // use std::right << set::precision(3) << std::fixed << f
	 // You may need to set std::left at the beginning of the line (typically for strings)
	 // To reset from std::fixed:
	 // std::cout.setf(std::ios::fixed, std::ios::floatfield);
	 
	 std::cout << std::left;
	 std::cout << std::setprecision(7);
	 std::cout.setf(std::ios::fixed, std::ios::floatfield);
	 
	 std::cout << std::setw(9) << stub
		   << "  "       << std::setw(6) << itv->first.bond_type
		   << "   N: "   << std::setw(4) << std::right << stats.size()
		   << "  mean: " << std::setw(8) << stats.mean()
		   << "   s: "   << std::setw(9) << sqrt(var)
		   << "   k: "   << std::setw(7) << std::right << std::setprecision(3)
		   << std::fixed << stats.kurtosis() << std::endl;
      }
   }

   std::cout << "Found " << n_found << " bond types with more than " << min_counts
	     << " entries"<< std::endl;

   return r;
} 

