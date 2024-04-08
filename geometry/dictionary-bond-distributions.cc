/* geometry/dictionary-bond-distributions.cc
 * 
 * Copyright 2014 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

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
   energy_type_key_t() { }
   bool operator<(const energy_type_key_t &k) const {
      if (k.bond_type < bond_type) { 
	 return true;
      } else {
	 if (k.bond_type == bond_type) { 
	    if (k.energy_type_1 < energy_type_1) {
	       return true;
	    } else {
	       if (k.energy_type_1 == energy_type_1) { 
		  if (k.energy_type_2 < energy_type_2)
		     return true;
	       }
	    }
	 }
      }
      return false;
   }
   bool operator==(const energy_type_key_t &k) const {
      if (k.bond_type == bond_type) {
	 if (k.energy_type_1 == energy_type_1) {
	    if (k.energy_type_2 == energy_type_2) {
	       return true;
	    }
	 }
      }
      return false;
   }
   
};

class bond_dist_info_t {

public:
   double dist;
   std::string comp_id;
   std::string atom_name_1;
   std::string atom_name_2;
   bond_dist_info_t(const double &d,
		    const std::string &c,
		    const std::string &a1,
		    const std::string &a2) {
      dist = d;
      comp_id = c;
      atom_name_1 = a1;
      atom_name_2 = a2;
   }
};

bool bonds_vec_sorter(const std::pair<energy_type_key_t, std::vector<bond_dist_info_t> > &k1,
		      const std::pair<energy_type_key_t, std::vector<bond_dist_info_t> > &k2) {


   return (k1.second.size() < k2.second.size()); // sort by counts

}

bool bonds_vec_k_sorter(const std::pair<energy_type_key_t, std::vector<bond_dist_info_t> > &k1,
			const std::pair<energy_type_key_t, std::vector<bond_dist_info_t> > &k2) {


   std::vector<double> v1(k1.second.size());
   std::vector<double> v2(k2.second.size());
   for (unsigned int i=0; i<k1.second.size(); i++) v1[i] = k1.second[i].dist;
   for (unsigned int i=0; i<k2.second.size(); i++) v2[i] = k2.second[i].dist;
   
   coot::stats::single s1(v1); 
   coot::stats::single s2(v2);

   if (0) 
      std::cout << "on sorting kurtosises: "
		<< " key sizes: " << k1.second.size() << " " << k2.second.size() << "   "
		<< k1.first.energy_type_1 << "-" << k1.first.energy_type_2
		<< " " << k1.first.bond_type
		<< " " << k1.second.size()
		<< " " << s1.mean()
		<< " " << s1.kurtosis() << " "
		<< "     ======     " << k2.first.energy_type_1 << "-" << k2.first.energy_type_2
		<< " " << k2.first.bond_type
		<< " " << k2.second.size()
		<< " " << s2.mean() 
		<< " " << s2.kurtosis()
		<< "\n";

   return (s1.kurtosis() < s2.kurtosis());
} 

void
write_bond_lengths(const std::string &file_name,
		   const std::vector<bond_dist_info_t> &d) {

   std::ofstream f(file_name.c_str());

   if (f) {
      for (unsigned int i=0; i<d.size(); i++) { 
	 f << d[i].dist << " " << d[i].comp_id << " " << d[i].atom_name_1 << " " << d[i].atom_name_2 << "\n";
      }
   }
   f.close();
}

int main(int argc, char **argv) {

   int r = 0;

   coot::protein_geometry geom;
   geom.set_verbose(false);

   int read_number = 0;
   for (int i=1; i<argc; i++) {
      std::string file_name = argv[i];
      coot::read_refmac_mon_lib_info_t rmi = geom.init_refmac_mon_lib(file_name, read_number);
      std::cout << "read " << file_name << " with status " << rmi.success << std::endl;
      read_number++;
   }

   std::map<energy_type_key_t, std::vector<bond_dist_info_t> > bonds_map;
   std::map<energy_type_key_t, std::vector<bond_dist_info_t> >::iterator it;

   if (1)
      std::cout << "Examining the bonds in the " << geom.size() << " entries."<< std::endl;
   
   for (unsigned int i=0; i<geom.size(); i++) {
      const coot::dictionary_residue_restraints_t &r = geom.get_monomer_restraints(i);
      std::cout << "adding to map bonds in " << r.residue_info.comp_id << std::endl;
      for (unsigned int ibond=0; ibond<r.bond_restraint.size(); ibond++) { 
	 const coot::dict_bond_restraint_t &br = r.bond_restraint[ibond];
	 std::string an1 = br.atom_id_1_4c();
	 std::string an2 = br.atom_id_2_4c();
	 if (!(r.is_hydrogen(an1)) &&
	     !(r.is_hydrogen(an2))) { 
	    energy_type_key_t k(r.type_energy(an1),
				r.type_energy(an2),
				br.type());
	    bond_dist_info_t bi(br.value_dist(), r.residue_info.comp_id, an1, an2);

	    // Doing this directly tickles a compiler bug?
	    // bonds[k].push_back(bi);
	    // so let's use the iterator:

	    it = bonds_map.find(k);
	    if (it != bonds_map.end()) { 
	       it->second.push_back(bi);
	    } else {
	       std::vector<bond_dist_info_t> v;
	       v.push_back(bi);
	       bonds_map[k] = v;
	    } 
	 }
      }
   }

   // convert to vector for sorting by N (counts)
   //
   std::vector<std::pair<energy_type_key_t, std::vector<bond_dist_info_t> > > bonds_vec;
   std::vector<std::pair<energy_type_key_t, std::vector<bond_dist_info_t> > >::const_iterator itv;

   for (it=bonds_map.begin(); it!=bonds_map.end(); it++) {
      std::pair<energy_type_key_t, std::vector<bond_dist_info_t> > p(it->first, it->second);
      bonds_vec.push_back(p);
   }

   // now sort bonds_vec how you like:
   std::sort(bonds_vec.begin(), bonds_vec.end(), bonds_vec_k_sorter);
   

   unsigned int min_counts = 0;
   unsigned int n_found = 0;

   for (itv=bonds_vec.begin(); itv!=bonds_vec.end(); itv++) {
      
      if (itv->second.size() >= min_counts) {

	 n_found++;
	 
	 std::string stub = itv->first.energy_type_1 + "-" + itv->first.energy_type_2;

	 std::string file_name = stub + "-" + itv->first.bond_type + ".tab";
	 write_bond_lengths(file_name, itv->second);

	 coot::stats::single stats;
	 for (unsigned int i=0; i<itv->second.size(); i++)
	    stats.add(itv->second[i].dist);
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

   if (0)
      std::cout << "Found " << n_found << " bond types with more than " << min_counts
		<< " entries"<< std::endl;

   return r;
} 

