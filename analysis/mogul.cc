
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>

#include "coot-utils.hh"
#include "mogul-interface.hh"

void
coot::mogul::parse(const std::string &file_name) {

   max_z_badness = 5.0;
   if (! file_exists(file_name)) {
      std::cout << "filename " << file_name << " does not exist " << std::endl;
   } else {
      std::ifstream f(file_name.c_str());
      if (!f) {
	 std::cout << "Failed to open " << file_name << std::endl;
      } else {

	 std::vector<std::string> lines;
	 std::string line;
	 while (std::getline(f, line)) { 
	    lines.push_back(line);
	 }
	 std::cout << lines.size() << " lines." << std::endl;

	 // 2 from the end because we do a lines+2 to try to read distributions
	 for (unsigned int iline=0; iline<(lines.size()); iline++) {
	    std::vector<std::string> bits = coot::util::split_string(lines[iline], ",");

	    std::cout << "considering ----------------------------------------- "
		      << lines[iline] << std::endl;
	    
	    if (bits[0] == "BOND") {
	       try {
		  mogul_item item = parse_item_line(bits, 2);
		  items.push_back(item);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "WARNING:: " << rte.what() << std::endl;
	       } 
	    }
	    if (bits[0] == "ANGLE") {
	       try {
		  mogul_item item = parse_item_line(bits, 3);
		  items.push_back(item);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "WARNING:: " << rte.what() << std::endl;
	       } 
	    }
	    if (bits[0] == "TORSION") {
	       try {
		  mogul_item item = parse_item_line(bits, 4);
		  std::cout << "ACCEPTED:::::::::::::::::::: " << lines[iline] << std::endl;
		  items.push_back(item);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "WARNING:: " << rte.what() << std::endl;
		  std::cout << "REJECTED:::::::::::::::::::: " << lines[iline] << std::endl;
	       } 
	    }
	 }
      }
   }
}


// can throw a runtime_error
// 
coot::mogul_item
coot::mogul::parse_item_line(const std::vector<std::string> &bits, int n_idx) const {

   // Fragment Type,Atom Indices,Query Value,Hits,Mean,Median,Standard Deviation,|z-score|,d(min),Distribution Minimum,Distribution Maximum,Bin Width,Number of Bins

   mogul_item r;
   bool debug = true;

   if (bits.size() > 6) {
      std::string indices_string = bits[1];
      std::vector<int> indices = get_indices(indices_string);

      float mean    = 0;
      float median  = 0;
      float std_dev = 0;
      float z       = 0;
      float dmin    = 0;
      float min     = 0;
      float max     = 0;

      std::cout << "got " << indices.size() << " indices " << std::endl;
      if (indices.size() > 1) {
	 float model_value = coot::util::string_to_float(bits[2]);
	 int counts    = coot::util::string_to_int(  bits[3]);

	 if (indices.size() == 2 || indices.size() == 3) { 
	    mean    = coot::util::string_to_float(bits[4]);
	    median  = coot::util::string_to_float(bits[5]);
	    std_dev = coot::util::string_to_float(bits[6]);
	    z       = coot::util::string_to_float(bits[7]);
	    min     = coot::util::string_to_float(bits[9]);
	    max     = coot::util::string_to_float(bits[10]);
	 } else {
	    dmin    = coot::util::string_to_float(bits[8]);
	    std::cout << "got dmin " << dmin << " from :" << dmin << ":" << std::endl;
	 } 
	 
	 if (n_idx == 2) {
	    if (indices.size() == 2) { 
	       r = mogul_item(indices[0], indices[1], model_value,
			      counts, mean, median, std_dev, z);
	    }
	 }
	 if (n_idx == 3) {
	    if (indices.size() == 3) { 
	       r = mogul_item(indices[0], indices[1], indices[2],
			      model_value, counts, mean, median, std_dev, z);
	    }
	 }
	 if (n_idx == 4) {
	    if (indices.size() == 4) { 
	       r = mogul_item(indices[0], indices[1], indices[2], indices[3],
			      model_value, counts, dmin);
	    }
	 }
	 r.set_max_z_badness(max_z_badness);

	 if (bits.size() > 11) {
	    std::vector<std::string> distribution_bits;
	    for (unsigned int i=9; i<bits.size(); i++)
	       distribution_bits.push_back(bits[i]);
	    r.distribution = mogul_distribution(distribution_bits);
	 }

	 if (debug)
	    std::cout << " item " << indices[0] << " " << indices[1] << " mean: " 
		      << mean << " min: " << min << " max: " << max << " median: "
		      << median << " std: " << std_dev
		      << "  Z-score:" << z << std::endl;
      }
   }
   return r;
}

std::vector<int>
coot::mogul::get_indices(const std::string &indices_string) const { 

   std::vector<int> v;
   std::vector<std::string> idx_bits = coot::util::split_string_no_blanks(indices_string, " ");
   for (unsigned int i=0; i<idx_bits.size(); i++)
      v.push_back(coot::util::string_to_int(idx_bits[i]));
   return v;
}

coot::mogul_distribution::mogul_distribution(const std::vector<std::string> &bits) {
   
   if (bits.size() > 6) {
      // lower-bottom and upper-top

      std::cout << "Geting bin stuff from here "
		<< bits[0] << "  "
		<< bits[1] << "  "
		<< bits[2] << "  "
		<< bits[3] << "  "
		<< bits[4] << "  "
		<< bits[5] << "  "
		<< std::endl;
      bin_start  = coot::util::string_to_float(bits[0]);
      bin_end    = coot::util::string_to_float(bits[1]);
      bin_width  = coot::util::string_to_float(bits[2]);
      n_bins     = coot::util::string_to_float(bits[3]);

      if (1) {
	 std::cout << "bin_start " << bin_start << std::endl;
	 std::cout << "bin_end   " << bin_end   << std::endl;
	 std::cout << "bin_width " << bin_width << std::endl;
	 std::cout << "n_bins    " << n_bins    << std::endl;
      }
      
      for (unsigned int ibin=0; ibin<n_bins; ibin++) {
	 int ibit = ibin + 5;
	 if (ibit < bits.size()) {
	    int v = coot::util::string_to_int(bits[ibit]);
	    counts.push_back(v);
	 }
      }
   }
} 

// can throw a runtime_error
// 
coot::mogul_item
coot::mogul::parse_angle_line(const std::vector<std::string> &bond_bits,
			      const std::vector<std::string> &stats_bits,
			      const std::vector<std::string> &distribution_bits) const {

   mogul_item r;
   if (bond_bits.size() > 6) {
      int atom_idx_1 = coot::util::string_to_int(bond_bits[1]);
      int atom_idx_2 = coot::util::string_to_int(bond_bits[2]);
      int atom_idx_3 = coot::util::string_to_int(bond_bits[3]);
      float model_value = coot::util::string_to_float(bond_bits[4]);
      if (stats_bits[0] == "STATS") {
	 int counts    =  coot::util::string_to_int(  stats_bits[1]);
	 float mean    =  coot::util::string_to_float(stats_bits[2]);
	 float min     =  coot::util::string_to_float(stats_bits[3]);
	 float max     =  coot::util::string_to_float(stats_bits[4]);
	 float median  =  coot::util::string_to_float(stats_bits[5]);
	 float std_dev =  coot::util::string_to_float(stats_bits[6]);
	 float z = fabsf((model_value - median)/std_dev);
	 r = mogul_item(atom_idx_1, atom_idx_2, atom_idx_3, model_value,
			counts, mean, median, std_dev, z);
	 r.set_max_z_badness(max_z_badness);
	 if (0)
	    std::cout << " angle " << atom_idx_1 << " " << atom_idx_2 << " mean: " 
		      << mean << " min: " << min << " max: " << max << " median: "
		      << median << " std: " << std_dev
		      << "  Z-score:" << z 
		      << std::endl;
	 r.distribution = mogul_distribution(distribution_bits);
      }
   }
   return r;
}


std::string
coot::mogul_item::colour() const {

   colour_holder ch;
   float max_z = max_badness;
   float min_z = 1.0;

   float this_z = z;
   float range = max_z - min_z;
   float f = (this_z-min_z)/range;
   if (f > 1.0) f = 1.0;
   if (f < 0.0) f = 0.0;

   ch.blue = 0.0;
   ch.red = f;
   ch.green = 1 - f; 

   return ch.hex();
} 


coot::dictionary_residue_restraints_t
coot::mogul::make_restraints(CResidue *residue_p,
			     const std::string &comp_id,
			     const coot::protein_geometry &geom) {

   coot::dictionary_residue_restraints_t r(comp_id, -1);
   if (residue_p) { 
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::pair<bool, dictionary_residue_restraints_t> current_restraints = 
	 geom.get_monomer_restraints(comp_id);

      for (unsigned int i=0; i<items.size(); i++) { 
	 if (items[i].type == mogul_item::BOND) {
	    int idx_1 = items[i].idx_1 - 1;
	    int idx_2 = items[i].idx_2 - 1;
	    if (idx_1 >= 0 && idx_1 < n_residue_atoms) { 
	       if (idx_2 >= 0 && idx_2 < n_residue_atoms) {
		  std::string name_1(residue_atoms[idx_1]->name);
		  std::string name_2(residue_atoms[idx_2]->name);
		  std::string type;
		  if (current_restraints.first) 
		     type = get_bond_type(current_restraints.second, name_1, name_2);
		  float dist = items[i].median;
		  float esd  = items[i].std_dev;
		  dict_bond_restraint_t rest(name_1, name_2, type, dist, esd);
		  r.bond_restraint.push_back(rest);
	       }
	    }
	 }


	 if (items[i].type == mogul_item::ANGLE) {
	    int idx_1 = items[i].idx_1 - 1;
	    int idx_2 = items[i].idx_2 - 1;
	    int idx_3 = items[i].idx_3 - 1;
	    if (idx_1 >= 0 && idx_1 < n_residue_atoms) { 
	       if (idx_2 >= 0 && idx_2 < n_residue_atoms) {
		  if (idx_3 >= 0 && idx_3 < n_residue_atoms) {
		     std::string name_1(residue_atoms[idx_1]->name);
		     std::string name_2(residue_atoms[idx_2]->name);
		     std::string name_3(residue_atoms[idx_3]->name);
		     float angle = items[i].median;
		     float esd  = items[i].std_dev;
		     dict_angle_restraint_t rest(name_1, name_2, name_3, angle, esd);
		     r.angle_restraint.push_back(rest);
		  }
	       }
	    }
	 }
      }
   }
   return r;
}

coot::dictionary_residue_restraints_t
coot::mogul::make_restraints(const std::string &comp_id,
			     const std::string &compound_name,
			     const std::vector<std::string> &atom_names,
			     int n_atoms_all, int n_atoms_non_H,
			     const coot::dictionary_residue_restraints_t &bond_types_dict) {

   coot::dictionary_residue_restraints_t r(comp_id, -1);

   r.residue_info.comp_id = comp_id;
   r.residue_info.three_letter_code = comp_id;
   r.residue_info.name = compound_name;
   r.residue_info.number_atoms_all = n_atoms_all;
   r.residue_info.number_atoms_nh = n_atoms_non_H;
   r.residue_info.group = "non-polymer";
   r.residue_info.description_level = "Partial";

   
   if (! atom_names.empty()) {
      for (unsigned int i=0; i<items.size(); i++) { 
	 if (items[i].type == mogul_item::BOND) {
	    int idx_1 = items[i].idx_1 - 1;
	    int idx_2 = items[i].idx_2 - 1;
	    if (idx_1 >= 0 && idx_1 < atom_names.size()) { 
	       if (idx_2 >= 0 && idx_2 < atom_names.size()) {
		  std::string name_1 = atom_names[idx_1];
		  std::string name_2 = atom_names[idx_2];
		  // bt is "" if not found.
		  std::string bt = bond_types_dict.get_bond_type(name_1, name_2);
		  float dist = items[i].median;
		  float esd  = items[i].std_dev;
		  dict_bond_restraint_t rest(name_1, name_2, bt, dist, esd);
		  r.bond_restraint.push_back(rest);
	       }
	    }
	 }


	 if (items[i].type == mogul_item::ANGLE) {
	    int idx_1 = items[i].idx_1 - 1;
	    int idx_2 = items[i].idx_2 - 1;
	    int idx_3 = items[i].idx_3 - 1;
	    if (idx_1 >= 0 && idx_1 < atom_names.size()) { 
	       if (idx_2 >= 0 && idx_2 < atom_names.size()) {
		  if (idx_3 >= 0 && idx_3 < atom_names.size()) {
		     std::string name_1 = atom_names[idx_1];
		     std::string name_2 = atom_names[idx_2];
		     std::string name_3 = atom_names[idx_3];
		     float angle = items[i].median;
		     float esd  = items[i].std_dev;
		     dict_angle_restraint_t rest(name_1, name_2, name_3, angle, esd);
		     r.angle_restraint.push_back(rest);
		  }
	       }
	    }
	 }
      }

   }
   return r;
} 



std::string
coot::mogul::get_bond_type(const coot::dictionary_residue_restraints_t &restraints,
			   const std::string &name_1,
			   const std::string &name_2) const {

   return restraints.get_bond_type(name_1, name_2);
} 


// can throw a runtime_exception
coot::mogul_item
coot::mogul::get_angle_item(const std::vector<int> &indices) const {

   if (indices.size() != 3) {
      throw(std::runtime_error("wrong size of indices"));
   } else {
      for (unsigned int j=0; j<items.size(); j++) { 
	 const mogul_item &item = items[j];
	 if (item.matches_indices(indices)) {
	    return item;
	 }
      }
   }
   throw(std::runtime_error("no such item"));
}

// can throw a runtime_exception
coot::mogul_item
coot::mogul::get_bond_item(const std::vector<int> &indices) const {

   if (indices.size() != 2) {
      throw(std::runtime_error("wrong size of indices"));
   } else {
      for (unsigned int j=0; j<items.size(); j++) { 
	 const mogul_item &item = items[j];
	 if (item.matches_indices(indices)) {
	    return item;
	 }
      }
   }
   throw(std::runtime_error("no such item"));
}

// can throw a runtime_exception
coot::mogul_item
coot::mogul::get_torsion_item(const std::vector<int> &indices) const {

   if (indices.size() != 4) {
      throw(std::runtime_error("wrong size of indices"));
   } else {
      for (unsigned int j=0; j<items.size(); j++) { 
	 const mogul_item &item = items[j];
	 if (item.matches_indices(indices)) {
	    return item;
	 }
      }
   }
   throw(std::runtime_error("no such item"));
}



bool
coot::mogul_item::matches_indices(const std::vector<int> &indices) const {

   if (indices.size() == 4) {
      if (indices[0] == idx_1) { 
	 if (indices[1] == idx_2) { 
	    if (indices[2] == idx_3) {
	       if (indices[3] == idx_4) {
		  return true;
	       }
	    }
	 }
      }
   } 
   if (indices.size() == 3) {
      if (indices[0] == idx_1) { 
	 if (indices[1] == idx_2) { 
	    if (indices[2] == idx_3) {
	       return true;
	    }
	 }
      }
   } 
   if (indices.size() == 2) {
      if (indices[0] == idx_1) { 
	 if (indices[1] == idx_2) {
	    return true;
	 }
      }
   }
   return false;
} 
