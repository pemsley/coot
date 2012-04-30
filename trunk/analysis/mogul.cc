
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

	 for (unsigned int iline=0; iline<lines.size(); iline++) { 
	    std::vector<std::string> bits = coot::util::split_string_no_blanks(lines[iline], " ");

	    if (bits[0] == "BOND") {
	       try {
		  std::vector<std::string> stats_bits =
		     coot::util::split_string_no_blanks(lines[iline+1], " ");
		  std::vector<std::string> distribution_bits =
		     coot::util::split_string_no_blanks(lines[iline+2], " ");
		  mogul_item item = parse_bond_line(bits, stats_bits, distribution_bits);
		  items.push_back(item);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "WARNING:: " << rte.what() << std::endl;
	       } 
	    }
	    if (bits[0] == "ANGLE") {
	       try {
		  std::vector<std::string> stats_bits =
		     coot::util::split_string_no_blanks(lines[iline+1], " ");
		  std::vector<std::string> distribution_bits =
		     coot::util::split_string_no_blanks(lines[iline+2], " ");
		  mogul_item item = parse_angle_line(bits, stats_bits, distribution_bits);
		  items.push_back(item);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "WARNING:: " << rte.what() << std::endl;
	       } 
	    }
	 }
      }
   } 
} 

// can throw a runtime_error
// 
coot::mogul_item
coot::mogul::parse_bond_line(const std::vector<std::string> &bond_bits,
			     const std::vector<std::string> &stats_bits,
			     const std::vector<std::string> &distribution_bits) const {
   
   mogul_item r;
   bool debug = false;

   if (bond_bits.size() > 6) {
      int atom_idx_1 = coot::util::string_to_int(bond_bits[1]);
      int atom_idx_2 = coot::util::string_to_int(bond_bits[2]);
      float model_value = coot::util::string_to_float(bond_bits[3]);
      if (stats_bits[0] == "STATS") {
	 int counts    =  coot::util::string_to_int(  stats_bits[1]);
	 float mean    =  coot::util::string_to_float(stats_bits[2]);
	 float min     =  coot::util::string_to_float(stats_bits[3]);
	 float max     =  coot::util::string_to_float(stats_bits[4]);
	 float median  =  coot::util::string_to_float(stats_bits[5]);
	 float std_dev =  coot::util::string_to_float(stats_bits[6]);
	 float z = fabsf((model_value - median)/std_dev);
	 r = mogul_item(atom_idx_1, atom_idx_2, model_value,
			counts, mean, median, std_dev, z);
	 r.set_max_z_badness(max_z_badness);
	 if (debug)
	    std::cout << " bond " << atom_idx_1 << " " << atom_idx_2 << " mean: " 
		      << mean << " min: " << min << " max: " << max << " median: "
		      << median << " std: " << std_dev
		      << "  Z-score:" << z << std::endl;
	 mogul_distribution md(distribution_bits);
	 // only add distributions if the numbers were sane
	 if (md.counts.size() == md.n_bins) {
	    r.add_distribution(md);
	 }
      }
   }
   return r;
}

coot::mogul_distribution::mogul_distribution(const std::vector<std::string> &bits) {
   
   if (bits[0] == "DISTRIBUTION") {
      if (bits.size() > 6) {
	 // lower-bottom and upper-top
	 bin_start  = coot::util::string_to_float(bits[1]);
	 bin_end    = coot::util::string_to_float(bits[2]);
	 bin_width  = coot::util::string_to_float(bits[3]);
	 n_bins     = coot::util::string_to_float(bits[4]);
	 // 5 is/should be a ":"
	 for (unsigned int ibin=0; ibin<n_bins; ibin++) {
	    int ibit = ibin + 6;
	    if (ibit < bits.size()) {
	       int v = coot::util::string_to_int(bits[ibit]);
	       counts.push_back(v);
	    }
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
		  std::cout << "calling get_bond_type() for :"
			    << name_1 << ": :" << name_2 << ":" << std::endl;
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
