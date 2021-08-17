/* analysis/mogul.cc
 * 
 * Copyright 2012, 2013, 2014 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <complex>
#include <math.h>

#include "utils/coot-utils.hh"
#include "lidia-core/lig-build.hh"

#include "mogul-interface.hh"


std::ostream &
coot::operator<<(std::ostream &s, const mogul_item &it) {

   if (it.type == mogul_item::NONE)
      s << "NONE";
   if (it.type == mogul_item::BOND)
      s << "BOND";
   if (it.type == mogul_item::ANGLE)
      s << "ANGLE";
   if (it.type == mogul_item::TORSION)
      s << "TORSION";
   if (it.type == mogul_item::PLANE)
      s << "PLANE";

   s << " ";
   if (it.type == mogul_item::BOND)
      s << it.idx_1 << " " << it.idx_2;
   if (it.type == mogul_item::ANGLE)
      s << it.idx_1 << " " << it.idx_2 << " " << it.idx_3;
   
   s << " counts: " << it.counts << " value: " << it.value << " mean: " << it.mean
     << " median: " << it.median << " sd: " << it.std_dev;
   return s;
} 


void
coot::mogul::parse(const std::string &file_name) {

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
	 // std::cout << "Read " << lines.size() << " lines." << std::endl;

	 // 2 from the end because we do a lines+2 to try to read distributions
	 for (unsigned int iline=0; iline<(lines.size()); iline++) {
	    std::vector<std::string> bits = coot::util::split_string(lines[iline], ",");

	    // mmmm std::cout << "parsing line " << lines[iline] << std::endl;

	    if (bits[0] == "BOND") {
	       try {
		  mogul_item item = parse_item_line(bits, 2);
		  items.push_back(item);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING:: problem reading bond info: " << lines[iline] << " "
			    << rte.what() << std::endl;
	       } 
	    }
	    if (bits[0] == "ANGLE") {
	       try {
		  mogul_item item = parse_item_line(bits, 3);
		  items.push_back(item);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING:: problem reading angle info: " << lines[iline] << " "
			    << rte.what() << std::endl;
	       } 
	    }
	    if (bits[0] == "TORSION") {
	       try {
		  mogul_item item = parse_item_line(bits, 4);
		  // std::cout << "ACCEPTED:::::::::::::::::::: " << lines[iline] << std::endl;
		  items.push_back(item);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING:: problem reading torsion info: " << lines[iline] << " "
			    << rte.what() << std::endl;
		  // std::cout << "REJECTED:::::::::::::::::::: " << lines[iline] << std::endl;
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

   //         0        1            2         3    4     5         6                 7       8       9                  10                    11         12
   // Fragment Type,Atom Indices,Query Value,Hits,Mean,Median,Standard Deviation,|z-score|,d(min),Distribution Minimum,Distribution Maximum,Bin Width,Number of Bins

   mogul_item r;
   bool debug = false;

   float min_acceptable_sigma_bond = 0.013;
   float min_acceptable_sigma_angle = 1.0; // degrees

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

      if (debug) {

	 // if query value is blank then that's because mogul didn't
	 // write one out because it thought that the input as in 2d.
	 
	 std::cout << "got " << bits.size() << " bits " << std::endl;
	 for (unsigned int i=0; i<bits.size() && i<13; i++) { 
	    std::cout << i << ":  " << bits[i];
	    if (i ==  0) std::cout << "    -> Type";
	    if (i ==  1) std::cout << "    -> Atom-indices";
	    if (i ==  2) std::cout << "    -> Model-Value";
	    if (i ==  3) std::cout << "    -> Hits";
	    if (i ==  4) std::cout << "    -> Mean";
	    if (i ==  5) std::cout << "    -> Median";
	    if (i ==  6) std::cout << "    -> St-dev";
	    if (i ==  7) std::cout << "    -> |z|";
	    if (i ==  8) std::cout << "    -> d(min)";
	    if (i ==  9) std::cout << "    -> min";
	    if (i == 10) std::cout << "    -> max";
	    if (i == 11) std::cout << "    -> bin_width";
	    if (i == 12) std::cout << "    -> n_bins";
	    std::cout << "\n";
	 }
	 std::cout << "got " << indices.size() << " indices and n_idx is  " << n_idx << std::endl;
      }

      
      
      if (indices.size() > 1) {
	 if (debug) { 
	    std::cout << "model  value from :" << bits[2] << ":" << std::endl;
	    std::cout << "counts value from :" << bits[3] << ":" << std::endl;
	 } 
	 float model_value = coot::util::string_to_float(bits[2]);
	 int counts    = coot::util::string_to_int(bits[3]);

	 if (debug) { 
	    std::cout << "   model_value: " << model_value << std::endl;
	    std::cout << "   counts: " << counts << std::endl;
	 }

	 if (counts > 0) {

	    if (indices.size() == 2 || indices.size() == 3) {
	       if (debug)
		  std::cout << "indices.size() is " << indices.size() << std::endl;
	       mean    = coot::util::string_to_float(bits[4]);
	       median  = coot::util::string_to_float(bits[5]);
	       min     = coot::util::string_to_float(bits[9]);
	       max     = coot::util::string_to_float(bits[10]);
	       if (counts > 1) {
		  std_dev = coot::util::string_to_float(bits[6]);
		  z       = coot::util::string_to_float(bits[7]);
	       }
	       
	    } else {
// 	       dmin    = coot::util::string_to_float(bits[8]);
// 	       if (debug)
// 		  std::cout << "got dmin " << dmin << " from :" << dmin << ":" << std::endl;
	    } 
	 
	    if (debug)
	       std::cout << "indices.size() is " << indices.size() << " and n_idx is "
			 << n_idx << std::endl;
	    
	    if (n_idx == 2) {
	       if (indices.size() == 2) {
		  if (apply_minimum_sigma_cap) { 
		     if (std_dev < min_acceptable_sigma_bond) { 
			std_dev = min_acceptable_sigma_bond;
			z = fabsf(model_value - median)/std_dev;			
		     }
		  }
		  r = mogul_item(indices[0], indices[1], model_value,
				 counts, mean, median, std_dev, z);
	       }
	    }
	    if (n_idx == 3) {
	       if (indices.size() == 3) { 
		  if (apply_minimum_sigma_cap) { 
		     if (std_dev < min_acceptable_sigma_angle) { 
			std_dev = min_acceptable_sigma_angle;
			z = fabsf(model_value - median)/std_dev;
		     }
		  }
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
	       r.add_distribution(mogul_distribution(distribution_bits));
	    }

	    if (debug)
	       std::cout << " item " << indices[0] << " " << indices[1] << " mean: " 
			 << mean << " min: " << min << " max: " << max << " median: "
			 << median << " std: " << std_dev
			 << "  Z-score:" << z << std::endl;
	 }
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

      if (0)
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

      if (0) {
	 std::cout << "bin_start " << bin_start << std::endl;
	 std::cout << "bin_end   " << bin_end   << std::endl;
	 std::cout << "bin_width " << bin_width << std::endl;
	 std::cout << "n_bins    " << n_bins    << std::endl;
      }
      
      for (unsigned int ibin=0; ibin<n_bins; ibin++) {
	 unsigned int ibit = ibin + 4;
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
   ch.blue = 0.25 - (f-0.5)*(f-0.5);
   ch.red = pow(f, 0.2);
   ch.green = pow(1 - f, 0.2);

   return ch.hex();
} 


coot::dictionary_residue_restraints_t
coot::mogul::make_restraints(mmdb::Residue *residue_p,
			     const std::string &comp_id,
			     int imol,
			     const coot::protein_geometry &geom) {

   coot::dictionary_residue_restraints_t r(comp_id, -1);
   if (residue_p) { 
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::pair<bool, dictionary_residue_restraints_t> current_restraints = 
	 geom.get_monomer_restraints(comp_id, imol);

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
		  dict_bond_restraint_t rest(name_1, name_2, type, dist, esd, 0.0, 0.0, false);
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
	    if (idx_1 >= 0 && idx_1 < int(atom_names.size())) { 
	       if (idx_2 >= 0 && idx_2 < int(atom_names.size())) {
		  std::string name_1 = atom_names[idx_1];
		  std::string name_2 = atom_names[idx_2];
		  // bt is "" if not found.
		  std::string bt = bond_types_dict.get_bond_type(name_1, name_2);
		  float dist = items[i].median;
		  float esd  = items[i].std_dev;
		  dict_bond_restraint_t rest(name_1, name_2, bt, dist, esd, 0.0, 0.0, false);
		  r.bond_restraint.push_back(rest);
	       }
	    }
	 }


	 if (items[i].type == mogul_item::ANGLE) {
	    int idx_1 = items[i].idx_1 - 1;
	    int idx_2 = items[i].idx_2 - 1;
	    int idx_3 = items[i].idx_3 - 1;
	    if (idx_1 >= 0 && idx_1 < int(atom_names.size())) { 
	       if (idx_2 >= 0 && idx_2 < int(atom_names.size())) {
		  if (idx_3 >= 0 && idx_3 < int(atom_names.size())) {
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


#ifdef HAVE_GSL
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#endif 

void 
coot::mogul_item::ft_model_torsion_distribution() {

#ifdef HAVE_GSL

   int i;
   const int n = 36;
   double data[2*n];
     
   gsl_fft_complex_wavetable *wavetable;
   gsl_fft_complex_workspace *workspace;
     
   for (i = 0; i < n; i++) {
      REAL(data,i) = 0.0;
      IMAG(data,i) = 0.0;
   }
     
// data[0] = 1.0;
     
//    for (i = 1; i <= 10; i++) {
//       REAL(data,i) = REAL(data,n-i) = 1.0;
//    }


   std::cout << "c.f. n " << n << " distribution counts: " << distribution.counts.size() << std::endl;
   
   for (unsigned int i=0; i<distribution.counts.size(); i++) { 
      REAL(data, i) = distribution.counts[i];
      REAL(data, n-i-1) = distribution.counts[i];
      // std::cout << i << " " << distribution.counts[i] << std::endl;
   }
     
   for (i = 0; i < n; i++) {
	 printf ("%d: %e %e\n", i, REAL(data,i), IMAG(data,i));
   }
   printf ("\n");
     
   wavetable = gsl_fft_complex_wavetable_alloc (n);
   workspace = gsl_fft_complex_workspace_alloc (n);
     
   for (unsigned i = 0; i < wavetable->nf; i++) {
      printf ("# factor %d: %ld\n", i, wavetable->factor[i]);
   }
     
   gsl_fft_complex_forward (data, 1, n, wavetable, workspace);
     
   for (i = 0; i < n; i++) {
      printf ("%d: %e %e\n", i, REAL(data,i), IMAG(data,i));
   }

   double model[2*n];
   for (i = 0; i < n; i++) {
      model[i] = 0.0;
      model[i] = -82;
   }
   
   for (i = 0; i < 14; i++) {
      std::complex<double> c(REAL(data,i), IMAG(data,i));
      double r = std::abs(c);
      double phi = std::arg(c);
      std::cout << "r: " << r << "  phi " << phi << " from " << c << std::endl;
      for (unsigned int j=0; j<n; j++) { 
	 model[j] += 2/double(n)*r*cos(phi + 2*M_PI*double(i*j)/double(n));
      }
   }

   for (i = 0; i < n; i++) {
      std::cout << "model: " << i << " " << model[i] << std::endl;
   }
   
   gsl_fft_complex_backward (data, 1, n, wavetable, workspace);
   for (i = 0; i < n; i++) {
      printf ("reversed: %d %e %e\n", i, 1/double(n) * REAL(data,i), IMAG(data,i));
   }
   printf ("\n");
     
   gsl_fft_complex_wavetable_free (wavetable);
   gsl_fft_complex_workspace_free (workspace);

#endif    
}

void
coot::mogul_item::spline_model_torsion_distribution() {

   return; 

   for (unsigned int i=0; i<distribution.counts.size()-3; i+=3) { 
      lig_build::pos_t p1(i,distribution.counts[i]);
      lig_build::pos_t p2(i+1,distribution.counts[i+1]);
      lig_build::pos_t p3(i+2,distribution.counts[i+2]);
      lig_build::pos_t p4(i+3,distribution.counts[i+3]);
      for (double t=0; t<1; t+=0.1) {
	 lig_build::pos_t p12   = lig_build::pos_t::fraction_point(p1,   p2,   t);
	 lig_build::pos_t p23   = lig_build::pos_t::fraction_point(p2,   p3,   t);
	 lig_build::pos_t p34   = lig_build::pos_t::fraction_point(p3,   p4,   t);
	 lig_build::pos_t p123  = lig_build::pos_t::fraction_point(p12,  p23,  t);
	 lig_build::pos_t p234  = lig_build::pos_t::fraction_point(p23,  p34,  t);
	 lig_build::pos_t p1234 = lig_build::pos_t::fraction_point(p123, p234, t);
	 std::cout << p1234.x << " " << p1234.y << std::endl;
      }
   }
}

// Dont consider items with less than 2 hits
// 
// return a negative value on unable for some reason.
std::pair<float, coot::mogul_item>
coot::mogul::get_max_z_badness(mogul_item::type_t t) const {

   float r = -1;
   coot::mogul_item worst_item;

   for (unsigned int i_item=0; i_item<items.size(); i_item++) {
      const mogul_item &item = items[i_item];
      if (item.counts > 1) { 
	 if (item.type == mogul_item::BOND) { 
	    if (t == mogul_item::BOND || t == mogul_item::BOND_OR_ANGLE) {
	       float s = item.std_dev;

	       // capped sigma?
	       if (apply_minimum_sigma_cap) 
		  if (s < 0.015) // where does this number come from?
		     s = 0.015;
	       float n_z = fabsf(item.value-item.median)/s;
	       if (n_z > r) { 
		  r = n_z;
		  worst_item = item;
	       } 
	    }
	 }
	 if (item.type == mogul_item::ANGLE) { 
	    if (t == mogul_item::ANGLE || t == mogul_item::BOND_OR_ANGLE) {
	       float s = item.std_dev;
	       // cap sigma?
	       if (apply_minimum_sigma_cap) 
		  if (s < 1.0) // where does this number come from?  It seems sane.
		     s = 1.0;
	       float n_z =  fabsf(item.value-item.median)/s;
	       if (n_z > r) { 
		  r = n_z;
		  worst_item = item;
	       }
	    }
	 }
      }
   }
   return std::pair<float, mogul_item> (r, worst_item);
}
