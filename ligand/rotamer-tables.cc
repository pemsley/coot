/* src/rotamer-tables.cc
 * 
 * Copyright 2008 The University of Oxford
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

#include <fstream>
#include <stdexcept>
#include <string.h>

#include "ccp4mg-utils/mgtree.h"

#include "primitive-chi-angles.hh"
#include "rotamer.hh"

#include "utils/coot-utils.hh"

#include "utils/logging.hh"
extern logging logger;

coot::a_rotamer_table::a_rotamer_table(const std::string &residue_name_in,
				       const std::string &file_name) {

   residue_name = residue_name_in;
   // std::cout << "DEBUG:: reading file for " << residue_name << "\n";
   n_chis = -1; // unset
   if (residue_name == "SER" || residue_name == "VAL" ||
       residue_name == "THR" || residue_name == "CYS" || residue_name == "PRO") {
      fill_chi_1(file_name);
      n_chis = 1;
   }
   if (residue_name == "ASN" || residue_name == "ASP" || residue_name == "PHE" ||
       residue_name == "TYR" || residue_name == "TRP" || residue_name == "HIS" ||
       residue_name == "ILE" || residue_name == "LEU") {
      fill_chi_1_2(file_name);
      n_chis = 2;
   }
   if (residue_name == "MET" || residue_name == "MSE" ||
       residue_name == "GLU" || residue_name == "GLN") {
      fill_chi_1_2_3(file_name);
      n_chis = 3;
   }
    if (residue_name == "LYS" || residue_name == "ARG") {
       fill_chi_1_2_3_4(file_name);
       n_chis = 4;
    }
}

// Throw a runtime_error on failure to read file.
void
coot::a_rotamer_table::fill_chi_1(const std::string& file_name) {

   std::ifstream f(file_name.c_str());
   char chars[1024];

   bool in_table = 0;
   float chi;
   float prob;
   float very_small = 0.000001;
   // pr_chi_1 = std::vector<float> (360, very_small);
   int chi_count = 0;
   long int n_x1_obj_count = -1;
   if (!f) {
      std::string mess = "Failed to open " + file_name;
      throw std::runtime_error(mess);
   } else { 
      while (!f.eof()) {
	 f >> chars;
	 if (! (f.eof())) {
	    // std::cout << chars << std::endl;
	    if (in_table) { 
	       chi = atof(chars);
	       f >> chars;
	       prob = atof(chars);
	       int chi_val = lrint(chi-0.4);
	       pr_chi_1[chi_val] = prob;
// 	       std::cout << ":" << chars << ": |" << chi << "| (" << chi_val << ") "
// 			 << float(chi_val + 0.5) << " " << pr_chi_1[chi_val]
// 			 << std::endl;
	    }
	    if (! in_table) {
	       // mark the start of the table
	       if (!strncmp("line.)", chars, 6)) { 
		  // std::cout << "DEBUG:: Table starts now..." << std::endl; 
		  in_table = 1;
	       }
	       // Let's get the block size of the sampling.  This
	       // gives us the vector index to chi conversion.
	       if (n_x1_obj_count == 0) {
		  std::string l(chars);
		  n_chi1_samples_per_360 = atoi(l.c_str());
// 		  std::cout << "DEBUG:: n_chi1_samples_per_360: " << n_chi1_samples_per_360
// 			    << std::endl;
		  n_x1_obj_count = -1; // done
		  pr_chi_1 = std::vector<float> (n_chi1_samples_per_360, very_small);
	       } 
	       if (n_x1_obj_count)
		  n_x1_obj_count--;
	       if (!strncmp("x1:", chars, 3)) {
		  n_x1_obj_count = 2;
	       }
	    }
	 } 
      }
   }
}

void
coot::a_rotamer_table::fill_chi_1_2(const std::string& file_name) {
   
   std::ifstream f(file_name.c_str());
   char chars[1024];
   
   bool in_table = 0;
   float chi_1;
   float chi_2;
   float prob;
   float very_small = 0.000001;
   int chi_count = 0;
   long int n_x1_obj_count = -1;
   long int n_x2_obj_count = -1;
   long int n_x3_obj_count = -1;
   long int n_x4_obj_count = -1;
   int samples_scale = 1;
   if (!f) {
      std::string mess = "Failed to open " + file_name;
      throw std::runtime_error(mess);
   } else { 
      while (!f.eof()) {
	 f >> chars;
	 if (! (f.eof())) {
	    // std::cout << chars << std::endl;
	    if (in_table) { 
	       chi_1 = atof(chars);
	       f >> chars;
	       chi_2 = atof(chars);
	       if (chi_1 < 0.0) chi_1 += 360.0;
	       if (chi_2 < 0.0) chi_2 += 360.0;
	       f >> chars;
	       prob = atof(chars);
	       int chi_val_1 = chi_angle_to_bin(chi_1, n_chi1_samples_per_360);
	       int chi_val_2 = chi_angle_to_bin(chi_2, n_chi2_samples_per_360);
//  	       std::cout << "   wrapped chis " << chi_1 << " " << chi_2 << std::endl;
// 	       std::cout << "   n_chi1 n_chi2 samples_per_360 " <<  n_chi1_samples_per_360 << " "
// 			 << n_chi2_samples_per_360 << std::endl;
// 	       std::cout << "    chi_val_1 and _2: " << chi_val_1 << " " << chi_val_2 << std::endl;
//  	       std::cout << "   setting bin " << chi_val_1 << " (of "
//  			 << pr_chi_1_2[chi_val_1].size() << ") "
//  			 << chi_val_2 << " (of " << pr_chi_1_2[chi_val_1].size()
//  			 << ") to " << prob << std::endl;
	       pr_chi_1_2[chi_val_1][chi_val_2] = prob;
//  	       std::cout << ":" << chars << ": |" << "| ("
// 			 << chi_val_1 << ") " << float(chi_val_1 + 0.5) << " "
// 			 << chi_val_2 << ") " << float(chi_val_2 + 0.5) << " "
// 			 << pr_chi_1_2[chi_val_1][chi_val_2]
//  			 << std::endl;
	    }
	    if (! in_table) {
	       // mark the start of the table
	       if (!strncmp("line.)", chars, 6)) { 
		  // std::cout << "Table starts now..." << std::endl;
		  in_table = 1;
	       }
	       // Let's get the block size of the sampling.  This
	       // gives us the vector index to chi conversion.
	       if (n_x1_obj_count == 0) {
		  std::string l(chars);
		  n_chi1_samples_per_360 = atoi(l.c_str());
//  		  std::cout << "n_chi1_samples_per_360: " << n_chi1_samples_per_360
//  			    << std::endl;
		  n_x1_obj_count = -1; // done
	       }
	       if (n_x2_obj_count == 1) {
		  std::string max_degrees_str(chars);
		  if (max_degrees_str == "180.0") {
		     samples_scale = 2;
		  }
	       }
	       if (n_x2_obj_count == 0) {
		  std::string l(chars);
		  n_chi2_samples_per_360 = atoi(l.c_str()) * samples_scale;
//  		  std::cout << "n_chi2_samples_per_360: " << n_chi2_samples_per_360
//  			    << std::endl;
		  n_x2_obj_count = -1; // done
		  std::vector<float> pr_chi_1 (n_chi2_samples_per_360, very_small);
		  pr_chi_1_2 = std::vector<std::vector<float> > (n_chi1_samples_per_360, pr_chi_1);
	       } 
	       if (n_x1_obj_count)
		  n_x1_obj_count--;
	       if (n_x2_obj_count)
		  n_x2_obj_count--;
	       if (!strncmp("x1:", chars, 3)) {
		  n_x1_obj_count = 2;
	       }
	       if (!strncmp("x2:", chars, 3)) {
		  n_x2_obj_count = 2;
	       }
	    }
	 } 
      }
   }
}

void
coot::a_rotamer_table::fill_chi_1_2_3(const std::string& file_name) {

   std::ifstream f(file_name.c_str());
   char chars[1024];
   
   bool in_table = 0;
   float chi_1;
   float chi_2;
   float chi_3;
   float prob;
   float very_small = 0.000001;
   int chi_count = 0;
   long int n_x1_obj_count = -1;
   long int n_x2_obj_count = -1;
   long int n_x3_obj_count = -1;
   int samples_scale = 1;
   if (!f) {
      std::string mess = "Failed to open " + file_name;
      throw std::runtime_error(mess);
   } else { 
      while (!f.eof()) {
	 f >> chars;
	 if (! (f.eof())) {
	    if (in_table) { 
	       chi_1 = atof(chars);
	       f >> chars;
	       chi_2 = atof(chars);
	       f >> chars;
	       chi_3 = atof(chars);
	       f >> chars;
	       prob = atof(chars);
	       
	       int chi_val_1 = chi_angle_to_bin(chi_1, n_chi1_samples_per_360);
	       int chi_val_2 = chi_angle_to_bin(chi_2, n_chi2_samples_per_360);
	       int chi_val_3 = chi_angle_to_bin(chi_3, n_chi3_samples_per_360);
	       
	       // 	       std::cout << "   wrapped chis " << chi_1 << " " << chi_2 << std::endl;
	       // 	       std::cout << "   setting bin " << chi_val_1 << " (of "
	       // 			 << pr_chi_1_2[chi_val_1].size() << ") "
	       // 			 << chi_val_2 << " (of " << pr_chi_1_2[chi_val_1].size()
	       // 			 << ") to " << prob << std::endl;
	       pr_chi_1_2_3[chi_val_1][chi_val_2][chi_val_3] = prob;
	       //  	       std::cout << ":" << chars << ": |" << "| ("
	       // 			 << chi_val_1 << ") " << float(chi_val_1 + 0.5) << " "
	       // 			 << chi_val_2 << ") " << float(chi_val_2 + 0.5) << " "
	       // 			 << pr_chi_1_2[chi_val_1][chi_val_2]
	       //  			 << std::endl;
	    }
	    if (! in_table) {
	       // mark the start of the table
	       if (!strncmp("line.)", chars, 6)) { 
		  // std::cout << "DEBUG:: Table starts now..." << std::endl;
		  in_table = 1;
	       }
	       // Let's get the block size of the sampling.  This
	       // gives us the vector index to chi conversion.
	       if (n_x1_obj_count == 0) {
		  std::string l(chars);
		  n_chi1_samples_per_360 = atoi(l.c_str());
//  		  std::cout << "n_chi1_samples_per_360: " << n_chi1_samples_per_360
//  			    << std::endl;
		  n_x1_obj_count = -1; // done
	       } 
	       if (n_x2_obj_count == 0) {
		  std::string l(chars);
		  n_chi2_samples_per_360 = atoi(l.c_str());
//  		  std::cout << "n_chi2_samples_per_360: " << n_chi2_samples_per_360
//  			    << std::endl;
		  n_x2_obj_count = -1; // done
	       } 
	       if (n_x3_obj_count == 1) {
		  std::string max_degrees_str(chars);
		  if (max_degrees_str == "180.0") {
		     samples_scale = 2;
		  }
	       }
	       if (n_x3_obj_count == 0) {
		  std::string l(chars);
		  n_chi3_samples_per_360 = atoi(l.c_str()) * samples_scale;
// 		  std::cout << "n_chi3_samples_per_360: " << n_chi3_samples_per_360
// 			    << std::endl;
		  n_x3_obj_count = -1; // done
		  std::vector<float> pr_chi_1 (n_chi3_samples_per_360, very_small);
		  std::vector<std::vector<float> > pr_chi_1_2(n_chi2_samples_per_360, pr_chi_1);
		  pr_chi_1_2_3 = std::vector<std::vector<std::vector<float> > > (n_chi1_samples_per_360, pr_chi_1_2);
	       } 
	       if (n_x1_obj_count)
		  n_x1_obj_count--;
	       if (n_x2_obj_count)
		  n_x2_obj_count--;
	       if (n_x3_obj_count)
		  n_x3_obj_count--;
	       if (!strncmp("x1:", chars, 3)) {
		  n_x1_obj_count = 2;
	       }
	       if (!strncmp("x2:", chars, 3)) {
		  n_x2_obj_count = 2;
	       }
	       if (!strncmp("x3:", chars, 3)) {
		  n_x3_obj_count = 2;
	       }
	    }
	 }
      }
   }
}

void
coot::a_rotamer_table::fill_chi_1_2_3_4(const std::string& file_name) {

   std::ifstream f(file_name.c_str());
   char chars[1024];
   
   bool in_table = 0;
   float chi_1;
   float chi_2;
   float chi_3;
   float chi_4;
   float prob;
   float very_small = 0.000001;
   int chi_count = 0;
   long int n_x1_obj_count = -1;
   long int n_x2_obj_count = -1;
   long int n_x3_obj_count = -1;
   long int n_x4_obj_count = -1;
   if (!f) {
      std::string mess = "Failed to open " + file_name;
      throw std::runtime_error(mess);
   } else { 
      while (!f.eof()) {
	 f >> chars;
	 if (! (f.eof())) {
	    // std::cout << chars << std::endl;
	    if (in_table) { 
	       chi_1 = atof(chars);
	       f >> chars;
	       chi_2 = atof(chars);
	       f >> chars;
	       chi_3 = atof(chars);
	       f >> chars;
	       chi_4 = atof(chars);
	       f >> chars;
	       prob = atof(chars);
	       	       
	       int chi_val_1 = chi_angle_to_bin(chi_1, n_chi1_samples_per_360);
	       int chi_val_2 = chi_angle_to_bin(chi_2, n_chi2_samples_per_360);
	       int chi_val_3 = chi_angle_to_bin(chi_3, n_chi3_samples_per_360);
	       int chi_val_4 = chi_angle_to_bin(chi_4, n_chi4_samples_per_360);
	       pr_chi_1_2_3_4[chi_val_1][chi_val_2][chi_val_3][chi_val_4] = prob;

	    }
	    if (! in_table) {
	       // mark the start of the table
	       if (!strncmp("line.)", chars, 6)) { 
		  // std::cout << "DEBUG:: Table starts now..." << std::endl;
		  in_table = 1;
	       }
	       // Let's get the block size of the sampling.  This
	       // gives us the vector index to chi conversion.
	       if (n_x1_obj_count == 0) {
		  std::string l(chars);
		  n_chi1_samples_per_360 = atoi(l.c_str());
// 		  std::cout << "n_chi1_samples_per_360: " << n_chi1_samples_per_360
// 			    << std::endl;
		  n_x1_obj_count = -1; // done
	       } 
	       if (n_x2_obj_count == 0) {
		  std::string l(chars);
		  n_chi2_samples_per_360 = atoi(l.c_str());
// 		  std::cout << "n_chi2_samples_per_360: " << n_chi2_samples_per_360
// 			    << std::endl;
		  n_x2_obj_count = -1; // done
	       } 
	       if (n_x3_obj_count == 0) {
		  std::string l(chars);
		  n_chi3_samples_per_360 = atoi(l.c_str());
// 		  std::cout << "n_chi3_samples_per_360: " << n_chi3_samples_per_360
// 			    << std::endl;
		  n_x3_obj_count = -1; // done
	       } 
	       if (n_x4_obj_count == 0) {
		  std::string l(chars);
		  n_chi4_samples_per_360 = atoi(l.c_str());
// 		  std::cout << "n_chi4_samples_per_360: " << n_chi4_samples_per_360
// 			    << std::endl;
		  n_x4_obj_count = -1; // done
		  std::vector<float> pr_chi_1 (n_chi4_samples_per_360, very_small);
		  std::vector<std::vector<float> > pr_chi_1_2(n_chi3_samples_per_360, pr_chi_1);
		  std::vector<std::vector<std::vector<float> > > pr_chi_1_2_3(n_chi2_samples_per_360, pr_chi_1_2);
		  
		  pr_chi_1_2_3_4 = std::vector<std::vector<std::vector<std::vector<float> > > > (n_chi1_samples_per_360, pr_chi_1_2_3);
		  
	       } 
	       if (n_x1_obj_count)
		  n_x1_obj_count--;
	       if (n_x2_obj_count)
		  n_x2_obj_count--;
	       if (n_x3_obj_count)
		  n_x3_obj_count--;
	       if (n_x4_obj_count)
		  n_x4_obj_count--;
	       if (!strncmp("x1:", chars, 3)) {
		  n_x1_obj_count = 2;
	       }
	       if (!strncmp("x2:", chars, 3)) {
		  n_x2_obj_count = 2;
	       }
	       if (!strncmp("x3:", chars, 3)) {
		  n_x3_obj_count = 2;
	       }
	       if (!strncmp("x4:", chars, 3)) {
		  n_x4_obj_count = 2;
	       }
	    }
	 }
      }
   }
}

// ignore_lys_and_arg_flag is a default arg (false, don't ignore)
void
coot::rotamer_probability_tables::fill_tables(const std::string &dir, bool ignore_lys_and_arg_flag) {

   // 20250802-PE the logic at the end of this function means that when this function is executed
   // a secondd time, either is_well_formatted_ is true or tried_and_failed_ is true.
   // Which means that this function can only the usefully run once.
   // That will need to obe reviewed when late-loading the LYS and ARG tables in Moorhen.

   if (is_well_formatted_)
      return; // filled already.
   if (tried_and_failed_)
      return; // tried before and failed.

   std::vector<std::pair<std::string, std::string> > res;
   res.push_back(std::pair<std::string, std::string> ("SER", "rota500-ser.data"));
   res.push_back(std::pair<std::string, std::string> ("VAL", "rota500-val.data"));
   res.push_back(std::pair<std::string, std::string> ("THR", "rota500-thr.data"));
   res.push_back(std::pair<std::string, std::string> ("CYS", "rota500-cys.data"));
   res.push_back(std::pair<std::string, std::string> ("PRO", "rota500-pro.data"));

   res.push_back(std::pair<std::string, std::string> ("ASN", "rota500-asn.data"));
   res.push_back(std::pair<std::string, std::string> ("ASP", "rota500-asp.data"));
   res.push_back(std::pair<std::string, std::string> ("PHE", "rota500-phetyr.data"));
   res.push_back(std::pair<std::string, std::string> ("TYR", "rota500-phetyr.data"));
   res.push_back(std::pair<std::string, std::string> ("TRP", "rota500-trp.data"));
   res.push_back(std::pair<std::string, std::string> ("HIS", "rota500-his.data"));
   res.push_back(std::pair<std::string, std::string> ("ILE", "rota500-ile.data"));
   res.push_back(std::pair<std::string, std::string> ("LEU", "rota500-leu.data"));

   res.push_back(std::pair<std::string, std::string> ("MET", "rota500-met.data"));
   res.push_back(std::pair<std::string, std::string> ("MSE", "rota500-met.data"));
   res.push_back(std::pair<std::string, std::string> ("GLU", "rota500-glu.data"));
   res.push_back(std::pair<std::string, std::string> ("GLN", "rota500-gln.data"));

   if (! ignore_lys_and_arg_flag) {
      res.push_back(std::pair<std::string, std::string> ("ARG", "rota500-arg.data"));
      res.push_back(std::pair<std::string, std::string> ("LYS", "rota500-lys.data"));
   }

   // xxx/share/coot/rama-data
   std::string file_name_stub = dir;
   bool bad_read = false;

   // std::cout << "INFO:: Reading rotamer probability tables...";
   logger.log(log_t::INFO, "Reading rotamer probability tables...");
   std::cout.flush();
   for (unsigned int i=0; i<res.size(); i++) {
      std::string file_name = file_name_stub + "/" + res[i].second;
      try {
         const std::string &res_name = res[i].first;
         coot::a_rotamer_table t = coot::a_rotamer_table(res_name, file_name);
         bool already_exists = false;

         for (const auto &t : tables) {
            if (t.residue_name == res_name) {
               already_exists = true;
               break;
            }
         }
         if (already_exists) {
            // skip the new read
         } else {
            tables.push_back(t);
         }
      }
      catch (const std::runtime_error &mess) {
         std::cout << "Failed to read rotamer probability table for " << res[i].first
                   << "\n" << mess.what() << std::endl;
         bad_read = 1;
      }
   }

   if (false)
      if (bad_read == false)
         bad_read = test_yourself();

   if (bad_read == false)
      is_well_formatted_ = true;
   else
      tried_and_failed_ = true;
}

// can throw an exception
const coot::a_rotamer_table &
coot::rotamer_probability_tables::operator[](unsigned int idx) const {

   if (idx < tables.size()) {
      return tables[idx];
   } else {
      std::string s("out-of-range rotamer (table)");
      throw std::runtime_error(s);
   }
}


// float, a state and a name:
coot::rotamer_probability_info_t
coot::rotamer_probability_tables::probability_this_rotamer(unsigned int i_table,
							   const std::vector<std::pair<int,float> > &chi_angles) const {

   if (false)
      std::cout << "rotamer_probability_tables::probability_this_rotamer()" << std::endl;

   std::vector<int> bins = chi_angles_to_bins(i_table, chi_angles);
   if (bins.size() != chi_angles.size()) {
      throw std::runtime_error("ERROR:: bin size and chi_angles size do not match");
   }
   unsigned int n = tables[i_table].n_chis;
   if (bins.size() != n) {
      std::string mess = "ERROR:: not enough chi angles found. bin.size() (";
      mess += coot::util::int_to_string(bins.size());
      mess += ") and n chis do not match ";
      mess += " (should be ";
      mess += coot::util::int_to_string(n);
      mess += ") for ";
      mess += tables[i_table].residue_name;
      throw std::runtime_error(mess);
   }

   if (false) { // debugging
      if (n == 1) std::cout << "debug:: " << tables[i_table].residue_name << "  bin: " << bins[0] << std::endl;
      if (n == 2) std::cout << "debug:: " << tables[i_table].residue_name << "angles " << chi_angles[0].second << " " << chi_angles[1].second
			    << "  bins: " << bins[0] << " " << bins[1] << std::endl;
      if (n == 3) std::cout << "debug:: " << tables[i_table].residue_name << "angles " << chi_angles[0].second << " " << chi_angles[1].second << " "
			    << chi_angles[2].second
			    << "  bins: " << bins[0] << " " << bins[1] << " " << bins[2] << std::endl;
      if (n == 4) std::cout << "debug:: " << tables[i_table].residue_name << "angles " << chi_angles[0].second << " " << chi_angles[1].second << " "
			    << chi_angles[2].second << " " << chi_angles[3].second
			    << "  bins: " << bins[0] << " " << bins[1] << " " << bins[2] << " " << bins[3]
			    << std::endl;
   }

   float pr = 0.0;
   if (n == 1) pr = tables[i_table].pr_chi_1[bins[0]];
   if (n == 2) pr = tables[i_table].pr_chi_1_2[bins[0]][bins[1]];
   if (n == 3) pr = tables[i_table].pr_chi_1_2_3[bins[0]][bins[1]][bins[2]];
   if (n == 4) pr = tables[i_table].pr_chi_1_2_3_4[bins[0]][bins[1]][bins[2]][bins[3]];

   // std::cout << "raw pr: " << pr << std::endl; 

   if (n<1 || n>4) {
      std::string mess = "ERROR: bad nbins chis " + coot::util::int_to_string(n);
      mess += " for i_table ";
      mess += coot::util::int_to_string(i_table);
      throw std::runtime_error(mess);
   }

   return rotamer_probability_info_t(rotamer_probability_info_t::OK, pr*100.0, tables[i_table].residue_name);

}

std::vector<coot::rotamer_probability_info_t>
coot::rotamer_probability_tables::probability_this_rotamer(mmdb::Residue *residue) const {

   if (false) { // debugging
      std::cout << "probability_this_rotamer(): getting probability_this_rotamer for "
		<< coot::residue_spec_t(residue) << " " << residue->GetResName() << std::endl;
   }

   std::string resname (residue->GetResName());
   if (resname == "GLY" || resname == "ALA") {
      coot::rotamer_probability_info_t pr(coot::rotamer_probability_info_t::RESIDUE_IS_GLY_OR_ALA,
					  1, resname);
      std::vector<coot::rotamer_probability_info_t> v;
      v.push_back(pr);
      return v;
   }

   coot::primitive_chi_angles pchis(residue);
   std::vector<coot::alt_confed_chi_angles> chis = pchis.get_chi_angles();

   // Note that primitive_chi_angles has only 1 chi angle for PRO
   // 
   // PRO only have one chi angle, (traditional/dunbrack
   // says that it has 2.  Molprobity probability tables say 1).

   int i_table = -1;

   if (false) { // debug
      std::cout << "probability_this_rotamer() tables.size() " << tables.size() << std::endl;
      for (unsigned int itab=0; itab<tables.size(); itab++) {
	 std::cout << "    " << tables[itab].residue_name << std::endl;
      }
   }

   for (unsigned int itab=0; itab<tables.size(); itab++) {
      if (tables[itab].residue_name == resname) {
	 i_table = itab;
	 break;
      }
   }

   if (i_table == -1) {
      std::string mess = "Fail to find rotamer table for residue type " + resname;
      throw std::runtime_error(mess);
   }
   
   std::vector<coot::rotamer_probability_info_t> v;
   for (unsigned int i_chi=0; i_chi<chis.size(); i_chi++) {
      //       std::cout << "DEBUG:: passing internal probability_this_rotamer chi_angles of size "
      //                 << chis[i_chi].chi_angles.size() << std::endl;
      coot::rotamer_probability_info_t pr = probability_this_rotamer(i_table, chis[i_chi].chi_angles);
      v.push_back(pr);
   }
   return v;
}

// throw an exception on failure to get a bin for the residue_type.
std::vector<int>
coot::rotamer_probability_tables::chi_angles_to_bins(unsigned int table_index, 
						     std::vector<std::pair<int,float> > chi_angles) const {

   // Note that for ASP, GLU, PHE and TYR we have a nomenclature
   // issue.  IUPAC rules say that chi2 (chi3 for GLU) should be -90<=chi2<=90.
   // However the tables go from (basically) 0->180.
   //
   // So, where we have ASP, GLU, PHE or TYR and chi2 (chi3 for GLU) <
   // 0, we should add 180 degrees to the value from the table.
   //
   float max_degrees = 360.0;
   if (tables[table_index].residue_name == "ASP" ||
       tables[table_index].residue_name == "TYR" ||
       tables[table_index].residue_name == "PHE") {
      if (chi_angles.size() > 1) { 
	 if (chi_angles[1].second < 0.0)
	    chi_angles[1].second += 180.0;
	 if (chi_angles[1].second > 180.0)
	    chi_angles[1].second -= 180.0;
      }
   }
   if (tables[table_index].residue_name == "GLU") { 
      if (chi_angles.size() > 2) { 
	 if (chi_angles[2].second < 0.0)
	    chi_angles[2].second += 180.0;
	 if (chi_angles[2].second > 180.0)
	    chi_angles[2].second -= 180.0;
      }
   }

   std::vector<int> r;
   for (unsigned int i=0; i<chi_angles.size(); i++) {
      float chi = chi_angles[i].second;
      if (chi < 0.0)
	 chi += 360.0;
      int n_chi_samples = tables[table_index].n_chi1_samples_per_360;
      if (i==1) n_chi_samples = tables[table_index].n_chi2_samples_per_360;
      if (i==2) n_chi_samples = tables[table_index].n_chi3_samples_per_360;
      if (i==3) n_chi_samples = tables[table_index].n_chi4_samples_per_360;

      float fbin = float (n_chi_samples) * chi/max_degrees;
      if ((fbin >= float(n_chi_samples)) || (fbin < 0)) {
	 std::string mess = "ERROR:: Bin failure! fbin is ";
	 mess += coot::util::float_to_string(fbin);
	 mess += " for chi=";
	 mess += coot::util::float_to_string(chi);
	 mess += " and n_chi_samples: ";
	 mess += coot::util::int_to_string(n_chi_samples);
	 throw std::runtime_error(mess);
      }
      int bin = lrintf(fbin-0.5);
      r.push_back(bin);
   } 
   return r;
}


// Note: return 0 on 'OK' state
bool
coot::rotamer_probability_tables::test_yourself() {


   // Write only? :-)
   
   std::vector<std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > > test_results;

   std::vector<std::pair< std::vector<std::pair<int,float> >, float> > results;
   std::vector<std::pair<int,float> > chi_local;
   float very_small = 0.000001;

   // VAL
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 0.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.00733626784));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 180.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.572362));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 62.3));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.12698412698412698));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 359.93));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.0082699));
   chi_local.clear();

   std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > residue("VAL", results);
   test_results.push_back(residue);


   // HIS
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 0.0));
   chi_local.push_back(std::pair<int,float>(2, 180.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, very_small));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 60.0));
   chi_local.push_back(std::pair<int,float>(2, 82.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.360434));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 60.0));
   chi_local.push_back(std::pair<int,float>(2, 282.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.565652));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 317.0));
   chi_local.push_back(std::pair<int,float>(2, 282.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.33));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 359.9));
   chi_local.push_back(std::pair<int,float>(2, 359.9));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, very_small));
   chi_local.clear();

   residue = std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > ("HIS", results);
   test_results.push_back(residue);
   
   // ASN
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 0.0));
   chi_local.push_back(std::pair<int,float>(2, 180.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, very_small));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 60.0));
   chi_local.push_back(std::pair<int,float>(2, 82.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.117189));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 60.0));
   chi_local.push_back(std::pair<int,float>(2, 282.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.0826276));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 317.0));
   chi_local.push_back(std::pair<int,float>(2, 282.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.074387));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 359.9));
   chi_local.push_back(std::pair<int,float>(2, 359.9));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, very_small));
   chi_local.clear();

   residue = std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > ("ASN", results);
   test_results.push_back(residue);

   // ASP
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 184.69));  // B 93 ASP rnase
   chi_local.push_back(std::pair<int,float>(2, 160.04));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.5956072351));
   chi_local.clear();

   residue = std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > ("ASP", results);
   test_results.push_back(residue);
   
   // TYR
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 0.0));
   chi_local.push_back(std::pair<int,float>(2, 180.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, very_small));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 60.1));
   chi_local.push_back(std::pair<int,float>(2, 82.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.413169));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 60.0));
   chi_local.push_back(std::pair<int,float>(2, 282.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.29954499));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 317.0));
   chi_local.push_back(std::pair<int,float>(2, 282.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.1205763397));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 359.9));
   chi_local.push_back(std::pair<int,float>(2, 359.9));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, very_small));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, -65.0));
   chi_local.push_back(std::pair<int,float>(2, -87.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.964358));
   chi_local.clear();

   residue = std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > ("TYR", results);
   test_results.push_back(residue);

   // GLN
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 4.0));
   chi_local.push_back(std::pair<int,float>(2, 180.0));
   chi_local.push_back(std::pair<int,float>(3, 4.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.0054330457));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 359.99));
   chi_local.push_back(std::pair<int,float>(2, 359.99));
   chi_local.push_back(std::pair<int,float>(3, 183.00));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 9.5877277085e-4));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 51.99));
   chi_local.push_back(std::pair<int,float>(2, 170.99));
   chi_local.push_back(std::pair<int,float>(3, 313.00));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.2144455097));
   chi_local.clear();

   residue = std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > ("GLN", results);
   test_results.push_back(residue);

   // GLU
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 348.0));
   chi_local.push_back(std::pair<int,float>(2, 172.0));
   chi_local.push_back(std::pair<int,float>(3,  3.90));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.0043528064146));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 356.0));
   chi_local.push_back(std::pair<int,float>(2, 300.0));
   chi_local.push_back(std::pair<int,float>(3, 293.0));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.00389461626));
   chi_local.clear();

   residue = std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > ("GLU", results);
   test_results.push_back(residue);

   
   // LYS
   results.clear();
   chi_local.push_back(std::pair<int,float>(1, 359.99));
   chi_local.push_back(std::pair<int,float>(2, 359.99));
   chi_local.push_back(std::pair<int,float>(3, 359.99));
   chi_local.push_back(std::pair<int,float>(4, 359.99));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, very_small));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1, 351.0));
   chi_local.push_back(std::pair<int,float>(2, 291.99));
   chi_local.push_back(std::pair<int,float>(3, 203.99));
   chi_local.push_back(std::pair<int,float>(4, 299.99));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 9.96512207274539E-4));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1,  35.99));
   chi_local.push_back(std::pair<int,float>(2, 152.99));
   chi_local.push_back(std::pair<int,float>(3,  63.99));
   chi_local.push_back(std::pair<int,float>(4, 173.99));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.10313901));
   chi_local.clear();

   chi_local.push_back(std::pair<int,float>(1,  43.99));
   chi_local.push_back(std::pair<int,float>(2, 163.99));
   chi_local.push_back(std::pair<int,float>(3, 173.99));
   chi_local.push_back(std::pair<int,float>(4, 168.99));
   results.push_back(std::pair< std::vector<std::pair<int,float> >, float>(chi_local, 0.3266068759342));
   chi_local.clear();

   residue = std::pair<std::string, std::vector<std::pair< std::vector<std::pair<int,float> >, float> > > ("LYS", results);
   test_results.push_back(residue);
   

   
   // ------------------------------------------------------------------------------
   

   // Now compare vs actual results:
   //
   bool fail = 0;
   for(unsigned int itres=0; itres<test_results.size(); itres++) {
      std::string residue_name = test_results[itres].first;
      for (unsigned int itab=0; itab<tables.size(); itab++) {
	 if (tables[itab].residue_name == residue_name) {
	    for (unsigned int i_res_test=0; i_res_test<test_results[itres].second.size(); i_res_test++) {
	       chi_local = test_results[itres].second[i_res_test].first;
	       float correct_result = test_results[itres].second[i_res_test].second;
	       try { 
		  coot::rotamer_probability_info_t actual_result = probability_this_rotamer(itab, chi_local);
		  if (fabsf(actual_result.probability - correct_result*100) > 0.001) {
		     std::cout << "Test chi " << residue_name << " " ;
		     for (unsigned int i_chi = 0; i_chi<chi_local.size(); i_chi++)
			std::cout << "[" << chi_local[i_chi].first << " " << chi_local[i_chi].second << "] ";
		     std::cout << "======= fail match: " << actual_result.probability << " should be "
			       << correct_result*100 << std::endl;
		     fail = 1;
		  }
	       }
	       catch (const std::runtime_error &mess) {
		  std::cout << "Oooops: " << mess.what() << " on testing " << residue_name << " ";
		  for (unsigned int i_chi = 0; i_chi<chi_local.size(); i_chi++)
		     std::cout << "[" << chi_local[i_chi].first << " " << chi_local[i_chi].second << "] ";
		  std::cout << std::endl;
		  fail = 1;
	       }
	    } 
	 }
      }
   }
   if (fail == 0)
      std::cout << "DEBUG:: rotamer_probability_tables self test passed" << std::endl;
   else 
      std::cout << "DEBUG:: rotamer_probability_tables self test failed" << std::endl;
   
   return fail;
}

