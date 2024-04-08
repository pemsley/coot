/*
 * analysis/bfkurt.hh
 *
 * Copyright 2009 by University of Oxford
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#ifndef HAVE_BFKURT_HH
#define HAVE_BFKURT_HH

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#include <mmdb2/mmdb_manager.h>

namespace coot_extras { 

   // a trivial helper class
   class my_stats_t {
   public:
      float mean;
      float std_dev;
      float skew;
      float kurtosis;
      int n; // atoms in the residue
      int resno;
      std::string inscode;
      std::string resname;
      std::string atom_name; // for the graph block click callback
      short int questionable_flag;
      my_stats_t() {
	 n = 0;
	 questionable_flag = 0;
	 atom_name = " CA ";
      }
   };

   class my_chain_of_stats_t {
   public:
      std::vector<my_stats_t> residue_properties; // one for each residue
      std::string chain_id;
   };

   class b_factor_analysis { 
      std::vector<std::pair<std::string, std::vector<my_stats_t> > > kurtoses;
      my_stats_t stats(mmdb::Residue *residue_p) const;
      void set_questionable_flags(float z);
      bool is_mol_from_shelx_flag;

   public:
      b_factor_analysis(const mmdb::Manager *mol, bool is_from_shelx_ins_flag_in);
      short int write_table(const std::string &filename,
			    const std::string &pdb_filename,
			    short int write_only_questionables_flag) const;

      std::vector<my_chain_of_stats_t> chain_details() const;

   };

} 


#endif // HAVE_BFKURT_HH
