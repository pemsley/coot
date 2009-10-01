/* coot-utils/coot-fasta.hh
 * 
 * Copyright 2006, by The University of York
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

#ifndef COOT_FASTA_HH
#define COOT_FASTA_HH

#ifndef HAVE_STRING
#include<string>
#define HAVE_STRING
#endif


namespace coot {

   class fasta {
   public:
      std::string name;
      std::string sequence;
      short int is_fasta_aa(const std::string &a) const;
      fasta(const std::string &combined_string); // decomposition happens in constructor
      fasta(const std::string &name_in, const std::string &seq);
      std::string format() const {
	 std::string s = "> ";
	 s += name;
	 s += "\n";
	 s += sequence;
	 return s;
      }
   };

} 

#endif // COOT_FASTA_HH
