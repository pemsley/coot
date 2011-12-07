/* ligand/ideal-rna.hh
 * 
 * Copyright 2006 by The University of York
 * Copyright 2009 by the University of Oxford
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
 * 02110-1301, USA.
 */

#ifndef IDEAL_RNA_HH
#define IDEAL_RNA_HH

#include <string>
#include <mmdb/mmdb_manager.h>

namespace coot {

   class ideal_rna {

   public:
      enum sense_direction_t {SENSE, ANTISENSE};
      enum form_t {A_FORM, B_FORM};

   private:
      std::string RNA_or_DNA_;
      std::string form_;
      std::string seq;
      short int single_stranged_flag;
      CResidue * get_standard_residue_instance(const std::string &residue_type,
					       CMMDBManager *standard_residues) const;
      clipper::RTop_orth n_turns(int nbase, int n_in_chain, coot::ideal_rna::form_t form_flag) const;
      CMMDBManager *standard_residues;
      bool is_valid_base(char base) const;
      char antisense_base(char base, bool is_dna_flag) const;
      int mutate_res(CResidue *res, char base, bool is_dna_flag) const;
      void delete_o2_prime(CResidue *res) const; // RNA -> DNA
      void add_o2_prime(CResidue *res) const;    // DNA -> RNAv

      // where "standard" means old-style:  "Ar", "Ad"...
      bool use_standard_refmac_names;

      // convert new names: "A", "DA" to old ones (the ones in the
      // standard residues file), "Ar", "Ad"
      // 
      std::string residue_name_from_type(const std::string &residue_type_in) const ;
      std::string residue_name_old_to_new(const std::string &residue_type_in, bool is_dna_flag) const;
      void fix_up_residue_and_atom_names(CResidue *res, bool is_dna_flag); // fiddle with res

   public:
      ideal_rna(const std::string &RNA_or_DNA, const std::string &form,
		short int single_stranged_flag,
		const std::string &sequence, CMMDBManager *standard_residues);
      void use_v3_names() { use_standard_refmac_names = 0; }

      CMMDBManager *make_molecule();

   };

} 

#endif // IDEAL_RNA_HH

