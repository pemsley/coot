
#ifndef IDEAL_RNA_HH
#define IDEAL_RNA_HH

#include <string>
#include "mmdb_manager.h"

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
      short int is_valid_base(char base) const;
      char antisense_base(char base, short int is_dna_flag) const;
      int mutate_res(CResidue *res, char base, short int is_dna_flag) const;
      void delete_o2_prime(CResidue *res) const;

   public:
      ideal_rna(const std::string &RNA_or_DNA, const std::string &form,
		short int single_stranged_flag,
		const std::string &sequence, CMMDBManager *standard_residues);

      CMMDBManager *make_molecule();

   };

} 

#endif // IDEAL_RNA_HH

