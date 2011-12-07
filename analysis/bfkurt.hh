
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

#include <mmdb/mmdb_manager.h>

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
      my_stats_t stats(CResidue *residue_p) const;
      void set_questionable_flags(float z);
      bool is_mol_from_shelx_flag;

   public:
      b_factor_analysis(const CMMDBManager *mol, bool is_from_shelx_ins_flag_in);
      short int write_table(const std::string &filename,
			    const std::string &pdb_filename,
			    short int write_only_questionables_flag) const;

      std::vector<my_chain_of_stats_t> chain_details() const;

   };

} 


#endif // HAVE_BFKURT_HH
