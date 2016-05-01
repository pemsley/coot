
#ifndef COD_ATOM_TYPE_T_HH
#define COD_ATOM_TYPE_T_HH

#include <list>

#include "third-neighbour-info-t.hh"


namespace cod {
   
   // a container for strings of COD atom types at various levels.
   // The first one was at the 4th level - most sophisticated and
   // with 3rd neighbour info
   //
   class atom_type_t {
   public:

      atom_type_t() {}

      atom_type_t(const std::string &l4) {
	 level_4 = l4;
	 level_3 = level_4_type_to_level_3_type(l4);
      }
      atom_type_t(const std::string &s1, const std::string &s2) {
	 level_3 = s1;
	 level_4 = s2;
      }
      std::string level_4;
      std::string level_3; // as 4 but without 3rd neighbour info
      std::string level_2;
      std::string hash_value;
      std::list<third_neighbour_info_t> tnil;

      // helper function
      static std::string level_4_type_to_level_3_type(const std::string &l4t);

      bool operator==(const atom_type_t &t_in) const {
	 return (t_in.level_4 == level_4);
      }

      bool operator<(const atom_type_t &t_in) const {
	 return (level_4 < t_in.level_4); // is this the right way round?
      }

      void set_hash_value(unsigned int hash_in) {
	 hash_value = hash_in;
      }
   };
}

#endif // COD_ATOM_TYPE_T_HH
