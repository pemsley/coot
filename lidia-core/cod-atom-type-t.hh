
#ifndef COD_ATOM_TYPE_T_HH
#define COD_ATOM_TYPE_T_HH

#include <list>

#include "third-neighbour-info-t.hh"


namespace cod {

   int hybridization_to_int(RDKit::Atom::HybridizationType h);

   // return the number of rings and the ringinfo string.
   std::pair<int, std::string> make_ring_info_string(RDKit::Atom *atom_p);

   class atom_level_2_type {
      std::string str;
      std::string element;
   public:

      atom_level_2_type() {}
      atom_level_2_type(const std::string &s) { str = s;} // read
      atom_level_2_type(RDKit::Atom *p, const RDKit::ROMol &rdkm);

      class atom_level_2_component_type {
      public:
	 std::string element;
	 unsigned int number_of_rings;
	 std::string ring_info_string;
	 std::vector<int> neighb_hybridizations;
	 atom_level_2_component_type(RDKit::Atom *at, const RDKit::ROMol &rdkm);
	 atom_level_2_component_type() {}
      };
      std::string string() const {return str; }
      static bool level_2_component_sorter(const atom_level_2_component_type &la,
					   const atom_level_2_component_type &lb);
      friend std::ostream &operator<<(std::ostream &s,
				      const atom_level_2_component_type &c);
   };
   std::ostream &operator<<(std::ostream &s,
			    const atom_level_2_type::atom_level_2_component_type &c);
   
   // a container for strings of COD atom types at various levels.
   // The first one was at the 4th level - most sophisticated and
   // with 3rd neighbour info
   //
   class atom_type_t {
   public:

      atom_type_t() {}

      atom_type_t(const std::string &s1, const atom_level_2_type &l2,
		  const std::string &s3, const std::string &s4);

      atom_type_t(const std::string &l4) {
	 level_4 = l4;
	 level_3 = level_4_type_to_level_3_type(l4);
	 hash_value = -1;
      }
      atom_type_t(const std::string &s1, const std::string &s2) {
	 level_3 = s1;
	 level_4 = s2;
	 hash_value = -1;
      }
      std::string level_4;
      std::string level_3; // as 4 but without 3rd neighbour info
      atom_level_2_type level_2;
      int hash_value; // can be zero (?), so use -1 for fail.
      std::list<third_neighbour_info_t> tnil;

      // helper function
      static std::string level_4_type_to_level_3_type(const std::string &l4t);

      bool operator==(const atom_type_t &t_in) const {
	 return (t_in.level_4 == level_4);
      }

      bool l3_match(const atom_type_t &t_in) const {
	 return (t_in.level_3 == level_3);
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
