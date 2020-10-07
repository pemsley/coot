
#ifndef PROTEIN_DONOR_ACCEPTORS_HH
#define PROTEIN_DONOR_ACCEPTORS_HH

#include <string>
#include <vector>
#include <map>
#include "hb-types.hh"

namespace coot {

   class quick_protein_donor_acceptors {
      void init();
   public:

      class key {
      public:
	 std::string res_name;
	 std::string atom_name;
	 key(const std::string &r, const std::string &a) : res_name(r), atom_name(a) {}
	 // this one is less than the key_in?
	 bool operator<(const key &key_in) const {
	    if (res_name < key_in.res_name) {
	       return true;
	    } else {
	       if (res_name > key_in.res_name) {
		  return false;
	       } else {
		  return atom_name < key_in.atom_name;
	       }
	    }
	 }
      };
      std::map<key, hb_t> hb_type_map;
      quick_protein_donor_acceptors() { init(); }
      hb_t get_type(const key &k) const;
      void test() const;
      std::pair<bool, bool> is_hydrogen_bond_by_types(const std::pair<key,key> &hbtp) const;
      std::pair<bool, bool> is_hydrogen_bond_by_types(const key &k1, const key &k2) const {
	 std::pair<key,key> kp(k1,k2);
	 return is_hydrogen_bond_by_types(kp);
      }
      std::vector<std::pair<bool, bool> > is_hydrogen_bond_by_types(std::vector<std::pair<key, key> > &hbtps) const;
   };

}

#endif // PROTEIN_DONOR_ACCEPTORS_HH
