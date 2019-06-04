
#ifndef NEW_LINKED_RESIDUE_T_HH
#define NEW_LINKED_RESIDUE_T_HH

namespace coot {

   class new_linked_residue_t {
   public:
      new_linked_residue_t(mmdb::Atom *at_1, mmdb::Atom *at_2,
			   mmdb::Residue *res_1_in, mmdb::Residue *res_2_in,
			   bool is_fixed_first_in, bool is_fixed_second_in,
			   const std::string &link_type,
			   bool order_switch_flag) : at_1(at_1), at_2(at_2), res_1(res_1_in), res_2(res_2_in),
						     link_type(link_type), order_switch_flag(order_switch_flag) {

	 is_fixed_first  = is_fixed_first_in;
	 is_fixed_second = is_fixed_second_in;

	 if (order_switch_flag) {
	    is_fixed_first = is_fixed_second_in;
	    is_fixed_second = is_fixed_first_in;
	    std::swap(res_1, res_2);
	    order_switch_flag = false;
	 }
      }
      mmdb::Atom *at_1, *at_2;
      mmdb::Residue *res_1, *res_2;
      bool is_fixed_first;
      bool is_fixed_second;
      std::string link_type;
      bool order_switch_flag;
      bool match_p(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 if (r1 == res_1 && r2 == res_2) {
	    return true;
	 } else {
	    return (r2 == res_1 && r1 == res_2);
	 }
      }
   };

   class new_linked_residue_list_t {
   public:
      new_linked_residue_list_t() {}
      std::vector<new_linked_residue_t> nlr_vec;
      void insert(mmdb::Atom *at_1, mmdb::Atom *at_2,
		  mmdb::Residue *res_1, mmdb::Residue *res_2,
		  bool is_fixed_first_in, bool is_fixed_second_in,
		  const std::string &link_type,
		  bool order_switch_flag) {
	 bool found = already_added_p(res_1, res_2);
	 if (! found) {
	    new_linked_residue_t nlr(at_1, at_2, res_1, res_2, is_fixed_first_in, is_fixed_second_in,
				     link_type, order_switch_flag);
	    nlr_vec.push_back(nlr);
	 }
      }
      bool already_added_p(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 bool found = false;
	 for (std::size_t i=0; i<nlr_vec.size(); i++) {
	    if (nlr_vec[i].match_p(r1, r2)) {
	       found = true;
	       break;
	    }
	 }
	 return found;
      }
      std::size_t size() const { return nlr_vec.size(); }
      bool empty() { return nlr_vec.empty(); }
      const new_linked_residue_t &operator[](std::size_t i) { return nlr_vec[i]; }
   };
}

#endif // NEW_LINED_RESIDUE_T_HH
