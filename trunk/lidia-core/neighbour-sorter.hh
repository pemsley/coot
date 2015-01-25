

namespace coot {

   class chiral_neighbour_info_t {

   public:
      int idx_mmcif;
      int idx_atom_list;
      mmdb::Atom *at;
      chiral_neighbour_info_t(mmdb::Atom *at_in, int idx_mmcif_in, int idx_atom_list_in) {
	 at = at_in;
	 idx_mmcif = idx_mmcif_in;
	 idx_atom_list = idx_atom_list_in;
      }
      bool operator==(const chiral_neighbour_info_t &test_v) const {
	 if (at == test_v.at)
	    if (idx_mmcif == test_v.idx_mmcif)
	       if (idx_atom_list == test_v.idx_atom_list)
		  return true;
	 return false;
      } 
      static bool neighbour_sorter(const chiral_neighbour_info_t &v1,
				   const chiral_neighbour_info_t &v2);
   };
}
