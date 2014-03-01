
#ifndef DICT_MISMATCHES_HH
#define DICT_MISMATCHES_HH

namespace coot { 
   class bond_mismatch_t {
   public:
      std::string atom_id_1;
      std::string atom_id_2;
      double dist_1;
      double dist_2;
      double abs_diff;
      double diff;
      bond_mismatch_t(const std::string &a1,
		      const std::string &a2,
		      double d1, double d2) {
	 atom_id_1 = a1;
	 atom_id_2 = a2;
	 dist_1 = d1;
	 dist_2 = d2;
	 diff = dist_2 - dist_1;
	 abs_diff = fabs(diff);
      }
      bool operator<(const bond_mismatch_t &m) const {
	 return m.abs_diff < abs_diff;
      } 
   };

   class angle_mismatch_t {
   public:
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      double angle_1, angle_2;
      double abs_diff;
      double diff;
      angle_mismatch_t(const std::string &a1,
		       const std::string &a2,
		       const std::string &a3,
		       double a_1, double a_2) {
	 atom_id_1 = a1;
	 atom_id_2 = a2;
	 atom_id_3 = a3;
	 angle_1 = a_1;
	 angle_2 = a_2;
	 diff = angle_2 - angle_1;
	 abs_diff = fabs(diff);
      }
      bool operator<(const angle_mismatch_t &m) const {
	 return m.abs_diff < abs_diff;
      } 
   }; 
   
}

#endif // DICT_MISMATCHES_HH
