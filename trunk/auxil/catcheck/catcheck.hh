
#include <string>
#include <vector>
#include "mmdb_manager.h"

CMMDBManager *get_mol(const std::string &filename);
void water_coordination_check(CMMDBManager *mol, float max_dist);

class contact_info_t {
public:
   int this_index;
   std::vector<int> contact_indices;
   contact_info_t(const int &this_index_in) {
      this_index = this_index_in;
   } 
   contact_info_t() {};
   static bool contacts_less(const contact_info_t &a, const contact_info_t &b);
};

