#ifndef GPHL_CHEM_COMP_INFO_HH
#define GPHL_CHEM_COMP_INFO_HH

#include <string>
#include <vector>

namespace coot {

   class gphl_chem_comp_info_t {
   public:
      gphl_chem_comp_info_t() {}
      std::vector<std::pair<std::string, std::string> > info;
      void add(const std::string &key, const std::string &data_string) {
         auto p = std::make_pair(key, data_string);
         info.push_back(p);
      }
      int get_index(const std::string &key) const {
         int idx = -1;
         for (unsigned int i=0; i<info.size(); i++) {
            if (info[i].first == key) {
               idx = i;
               break;
            }
         }
         return idx;
      }
      // caller ensures that the idx is in range (using the above function)
      std::pair<std::string, std::string> & operator[](int idx) {
         return info[idx];
      }
   };

}

#endif // GPHL_CHEM_COMP_INFO_HH
