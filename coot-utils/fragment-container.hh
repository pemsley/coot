#ifndef FRAGMENT_INFO_HH
#define FRAGMENT_INFO_HH

#include <vector>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class fragment_container_t {
   public:

      class fragment_range_t {
      public:
         std::vector<mmdb::Residue *> residues; // ordered residues in range
         std::string chain_id;
	 residue_spec_t start_res;
	 residue_spec_t end_res;
	 fragment_range_t(const std::string &chain_id_in,
                          const residue_spec_t &r1, const residue_spec_t &r2) : chain_id(chain_id_in), start_res(r1), end_res(r2) {}
         // friend std::ostream &operator<<(std::ostream &s, const fragment_range_t &fc); // see below
      };
      // fix this another time!
      // std::ostream& operator<<(std::ostream &s, const fragment_range_t &fc);

      std::vector<fragment_range_t> ranges;
      fragment_container_t() {}
      void add(const fragment_range_t &r) {
         ranges.push_back(r);
      }
      void print_fragments() const;
   };

   fragment_container_t make_fragments(mmdb::Manager *mol);
}



#endif // FRAGMENT_INFO_HH
