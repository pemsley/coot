
#include "mini-mol.hh"

namespace coot {

   class db_strands {

   public:
      db_strands() {}
      std::vector<minimol::molecule> get_reference_strands(int n_strands, int strand_length);
   };
} 
