
#include <vector>

class lsq_results_t {
   public:
   //! the 9-element rotation matrix
   //! Note to self: what is the ordering?
   std::vector<double> rotation_matrix;
   //! the translation vector
   std::vector<double> translation;
   bool empty() { return rotation_matrix.empty(); }
};
