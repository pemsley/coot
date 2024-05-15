
class lsq_results_t {
   public:
   std::vector<double> rotation_matrix;
   std::vector<double> translation;
   bool empty() { return rotation_matrix.empty(); }
};
