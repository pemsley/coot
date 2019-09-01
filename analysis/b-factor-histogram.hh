
#include <vector>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   class b_factor_histogram {
      int n_bins;
      int n_atoms;
      float b_max;
      void init();
      int get_n_atoms() const;
      int get_n_bins() const;
      int b_to_bin(const float &b) const;
      std::vector<std::vector<float> > b_vector;
      float alpha_estimate;
      float beta_estimate;
      float shift_estimate;
      double ig(const double &x) const;
      double Gamma(const double &b) const;
      void optimize_estimates();
      std::vector<double> select_from_model() const;
   public:
      b_factor_histogram(mmdb::Manager *mol);
      b_factor_histogram(mmdb::Manager *mol, int atom_selection_handle);
      void model();
      std::vector<std::pair<double, double> > get_data() const;
      std::vector<std::pair<double, double> > get_model() const; // call model() first
   };
}

