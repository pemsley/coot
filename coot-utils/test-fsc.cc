
// Tests for coot::util::fsc(), the Fourier Shell Correlation between two maps,
// with and without a mask.
//
// The maps are made from the tutorial data, so set COOT_TEST_DATA_DIR to the
// top of the coot source tree to run this before coot has been installed.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "utils/coot-utils.hh"
#include "coot-map-utils.hh"

typedef std::vector<std::pair<clipper::Resolution, double> > fsc_curve_t;

// the maps for the tests, filled by read_test_maps()
clipper::Xmap<float> xmap_2fofc;
clipper::Xmap<float> xmap_fofc;
clipper::Xmap<float> xmap_other_grid; // not on the same grid as the others

int n_tests = 0;

// Resolve a coot data fixture across the two places it can live: the source tree at
// build/test time (pointed at by COOT_TEST_DATA_DIR), and $prefix/share/coot after
// installation (coot::package_data_dir()). Returns the first candidate that exists,
// else the bare filename.
std::string coot_test_data(const std::string &file) {
   std::vector<std::string> dirs;
   if (const char *e = getenv("COOT_TEST_DATA_DIR")) dirs.push_back(e);
   dirs.push_back(coot::package_data_dir());
   dirs.push_back(".");
   for (const auto &d : dirs) {
      std::string p = coot::util::append_dir_file(d, file);
      if (coot::file_exists(p))
         return p;
   }
   return file;
}

bool read_test_maps() {
   std::string mtz_file_name = coot_test_data("data/rnasa-1.8-all_refmac1.mtz");
   bool status_1 = coot::util::map_fill_from_mtz(&xmap_2fofc,      mtz_file_name, "FWT",    "PHWT",    "", 0);
   bool status_2 = coot::util::map_fill_from_mtz(&xmap_fofc,       mtz_file_name, "DELFWT", "PHDELWT", "", 0);
   bool status_3 = coot::util::map_fill_from_mtz(&xmap_other_grid, mtz_file_name, "FWT",    "PHWT",    "", 0, 2.5);
   if (! (status_1 && status_2 && status_3))
      std::cout << "Failed to make the test maps from " << mtz_file_name << std::endl;
   return (status_1 && status_2 && status_3);
}

// a map on the same grid as xmap_2fofc, with every grid point set to value
clipper::Xmap<float> make_flat_mask(float value) {
   clipper::Xmap<float> xmap_mask = xmap_2fofc;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_mask.first(); !ix.last(); ix.next())
      xmap_mask[ix] = value;
   return xmap_mask;
}

// The correlation in each of 100 shells, calculated here from the amplitudes and the
// phase differences, as a check on coot::util::fsc(). The shells are those of
// coot::util::fsc().
std::vector<double>
reference_fsc(const clipper::Xmap<float> &xmap_1, const clipper::Xmap<float> &xmap_2) {

   unsigned int n_bins = 100;
   float mg = coot::util::max_gridding(xmap_1);
   clipper::Resolution reso(2.0 * mg);
   clipper::HKL_info hkl_info(xmap_1.spacegroup(), xmap_1.cell(), reso, true);
   clipper::HKL_data<clipper::datatypes::F_phi<float> > fphis_1(hkl_info);
   clipper::HKL_data<clipper::datatypes::F_phi<float> > fphis_2(hkl_info);
   xmap_1.fft_to(fphis_1);
   xmap_2.fft_to(fphis_2);
   clipper::HKL_info::HKL_reference_index hri;
   double max_reso = 0.0;
   for (hri = fphis_1.first(); !hri.last(); hri.next()) {
      float ir = sqrt(hri.invresolsq());
      if (ir > max_reso)
         max_reso = ir;
   }
   std::vector<double> ff_sum(n_bins, 0.0);
   std::vector<double> f1_sqrd_sum(n_bins, 0.0);
   std::vector<double> f2_sqrd_sum(n_bins, 0.0);
   for (hri = fphis_1.first(); !hri.last(); hri.next()) {
      double f_1 = fphis_1[hri].f();
      double f_2 = fphis_2[hri].f();
      float ir = sqrtf(hri.invresolsq());
      unsigned int bin_no = static_cast<unsigned int>(static_cast<double>(n_bins) * ir/max_reso);
      if (bin_no >= n_bins) bin_no = n_bins - 1;
      // |F1||F2|cos(phi_1 - phi_2) is Re(F1 conj(F2))
      ff_sum[bin_no] += f_1 * f_2 * cos(fphis_1[hri].phi() - fphis_2[hri].phi());
      f1_sqrd_sum[bin_no] += f_1 * f_1;
      f2_sqrd_sum[bin_no] += f_2 * f_2;
   }
   std::vector<double> v(n_bins, 0.0);
   for (unsigned int i=0; i<n_bins; i++) {
      double d = f1_sqrd_sum[i] * f2_sqrd_sum[i];
      if (d > 0.0)
         v[i] = ff_sum[i]/sqrt(d);
   }
   return v;
}

// every shell that has anything in it should have a correlation of 1
int test_fsc_of_a_map_with_itself() {
   int status = 1;
   fsc_curve_t curve = coot::util::fsc(xmap_2fofc, xmap_2fofc);
   unsigned int n_ones = 0;
   for (unsigned int i=0; i<curve.size(); i++) {
      const double &c = curve[i].second;
      if (fabs(c - 1.0) < 1.0e-6) {
         n_ones++;
      } else {
         if (c != 0.0) {
            std::cout << "   shell " << i << " has correlation " << c << std::endl;
            status = 0;
         }
      }
   }
   if (curve.size() != 100) status = 0;
   if (n_ones < 50) status = 0;
   return status;
}

int test_fsc_against_reference() {
   int status = 1;
   fsc_curve_t curve = coot::util::fsc(xmap_2fofc, xmap_fofc);
   std::vector<double> reference = reference_fsc(xmap_2fofc, xmap_fofc);
   unsigned int n_informative = 0; // shells that are neither empty nor perfectly correlated
   if (curve.size() == reference.size()) {
      for (unsigned int i=0; i<curve.size(); i++) {
         const double &c = curve[i].second;
         if (fabs(c - reference[i]) > 1.0e-4) {
            std::cout << "   shell " << i << " has correlation " << c
                      << " but the reference is " << reference[i] << std::endl;
            status = 0;
         }
         if (fabs(c) > 0.01 && fabs(c) < 0.99)
            n_informative++;
      }
   } else {
      status = 0;
   }
   if (n_informative < 10) status = 0;
   return status;
}

// multiplying by 1 changes nothing
int test_masked_fsc_with_a_mask_of_ones() {
   int status = 1;
   clipper::Xmap<float> xmap_mask = make_flat_mask(1.0f);
   fsc_curve_t curve_masked   = coot::util::fsc(xmap_2fofc, xmap_fofc, xmap_mask);
   fsc_curve_t curve_unmasked = coot::util::fsc(xmap_2fofc, xmap_fofc);
   if (curve_masked.size() == curve_unmasked.size() && ! curve_masked.empty()) {
      for (unsigned int i=0; i<curve_masked.size(); i++) {
         double delta_c = curve_masked[i].second - curve_unmasked[i].second;
         double delta_r = curve_masked[i].first.limit() - curve_unmasked[i].first.limit();
         if (fabs(delta_c) > 1.0e-9 || fabs(delta_r) > 1.0e-9) {
            std::cout << "   shell " << i << " differs: " << curve_masked[i].second
                      << " vs " << curve_unmasked[i].second << std::endl;
            status = 0;
         }
      }
   } else {
      status = 0;
   }
   return status;
}

// nothing is left of either map, so every shell is without signal: 0, not NaN
int test_masked_fsc_with_a_mask_of_zeros() {
   int status = 1;
   clipper::Xmap<float> xmap_mask = make_flat_mask(0.0f);
   fsc_curve_t curve = coot::util::fsc(xmap_2fofc, xmap_fofc, xmap_mask);
   for (unsigned int i=0; i<curve.size(); i++) {
      if (curve[i].second != 0.0) {
         std::cout << "   shell " << i << " has correlation " << curve[i].second << std::endl;
         status = 0;
      }
   }
   if (curve.size() != 100) status = 0;
   return status;
}

int test_masked_fsc_with_a_mask_on_a_different_grid() {
   int status = 0;
   fsc_curve_t curve = coot::util::fsc(xmap_2fofc, xmap_fofc, xmap_other_grid);
   if (curve.empty())
      status = 1;
   return status;
}

int run_test(int (*test_func)(), const std::string &test_name) {
   n_tests++;
   int status = test_func();
   std::string status_string = "FAIL: ";
   if (status == 1)
      status_string = "PASS: ";
   std::cout << status_string << test_name << std::endl;
   return status;
}

int main(int argc, char **argv) {

   int all_tests_status = 1; // fail

   if (read_test_maps()) {
      int status = 0;
      status += run_test(test_fsc_of_a_map_with_itself,                   "fsc of a map with itself");
      status += run_test(test_fsc_against_reference,                      "fsc against a reference calculation");
      status += run_test(test_masked_fsc_with_a_mask_of_ones,             "masked fsc with a mask of ones");
      status += run_test(test_masked_fsc_with_a_mask_of_zeros,            "masked fsc with a mask of zeros");
      status += run_test(test_masked_fsc_with_a_mask_on_a_different_grid, "masked fsc with a mask on a different grid");
      if (status == n_tests) all_tests_status = 0;
   }
   return all_tests_status;
}
