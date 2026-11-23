
// Tests for coot::local_resolution(), the local resolution of a reconstruction from
// its two half-maps.
//
// The maps are made here, from random numbers, so this needs no data files. A map of
// noise has signal at every resolution, so two copies of it are perfectly correlated to
// the Nyquist limit, and two different ones are not correlated at all.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>

#include "coot-locres.hh"

const float grid_spacing = 2.0f; // Angstroms
const float nyquist = 2.0f * grid_spacing;

// the maps for the tests, filled by make_test_maps()
clipper::Xmap<float> xmap_signal;        // noise, used as the signal that the half-maps share
clipper::Xmap<float> xmap_other_signal;  // different noise
clipper::Xmap<float> xmap_half_1;        // the signal, with a little noise added where u is low and a lot where u is high
clipper::Xmap<float> xmap_half_2;        // the same, with different noise

int n_tests = 0;

// windows that are small enough for the test maps: 8 grid points apart, a mask radius
// of 4 and an edge of 4
coot::locres_params_t test_params() {
   coot::locres_params_t params;
   params.sampling_A = 16.0f;
   params.mask_radius_A = 8.0f;
   params.edge_width_A = 8.0f;
   return params;
}

// a P1 map of noise, with nu, nv and nw grid points along its axes
clipper::Xmap<float> make_noise_map(int nu, int nv, int nw, unsigned int seed) {
   clipper::Spacegroup spacegroup(clipper::Spacegroup::P1);
   clipper::Cell cell(clipper::Cell_descr(nu * grid_spacing, nv * grid_spacing, nw * grid_spacing, 90.0, 90.0, 90.0));
   clipper::Grid_sampling grid(nu, nv, nw);
   clipper::Xmap<float> xmap(spacegroup, cell, grid);
   std::mt19937 generator(seed);
   std::normal_distribution<float> distribution(0.0f, 1.0f);
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next())
      xmap[ix] = distribution(generator);
   return xmap;
}

void make_test_maps() {
   int n = 64;
   xmap_signal       = make_noise_map(n, n, n, 1);
   xmap_other_signal = make_noise_map(n, n, n, 2);
   clipper::Xmap<float> xmap_noise_1 = make_noise_map(n, n, n, 3);
   clipper::Xmap<float> xmap_noise_2 = make_noise_map(n, n, n, 4);
   xmap_half_1 = xmap_signal;
   xmap_half_2 = xmap_signal;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_signal.first(); !ix.last(); ix.next()) {
      float noise_level = (ix.coord().u() < n/2) ? 0.1f : 10.0f;
      xmap_half_1[ix] += noise_level * xmap_noise_1[ix];
      xmap_half_2[ix] += noise_level * xmap_noise_2[ix];
   }
}

// the values of this map that are not 0, sorted
std::vector<float> estimates(const clipper::Xmap<float> &xmap) {
   std::vector<float> v;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next())
      if (xmap[ix] != 0.0f)
         v.push_back(xmap[ix]);
   std::sort(v.begin(), v.end());
   return v;
}

// the values of this map in a cube of 9 grid points centred here, sorted
std::vector<float> values_around(const clipper::Xmap<float> &xmap, int u, int v, int w) {
   std::vector<float> values;
   for (int i=-4; i<=4; i++)
      for (int j=-4; j<=4; j++)
         for (int k=-4; k<=4; k++)
            values.push_back(xmap.get_data(clipper::Coord_grid(u + i, v + j, w + k)));
   std::sort(values.begin(), values.end());
   return values;
}

bool same_values(const clipper::Xmap<float> &xmap_a, const clipper::Xmap<float> &xmap_b) {
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_a.first(); !ix.last(); ix.next())
      if (xmap_a[ix] != xmap_b[ix])
         return false;
   return true;
}

// identical half-maps are correlated all the way to the Nyquist limit
int test_identical_half_maps() {
   int status = 1;
   clipper::Xmap<float> xmap = coot::local_resolution(xmap_signal, xmap_signal, nullptr, test_params());
   std::vector<float> v = estimates(xmap);
   std::size_t n_grid = 64 * 64 * 64;
   // the windows do not reach the corners of the map, so some of it has no estimate
   if (v.size() < n_grid / 2 || v.size() == n_grid) {
      std::cout << "   " << v.size() << " of " << n_grid << " grid points have an estimate" << std::endl;
      status = 0;
   }
   if (! v.empty()) {
      if (fabsf(v.front() - nyquist) > 1.0e-4f || fabsf(v.back() - nyquist) > 1.0e-4f) {
         std::cout << "   the estimates range from " << v.front() << " to " << v.back() << std::endl;
         status = 0;
      }
   }
   return status;
}

// the same, in a map that is not a cube
int test_identical_half_maps_in_a_box_that_is_not_a_cube() {
   int status = 1;
   clipper::Xmap<float> xmap_box = make_noise_map(64, 48, 40, 5);
   clipper::Xmap<float> xmap = coot::local_resolution(xmap_box, xmap_box, nullptr, test_params());
   std::vector<float> v = estimates(xmap);
   if (v.size() < 64 * 48 * 40 / 2) status = 0;
   if (! v.empty()) {
      if (fabsf(v.front() - nyquist) > 1.0e-4f || fabsf(v.back() - nyquist) > 1.0e-4f) {
         std::cout << "   the estimates range from " << v.front() << " to " << v.back() << std::endl;
         status = 0;
      }
   }
   return status;
}

// half-maps that have nothing in common have no resolution to speak of
int test_unrelated_half_maps() {
   int status = 0;
   clipper::Xmap<float> xmap = coot::local_resolution(xmap_signal, xmap_other_signal, nullptr, test_params());
   std::vector<float> v = estimates(xmap);
   if (! v.empty()) {
      float median = v[v.size()/2];
      // nothing is worse than the low resolution limit, give or take the rounding of the means
      if (median > 5.0f * nyquist && v.back() < 50.001f)
         status = 1;
      else
         std::cout << "   median " << median << " worst " << v.back() << std::endl;
   }
   return status;
}

// the estimate is good where the noise is low and bad where the noise is high
int test_resolution_follows_the_noise() {
   int status = 0;
   clipper::Xmap<float> xmap = coot::local_resolution(xmap_half_1, xmap_half_2, nullptr, test_params());
   // no window that reaches these places straddles the change in noise at u = 32
   std::vector<float> low_noise  = values_around(xmap, 16, 32, 32);
   std::vector<float> high_noise = values_around(xmap, 48, 32, 32);
   float low_noise_worst    = low_noise.back();
   float high_noise_median  = high_noise[high_noise.size()/2];
   if (low_noise.front() > 0.0f && low_noise_worst < 1.2f * nyquist && high_noise_median > 5.0f * nyquist)
      status = 1;
   else
      std::cout << "   low noise: " << low_noise.front() << " to " << low_noise_worst
                << ", high noise: median " << high_noise_median << std::endl;
   return status;
}

// a mask sets the outside to 0 and leaves the inside as it is without a mask
int test_mask() {
   int status = 1;
   clipper::Xmap<float> xmap_mask = xmap_signal;
   clipper::Xmap_base::Map_reference_index ix;
   unsigned int n_inside = 0;
   for (ix = xmap_mask.first(); !ix.last(); ix.next()) {
      clipper::Coord_grid c = ix.coord();
      int du = c.u() - 32;
      int dv = c.v() - 28;
      int dw = c.w() - 30;
      bool inside = (du * du + dv * dv + dw * dw < 12 * 12);
      xmap_mask[ix] = inside ? 1.0f : 0.0f;
      if (inside) n_inside++;
   }
   clipper::Xmap<float> xmap_without = coot::local_resolution(xmap_half_1, xmap_half_2, nullptr, test_params());
   clipper::Xmap<float> xmap_with    = coot::local_resolution(xmap_half_1, xmap_half_2, &xmap_mask, test_params());
   unsigned int n_wrong = 0;
   for (ix = xmap_mask.first(); !ix.last(); ix.next()) {
      if (xmap_mask[ix] > 0.5f) {
         if (xmap_with[ix] == 0.0f || fabsf(xmap_with[ix] - xmap_without[ix]) > 1.0e-5f)
            n_wrong++;
      } else {
         if (xmap_with[ix] != 0.0f)
            n_wrong++;
      }
   }
   if (n_wrong > 0) {
      std::cout << "   " << n_wrong << " grid points are wrong" << std::endl;
      status = 0;
   }
   if (n_inside == 0) status = 0;
   return status;
}

// the result does not depend on how the work is shared out
int test_threads() {
   int status = 1;
   coot::locres_params_t params = test_params();
   params.n_threads = 1;
   clipper::Xmap<float> xmap_1_thread = coot::local_resolution(xmap_half_1, xmap_half_2, nullptr, params);
   params.n_threads = 3;
   clipper::Xmap<float> xmap_3_threads = coot::local_resolution(xmap_half_1, xmap_half_2, nullptr, params);
   std::size_t n_calls = 0;
   params.parallel_for = [&n_calls] (std::size_t n, const std::function<void(std::size_t)> &task) {
      for (std::size_t i=n; i>0; i--) { // backwards, to show that the order does not matter
         task(i-1);
         n_calls++;
      }
   };
   clipper::Xmap<float> xmap_parallel_for = coot::local_resolution(xmap_half_1, xmap_half_2, nullptr, params);
   if (! same_values(xmap_1_thread, xmap_3_threads)) status = 0;
   if (! same_values(xmap_1_thread, xmap_parallel_for)) status = 0;
   if (n_calls == 0) status = 0;
   return status;
}

int test_maps_that_cannot_be_used() {
   int status = 1;
   // half-maps on different grids
   clipper::Xmap<float> xmap_small = make_noise_map(48, 48, 48, 6);
   try {
      coot::local_resolution(xmap_signal, xmap_small, nullptr, test_params());
      status = 0;
   }
   catch (const std::runtime_error &e) {
   }
   // a mask on a different grid
   try {
      coot::local_resolution(xmap_signal, xmap_signal, &xmap_small, test_params());
      status = 0;
   }
   catch (const std::runtime_error &e) {
   }
   // a window that is bigger than the map: the default window is 75 Angstroms across
   clipper::Xmap<float> xmap_tiny = make_noise_map(16, 16, 16, 7);
   try {
      coot::local_resolution(xmap_tiny, xmap_tiny);
      status = 0;
   }
   catch (const std::runtime_error &e) {
   }
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

int main() {

   int all_tests_status = 1; // fail

   make_test_maps();
   int status = 0;
   status += run_test(test_identical_half_maps,                             "identical half-maps");
   status += run_test(test_identical_half_maps_in_a_box_that_is_not_a_cube, "identical half-maps in a box that is not a cube");
   status += run_test(test_unrelated_half_maps,                             "unrelated half-maps");
   status += run_test(test_resolution_follows_the_noise,                    "resolution follows the noise");
   status += run_test(test_mask,                                            "mask");
   status += run_test(test_threads,                                         "threads");
   status += run_test(test_maps_that_cannot_be_used,                        "maps that cannot be used");
   if (status == n_tests) all_tests_status = 0;

   return all_tests_status;
}
