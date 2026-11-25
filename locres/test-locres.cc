
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

// a mask sets the outside to 99.9 and leaves the inside as it is without a mask
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
         if (xmap_with[ix] != 99.9f)
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

// the mean of each grid point and its 26 neighbours: a map with less in it the higher the resolution
clipper::Xmap<float> smoothed(const clipper::Xmap<float> &xmap) {
   clipper::Xmap<float> xmap_smooth = xmap;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      clipper::Coord_grid c = ix.coord();
      float sum = 0.0f;
      for (int i=-1; i<=1; i++)
         for (int j=-1; j<=1; j++)
            for (int k=-1; k<=1; k++)
               sum += xmap.get_data(c + clipper::Coord_grid(i, j, k));
      xmap_smooth[ix] = sum / 27.0f;
   }
   return xmap_smooth;
}

// with match_relion, the resolution of a window is that of the last whole shell before its
// FSC falls below the threshold, where shell i is at (the length of the map) / i
int test_match_relion() {
   int status = 1;
   // a smooth signal under noise: the FSC is high at low resolution, falls with resolution and
   // crosses the threshold between the Nyquist limit and the low resolution limit
   clipper::Xmap<float> xmap_a = smoothed(xmap_signal);
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_a.first(); !ix.last(); ix.next())
      xmap_a[ix] *= 3.0f;
   clipper::Xmap<float> xmap_b = xmap_a;
   clipper::Xmap<float> xmap_noise_a = make_noise_map(64, 64, 64, 8);
   clipper::Xmap<float> xmap_noise_b = make_noise_map(64, 64, 64, 9);
   for (ix = xmap_a.first(); !ix.last(); ix.next()) {
      xmap_a[ix] += xmap_noise_a[ix];
      xmap_b[ix] += xmap_noise_b[ix];
   }
   coot::locres_params_t params = test_params();
   clipper::Xmap<float> xmap_interpolated = coot::local_resolution(xmap_a, xmap_b, nullptr, params);
   params.match_relion = true;
   clipper::Xmap<float> xmap_whole_shells = coot::local_resolution(xmap_a, xmap_b, nullptr, params);

   // only one window reaches each of these places, so the value there is the resolution of that window
   float map_length = 64.0f * grid_spacing;
   clipper::Coord_grid places[3] = { clipper::Coord_grid(2, 32, 32), clipper::Coord_grid(62, 32, 32),
                                     clipper::Coord_grid(32, 2, 32) };
   for (int i=0; i<3; i++) {
      float r_interpolated = xmap_interpolated.get_data(places[i]);
      float r_whole_shells = xmap_whole_shells.get_data(places[i]);
      float r_expected = map_length / floorf(map_length / r_interpolated + 1.0e-6f);
      bool between_the_limits = (r_interpolated > 1.2f * nyquist && r_interpolated < 30.0f);
      if (! between_the_limits || fabsf(r_whole_shells - r_expected) > 1.0e-3f) {
         std::cout << "   at " << places[i].format() << " interpolated " << r_interpolated << " whole shells "
                   << r_whole_shells << " expected " << r_expected << std::endl;
         status = 0;
      }
   }
   // a whole shell is never at a better resolution than the crossing
   unsigned int n_better = 0;
   unsigned int n_worse = 0;
   for (ix = xmap_interpolated.first(); !ix.last(); ix.next()) {
      if (xmap_whole_shells[ix] < xmap_interpolated[ix] - 1.0e-4f) n_better++;
      if (xmap_whole_shells[ix] > xmap_interpolated[ix] + 1.0e-4f) n_worse++;
   }
   if (n_better > 0 || n_worse == 0) {
      std::cout << "   whole shells better at " << n_better << " grid points, worse at " << n_worse << std::endl;
      status = 0;
   }
   return status;
}

// The phase randomization takes out of the FSC of a window what the window spreads from low
// resolution to higher resolution
int test_phase_randomization() {
   int status = 1;
   // Half-maps that have a few strong waves in common, all of them longer than 25 A, and
   // nothing else: the rest is noise. A window is not much bigger than one of the waves, so
   // it spreads them over many shells.
   const int n = 64;
   const int n_waves = 6;
   const int waves[n_waves][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4}, {3, 2, 1}, {1, 3, 2}, {2, 1, 3} };
   clipper::Xmap<float> xmap_a = make_noise_map(n, n, n, 10);
   clipper::Xmap<float> xmap_b = make_noise_map(n, n, n, 11);
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_a.first(); !ix.last(); ix.next()) {
      clipper::Coord_grid c = ix.coord();
      float sum = 0.0f;
      for (int i=0; i<n_waves; i++) {
         int phase = waves[i][0] * c.u() + waves[i][1] * c.v() + waves[i][2] * c.w();
         sum += cosf(6.2831853f * static_cast<float>(phase) / static_cast<float>(n) + static_cast<float>(i));
      }
      xmap_a[ix] += 10.0f * sum;
      xmap_b[ix] += 10.0f * sum;
   }
   coot::locres_params_t params = test_params();
   clipper::Xmap<float> xmap_plain = coot::local_resolution(xmap_a, xmap_b, nullptr, params);
   params.phase_randomization = true;
   clipper::Xmap<float> xmap_corrected = coot::local_resolution(xmap_a, xmap_b, nullptr, params);
   std::vector<float> v_plain = estimates(xmap_plain);
   std::vector<float> v_corrected = estimates(xmap_corrected);
   if (v_plain.empty() || v_corrected.empty()) return 0;
   float median_plain = v_plain[v_plain.size()/2];
   float median_corrected = v_corrected[v_corrected.size()/2];
   // Without the correction the waves look like signal to the Nyquist limit or near it. With it,
   // there is none beyond the first shell that is corrected, which is at 16 A in these windows.
   if (median_plain > 2.0f * nyquist || median_corrected < 3.5f * nyquist || median_corrected > 25.0f) {
      std::cout << "   median without the correction " << median_plain << ", with it " << median_corrected << std::endl;
      status = 0;
   }

   // the random numbers are the same every time, on any number of threads
   clipper::Xmap<float> xmap_again = coot::local_resolution(xmap_a, xmap_b, nullptr, params);
   params.n_threads = 1;
   clipper::Xmap<float> xmap_1_thread = coot::local_resolution(xmap_a, xmap_b, nullptr, params);
   if (! same_values(xmap_corrected, xmap_again)) status = 0;
   if (! same_values(xmap_corrected, xmap_1_thread)) status = 0;

   // identical half-maps are still correlated to the Nyquist limit
   params.n_threads = 0;
   clipper::Xmap<float> xmap_identical = coot::local_resolution(xmap_signal, xmap_signal, nullptr, params);
   std::vector<float> v_identical = estimates(xmap_identical);
   if (v_identical.empty() || fabsf(v_identical.front() - nyquist) > 1.0e-4f || fabsf(v_identical.back() - nyquist) > 1.0e-4f) {
      std::cout << "   identical half-maps: " << v_identical.size() << " estimates" << std::endl;
      status = 0;
   }
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

// the function that a program finds by name when it has loaded this library
int test_c_linkage_function() {
   int status = 1;
   coot::locres_params_t params = test_params();
   clipper::Xmap<float> xmap_expected = coot::local_resolution(xmap_half_1, xmap_half_2, nullptr, params);
   clipper::Xmap<float> xmap_result;
   std::string message;
   if (coot_local_resolution(&xmap_half_1, &xmap_half_2, nullptr, &params, &xmap_result, &message) != 1) status = 0;
   if (! message.empty()) status = 0;
   if (xmap_result.is_null()) {
      status = 0;
   } else {
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = xmap_expected.first(); !ix.last(); ix.next())
         if (xmap_result[ix] != xmap_expected[ix])
            status = 0;
   }
   // half-maps on different grids: it fails, says why, and leaves the result alone
   clipper::Xmap<float> xmap_small = make_noise_map(48, 48, 48, 6);
   clipper::Xmap<float> xmap_untouched;
   if (coot_local_resolution(&xmap_signal, &xmap_small, nullptr, &params, &xmap_untouched, &message) != 0) status = 0;
   if (message.empty()) status = 0;
   if (! xmap_untouched.is_null()) status = 0;
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
   status += run_test(test_match_relion,                                    "match RELION");
   status += run_test(test_phase_randomization,                             "phase randomization");
   status += run_test(test_maps_that_cannot_be_used,                        "maps that cannot be used");
   status += run_test(test_c_linkage_function,                              "the function with C linkage");
   if (status == n_tests) all_tests_status = 0;

   return all_tests_status;
}
