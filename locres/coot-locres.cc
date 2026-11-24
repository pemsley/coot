
#include <cmath>
#include <vector>
#include <complex>
#include <string>
#include <random>
#include <thread>
#include <atomic>
#include <algorithm>
#include <stdexcept>

#include <clipper/core/fftmap.h>

#include "fsc-core.hh"
#include "coot-locres.hh"

namespace {

   // The sizes of a window, in grid points
   class window_geometry_t {
   public:
      int step;            // the distance between the centres of neighbouring windows
      int mask_radius;     // a window has its full weight within this radius
      int reach;           // a window has no weight at this radius or beyond: the mask radius plus the edge width
      int n_sub;           // the side of the cube of grid points that holds a window
      int n_fft;           // the side of the cube that is Fourier transformed: that cube, padded with zeros
      double grid_spacing; // in Angstroms
      double map_length;   // the longest side of the map, in Angstroms
      // for the phase randomization: shell i of the transform of the map is at a resolution of map_length / i
      int randomize_shell;       // phases are random from this shell of the transform of the map
      int first_corrected_shell; // the FSC of a window is corrected from this shell of the transform of its cube
   };

   // the index of a grid point in a map that has been copied to a vector by dense_map()
   inline std::size_t dense_index(const clipper::Grid_sampling &grid, int u, int v, int w) {
      return (static_cast<std::size_t>(u) * grid.nv() + v) * grid.nw() + w;
   }

   bool same_grid(const clipper::Xmap<float> &xmap_a, const clipper::Xmap<float> &xmap_b) {
      const clipper::Grid_sampling &gs_a = xmap_a.grid_sampling();
      const clipper::Grid_sampling &gs_b = xmap_b.grid_sampling();
      return (gs_a.nu() == gs_b.nu() && gs_a.nv() == gs_b.nv() && gs_a.nw() == gs_b.nw() &&
              xmap_a.spacegroup().hash() == xmap_b.spacegroup().hash() &&
              xmap_a.cell().equals(xmap_b.cell()));
   }

   // throw if these maps cannot be used
   void check_maps(const clipper::Xmap<float> &half_1, const clipper::Xmap<float> &half_2,
                   const clipper::Xmap<float> *xmap_mask_p) {

      std::string fn = "coot::local_resolution(): ";
      if (half_1.is_null() || half_2.is_null())
         throw std::runtime_error(fn + "a half-map is empty");
      if (! same_grid(half_1, half_2))
         throw std::runtime_error(fn + "the half-maps do not have the same cell and grid");
      if (xmap_mask_p) {
         if (xmap_mask_p->is_null())
            throw std::runtime_error(fn + "the mask is empty");
         if (! same_grid(half_1, *xmap_mask_p))
            throw std::runtime_error(fn + "the mask does not have the same cell and grid as the half-maps");
      }
      if (half_1.spacegroup().num_symops() != 1)
         throw std::runtime_error(fn + "the half-maps are not in P1");
      const clipper::Cell &cell = half_1.cell();
      double right_angle = 0.5 * clipper::Util::pi();
      double angle_tolerance = 1.0e-4; // radians
      if (fabs(cell.alpha() - right_angle) > angle_tolerance ||
          fabs(cell.beta()  - right_angle) > angle_tolerance ||
          fabs(cell.gamma() - right_angle) > angle_tolerance)
         throw std::runtime_error(fn + "the cell of the half-maps is not orthogonal");
      const clipper::Grid_sampling &grid = half_1.grid_sampling();
      double spacing_a = cell.a() / static_cast<double>(grid.nu());
      double spacing_b = cell.b() / static_cast<double>(grid.nv());
      double spacing_c = cell.c() / static_cast<double>(grid.nw());
      if (fabs(spacing_b - spacing_a) > 1.0e-3 * spacing_a || fabs(spacing_c - spacing_a) > 1.0e-3 * spacing_a)
         throw std::runtime_error(fn + "the grid spacing of the half-maps is not the same along each axis");
   }

   // The smallest size that is not less than n and that has no prime factor greater than 7.
   // Nothing here is Fourier transformed on a grid of any other size.
   int fft_size(int n) {
      while (true) {
         int m = n;
         for (int p : {2, 3, 5, 7})
            while (m % p == 0)
               m /= p;
         if (m == 1)
            return n;
         n++;
      }
   }

   // throw if the windows cannot be used with this map
   window_geometry_t make_window_geometry(const clipper::Xmap<float> &xmap, const coot::locres_params_t &params) {

      std::string fn = "coot::local_resolution(): ";
      if (! (params.sampling_A > 0.0f))
         throw std::runtime_error(fn + "the sampling must be greater than 0");
      if (params.mask_radius_A < 0.0f || params.edge_width_A < 0.0f)
         throw std::runtime_error(fn + "the mask radius and the edge width must not be negative");
      if (! (params.low_resolution_limit_A > 0.0f))
         throw std::runtime_error(fn + "the low resolution limit must be greater than 0");

      const clipper::Grid_sampling &grid = xmap.grid_sampling();
      window_geometry_t geom;
      geom.grid_spacing = xmap.cell().a() / static_cast<double>(grid.nu());
      geom.map_length = geom.grid_spacing * static_cast<double>(std::max(grid.nu(), std::max(grid.nv(), grid.nw())));
      // RELION rounds each of these to a whole number of grid points too
      geom.step        = static_cast<int>(std::lround(params.sampling_A    / geom.grid_spacing));
      geom.mask_radius = static_cast<int>(std::lround(params.mask_radius_A / geom.grid_spacing));
      int edge_width   = static_cast<int>(std::lround(params.edge_width_A  / geom.grid_spacing));
      geom.reach = geom.mask_radius + edge_width;
      geom.n_sub = 2 * geom.reach;
      geom.n_fft = fft_size(2 * geom.n_sub);
      // RELION randomizes phases beyond 25 A (--locres_randomize_at), from a whole shell, and
      // corrects the FSC from two shells after that. first_corrected_shell is at that resolution.
      geom.randomize_shell = std::max(1, static_cast<int>(geom.map_length / 25.0));
      double cube_length = geom.n_fft * geom.grid_spacing;
      double first_corrected = (geom.randomize_shell + 2) * cube_length / geom.map_length;
      geom.first_corrected_shell = static_cast<int>(std::ceil(first_corrected - 1.0e-9));
      if (geom.step < 1)
         throw std::runtime_error(fn + "the sampling is less than half of the grid spacing");
      if (geom.reach < 1)
         throw std::runtime_error(fn + "the window is smaller than the grid spacing");
      int n_min = std::min(grid.nu(), std::min(grid.nv(), grid.nw()));
      if (geom.n_sub > n_min || 2 * geom.mask_radius >= n_min)
         throw std::runtime_error(fn + "the window is too big for the map");
      return geom;
   }

   // Copy a P1 map to a vector that can be indexed by dense_index()
   std::vector<float> dense_map(const clipper::Xmap<float> &xmap) {

      const clipper::Grid_sampling &grid = xmap.grid_sampling();
      std::vector<float> dense(static_cast<std::size_t>(grid.nu()) * grid.nv() * grid.nw());
      std::size_t idx = 0;
      clipper::Xmap_base::Map_reference_coord ix(xmap);
      for (int u=0; u<grid.nu(); u++) {
         for (int v=0; v<grid.nv(); v++) {
            ix.set_coord(clipper::Coord_grid(u, v, 0));
            for (int w=0; w<grid.nw(); w++) {
               dense[idx++] = xmap[ix];
               ix.next_w();
            }
         }
      }
      return dense;
   }

   // A copy of this map (a vector made by dense_map()) in which the phase of each Fourier
   // coefficient from the shell randomize_shell on is a random number. The amplitudes are
   // kept. Like RELION, this gives every stored coefficient a phase of its own, also in the
   // planes of the half transform that hold a coefficient and its Friedel mate.
   //
   // Unlike RELION, this does not transform the map on its own grid if a side of that has a
   // prime factor greater than 7 (314 is 2 x 157, and slow). The map is then transformed in
   // the smallest box that has none, the rest of which is 0, and cut out of it afterwards.
   std::vector<float> phase_randomized(const std::vector<float> &map, const clipper::Grid_sampling &grid,
                                       const window_geometry_t &geom, std::mt19937 *generator_p) {

      clipper::Grid_sampling fft_grid(fft_size(grid.nu()), fft_size(grid.nv()), fft_size(grid.nw()));
      clipper::FFTmap_p1 fft(fft_grid);
      std::size_t idx = 0;
      for (int u=0; u<grid.nu(); u++)
         for (int v=0; v<grid.nv(); v++)
            for (int w=0; w<grid.nw(); w++)
               fft.real_data(clipper::Coord_grid(u, v, w)) = map[idx++];
      fft.fft_x_to_h(1.0);

      // The indices of a coefficient, in units of 1 / map_length. Those are its own indices if
      // the map is a cube and is transformed on its own grid, and then this is what RELION does.
      double scale_h = geom.map_length / (fft_grid.nu() * geom.grid_spacing);
      double scale_k = geom.map_length / (fft_grid.nv() * geom.grid_spacing);
      double scale_l = geom.map_length / (fft_grid.nw() * geom.grid_spacing);
      double limit_sq = static_cast<double>(geom.randomize_shell) * geom.randomize_shell;
      const double two_pi = 6.283185307179586;
      for (int h=0; h<fft_grid.nu(); h++) {
         double hs = scale_h * ((2 * h <= fft_grid.nu()) ? h : h - fft_grid.nu());
         for (int k=0; k<fft_grid.nv(); k++) {
            double ks = scale_k * ((2 * k <= fft_grid.nv()) ? k : k - fft_grid.nv());
            for (int l=0; l<=fft_grid.nw()/2; l++) {
               double ls = scale_l * l;
               if (hs * hs + ks * ks + ls * ls >= limit_sq) {
                  clipper::Coord_grid c(h, k, l);
                  float amplitude = std::abs(fft.cplx_data(c));
                  // the 32 bits of the generator, so that the phases are the same on every platform
                  float phase = static_cast<float>(two_pi * ((*generator_p)() / 4294967296.0));
                  fft.cplx_data(c) = std::polar(amplitude, phase);
               }
            }
         }
      }
      fft.fft_h_to_x(1.0);

      std::vector<float> randomized(map.size());
      idx = 0;
      for (int u=0; u<grid.nu(); u++)
         for (int v=0; v<grid.nv(); v++)
            for (int w=0; w<grid.nw(); w++)
               randomized[idx++] = fft.real_data(clipper::Coord_grid(u, v, w));
      return randomized;
   }

   // The weight of a window at each grid point of the cube that holds it. The centre of the
   // window is at index "reach" along each axis of the cube. This is RELION's raised cosine:
   // 1 within the mask radius, 0 at the reach and beyond, a half cycle of a cosine between.
   std::vector<float> window_weights(const window_geometry_t &geom) {

      int n = geom.n_sub;
      std::vector<float> weights(static_cast<std::size_t>(n) * n * n);
      double r_in  = geom.mask_radius;
      double r_out = geom.reach;
      std::size_t idx = 0;
      for (int i=0; i<n; i++) {
         for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
               double di = i - geom.reach;
               double dj = j - geom.reach;
               double dk = k - geom.reach;
               double d = sqrt(di * di + dj * dj + dk * dk);
               double weight = 0.0;
               if (d <= r_in)
                  weight = 1.0;
               else if (d < r_out)
                  weight = 0.5 - 0.5 * cos(clipper::Util::pi() * (r_out - d) / (r_out - r_in));
               weights[idx++] = static_cast<float>(weight);
            }
         }
      }
      return weights;
   }

   // Call func(map_idx, weight, i, j, k) for each grid point of the map at which the window
   // centred at "centre" has some weight. map_idx is the dense_index() of that grid point and
   // i, j, k are its indices in the cube that holds the window. A window does not wrap around
   // to the other side of the map: the part of it that is outside the map is left out, as it
   // is by RELION (a cryo-EM map is not periodic).
   template <typename T>
   void for_each_window_point(const clipper::Grid_sampling &grid, const clipper::Coord_grid &centre,
                              const std::vector<float> &weights, const window_geometry_t &geom, T func) {

      int n = geom.n_sub;
      for (int i=0; i<n; i++) {
         int u = centre.u() + i - geom.reach;
         if (u < 0 || u >= grid.nu()) continue;
         for (int j=0; j<n; j++) {
            int v = centre.v() + j - geom.reach;
            if (v < 0 || v >= grid.nv()) continue;
            std::size_t idx = (static_cast<std::size_t>(i) * n + j) * n;
            for (int k=0; k<n; k++) {
               int w = centre.w() + k - geom.reach;
               if (w < 0 || w >= grid.nw()) continue;
               float weight = weights[idx + k];
               if (weight > 0.0f)
                  func(dense_index(grid, u, v, w), weight, i, j, k);
            }
         }
      }
   }

   // Does the window centred here have any weight at a grid point that is inside the mask?
   bool window_reaches_mask(const std::vector<float> &mask, const clipper::Grid_sampling &grid,
                            const clipper::Coord_grid &centre, const std::vector<float> &weights,
                            const window_geometry_t &geom, float mask_threshold) {

      bool reaches = false;
      auto check = [&mask, &reaches, mask_threshold] (std::size_t map_idx, float, int, int, int) {
         if (mask[map_idx] > mask_threshold)
            reaches = true;
      };
      for_each_window_point(grid, centre, weights, geom, check);
      return reaches;
   }

   // The centres of the windows: the grid points whose indices are multiples of the step,
   // inside the ellipsoid that fits in the cell less the mask radius, and that reach the
   // mask if there is one (mask is empty if there is not). For a cubic cell these are
   // RELION's windows, in RELION's order.
   std::vector<clipper::Coord_grid>
   window_centres(const clipper::Grid_sampling &grid, const window_geometry_t &geom,
                  const std::vector<float> &weights, const std::vector<float> &mask, float mask_threshold) {

      std::vector<clipper::Coord_grid> centres;
      // the semi-axes of the ellipsoid
      long long a = grid.nu()/2 - geom.mask_radius;
      long long b = grid.nv()/2 - geom.mask_radius;
      long long c = grid.nw()/2 - geom.mask_radius;
      for (int u=0; u<grid.nu(); u+=geom.step) {
         for (int v=0; v<grid.nv(); v+=geom.step) {
            for (int w=0; w<grid.nw(); w+=geom.step) {
               long long x = u - grid.nu()/2;
               long long y = v - grid.nv()/2;
               long long z = w - grid.nw()/2;
               // (x/a)^2 + (y/b)^2 + (z/c)^2 < 1, in whole numbers
               if (x*x * b*b * c*c + y*y * a*a * c*c + z*z * a*a * b*b < a*a * b*b * c*c) {
                  clipper::Coord_grid centre(u, v, w);
                  if (mask.empty() || window_reaches_mask(mask, grid, centre, weights, geom, mask_threshold))
                     centres.push_back(centre);
               }
            }
         }
      }
      return centres;
   }

   // The shell of each Fourier coefficient of a transform of side n, in the order that
   // window_resolution() reads them. A shell is the distance from the origin in grid
   // points of the transform, rounded. Coefficients beyond shell n/2 get n/2 + 1.
   std::vector<unsigned short> make_shell_table(int n) {

      int n_l = n/2 + 1;
      std::vector<unsigned short> table(static_cast<std::size_t>(n) * n * n_l);
      std::size_t idx = 0;
      for (int h=0; h<n; h++) {
         double hs = (2 * h <= n) ? h : h - n;
         for (int k=0; k<n; k++) {
            double ks = (2 * k <= n) ? k : k - n;
            for (int l=0; l<n_l; l++) {
               double r = sqrt(hs * hs + ks * ks + static_cast<double>(l) * l);
               int shell = static_cast<int>(r + 0.5);
               if (shell > n_l) shell = n_l;
               table[idx++] = static_cast<unsigned short>(shell);
            }
         }
      }
      return table;
   }

   // Where the FSC first falls below the threshold, between the shells on either side of
   // that. fsc[i] is the FSC of shell i, which is at a resolution of box_length_A / i.
   double resolution_from_fsc(const std::vector<double> &fsc, double box_length_A, const coot::locres_params_t &params) {

      double threshold = params.fsc_threshold;
      double limit = params.low_resolution_limit_A;
      for (std::size_t i=1; i<fsc.size(); i++) {
         if (fsc[i] < threshold) {
            if (i == 1)
               return limit; // no signal
            double fraction = (fsc[i-1] - threshold) / (fsc[i-1] - fsc[i]);
            double resolution = box_length_A / (static_cast<double>(i-1) + fraction);
            return std::min(resolution, limit);
         }
      }
      // the FSC is above the threshold all the way to the last shell
      return box_length_A / static_cast<double>(fsc.size() - 1);
   }

   // The resolution that RELION gives to a window whose FSC crosses the threshold at this
   // resolution: that of the last whole shell before the crossing. RELION transforms the
   // whole map, so its shell i is at a resolution of map_length_A / i.
   double whole_shell_resolution(double resolution, double map_length_A, double limit) {

      double shell = std::floor(map_length_A / resolution + 1.0e-6); // a crossing on a shell is in that shell
      if (shell < 1.0)
         return limit;
      return std::min(map_length_A / shell, limit);
   }

   // The FSC of two maps in the window centred here, for each shell of the transform of the
   // cube that holds the window: multiply both maps by the window and Fourier transform them.
   std::vector<double> window_fsc(const std::vector<float> &map_1, const std::vector<float> &map_2,
                                  const clipper::Grid_sampling &grid, const clipper::Coord_grid &centre,
                                  const std::vector<float> &weights, const std::vector<unsigned short> &shell_table,
                                  const window_geometry_t &geom) {

      int n = geom.n_fft;
      clipper::Grid_sampling fft_grid(n, n, n);
      // clipper locks around the making of FFTW plans, so these are safe to use in threads
      clipper::FFTmap_p1 fft_1(fft_grid);
      clipper::FFTmap_p1 fft_2(fft_grid);

      auto fill = [&fft_1, &fft_2, &map_1, &map_2] (std::size_t map_idx, float weight, int i, int j, int k) {
         clipper::Coord_grid c(i, j, k);
         fft_1.real_data(c) = weight * map_1[map_idx];
         fft_2.real_data(c) = weight * map_2[map_idx];
      };
      for_each_window_point(grid, centre, weights, geom, fill);
      fft_1.fft_x_to_h(1.0);
      fft_2.fft_x_to_h(1.0);

      // Only half of the transform is stored. A coefficient with l between 0 and n/2 stands
      // for its Friedel mate too. Both mates are stored where l is 0 or n/2.
      int n_l = n/2 + 1;
      coot::fsc_shells_t shells(n_l);
      std::size_t idx = 0;
      for (int h=0; h<n; h++) {
         for (int k=0; k<n; k++) {
            for (int l=0; l<n_l; l++) {
               unsigned short shell = shell_table[idx++];
               if (shell < n_l) {
                  unsigned int n_times = (l == 0 || 2 * l == n) ? 1 : 2;
                  clipper::Coord_grid c(h, k, l);
                  shells.add(shell, fft_1.cplx_data(c), fft_2.cplx_data(c), n_times);
               }
            }
         }
      }
      std::vector<double> fsc(n_l);
      for (int i=0; i<n_l; i++)
         fsc[i] = shells.fsc(i);
      return fsc;
   }

   // The resolution of the window centred here: where the FSC of the half-maps in the window
   // falls below the threshold. If there are phase-randomized half-maps, the FSC is first
   // corrected for what the window alone gives, which is the FSC of those maps in the window.
   double window_resolution(const std::vector<float> &map_1, const std::vector<float> &map_2,
                            const std::vector<float> *random_1_p, const std::vector<float> *random_2_p,
                            const clipper::Grid_sampling &grid, const clipper::Coord_grid &centre,
                            const std::vector<float> &weights, const std::vector<unsigned short> &shell_table,
                            const window_geometry_t &geom, const coot::locres_params_t &params) {

      std::vector<double> fsc = window_fsc(map_1, map_2, grid, centre, weights, shell_table, geom);
      if (random_1_p && random_2_p) {
         std::vector<double> fsc_random = window_fsc(*random_1_p, *random_2_p, grid, centre, weights, shell_table, geom);
         for (std::size_t i=geom.first_corrected_shell; i<fsc.size(); i++) {
            double d = 1.0 - fsc_random[i];
            fsc[i] = (d > 0.0) ? (fsc[i] - fsc_random[i]) / d : 0.0;
         }
      }
      double resolution = resolution_from_fsc(fsc, geom.n_fft * geom.grid_spacing, params);
      if (params.match_relion)
         resolution = whole_shell_resolution(resolution, geom.map_length, params.low_resolution_limit_A);
      return resolution;
   }

   // Call task(i) for each i from 0 to n_tasks - 1, on n_threads threads (0 means one for each core)
   void run_tasks(std::size_t n_tasks, const std::function<void(std::size_t)> &task, int n_threads_in) {

      std::size_t n_threads = (n_threads_in > 0) ? n_threads_in : std::thread::hardware_concurrency();
      if (n_threads > n_tasks) n_threads = n_tasks;
      if (n_threads <= 1) {
         for (std::size_t i=0; i<n_tasks; i++)
            task(i);
      } else {
         std::atomic<std::size_t> next(0);
         auto worker = [&next, &task, n_tasks] () {
            for (std::size_t i=next++; i<n_tasks; i=next++)
               task(i);
         };
         std::vector<std::thread> threads;
         for (std::size_t i=0; i<n_threads; i++)
            threads.push_back(std::thread(worker));
         for (std::size_t i=0; i<n_threads; i++)
            threads[i].join();
      }
   }
}

clipper::Xmap<float>
coot::local_resolution(const clipper::Xmap<float> &half_1, const clipper::Xmap<float> &half_2,
                       const clipper::Xmap<float> *xmap_mask_p, const coot::locres_params_t &params) {

   check_maps(half_1, half_2, xmap_mask_p);
   const clipper::Grid_sampling &grid = half_1.grid_sampling();
   window_geometry_t geom = make_window_geometry(half_1, params);
   std::vector<float> weights = window_weights(geom);
   std::vector<float> mask;
   if (xmap_mask_p)
      mask = dense_map(*xmap_mask_p);
   std::vector<clipper::Coord_grid> centres = window_centres(grid, geom, weights, mask, params.mask_threshold);

   // the resolution of each window
   std::vector<double> resolutions(centres.size(), 0.0);
   {
      std::vector<float> map_1 = dense_map(half_1);
      std::vector<float> map_2 = dense_map(half_2);
      std::vector<float> random_1;
      std::vector<float> random_2;
      if (params.phase_randomization) {
         std::mt19937 generator(1); // the same random numbers every time, so that a run can be repeated
         random_1 = phase_randomized(map_1, grid, geom, &generator);
         random_2 = phase_randomized(map_2, grid, geom, &generator);
      }
      const std::vector<float> *random_1_p = params.phase_randomization ? &random_1 : nullptr;
      const std::vector<float> *random_2_p = params.phase_randomization ? &random_2 : nullptr;
      std::vector<unsigned short> shell_table = make_shell_table(geom.n_fft);
      auto task = [&] (std::size_t i) {
         resolutions[i] = window_resolution(map_1, map_2, random_1_p, random_2_p, grid, centres[i], weights, shell_table,
                                            geom, params);
      };
      if (params.parallel_for)
         params.parallel_for(centres.size(), task);
      else
         run_tasks(centres.size(), task, params.n_threads);
   }

   // each window adds its weight, and its weight divided by its resolution, to the grid points that it reaches
   std::size_t n_grid = static_cast<std::size_t>(grid.nu()) * grid.nv() * grid.nw();
   std::vector<float> sum_weights(n_grid, 0.0f);
   std::vector<float> sum_weights_over_resolution(n_grid, 0.0f);
   for (std::size_t iw=0; iw<centres.size(); iw++) {
      float reciprocal_resolution = static_cast<float>(1.0 / resolutions[iw]);
      auto add = [&sum_weights, &sum_weights_over_resolution, reciprocal_resolution]
         (std::size_t map_idx, float weight, int, int, int) {
         sum_weights[map_idx] += weight;
         sum_weights_over_resolution[map_idx] += weight * reciprocal_resolution;
      };
      for_each_window_point(grid, centres[iw], weights, geom, add);
   }

   // What a grid point outside the mask is given. It is not a resolution that a grid point
   // inside the mask can have, and it is on the poor side of them all, not the good side.
   const float outside_the_mask = 99.9f;
   clipper::Xmap<float> xmap_resolution(half_1.spacegroup(), half_1.cell(), grid);
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_resolution.first(); !ix.last(); ix.next()) {
      clipper::Coord_grid c = ix.coord();
      std::size_t map_idx = dense_index(grid, c.u(), c.v(), c.w());
      float resolution = 0.0f; // no window reaches this grid point
      bool inside_mask = mask.empty() || mask[map_idx] > params.mask_threshold;
      if (! inside_mask)
         resolution = outside_the_mask;
      else if (sum_weights_over_resolution[map_idx] > 0.0f)
         resolution = sum_weights[map_idx] / sum_weights_over_resolution[map_idx];
      xmap_resolution[ix] = resolution;
   }
   return xmap_resolution;
}
