
#ifndef COOT_LOCRES_HH
#define COOT_LOCRES_HH

// The local resolution of a cryo-EM reconstruction, from its two half-maps.
//
// The Fourier Shell Correlation (FSC) of the half-maps is calculated in small
// spherical windows with soft edges, on a regular lattice across the map, as
// "relion_postprocess --locres" does.
//
// Beware: multiplying both half-maps by the same window makes their FSC too high
// at high resolution, so the resolutions from this library are optimistic. This
// is not corrected (RELION corrects it by randomising phases). A smaller window,
// or one with a sharper edge, makes it worse.
//
// This library uses clipper and the standard library only.

#include <cstddef>
#include <functional>

#include <clipper/core/xmap.h>

namespace coot {

   //! \brief The parameters of coot::local_resolution()
   //!
   //! The defaults are those of RELION with "--locres_sampling 25", where the mask
   //! radius is half of the sampling and the edge width is equal to it.
   class locres_params_t {
   public:
      //! the distance between the centres of neighbouring windows, in Angstroms
      float sampling_A;
      //! the radius within which a window has its full weight, in Angstroms
      float mask_radius_A;
      //! the width of the soft edge of a window, outside the mask radius, in Angstroms
      float edge_width_A;
      //! the resolution of a window is where its FSC first falls below this
      float fsc_threshold;
      //! the resolution, in Angstroms, of a window that has no signal, and the worst resolution that is reported
      float low_resolution_limit_A;
      //! a grid point is inside the mask if the mask is greater than this there
      float mask_threshold;
      //! the number of threads to use, 0 means one for each core
      int n_threads;
      //! \brief A function that runs the tasks, to use in place of the threads of this library
      //!
      //! If this is set, it must call task(i) once for each i from 0 to n - 1, in any order
      //! and on any threads, and return when every call has returned. n_threads is then not
      //! used. Use it to run the calculation on a thread pool that already exists.
      std::function<void(std::size_t n, const std::function<void(std::size_t)> &task)> parallel_for;

      locres_params_t() :
         sampling_A(25.0f), mask_radius_A(12.5f), edge_width_A(25.0f), fsc_threshold(0.143f),
         low_resolution_limit_A(50.0f), mask_threshold(0.5f), n_threads(0) {}
   };

   //! \brief Estimate the local resolution of a reconstruction from its two half-maps
   //!
   //! The windows are centred on the grid points whose indices are multiples of the
   //! sampling, inside the ellipsoid that fits in the unit cell less the mask radius
   //! (a sphere, for a cubic cell). The resolution of a window is where the FSC of the
   //! two half-maps in that window first falls below the threshold. The resolution at a
   //! grid point is the reciprocal of the mean of 1/resolution of the windows that reach
   //! it, each weighted by the weight of the window at that grid point.
   //!
   //! Beware: these resolutions are optimistic, because the half-maps are multiplied by
   //! the same window and that is not corrected. See the top of this file.
   //!
   //! A mask does not change the values inside it. It saves time, because a window that
   //! reaches no grid point inside the mask is not calculated, and the grid points
   //! outside it are set to 0.
   //!
   //! This function throws a std::runtime_error if the maps or the parameters cannot be
   //! used. It does not write to standard output.
   //!
   //! \param half_1 the first half-map, in P1, with an orthogonal cell and the same grid
   //!        spacing along each axis (as cryo-EM maps have)
   //! \param half_2 the second half-map, with the same cell and grid as the first
   //! \param xmap_mask_p a mask with the same cell and grid as the half-maps, or nullptr
   //!        for no mask
   //! \param params the size and spacing of the windows, the thresholds and the threads
   //! \return the resolution in Angstroms at each grid point, on the grid of the half-maps.
   //!         0 means that there is no estimate: no window reaches that grid point, or it
   //!         is outside the mask. Other values are between twice the grid spacing (the FSC
   //!         did not fall below the threshold) and the low resolution limit.
   clipper::Xmap<float> local_resolution(const clipper::Xmap<float> &half_1,
                                         const clipper::Xmap<float> &half_2,
                                         const clipper::Xmap<float> *xmap_mask_p = nullptr,
                                         const locres_params_t &params = locres_params_t());
}

#endif // COOT_LOCRES_HH
