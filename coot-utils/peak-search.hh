
#ifndef HAVE_VECTOR
#include <vector>
#endif

#include "mmdb_manager.h"
#include "clipper/core/coords.h"
#include "clipper/core/xmap.h"
#include "clipper/contrib/skeleton.h"

namespace coot { 

   class peak_search { 
      float map_rms;
      clipper::Coord_orth move_grid_to_peak(const clipper::Xmap<float> &xmap,
					    const clipper::Coord_grid &c_g); 
      void peak_search_0(const clipper::Xmap<float> &xmap,
			 clipper::Xmap<short int> *marked_map_p,
			 float n_sigma) const;
      void peak_search_0_negative(const clipper::Xmap<float> &xmap,
				  clipper::Xmap<short int> *marked_map_p,
				  float n_sigma);
      void peak_search_1(const clipper::Xmap<float> &xmap,
			 clipper::Xmap<short int> *marked_map_p);
      void peak_search_2(const clipper::Xmap<float> &xmap,
			 clipper::Xmap<short int> *marked_map_p);

      void peak_search_trace_along(const clipper::Coord_grid &c_g_start, 
				   const clipper::Skeleton_basic::Neighbours &neighb,
				   const clipper::Xmap<float> &xmap,
				   clipper::Xmap<short int> *marked_map);
      void mask_around_coord(clipper::Xmap<float> *xmap,
			     const clipper::Coord_orth &co, float atom_radius) const;
      std::vector<clipper::Coord_orth>
      make_sample_protein_coords(CMMDBManager *mol) const;
      clipper::Coord_orth
      move_point_close_to_protein(const clipper::Coord_orth &pt,
				  const std::vector<clipper::Coord_orth> &protein,
				  const std::vector<int> &itrans,
				  const clipper::Xmap<float> &xmap) const;

      const std::vector<int>
      find_protein_to_origin_translations(const std::vector<clipper::Coord_orth> &sampled_protein_coords, const clipper::Xmap<float> &xmap) const;
      
      double
      min_dist_to_protein(const clipper::Coord_orth &point,
			  const std::vector<clipper::Coord_orth> &sampled_protein_coords) const;

      static bool compare_ps_peaks(const std::pair<clipper::Coord_orth, float> &a,
				   const std::pair<clipper::Coord_orth, float> &b);
      static bool compare_ps_peaks_mri(const std::pair<clipper::Xmap<float>::Map_reference_index, float> &a,
				       const std::pair<clipper::Xmap<float>::Map_reference_index, float> &b);
      static bool compare_ps_peaks_cg(const std::pair<clipper::Coord_grid, float> &a,
				      const std::pair<clipper::Coord_grid, float> &b);

   public:
      peak_search(const clipper::Xmap<float> &xmap);
      std::vector<clipper::Coord_orth>
      get_peaks(const clipper::Xmap<float> &xmap,
		float n_sigma);
      std::vector<std::pair<clipper::Xmap<float>::Map_reference_index, float> >
      get_peak_map_indices(const clipper::Xmap<float> &xmap,
			   float n_sigma) const;
      std::vector<std::pair<clipper::Coord_grid, float> >
      get_peak_grid_points(const clipper::Xmap<float> &xmap,
			   float n_sigma) const;
      std::vector<std::pair<clipper::Coord_orth, float> >
      get_peaks(const clipper::Xmap<float> &xmap,
		CMMDBManager *mol, 
		float n_sigma,
		int do_positive_levels_flag,
		int also_negative_levels_flag);
      std::vector<std::pair<clipper::Coord_orth, float> >
      get_peaks(const clipper::Xmap<float> &xmap,
		float n_sigma,
		int do_positive_levels_flag,
		int also_negative_levels_flag);
      
      void mask_map(clipper::Xmap<float> *xmap,
		    const std::vector<clipper::Coord_orth> &ps_peaks) const;

      void add_peak_vectors(std::vector<clipper::Coord_orth> *in,
			    const std::vector<clipper::Coord_orth> &extras) const;

   };
}

