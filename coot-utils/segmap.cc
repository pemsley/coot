
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <queue>
#include <set>

#include "peak-search.hh"
#include "xmap-stats.hh"
#include "segmap.hh"

void
coot::segmap::proc() {

   // std::pair<float, float> mv = mean_and_variance(xmap);
   mean_and_variance<float> map_stats = map_density_distribution(xmap, 1000, false, true);

   // map stats should contain only those parts of the map that are unmasked.

   {
      std::ofstream f("map.hist");
      for (std::size_t i=0; i<map_stats.size(); i++) {
	 float this_level = static_cast<float>(i) * map_stats.bin_width + map_stats.min_density;
	 f << this_level << " " << map_stats.bins[i] << "\n";
      }
   }

   // what contour level do we want? The level that selects 20% of the map (20% is above
   // the given contour level (Terwilliger says 20%)

   float mean = map_stats.mean;
   unsigned int n_count = 0;
   unsigned int n_count_above_mean = 0;
   clipper::Xmap_base::Map_reference_index ix;
   float plausibly_protein_level = mean + sqrt(map_stats.variance);
   for (ix=xmap.first(); !ix.last(); ix.next()) {
      n_count++;
      if (xmap[ix] > plausibly_protein_level)
	 n_count_above_mean++;
   }


   float frac_plausibly_protein = static_cast<float>(n_count_above_mean)/static_cast<float>(n_count);

   std::cout << " Of " << n_count << " points " << frac_plausibly_protein
	     << " were above the plausibly_protein_level" << std::endl;

   float level_frac = 0.2; // was 0.2
   // level_frac = some_func(plausibly_protein_level) // this function will need turning for other proteins
   level_frac = 0.003; // OK for EMD-3908
   level_frac = 0.01;

   int i_top = map_stats.size() -1;
   unsigned int n_in_histogram = 0;
   for (int i=i_top; i>=0; i--) {
      float this_level = static_cast<float>(i) * map_stats.bin_width + map_stats.min_density;
      n_in_histogram += map_stats.bins[i];
   }

   float contour_level = 999.9;
   float v_cumulation = 0.0;
   for (int i=i_top; i>=0; i--) {
      float this_level = static_cast<float>(i) * map_stats.bin_width + map_stats.min_density;
      float v = map_stats.bins[i];
      v_cumulation += v/static_cast<float>(n_in_histogram);
      std::cout << "i " << i << " in this bin: " << v << " this_level " << this_level
		<< " running_sum frac: " << v_cumulation << std::endl;
      if (v_cumulation > level_frac) {
	 contour_level = this_level;
	 break;
      }
   }

   std::cout << "contour-level: " << contour_level << std::endl;

   // This might need to be also (or instead) a function of mean and standard deviation
   //
   float cut_off_for_peak_search = 0.4 * map_stats.max_density;

   std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > > peaks = find_peaks(cut_off_for_peak_search);
   clipper::Xmap<float> m = flood_from_peaks(peaks, contour_level);

   // write results
   std::string map_file_name("segmented.map");
   clipper::CCP4MAPfile file;
   try {
      file.open_write(map_file_name);
      file.export_xmap(m);
   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
   }
}

std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > >
coot::segmap::find_peaks(float cut_off) const {

   clipper::Xmap_base::Map_reference_index ix;
   std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > > peaks;

   clipper::Skeleton_basic::Neighbours neighb(xmap, 0.25, 1.75); // 3x3x3 cube, not centre

   for (ix = xmap.first(); !ix.last(); ix.next())  {
      const float &v = xmap[ix];
      if (v > cut_off) {
	 bool is_peak = true;
	 for (int i=0; i<neighb.size(); i++) {
	    clipper::Coord_grid c_g(ix.coord() + neighb[i]);
	    if (v < xmap.get_data(c_g)) {
	       is_peak = false;
	       break;
	    }
	 }
	 if (is_peak) {
	    std::pair<clipper::Xmap_base::Map_reference_index, float> p(ix, v);
	    peaks.push_back(p);
	 }
      }
   }

   auto sorter = [] (const std::pair<clipper::Xmap_base::Map_reference_index, float> &p1,
		     const std::pair<clipper::Xmap_base::Map_reference_index, float> &p2) {
      return (p1.second > p2.second);
   };

   std::sort(peaks.begin(), peaks.end(), sorter);

   std::cout << "peaks: size " << peaks.size() << std::endl;

   int top_n_peaks_for_show = 5;
   if (int(peaks.size()) < top_n_peaks_for_show) top_n_peaks_for_show = peaks.size();

   for (int i=0; i<top_n_peaks_for_show; i++)
      std::cout << "   " << peaks[i].first.coord().format() << " " << peaks[i].second << std::endl;

   return peaks;
}

clipper::Xmap<float>
coot::segmap::flood_from_peaks(const std::vector<std::pair<clipper::Xmap_base::Map_reference_index, float > > &peaks,
			       float cut_off_for_flood) {

   int top_n_peaks = 1;
   float low_level = 0.0; // if density values are below the cut_off_for_flood and more than zero, make them this

   if (int(peaks.size()) > top_n_peaks) top_n_peaks = peaks.size();
   std::cout << "debug:: in flood_from_peaks top_n_peaks is " << top_n_peaks << std::endl;

   // std::vector<clipper::Coord_grid> inside_points; // use a set of sortable coord grids
   std::queue<clipper::Coord_grid> q;
   clipper::Xmap<int> queued; // bool type doesn't work with the below for-loop
   clipper::Xmap<int> considered;
   clipper::Xmap<int> inside_points;
   inside_points.init(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   considered.init(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   queued.init(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix=considered.first(); !ix.last(); ix.next())
      considered[ix] = 0;
   for (ix=queued.first(); !ix.last(); ix.next())
      queued[ix] = 0;
   for (ix=inside_points.first(); !ix.last(); ix.next())
      inside_points[ix] = 0;
   for (int j=0; j<top_n_peaks; j++) {
      const clipper::Xmap_base::Map_reference_index &ix = peaks[j].first;
      q.push(ix.coord());
      while (! q.empty()) {
	 clipper::Coord_grid gp_new_centre = q.front();
	 considered.set_data(gp_new_centre, 1);
	 q.pop();
	 clipper::Skeleton_basic::Neighbours neighb(xmap, 0.25, 1.75); // 3x3x3 cube, not centre
	 for (int i=0; i<neighb.size(); i++) {
	    clipper::Coord_grid c_g(gp_new_centre + neighb[i]);
	    if (xmap.get_data(c_g) > cut_off_for_flood) {
	       if (considered.get_data(c_g) == 0) {
		  // std::cout << "pushing back " << c_g.format() << " " << std::endl;
		  if (queued.get_data(c_g) == 0) {
		     q.push(c_g);
		     inside_points.set_data(c_g, 1);
		     queued.set_data(c_g, 1);
		  }
	       }
	    }
	 }
      }
   }

   clipper::Xmap<float> segmented = xmap;
   for (ix=inside_points.first(); !ix.last(); ix.next()) {
      if (inside_points[ix] == 0) {
	 if (segmented[ix] > 0.0)
	    segmented[ix] = low_level;
      }
   }

   return segmented;
}
