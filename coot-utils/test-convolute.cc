
#include <fstream>
#include <iomanip>
#include <thread>

#include <clipper/ccp4/ccp4_map_io.h>

#include "coot-coord-utils.hh"
#include "coot-map-utils.hh"
#include "coot-map-heavy.hh"
#include "atom-selection-container.hh"

void thread_fill_a_map() {

   int n_threads = 5;

   clipper::Xmap<float> xmap; // not blank in real code

   std::string file_name("test.map");
   clipper::CCP4MAPfile file;
   try {
      file.open_read(file_name);
      file.import_xmap(xmap);
   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << file_name << std::endl;
   }


   typedef clipper::Xmap<float>::Map_reference_index MRI;

   // this function does't work - needs some debugging.
   std::vector<std::pair<MRI, MRI> > map_ref_start_stops =
      coot::make_map_reference_index_start_stops(xmap, n_threads);

   for (unsigned int i=0; i<map_ref_start_stops.size(); i++) {
      std::cout << "checking returned start stops:    " << i << " "
		<< map_ref_start_stops[i].first.coord().format() << " "
		<< map_ref_start_stops[i].first.index() << " "
		<< map_ref_start_stops[i].second.coord().format() << " "
		<< map_ref_start_stops[i].second.index() << "\n";
   }

   std::vector<std::thread> threads;

   unsigned int i_count_basic = 0;
   unsigned int i_count_from_start_stops = 0;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next())
      xmap[ix] = 0.0;

   for (ix = xmap.first(); !ix.last(); ix.next())
      i_count_basic += 1;
   for (unsigned int i=0; i<map_ref_start_stops.size(); i++) {
      const std::pair<MRI, MRI> &ss = map_ref_start_stops[i];      
      for (MRI ix = ss.first; ix.index() != ss.second.index(); ix.next())
	 i_count_from_start_stops += 1;
   }

   unsigned int n_zero = 0;
   unsigned int n_one  = 0;
   unsigned int n_greater = 0;
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      if (xmap[ix] < 0.5)
	 n_zero++;
      else
	 if (xmap[ix] < 1.001) {
	    n_greater++;
	 } else {
	    n_one++;
	 }
   }
   std::cout << "DEBUG:: counts n_zero: " << n_zero << " " << n_one << " " << n_one << std::endl;
   

   std::cout << "compare counts: basic: " << i_count_basic << " vs " << i_count_from_start_stops
	     << std::endl;

   // this is for editing the xmap. If you don't want to do that,
   // pass it a const and use std::cref() in the calling function.
   // This (l) is a lambda function.
   auto l = [](clipper::Xmap<float> &xmap, const std::pair<MRI, MRI> &ss) {
      std::cout << "info:: start: " << ss.first.coord().format() << " " << ss.first.index()
      << "   end " << ss.second.coord().format() << " " << ss.second.index() << std::endl;
      for (MRI ix = ss.first; ix.index() != ss.second.index(); ix.next()) {
         xmap[ix] += 1.0; // or whatever
      }
   };
   for (unsigned int i=0; i<map_ref_start_stops.size(); i++) {
      std::cout << "making a thread with start stops: "
		<< map_ref_start_stops[i].first.coord().format() << " "
		<< map_ref_start_stops[i].first.index() << " "
		<< map_ref_start_stops[i].second.coord().format() << " "
		<< map_ref_start_stops[i].second.index() << "\n";
      threads.push_back(std::thread(l, std::ref(xmap), std::cref(map_ref_start_stops[i])));
   }
   for (unsigned int i=0; i<map_ref_start_stops.size(); i++)
      threads[i].join();

   // if needed.
   unsigned int max_count = 0;
   if (!map_ref_start_stops.empty())
      max_count = map_ref_start_stops.back().second.index();

   clipper::CCP4MAPfile outmapfile;
   outmapfile.open_write("thread-filled-testing-out.map");
   outmapfile.export_xmap(xmap);
   outmapfile.close_write();

}

float get_spherically_averaged_density_value(const std::vector<std::pair<double, double> > &sam,
					     float angstroms_per_bin,
					     const float &dist);

clipper::Xmap<float>
convolute_map(const std::vector<std::pair<double, double> > &sam,
	      float angstroms_per_bin,
	      const clipper::Xmap<float> &xmap) {

   clipper::Xmap<float> r = xmap;
   double mol_radius = 26;
   double dd_crit = mol_radius * mol_radius;

   unsigned int count = 0;
   clipper::Xmap_base::Map_reference_index ix_1, ix_2;;
   for (ix_1 = xmap.first(); !ix_1.last(); ix_1.next()) {
      r[ix_1] = 0.0;
      count++;
   }

   std::cout << "Map count " << count << std::endl;
#if 0 // single threaded
   count = 0;
   for (ix_1 = xmap.first(); !ix_1.last(); ix_1.next()) {
      clipper::Coord_map  cm_1 = ix_1.coord().coord_map();
      clipper::Coord_frac cf_1 = cm_1.coord_frac(xmap.grid_sampling());
      clipper::Coord_orth co_1 = cf_1.coord_orth(xmap.cell());
      for (ix_2 = xmap.first(); !ix_2.last(); ix_2.next()) {
	 clipper::Coord_map  cm_2 = ix_2.coord().coord_map();
	 clipper::Coord_frac cf_2 = cm_2.coord_frac(xmap.grid_sampling());
	 clipper::Coord_orth co_2 = cf_2.coord_orth(xmap.cell());
	 double dd = (co_2-co_1).lengthsq();
	 if (dd < dd_crit) {
	    float d = static_cast<float>(std::sqrt(dd));
	    float cf = get_spherically_averaged_density_value(sam, angstroms_per_bin, d);
	    r[ix_2] += cf;
	 }
      }
      count++;
      if (count%100 == 0) {
	 std::cout << count << " ";
	 std::cout.flush();
      }
   }
#endif

   unsigned int n_threads = 4;
   typedef clipper::Xmap<float>::Map_reference_index MRI;
   std::vector<std::pair<MRI, MRI> > map_ref_start_stops =
      coot::make_map_reference_index_start_stops(xmap, n_threads);
   auto l = [](const clipper::Xmap<float> &xmap,
	       const std::pair<MRI, MRI> &ss,
	       const std::vector<std::pair<double, double> > &sam,
	       float angstroms_per_bin, clipper::Xmap<float> *r_xmap_p) {
      std::cout << "info:: start: " << ss.first.coord().format() << " " << ss.first.index()
                << "   end " << ss.second.coord().format() << " " << ss.second.index() << std::endl;
      double mol_radius = 26;
      double dd_crit = mol_radius * mol_radius;
      for (MRI ix_1 = ss.first; ix_1.index() != ss.second.index(); ix_1.next()) {
	 clipper::Coord_map  cm_1 = ix_1.coord().coord_map();
	 clipper::Coord_frac cf_1 = cm_1.coord_frac(xmap.grid_sampling());
	 clipper::Coord_orth co_1 = cf_1.coord_orth(xmap.cell());
	 for (MRI ix_2 = xmap.first(); !ix_2.last(); ix_2.next()) {
	    clipper::Coord_map  cm_2 = ix_2.coord().coord_map();
	    clipper::Coord_frac cf_2 = cm_2.coord_frac(xmap.grid_sampling());
	    clipper::Coord_orth co_2 = cf_2.coord_orth(xmap.cell());
	    double dd = (co_2-co_1).lengthsq();
	    if (dd < dd_crit) {
	       float d = static_cast<float>(std::sqrt(dd));
	       float cf = get_spherically_averaged_density_value(sam, angstroms_per_bin, d);
	       (*r_xmap_p)[ix_1] += cf;
	    }
	 }
      }
   };
   std::vector<std::thread> threads;
   for (unsigned int i=0; i<map_ref_start_stops.size(); i++) {
      const std::pair<MRI, MRI> &ss = map_ref_start_stops[i];
      std::cout << "making a thread with start stops: "
		<< map_ref_start_stops[i].first.coord().format() << " "
		<< map_ref_start_stops[i].first.index() << " "
		<< map_ref_start_stops[i].second.coord().format() << " "
		<< map_ref_start_stops[i].second.index() << "\n";
      threads.push_back(std::thread(l, std::cref(xmap), std::cref(map_ref_start_stops[i]),
				    std::cref(sam), angstroms_per_bin, &r));
   }
   for (unsigned int i=0; i<threads.size(); i++)
      threads[i].join();

   std::cout << std::endl;

   return r;
}

float
get_spherically_averaged_density_value(const std::vector<std::pair<double, double> > &sam,
				       float angstroms_per_bin,
				       const float &dist) {

   float v = 0.0;

   double radius_max = sam.back().first + 0.5; // that should depend on angstroms_per_bin
   if (dist < radius_max) {
      int n_bins = sam.size();
      int bin_nearest = static_cast<int>((dist/radius_max) * n_bins +0.5);
      if (bin_nearest < n_bins)
	 v = sam[bin_nearest].second;
      else
	 std::cout << "error " << bin_nearest << std::endl;
   }
   return v;
}


int main(int argc, char **argv) {

   int status = 0;

   // thread_fill_a_map();
   // return 0;

   if (argc > 3) {

      std::string pdb_file_name   = argv[1];
      std::string table_file_name = argv[2];
      std::string xmap_file_name  = argv[3];
      float border = 40.0;

      atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
      float angstroms_per_bin = 1.0;

      std::cout << "asc n_selected_atoms " << asc.n_selected_atoms << std::endl;
      std::pair<bool, clipper::Coord_orth> centre = coot::centre_of_molecule(asc.mol);
      if (centre.first) {

	 std::cout << "DEBUG:: molecule  centre after recentering: "
		   << centre.first << " " << centre.second.format() << std::endl;

	 // move to origin (!)
	 asc.apply_shift(centre.second.x(), centre.second.y(), centre.second.z());

	 std::pair<clipper::Coord_orth, clipper::Coord_orth> e =
	    coot::util::extents(asc.mol, asc.SelectionHandle);
	 double x_range = e.second.x() - e.first.x();
	 double y_range = e.second.y() - e.first.y();
	 double z_range = e.second.z() - e.first.z();

	 double nr = clipper::Util::d2rad(90.0);
	 clipper::Cell_descr cell_descr(x_range + 2*border,
					y_range + 2*border,
					z_range + 2*border, nr, nr, nr);
	 clipper::Cell cell = clipper::Cell(cell_descr);
	 clipper::Spacegroup spacegroup = clipper::Spacegroup::p1();
	 clipper::Resolution reso = clipper::Resolution(3.0);
	 clipper::Grid_sampling gs(spacegroup, cell, reso);
	 clipper::Xmap<float> xmap_calc =
	    coot::util::calc_atom_map(asc.mol, asc.SelectionHandle, cell, spacegroup, gs);

	 if (true) {
	    clipper::CCP4MAPfile outmapfile;
	    outmapfile.open_write("atom_calc.map");
	    outmapfile.export_xmap(xmap_calc);
	    outmapfile.close_write();
	 }

	 // New method:

	 unsigned int n_bins = 30; // more than 6
	 unsigned int n_sphere_points = 10000;
	 double molecule_radius = 30.0; // was 26
	 double one_over_n_bins = 1.0/static_cast<double>(n_bins);

	 std::vector<std::pair<double, double> > vp(n_bins);
	 std::vector<coot::util::phitheta> phithetas = coot::util::make_phi_thetas(n_sphere_points);
	 for (unsigned int i=2; i<n_bins; i++) {
	    float radius = molecule_radius * one_over_n_bins * (static_cast<float>(i) + 0.5);
	    float average = coot::util::average_of_sample_map_at_sphere_points(centre.second,
									       radius, phithetas,
									       xmap_calc);
	    vp[i].first = radius;
	    vp[i].second = average;
	 }

	 // smooth out the inner 2 bins - or fft -> low pass filter -> inv-fft

	 // first give them x values
	 for (unsigned int i=0; i<2; i++) {
	    float radius = molecule_radius * one_over_n_bins * (static_cast<float>(i) + 0.5);
	    vp[i].first = radius;
	 }

	 double sum_for_smooth = 0.0;
	 for (unsigned int i=2; i<6; i++) {
	    sum_for_smooth += vp[i].second;
	 }
	 vp[0].second = 0.25 * sum_for_smooth;
	 vp[1].second = 0.25 * sum_for_smooth;

	 std::vector<std::complex<double> > vp_fft(vp.size());

	 std::ofstream f(table_file_name);
	 for (std::size_t i=0; i<vp.size(); i++) {
	    f << "   " << vp[i].first << " " << vp[i].second << std::endl;
	 }
	 f.close();

	 // Read in the reference map to convolute
	 clipper::Xmap<float> xmap;
	 try {
	    clipper::CCP4MAPfile file;
	    file.open_read(xmap_file_name);
	    clipper::Grid_sampling fgs = file.grid_sampling();
	    file.import_xmap(xmap);
	 }
	 catch (const clipper::Message_base &exc) {
	    std::cout << "WARNING:: failed to open " << xmap_file_name << std::endl;
	 }
	 clipper::Xmap<float> convoluted_map = convolute_map(vp, angstroms_per_bin, xmap);

	 clipper::CCP4MAPfile outmapfile;
	 outmapfile.open_write("slow-convolute.map");
	 outmapfile.export_xmap(convoluted_map);
	 outmapfile.close_write();
      }
   }
   return status;

}
