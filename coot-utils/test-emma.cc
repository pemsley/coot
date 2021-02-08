
#include <fstream>
#include <iomanip>

#include <fftw.h>

#include <mmdb2/mmdb_manager.h>
#include <clipper/ccp4/ccp4_map_io.h>

#include "coot-coord-utils.hh"
#include "coot-map-utils.hh"
#include "atom-selection-container.hh"
#include "emma.hh"

void overlap_map(const std::vector<std::complex<double> > &vp_fft,
                 float angstroms_per_bin,
                 const std::string &xmap_file_name) {
   
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(xmap_file_name);
      file.import_xmap(xmap);

      clipper::Resolution reso(4.16);
      clipper::HKL_info myhkl(xmap.spacegroup(), xmap.cell(), reso, true);
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphi_map(myhkl);

      // perhaps I want to use 
      xmap.fft_to(fphi_map);

      const unsigned int max_bin_no = vp_fft.size()/2;

      // a reflection at 1A will select the bin at max_bin_no
      // if angstroms_per_bin is 1.
      //

      clipper::HKL_info::HKL_reference_index hri;
      for (hri = fphi_map.first(); !hri.last(); hri.next()) {
	      float irs = hri.invresolsq();
	      float ir = sqrt(irs);
	      unsigned int bin_idx = static_cast<unsigned int>(ir * max_bin_no / angstroms_per_bin);
	      float bin_res = bin_idx * (1.0/static_cast<float>(max_bin_no));
	      std::cout << hri.hkl().format() << " bin_idx " << bin_idx
		             << " bin resolution " << bin_res
		             << " vs " << ir << " ir(Angstroms) " << 1.0/ir << "\n";

	      // We only need to look at the "low" resolution reflections
	      //
	      if (bin_idx < max_bin_no) {
	         float a = fphi_map[hri].a();
	         float b = fphi_map[hri].b();
	         std::complex<float> fphi(a,b);
	         std::complex<float> saft(vp_fft[bin_idx]);
	         std::complex<float> p = fphi * std::conj(-saft); // or the other way round or negative.
	         float new_f   = std::abs(p);
	         float new_phi = std::arg(p);
	         fphi_map[hri].f()   = new_f;
	         fphi_map[hri].phi() = new_phi;
	      } else {
	         fphi_map[hri].f()   = 0.0;
	         fphi_map[hri].phi() = 0.0;
	      }
      }

      xmap.fft_from(fphi_map);
      clipper::CCP4MAPfile outmapfile;
      outmapfile.open_write("FFT-convoluted.map");
      outmapfile.export_xmap(xmap);
      outmapfile.close_write();

}


int main(int argc, char **argv) {

   int status = 0;

   if (argc > 3) {

      std::string pdb_file_name   = argv[1];
      std::string table_file_name = argv[2];
      std::string xmap_file_name  = argv[3];
      float border = 5.0;

      atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, true, false);
      float angstroms_per_bin = 1.0;

      // I don't like spherically_averaged_molecule now.  Better do this:
      //
      // o Move the molecule to the origin,
      // o find the extents of the molecule and hence the cell
      // o Use clipper::Xmap<float> calc_atom_map(mmdb::Manager *mol,..);
      // o calculate N phi_thetas
      // o sample the points on a sphere: average_of_sample_map_at_sphere_points()
      //

      std::cout << "----------- here with asc n_selected_atoms " << asc.n_selected_atoms << std::endl;
      std::pair<bool, clipper::Coord_orth> centre = coot::centre_of_molecule(asc.mol);
      if (centre.first) {

         std::cout << "DEBUG:: molecule  centre after recentering: "
		             << centre.first << " " << centre.second.format() << std::endl;

	      // move to origin (!)
	      asc.apply_shift(centre.second.x(), centre.second.y(), centre.second.z());

	      std::pair<clipper::Coord_orth, clipper::Coord_orth> e = coot::util::extents(asc.mol, asc.SelectionHandle);
	      double x_range = e.second.x() - e.first.x();
	      double y_range = e.second.y() - e.first.y();
	      double z_range = e.second.z() - e.first.z();

	      float border = 5.0;
	      double nr = clipper::Util::d2rad(90.0);
	      clipper::Cell_descr cell_descr(x_range + 2*border,
					     y_range + 2*border,
					     z_range + 2*border, nr, nr, nr);
	      clipper::Cell cell = clipper::Cell(cell_descr);
	      clipper::Spacegroup spacegroup = clipper::Spacegroup::p1();
	      clipper::Resolution reso = clipper::Resolution(3.0);
	      clipper::Grid_sampling gs(spacegroup, cell, reso);
	      clipper::Xmap<float> xmap =
	         coot::util::calc_atom_map(asc.mol, asc.SelectionHandle, cell, spacegroup, gs);

	      if (true) {
	         clipper::CCP4MAPfile outmapfile;
	         outmapfile.open_write("atom_calc.map");
	         outmapfile.export_xmap(xmap);
	         outmapfile.close_write();
	      }

	 // std::vector<std::pair<double, double> > vp =
	 // coot::util::spherically_averaged_molecule(asc, angstroms_per_bin);
	 // std::vector<std::complex<double> > vp_fft(vp.size());

	 // New method:

	 unsigned int n_bins = 32; // more than 6
	 unsigned int n_sphere_points = 100;

	 std::vector<std::pair<double, double> > vp(n_bins);
	 std::vector<coot::util::phitheta> phithetas = coot::util::make_phi_thetas(n_sphere_points);
	 for (unsigned int i=2; i<n_bins; i++) {
	    float radius = static_cast<float>(i) + 0.5;
	    float average = coot::util::average_of_sample_map_at_sphere_points(centre.second,
									       radius, phithetas, xmap);
	    std::cout << "radius: " << radius << " average " << average << std::endl;
	    vp[i].first = radius;
	    vp[i].second = average;
	 }
	 // smooth out the inner 2 bins - or fft -> low pass filter -> inv-fft
	 double sum_for_smooth = 0.0;
	 for (unsigned int i=2; i<6; i++)
	    sum_for_smooth += vp[i].second;
	 vp[0].second = 0.25 * sum_for_smooth;
	 vp[1].second = 0.25 * sum_for_smooth;

	 std::vector<std::complex<double> > vp_fft(vp.size());

	 std::ofstream f(table_file_name);
	 for (std::size_t i=0; i<vp.size(); i++) {
	    f << "   " << vp[i].first << " " << vp[i].second << std::endl;
	 }
	 f.close();

	 int N = vp.size();

	 {
	    fftw_complex in[N], out[N];
	    fftw_plan p = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE);

	    for (std::size_t i=0; i<vp.size(); i++) {
	       // with 42 bins 3 * i is about as high a frequency as I can get.
	       // in[i].re = sin(1.0*i); // vp[i].first;
	       in[i].re = vp[i].second -26;
	       in[i].im = 0.0;
	    }

	    fftw_one(p, in, out);
	    for (std::size_t i=0; i<vp.size(); i++)
	       std::cout << "  in " << i << " " << in[i].re << " " << in[i].im << std::endl;
	    for (std::size_t i=0; i<vp.size(); i++) {
	       double theta = atan2(out[i].re, out[i].im);
	       double amp = sqrt(out[i].re*out[i].re + out[i].im*out[i].im);
	       std::cout << "  out: bin: " << i
			 << " re " << std::setw(8) << std::setprecision(3) << std::fixed << std::right << out[i].re
			 << " im " << std::setw(8) << std::setprecision(3) << std::fixed << std::right << out[i].im
			 << " ampl " << amp << " phi " << 180.0 * theta/M_PI << std::endl;
	    }

	    for (std::size_t i=0; i<vp.size(); i++)
	       vp_fft[i] = std::complex<double>(out[i].re, out[i].im);

	    overlap_map(vp_fft, angstroms_per_bin, xmap_file_name);

	    fftw_destroy_plan(p);
	 }
      }
   } else {
      std::cout << "Usage: pdb-file-name hist-spherical-table-out map-file-name" << std::endl;
   }
   return status;

}

