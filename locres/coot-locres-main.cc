
// coot-locres: the local resolution of a cryo-EM reconstruction, from its two half-maps.
//
// This reads the maps, calls coot::local_resolution() and writes the result. The
// calculation itself is in libcoot-locres.

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include <clipper/ccp4/ccp4_map_io.h>

#include "coot-locres.hh"

void print_usage() {
   std::cout << "Usage: coot-locres half-map-1 half-map-2 [options]\n"
             << "\n"
             << "Estimate the local resolution of a cryo-EM reconstruction from its two half-maps and\n"
             << "write it as a map: the resolution in Angstroms at each grid point, or 0 where there is\n"
             << "no estimate. The FSC of the half-maps is calculated in spherical windows with soft\n"
             << "edges, across the map, as \"relion_postprocess --locres\" does. The resolution of a\n"
             << "window is where its FSC first falls below the threshold.\n"
             << "\n"
             << "Options:\n"
             << "  --mask file-name      a mask on the grid of the half-maps. The result is 0 where the\n"
             << "                        mask is 0.5 or less, and a window that reaches nowhere else is\n"
             << "                        not calculated. The mask does not change the other values\n"
             << "  --sampling d          the distance between the centres of neighbouring windows, in\n"
             << "                        Angstroms (default: 25)\n"
             << "  --mask-radius r       the radius within which a window has its full weight, in\n"
             << "                        Angstroms (default: half of the sampling)\n"
             << "  --edge-width w        the width of the soft edge of a window, outside that radius,\n"
             << "                        in Angstroms (default: the sampling)\n"
             << "  --fsc-threshold t     the FSC threshold (default: 0.143)\n"
             << "  --threads n           the number of threads (default: one for each core)\n"
             << "  --output file-name    the map to write (default: locres.mrc, beside half-map-1)\n"
             << "  --help                print this and exit\n"
             << "\n"
             << "The half-maps must be in P1, with an orthogonal cell and the same grid spacing along\n"
             << "each axis, as cryo-EM maps are. The grid spacing is taken from the header of the maps.\n"
             << "\n"
             << "Beware: these resolutions are optimistic. Each window multiplies both half-maps by the\n"
             << "same sphere, which raises their FSC, and that is not corrected (RELION corrects it by\n"
             << "randomising phases). It matters little where the density is good and most in weak\n"
             << "density and solvent. A smaller window, or one with a sharper edge, makes it worse.\n";
}

// if this option was given, set *value_p to the number that it was given
template <typename T>
void set_from_option(const std::map<std::string, std::string> &options, const std::string &name, T *value_p) {
   const std::string &text = options.at(name);
   if (text.empty()) return;
   std::istringstream stream(text);
   T value;
   char c;
   if (!(stream >> value) || (stream >> c))
      throw std::runtime_error(name + " cannot be \"" + text + "\"");
   *value_p = value;
}

clipper::Xmap<float> read_map(const std::string &file_name) {
   clipper::Xmap<float> xmap;
   try {
      clipper::CCP4MAPfile file;
      file.open_read(file_name);
      file.import_xmap(xmap);
      file.close_read();
   }
   catch (const clipper::Message_base &exc) {
      throw std::runtime_error("cannot read the map " + file_name);
   }
   return xmap;
}

void write_map(const clipper::Xmap<float> &xmap, const std::string &file_name) {
   try {
      clipper::CCP4MAPfile file;
      file.open_write(file_name);
      file.export_xmap(xmap);
      file.close_write();
   }
   catch (const clipper::Message_base &exc) {
      throw std::runtime_error("cannot write the map " + file_name);
   }
}

// the statistics that RELION prints for its histogram of local resolutions
void print_statistics(const clipper::Xmap<float> &xmap) {
   std::size_t n_grid_points = 0;
   std::size_t n = 0;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      n_grid_points++;
      if (xmap[ix] > 0.0f) n++;
   }
   if (n == 0) {
      std::cout << "None of the " << n_grid_points << " grid points has an estimate" << std::endl;
      return;
   }
   std::vector<float> v;
   v.reserve(n);
   double sum = 0.0;
   double sum_sq = 0.0;
   for (ix = xmap.first(); !ix.last(); ix.next()) {
      if (xmap[ix] > 0.0f) {
         v.push_back(xmap[ix]);
         sum += xmap[ix];
         sum_sq += static_cast<double>(xmap[ix]) * xmap[ix];
      }
   }
   // the value that this fraction of the others are below
   auto quantile = [&v, n] (double fraction) {
      std::size_t idx = static_cast<std::size_t>(fraction * static_cast<double>(n - 1));
      std::nth_element(v.begin(), v.begin() + idx, v.end());
      return v[idx];
   };
   double mean = sum / static_cast<double>(n);
   double variance = std::max(0.0, sum_sq / static_cast<double>(n) - mean * mean);
   std::cout << "Local resolution in Angstroms, at the " << n << " of " << n_grid_points
             << " grid points that have an estimate:" << std::endl;
   std::cout << "   Min: " << quantile(0.0) << " Q1: " << quantile(0.25) << " Median: " << quantile(0.5)
             << " Q3: " << quantile(0.75) << " Max: " << quantile(1.0) << std::endl;
   std::cout << "   Mean: " << mean << " Std: " << std::sqrt(variance) << std::endl;
}

int main(int argc, char **argv) {

   int status = 0;

   try {
      // the options, and the value that each was given
      std::map<std::string, std::string> options = {
         {"--mask", ""}, {"--sampling", ""}, {"--mask-radius", ""}, {"--edge-width", ""},
         {"--fsc-threshold", ""}, {"--threads", ""}, {"--output", ""}};
      std::vector<std::string> file_names;
      for (int i=1; i<argc; i++) {
         std::string arg = argv[i];
         if (arg == "--help" || arg == "-h") {
            print_usage();
            return 0;
         }
         if (arg.empty() || arg[0] != '-') {
            file_names.push_back(arg);
            continue;
         }
         // "--option value" or "--option=value"
         std::string::size_type pos = arg.find('=');
         std::string name = arg.substr(0, pos);
         if (options.find(name) == options.end())
            throw std::runtime_error("unknown option " + name + " (--help lists the options)");
         if (pos != std::string::npos)
            options[name] = arg.substr(pos + 1);
         else if (i + 1 < argc)
            options[name] = argv[++i];
         if (options[name].empty())
            throw std::runtime_error(name + " needs a value");
      }
      if (file_names.size() != 2)
         throw std::runtime_error("two half-maps are needed (--help says how to use this)");

      // as in RELION, the size of a window follows the sampling unless it is given
      coot::locres_params_t params;
      set_from_option(options, "--sampling", &params.sampling_A);
      params.mask_radius_A = 0.5f * params.sampling_A;
      params.edge_width_A = params.sampling_A;
      set_from_option(options, "--mask-radius", &params.mask_radius_A);
      set_from_option(options, "--edge-width", &params.edge_width_A);
      set_from_option(options, "--fsc-threshold", &params.fsc_threshold);
      set_from_option(options, "--threads", &params.n_threads);
      if (params.n_threads < 0)
         throw std::runtime_error("--threads cannot be negative");

      std::string output_file_name = options["--output"];
      if (output_file_name.empty()) {
         std::string::size_type pos = file_names[0].find_last_of("/\\");
         if (pos != std::string::npos)
            output_file_name = file_names[0].substr(0, pos + 1);
         output_file_name += "locres.mrc";
      }

      clipper::Xmap<float> half_1 = read_map(file_names[0]);
      clipper::Xmap<float> half_2 = read_map(file_names[1]);
      clipper::Xmap<float> xmap_mask;
      bool have_mask = ! options["--mask"].empty();
      if (have_mask)
         xmap_mask = read_map(options["--mask"]);

      clipper::Grid_sampling grid = half_1.grid_sampling();
      std::cout << "Half-maps: " << grid.nu() << " x " << grid.nv() << " x " << grid.nw() << " grid points, "
                << half_1.cell().a() / static_cast<double>(grid.nu()) << " A apart" << std::endl;
      std::cout << "Windows: " << params.sampling_A << " A apart, full weight within " << params.mask_radius_A
                << " A, soft edge " << params.edge_width_A << " A wide" << std::endl;

      clipper::Xmap<float> xmap = coot::local_resolution(half_1, half_2, have_mask ? &xmap_mask : nullptr, params);
      print_statistics(xmap);
      write_map(xmap, output_file_name);
      std::cout << "Wrote " << output_file_name << std::endl;
   }
   catch (const std::exception &e) {
      std::cerr << "coot-locres: " << e.what() << std::endl;
      status = 1;
   }
   return status;
}
