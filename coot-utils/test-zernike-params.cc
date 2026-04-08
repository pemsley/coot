/*
 * coot-utils/test-zernike-params.cc
 *
 * Parameter exploration for Zernike local scan
 * Focuses on D_chain only with configurable parameters
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <clipper/ccp4/ccp4_map_io.h>
#include <mmdb2/mmdb_manager.h>
#include "zernike.hh"

void print_usage(const char *prog) {
   std::cout << "Usage: " << prog << " <map> <pdb> [options]\n"
             << "Options:\n"
             << "  --n_grid N     Grid points per dimension (default: 32)\n"
             << "  --order N      Zernike order (default: 24)\n"
             << "  --radius R     Sphere radius in Angstroms (default: 23)\n"
             << "  --range R      Scan range +/- in Angstroms (default: 5)\n"
             << "  --step S       Scan step size in Angstroms (default: 1)\n"
             << "  --mask_radius R  CA mask radius (default: 5)\n"
             << "  --scan_centre X,Y,Z  Centre for scanning (default: use PDB CA centroid)\n"
             << "  --output FILE  Output file (default: stdout)\n"
             << std::endl;
}

int main(int argc, char **argv) {

   if (argc < 3) {
      print_usage(argv[0]);
      return 1;
   }

   std::string map_file = argv[1];
   std::string pdb_file = argv[2];

   // Default parameters
   int n_grid = 32;
   int order = 24;
   float radius = 23.0f;
   float scan_range = 5.0f;
   float scan_step = 1.0f;
   float mask_radius = 5.0f;
   std::string output_file;
   bool use_custom_centre = false;
   double scan_cx = 0, scan_cy = 0, scan_cz = 0;

   // Parse command line options
   for (int i = 3; i < argc; i++) {
      std::string arg = argv[i];
      if (arg == "--n_grid" && i + 1 < argc) {
         n_grid = std::atoi(argv[++i]);
      } else if (arg == "--order" && i + 1 < argc) {
         order = std::atoi(argv[++i]);
      } else if (arg == "--radius" && i + 1 < argc) {
         radius = std::atof(argv[++i]);
      } else if (arg == "--range" && i + 1 < argc) {
         scan_range = std::atof(argv[++i]);
      } else if (arg == "--step" && i + 1 < argc) {
         scan_step = std::atof(argv[++i]);
      } else if (arg == "--mask_radius" && i + 1 < argc) {
         mask_radius = std::atof(argv[++i]);
      } else if (arg == "--output" && i + 1 < argc) {
         output_file = argv[++i];
      } else if (arg == "--scan_centre" && i + 1 < argc) {
         std::string coords = argv[++i];
         if (sscanf(coords.c_str(), "%lf,%lf,%lf", &scan_cx, &scan_cy, &scan_cz) == 3) {
            use_custom_centre = true;
         } else {
            std::cerr << "Error: --scan_centre requires X,Y,Z format\n";
            return 1;
         }
      }
   }

   // Output stream
   std::ostream *out = &std::cout;
   std::ofstream file_out;
   if (!output_file.empty()) {
      file_out.open(output_file);
      if (file_out.is_open()) {
         out = &file_out;
      }
   }

   *out << "# Zernike Parameter Scan - D_chain only\n";
   *out << "# Parameters:\n";
   *out << "#   n_grid: " << n_grid << "\n";
   *out << "#   order: " << order << "\n";
   *out << "#   radius: " << radius << " A\n";
   *out << "#   scan_range: +/- " << scan_range << " A\n";
   *out << "#   scan_step: " << scan_step << " A\n";
   *out << "#   mask_radius: " << mask_radius << " A\n";

   try {
      // Load map
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file);
      file.import_xmap(xmap);
      file.close_read();

      // Load PDB
      mmdb::Manager *mol = new mmdb::Manager();
      mmdb::ERROR_CODE rc = mol->ReadCoorFile(pdb_file.c_str());
      if (rc != mmdb::Error_NoError) {
         std::cerr << "Error reading PDB file: " << pdb_file << std::endl;
         delete mol;
         return 1;
      }

      // Calculate reference centre from CA atom centroid
      int sel_hnd = mol->NewSelection();
      mol->SelectAtoms(sel_hnd, 0, "*",
                       mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                       "*", " CA ", "*", "*");
      mmdb::PPAtom ca_atoms = nullptr;
      int n_ca = 0;
      mol->GetSelIndex(sel_hnd, ca_atoms, n_ca);

      double sum_x = 0, sum_y = 0, sum_z = 0;
      for (int i = 0; i < n_ca; ++i) {
         sum_x += ca_atoms[i]->x;
         sum_y += ca_atoms[i]->y;
         sum_z += ca_atoms[i]->z;
      }
      mol->DeleteSelection(sel_hnd);

      clipper::Coord_orth ref_centre(sum_x / n_ca, sum_y / n_ca, sum_z / n_ca);
      *out << "# Reference centre (from " << n_ca << " CA atoms): "
           << ref_centre.x() << ", " << ref_centre.y() << ", " << ref_centre.z() << "\n";

      // Compute relative CA positions for consistent masking
      std::vector<clipper::Coord_orth> relative_ca_positions;
      sel_hnd = mol->NewSelection();
      mol->SelectAtoms(sel_hnd, 0, "*",
                       mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                       "*", " CA ", "*", "*");
      mol->GetSelIndex(sel_hnd, ca_atoms, n_ca);
      for (int i = 0; i < n_ca; ++i) {
         relative_ca_positions.push_back(clipper::Coord_orth(
            ca_atoms[i]->x - ref_centre.x(),
            ca_atoms[i]->y - ref_centre.y(),
            ca_atoms[i]->z - ref_centre.z()));
      }
      mol->DeleteSelection(sel_hnd);
      *out << "# Computed " << relative_ca_positions.size() << " relative CA positions for masking\n";

      // Scan centre: use custom if specified, otherwise use reference centre
      clipper::Coord_orth scan_centre;
      if (use_custom_centre) {
         scan_centre = clipper::Coord_orth(scan_cx, scan_cy, scan_cz);
         *out << "# Scan centre (custom): " << scan_centre.x() << ", " << scan_centre.y() << ", " << scan_centre.z() << "\n";
      } else {
         scan_centre = ref_centre;
         *out << "# Scan centre (using ref centre): " << scan_centre.x() << ", " << scan_centre.y() << ", " << scan_centre.z() << "\n";
      }

      // Compute reference descriptor WITH protein mask
      *out << "# Computing reference descriptor...\n";
      auto t_start = std::chrono::high_resolution_clock::now();

      coot::ZernikeDescriptor zd_ref(xmap, ref_centre, radius, relative_ca_positions, mask_radius, n_grid, order);

      auto t_ref = std::chrono::high_resolution_clock::now();
      double ref_time = std::chrono::duration<double>(t_ref - t_start).count();
      *out << "# Reference computation time: " << ref_time << " s\n";

      auto diag = zd_ref.get_diagnostics();
      *out << "# Reference diagnostics:\n";
      *out << "#   Points in mask: " << diag.n_points_in_mask << " / " << diag.n_points_in_sphere << "\n";
      *out << "#   Percent in mask: " << (100.0 * diag.n_points_in_mask / diag.n_points_in_sphere) << "%\n";

      // Scan parameters
      int n_steps = static_cast<int>(2 * scan_range / scan_step) + 1;
      int total_points = n_steps * n_steps * n_steps;
      *out << "# Scan grid: " << n_steps << "^3 = " << total_points << " points\n";
      *out << "#\n";
      *out << "# x y z L2_norm cosine_sim dx dy dz dist_from_centre\n";
      out->flush();

      // Scan around scan_centre
      int count = 0;
      auto t_scan_start = std::chrono::high_resolution_clock::now();

      for (float dx = -scan_range; dx <= scan_range + 0.001; dx += scan_step) {
         for (float dy = -scan_range; dy <= scan_range + 0.001; dy += scan_step) {
            for (float dz = -scan_range; dz <= scan_range + 0.001; dz += scan_step) {
               double x = scan_centre.x() + dx;
               double y = scan_centre.y() + dy;
               double z = scan_centre.z() + dz;
               clipper::Coord_orth pt(x, y, z);

               coot::ZernikeDescriptor zd(xmap, pt, radius, relative_ca_positions, mask_radius, n_grid, order);
               float dist = coot::ZernikeDescriptor::distance(zd_ref, zd);
               float cos_sim = coot::ZernikeDescriptor::cosine_similarity(zd_ref, zd);

               float dist_from_centre = std::sqrt(dx*dx + dy*dy + dz*dz);

               *out << x << " " << y << " " << z << " "
                    << dist << " " << cos_sim << " "
                    << dx << " " << dy << " " << dz << " "
                    << dist_from_centre << "\n";

               count++;
               if (count % 100 == 0) {
                  out->flush();
                  auto t_now = std::chrono::high_resolution_clock::now();
                  double elapsed = std::chrono::duration<double>(t_now - t_scan_start).count();
                  double rate = count / elapsed;
                  double remaining = (total_points - count) / rate;
                  std::cerr << "\rProgress: " << count << "/" << total_points
                            << " (" << (100*count/total_points) << "%) "
                            << "ETA: " << remaining << "s   " << std::flush;
               }
            }
         }
      }

      auto t_end = std::chrono::high_resolution_clock::now();
      double scan_time = std::chrono::duration<double>(t_end - t_scan_start).count();

      *out << "#\n";
      *out << "# Scan completed: " << count << " points in " << scan_time << " s\n";
      *out << "# Rate: " << (count / scan_time) << " points/s\n";

      std::cerr << "\nDone. " << count << " points in " << scan_time << " s\n";

      delete mol;

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;
   }

   return 0;
}
