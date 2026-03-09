/*
 * coot-utils/test-zernike-search.cc
 *
 * Zernike-guided domain search with local fitting refinement
 *
 * 1. Compute Zernike descriptor from calculated density (probe)
 * 2. Scan experimental map for best matches
 * 3. Take top N peaks
 * 4. For each peak, translate model and do local rigid body fit
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/contrib/edcalc.h>
#include <mmdb2/mmdb_manager.h>
#include "zernike.hh"
#include "coot-map-heavy.hh"
#include "mini-mol/mini-mol.hh"
#include "ligand/ligand.hh"

struct ScanResult {
   double x, y, z;
   float l2_norm;
   float cosine_sim;
   float dx, dy, dz;
   float dist_from_centre;
};

bool compare_by_l2(const ScanResult &a, const ScanResult &b) {
   return a.l2_norm < b.l2_norm;
}

void print_usage(const char *prog) {
   std::cout << "Usage: " << prog << " <pdb-file> --exp_map <map-file> [options]\n"
             << "\n"
             << "Zernike-guided domain search with local fitting\n"
             << "\n"
             << "Options:\n"
             << "  --n_grid N       Grid points per dimension (default: 32)\n"
             << "  --order N        Zernike order (default: 20)\n"
             << "  --radius R       Sphere radius in Angstroms (default: 23)\n"
             << "  --grid_spacing S Map grid spacing in Angstroms (default: 1.0)\n"
             << "  --exp_map FILE   Experimental map (required)\n"
             << "  --scan_centre X,Y,Z  Centre for scanning (default: map centre)\n"
             << "  --range R        Scan range +/- in Angstroms (default: 45)\n"
             << "  --step S         Scan step size in Angstroms (default: 3)\n"
             << "  --n_peaks N      Number of top peaks to refine (default: 20)\n"
             << "  --n_trials N     Jiggle fit trials per peak (default: 100)\n"
             << "  --jiggle_scale F Jiggle scale factor (default: 1.0)\n"
             << "  --output FILE    Output file for results\n"
             << std::endl;
}

int main(int argc, char **argv) {

   if (argc < 2) {
      print_usage(argv[0]);
      return 1;
   }

   std::string pdb_file = argv[1];

   // Default parameters
   int n_grid = 32;
   int order = 20;
   float radius = 23.0f;
   float grid_spacing = 1.0f;
   std::string exp_map_file;
   bool use_custom_centre = false;
   double scan_cx = 0, scan_cy = 0, scan_cz = 0;
   float scan_range = 45.0f;
   float scan_step = 3.0f;
   int n_peaks = 20;
   int n_trials = 100;
   float jiggle_scale = 1.0f;
   std::string output_file;

   // Parse options
   for (int i = 2; i < argc; i++) {
      std::string arg = argv[i];
      if (arg == "--n_grid" && i + 1 < argc) {
         n_grid = std::atoi(argv[++i]);
      } else if (arg == "--order" && i + 1 < argc) {
         order = std::atoi(argv[++i]);
      } else if (arg == "--radius" && i + 1 < argc) {
         radius = std::atof(argv[++i]);
      } else if (arg == "--grid_spacing" && i + 1 < argc) {
         grid_spacing = std::atof(argv[++i]);
      } else if (arg == "--exp_map" && i + 1 < argc) {
         exp_map_file = argv[++i];
      } else if (arg == "--scan_centre" && i + 1 < argc) {
         std::string coords = argv[++i];
         if (sscanf(coords.c_str(), "%lf,%lf,%lf", &scan_cx, &scan_cy, &scan_cz) == 3) {
            use_custom_centre = true;
         } else {
            std::cerr << "Error: --scan_centre requires X,Y,Z format\n";
            return 1;
         }
      } else if (arg == "--range" && i + 1 < argc) {
         scan_range = std::atof(argv[++i]);
      } else if (arg == "--step" && i + 1 < argc) {
         scan_step = std::atof(argv[++i]);
      } else if (arg == "--n_peaks" && i + 1 < argc) {
         n_peaks = std::atoi(argv[++i]);
      } else if (arg == "--n_trials" && i + 1 < argc) {
         n_trials = std::atoi(argv[++i]);
      } else if (arg == "--jiggle_scale" && i + 1 < argc) {
         jiggle_scale = std::atof(argv[++i]);
      } else if (arg == "--output" && i + 1 < argc) {
         output_file = argv[++i];
      }
   }

   if (exp_map_file.empty()) {
      std::cerr << "Error: --exp_map is required\n";
      print_usage(argv[0]);
      return 1;
   }

   std::cout << "=== Zernike Domain Search ===" << std::endl;
   std::cout << "Parameters:\n"
             << "  n_grid: " << n_grid << "\n"
             << "  order: " << order << "\n"
             << "  radius: " << radius << " A\n"
             << "  scan_range: +/- " << scan_range << " A\n"
             << "  scan_step: " << scan_step << " A\n"
             << "  n_peaks: " << n_peaks << "\n"
             << "  n_trials: " << n_trials << "\n"
             << "  jiggle_scale: " << jiggle_scale << "\n"
             << std::endl;

   // ============ Load PDB and compute calculated density ============
   mmdb::Manager *mol = new mmdb::Manager();
   mmdb::ERROR_CODE rc = mol->ReadCoorFile(pdb_file.c_str());
   if (rc != mmdb::Error_NoError) {
      std::cerr << "Error reading PDB file: " << pdb_file << std::endl;
      delete mol;
      return 1;
   }
   std::cout << "Loaded PDB: " << pdb_file << std::endl;

   // Select all atoms and compute centre of mass
   int sel_hnd = mol->NewSelection();
   mol->SelectAtoms(sel_hnd, 0, "*",
                    mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "*", "*", "*", "*");

   mmdb::PPAtom atoms = nullptr;
   int n_atoms = 0;
   mol->GetSelIndex(sel_hnd, atoms, n_atoms);
   std::cout << "Selected " << n_atoms << " atoms\n";

   double sum_x = 0, sum_y = 0, sum_z = 0;
   double min_x = 1e9, max_x = -1e9;
   double min_y = 1e9, max_y = -1e9;
   double min_z = 1e9, max_z = -1e9;

   for (int i = 0; i < n_atoms; ++i) {
      double x = atoms[i]->x;
      double y = atoms[i]->y;
      double z = atoms[i]->z;
      sum_x += x;
      sum_y += y;
      sum_z += z;
      if (x < min_x) min_x = x;
      if (x > max_x) max_x = x;
      if (y < min_y) min_y = y;
      if (y > max_y) max_y = y;
      if (z < min_z) min_z = z;
      if (z > max_z) max_z = z;
   }

   clipper::Coord_orth model_centre(sum_x / n_atoms, sum_y / n_atoms, sum_z / n_atoms);
   std::cout << "Model centre of mass: " << model_centre.x() << ", "
             << model_centre.y() << ", " << model_centre.z() << std::endl;

   mol->DeleteSelection(sel_hnd);

   // Create NXmap for calculated density
   float border = 5.0f;
   float map_min_x = min_x - border;
   float map_min_y = min_y - border;
   float map_min_z = min_z - border;
   float map_max_x = max_x + border;
   float map_max_y = max_y + border;
   float map_max_z = max_z + border;

   int nx = static_cast<int>((map_max_x - map_min_x) / grid_spacing) + 1;
   int ny = static_cast<int>((map_max_y - map_min_y) / grid_spacing) + 1;
   int nz = static_cast<int>((map_max_z - map_min_z) / grid_spacing) + 1;

   std::cout << "Creating NXmap: " << nx << " x " << ny << " x " << nz << " points\n";

   clipper::Grid grid(nx, ny, nz);
   clipper::RTop_orth rtop(clipper::Mat33<>(grid_spacing, 0, 0,
                                            0, grid_spacing, 0,
                                            0, 0, grid_spacing),
                           clipper::Vec3<>(-map_min_x, -map_min_y, -map_min_z));

   clipper::NXmap<float> nxmap(grid, rtop);

   // Build atom list for EDcalc
   sel_hnd = mol->NewSelection();
   mol->SelectAtoms(sel_hnd, 0, "*",
                    mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "*", "*", "*", "*");
   mol->GetSelIndex(sel_hnd, atoms, n_atoms);

   std::vector<clipper::Atom> atom_list;
   for (int i = 0; i < n_atoms; ++i) {
      std::string ele(atoms[i]->element);
      size_t start = ele.find_first_not_of(" ");
      size_t end = ele.find_last_not_of(" ");
      if (start != std::string::npos)
         ele = ele.substr(start, end - start + 1);
      else
         ele = "C";

      clipper::Atom cat;
      cat.set_element(ele);
      cat.set_coord_orth(clipper::Coord_orth(atoms[i]->x, atoms[i]->y, atoms[i]->z));
      float u_iso = (atoms[i]->tempFactor > 0) ? atoms[i]->tempFactor * 0.0125f : 0.5f;
      cat.set_u_iso(u_iso);
      cat.set_occupancy(1.0);
      atom_list.push_back(cat);
   }
   mol->DeleteSelection(sel_hnd);

   std::cout << "Computing calculated density...\n";
   clipper::EDcalc_iso<float> edc(3.0);
   clipper::Atom_list al(atom_list);
   edc(nxmap, al);

   // Compute Zernike descriptor from calculated density (probe)
   std::cout << "Computing Zernike descriptor from calculated density...\n";
   coot::ZernikeDescriptor zd_probe(nxmap, model_centre, radius, n_grid, order);
   zd_probe.print_diagnostics();

   // ============ Load experimental map ============
   std::cout << "\nLoading experimental map: " << exp_map_file << std::endl;
   clipper::CCP4MAPfile mapfile;
   clipper::Xmap<float> xmap;
   mapfile.open_read(exp_map_file);
   mapfile.import_xmap(xmap);
   mapfile.close_read();

   // Compute map sigma for fitting
   float map_sum = 0, map_sum_sq = 0;
   int map_count = 0;
   for (clipper::Xmap<float>::Map_reference_index ix = xmap.first(); !ix.last(); ix.next()) {
      float v = xmap[ix];
      map_sum += v;
      map_sum_sq += v * v;
      map_count++;
   }
   float map_mean = map_sum / map_count;
   float map_sigma = std::sqrt(map_sum_sq / map_count - map_mean * map_mean);
   std::cout << "Map sigma: " << map_sigma << std::endl;

   // Scan centre
   clipper::Coord_orth scan_centre;
   if (use_custom_centre) {
      scan_centre = clipper::Coord_orth(scan_cx, scan_cy, scan_cz);
   } else {
      clipper::Cell cell = xmap.cell();
      scan_centre = clipper::Coord_orth(cell.a()/2, cell.b()/2, cell.c()/2);
   }
   std::cout << "Scan centre: " << scan_centre.x() << ", "
             << scan_centre.y() << ", " << scan_centre.z() << std::endl;

   // ============ Grid scan ============
   std::cout << "\n=== Grid Scan ===" << std::endl;
   int n_steps = static_cast<int>(2 * scan_range / scan_step) + 1;
   int total_points = n_steps * n_steps * n_steps;
   std::cout << "Scan grid: " << n_steps << "^3 = " << total_points << " points\n";

   std::vector<ScanResult> results;
   results.reserve(total_points);

   auto t_scan_start = std::chrono::high_resolution_clock::now();
   int count = 0;

   for (float dx = -scan_range; dx <= scan_range + 0.001; dx += scan_step) {
      for (float dy = -scan_range; dy <= scan_range + 0.001; dy += scan_step) {
         for (float dz = -scan_range; dz <= scan_range + 0.001; dz += scan_step) {
            double x = scan_centre.x() + dx;
            double y = scan_centre.y() + dy;
            double z = scan_centre.z() + dz;
            clipper::Coord_orth pt(x, y, z);

            coot::ZernikeDescriptor zd_target(xmap, pt, radius, n_grid, order);
            float l2 = coot::ZernikeDescriptor::distance(zd_probe, zd_target);
            float cos_sim = coot::ZernikeDescriptor::cosine_similarity(zd_probe, zd_target);
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);

            ScanResult r;
            r.x = x; r.y = y; r.z = z;
            r.l2_norm = l2;
            r.cosine_sim = cos_sim;
            r.dx = dx; r.dy = dy; r.dz = dz;
            r.dist_from_centre = dist;
            results.push_back(r);

            count++;
            if (count % 100 == 0) {
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

   auto t_scan_end = std::chrono::high_resolution_clock::now();
   double scan_time = std::chrono::duration<double>(t_scan_end - t_scan_start).count();
   std::cerr << "\nScan completed: " << count << " points in " << scan_time << " s\n";

   // ============ Sort and get top peaks ============
   std::sort(results.begin(), results.end(), compare_by_l2);

   int n_to_refine = std::min(n_peaks, static_cast<int>(results.size()));
   std::cout << "\n=== Top " << n_to_refine << " Peaks ===" << std::endl;
   std::cout << "Rank  Position             L2_norm    Cosine_sim  Dist_from_centre\n";
   for (int i = 0; i < n_to_refine; ++i) {
      const ScanResult &r = results[i];
      printf("%3d   (%7.1f,%7.1f,%7.1f)  %10.6f  %10.6f  %8.2f\n",
             i+1, r.x, r.y, r.z, r.l2_norm, r.cosine_sim, r.dist_from_centre);
   }

   // ============ Local fitting at each peak ============
   std::cout << "\n=== Local Fitting ===" << std::endl;

   // Create minimol from the model
   coot::minimol::molecule mol_orig(mol);

   // Output stream
   std::ostream *out = &std::cout;
   std::ofstream file_out;
   if (!output_file.empty()) {
      file_out.open(output_file);
      if (file_out.is_open()) {
         out = &file_out;
      }
   }

   *out << "# Zernike Domain Search Results\n";
   *out << "# PDB: " << pdb_file << "\n";
   *out << "# Map: " << exp_map_file << "\n";
   *out << "# Order: " << order << ", Radius: " << radius << "\n";
   *out << "#\n";
   *out << "# Rank  Peak_X  Peak_Y  Peak_Z  L2_norm  Cosine_sim  Fit_score  Dist_from_model_centre\n";

   for (int i = 0; i < n_to_refine; ++i) {
      const ScanResult &r = results[i];

      // Calculate translation vector (from model centre to peak position)
      clipper::Coord_orth translation(r.x - model_centre.x(),
                                       r.y - model_centre.y(),
                                       r.z - model_centre.z());

      // Make a copy of the molecule and translate it
      coot::minimol::molecule mol_copy = mol_orig;
      mol_copy.translate(translation);

      // Rigid body fit using coot::ligand
      coot::ligand lig;
      lig.import_map_from(xmap, map_sigma);
      lig.install_ligand(mol_copy);
      lig.find_centre_by_ligand(0);
      lig.set_map_atom_mask_radius(0.5);
      lig.set_dont_write_solutions();
      lig.set_dont_test_rotations();
      lig.set_acceptable_fit_fraction(0.1);
      lig.fit_ligands_to_clusters(1);

      // Get the fitted solution and score
      coot::minimol::molecule fitted_mol = lig.get_solution(0, 0);
      coot::ligand_score_card score_card = lig.get_solution_score(0, 0);
      float fit_score = score_card.get_score();

      // Calculate distance from original model centre
      float dist_from_model = std::sqrt(
         (r.x - model_centre.x()) * (r.x - model_centre.x()) +
         (r.y - model_centre.y()) * (r.y - model_centre.y()) +
         (r.z - model_centre.z()) * (r.z - model_centre.z()));

      *out << i+1 << "  " << r.x << "  " << r.y << "  " << r.z << "  "
           << r.l2_norm << "  " << r.cosine_sim << "  "
           << fit_score << "  " << dist_from_model << "\n";

      std::cout << "Peak " << (i+1) << ": fit_score = " << fit_score
                << " (dist from model: " << dist_from_model << " A)\n";
   }

   if (file_out.is_open()) {
      file_out.close();
      std::cout << "\nResults written to: " << output_file << std::endl;
   }

   delete mol;

   return 0;
}
