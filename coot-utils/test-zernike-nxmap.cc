/*
 * coot-utils/test-zernike-nxmap.cc
 *
 * Test Zernike descriptors with calculated density (NXmap)
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/contrib/edcalc.h>
#include <mmdb2/mmdb_manager.h>
#include "zernike.hh"
#include "coot-map-heavy.hh"
#include "coot-map-utils.hh"

// Convert NXmap to Xmap in P1 cell for blurring
clipper::Xmap<float> nxmap_to_xmap(const clipper::NXmap<float> &nxmap,
                                    float grid_spacing,
                                    float min_x, float min_y, float min_z,
                                    float max_x, float max_y, float max_z) {

   // Create P1 cell covering the NXmap region with some padding
   float padding = 5.0f;
   float cell_a = (max_x - min_x) + 2 * padding;
   float cell_b = (max_y - min_y) + 2 * padding;
   float cell_c = (max_z - min_z) + 2 * padding;

   clipper::Cell cell(clipper::Cell_descr(cell_a, cell_b, cell_c, 90, 90, 90));
   clipper::Spacegroup spgr(clipper::Spacegroup::P1);

   // Grid sampling to match NXmap spacing
   int nu = static_cast<int>(cell_a / grid_spacing) + 1;
   int nv = static_cast<int>(cell_b / grid_spacing) + 1;
   int nw = static_cast<int>(cell_c / grid_spacing) + 1;

   clipper::Grid_sampling gs(nu, nv, nw);
   clipper::Xmap<float> xmap(spgr, cell, gs);

   // Origin offset: Xmap fractional (0,0,0) corresponds to orthogonal (0,0,0)
   // We want to map NXmap coordinates to Xmap coordinates
   // NXmap origin in orthogonal coords is (min_x, min_y, min_z)
   // Xmap origin in orthogonal coords is (0, 0, 0)
   // So we offset by (-min_x + padding, -min_y + padding, -min_z + padding)

   float offset_x = -min_x + padding;
   float offset_y = -min_y + padding;
   float offset_z = -min_z + padding;

   // Initialize xmap to zero
   for (clipper::Xmap<float>::Map_reference_index ix = xmap.first(); !ix.last(); ix.next()) {
      xmap[ix] = 0.0f;
   }

   // Copy data from NXmap to Xmap
   for (clipper::NXmap<float>::Map_reference_index ix = nxmap.first(); !ix.last(); ix.next()) {
      clipper::Coord_grid nxmap_cg = ix.coord();
      clipper::Coord_map nxmap_cm(nxmap_cg.u(), nxmap_cg.v(), nxmap_cg.w());
      clipper::Coord_orth co = nxmap.coord_orth(nxmap_cm);
      // Transform to Xmap coordinates
      clipper::Coord_orth co_xmap(co.x() + offset_x, co.y() + offset_y, co.z() + offset_z);
      clipper::Coord_frac cf = co_xmap.coord_frac(cell);
      clipper::Coord_map cm = cf.coord_map(gs);
      clipper::Coord_grid cg = cm.coord_grid();

      if (cg.u() >= 0 && cg.u() < nu &&
          cg.v() >= 0 && cg.v() < nv &&
          cg.w() >= 0 && cg.w() < nw) {
         xmap.set_data(cg, nxmap[ix]);
      }
   }

   return xmap;
}

void print_usage(const char *prog) {
   std::cout << "Usage: " << prog << " <pdb-file> [options]\n"
             << "\n"
             << "Computes Zernike descriptor from calculated density\n"
             << "\n"
             << "Options:\n"
             << "  --n_grid N     Grid points per dimension (default: 32)\n"
             << "  --order N      Zernike order (default: 20)\n"
             << "  --radius R     Sphere radius in Angstroms (default: 23)\n"
             << "  --grid_spacing S  Map grid spacing in Angstroms (default: 1.0)\n"
             << "  --exp_map FILE   Experimental map for comparison/scanning\n"
             << "  --scan           Enable grid search mode\n"
             << "  --scan_centre X,Y,Z  Centre for scanning (default: map centre)\n"
             << "  --range R      Scan range +/- in Angstroms (default: 45)\n"
             << "  --step S       Scan step size in Angstroms (default: 3)\n"
             << "  --output FILE  Output file for scan results\n"
             << "  --mainchain      Use mainchain-weighted density (MC+CB=1.0, CG=0.5, others ignored)\n"
             << "  --blur B         Blur both maps by B (A^2), e.g. 80\n"
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
   bool do_scan = false;
   bool use_custom_centre = false;
   double scan_cx = 0, scan_cy = 0, scan_cz = 0;
   float scan_range = 45.0f;
   float scan_step = 3.0f;
   std::string output_file;
   bool mainchain_weighted = false;
   float blur_b_factor = 0.0f;

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
      } else if (arg == "--scan") {
         do_scan = true;
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
      } else if (arg == "--output" && i + 1 < argc) {
         output_file = argv[++i];
      } else if (arg == "--mainchain") {
         mainchain_weighted = true;
      } else if (arg == "--blur" && i + 1 < argc) {
         blur_b_factor = std::atof(argv[++i]);
      }
   }

   std::cout << "Parameters:\n"
             << "  n_grid: " << n_grid << "\n"
             << "  order: " << order << "\n"
             << "  radius: " << radius << " A\n"
             << "  grid_spacing: " << grid_spacing << " A\n"
             << "  mainchain_weighted: " << (mainchain_weighted ? "yes" : "no") << "\n"
             << "  blur_b_factor: " << blur_b_factor << " A^2\n"
             << std::endl;

   // Load PDB
   mmdb::Manager *mol = new mmdb::Manager();
   mmdb::ERROR_CODE rc = mol->ReadCoorFile(pdb_file.c_str());
   if (rc != mmdb::Error_NoError) {
      std::cerr << "Error reading PDB file: " << pdb_file << std::endl;
      delete mol;
      return 1;
   }

   std::cout << "Loaded PDB: " << pdb_file << std::endl;

   // Select all atoms
   int sel_hnd = mol->NewSelection();
   mol->SelectAtoms(sel_hnd, 0, "*",
                    mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "*", "*", "*", "*");

   mmdb::PPAtom atoms = nullptr;
   int n_atoms = 0;
   mol->GetSelIndex(sel_hnd, atoms, n_atoms);
   std::cout << "Selected " << n_atoms << " atoms\n";

   // Calculate centre of mass
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

   clipper::Coord_orth centre(sum_x / n_atoms, sum_y / n_atoms, sum_z / n_atoms);
   std::cout << "Centre of mass: " << centre.x() << ", " << centre.y() << ", " << centre.z() << std::endl;
   std::cout << "Bounding box: [" << min_x << ", " << max_x << "] x ["
             << min_y << ", " << max_y << "] x ["
             << min_z << ", " << max_z << "]\n";

   // Create NXmap for calculated density
   // The map needs to cover the molecule plus some border
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

   std::cout << "Creating NXmap: " << nx << " x " << ny << " x " << nz
             << " = " << (nx * ny * nz) << " points\n";

   // Create grid and transformation
   // NXmap coord_orth computes: orth = rot * grid - trn
   // So to have grid (0,0,0) -> (map_min_x, map_min_y, map_min_z)
   // we need trn = -map_min
   clipper::Grid grid(nx, ny, nz);
   clipper::RTop_orth rtop(clipper::Mat33<>(grid_spacing, 0, 0,
                                            0, grid_spacing, 0,
                                            0, 0, grid_spacing),
                           clipper::Vec3<>(-map_min_x, -map_min_y, -map_min_z));

   std::cout << "RTop translation: " << -map_min_x << ", " << -map_min_y << ", " << -map_min_z << std::endl;
   std::cout << "RTop scale: " << grid_spacing << std::endl;

   clipper::NXmap<float> nxmap_ref(grid, rtop);

   std::cout << "Computing calculated density...\n";

   // Debug: check NXmap ref properties
   std::cout << "NXmap ref grid: " << nxmap_ref.grid().nu() << " x "
             << nxmap_ref.grid().nv() << " x " << nxmap_ref.grid().nw() << std::endl;

   // Use EDcalc_iso directly with trimmed element strings
   clipper::NXmap<float> nxmap(nxmap_ref.grid(), nxmap_ref.operator_orth_grid());

   clipper::EDcalc_iso<float> edc(3.0);
   std::vector<clipper::Atom> atom_list;

   int n_included = 0;
   int n_skipped = 0;

   for (int i = 0; i < n_atoms; ++i) {
      std::string ele(atoms[i]->element);
      // Trim whitespace from element
      size_t start = ele.find_first_not_of(" ");
      size_t end = ele.find_last_not_of(" ");
      if (start != std::string::npos)
         ele = ele.substr(start, end - start + 1);
      else
         ele = "C";

      // Get atom name and trim whitespace
      std::string atom_name(atoms[i]->name);
      start = atom_name.find_first_not_of(" ");
      end = atom_name.find_last_not_of(" ");
      if (start != std::string::npos)
         atom_name = atom_name.substr(start, end - start + 1);

      float occupancy = 1.0f;

      if (mainchain_weighted) {
         // Mainchain atoms (N, CA, C, O) and CB: weight 1.0
         // CG atoms (CG, CG1, CG2): weight 0.5
         // All others: skip
         if (atom_name == "N" || atom_name == "CA" || atom_name == "C" || atom_name == "O" ||
             atom_name == "CB") {
            occupancy = 1.0f;
         } else if (atom_name == "CG" || atom_name == "CG1" || atom_name == "CG2") {
            occupancy = 0.5f;
         } else {
            n_skipped++;
            continue;  // Skip this atom
         }
      }

      clipper::Atom cat;
      cat.set_element(ele);
      cat.set_coord_orth(clipper::Coord_orth(atoms[i]->x, atoms[i]->y, atoms[i]->z));
      float u_iso = (atoms[i]->tempFactor > 0) ? atoms[i]->tempFactor * 0.0125f : 0.5f;
      cat.set_u_iso(u_iso);
      cat.set_occupancy(occupancy);
      atom_list.push_back(cat);
      n_included++;
   }

   if (mainchain_weighted) {
      std::cout << "Mainchain weighting: included " << n_included << " atoms, skipped " << n_skipped << std::endl;
   }

   // Debug: print first few processed atoms
   for (int i = 0; i < 5 && i < static_cast<int>(atom_list.size()); ++i) {
      std::cout << "Clipper Atom " << i << ": ele='" << atom_list[i].element() << "'"
                << " pos=" << atom_list[i].coord_orth().format()
                << " u_iso=" << atom_list[i].u_iso()
                << " occ=" << atom_list[i].occupancy() << std::endl;
   }

   // Debug: check NXmap coordinate transformation
   clipper::RTop<> actual_rtop = nxmap.operator_orth_grid();
   std::cout << "Actual NXmap RTop translation: " << actual_rtop.trn()[0] << ", "
             << actual_rtop.trn()[1] << ", " << actual_rtop.trn()[2] << std::endl;

   // Test: what does grid (0,0,0) map to?
   clipper::Coord_map cm_origin(0, 0, 0);
   clipper::Coord_orth origin_orth = nxmap.coord_orth(cm_origin);
   std::cout << "Grid (0,0,0) -> orth: " << origin_orth.format() << std::endl;

   clipper::Coord_orth test_pt(centre.x(), centre.y(), centre.z());
   clipper::Coord_map test_cm = nxmap.coord_map(test_pt);
   std::cout << "Centre in orth: " << test_pt.format() << std::endl;
   std::cout << "Centre in map coords: " << test_cm.format() << std::endl;
   std::cout << "Grid range: [0," << nx-1 << "] x [0," << ny-1 << "] x [0," << nz-1 << "]" << std::endl;

   clipper::Atom_list al(atom_list);
   std::cout << "Calling EDcalc_iso with " << atom_list.size() << " atoms...\n";
   edc(nxmap, al);
   std::cout << "EDcalc_iso done\n";

   // Debug: check a few grid points
   std::cout << "Checking grid values around centre...\n";
   int cx = static_cast<int>(test_cm.u());
   int cy = static_cast<int>(test_cm.v());
   int cz = static_cast<int>(test_cm.w());
   for (int dx = -2; dx <= 2; ++dx) {
      int gx = cx + dx;
      int gy = cy;
      int gz = cz;
      if (gx >= 0 && gx < nx && gy >= 0 && gy < ny && gz >= 0 && gz < nz) {
         clipper::Coord_grid cg(gx, gy, gz);
         float val = nxmap.get_data(cg);
         std::cout << "  Grid[" << gx << "," << gy << "," << gz << "] = " << val << std::endl;
      }
   }

   mol->DeleteSelection(sel_hnd);

   // Check density statistics
   float density_min = 1e9, density_max = -1e9, density_sum = 0;
   int n_points = 0;
   for (clipper::NXmap<float>::Map_reference_index ix = nxmap.first(); !ix.last(); ix.next()) {
      float d = nxmap[ix];
      if (d < density_min) density_min = d;
      if (d > density_max) density_max = d;
      density_sum += d;
      n_points++;
   }
   std::cout << "Density stats: min=" << density_min << " max=" << density_max
             << " mean=" << (density_sum / n_points) << std::endl;

   // Compute Zernike descriptor
   // If blurring is requested, convert NXmap to Xmap and blur it
   coot::ZernikeDescriptor zd;
   if (blur_b_factor > 0.001f) {
      std::cout << "\nConverting calculated density to Xmap for blurring...\n";
      clipper::Xmap<float> calc_xmap = nxmap_to_xmap(nxmap, grid_spacing,
                                                      map_min_x, map_min_y, map_min_z,
                                                      map_max_x, map_max_y, map_max_z);
      std::cout << "Blurring calculated density by " << blur_b_factor << " A^2...\n";
      calc_xmap = coot::util::sharpen_blur_map(calc_xmap, blur_b_factor);
      std::cout << "Computing Zernike descriptor from blurred Xmap...\n";
      zd = coot::ZernikeDescriptor(calc_xmap, centre, radius, n_grid, order);
   } else {
      std::cout << "\nComputing Zernike descriptor from NXmap...\n";
      zd = coot::ZernikeDescriptor(nxmap, centre, radius, n_grid, order);
   }

   zd.print_diagnostics();

   // Print first few invariants
   std::vector<float> inv = zd.get_invariants();
   std::cout << "\nFirst 20 invariants:\n";
   for (size_t i = 0; i < std::min(inv.size(), size_t(20)); ++i) {
      std::cout << "  F[" << i << "] = " << inv[i] << "\n";
   }

   std::cout << "\nTotal invariants: " << inv.size() << std::endl;

   // Compare with experimental map if provided
   if (!exp_map_file.empty()) {
      std::cout << "\n=== Loading experimental map ===\n";
      std::cout << "Loading: " << exp_map_file << std::endl;

      clipper::CCP4MAPfile mapfile;
      clipper::Xmap<float> xmap;
      mapfile.open_read(exp_map_file);
      mapfile.import_xmap(xmap);
      mapfile.close_read();

      // Apply blurring if requested
      if (blur_b_factor > 0.001f) {
         std::cout << "Blurring map by " << blur_b_factor << " A^2...\n";
         xmap = coot::util::sharpen_blur_map(xmap, blur_b_factor);
      }

      if (do_scan) {
         // Grid search mode: use calculated density descriptor as probe
         std::cout << "\n=== Grid Search Mode ===\n";
         std::cout << "Probe: Zernike descriptor from calculated density (NXmap)\n";
         std::cout << "Target: Experimental map\n";

         // Output stream
         std::ostream *out = &std::cout;
         std::ofstream file_out;
         if (!output_file.empty()) {
            file_out.open(output_file);
            if (file_out.is_open()) {
               out = &file_out;
               std::cout << "Writing output to: " << output_file << std::endl;
            }
         }

         // Scan centre
         clipper::Coord_orth scan_centre;
         if (use_custom_centre) {
            scan_centre = clipper::Coord_orth(scan_cx, scan_cy, scan_cz);
         } else {
            // Use map centre as default scan centre
            clipper::Cell cell = xmap.cell();
            scan_centre = clipper::Coord_orth(cell.a()/2, cell.b()/2, cell.c()/2);
         }

         *out << "# Zernike Grid Search - Calculated density probe\n";
         *out << "# Parameters:\n";
         *out << "#   n_grid: " << n_grid << "\n";
         *out << "#   order: " << order << "\n";
         *out << "#   radius: " << radius << " A\n";
         *out << "#   scan_range: +/- " << scan_range << " A\n";
         *out << "#   scan_step: " << scan_step << " A\n";
         *out << "#   PDB file: " << pdb_file << "\n";
         *out << "#   Exp map: " << exp_map_file << "\n";
         *out << "# Probe centre (model COM): " << centre.x() << ", " << centre.y() << ", " << centre.z() << "\n";
         *out << "# Scan centre: " << scan_centre.x() << ", " << scan_centre.y() << ", " << scan_centre.z() << "\n";

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

                  // Compute descriptor from experimental map at this position
                  coot::ZernikeDescriptor zd_exp(xmap, pt, radius, n_grid, order);

                  // Compare with calculated density probe
                  float dist = coot::ZernikeDescriptor::distance(zd, zd_exp);
                  float cos_sim = coot::ZernikeDescriptor::cosine_similarity(zd, zd_exp);

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

      } else {
         // Single comparison mode
         std::cout << "Computing Zernike descriptor from experimental map at same centre...\n";
         coot::ZernikeDescriptor zd_exp(xmap, centre, radius, n_grid, order);
         zd_exp.print_diagnostics();

         std::vector<float> inv_exp = zd_exp.get_invariants();

         // Compare descriptors
         float l2_dist = coot::ZernikeDescriptor::distance(zd, zd_exp);
         float cos_sim = coot::ZernikeDescriptor::cosine_similarity(zd, zd_exp);
         float corr = coot::ZernikeDescriptor::correlation(zd, zd_exp);

         std::cout << "\n=== Comparison Results ===\n";
         std::cout << "L2 distance:        " << l2_dist << std::endl;
         std::cout << "Cosine similarity:  " << cos_sim << std::endl;
         std::cout << "Correlation:        " << corr << std::endl;

         // Print side-by-side invariants
         std::cout << "\nFirst 20 invariants (Calc vs Exp):\n";
         for (size_t i = 0; i < std::min(inv.size(), size_t(20)); ++i) {
            std::cout << "  F[" << i << "]: " << inv[i] << " vs " << inv_exp[i]
                      << " (ratio: " << (inv_exp[i] > 1e-10 ? inv[i]/inv_exp[i] : 0) << ")\n";
         }
      }
   }

   delete mol;

   return 0;
}
