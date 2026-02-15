/*
 * coot-utils/test-zernike.cc
 *
 * Copyright 2025 by Medical Research Council
 * Author: Paul Emsley
 *
 * Test program for 3D Zernike moment computation
 */

#include <iostream>
#include <cstdlib>
#include <clipper/ccp4/ccp4_map_io.h>
#include <mmdb2/mmdb_manager.h>
#include "zernike.hh"

void test_invariant_count() {
   std::cout << "=== Testing invariant count ===" << std::endl;

   // For order N, count should be sum_{n=0}^{N} floor((n+2)/2)
   // order 0: 1 pair (0,0)
   // order 1: 1 pair (1,1)
   // order 2: 2 pairs (2,2), (2,0)
   // order 3: 2 pairs (3,3), (3,1)
   // order 4: 3 pairs (4,4), (4,2), (4,0)
   // etc.

   for (int order = 0; order <= 5; ++order) {
      size_t count = coot::zernike_invariant_count(order);
      std::cout << "Order " << order << ": " << count << " invariants" << std::endl;
   }

   // Order 20 should have: sum_{n=0}^{20} floor((n+2)/2)
   // = 1+1+2+2+3+3+4+4+5+5+6+6+7+7+8+8+9+9+10+10+11 = 121
   size_t count20 = coot::zernike_invariant_count(20);
   std::cout << "Order 20: " << count20 << " invariants (expected 121)" << std::endl;
}

void test_with_map(const std::string &map_file) {
   std::cout << "\n=== Testing with map: " << map_file << " ===" << std::endl;

   try {
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file);
      file.import_xmap(xmap);
      file.close_read();

      // Get cell parameters
      const clipper::Cell &cell = xmap.cell();
      std::cout << "Cell: " << cell.a() << " " << cell.b() << " " << cell.c() << std::endl;

      // Find a high-density point to use as centre
      clipper::Coord_orth centre(cell.a()/2, cell.b()/2, cell.c()/2);

      // Compute Zernike descriptor
      float radius = 15.0;  // 15 Angstrom sphere
      int n_grid = 32;      // 32^3 grid (faster for testing)
      int order = 10;       // order 10 for testing

      std::cout << "Computing Zernike descriptor..." << std::endl;
      std::cout << "  Centre: " << centre.x() << ", " << centre.y() << ", " << centre.z() << std::endl;
      std::cout << "  Radius: " << radius << " A" << std::endl;
      std::cout << "  Grid: " << n_grid << "^3" << std::endl;
      std::cout << "  Order: " << order << std::endl;

      coot::ZernikeDescriptor zd(xmap, centre, radius, n_grid, order);

      // Print diagnostics
      std::cout << "\n--- Descriptor 1 (centre) ---" << std::endl;
      zd.print_diagnostics();
      zd.print_moments(2);  // Print moments up to n=2

      std::vector<float> invariants = zd.get_invariants();
      std::cout << "Computed " << invariants.size() << " invariants" << std::endl;

      // Print first few invariants
      std::cout << "First 10 invariants: ";
      for (size_t i = 0; i < std::min(size_t(10), invariants.size()); ++i) {
         std::cout << invariants[i] << " ";
      }
      std::cout << std::endl;

      // Test rotation invariance by computing at slightly shifted position
      clipper::Coord_orth centre2(5.0 + centre.x(), centre.y(), centre.z());
      // clipper::Coord_orth centre2(107.0, 124.0, 82.0); // centre of the D chain WD40
      coot::ZernikeDescriptor zd2(xmap, centre2, radius, n_grid, order);

      std::cout << "\n--- Descriptor 2 (shifted) ---" << std::endl;
      zd2.print_diagnostics();

      float dist = coot::ZernikeDescriptor::distance(zd, zd2);
      float corr = coot::ZernikeDescriptor::correlation(zd, zd2);
      std::cout << "Distance between centre and shifted centre: " << dist << std::endl;
      std::cout << "Correlation: " << corr << std::endl;

      // Test self-comparison
      float self_dist = coot::ZernikeDescriptor::distance(zd, zd);
      float self_corr = coot::ZernikeDescriptor::correlation(zd, zd);
      std::cout << "Self-distance (should be 0): " << self_dist << std::endl;
      std::cout << "Self-correlation (should be 1): " << self_corr << std::endl;

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }
}

void test_synthetic() {
   std::cout << "\n=== Testing with synthetic data ===" << std::endl;

   // Create a simple test with a synthetic map
   // For now, just test the basic infrastructure

   // Test that default constructor creates empty descriptor
   coot::ZernikeDescriptor empty_zd;
   if (!empty_zd.is_valid()) {
      std::cout << "Empty descriptor is correctly marked as invalid" << std::endl;
   }
}

void sweep_test(const std::string &map_file) {
   std::cout << "\n=== Sweep test for capture radius ===" << std::endl;

   try {
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file);
      file.import_xmap(xmap);
      file.close_read();

      // Reference centre: WD40 domain centre
      clipper::Coord_orth centre(127.01811981201172, 101.62623596191406, 105.45);

      // Parameters
      float radius = 15.0;
      int n_grid = 32;
      int order = 10;

      std::cout << "Reference centre: " << centre.x() << ", " << centre.y() << ", " << centre.z() << std::endl;
      std::cout << "Radius: " << radius << " A, Grid: " << n_grid << "^3, Order: " << order << std::endl;

      // Compute reference descriptor
      std::cout << "Computing reference descriptor..." << std::endl;
      coot::ZernikeDescriptor zd_ref(xmap, centre, radius, n_grid, order);
      zd_ref.print_diagnostics();

      // Output header
      std::cout << "\n# Sweep results: offset_x, offset_y, offset_z, distance" << std::endl;
      std::cout << "# X sweep (Y=0, Z=0)" << std::endl;

      // Sweep in X
      for (int dx = -30; dx <= 30; dx += 2) {
         clipper::Coord_orth pt(centre.x() + dx, centre.y(), centre.z());
         coot::ZernikeDescriptor zd(xmap, pt, radius, n_grid, order);
         float dist = coot::ZernikeDescriptor::distance(zd_ref, zd);
         std::cout << dx << " 0 0 " << dist << std::endl;
      }

      std::cout << "\n# Y sweep (X=0, Z=0)" << std::endl;

      // Sweep in Y
      for (int dy = -30; dy <= 30; dy += 2) {
         clipper::Coord_orth pt(centre.x(), centre.y() + dy, centre.z());
         coot::ZernikeDescriptor zd(xmap, pt, radius, n_grid, order);
         float dist = coot::ZernikeDescriptor::distance(zd_ref, zd);
         std::cout << "0 " << dy << " 0 " << dist << std::endl;
      }

      std::cout << "\n# Z sweep (X=0, Y=0)" << std::endl;

      // Sweep in Z
      for (int dz = -30; dz <= 30; dz += 2) {
         clipper::Coord_orth pt(centre.x(), centre.y(), centre.z() + dz);
         coot::ZernikeDescriptor zd(xmap, pt, radius, n_grid, order);
         float dist = coot::ZernikeDescriptor::distance(zd_ref, zd);
         std::cout << "0 0 " << dz << " " << dist << std::endl;
      }

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }
}

void sweep_3d_test(const std::string &map_file) {

   std::cout << "\n=== 3D Sweep test - no masking ===" << std::endl;

   try {
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file);
      file.import_xmap(xmap);
      file.close_read();

      // Reference centre: WD40 domain centre
      clipper::Coord_orth ref_centre(127.01811981201172, 101.62623596191406, 105.45);

      // Sweep centre: middle of map
      clipper::Coord_orth sweep_centre(112.0, 112.0, 112.0);

      // Parameters
      float radius = 22.0;
      int n_grid   = 32;
      int order    = 16;
      int step     =  5;
      int range    = 40;

      std::cout << "Reference centre (WD40): " << ref_centre.x() << ", "
                << ref_centre.y() << ", " << ref_centre.z() << std::endl;
      std::cout << "Sweep centre: " << sweep_centre.x() << ", "
                << sweep_centre.y() << ", " << sweep_centre.z() << std::endl;
      std::cout << "Radius: " << radius << " A, Grid: " << n_grid
                << "^3, Order: " << order << std::endl;
      std::cout << "Sweep range: +/-" << range << " A, Step: " << step << " A" << std::endl;

      // Compute reference descriptor
      std::cout << "Computing reference descriptor..." << std::endl;
      coot::ZernikeDescriptor zd_ref(xmap, ref_centre, radius, n_grid, order);
      zd_ref.print_diagnostics();

      // Output header
      std::cout << "\n# 3D Sweep results: x, y, z, distance" << std::endl;

      // 3D sweep
      for (int dx = -range; dx <= range; dx += step) {
         for (int dy = -range; dy <= range; dy += step) {
            for (int dz = -range; dz <= range; dz += step) {
               double x = sweep_centre.x() + dx;
               double y = sweep_centre.y() + dy;
               double z = sweep_centre.z() + dz;
               clipper::Coord_orth pt(x, y, z);
               coot::ZernikeDescriptor zd(xmap, pt, radius, n_grid, order);
               float dist = coot::ZernikeDescriptor::distance(zd_ref, zd);
               std::cout << x << " " << y << " " << z << " " << dist << std::endl;
            }
         }
      }

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }
}

void sweep_3d_masked_test(const std::string &map_file, const std::string &pdb_file) {

   std::cout << "\n=== 3D Sweep test with protein mask ===" << std::endl;

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
         return;
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
      std::cout << "Reference centre (from " << n_ca << " CA atoms): "
                << ref_centre.x() << ", " << ref_centre.y() << ", " << ref_centre.z() << std::endl;

      // Sweep centre: middle of map
      clipper::Coord_orth sweep_centre(112.0, 112.0, 112.0);

      // Parameters
      float radius = 23.0;
      int n_grid = 32;
      int order  = 10;
      int step   = 5;
      int range  = 44;
      float mask_radius = 5.0;  // 5A from CA atoms

      std::cout << "Reference centre (WD40): " << ref_centre.x() << ", "
                << ref_centre.y() << ", " << ref_centre.z() << std::endl;
      std::cout << "Sweep centre: " << sweep_centre.x() << ", "
                << sweep_centre.y() << ", " << sweep_centre.z() << std::endl;
      std::cout << "Radius: " << radius << " A, Grid: " << n_grid
                << "^3, Order: " << order << std::endl;
      std::cout << "Sweep range: +/-" << range << " A, Step: " << step << " A" << std::endl;
      std::cout << "Mask radius: " << mask_radius << " A from CA atoms" << std::endl;
      std::cout << "PDB file: " << pdb_file << std::endl;

      // Compute reference descriptor WITH protein mask
      std::cout << "Computing reference descriptor (protein-masked)..." << std::endl;
      coot::ZernikeDescriptor zd_ref(xmap, ref_centre, radius, mol, mask_radius, n_grid, order);
      zd_ref.print_diagnostics();

      // Output header
      std::cout << "\n# 3D Sweep results (masked reference): x, y, z, distance, cosine_similarity" << std::endl;

      // 3D sweep - search points are NOT masked, only reference is masked
      for (int dx = -range; dx <= range; dx += step) {
         for (int dy = -range; dy <= range; dy += step) {
            for (int dz = -range; dz <= range; dz += step) {
               double x = sweep_centre.x() + dx;
               double y = sweep_centre.y() + dy;
               double z = sweep_centre.z() + dz;
               clipper::Coord_orth pt(x, y, z);
               // Search descriptor is NOT masked (we don't know where protein is)
               coot::ZernikeDescriptor zd(xmap, pt, radius, n_grid, order);
               float dist = coot::ZernikeDescriptor::distance(zd_ref, zd);
               float cos_sim = coot::ZernikeDescriptor::cosine_similarity(zd_ref, zd);
               std::cout << x << " " << y << " " << z << " " << dist << " " << cos_sim << std::endl;
            }
         }
      }

      delete mol;

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }
}

void local_scan_test(const std::string &map_file, const std::string &pdb_file, int position = 0) {
   // position: 0 = all positions, 1-4 = specific position
   std::cout << "\n=== Local scan around true WD40 positions ===" << std::endl;
   if (position > 0)
      std::cout << "Scanning position " << position << " only" << std::endl;

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
         return;
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
      std::cout << "Reference centre (from " << n_ca << " CA atoms): "
                << ref_centre.x() << ", " << ref_centre.y() << ", " << ref_centre.z() << std::endl;

      // Parameters
      float radius = 23.0;
      int n_grid = 64;   // was 32, increased for smoother sampling (0.72A step vs 1.44A)
      int order = 24;
      float mask_radius = 5.0;

      std::cout << "Radius: " << radius << " A, Grid: " << n_grid
                << "^3, Order: " << order << std::endl;

      // Compute reference descriptor WITH protein mask
      std::cout << "Computing reference descriptor (protein-masked)..." << std::endl;
      coot::ZernikeDescriptor zd_ref(xmap, ref_centre, radius, mol, mask_radius, n_grid, order);
      zd_ref.print_diagnostics();

      // True WD40 positions
      std::vector<clipper::Coord_orth> true_positions = {
         clipper::Coord_orth(107.08999633789062, 124.80999755859375,  82.05999755859375),   // D chain
         clipper::Coord_orth(107.47451782226562, 124.73106384277344, 126.4388656616211),
         clipper::Coord_orth(128.3690185546875,  101.39816284179688, 107.40684509277344),
         clipper::Coord_orth(117.29839324951172,  79.96355438232422, 146.205322265625)
      };

      std::vector<std::string> position_names = {
         "D_chain", "position_2", "position_3", "position_4"
      };

      // Local scan parameters
      int scan_range = 10;  // +/- 5A = 10A total
      int scan_step = 1;

      // Scan around each true position (or just the specified one)
      size_t start_idx = (position > 0) ? position - 1 : 0;
      size_t end_idx = (position > 0) ? position : true_positions.size();
      for (size_t i = start_idx; i < end_idx; ++i) {
         const auto &true_pos = true_positions[i];
         std::cout << "\n# Scanning around " << position_names[i] << ": "
                   << true_pos.x() << ", " << true_pos.y() << ", " << true_pos.z() << std::endl;
         std::cout << "# x y z L2_norm cosine_similarity" << std::endl;

         for (int dx = -scan_range; dx <= scan_range; dx += scan_step) {
            for (int dy = -scan_range; dy <= scan_range; dy += scan_step) {
               for (int dz = -scan_range; dz <= scan_range; dz += scan_step) {
                  double x = true_pos.x() + dx;
                  double y = true_pos.y() + dy;
                  double z = true_pos.z() + dz;
                  clipper::Coord_orth pt(x, y, z);
                  coot::ZernikeDescriptor zd(xmap, pt, radius, n_grid, order);
                  float dist = coot::ZernikeDescriptor::distance(zd_ref, zd);
                  float cos_sim = coot::ZernikeDescriptor::cosine_similarity(zd_ref, zd);
                  std::cout << x << " " << y << " " << z << " " << dist << " " << cos_sim << std::endl;
               }
            }
         }
      }

      delete mol;

   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }
}

int main(int argc, char **argv) {

   test_invariant_count();
   test_synthetic();

   if (argc > 1) {
      std::string arg1 = argv[1];
      if (arg1 == "--sweep" && argc > 2) {
         sweep_test(argv[2]);
      } else if (arg1 == "--sweep3d" && argc > 2) {
         sweep_3d_test(argv[2]);
      } else if (arg1 == "--sweep3d-masked" && argc > 3) {
         sweep_3d_masked_test(argv[2], argv[3]);
      } else if (arg1 == "--local-scan" && argc > 3) {
         int position = 0;  // 0 = all positions
         if (argc > 4) {
            position = std::atoi(argv[4]);
            if (position < 1 || position > 4) {
               std::cerr << "Position must be 1-4, got: " << argv[4] << std::endl;
               return 1;
            }
         }
         local_scan_test(argv[2], argv[3], position);
      } else {
         test_with_map(argv[1]);
      }
   } else {
      std::cout << "\nUsage:" << std::endl;
      std::cout << "  test-zernike <map.ccp4>                        - run basic tests" << std::endl;
      std::cout << "  test-zernike --sweep <map.ccp4>                - run 1D sweep tests" << std::endl;
      std::cout << "  test-zernike --sweep3d <map.ccp4>              - run 3D sweep test" << std::endl;
      std::cout << "  test-zernike --sweep3d-masked <map> <pdb>      - 3D sweep with protein mask" << std::endl;
      std::cout << "  test-zernike --local-scan <map> <pdb> [1-4]    - scan around WD40 position(s)" << std::endl;
   }

   return 0;
}
