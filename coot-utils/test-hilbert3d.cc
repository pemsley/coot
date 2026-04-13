#include "hilbert3d.hh"
#include <iostream>
#include <set>
#include <cstdlib>

int main() {

   for (int order=1; order<=5; order++) {
      int n = 1 << order;
      int total = n * n * n;
      std::set<uint64_t> seen;

      // Round-trip: xyz -> d -> xyz
      for (int x=0; x<n; x++) {
         for (int y=0; y<n; y++) {
            for (int z=0; z<n; z++) {
               uint64_t d = coot::xyz_to_hilbert_3d(x, y, z, order);
               if (d < 0 || static_cast<int>(d) >= total) {
                  std::cout << "FAIL: d=" << d << " out of range, order=" << order << std::endl;
                  return 1;
               }
               if (seen.count(d)) {
                  std::cout << "FAIL: duplicate d=" << d << " at (" << x << "," << y << "," << z
                            << "), order=" << order << std::endl;
                  return 1;
               }
               seen.insert(d);
               auto [x2, y2, z2] = coot::hilbert_to_xyz_3d(d, order);
               if (x2 != x || y2 != y || z2 != z) {
                  std::cout << "FAIL: round-trip (" << x << "," << y << "," << z << ") -> "
                            << d << " -> (" << x2 << "," << y2 << "," << z2 << ")" << std::endl;
                  return 1;
               }
            }
         }
      }

      // Round-trip: d -> xyz -> d
      for (int d=0; d<total; d++) {
         auto [x, y, z] = coot::hilbert_to_xyz_3d(d, order);
         uint64_t d2 = coot::xyz_to_hilbert_3d(x, y, z, order);
         if (static_cast<int>(d2) != d) {
            std::cout << "FAIL: inverse round-trip d=" << d << " -> (" << x << "," << y << "," << z
                      << ") -> " << d2 << std::endl;
            return 1;
         }
      }

      // Adjacency: consecutive indices must be Manhattan distance 1
      auto [px, py, pz] = coot::hilbert_to_xyz_3d(0, order);
      for (int d=1; d<total; d++) {
         auto [cx, cy, cz] = coot::hilbert_to_xyz_3d(d, order);
         int dist = std::abs(cx-px) + std::abs(cy-py) + std::abs(cz-pz);
         if (dist != 1) {
            std::cout << "FAIL: order=" << order << " d=" << d << " dist=" << dist
                      << " from (" << px << "," << py << "," << pz << ") to ("
                      << cx << "," << cy << "," << cz << ")" << std::endl;
            return 1;
         }
         px = cx; py = cy; pz = cz;
      }

      std::cout << "order " << order << ": all " << total << " points OK" << std::endl;
   }

   std::cout << "All tests passed." << std::endl;
   return 0;
}
