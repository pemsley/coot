#ifndef HILBERT3D_HH
#define HILBERT3D_HH

#include <array>
#include <cstdint>
#include <tuple>

namespace coot {

   // 3D Hilbert curve conversions for a 2^order cube.
   // Based on Skilling, "Programming the Hilbert curve",
   // AIP Conf. Proc. 707, 381 (2004).

   inline uint64_t xyz_to_hilbert_3d(int x, int y, int z, int order) {

      std::array<int, 3> X = {x, y, z};
      int M = 1 << (order - 1);

      // Inverse undo
      int Q = M;
      while (Q > 1) {
         int P = Q - 1;
         for (int i=0; i<3; i++) {
            if (X[i] & Q) {
               X[0] ^= P;
            } else {
               int t = (X[0] ^ X[i]) & P;
               X[0] ^= t;
               X[i] ^= t;
            }
         }
         Q >>= 1;
      }

      // Gray encode
      for (int i=1; i<3; i++)
         X[i] ^= X[i - 1];
      int t = 0;
      Q = M;
      while (Q > 1) {
         if (X[2] & Q)
            t ^= Q - 1;
         Q >>= 1;
      }
      for (int i=0; i<3; i++)
         X[i] ^= t;

      // Interleave bits into scalar index
      uint64_t d = 0;
      for (int i=order-1; i>=0; i--)
         for (int dim=0; dim<3; dim++)
            d = (d << 1) | ((X[dim] >> i) & 1);

      return d;
   }

   inline std::tuple<int, int, int> hilbert_to_xyz_3d(uint64_t d, int order) {

      // De-interleave scalar index into transposed form
      std::array<int, 3> X = {0, 0, 0};
      for (int i=0; i<order; i++)
         for (int dim=2; dim>=0; dim--) {
            X[dim] |= (static_cast<int>(d) & 1) << i;
            d >>= 1;
         }

      int N = 1 << order;

      // Gray decode
      int t = X[2] >> 1;
      for (int i=2; i>0; i--)
         X[i] ^= X[i - 1];
      X[0] ^= t;

      // Undo excess work
      int Q = 2;
      while (Q != N) {
         int P = Q - 1;
         for (int i=2; i>=0; i--) {
            if (X[i] & Q) {
               X[0] ^= P;
            } else {
               t = (X[0] ^ X[i]) & P;
               X[0] ^= t;
               X[i] ^= t;
            }
         }
         Q <<= 1;
      }

      return {X[0], X[1], X[2]};
   }

} // namespace coot

#endif // HILBERT3D_HH
