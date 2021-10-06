

#ifndef SIMPLE_DISTANCE_OBJECT_T_HH
#define SIMPLE_DISTANCE_OBJECT_T_HH

#include <clipper/core/coords.h>

namespace coot {

   // this is a copy of what's in graphics-info.h and ideally graphics-info.h should
   // include this header.
   //
   // 20211006-PE it does now

   class simple_distance_object_t {
   public:
      clipper::Coord_orth start_pos;
      clipper::Coord_orth end_pos;
      int imol_start;
      int imol_end;
      simple_distance_object_t(int imol1,
                               const clipper::Coord_orth &start,
                               int imol2,
                               const clipper::Coord_orth &end) {
         start_pos = start;
         end_pos = end;
         imol_start = imol1;
         imol_end = imol2;
      }
      double length() const {
         double len_sqrd = (end_pos-start_pos).lengthsq();
         return std::sqrt(len_sqrd);
      }
      clipper::Coord_orth uv() const {
         clipper::Coord_orth delta = end_pos - start_pos;
         double l = length();
         if (l > 0.0) {
            clipper::Coord_orth delta_uv = (1.0/l) * delta;
            return delta_uv;
         } else {
            // make one up
            return clipper::Coord_orth(0,0,1);
         }
      }
      friend std::ostream& operator<<(std::ostream &s, simple_distance_object_t o);
   };
   std::ostream& operator<<(std::ostream &s, simple_distance_object_t o);

}

#endif // SIMPLE_DISTANCE_OBJECT_T_HH

