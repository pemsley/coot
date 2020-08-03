

#ifndef SIMPLE_DISTANCE_OBJECT_T_HH
#define SIMPLE_DISTANCE_OBJECT_T_HH

namespace coot {

   // this is a copy of what's in graphics-info.h and ideally graphics-info.h should
   // include this header

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
     friend std::ostream& operator<<(std::ostream &s, simple_distance_object_t o);
   };
   std::ostream& operator<<(std::ostream &s, simple_distance_object_t o);

}

#endif // SIMPLE_DISTANCE_OBJECT_T_HH

