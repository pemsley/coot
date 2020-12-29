
#ifndef ARC_INFO_TYPE_HH
#define ARC_INFO_TYPE_HH

namespace coot {
// can throw an exception (e.g. null pointers, overlapping atoms)
   //
   class arc_info_type {
   public:
      float delta; // the difference between the starting angle delta
                   // and the end - using the orientation_matrix, means
                   // that the start angle is 0.
      clipper::Coord_orth start_point;
      clipper::Coord_orth start_dir;
      clipper::Coord_orth normal;
      clipper::Mat33<double> orientation_matrix;
      arc_info_type(mmdb::Atom *at_1, mmdb::Atom *at_2, mmdb::Atom *at_3);
   };

}

#endif // ARC_INFO_TYPE_HH
