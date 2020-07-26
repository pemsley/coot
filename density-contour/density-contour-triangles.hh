

#ifndef DENSITY_CONTOUR_TRIANGLES_HH
#define DENSITY_CONTOUR_TRIANGLES_HH

struct POINT3DID {
   unsigned int newID;
   float x, y, z;
};

typedef std::map<unsigned int, POINT3DID> ID2POINT3DID;

class TRIANGLE {
public:
   unsigned int pointID[3];
   bool reject_this; // set in post-processing
   clipper::Coord_orth mid_point;
   double back_front_projection_distance;
   clipper::Coord_orth normal_for_flat_shading;
   bool operator<(const TRIANGLE &t) const {
      return (back_front_projection_distance < t.back_front_projection_distance);
   }
   TRIANGLE() {
      reject_this = false;
   }
};

namespace coot {


   class density_contour_triangles_container_t {
   public:
      std::vector<clipper::Coord_orth> points;
      std::vector<clipper::Coord_orth> normals; // for Gouraud shading
      std::vector<TRIANGLE> point_indices;
      void depth_sort(const clipper::Coord_orth &back_plane_point,
                      const clipper::Coord_orth &front_plane_point);
      void calculate_normals(); // average normals on shared points
      void remove_small_triangles();
      bool empty() const { return (points.empty()); }
      void clear() {
         points.clear();
         normals.clear();
         point_indices.clear();
      }
   };

}

#endif // DENSITY_CONTOUR_TRIANGLES_HH
