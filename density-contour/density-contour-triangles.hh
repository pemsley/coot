

#ifndef DENSITY_CONTOUR_TRIANGLES_HH
#define DENSITY_CONTOUR_TRIANGLES_HH

struct POINT3DID {
	unsigned int newID;
	float x, y, z;
};

typedef std::map<unsigned int, POINT3DID> ID2POINT3DID;

struct TRIANGLE {
	unsigned int pointID[3];
};

namespace coot { 

  class density_contour_triangles_container_t { 
    
  public:
    
    std::vector<clipper::Coord_orth> points;
    std::vector<TRIANGLE> point_indices;

  };

}

#endif // DENSITY_CONTOUR_TRIANGLES_HH
