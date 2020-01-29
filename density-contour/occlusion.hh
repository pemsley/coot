
#include <set>

#include <clipper/core/coords.h>

namespace coot {

	class occlusion_triangle {
		bool occlusion_is_set;
           public:
		int idx_1, idx_2, idx3;
		clipper::Coord_orth mid_point;
                float surface_area;
                clipper::Coord_orth normal; // normalized normal
		float occlusion_factor;
	};

         class augmented_position {
         public:
            clipper::Coord_orth position;
            clipper::Coord_orth normal; // normalized normal
            float occlusion_factor;
            augmented_position(const clipper::Coord_orth &p, const clipper::Coord_orth &n) : position(p), normal(n) {
              occlusion_factor = 0.0;
            }
            augmented_position() {}
         };

        // this presumes that the surface area and the normal have been calculated/set for the triangles
	void set_occlusions(std::vector<occlusion_triangle> &tris, const std::vector<clipper::Coord_orth> &positions);
        void set_lower_left_and_range(const std::vector<occlusion_triangle> &tris, const std::vector<clipper::Coord_orth> &positions,
                                      float brick_size,
                                      clipper::Coord_orth *lower_left_p, int *brick_range_p);
        // set the occlusion factor on the positions - we don't core about triangles
        void set_occlusions(std::vector<augmented_position> &positions);
        void set_lower_left_and_range(const std::vector<augmented_position> &positions,
                                    float brick_size,
                                    clipper::Coord_orth *lower_left_p, int *brick_range);
       // fill positions_in_bricks
       void fill_the_bricks(const std::vector<augmented_position> &positions, float brick_size, int *brick_range_p, 
                            const clipper::Coord_orth &lower_left,
                            std::vector<std::set<unsigned int> > *positions_in_bricks_p);

       unsigned int occlusion_idx_3d_to_idx_1d(const int idx_3d[3], const int *brick_range);

       void occlusion_of_positions_within_bricks(const std::vector<std::set<unsigned int> > &indices,
                                                 std::vector<augmented_position> &positions, // change positions by adding occlusion
                                                 double occlusion_limit);
       void occlusion_of_positions_between_bricks(const std::vector<std::set<unsigned int> > &bricks,
                                                  std::vector<augmented_position> &positions,  // change positions y adding occlusion
                                                  double occlusion_limit, const int *brick_range_p);




}

