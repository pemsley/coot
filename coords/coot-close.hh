
#include "clipper/core/coords.h"
#include "mmdb_manager.h"

// clipper::Coord_orth
// closest_approach(const clipper::Coord_orth &moving_point, 
// 		 const clipper::Coord_orth &reference_point,
// 		 const atom_selection_container_t &asc);

clipper::Coord_orth
closest_approach(const clipper::Coord_orth &moving_point, 
		 const clipper::Coord_orth &reference_point,
		 CMMDBManager *mol);

clipper::RTop_orth 
closest_approach_transformation(const clipper::Coord_orth &moving_point, 
				const clipper::Coord_orth &reference_point,
				CMMDBManager *mol);



