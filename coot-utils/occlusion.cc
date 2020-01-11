

#include <set>
#include "occlusion.hh"

void coot::set_lower_left_and_range(const std::vector<occlusion_triangle> &tris, const std::vector<clipper::Coord_orth> &positions,
                                    float brick_size,
                                    clipper::Coord_orth *lower_left_p, int *brick_range) {

   *lower_left_p = clipper::Coord_orth(0,0,0);
   if (tris.size() > 0) {
      clipper::Coord_orth ll( 99990,  99990,  99990); // lower left
      clipper::Coord_orth tr(-99990, -99990, -99990); // top right
      for (unsigned int i=0; i<tris.size(); i++) {
         const occlusion_triangle &t = tris[i];
         if (t.mid_point.x() < ll.x()) ll = clipper::Coord_orth(t.mid_point.x(), ll.y(), ll.z());
         if (t.mid_point.y() < ll.y()) ll = clipper::Coord_orth(ll.x(), t.mid_point.y(), ll.z());
         if (t.mid_point.z() < ll.z()) ll = clipper::Coord_orth(ll.x(), ll.y(), t.mid_point.z());
         if (t.mid_point.x() > tr.x()) tr = clipper::Coord_orth(t.mid_point.x(), tr.y(), tr.z());
         if (t.mid_point.y() > tr.y()) tr = clipper::Coord_orth(tr.x(), t.mid_point.y(), tr.z());
         if (t.mid_point.z() > tr.z()) tr = clipper::Coord_orth(tr.x(), tr.y(), t.mid_point.z());
      }
      *lower_left_p = ll;
      brick_range[0] = static_cast<int>((tr.x() - ll.x())/brick_size) + 1;
      brick_range[1] = static_cast<int>((tr.y() - ll.y())/brick_size) + 1;
      brick_range[2] = static_cast<int>((tr.z() - ll.x())/brick_size) + 1;
   }
}

void coot::set_lower_left_and_range(const std::vector<augmented_position> &positions,
                                    float brick_size,
                                    clipper::Coord_orth *lower_left_p, int *brick_range) {

   *lower_left_p = clipper::Coord_orth(0,0,0);
   if (positions.size() > 0) {
      clipper::Coord_orth ll( 99990,  99990,  99990); // lower left
      clipper::Coord_orth tr(-99990, -99990, -99990); // top right
      for (unsigned int i=0; i<positions.size(); i++) {
         const clipper::Coord_orth &p = positions[i].position;
         if (p.x() < ll.x()) ll = clipper::Coord_orth(p.x(), ll.y(), ll.z());
         if (p.y() < ll.y()) ll = clipper::Coord_orth(ll.x(), p.y(), ll.z());
         if (p.z() < ll.z()) ll = clipper::Coord_orth(ll.x(), ll.y(), p.z());
         if (p.x() > tr.x()) tr = clipper::Coord_orth(p.x(), tr.y(), tr.z());
         if (p.y() > tr.y()) tr = clipper::Coord_orth(tr.x(), p.y(), tr.z());
         if (p.z() > tr.z()) tr = clipper::Coord_orth(tr.x(), tr.y(), p.z());
      }
   }
}


void coot::set_occlusions(std::vector<occlusion_triangle> &tris, const std::vector<clipper::Coord_orth> &positions) {

   // each triangle in tris has the centre position, area and normal filled and needs to have the occlusion factor filled.
   // and that is what this function does.

   float occlusion_limit = 8.8; // only consider triangles that are closer than this distance
   float brick_size = 10.0;

   clipper::Coord_orth lower_left;
   int brick_range[3]; // how many bricks in each dimension

   set_lower_left_and_range(tris, positions, brick_size, &lower_left, brick_range);
   std::cout << "brick ranges: " << brick_range[0] << " " << brick_range[1] << " " << brick_range[2] << std::endl;

   std::vector<std::set<unsigned int> > triangles_in_bricks;

   float inv_brick_size = 1.0/brick_size;
   for (unsigned int i=0; i<tris.size(); i++) {
      const occlusion_triangle &t = tris[i];
      int idx_3d[3];
      idx_3d[0] = static_cast<int> ((t.mid_point.x() - lower_left[0]) * inv_brick_size);
      idx_3d[1] = static_cast<int> ((t.mid_point.y() - lower_left[1]) * inv_brick_size);
      idx_3d[2] = static_cast<int> ((t.mid_point.z() - lower_left[2]) * inv_brick_size);
   }


}


// set the occlusion factor on the positions - we don't care about triangles
//
void coot::set_occlusions(std::vector<augmented_position>  &positions) {

   // double because we are using lengths of clipper Coord_orths
   double occlusion_limit = 8.8; // only consider triangles that are closer than this distance
   float brick_size = 10.0;

   clipper::Coord_orth lower_left;
   int brick_range[3]; // how many bricks in each dimension

   set_lower_left_and_range(positions, brick_size, &lower_left, brick_range);
   std::cout << "brick ranges: " << brick_range[0] << " " << brick_range[1] << " " << brick_range[2] << std::endl;
   std::vector<std::set<unsigned int> > positions_in_bricks;
   positions_in_bricks.resize(brick_range[0] * brick_range[1] * brick_range[2]);

   float inv_brick_size = 1.0/brick_size;
   fill_the_bricks(positions, brick_size, brick_range, lower_left, &positions_in_bricks);
   occlusion_of_positions_within_bricks(positions_in_bricks, positions, occlusion_limit);
   occlusion_of_positions_between_bricks(positions_in_bricks, positions, occlusion_limit, brick_range);
}

unsigned int
coot::occlusion_idx_3d_to_idx_1d(const int idx_3d[3], const int *range) {

   unsigned int idx = range[0] * range[1] * idx_3d[2] + range[1] * idx_3d[1] + idx_3d[0];
   return idx;

}

// fill positions_in_bricks
void coot::fill_the_bricks(const std::vector<coot::augmented_position> &positions, float brick_size, int *brick_range_p,
                           const clipper::Coord_orth &lower_left,
                           std::vector<std::set<unsigned int> > *positions_in_bricks_p) {

   float inv_brick_size = 1.0/brick_size;
   for (unsigned int i=0; i<positions.size(); i++) {
      const clipper::Coord_orth &p = positions[i].position;
      int idx_3d[3];
      idx_3d[0] = static_cast<int> ((p.x() - lower_left[0]) * inv_brick_size);
      idx_3d[1] = static_cast<int> ((p.y() - lower_left[1]) * inv_brick_size);
      idx_3d[2] = static_cast<int> ((p.z() - lower_left[2]) * inv_brick_size);
      int idx_1d = occlusion_idx_3d_to_idx_1d(idx_3d, brick_range_p);
      positions_in_bricks_p->at(idx_1d).insert(i);
   }

}

void 
coot::occlusion_of_positions_between_bricks(const std::vector<std::set<unsigned int> > &bricks,
                                            std::vector<augmented_position> &positions,
                                            double occlusion_limit, const int *brick_range_p) {

   int brick_index_max = brick_range_p[0] * brick_range_p[1] * brick_range_p[2];

   double occ_lim_sqrd = occlusion_limit * occlusion_limit;
   int bricks_size = bricks.size();
   for (int ib=0; ib<bricks_size; ib++) {
      const std::set<unsigned int> &brick_base = bricks[ib];
      for (int iz=-1; iz<2; iz++) {
         for (int iy= -1; iy<2; iy++) {
            for (int ix= -1; ix<2; ix++) {
               int ib_neighb = ib + ix + iy * brick_range_p[0] + iz * brick_range_p[0] * brick_range_p[1];
               if ((ib_neighb >= 0) && (ib_neighb != ib)) {
                  if (ib_neighb < brick_index_max) {
                     std::set<unsigned int>::const_iterator it_base;
                     std::set<unsigned int>::const_iterator it_neighb;
                     const std::set<unsigned int> &brick_neighb = bricks[ib];
                     for (it_base=brick_base.begin(); it_base!=brick_base.end(); it_base++) {
                        const clipper::Coord_orth &pt_1 = positions[*it_base].position;
                        for (it_neighb=brick_neighb.begin(); it_neighb!=brick_neighb.end(); it_neighb++) {
                           const clipper::Coord_orth &pt_2 = positions[*it_neighb].position;
                           clipper::Coord_orth delta_vector(pt_2-pt_1);
                           double dd = delta_vector.lengthsq();
                           if (dd < occ_lim_sqrd) {
                              double dp_1 = clipper::Coord_orth::dot(positions[*it_base].normal, delta_vector);
                              if (dp_1 > 0.0) {
                                // surface is concave, therefore occluding...
                                double dp_2 = clipper::Coord_orth::dot(positions[*it_base].normal, positions[*it_neighb].normal);
                                double d = sqrt(dd);
                                if (d < 1.0) d = 1.0;
                                double occlusion_factor = 0.5 * (1.0 + dp_2)/d;
                                // now add occlusion_factor to positions[i].
                                positions[*it_base].occlusion_factor += occlusion_factor;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

}


void 
coot::occlusion_of_positions_within_bricks(const std::vector<std::set<unsigned int> > &bricks,
                                           std::vector<augmented_position> &positions,
                                           double occlusion_limit) {

   double occ_lim_sqrd = occlusion_limit * occlusion_limit;
   for (unsigned int i=0; i<bricks.size(); i++) {
      const std::set<unsigned int> &indices = bricks[i];
      std::set<unsigned int>::const_iterator it_1;
      for (it_1=indices.begin(); it_1 != indices.end(); it_1++) {
         const clipper::Coord_orth &pt_1 = positions[*it_1].position;
         std::set<unsigned int>::const_iterator it_2;
         for (it_2=indices.begin(); it_2 != indices.end(); it_2++) {
            if (it_1 != it_2) {
               const clipper::Coord_orth &pt_2 = positions[*it_2].position;
               clipper::Coord_orth delta_vector(pt_2-pt_1);
               double dd = delta_vector.lengthsq();
               if (dd < occ_lim_sqrd) {
                  double dp_1 = clipper::Coord_orth::dot(positions[*it_1].normal, delta_vector);
                  if (dp_1 > 0.0) {
                     // surface is concave, therefore occluding...
                     double dp_2 = clipper::Coord_orth::dot(positions[*it_2].normal, positions[*it_2].normal);
                     double d = sqrt(dd);
                     if (d < 1.0) d = 1.0;
                     double occlusion_factor = 0.5 * (1.0 + dp_2)/d;
                     // now add occlusion_factor to positions[i].
                     positions[*it_1].occlusion_factor += occlusion_factor;
                     positions[*it_2].occlusion_factor += occlusion_factor;
                  }
               }
            }
         }
      }
   }
}

