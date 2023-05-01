
#include "coot-map-utils.hh"
#include "blob-line.hh"

std::pair<bool, clipper::Coord_orth>
coot::find_peak_along_line_favour_front(const clipper::Coord_orth &p1,
                                        const clipper::Coord_orth &p2,
                                        float contour_level,
                                        const clipper::Xmap<float> &xmap) {

   float best_score = -9999999.9;
   clipper::Coord_orth pbest;
   int istep_max = 500;
   bool point_set = false;

   for (int istep=0; istep<=istep_max; istep++) {
      float fr = float(istep)/float(istep_max);
      clipper::Coord_orth pc = p1 + fr*(p2-p1);
      float d = util::density_at_point(xmap, pc);
      if (d > contour_level) {
         // OK, so the point we want is somewhere in this peak
         for (int jstep=istep; jstep<=istep_max; jstep++) {
            fr = float(jstep)/float(istep_max);
            pc = p1 + fr*(p2-p1);
            d = util::density_at_point(xmap, pc);
            if (d > contour_level) {
               if (d> best_score) {
                  best_score = d;
                  pbest = pc;
                  point_set = true;
               }
            } else {
               // the front peak is over (now below the contour level), we have the pbest.
               break;
            }
         }
         break;
      }
   }
   return std::make_pair(point_set, pbest);
}
