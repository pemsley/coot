
#include "loop-path.hh"
#include "coot-utils/coot-coord-utils.hh"

std::pair<bool, std::vector<coot::CartesianPair> >
coot::loop_path(mmdb::Atom *start_back_2,
		mmdb::Atom *start,
		mmdb::Atom *end,
		mmdb::Atom *end_plus_2,
		unsigned int n_line_segments) {

   std::vector<CartesianPair> loop_line_segments;
   bool needs_CA_CA_badness_spots = false;
   std::pair<bool, std::vector<coot::CartesianPair> > blank(false, loop_line_segments);

   if (! start_back_2) return blank;
   if (! start) return blank;
   if (! end) return blank;
   if (! end_plus_2) return blank;

   // sane input


   // 20190118-PE, we want to add the John Berrisford request: that residues
   // that do not have enough of a gap in the residue numbering for the distance
   // between them should be flagged with something big and red.
   // Also, they will be drawn as a straight line, not a spline.

   int res_no_start = start->residue->GetSeqNum();
   int res_no_end   =   end->residue->GetSeqNum();
   int res_no_delta = res_no_end - res_no_start;
   n_line_segments = 2 * res_no_delta; // don't listen to call parameter!
   if (n_line_segments < 4) n_line_segments = 4; // sanitize

   clipper::Coord_orth P0 = co(start_back_2);
   clipper::Coord_orth P1 = co(start);
   clipper::Coord_orth P4 = co(end);
   clipper::Coord_orth P5 = co(end_plus_2);

   // Now make P2 = P1 + s(P1 - P0);
   // Now make P3 = P4 + s(P4 - P5);

   double d2 = sqrt(clipper::Coord_orth(P1-P4).lengthsq());

   bool sird = is_sane_inter_residue_distance(d2, res_no_delta);

   if (sird) {
      double d = 0.4 * (d2); // this number could be optimized
      if (d < 0.50) d = 0.50; // and this one
      if (d > 30.0) d = 10.0; // and this one

      double s = sqrt(d);
      clipper::Coord_orth P2 = P1 + s * P1 - s * P0;
      clipper::Coord_orth P3 = P4 + s * P4 - s * P5;

      if (false) {
	 std::cout << "oo = new_generic_object_number('points')" << std::endl;
	 std::cout << "to_generic_object_add_point(oo, \"green\", 4, " << P1.x() << ", " << P1.y() << ", " << P1.z() << ")" << std::endl;
	 std::cout << "to_generic_object_add_point(oo, \"green\", 4, " << P2.x() << ", " << P2.y() << ", " << P2.z() << ")" << std::endl;
	 std::cout << "to_generic_object_add_point(oo, \"green\", 4, " << P3.x() << ", " << P3.y() << ", " << P3.z() << ")" << std::endl;
	 std::cout << "to_generic_object_add_point(oo, \"green\", 4, " << P4.x() << ", " << P4.y() << ", " << P4.z() << ")" << std::endl;
	 std::cout << "set_display_generic_object(oo, 1)" << std::endl;
      }

      unsigned int n_pts = 2 * n_line_segments;

      for (unsigned int i=0; i<n_pts; i+=2) {
	 double t = static_cast<float>(i)/static_cast<float>(n_pts);
	 clipper::Coord_orth comp_1 = (1.0-t)*(1.0-t)*(1.0-t)*P1;
	 clipper::Coord_orth comp_2 = 3.0*(1.0-t)*(1.0-t)*t*P2;
	 clipper::Coord_orth comp_3 = 3.0*(1.0-t)*t*t*P3;
	 clipper::Coord_orth comp_4 = t*t*t*P4;
	 clipper::Coord_orth ls_start = comp_1 + comp_2 + comp_3 + comp_4;
	 t = static_cast<float>(i+1)/static_cast<float>(n_pts);
	 comp_1 = (1.0-t)*(1.0-t)*(1.0-t)*P1;
	 comp_2 = 3.0*(1.0-t)*(1.0-t)*t*P2;
	 comp_3 = 3.0*(1.0-t)*t*t*P3;
	 comp_4 = t*t*t*P4;
	 clipper::Coord_orth ls_end = comp_1 + comp_2 + comp_3 + comp_4;
	 Cartesian ls_start_c(ls_start);
	 Cartesian ls_end_c(ls_end);
	 loop_line_segments.push_back(CartesianPair(ls_start_c, ls_end_c));
      }
   } else {

      // non-sane distance between residues (i.e. the distance between residues is too short for the
      // residue number distance)

      needs_CA_CA_badness_spots = true;
      unsigned int n_pts = 4 * n_line_segments;
      double recip_res_no_delta = 1.0/static_cast<double>(n_pts);
      clipper::Coord_orth ls_delta((P4-P1) * recip_res_no_delta);
      for (unsigned int i=0; i<n_pts; i+=2) {

	 clipper::Coord_orth ls_start = P1 + i * ls_delta;
	 clipper::Coord_orth ls_end   = ls_start + ls_delta;
	 Cartesian ls_start_c(ls_start);
	 Cartesian ls_end_c(ls_end);

	 // these loop line segments have points added by calling function
	 // (currently rendered as orange discs)
	 //
	 loop_line_segments.push_back(CartesianPair(ls_start_c, ls_end_c));
      }
   }
   return std::pair<bool, std::vector<coot::CartesianPair> >(needs_CA_CA_badness_spots, loop_line_segments);
}

// needs extra arg for P-P vs CA-CA
bool
coot::is_sane_inter_residue_distance(double dist_between_residues, int res_no_delta) {

   bool status = true;

   // return false if the residue number difference is too small for the
   // position difference of the loop residues

   // old compiler
   double dist_crit = static_cast<double>(std::abs(static_cast<double>(res_no_delta))) * 3.7;

   if (dist_between_residues > dist_crit)
      status = false;

   if (false)
      std::cout << "debug:: is_sane_inter_residue_distance() returns " << status
		<< " based on " << res_no_delta << " -> " << dist_crit << " vs "
		<< dist_between_residues << std::endl;

   return status;
}
