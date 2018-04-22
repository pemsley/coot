
#include "ca-torsion-info.hh"
#include "coot-utils/coot-coord-utils.hh"

std::pair<bool, double>
coot::CA_torsion_info_t::score_fragment(const coot::minimol::fragment &frag,
					mmdb::Residue *residue_p,
					mmdb::Residue *res_prev_p,
					int seqnum, int offset) {
   
   std::pair<bool, double> score(false,0);

   // extract the 4 CA positions!

   // we can't catch "bad" access attempts for residues and atoms, baah.
   
   std::pair<bool, clipper::Coord_orth> CA_0;
   std::pair<bool, clipper::Coord_orth> CA_1;
   std::pair<bool, clipper::Coord_orth> CA_2;
   std::pair<bool, clipper::Coord_orth> CA_3;
   std::pair<bool, clipper::Coord_orth> CA_4;

   std::string CA = " CA "; // PDBv3 FIXME

   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 std::string atom_name(at->GetAtomName());
	 if (atom_name == " CA ") {
	    clipper::Coord_orth pt = co(at);
	    CA_2 = std::pair<bool, clipper::Coord_orth> (true, pt);
	 }
      }
   }

   if (offset == 1) {
      if (res_prev_p) {
	 mmdb::Atom **residue_atoms = 0;
	 int n_residue_atoms;
	 res_prev_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    mmdb::Atom *at = residue_atoms[iat];
	    std::string atom_name(at->GetAtomName());
	    if (atom_name == " CA ") {
	       clipper::Coord_orth pt = co(at);
	       CA_1 = std::pair<bool, clipper::Coord_orth> (true, pt);
	    }
	 }
      }
   }

   for (int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
      const minimol::residue &res = frag[ires];
      for (std::size_t iat=0; iat<res.atoms.size(); iat++) {
	 if (res.atoms[iat].name == CA) {
	    if (ires == (seqnum - 3)) CA_0 = std::pair<bool, clipper::Coord_orth> (true, res.atoms[iat].pos);
	    if (ires == (seqnum - 2)) CA_1 = std::pair<bool, clipper::Coord_orth> (true, res.atoms[iat].pos);
	    if (ires == (seqnum - 1)) CA_2 = std::pair<bool, clipper::Coord_orth> (true, res.atoms[iat].pos);
	    if (ires == (seqnum    )) CA_3 = std::pair<bool, clipper::Coord_orth> (true, res.atoms[iat].pos);
	    if (ires == (seqnum + 1)) CA_4 = std::pair<bool, clipper::Coord_orth> (true, res.atoms[iat].pos);
	 }
      }
   }

   if (false)
      std::cout << "debug:: founds: "
		<< CA_0.first << " " << CA_1.first << " "
		<< CA_2.first << " " << CA_3.first << " "
		<< CA_4.first << " " << std::endl;

   if (CA_1.first && CA_2.first && CA_3.first && CA_4.first) {

      double a = clipper::Util::rad2d(clipper::Coord_orth::angle(CA_2.second, CA_3.second, CA_4.second));
      double t = clipper::Util::rad2d(clipper::Coord_orth::torsion(CA_1.second, CA_2.second, CA_3.second, CA_4.second));

      float pr = ai.prob_angle_torsion(a,t);
      // std::cout << "here with angle " << a << " torsion " << t << " pr " << pr << std::endl;
      score.first = true;
      score.second = pr;
   }

   return score;

}
