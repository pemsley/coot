

#include <algorithm> // for std::find

#include "plane-utils.hh"
#include "coot-coord-utils.hh"

#include <compat/coot-sysdep.h>

// the vector points towards the ring
std::pair<bool, double>
coot::angle_betwen_plane_and_vector(mmdb::Residue *residue_p,
				    const std::vector<std::string> &ring_atom_names,
				    const std::string &altconf_in,
				    const clipper::Coord_orth &vector) {

   std::pair<bool, double> r(false,0);

   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<clipper::Coord_orth> ring_atom_positions;
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 std::string atom_name(at->GetAtomName());
	 std::string alt_conf(at->altLoc);
	 std::vector<std::string>::const_iterator it =
	    std::find(ring_atom_names.begin(), ring_atom_names.end(), atom_name);
	 if (it != ring_atom_names.end()) {
	    if (alt_conf == altconf_in) {
	       clipper::Coord_orth pos = co(at);
	       ring_atom_positions.push_back(pos);
	    }
	 }

	 if (ring_atom_positions.size() > 4) {

	    lsq_plane_info_t l(ring_atom_positions);
	    double angle = l.angle(vector); // degrees
	    r.first = true;
	    r.second = angle;
	 }
      }
   }

   return r;

}

std::pair<bool, double>
coot::angle_betwen_plane_and_vector(mmdb::Residue *residue_p,
				    mmdb::Atom *atom_in_ring,
				    mmdb::Atom *bonding_atom) {

   std::pair<bool, double> r(false,0);

   clipper::Coord_orth pt_1 = co(atom_in_ring);
   clipper::Coord_orth pt_2 = co(bonding_atom);
   clipper::Coord_orth dv = pt_1 - pt_2;

   std::string alt_conf(bonding_atom->altLoc);
   std::string res_name(residue_p->GetResName());

   std::vector<std::string> ring_atom_names;

   // maybe look up the group of the comp_id and use that for a mor general
   // test.

   if (res_name == "NAG") {
      ring_atom_names.push_back(" C1 ");
      ring_atom_names.push_back(" C2 ");
      ring_atom_names.push_back(" C3 ");
      ring_atom_names.push_back(" C4 ");
      ring_atom_names.push_back(" C5 ");
      ring_atom_names.push_back(" O5 ");
   }
   if (res_name == "MAN" || res_name == "BMA" || res_name == "GLC" || res_name == "GAL") {
      ring_atom_names.push_back(" C1 ");
      ring_atom_names.push_back(" C2 ");
      ring_atom_names.push_back(" C3 ");
      ring_atom_names.push_back(" C4 ");
      ring_atom_names.push_back(" C5 ");
      ring_atom_names.push_back(" O5 ");
   }
   if (res_name == "SIA") {
      ring_atom_names.push_back(" C2 ");
      ring_atom_names.push_back(" C3 ");
      ring_atom_names.push_back(" C4 ");
      ring_atom_names.push_back(" C5 ");
      ring_atom_names.push_back(" C6 ");
      ring_atom_names.push_back(" O6 ");
   }
   r = angle_betwen_plane_and_vector(residue_p, ring_atom_names, alt_conf, dv);
   return r;
}
