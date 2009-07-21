

#include <stdexcept>
#include <iostream>
#include <vector>

#include "richardson-rotamer.hh"
#include "coot-map-utils.hh"
#include "backrub-rotamer.hh"

// thow an exception on failure to get a good result.
std::pair<coot::minimol::molecule,float>
coot::backrub::search() {

   coot::minimol::molecule mol;
   int n_vr = 100;                     // totat is 2*n_vr+2
   double vector_rotation_range = 120; // degrees either side of the starting position.

   float best_score = -9999;
   coot::minimol::fragment best_frag;
   coot::richardson_rotamer r_rotamer(orig_this_residue, stored_mol, 0.1, 0);
   std::vector<float> pr = r_rotamer.probabilities();
   int n_rotatmers = pr.size();

   for (unsigned int irot=0; irot<n_rotatmers; irot++) {
      CResidue *r = r_rotamer.GetResidue(irot);
      for (int ivr=-n_vr; ivr<=n_vr; ivr++) {
	 double rotation_angle = vector_rotation_range * double(ivr)/double(n_vr);
	 coot::minimol::fragment frag = make_test_fragment(r, rotation_angle);
	 float f = score_fragment(frag);
	 if (f > best_score) {
	    best_score = f;
	    best_frag  = frag;
	 }
      }
      delete r;
   }

   if (best_frag.n_filled_residues() == 3) {
      best_frag.fragment_id = chain_id;
      mol.fragments.push_back(best_frag);
   } else {
      std::string mess = "  Failed to get a good fitting result";
      throw std::runtime_error(mess);
   } 
   return std::pair<coot::minimol::molecule, float> (mol, best_score);
}


// rotation angles in degrees
// 
coot::minimol::fragment
coot::backrub::make_test_fragment(CResidue *r, double rotation_angle) const {

   coot::minimol::fragment f;
   std::vector<std::string> prev_res_atoms;
   std::vector<std::string> next_res_atoms;

   prev_res_atoms.push_back(" C  ");
   prev_res_atoms.push_back(" O  ");
   next_res_atoms.push_back(" N  ");
   next_res_atoms.push_back(" H  ");

   coot::minimol::residue this_residue(r);
   coot::minimol::residue prev_residue =
      make_residue_include_only(orig_prev_residue, prev_res_atoms);
   coot::minimol::residue next_residue =
      make_residue_include_only(orig_next_residue, next_res_atoms);
   f.addresidue(prev_residue, 0);
   f.addresidue(r, 0);
   f.addresidue(next_residue, 0);

   // now rotate fragment around the ca_prev -> ca_next vector

   for (unsigned int ires=f.min_res_no(); ires<=f.max_residue_number(); ires++) {
      for (unsigned int iat=0; iat<f[ires].n_atoms(); iat++) {
	 clipper::Coord_orth pt(f[ires][iat].pos);
	 // rotate pt
	 clipper::Coord_orth dir = ca_next - ca_prev;
	 double ra = M_PI*rotation_angle/180.0;
	 clipper::Coord_orth pt_new =
	    coot::util::rotate_round_vector(dir, pt, ca_prev, ra);
	 f[ires][iat].pos = pt_new;
      } 
   }

   return f;
}


float
coot::backrub::score_fragment(minimol::fragment &frag) const {

   float d_score = 0;
   for (unsigned int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
      for (unsigned int iat=0; iat<frag[ires].n_atoms(); iat++) {
	 float d = coot::util::density_at_point(xmap, frag[ires][iat].pos);
	 d_score += d;
      } 
   }
   return d_score;
} 

coot::minimol::residue
coot::backrub::make_residue_include_only(CResidue *orig_prev_residue,
					 const std::vector<std::string> &prev_res_atoms) const {

   coot::minimol::residue r(orig_prev_residue, prev_res_atoms);
   return r;
} 



// Throw an exception of previous or next residues are null
// 
// throw an exception on failure to find CA of prev or next.
void
coot::backrub::setup_prev_next_ca_positions() {

   // run through the atoms of both orig_prev_residue and
   // orig_next_residue looking for the CA atom in the same alt conf
   // with which this object was constructed.
   //
   // That sets ca_prev and ca_next.
   //
   // Throw an exception if we can't find either CA atoms.

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;

   short int found = 0;

   if (! orig_this_residue) {
      std::string mess(" Null this residue ");
      throw std::runtime_error(mess);
   } 
   if (! orig_prev_residue) {
      std::string mess(" Null previous residue ");
      throw std::runtime_error(mess);
   } 
   if (! orig_next_residue) {
      std::string mess(" Null next residue ");
      throw std::runtime_error(mess);
   } 
   
   orig_prev_residue->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::string atom_alt_conf(residue_atoms[iat]->altLoc);
      if (atom_name == " CA " ) {
	 if (atom_alt_conf == alt_conf) {
	    found = 1;
	    ca_prev = clipper::Coord_orth(residue_atoms[iat]->x,
					  residue_atoms[iat]->y,
					  residue_atoms[iat]->z);
	 }
      }
   }

   if (! found) {
      std::string mess(" CA atom of previous residue in alt conf ");
      mess += alt_conf;
      mess += " not found";
      throw std::runtime_error(mess);
   }

   found = 0;
   residue_atoms = 0;
   orig_next_residue->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::string atom_alt_conf(residue_atoms[iat]->altLoc);
      if (atom_name == " CA " ) {
	 if (atom_alt_conf == alt_conf) {
	    found = 1;
	    ca_next = clipper::Coord_orth(residue_atoms[iat]->x,
					  residue_atoms[iat]->y,
					  residue_atoms[iat]->z);
	 }
      }
   }

   if (! found) {
      std::string mess(" CA atom of next residue in alt conf ");
      mess += alt_conf;
      mess += " not found";
      throw std::runtime_error(mess);
   } 
} 
