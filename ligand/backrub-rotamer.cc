/* ligand/backrub-rotamer.cc
 * 
 * Copyright 2009 by The University of Oxford.
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */


#include <stdexcept>
#include <iostream>
#include <vector>

#include "richardson-rotamer.hh"
#include "utils/coot-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "backrub-rotamer.hh"

// thow an exception on failure to get a good result.
std::pair<coot::minimol::molecule,float>
coot::backrub::search(const coot::dictionary_residue_restraints_t &rest) {


    // waters either can clash and H-bond with the side chain  (mode 1)
    // or they should be ignored in the clash analysis and close waters
    // should be marked for deletion (mode 0).
   int water_interaction_mode = 0;
   coot::minimol::molecule mol;
   int n_vr = 15;                     // total is 2*n_vr+1
   double vector_rotation_range = 15; // degrees either side of the starting position.

   float best_score = -9999;
   std::vector<mmdb::Atom *> clashing_waters_for_best_score_local;
   coot::minimol::fragment best_frag;
   coot::richardson_rotamer r_rotamer(orig_this_residue, alt_conf, stored_mol, 0.1, 0);
   std::vector<float> pr = r_rotamer.probabilities();
   unsigned int n_rotatmers = pr.size();

   // First, where about in space is this residue centred?  Let's use
   // that to do an atom selection and a generation of the molecule
   // for a sphere of residues/atoms.  We can then use these sphere atoms to score the clashes

   clipper::Coord_orth rc = rotamer_residue_centre();
   float rr = residue_radius(rc);
   float rrr = rr + 6.0; // guess
   int SelectionHandle = stored_mol->NewSelection();

   stored_mol->SelectSphere(SelectionHandle, mmdb::STYPE_ATOM, rc.x(), rc.y(), rc.z(), rrr, mmdb::SKEY_OR);
   mmdb::PPAtom sphere_atoms = 0;
   int n_sphere_atoms = 0;
   stored_mol->GetSelIndex(SelectionHandle, sphere_atoms, n_sphere_atoms);

   atom_selection_container_t stored_mol_asc = make_asc(stored_mol);

   for (unsigned int irot=0; irot<n_rotatmers; irot++) {
      mmdb::Residue *r = r_rotamer.GetResidue(rest, irot);

      for (int ivr=-n_vr; ivr<=n_vr; ivr++) {

         double rotation_angle = vector_rotation_range * double(ivr)/double(n_vr);
         if (0) {
            std::cout << "DEBUG:: rotamer " << irot << " rotation_angle " << rotation_angle
                      << " from "
                      << vector_rotation_range << "*"  << double(ivr) << "/" << double(n_vr)
                      << std::endl;
         }
         coot::minimol::fragment frag = make_test_fragment(r, rotation_angle);
         float f = score_fragment(frag);
         std::pair<float, std::vector<mmdb::Atom *> > cs =
            get_clash_score(coot::minimol::molecule(frag), sphere_atoms, n_sphere_atoms,
                            water_interaction_mode);

         // the clash score is large and positive for big clashes.
         // clash score on to the same scale as density fit.
         float total_score =  -0.02 * cs.first + f;
         if (false) {
            std::cout << "   angle: " << rotation_angle << " density score " << f
                      << "  clash score " << cs.first << " total " << total_score << std::endl;
         }
         if (total_score > best_score) {
            best_score = total_score;
            best_frag  = frag;
            clashing_waters_for_best_score_local = cs.second;
         }
      }
      delete r;
   }
   stored_mol_asc.delete_atom_selection();

   if (best_frag.n_filled_residues() == 3) {
      best_frag.fragment_id = chain_id;
      mol.fragments.push_back(best_frag);
   } else {
      std::string mess = "  Failed to get a good fitting result";
      throw std::runtime_error(mess);
   }

   clashing_waters = clashing_waters_for_best_score_local;
   return std::pair<coot::minimol::molecule, float> (mol, best_score);
}


// rotation angles in degrees
//
// Throws an exception if fragment does not have 3 residues.
// 
coot::minimol::fragment
coot::backrub::make_test_fragment(mmdb::Residue *r, double rotation_angle) const {

   coot::minimol::fragment f(chain_id);
   std::vector<std::string> prev_res_atoms;
   std::vector<std::string> next_res_atoms;

   prev_res_atoms.push_back(" C  ");
   prev_res_atoms.push_back(" O  ");
   next_res_atoms.push_back(" N  ");
   next_res_atoms.push_back(" H  ");

   if (0) { // debug 
      std::cout << "DEBUG:: orig_prev_residue " << orig_prev_residue << std::endl;
      std::cout << "DEBUG:: orig_next_residue " << orig_next_residue << std::endl;
      if (orig_prev_residue)
         std::cout << "DEBUG:: orig_prev_residue "
                   << orig_prev_residue->GetChainID() << " " 
                   << orig_prev_residue->GetSeqNum() << " :" 
                   << orig_prev_residue->GetInsCode() << ":" 
                   << std::endl;
      if (orig_next_residue)
         std::cout << "DEBUG:: orig_prev_residue "
                   << orig_next_residue->GetChainID() << " " 
                   << orig_next_residue->GetSeqNum() << " :" 
                   << orig_next_residue->GetInsCode() << ":" 
                   << std::endl;
   }


   coot::minimol::residue this_residue(r);
   coot::minimol::residue prev_residue =
      make_residue_include_only(orig_prev_residue, prev_res_atoms);
   coot::minimol::residue next_residue =
      make_residue_include_only(orig_next_residue, next_res_atoms);

   // addresidue() fails when adding residues with the same residue
   // number (the ins code is ignored).  Urgh.
   //
   try {
      f.addresidue(prev_residue, 0);
      f.addresidue(coot::minimol::residue(r), 0);
      f.addresidue(next_residue, 0);
   }
   catch (const std::runtime_error &rte) {
      std::cout << "ERROR:: make_test_fragment() " << rte.what() << std::endl;
   }


   // now rotate fragment around the ca_prev -> ca_next vector
   clipper::Coord_orth dir = ca_next - ca_prev;

   for (int ires=f.min_res_no(); ires<=f.max_residue_number(); ires++) {
      for (unsigned int iat=0; iat<f[ires].n_atoms(); iat++) {
         clipper::Coord_orth pt(f[ires][iat].pos);
         // rotate pt
         double ra = M_PI*rotation_angle/180.0;
         clipper::Coord_orth pt_new =
            coot::util::rotate_around_vector(dir, pt, ca_prev, ra);
         f[ires][iat].pos = pt_new;
      }
   }

   rotate_individual_peptides_back_best(r, rotation_angle, &f);

   if (f.n_filled_residues() != 3) {
      std::string mess = "  Failed to get 3 residues with atoms in test fragment. Got ";
      mess += coot::util::int_to_string(f.n_filled_residues());
      if (false) {
         // debug frag
         coot::minimol::molecule mol_for_frag;
         mol_for_frag.fragments.push_back(f);
         mol_for_frag.write_file("test-frag.pdb", 0);
      }
      throw std::runtime_error(mess);
   }
   return f;
}

clipper::Coord_orth
coot::backrub::rotamer_residue_centre() const {

   mmdb::PPAtom residue_atoms;
   int n_residue_atoms;
   orig_this_residue->GetAtomTable(residue_atoms, n_residue_atoms);
   float sum_x=0, sum_y=0, sum_z=0;
   for (int iat=0; iat<n_residue_atoms; iat++) {
      sum_x += residue_atoms[iat]->x;
      sum_y += residue_atoms[iat]->y;
      sum_z += residue_atoms[iat]->z;
   }
   if (n_residue_atoms > 0) {
      float inv = 1.0/float(n_residue_atoms);
      clipper::Coord_orth pt(sum_x*inv, sum_y*inv, sum_z*inv);
      return pt;
   } else {
      return clipper::Coord_orth(0,0,0);
   }
} 

float
coot::backrub::residue_radius(const clipper::Coord_orth &rc) {

   float r = 0; 

   mmdb::PPAtom residue_atoms;
   int n_residue_atoms;
   orig_this_residue->GetAtomTable(residue_atoms, n_residue_atoms);
   float longest_length = 0.0;
   for (int iat=0; iat<n_residue_atoms; iat++) {
      clipper::Coord_orth pt(residue_atoms[iat]->x - rc.x(),
                             residue_atoms[iat]->y - rc.y(),
                             residue_atoms[iat]->z - rc.z());
      float this_length_sq = pt.lengthsq();
      if (this_length_sq > longest_length) {
         longest_length = this_length_sq;
         r = sqrt(longest_length);
      } 
   }

   return r;
} 



float
coot::backrub::score_fragment(minimol::fragment &frag) const {

   float d_score = 0;
   for (int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
      for (unsigned int iat=0; iat<frag[ires].n_atoms(); iat++) {
         float d = util::density_at_point(*xmap_p, frag[ires][iat].pos);
         d_score += d;
      }
   }
   return d_score;
} 

coot::minimol::residue
coot::backrub::make_residue_include_only(mmdb::Residue *orig_prev_residue,
                                         const std::vector<std::string> &prev_res_atoms) const {

   coot::minimol::residue r(orig_prev_residue, prev_res_atoms);
   return r;
} 



// Throw an exception of previous or next residues are null
//
// throw an exception on failure to find CA of prev or next.
void
coot::backrub::setup_this_and_prev_next_ca_positions() {

   // run through the atoms of both orig_prev_residue and
   // orig_next_residue looking for the CA atom in the same alt conf
   // with which this object was constructed.
   //
   // That sets ca_prev and ca_next.
   //
   // Throw an exception if we can't find either CA atoms.

   mmdb::PPAtom residue_atoms = 0;
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

   orig_this_residue->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::string atom_alt_conf(residue_atoms[iat]->altLoc);
      if (atom_name == " CA " ) {
         if (atom_alt_conf == alt_conf) {
            found = 1;
            ca_this = clipper::Coord_orth(residue_atoms[iat]->x,
                                          residue_atoms[iat]->y,
                                          residue_atoms[iat]->z);
         }
      }
   }

   if (! found) {
      std::string mess(" CA atom of this residue in alt conf \"");
      mess += alt_conf;
      mess += "\" not found";
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
      std::string mess(" CA atom of previous residue in alt conf \"");
      mess += alt_conf;
      mess += "\" not found";
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
      std::string mess(" CA atom of next residue in alt conf \"");
      mess += alt_conf;
      mess += "\" not found";
      throw std::runtime_error(mess);
   } 
}


// const because we are fiddling with f, not data members.
// 
void
coot::backrub::rotate_individual_peptides_back_best(mmdb::Residue *r, double rotation_angle,
                                                    coot::minimol::fragment *f) const {

   // Simple-minded optimisation.  Do a grid-search.  Pick the
   // smallest.  Could be speeded up with better optimisation.  The
   // function is 1-D and smooth.
   // 
   double best_back_rotation_leading =
      sample_individual_peptide(r, rotation_angle, f, orig_prev_residue, orig_this_residue, 1);
   double best_back_rotation_trailing = 
      sample_individual_peptide(r, rotation_angle, f, orig_this_residue, orig_next_residue, 0);

//    std::cout << "best_back_rotations: main: " << rotation_angle << " leading "
//              << best_back_rotation_leading << " trailing "
//              << best_back_rotation_trailing << std::endl;

   apply_back_rotation(f, 1, best_back_rotation_leading);
   apply_back_rotation(f, 0, best_back_rotation_trailing);
}

// This function is just to test that the back rotation of the
// individual peptide atoms.  The rotation angle is the angle around
// the major backrub vector (not the one to be applied to individual
// peptides).
//
// 
double
coot::backrub::sample_individual_peptide(mmdb::Residue *r, double rotation_angle,
                                         coot::minimol::fragment *f,
                                         mmdb::Residue *residue_front,
                                         mmdb::Residue *residue_back,
                                         bool is_leading_peptide_flag) const {

   // So f has had the big rotation applied.  Now we want to apply
   // compensating rotation of previous->this peptide.  What is the
   // best rotation given the rotation_angle?
   // 
   // Let's only consider the position of the O when doing that.
   //

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   clipper::Coord_orth O_pos;
   bool found_O_pos = 0;
   double best_back_rotation_angle = 0.0;
   
   residue_front->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string atom_name(residue_atoms[iat]->name);
      std::string atom_alt_conf(residue_atoms[iat]->altLoc);
      if (atom_name == " O  " ) {
         O_pos = clipper::Coord_orth(residue_atoms[iat]->x, 
                                     residue_atoms[iat]->y,
                                     residue_atoms[iat]->z);
         found_O_pos = 1;
      }
   }
   
   if (found_O_pos) { 

      int prev_resno = orig_prev_residue->GetSeqNum();
      int this_resno = orig_this_residue->GetSeqNum();
      int next_resno = orig_next_residue->GetSeqNum();

      // trailing peptide initially:
      clipper::Coord_orth ca_ref_rot = ca_this;
      clipper::Coord_orth dir = ca_next - ca_this;
      int res_offset = 1;
      int resno_for_CO = this_resno;
      int resno_for_N  = next_resno;
      if (is_leading_peptide_flag) { 
         ca_ref_rot = ca_prev;
         dir = ca_this - ca_prev;
         res_offset = 2;
         resno_for_CO = prev_resno;
         resno_for_N  = this_resno;
      } 
      
      double rar = 2 * fabs(rotation_angle); // rotation angle range
      double best_sum_dist = 9999999.9;
      for (double sample_rotation_angle=-rar;
           sample_rotation_angle<=rar;
           sample_rotation_angle += (rar*0.02 + 0.0001)) {
         double sum_dist = 0.0;
         for (int ires=(f->max_residue_number()-res_offset);
              ires<=(f->max_residue_number()+1-res_offset);
              ires++) {
            for (unsigned int iat=0; iat<(*f)[ires].n_atoms(); iat++) {
               if (ires == resno_for_CO) {
                  if ((*f)[ires][iat].name == " O  ") {
                     clipper::Coord_orth pt((*f)[ires][iat].pos);
                     // rotate pt
                     double ra = M_PI*sample_rotation_angle/180.0;
                     clipper::Coord_orth pt_new =
                        coot::util::rotate_around_vector(dir, pt, ca_ref_rot, ra);
                     double d =  clipper::Coord_orth::length(O_pos, pt_new);
                     sum_dist += d;
                  }
               }
            }
         }

         if (sum_dist < best_sum_dist) {
            best_sum_dist = sum_dist;
            best_back_rotation_angle = sample_rotation_angle;
         }
      }

      if (0) {
         // This is for generation of the correlation of the out of
         // plane distance and the best back rotation.  There is no
         // useful correlation, it turns out.
         
         // What is the distance of CA-this to the CA-prev->CA-next
         // vector, and how does it correlate with the gradient
         // of the best peptide rotate_angle vs major (backrub) rotation
         // angle?
         //
         // Let's call the CA vector distance d_out.
         // A = CA_this - CA_prev
         // B = CA_next - CA_prev
         // cos(theta) = A.B/(|A||B|)
         // d = |A|sin(theta)

         clipper::Coord_orth A = ca_this - ca_prev;
         clipper::Coord_orth B = ca_next - ca_prev;
         double A_length = sqrt(A.lengthsq());
         double B_length = sqrt(B.lengthsq());
         double cos_theta = clipper::Coord_orth::dot(A,B)/(A_length*B_length);
         double theta = acos(cos_theta);
         double d_out = A_length * sin(theta);
                     
         std::cout << "DEBUG:: best_sum_dist " << best_sum_dist << " with d_out "
                   << d_out << std::endl;
      }
   }
   return best_back_rotation_angle;
}


void
coot::backrub::apply_back_rotation(coot::minimol::fragment *f,
                                   bool is_leading_peptide_flag,
                                   double best_back_rotation_angle) const {

   int prev_resno = orig_prev_residue->GetSeqNum();
   int this_resno = orig_this_residue->GetSeqNum();
   int next_resno = orig_next_residue->GetSeqNum();
   
   // trailing peptide initially:
   clipper::Coord_orth ca_ref_rot = ca_this;
   clipper::Coord_orth dir = ca_next - ca_this;
   int res_offset = 1;
   int resno_for_CO = this_resno;
   int resno_for_N  = next_resno;
   if (is_leading_peptide_flag) { 
      ca_ref_rot = ca_prev;
      dir = ca_this - ca_prev;
      res_offset = 2;
      resno_for_CO = prev_resno;
      resno_for_N  = this_resno;
   }
   
   for (int ires=(f->max_residue_number()-res_offset);
        ires<=(f->max_residue_number()+1-res_offset);
        ires++) {
      for (unsigned int iat=0; iat<(*f)[ires].n_atoms(); iat++) {
         bool rotate_it = 0; 
         if (ires == resno_for_CO) {
            if (((*f)[ires][iat].name == " C  ") ||
                ((*f)[ires][iat].name == " O  ")) {
               rotate_it = 1;
            }
         }
         if (ires == resno_for_N) {
            if (((*f)[ires][iat].name == " N  ") ||
                ((*f)[ires][iat].name == " H  ")) {
               rotate_it = 1;
            }
         }
         if (rotate_it) {
            clipper::Coord_orth pt((*f)[ires][iat].pos);
            // rotate pt
            double ra = M_PI*best_back_rotation_angle/180.0;
            clipper::Coord_orth pt_new =
               coot::util::rotate_around_vector(dir, pt, ca_ref_rot, ra);
            (*f)[ires][iat].pos = pt_new;
         } 
      }
   }
}



std::vector<coot::atom_spec_t>
coot::backrub::waters_for_deletion() const {

   std::vector<atom_spec_t> baddies;
   for (unsigned int i=0; i<clashing_waters.size(); i++) {
      baddies.push_back(clashing_waters[i]);
   }
   return baddies;
}

// How are clashes scored?  Hmmm... well, I think no clashes at all should have a score
// of 0.0 (c.f. auto_best_fit_rotamer()).  I think a badly crashing residue should have a
// score of around 1000.  A single 2.0A crash will have a score of 16.7 and a 1.0A crash
// 66.7.
//
// The second is the clashing waters if water_interaction_mode is 0.
std::pair<float, std::vector<mmdb::Atom *> >
coot::get_clash_score(const coot::minimol::molecule &a_rotamer,
                      atom_selection_container_t asc, int water_interaction_mode) {

   // Is this the clash score function that you want?

   float score = 0;
   // float dist_crit = 2.7;
   float dist_crit = 2.1;
   std::vector<mmdb::Atom *> clashing_waters;

   // First, where is the middle of the rotamer residue atoms and what
   // is the mean and maximum distance of coordinates from that point?

   // double std_dev_residue_pos;

   std::pair<double, clipper::Coord_orth> rotamer_info = a_rotamer.get_pos();
   double max_dev_residue_pos = rotamer_info.first;
   clipper::Coord_orth mean_residue_pos = rotamer_info.second;
   if (rotamer_info.first < 0.0) {
      // there were no atoms then
      std::cout << "ERROR: clash score: there are no atoms in the residue" << std::endl;
   } else {

      // So now we know the centre of the residue and the maximum distance of one of its
      // atoms from the centre.  Now let's run over the atoms of the atom_sel.
      //
      // When we find a distance between the middle of the residue and an atom_sel atom
      // that is less than (max_dev_residue_pos + distance), then we have found
      // potentially clashing atoms, so check for a clash of that atom_sel atom with all
      // the atoms of the residue.
      double d;
      double d_atom;
      float badness;
      for (int i=0; i<asc.n_selected_atoms; i++) {
         mmdb::Atom *at = asc.atom_selection[i];
         clipper::Coord_orth atom_sel_atom(asc.atom_selection[i]->x,
                                           asc.atom_selection[i]->y,
                                           asc.atom_selection[i]->z);
         d = clipper::Coord_orth::length(atom_sel_atom, mean_residue_pos);
         if (d < (max_dev_residue_pos + dist_crit)) {
            for (unsigned int ifrag=0; ifrag<a_rotamer.fragments.size(); ifrag++) {
               for (int ires=a_rotamer[ifrag].min_res_no(); ires<=a_rotamer[ifrag].max_residue_number(); ires++) {
                  for (unsigned int iat=0; iat<a_rotamer[ifrag][ires].n_atoms(); iat++) {
                     d_atom = clipper::Coord_orth::length(a_rotamer[ifrag][ires][iat].pos,atom_sel_atom);
                     std::cout << "  d_atom " << d_atom << "\n";
                     if (d_atom < dist_crit) {
                        int atom_sel_atom_resno = asc.atom_selection[i]->GetSeqNum();
                        std::string atom_sel_atom_chain(asc.atom_selection[i]->GetChainID());

                        std::cout << "comparing rotamer chain :" << a_rotamer[ifrag].fragment_id << ": this res chain "
                                  << atom_sel_atom_chain << " and resnos "
                                  << ires << " with this resno " << atom_sel_atom_resno << std::endl;
                        if (! ((ires == atom_sel_atom_resno) &&
                               (a_rotamer[ifrag].fragment_id == atom_sel_atom_chain)) ) {
                           if ( (a_rotamer[ifrag][ires][iat].name != " N  ") &&
                                (a_rotamer[ifrag][ires][iat].name != " C  ") &&
                                (a_rotamer[ifrag][ires][iat].name != " CA ") &&
                                (a_rotamer[ifrag][ires][iat].name != " O  ") &&
                                (a_rotamer[ifrag][ires][iat].name != " H  ") ) {
                              bool is_water = false;
                              if (std::string(at->GetResName()) == "HOH") is_water = true;
                              if (! is_water || water_interaction_mode == 1) {
                                 badness = 100.0 * (1.0/d_atom - 1.0/dist_crit);
                                 std::cout << "adding badness " << badness << " for "
                                           <<  asc.atom_selection[i] << "\n";
                                  if (badness > 100.0)
                                     badness = 100.0;
                                 score += badness;
                              } else {
                                 // we don't count waters in the baddie score, just accumulate them
                                 // for deletion.
                                 clashing_waters.push_back(at);
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
   std::pair<float, std::vector<mmdb::Atom *> > p(score, clashing_waters);
   return p;
}


// How are clashes scored?  Hmmm... well, I think no clashes at all should have a score
// of 0.0 (c.f. auto_best_fit_rotamer()).  I think a badly crashing residue should have a
// score of around 1000.  A single 2.0A crash will have a score of 16.7 and a 1.0A crash
// 66.7.
//
//
// Uses class data orig_next_residue etc.
//
std::pair<float, std::vector<mmdb::Atom *> >
coot::backrub::get_clash_score(const coot::minimol::molecule &a_rotamer,
                               mmdb::PPAtom atom_selection, int n_sphere_atoms,
                               int water_interaction_mode) const {

   // Clashes to that are to H bonders have a different dist_crit
   // (smaller).  Tested on both atoms not being carbon (for now).
   // 
   std::vector<mmdb::Atom *> clashing_waters;
   float clash_score = 0.0;
   double dist_crit = 3.2;
   double dist_crit_H_bonder = 2.5; 
   double dist_crit_sq = dist_crit * dist_crit;
   double dist_crit_H_bonder_sq = dist_crit_H_bonder * dist_crit_H_bonder;
   std::pair<double, clipper::Coord_orth> rotamer_info = a_rotamer.get_pos();
   double max_dev_residue_pos = rotamer_info.first;
   clipper::Coord_orth mean_residue_pos = rotamer_info.second;
   int resno_1 = orig_prev_residue->GetSeqNum();
   int resno_2 = orig_this_residue->GetSeqNum();
   int resno_3 = orig_next_residue->GetSeqNum();
   if (rotamer_info.first < 0.0) {
      // there were no atoms then
      std::cout << "ERROR: clash score: there are no atoms in the residue" << std::endl;
   } else {

      for (int i=0; i<n_sphere_atoms; i++) {
         mmdb::Atom *at = atom_selection[i];
         clipper::Coord_orth atom_sel_atom_pos(atom_selection[i]->x,
                                               atom_selection[i]->y,
                                               atom_selection[i]->z);
         int atom_sel_resno = atom_selection[i]->GetSeqNum();
         std::string atom_sel_atom_chain(atom_selection[i]->GetChainID());
         std::string atom_sel_ele = atom_selection[i]->element;
         bool count_it = 1;
         if (chain_id == atom_sel_atom_chain) {
            if (atom_sel_resno==resno_1 || atom_sel_resno==resno_2 || atom_sel_resno==resno_3) {
               count_it = 0;
            }
         }
         if (count_it) {
            for (unsigned int ifrag=0; ifrag<a_rotamer.fragments.size(); ifrag++) {
               for (int ires=a_rotamer[ifrag].min_res_no(); ires<=a_rotamer[ifrag].max_residue_number(); ires++) {
                  for (unsigned int iat=0; iat<a_rotamer[ifrag][ires].n_atoms(); iat++) {
                     double dlsq = (a_rotamer[ifrag][ires][iat].pos-atom_sel_atom_pos).lengthsq();
                     if (dlsq <= 0.001)
                        dlsq = 0.001;
                     if (dlsq < dist_crit_sq) {
                        bool is_water = false;
                        if (std::string(at->GetResName()) == "HOH") is_water = true;
                        if (!is_water || water_interaction_mode == 1) {
                           // We take off  1.0/dist_crit_sq so that the clash score goes to 0 at dist_crit.
                           float extra_clash_score = 100 * (1.0/dlsq - 1.0/dist_crit_sq);
                           if (atom_sel_ele != " C" && a_rotamer[ifrag][ires][iat].element != " C") {
                              if (dlsq > dist_crit_H_bonder_sq)
                                 extra_clash_score = 0.0;
                              else
                                 extra_clash_score = 100 * (1.0/dlsq - 1.0/dist_crit_H_bonder_sq);
                           }
                           clash_score += extra_clash_score;
                           if (false)
                              std::cout << "   really adding badness " << extra_clash_score << " for distance "
                                        << sqrt(dlsq) << " \t" << atom_selection[i] << " to "
                                        << ires << " " << a_rotamer[ifrag][ires][iat].name
                                        << std::endl;
                         } else {
                            double d = sqrt(dlsq);
                            if (d < dist_crit_H_bonder) {
                               clashing_waters.push_back(at);
                            }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   std::pair<float, std::vector<mmdb::Atom *> > p(clash_score, clashing_waters);
   return p;
}
