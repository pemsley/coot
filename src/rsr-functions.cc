
#include "graphics-info.h"
#include "rsr-functions.hh"

// this is not in a header, it seems.
coot::refinement_results_t
refine_residues_with_alt_conf(int imol, const std::vector<coot::residue_spec_t> &residue_specs,
			      const std::string &alt_conf);


void regularize_residue() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         std::vector<mmdb::Residue *> rv;
         rv.push_back(at->residue);
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = g.regularize_residues_vec(imol, rv, alt_conf, mol);
      }
   }
}

void regularize_tandem_3() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<mmdb::Residue *> rv;
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         mmdb::Residue *r_p_1 = coot::util::get_following_residue(rspec, mol);
         mmdb::Residue *r_p_2 = coot::util::get_following_residue(coot::residue_spec_t(r_p_1), mol);
         mmdb::Residue *r_p_3 = coot::util::get_following_residue(coot::residue_spec_t(r_p_2), mol);
         mmdb::Residue *r_m_1 = coot::util::get_previous_residue(rspec, mol);
         mmdb::Residue *r_m_2 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_1), mol);
         mmdb::Residue *r_m_3 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_2), mol);
         rv.push_back(r_m_3);
         rv.push_back(r_m_2);
         rv.push_back(r_m_1);
         rv.push_back(at->residue);
         rv.push_back(r_p_1);
         rv.push_back(r_p_2);
         rv.push_back(r_p_3);
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = g.regularize_residues_vec(imol, rv, alt_conf, mol);
      }
   }
}

void regularize_sphere() {

   float radius = 6.6;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         std::vector<coot::residue_spec_t> v = graphics_info_t::molecules[imol].residues_near_residue(rspec, radius);
         v.push_back(rspec);
         std::vector<mmdb::Residue *> rv;
         for (unsigned int i=0; i<v.size(); i++) {
            mmdb::Residue *r = coot::util::get_residue(v[i], mol);
            if (r) rv.push_back(r);
         }
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = g.regularize_residues_vec(imol, rv, alt_conf, mol);
      }
   }
}

void rsr_refine_residue() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<coot::residue_spec_t> v;
         v.push_back(rspec);
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = refine_residues_with_alt_conf(imol, v, alt_conf);
      }
   }
}

void rsr_refine_tandem_5() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<coot::residue_spec_t> v;
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         mmdb::Residue *r_p_1 = coot::util::get_following_residue(rspec, mol);
         mmdb::Residue *r_p_2 = coot::util::get_following_residue(coot::residue_spec_t(r_p_1), mol);
         mmdb::Residue *r_p_3 = coot::util::get_following_residue(coot::residue_spec_t(r_p_2), mol);
         mmdb::Residue *r_p_4 = coot::util::get_following_residue(coot::residue_spec_t(r_p_3), mol);
         mmdb::Residue *r_p_5 = coot::util::get_following_residue(coot::residue_spec_t(r_p_4), mol);
         mmdb::Residue *r_m_1 = coot::util::get_previous_residue(rspec, mol);
         mmdb::Residue *r_m_2 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_1), mol);
         mmdb::Residue *r_m_3 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_2), mol);
         mmdb::Residue *r_m_4 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_3), mol);
         mmdb::Residue *r_m_5 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_4), mol);
         v.push_back(coot::residue_spec_t(r_m_5));
         v.push_back(coot::residue_spec_t(r_m_4));
         v.push_back(coot::residue_spec_t(r_m_3));
         v.push_back(coot::residue_spec_t(r_m_2));
         v.push_back(coot::residue_spec_t(r_m_1));
         v.push_back(rspec);
         v.push_back(coot::residue_spec_t(r_p_1));
         v.push_back(coot::residue_spec_t(r_p_2));
         v.push_back(coot::residue_spec_t(r_p_3));
         v.push_back(coot::residue_spec_t(r_p_4));
         v.push_back(coot::residue_spec_t(r_p_5));
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = refine_residues_with_alt_conf(imol, v, alt_conf);
      }
   }
}

void rsr_refine_tandem_3() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<coot::residue_spec_t> v;
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         mmdb::Residue *r_p_1 = coot::util::get_following_residue(rspec, mol);
         mmdb::Residue *r_p_2 = coot::util::get_following_residue(coot::residue_spec_t(r_p_1), mol);
         mmdb::Residue *r_p_3 = coot::util::get_following_residue(coot::residue_spec_t(r_p_2), mol);
         mmdb::Residue *r_m_1 = coot::util::get_previous_residue(rspec, mol);
         mmdb::Residue *r_m_2 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_1), mol);
         mmdb::Residue *r_m_3 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_2), mol);
         v.push_back(coot::residue_spec_t(r_m_3));
         v.push_back(coot::residue_spec_t(r_m_2));
         v.push_back(coot::residue_spec_t(r_m_1));
         v.push_back(rspec);
         v.push_back(coot::residue_spec_t(r_p_1));
         v.push_back(coot::residue_spec_t(r_p_2));
         v.push_back(coot::residue_spec_t(r_p_3));
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = refine_residues_with_alt_conf(imol, v, alt_conf);
      }
   }
   
}

void rsr_refine_chain() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         mmdb::Chain *chain_p = at->residue->chain;
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<mmdb::Residue *> residues = coot::util::residues_in_chain(chain_p);
         std::vector<coot::residue_spec_t> v;
         for (unsigned int i=0; i<residues.size(); i++) {
            coot::residue_spec_t spec(residues[i]);
            v.push_back(spec);
         }
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = refine_residues_with_alt_conf(imol, v, alt_conf);
      }
   }
}


void rsr_refine_all_atoms() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {

         // add a fix for all chains - needs reworking.

         mmdb::Chain *chain_p = at->residue->chain;
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<mmdb::Residue *> residues = coot::util::residues_in_chain(chain_p);
         std::vector<coot::residue_spec_t> v;
         for (unsigned int i=0; i<residues.size(); i++) {
            coot::residue_spec_t spec(residues[i]);
            v.push_back(spec);
         }
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = refine_residues_with_alt_conf(imol, v, alt_conf);
      }
   }

}

void rsr_sphere_refine_plus() {

   float radius = 6.6;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<coot::residue_spec_t> v = graphics_info_t::molecules[imol].residues_near_residue(rspec, radius);
         v.push_back(rspec);
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = refine_residues_with_alt_conf(imol, v, alt_conf);
      }
   }
}

void rsr_sphere_refine() {

   float radius = 4.3;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      auto atom_spec = active_atom.second.second;
      mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
      if (at) {
         std::string alt_conf = at->altLoc;
         coot::residue_spec_t rspec(atom_spec);
         std::vector<coot::residue_spec_t> v = graphics_info_t::molecules[imol].residues_near_residue(rspec, radius);
         v.push_back(rspec);
         g.residue_type_selection_was_user_picked_residue_range = false;
         coot::refinement_results_t rr = refine_residues_with_alt_conf(imol, v, alt_conf);
      }
   }
}

