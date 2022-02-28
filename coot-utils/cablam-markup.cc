
#include <compat/coot-sysdep.h>
#include <fstream>
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

#include "compat/coot-sysdep.h"
#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "cablam-markup.hh"

coot::cablam_markup_t::cablam_markup_t(mmdb::Atom *O_prev_at,
                                       mmdb::Atom *O_this_at,
                                       mmdb::Atom *O_next_at,
                                       mmdb::Atom *CA_prev_at,
                                       mmdb::Atom *CA_this_at,
                                       mmdb::Atom *CA_next_at,
                                       mmdb::Atom *CA_next_next_at) {
   residue = 0;
   if (! O_prev_at)  return;
   if (! O_this_at)  return;
   if (! O_next_at)  return;
   if (! CA_prev_at) return;
   if (! CA_this_at) return;
   if (! CA_next_at) return;
   score = -1;
   residue = O_this_at->residue;
   O_prev_pos = co(O_prev_at);
   O_this_pos = co(O_this_at);
   O_next_pos = co(O_next_at);
   clipper::Coord_orth CA_prev_pos = co(CA_prev_at);
   clipper::Coord_orth CA_this_pos = co(CA_this_at);
   clipper::Coord_orth CA_next_pos = co(CA_next_at);
   clipper::Coord_orth CA_next_next_pos = co(CA_next_next_at);

   // scrap of paper algebra
   clipper::Coord_orth d(O_this_pos - CA_this_pos);
   clipper::Coord_orth r(CA_next_pos - CA_this_pos);
   clipper::Coord_orth r_uv(r.unit());
   double dp = clipper::Coord_orth::dot(r,d);
   double r_len = std::sqrt(r.lengthsq());
   double r1_len = dp/r_len;

   CA_proj_point_this = CA_this_pos + r1_len * r_uv;

   // now for the other pair
   d = clipper::Coord_orth(O_prev_pos - CA_prev_pos);
   r = clipper::Coord_orth(CA_this_pos - CA_prev_pos);
   r_uv = clipper::Coord_orth(r.unit());
   dp = clipper::Coord_orth::dot(r,d);
   r_len = std::sqrt(r.lengthsq());
   r1_len = dp/r_len;

   CA_proj_point_prev = CA_prev_pos + r1_len * r_uv;

   // and the "next" CO which needs the next_next CA:
   d = clipper::Coord_orth(O_next_pos - CA_next_pos);
   r = clipper::Coord_orth(CA_next_next_pos - CA_next_pos);
   r_uv = clipper::Coord_orth(r.unit());
   dp = clipper::Coord_orth::dot(r,d);
   r_len = std::sqrt(r.lengthsq());
   r1_len = dp/r_len;

   CA_proj_point_next = CA_next_pos + r1_len * r_uv;
}


coot::cablam_markup_t
coot::calc_cablam(mmdb::Chain *chain_p, mmdb::Residue *residue_this_p,
                  int ires, double score) {

   cablam_markup_t cm; // null residue in this constructor
   mmdb::Residue *residue_prev_p = chain_p->GetResidue(ires-1);
   mmdb::Residue *residue_next_p = chain_p->GetResidue(ires+1);
   mmdb::Residue *residue_next_next_p = chain_p->GetResidue(ires+2);

   if (residue_prev_p->GetSeqNum() + 1 == residue_this_p->GetSeqNum()) {
   } else {
      return cm; // fail on tandem residues test
   }
   if (residue_next_p->GetSeqNum() - 1 == residue_this_p->GetSeqNum()) {
   } else {
      return cm; // fail on tandem residues test
   }
#if 0
   // I am not yet that I want this test.
   if (residue_next_next_p->GetSeqNum() - 2 == residue_this_p->GetSeqNum()) {
   } else {
      return cm; // fail on tandem residues test
   }
#endif
   mmdb::Atom *O_this = 0;
   mmdb::Atom *O_prev = 0;
   mmdb::Atom *O_next = 0;
   mmdb::Atom *CA_this = 0;
   mmdb::Atom *CA_prev = 0;
   mmdb::Atom *CA_next = 0;
   mmdb::Atom *CA_next_next = 0;
   int n_atoms_prev = residue_prev_p->GetNumberOfAtoms();
   int n_atoms_this = residue_this_p->GetNumberOfAtoms();
   int n_atoms_next = residue_next_p->GetNumberOfAtoms();
   int n_atoms_next_next = residue_next_next_p->GetNumberOfAtoms();
   for (int iat=0; iat<n_atoms_prev; iat++) {
      mmdb::Atom *at = residue_prev_p->GetAtom(iat);
      std::string alt_loc(at->altLoc);
      if (alt_loc.empty()) { // no cablams for altconfed atoms
         std::string atom_name(at->GetAtomName());
         if (atom_name == " O  ") {
            O_prev = at;
            continue;
         }
         if (atom_name == " CA ") {
            CA_prev = at;
            continue;
         }
      }
   }
   for (int iat=0; iat<n_atoms_this; iat++) {
      mmdb::Atom *at = residue_this_p->GetAtom(iat);
      std::string alt_loc(at->altLoc);
      if (alt_loc.empty()) { // no cablams for altconfed atoms
         std::string atom_name(at->GetAtomName());
         if (atom_name == " O  ") {
            O_this = at;
            continue;
         }
         if (atom_name == " CA ") {
            CA_this = at;
            continue;
         }
      }
   }
   for (int iat=0; iat<n_atoms_next; iat++) {
      mmdb::Atom *at = residue_next_p->GetAtom(iat);
      std::string alt_loc(at->altLoc);
      if (alt_loc.empty()) { // no cablams for altconfed atoms
         std::string atom_name(at->GetAtomName());
         if (atom_name == " O  ") {
            O_next = at;
            continue;
         }
         if (atom_name == " CA ") {
            CA_next = at;
            continue;
         }
      }
   }
   for (int iat=0; iat<n_atoms_next_next; iat++) {
      mmdb::Atom *at = residue_next_next_p->GetAtom(iat);
      std::string alt_loc(at->altLoc);
      if (alt_loc.empty()) { // no cablams for altconfed atoms
         std::string atom_name(at->GetAtomName());
         if (atom_name == " CA ") {
            CA_next_next = at;
            continue;
         }
      }
   }
   if (O_prev && O_this && O_next) {
      if (CA_prev && CA_this && CA_next && CA_next_next) {
         clipper::Coord_orth p1 = co(CA_prev);
         clipper::Coord_orth p2 = co(CA_this);
         clipper::Coord_orth p3 = co(CA_next);
         double v1_sqrd = (p2 - p1).lengthsq();
         double v2_sqrd = (p2 - p3).lengthsq();
         double v1 = std::sqrt(v1_sqrd);
         double v2 = std::sqrt(v2_sqrd);
         if (v1 < 3.9 && v2 < 3.9) {
            cablam_markup_t cm_local(O_prev, O_this, O_next, CA_prev, CA_this, CA_next, CA_next_next);
            cm_local.score = score;
            cm = cm_local;
         }
      }
   }
   return cm;
}


std::vector<coot::cablam_markup_t>
coot::make_cablam_markups(const std::vector<std::pair<residue_spec_t, double> > &residues,
                          mmdb::Manager *mol) {

   std::vector<cablam_markup_t> v;
   std::vector<std::pair<residue_spec_t, double> >::const_iterator it;
   for (it=residues.begin(); it!=residues.end(); ++it) {
      const residue_spec_t &cablam_res_spec(it->first);
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            int ires_max = nres - 2;
            for (int ires=1; ires<ires_max; ires++) {
               mmdb::Residue *residue_this_p = chain_p->GetResidue(ires);
               residue_spec_t spec(residue_this_p);
               if (cablam_res_spec == spec) {
                  auto cm = calc_cablam(chain_p, residue_this_p, ires, it->second);
                  v.push_back(cm);
               }
            }
         }
      }
   }
   return v;
}

// parse this log file and call the above function for each cablam outlier residue
std::vector<coot::cablam_markup_t>
coot::make_cablam_markups(mmdb::Manager *mol, const std::string &cablam_output_file_name) {

   std::vector<cablam_markup_t> v;

   std::vector<std::pair<residue_spec_t, double> > scored_baddie_specs;
   if (file_exists(cablam_output_file_name)) {
      std::ifstream f(cablam_output_file_name.c_str());
      if (f) {
         std::cout << "debug:: opened " << cablam_output_file_name << std::endl;
         std::string line;
         while (std::getline(f, line)) {
            // parse by position, not field!
            if (line.length() == 90) {
               std::string chain_id = line.substr(0,2); // need to try other files
               if (chain_id[0] == ' ') chain_id = line.substr(1,1);
               std::string resno_string = line.substr(2,4);
               try {
                  int res_no = util::string_to_int(resno_string);
                  std::string residue_type  = line.substr( 8,3);
                  std::string cablam_string = line.substr(13,6);
                  std::string cablam_type_string = line.substr(20,7);
                  if (false)
                     std::cout << "debug:: " << chain_id << " " << res_no << " "
                              << residue_type << " "  // << level_string << " " << level
                              << " cablam_string \"" << cablam_string << "\""
                              << " type \"" << cablam_type_string << "\"" << std::endl;
                  if (cablam_string == "CaBLAM") {
                     // either Disfavoured or an Outlier
                     if (cablam_type_string == "Outlier") {
                        std::string level_string = line.substr(34,6);
                        double level = util::string_to_float(level_string);

                        residue_spec_t res_spec(chain_id, res_no, "");
                        std::pair<residue_spec_t, double> p(res_spec, level);
                        scored_baddie_specs.push_back(p);
                     }
                  }
               }
               catch (const std::runtime_error &rte) {
                  std::cout << rte.what() << std::endl;
               }
            }
         }
      }
   } else {
      std::cout << "WARNING:: file not found " << cablam_output_file_name << std::endl;
   }
   v = make_cablam_markups(scored_baddie_specs, mol);
   return v;
}


coot::cablam_like_geometry_stats_t::cablam_like_geometry_stats_t(const coot::cablam_markup_t &cm) {

   auto v_prev = cm.O_prev_pos - cm.CA_proj_point_prev;
   auto v_this = cm.O_this_pos - cm.CA_proj_point_this;
   auto v_next = cm.O_next_pos - cm.CA_proj_point_next;

   double dd = (cm.CA_proj_point_next - cm.CA_proj_point_prev).lengthsq();
   double d = std::sqrt(dd);
   double a1 = clipper::Coord_orth::dot(v_prev, v_this);
   double a2 = clipper::Coord_orth::dot(v_next, v_this);

   residue = cm.residue;
   dp_prev_to_mid = a1;
   dp_next_to_mid = a2;
   dist_proj_point_prev_to_next = d;

}

std::vector<coot::cablam_like_geometry_stats_t>
coot::get_cablam_like_geometry_stats(mmdb::Manager *mol) {

   std::vector<cablam_like_geometry_stats_t> v;

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         int n_res_max = n_res - 2;
         for (int ires=1; ires<n_res_max; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               auto cm = calc_cablam(chain_p, residue_p, ires);
               if (cm.residue) {
                  v.push_back(cablam_like_geometry_stats_t(cm));
               }
            }
         }
      }
   }

   return v;
}
