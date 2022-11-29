
#include "ramachandran-validation.hh"

std::vector<coot::phi_psi_prob_t>
coot::ramachandran_validation(mmdb::Manager *mol, const ramachandrans_container_t &rc) {

   auto have_close_peptide_bond = [] (mmdb::Residue *residue_1, mmdb::Residue *residue_2) {

      bool status = false;
      mmdb::Atom **residue_atoms_1 = 0;
      int n_residue_atoms_1 = 0;
      // I should iterate over all alt confs
      residue_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
      for (int iat=0; iat<n_residue_atoms_1; iat++) {
         mmdb::Atom *at_1 = residue_atoms_1[iat];
         if (! at_1->isTer()) {
            std::string atom_name_1(at_1->GetAtomName());
            if (atom_name_1 == " C  ") {
               coot::Cartesian pt_c(at_1->x, at_1->y, at_1->z);
               mmdb::Atom **residue_atoms_2 = 0;
               int n_residue_atoms_2 = 0;
               // I should iterate over all alt confs
               residue_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
               for (int iat_inner=0; iat_inner<n_residue_atoms_2; iat_inner++) {
                  mmdb::Atom *at_2 = residue_atoms_2[iat_inner];
                  if (! at_2->isTer()) {
                     std::string atom_name_2(at_2->GetAtomName());
                     if (atom_name_2 == " N  ") {
                        coot::Cartesian pt_n(at_1->x, at_1->y, at_1->z);
                        double dd = coot::Cartesian::lengthsq(pt_c, pt_n);
                        double d = std::sqrt(dd);
                        if (d < 3.0) {
                           status = true;
                           break;
                        }
                     }
                  }
               }
            }
         }
      }
      return status;
   };

   auto get_HA_unit_vector = [] (mmdb::Residue *r) {
      bool status = false;
      coot::Cartesian dir;
      mmdb::Atom *CA = r->GetAtom(" CA ");
      mmdb::Atom *C  = r->GetAtom(" C  ");
      mmdb::Atom *N  = r->GetAtom(" N  ");
      mmdb::Atom *CB = r->GetAtom(" CB ");

      if (CA && C && N && CB) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian dir_3 = ca_pos - cb_pos;
         coot::Cartesian r = dir_1 + dir_2 + dir_3;
         dir = r.unit();
         status = true;
      } else {
         if (CA && C && N) {
            coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
            coot::Cartesian  c_pos( C->x,  C->y,  C->z);
            coot::Cartesian  n_pos( N->x,  N->y,  N->z);
            coot::Cartesian dir_1 = ca_pos - c_pos;
            coot::Cartesian dir_2 = ca_pos - n_pos;
            coot::Cartesian r = dir_1 + dir_2;
            dir = r.unit();
            status = true;
         }
      }
      return std::make_pair(status, dir);
   };

   std::vector<coot::phi_psi_prob_t> v;

   float rama_ball_pos_offset_scale = 0.6;

   rama_plot::phi_psis_for_model_t ppm(mol);
   // This: std::map<coot::residue_spec_t, phi_psi_t> ppm.phi_psi  is now filled

   std::map<coot::residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
   for (it=ppm.phi_psi.begin(); it!=ppm.phi_psi.end(); ++it) {
      const auto &phi_psi(it->second);
      mmdb::Residue *rp = phi_psi.residue_prev;
      mmdb::Residue *rt = phi_psi.residue_this;
      mmdb::Residue *rn = phi_psi.residue_next;
      if (rp && rt && rn) {
         if (have_close_peptide_bond(rp, rt)) {
            if (have_close_peptide_bond(rt, rn)) {
               mmdb::Atom *at = rt->GetAtom(" CA "); // 20221006-PE alt-confs another day
               if (at) {
                  coot::Cartesian pos(at->x, at->y, at->z);
                  std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(rt);
                  coot::Cartesian offset(0,0,rama_ball_pos_offset_scale);
                  if (hav.first) offset = hav.second * rama_ball_pos_offset_scale;
                  coot::util::phi_psi_t cupp(rp, rt, rn);
                  coot::phi_psi_prob_t ppp(cupp, pos + offset, rc);
                  v.push_back(ppp);
               }
            }
         }
      }
   }
   return v;
}

