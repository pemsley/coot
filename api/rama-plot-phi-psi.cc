
#include <clipper/core/coords.h>
#include "rama-plot-phi-psi.hh"

rama_plot::phi_psi_t::phi_psi_t(double phi_in, double psi_in, const std::string &res_name_in, const std::string &label_in,
                int res_no_in, const std::string &ins_code_in, const std::string &chain_id_in, bool is_pre_pro_in) : phi(phi_in), psi(psi_in) {
   init();
   label = label_in;
   residue_name = res_name_in;
   res_no = res_no_in;
   ins_code = ins_code_in;
   chain_id = chain_id_in;
   is_pre_pro = is_pre_pro_in;
}

// each residue needs to be non-null (there is no protection
// for that in this function).
std::pair<bool, rama_plot::phi_psi_t>
rama_plot::util::get_phi_psi(mmdb::Residue *residue_0, mmdb::Residue *residue_1, mmdb::Residue *residue_2) {

   bool is_valid_flag = 0;
   bool is_pre_pro = 0;
   rama_plot::phi_psi_t phi_psi; // part of the returned value
   int nResidueAtoms;
   mmdb::PPAtom res_selection;
   int natom = 0;
   int ires = residue_1->GetSeqNum();
   clipper::Coord_orth c_prev, n_this, ca_this, c_this, n_next;

   residue_0->GetAtomTable(res_selection, nResidueAtoms);
   if (nResidueAtoms > 0) {
      for (int j=0; j<nResidueAtoms; j++) {
         std::string atom_name = res_selection[j]->name;
         if (atom_name == " C  ") {
            c_prev = clipper::Coord_orth(res_selection[j]->x,
                                         res_selection[j]->y,
                                         res_selection[j]->z);
            natom++;
         }
      }
   }

   std::string res_name_1(residue_1->GetResName());
   residue_1->GetAtomTable(res_selection, nResidueAtoms);
   if (nResidueAtoms > 0) {
      for (int j=0; j<nResidueAtoms; j++) {
         std::string atom_name = res_selection[j]->name;
         if (atom_name == " C  ") {
            c_this = clipper::Coord_orth(res_selection[j]->x,
                                         res_selection[j]->y,
                                         res_selection[j]->z);
            natom++;
         }
         if (atom_name == " CA ") {
            ca_this = clipper::Coord_orth(res_selection[j]->x,
                                          res_selection[j]->y,
                                          res_selection[j]->z);
            natom++;
         }
         if (atom_name == " N  ") {
            n_this = clipper::Coord_orth(res_selection[j]->x,
                                         res_selection[j]->y,
                                         res_selection[j]->z);
            natom++;
         }
      }
   }

   residue_2->GetAtomTable(res_selection, nResidueAtoms);
   if (std::string(residue_2->GetResName()) == "PRO")
      is_pre_pro = 1;
   if (nResidueAtoms > 0) {
      for (int j=0; j<nResidueAtoms; j++) {
         std::string atom_name = res_selection[j]->name;
         if (atom_name == " N  ") {
            n_next = clipper::Coord_orth(res_selection[j]->x,
                                         res_selection[j]->y,
                                         res_selection[j]->z);
            natom++;
         }
      }
   }

   if (natom == 5) {
      char num[30];
      snprintf(num,20,"%d",ires);
      std::string label(num);
      std::string segid = residue_1->GetChainID();
      std::string inscode = residue_1->GetInsCode();
      label += inscode;
      label += " ";
      label += segid;
      label += " ";
      label += residue_1->name;

      double phi = clipper::Util::rad2d(ca_this.torsion(c_prev, n_this, ca_this, c_this));
      double psi = clipper::Util::rad2d(ca_this.torsion(n_this, ca_this, c_this, n_next));

      phi_psi = rama_plot::phi_psi_t(phi, psi,
                                     residue_1->name,
                                     label.c_str(),
                                     ires,
                                     inscode,
                                     segid,
                                     is_pre_pro);
      // peptide bonding atoms have to be within 2.0A, or this is not
      // a valid peptide.
      //
      double dist_1 = clipper::Coord_orth::length(c_prev, n_this);
      double dist_2 = clipper::Coord_orth::length(c_this, n_next);

      if (dist_1 < 2.0)
         if (dist_2 < 2.0)
            is_valid_flag = 1;

   } else {
      // std::cout << "only found " << natom << " atoms " << std::endl;
   }

   rama_plot::phi_psi_t phi_psi_with_residues(phi_psi);
   phi_psi_with_residues.residue_prev = residue_0;
   phi_psi_with_residues.residue_this = residue_1;
   phi_psi_with_residues.residue_next = residue_2;

   // phi_psi_with_residues.type = clipper::Ramachandran::All;
   if (res_name_1 == "GLY") {
      phi_psi_with_residues.type = clipper::Ramachandran::Gly2;
   } else {
      if (res_name_1 == "PRO") {
         phi_psi_with_residues.type = clipper::Ramachandran::Pro2;
      } else {
         if (res_name_1 == "ILE" || res_name_1 == "VAL") {
            phi_psi_with_residues.type = clipper::Ramachandran::IleVal2;
         } else {
            phi_psi_with_residues.type = clipper::Ramachandran::NoGPIVpreP2;
         }
      }
   }
   if (is_pre_pro) {
      phi_psi_with_residues.type = clipper::Ramachandran::PrePro2;
   }

   return std::pair<bool, rama_plot::phi_psi_t> (is_valid_flag, phi_psi_with_residues);

}

// from ../src/rama-plot.cc
//
// this can throw an exception (e.g. bonding atoms too far
// apart).
rama_plot::phi_psi_t::phi_psi_t(mmdb::Residue *prev_res, mmdb::Residue *this_res, mmdb::Residue *next_res) {

   if (prev_res && this_res && next_res) {
      std::pair<bool, rama_plot::phi_psi_t> bpp = rama_plot::util::get_phi_psi(prev_res, this_res, next_res);

      if (! bpp.first) {
         std::string mess = "bad residues for phi,psi calculation";
         throw std::runtime_error(mess);
      } else {
         *this = bpp.second;
      }
   }
}


// from ../src/rama-plot.cc
//
// fill phi_psi vector
//
void
rama_plot::phi_psis_for_model_t::generate_phi_psis(mmdb::Manager *mol_in) {

   if (! mol_in) return;

   int n_models = mol_in->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {
      mmdb::Model *model_p = mol_in->GetModel(imod);
      if (model_p) {
         mmdb::Chain *chain_p;
         int nchains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<nchains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            if (nres > 2) {
               for (int ires=1; ires<(nres-1); ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);

                  // this could be improved
                  mmdb::Residue *res_prev = chain_p->GetResidue(ires-1);
                  mmdb::Residue *res_next = chain_p->GetResidue(ires+1);

                  if (res_prev && residue_p && res_next) {
                     try {
                        // coot::phi_psi_t constructor can throw an error
                        // (e.g. bonding atoms too far apart).
                        coot::residue_spec_t spec(residue_p);
                        rama_plot::phi_psi_t pp(res_prev, residue_p, res_next);
                        add_phi_psi(spec, pp);
                     }
                     catch (const std::runtime_error &rte) {
                        // nothing too bad, just don't add that residue
                        // to the plot
                     }
                  }
               }
            }
         }
      }
   }
}

