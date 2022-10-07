
#include "rama-plot-phi-psi.hh"


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

