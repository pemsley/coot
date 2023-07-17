
#include "coot-phi-psi.hh"

void
coot::phi_psis_for_model_t::process(mmdb::Manager *mol, int imod_in) {

   coot::phi_psis_for_model_t ppm(imod_in);
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      if (imod == imod_in) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int nchains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<nchains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               if (nres > 2) { 
                  for (int ires=1; ires<(nres-1); ires++) { 
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     
                     // this could be improved - using serial numbers.
                     mmdb::Residue *res_prev = chain_p->GetResidue(ires-1);
                     mmdb::Residue *res_next = chain_p->GetResidue(ires+1);

                     if (res_prev && residue_p && res_next) {
                        try {
                           // coot::phi_psi_t constructor can throw an error
                           // (e.g. bonding atoms too far apart).
                           coot::residue_spec_t spec(residue_p);
                           coot::util::phi_psi_with_residues_t pp(res_prev, residue_p, res_next);
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

}
