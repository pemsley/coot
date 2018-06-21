
#include "get-residue.hh"

// Return NULL on residue not found in this molecule. Only look in
// MODEL 1.
// 
mmdb::Residue *
coot::get_residue(const coot::residue_spec_t &res_spec, mmdb::Manager *mol) {

   mmdb::Residue *res = nullptr;
   bool found_res = false;

   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
	 mmdb::Chain *chain_p;
	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    chain_p = model_p->GetChain(i_chain);
	    std::string mol_chain(chain_p->GetChainID());
	    if (mol_chain == res_spec.chain_id) {
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *residue_p;
	       for (int ires=0; ires<nres; ires++) {
		  residue_p = chain_p->GetResidue(ires);
		  if (residue_p->GetSeqNum() == res_spec.res_no) {
		     std::string ins_code(residue_p->GetInsCode());
		     if (ins_code == res_spec.ins_code) {
			res = residue_p;
			found_res = true;
			break;
		     }
		  }
		  if (found_res) break;
	       }
	    }
	    if (found_res) break;
	 }
      }
   }
   return res;
}

clipper::Coord_orth
coot::co(mmdb::Atom *at) {
   return clipper::Coord_orth(at->x, at->y, at->z);
}
