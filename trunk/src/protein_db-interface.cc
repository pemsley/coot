
#include "protein_db-interface.hh"
#include "coot-coord-utils.hh"
#include "mini-mol-utils.hh"

CMMDBManager *make_mol(const std::vector<ProteinDB::Chain> &chains) {

   CMMDBManager *mol = new CMMDBManager;
   std::vector<CResidue *> needs_cb_and_o;

   for (unsigned int ich=0; ich<chains.size(); ich++) { 
      CModel *model_p = new CModel;
      CChain *chain_p = new CChain;
      mol->AddModel(model_p);
      model_p->AddChain(chain_p);
      chain_p->SetChainID("Z");
      for (unsigned int ires=0; ires<chains[ich].size(); ires++) { 
	 if (chains[ich][ires].flag() == ProteinDB::Residue::NONE) {
	    // do nothing
	 } else {
	    // we have a CA at least, and possibly a N and C too.
	    CResidue *residue_p = new CResidue;
	    chain_p->AddResidue(residue_p);
	    residue_p->seqNum = ires+1;
	    CAtom *at_p = new CAtom;
	    residue_p->AddAtom(at_p);
	    residue_p->SetResName("UNK");
	    clipper::Coord_orth ca_pos = chains[ich][ires].coord_ca();
	    at_p->SetCoordinates(ca_pos.x(), ca_pos.y(), ca_pos.z(), 1.0, 30.0);
	    at_p->SetElementName(" C");
	    at_p->SetAtomName(" CA ");
	    if (chains[ich][ires].flag() == ProteinDB::Residue::NORMAL) {
	       clipper::Coord_orth n_pos = chains[ich][ires].coord_n();
	       clipper::Coord_orth c_pos = chains[ich][ires].coord_c();
	       CAtom *at_p_1 = new CAtom;
	       CAtom *at_p_2 = new CAtom;
	       residue_p->AddAtom(at_p_1);
	       residue_p->AddAtom(at_p_2);
	       at_p_1->SetCoordinates(n_pos.x(), n_pos.y(), n_pos.z(), 1.0, 30.0);
	       at_p_2->SetCoordinates(c_pos.x(), c_pos.y(), c_pos.z(), 1.0, 30.0);
	       at_p_1->SetElementName(" N");
	       at_p_2->SetElementName(" C");
	       at_p_1->SetAtomName(" N  ");
	       at_p_2->SetAtomName(" C  ");
	       needs_cb_and_o.push_back(residue_p);
	    } 
	 } 
      }
   }
   mol->FinishStructEdit();

   if (needs_cb_and_o.size()) {

      for (unsigned int ires=0; ires<needs_cb_and_o.size(); ires++) { 

	 // We can position the CB of all residues
	 //
	 // We can position the O of all except the last one.
	 // 
	 CResidue *residue_p = needs_cb_and_o[ires];

	 coot::minimol::residue mini_res_this(residue_p);
	 std::pair<bool, clipper::Coord_orth> cb = coot::cbeta_position(mini_res_this);
	 if (! cb.first) {
	    std::cout << "failed to get CB pos " << std::endl;
	 } else { 
	    CAtom *cb_p = new CAtom;
	    residue_p->AddAtom(cb_p);
	    cb_p->SetElementName(" C");
	    cb_p->SetAtomName(" CB ");
	    cb_p->SetCoordinates(cb.second.x(), cb.second.y(), cb.second.z(), 1.0, 30.0);
	 }

	 if (ires<(needs_cb_and_o.size()-1)) {
	    CResidue *residue_p_next=NULL;
	    coot::residue_spec_t next_res_spec(coot::residue_spec_t(residue_p).next());

	    // is residue_p_next in the vector of residues needs_cb_and_o?
	    for (unsigned int iloop=0; iloop<needs_cb_and_o.size(); iloop++) {
	       // pointer comparison
	       if (needs_cb_and_o[iloop]->chain == residue_p->chain) { 
		  if (coot::residue_spec_t(needs_cb_and_o[iloop]) == next_res_spec) {
		     residue_p_next = needs_cb_and_o[iloop];
		     break;
		  }
	       }
	    }

	    if (!residue_p_next) {
	       // std::cout << "   spec for next residue of " << coot::residue_spec_t(residue_p)
	       // << " was not found " << std::endl;
	    } else { 
	       coot::minimol::residue mini_res_next(residue_p_next);
	       std::pair<bool, clipper::Coord_orth> o = coot::o_position(mini_res_this, mini_res_next);
	       if (! o.first) {
		  std::cout << "   failed to get O pos " << std::endl;
	       } else { 
		  CAtom *o_p = new CAtom;
		  residue_p->AddAtom(o_p);
		  o_p->SetElementName(" O");
		  o_p->SetAtomName(" O  ");
		  o_p->SetCoordinates(o.second.x(), o.second.y(), o.second.z(), 1.0, 30.0);
	       }
	    }
	 }
      }
   } 
   mol->FinishStructEdit();
   return mol;
} 
