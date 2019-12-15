
#include "strand-fragments.hh"
#include "coot-coord-utils.hh"

std::vector<std::vector<clipper::Coord_orth> >
coot::mol_to_5_residue_strand_fragments(mmdb::Manager *mol) {

   std::vector<std::vector<clipper::Coord_orth> > v;
   std::cout << "               Sheet info: " << std::endl;
   std::cout << "------------------------------------------------\n";

   mmdb::Model *model_p = mol->GetModel(1);
   if (!model_p) return v;

   int nsheet = model_p->GetNumberOfSheets();

   for (int is=1; is<=nsheet; is++) {
      mmdb::Sheet *sheet_p = model_p->GetSheet(is);

      int nstrand = sheet_p->nStrands;
      for (int istrand=0; istrand<nstrand; istrand++) {
         mmdb::Strand *strand_p = sheet_p->strand[istrand];
         if (strand_p) {
            std::cout << strand_p->sheetID << " " << strand_p->strandNo << " "
                      << strand_p->initChainID << " " << strand_p->initSeqNum
                      << " " << strand_p->endChainID << " " << strand_p->endSeqNum
                      << std::endl;

            std::vector<clipper::Coord_orth> frag;	    
            for (int ires=strand_p->initSeqNum; ires<=(strand_p->endSeqNum-4); ires++) {
               for (int i=0; i<5; i++) {
                  if ((ires+i) > strand_p->endSeqNum) break;
                  mmdb::Residue *residue_p = util::get_residue(strand_p->initChainID, ires+i, "", mol);
                  if (residue_p) {
                     // std::cout << "   found residue " << residue_spec_t(residue_p) << std::endl;
                     mmdb::Atom *ca_at = residue_p->GetAtom(" CA ", 0, "");
                     mmdb::Atom *n_at  = residue_p->GetAtom(" N  ", 0, "");
                     mmdb::Atom *o_at  = residue_p->GetAtom(" O  ", 0, "");
                     mmdb::Atom *c_at  = residue_p->GetAtom(" C  ", 0, "");
                     if (ca_at && n_at && n_at && c_at) {
                        // std::cout << "   found atom " << atom_spec_t(at) << std::endl;
                        frag.push_back(co(ca_at));
                     }
                  }
               }
	            if (frag.size() == 5)
                  v.push_back(frag);
            }
         }
      }
   }
   std::cout << "------------------------------------------------\n";


   return v;
}
