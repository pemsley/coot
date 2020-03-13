
#include "strand-fragments.hh"
#include "coot-coord-utils.hh"

// now return the mainchain fragments after the CA coordinates
std::pair<std::vector<std::vector<clipper::Coord_orth> >, std::vector<std::vector<clipper::Coord_orth> > >
coot::mol_to_5_residue_strand_fragments(mmdb::Manager *mol) {

   std::vector<std::vector<clipper::Coord_orth> > v;
   std::vector<std::vector<clipper::Coord_orth> > mc_v; // mainchain fragments

   std::cout << "               Sheet info: " << std::endl;
   std::cout << "------------------------------------------------\n";

   mmdb::Model *model_p = mol->GetModel(1);
   if (!model_p) return std::pair<std::vector<std::vector<clipper::Coord_orth> >, std::vector<std::vector<clipper::Coord_orth> > >(v,v);

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

            std::vector<clipper::Coord_orth> ca_frag;
            std::vector<clipper::Coord_orth> mc_frag; // atom order: N CA C O
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
                     if (ca_at && n_at && c_at && o_at) {
                        // std::cout << "   found atom " << atom_spec_t(at) << std::endl;
                        ca_frag.push_back(co(ca_at));
                        mc_frag.push_back(co(n_at));
                        mc_frag.push_back(co(ca_at));
                        mc_frag.push_back(co(c_at));
                        mc_frag.push_back(co(o_at));
                     }
                  }
               }
	            if (ca_frag.size() == 5) {
                  v.push_back(ca_frag);
                  mc_v.push_back(mc_frag);
               }
            }
         }
      }
   }
   std::cout << "------------------------------------------------\n";
   std::pair<std::vector<std::vector<clipper::Coord_orth> >, std::vector<std::vector<clipper::Coord_orth> > > p(v, mc_v);
   return p;
}
