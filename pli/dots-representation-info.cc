
#include "coot-utils/coot-coord-utils.hh"
#include "dots-representation-info.hh"

// create (and later delete, of course) a new molecule by deep copying
// and assembling the passed residues.  Use that to make an atom
// selection which gets passed to
// dots_representation_info_t::solvent_exposure().  Notice that we
// pass back atom specs.
//
std::vector<std::pair<coot::atom_spec_t, float> >
pli::dots_representation_info_t::solvent_accessibilities(mmdb::Residue *res_ref,
                                                         const std::vector<mmdb::Residue *> &near_residues) {

   std::vector<std::pair<coot::atom_spec_t, float> > v;

   // no solvent exposures if the ligand does not have *any*
   // neighbours (i.e. it's floating in space (probably the only
   // residue in the molecule).
   if (near_residues.size() == 0)
      return v;

   std::vector<mmdb::Residue *> residues = near_residues;
   residues.push_back(res_ref);

   std::pair<bool, mmdb::Manager *> mol =
      coot::util::create_mmdbmanager_from_residue_vector(residues, 0);

   if (mol.first) {

      int SelHnd = mol.second->NewSelection();
      mol.second->SelectAtoms(SelHnd, 0, res_ref->GetChainID(),
                              res_ref->GetSeqNum(), res_ref->GetInsCode(),
                              res_ref->GetSeqNum(), res_ref->GetInsCode(),
                              "*", "*", "*", "*");

      std::vector<std::pair<mmdb::Atom *, float> > se = solvent_exposure(SelHnd, mol.second);
      v.resize(se.size());
      for (unsigned int i=0; i<se.size(); i++) {
         v[i] = std::pair<coot::atom_spec_t, float> (coot::atom_spec_t(se[i].first), se[i].second);
      }

      mol.second->DeleteSelection(SelHnd);
      delete mol.second;
   }
   return v;
}


std::vector<pli::solvent_exposure_difference_helper_t>
pli::dots_representation_info_t::solvent_exposure_differences(mmdb::Residue *res_ref,
                                                              const std::vector<mmdb::Residue *> &near_residues) const {

   std::vector<pli::solvent_exposure_difference_helper_t> v;

   std::vector<mmdb::Residue *> residues = near_residues;
   residues.push_back(res_ref);

   std::pair<bool, mmdb::Manager *> mol_holo =
      coot::util::create_mmdbmanager_from_residue_vector(residues, 0);

   std::pair<bool, mmdb::Manager *> mol_apo =
      coot::util::create_mmdbmanager_from_residue_vector(near_residues, 0);

   if (mol_holo.first) {
      if (mol_apo.first) {

         for (unsigned int ir=0; ir<near_residues.size(); ir++) {
            std::string res_name = near_residues[ir]->GetResName();
            if (res_name != "HOH") {
               int SelHnd_holo = mol_holo.second->NewSelection();
               int SelHnd_apo  =  mol_apo.second->NewSelection();
               mol_holo.second->SelectAtoms(SelHnd_holo, 0, near_residues[ir]->GetChainID(),
                                            near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                            near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                            "*", "*", "!H", "*");
               mol_apo.second->SelectAtoms(SelHnd_apo, 0, near_residues[ir]->GetChainID(),
                                           near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                           near_residues[ir]->GetSeqNum(), near_residues[ir]->GetInsCode(),
                                           "*", "*", "!H", "*");
               std::vector<std::pair<mmdb::Atom *, float> > se_holo = solvent_exposure(SelHnd_holo, mol_holo.second);
               std::vector<std::pair<mmdb::Atom *, float> > se_apo  = solvent_exposure(SelHnd_apo,   mol_apo.second);

               double se_frac_holo = 0.0;
               double se_frac_apo  = 0.0;

               for (unsigned int iah=0; iah<se_holo.size(); iah++) {
                  se_frac_holo += se_holo[iah].second;
               }
               for (unsigned int iaa=0; iaa<se_apo.size(); iaa++) {
                  se_frac_apo += se_apo[iaa].second;
               }

               if (0)
                  std::cout << "storing " << coot::residue_spec_t(near_residues[ir]) << " "
                            << near_residues[ir]->GetResName() << " "
                            << se_frac_holo << " " << se_frac_apo << std::endl;
               coot::residue_spec_t res_spec(near_residues[ir]);
               pli::solvent_exposure_difference_helper_t sed(res_spec, se_frac_holo, se_frac_apo);
               v.push_back(sed);

               mol_holo.second->DeleteSelection(SelHnd_holo);
               mol_apo.second->DeleteSelection(SelHnd_apo);
            }
         }

         delete mol_apo.second;
      }
      delete mol_holo.second;
   }
   return v;
}



// simply transfer the atoms of mol to the points vector
//
void
pli::dots_representation_info_t::pure_points(mmdb::Manager *mol) {

   is_closed = 0;
   int imod = 1;
   std::vector<clipper::Coord_orth> local_points;

   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::Residue *residue_p;
      mmdb::Atom *at;
      for (int ires=0; ires<nres; ires++) {
         residue_p = chain_p->GetResidue(ires);
         int n_atoms = residue_p->GetNumberOfAtoms();

         for (int iat=0; iat<n_atoms; iat++) {
            at = residue_p->GetAtom(iat);
            local_points.push_back(clipper::Coord_orth(at->x, at->y, at->z));
         }
      }
   }
   coot::colour_holder col(0.3, 0.4, 0.5);
   std::pair<coot::colour_holder, std::vector<clipper::Coord_orth> > p(col, local_points);
   points.push_back(p);
}
