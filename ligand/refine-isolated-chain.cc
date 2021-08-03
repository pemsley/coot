
#include "refine-isolated-chain.hh"
#include "ideal/simple-restraint.hh"

void coot::refine_isolated_chain(mmdb::Chain *chain_p, mmdb::Manager *mol_for_this_chain, const coot::protein_geometry &geom,
                                 ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight,
                                 const clipper::Xmap<float> &xmap) {
   std::vector<std::pair<bool, mmdb::Residue *> > residues;
   int n_res = chain_p->GetNumberOfResidues();
   unsigned int n_atoms_for_refinement = 0; // debug crash
   for (int ires=0; ires<n_res; ires++) {
      mmdb::Residue *residue_p = chain_p->GetResidue(ires);
      if (residue_p) {
         residues.push_back(std::make_pair(false, residue_p));

         if (true) {
            int n_atoms = residue_p->GetNumberOfAtoms();
            for (int iat=0; iat<n_atoms; iat++) {
               mmdb::Atom *at = residue_p->GetAtom(iat);
               if (! at->isTer()) {
                  n_atoms_for_refinement++;
               }
            }
         }
      }
   }

   if (true) {
      std::string chain_id(chain_p->GetChainID());
      std::cout << "in refine_isolated_chain(): chain " << chain_id
                << " n_atoms_for_refinement: " << n_atoms_for_refinement << std::endl;
   }

   std::vector<mmdb::Link> links;
   std::vector<atom_spec_t> fixed_atom_specs;
   coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;
   coot::restraints_container_t restraints(residues, links, geom, mol_for_this_chain, fixed_atom_specs, &xmap);
   restraints.thread_pool(thread_pool_p, n_threads);
   restraints.set_quiet_reporting();
   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
   bool do_internal_torsions = false;
   restraints.add_map(weight);
   int imol = 0;
   restraints.make_restraints(imol, geom, flags, do_internal_torsions, false, 0, 0, true, true, false, pseudos);
   restraints.minimize(flags);
}

