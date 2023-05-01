
#include "simple-restraint.hh"
#include "geometry/mol-utils-2.hh"
#include "analysis/stats.hh"
#include "model-bond-deltas.hh"

coot::model_bond_deltas::model_bond_deltas(mmdb::Manager *mol_in, int imol_in,
					   protein_geometry *geom_p_in) {

   imol = imol_in;
   mol = mol_in;
   geom_p = geom_p_in;

}

void
coot::model_bond_deltas::resolve() {

   if (mol) {

      double mean = -1;
      std::vector<model_bond_deltas> resultants;
	 
      std::map<std::string, std::pair<int, int> > limits = get_residue_number_limits(mol);
      std::map<std::string, std::pair<int, int> >::const_iterator it;

      for (it=limits.begin(); it!=limits.end(); it++) {
	 short int have_flanking_residue_at_start = 0;
	 short int have_flanking_residue_at_end   = 0;
	 short int have_disulfide_residues = 0;
	 int istart_res = it->second.first;
	 int iend_res   = it->second.second;
	 bool do_link_restraints = false;
	 bool do_flank_restraints = false;
	 const std::string &chain_id = it->first;
	 std::string altloc;
	 std::vector<coot::atom_spec_t> fixed_atom_specs;
	 clipper::Xmap<float> dummy_xmap;

	 // Old - and n_atoms_limit_for_nbc is not set
// 	 restraints_container_t restraints(istart_res,
// 					   iend_res,
// 					   have_flanking_residue_at_start,
// 					   have_flanking_residue_at_end,
// 					   have_disulfide_residues,
// 					   altloc,
// 					   chain_id,
// 					   mol,
// 					   fixed_atom_specs, &dummy_xmap);

	 std::vector<std::pair<bool,mmdb::Residue *> > residues;
	 // now fill residues
	 int imod = 1;
	 mmdb::Model *model_p = mol->GetModel(imod);
	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
	       std::string this_chain_id = chain_p->GetChainID();
	       if (this_chain_id == chain_id) {
		  int nres = chain_p->GetNumberOfResidues();
		  residues.reserve(nres);
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     if (residue_p) {
			std::pair<bool, mmdb::Residue *> p(false, residue_p);
			residues.push_back(p);
		     }
		  }
	       }
	    }
	 }

	 restraints_container_t restraints(residues, *geom_p, mol, &dummy_xmap);

	 int n_threads = 2;
	 ctpl::thread_pool thread_pool(n_threads);
	 restraints.thread_pool(&thread_pool, n_threads);

	 restraint_usage_Flags flags = coot::BONDS;
	 pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	 bool do_trans_peptide_restraints = false;
	 restraints.make_restraints(imol, *geom_p, flags, 1, do_trans_peptide_restraints,
				    0.0, 0, false, false, false, pseudos,
				    do_link_restraints,
				    do_flank_restraints);
	 model_bond_deltas resultant = restraints.resolve_bonds();
	 if (resultant.size() > 0) {
	    std::cout << "INFO:: resultant for chain " << chain_id << " n_bonds: "
		      << resultant.size() << std::endl;
	    resultants.push_back(resultant);
	 }
      }

      // loop over chains
      std::vector<double> data;
      for (unsigned int ich=0; ich<resultants.size(); ich++) {
	 // const std::vector<double> &d = resultants[ich];
	 // data.insert(data.end(), d.begin(), d.end());
	 const model_bond_deltas &mbd = resultants[ich];
	 for (std::size_t j=0; j<mbd.size(); j++) {
	    data.push_back(mbd.xyzd.deltas[j]);
	 }
      }
      coot::stats::single s(data);
      mean = s.mean();
      double sd   = sqrt(s.variance());

      // now filter
      double lim_1 = mean - 3.0 * sd;
      double lim_2 = mean + 3.0 * sd;
      std::vector<double> filtered_data;
      for (unsigned int i=0; i<data.size(); i++) {
	 const double &d = data[i];
	 if (d > lim_1 && d < lim_2) {
	    filtered_data.push_back(d);
	 }
      }
      coot::stats::single sf(filtered_data);
      mean = sf.mean();
      sd = sqrt(sf.variance());

      std::cout << "Mean Bond length deltas: " << mean << " +/- " << sd << std::endl;
   }
}
