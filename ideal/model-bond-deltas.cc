
#include "geometry/mol-utils-2.hh"
#include "analysis/stats.hh"
#include "model-bond-deltas.hh"
#include "simple-restraint.hh"

coot::model_bond_deltas::model_bond_deltas(mmdb::Manager *mol_in, int imol_in,
					   protein_geometry *geom_p_in) {

   imol = imol_in;
   mol = mol_in;
   geom_p = geom_p_in;

}

void
coot::model_bond_deltas::resolve() {

   if (mol) {

      std::vector<double> means(3, 0.0);
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
	 restraints_container_t restraints(istart_res,
					   iend_res,
					   have_flanking_residue_at_start,
					   have_flanking_residue_at_end,
					   have_disulfide_residues,
					   altloc,
					   chain_id,
					   mol,
					   fixed_atom_specs, &dummy_xmap);
 

	 restraint_usage_Flags flags = coot::BONDS;
	 pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	 bool do_trans_peptide_restraints = false;
	 restraints.make_restraints(imol, *geom_p, flags, 1, do_trans_peptide_restraints,
				    0.0, 0, false, false, pseudos,
				    do_link_restraints,
				    do_flank_restraints);
	 model_bond_deltas resultant = restraints.resolve_bonds();
	 if (resultant.size() > 0) {
	    std::cout << "INFO:: resultant for chain " << chain_id << " n_bonds: "
		      << resultant.size() << std::endl;
	    resultants.push_back(resultant);
	 }
      }

      // loop over x y z
      for (unsigned int j=0; j<3; j++) {

	 // loop over chains
	 std::vector<double> data;
	 for (unsigned int ich=0; ich<resultants.size(); ich++) {
	    const std::vector<double> &d = resultants[ich].xyz.x[j];
	    data.insert(data.end(), d.begin(), d.end());
	 }
	 coot::stats::single s(data);
	 double mean = s.mean();
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
	 means[j] = sf.mean();
      }
      std::cout << "means: x: " << means[0] << " y: " << means[1] << " z: " << means[2]
		<< std::endl;

   }
}
