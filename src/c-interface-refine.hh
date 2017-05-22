
void    add_initial_position_restraints(int imol, const std::vector<coot::residue_spec_t> &residue_specs, double weight);
// this removed all initial position restraints, not just those listed.
void remove_initial_position_restraints(int imol, const std::vector<coot::residue_spec_t> &residue_specs);

void use_unimodal_ring_torsion_restraints(const std::string &res_name);

void set_refinement_geman_mcclure_alpha(float alpha);
