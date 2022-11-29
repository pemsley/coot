#include "residue-validation-information.hh"

coot::residue_validation_information_t::residue_validation_information_t(
    const coot::residue_spec_t &rs,
    const coot::atom_spec_t &atom_spec_in,
    double distortion_in, const std::string &l) 
    :residue_spec(rs), 
    atom_spec(atom_spec_in), 
    distortion(distortion_in), 
    block_colour("#904040"),
    label(l){
        
}
