

#ifndef C_INTERFACE_LIGANDS_SWIG_HH
#define C_INTERFACE_LIGANDS_SWIG_HH

#include "probe-clash-score.hh"
#include "ligand-check.hh"

/*! \file
  \brief Coot Scripting Interface - Ligands interface
*/

// We don't need to SWIG this one...
std::pair<mmdb::Residue *, int>
new_molecule_sans_biggest_ligand(int imol);


// return a new molecule number
int get_monomer_molecule_by_network_and_dict_gen(const std::string &text);


#ifdef USE_GUILE
// Create a new molecule, which is a copy of this molecule without the
// biggest hetgroup.
// 
// return a list of the new molecule number and the spec of the removed residue 
// (or scheme false).
// 
SCM new_molecule_sans_biggest_ligand_scm(int imol);
void gui_ligand_metrics_scm(SCM ligand_spec, SCM ligand_metrics, double percentile_limit);

#endif

#ifdef USE_PYTHON
PyObject *new_molecule_sans_biggest_ligand_py(int imol);
void gui_ligand_metrics_py(PyObject *ligand_spec, PyObject *ligand_metrics, double percentile_limit);
#endif

#ifdef USE_PYTHON
// this is not a ligands function (although it can be used for ligands)
// it doesn't belong here
PyObject *residues_distortions_py(int imol, PyObject *residue_spec_list);
PyObject *get_intermediate_atoms_distortions_py();
#endif

// This don't call graphics_draw(), so the caller needs to do so.
//
coot::probe_clash_score_t
probe_clash_score(const std::string &dots_file_name);

#ifdef USE_GUILE
// internal bumps scoring, sphere overlap
SCM ligand_atom_overlaps_scm(int imol, SCM ligand_spec, double neighb_radius);
#endif

#ifdef USE_PYTHON
// internal bumps scoring, sphere overlap
PyObject *ligand_atom_overlaps_py(int imol, PyObject *ligand_spec, double neighb_radius);
#endif



#ifdef USE_GUILE
bool
residues_torsions_match_scm(int imol_1, SCM res_1,
			    int imol_2, SCM res_2,
			    float tolerance); // in degrees
#endif // USE_GUILE

#ifdef USE_PYTHON
bool
residues_torsions_match_py(int imol_1, PyObject *res_1,
			   int imol_2, PyObject *res_2,
			   float tolerance); // in degrees
#endif // USE_PYTHON

#ifdef USE_GUILE
double kolmogorov_smirnov_scm(SCM l1, SCM l2);
double kolmogorov_smirnov_vs_normal_scm(SCM l1, double mean, double std_dev);
#endif

#ifdef USE_GUILE
SCM kullback_liebler_scm(SCM l1, SCM l2);
#endif

#ifdef USE_PYTHON
double kolmogorov_smirnov_py(PyObject *l1, PyObject *l2);
double kolmogorov_smirnov_vs_normal_py(PyObject *l1, double mean, double std_dev);
#endif

#ifdef USE_PYTHON
PyObject *kullback_liebler_py(PyObject *l1, PyObject *l2);
#endif

// Returning void ATM.  We shoud return an interesting object at some
// stage. Perhaps a coot::geometry_distortion_info_container_t?
//
double
print_residue_distortions(int imol, std::string chain_id, int res_no, std::string ins_code);
void
display_residue_distortions(int imol, std::string chain_id, int res_no, std::string ins_code);

void display_residue_hydrogen_bond_atom_status_using_dictionary(int imol, std::string chain_id, int res_no,
								std::string ins_code);
void
write_dictionary_from_residue(int imol, std::string chain_id, int res_no, std::string ins_code, std::string cif_file_name);
void 
add_dictionary_from_residue(int imol, std::string chain_id, int res_no, std::string ins_code);

void
invert_chiral_centre(int imol, std::string chain_id, int res_no, std::string ins_code, std::string atom_name);

// Read cif_dict_in, match the atoms there-in to those of the dictionary of reference_comp_id.
// Write a new dictionary to cif_dict_out.
// 
// Return 1 if was successful in doing the atom matching and writing the cif file.
//
// e.g: match_residue_and_dictionary(0, "A", 1, "", "DRG-pyrogen.cif", "DRG-renamed.cif", "DRG", "LYS")
int
match_residue_and_dictionary(int imol, std::string chain_id, int res_no, std::string ins_code,
			     std::string cif_dict_in,
			     std::string cif_dict_out,
			     std::string cif_dict_comp_id, // comp-id in the input file
			     std::string reference_comp_id,
			     std::string output_comp_id,
			     std::string output_compound_name);

// This is the GUI interface: User sits over on top of the residue
// they want to change and provides an output dictionary cif file
// name, a reference-comp-id and a compd-id for the new residue (type).
int
match_this_residue_and_dictionary(int imol, std::string chain_id, int res_no, std::string ins_code,
				  std::string cif_dict_out,
				  std::string reference_comp_id,
				  std::string output_comp_id);
				  
// return False if unknown
bool comprised_of_organic_set_p(const std::string &rn);

// all-atom contact dots.  This is not the place for this declaration (not a ligand function)
//
//! \brief calculate all-atom contact dots
//!
//! remove contact dots objects using the Generic Display Objects dialog
void coot_all_atom_contact_dots(int imol);

#ifdef USE_PYTHON
void coot_contact_dots_for_ligand_py(int imol, PyObject *ligand_spec);
// change HE2 to HD1 and vice versa
void switch_HIS_protonation_py(int imol, PyObject *residue_spec);
#endif

// this is not a ligand function - it does not belong here.
void coot_reduce(int imol);


#ifdef USE_GUILE
void coot_contact_dots_for_ligand_scm(int imol, SCM residue_spec_scm);
// change HE2 to HD1 and vice versa
void switch_HIS_protonation_scm(int imol, SCM residue_spec_scm);
#endif



// we want to read in the built-in database to convert these scores to percentiles
// return -1 (test for negative) on failure
// metric_name examples: density_correlation coot_diff_map_KS_2 mogul_z_worst bumps_1
// range: 0->100
// 
// reverse order for mogul and bumps (for example) because low-is-good
// 
double get_ligand_percentile(std::string metric_name, double metric_value, short int reverse_order);


// is enhanced ligand version
bool enhanced_ligand_coot_p();


#endif // C_INTERFACE_LIGANDS_SWIG_HH
