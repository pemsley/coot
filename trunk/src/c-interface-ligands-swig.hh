

#ifndef C_INTERFACE_LIGANDS_SWIG_HH
#define C_INTERFACE_LIGANDS_SWIG_HH

#include "probe-clash-score.hh"

// We don't need to SWIG this one...
std::pair<CResidue *, int>
new_molecule_sans_biggest_ligand(int imol);


#ifdef USE_GUILE
// Create a new molecule, which is a copy of this molecule without the
// biggest hetgroup.
// 
// return a list of the new molecule number and the spec of the removed residue 
// (or scheme false).
// 
SCM new_molecule_sans_biggest_ligand_scm(int imol);
#endif


coot::probe_clash_score_t
probe_clash_score(const std::string &dots_file_name);


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
#endif

#ifdef USE_GUILE
SCM kullback_liebler_scm(SCM l1, SCM l2);
#endif

#endif // C_INTERFACE_LIGANDS_SWIG_HH
