

#ifdef USE_MOLECULES_TO_TRIANGLES
#ifdef USE_PYTHON

// Martin's Triangles

// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation_py(int imol, PyObject *atom_selection_py, PyObject *ColorScheme_py, PyObject *style_py);
#endif // USE_PYTHON

#ifdef USE_GUILE
int add_molecular_representation_scm(int imol, SCM atom_selection_scm, SCM ColorScheme_scm, SCM style_scm);
#endif // USE_GUILE

// not dependent on scm or python
void remove_molecular_representation(int imol, int rep_no);

#endif // USE_MOLECULES_TO_TRIANGLES
