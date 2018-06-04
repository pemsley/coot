

#ifdef USE_PYTHON

// Martin's Triangles

// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation(int imol, PyObject *atom_selection_py, PyObject *ColorScheme_py, PyObject *style_py);
void remove_molecular_represenation(int imol, int rep_no);
#endif // USE_PYTHON
