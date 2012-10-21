
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
