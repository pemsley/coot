
// needs to have included mmdb_manager.h"
// and "mmdb-extras.h" for atom_selection_container_t
#include "Cartesian.h"


coot::Cartesian
centre_of_molecule(atom_selection_container_t SelAtom); 


atom_selection_container_t get_atom_selection(std::string t);
int fix_nucleic_acid_residue_names(atom_selection_container_t asc);
int fix_nucleic_acid_residue_name(CResidue *r); // return whether it was changed or not.

// return the number of fixed atoms
int fix_away_atoms(atom_selection_container_t asc);

// return the number of changed hydrogen names
int 
fix_hydrogen_names(atom_selection_container_t asc);

int write_atom_selection_file(atom_selection_container_t asc,
			      const std::string &filename);


// needs <iostream.h>
// 
ostream& operator<<(ostream& s, CAtom &atom);

ostream& operator<<(ostream& s, PCAtom atom); 
