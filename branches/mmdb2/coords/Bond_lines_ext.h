
/* #include <string> */
/* #include <fstream> */
/* #include <vector> */

/* #include "Cartesian.h" */
/* #include "mmdb_manager.h" */
/* #include "mmdb-extras.h" */
/* #include "mmdb.h" */
/* #include "mmdb-crystal.h"  // should be merged with extras */


// #include "Bond_lines.h"

class Bond_lines_ext : public Bond_lines_container {

   // Add clipper facilities here with which we don't want to complicate
   // Bond_lines
 public:

   Bond_lines_ext() {}; // for use with find_molecule_middle(); 

   void find_skel_atom_bonds(atom_selection_container_t SelAtom);

   Bond_lines_ext(atom_selection_container_t SelAtom) {
      find_skel_atom_bonds(SelAtom); }

   coot::Cartesian find_molecule_middle(atom_selection_container_t SelAtom,
					float max_neighbour_dist); 

};
