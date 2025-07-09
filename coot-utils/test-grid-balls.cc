
#include "grid-balls.hh"
#include "atom-selection-container.hh"

int main(int argc, char **argv) {

   int status = 0;
   bool use_gemmi = false;
   if (argc > 1) {
      std::string pdb_file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, use_gemmi);
      if (asc.read_success) {
         coot::grid_balls_t gb(asc.mol, 1.4, 4.5);
      }
   }
   return status;

}
