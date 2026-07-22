
#include "grid-balls.hh"
#include "atom-selection-container.hh"

int main(int argc, char **argv) {

   int status = 0;
   bool use_gemmi = false;
   if (argc > 1) {
      std::string pdb_file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, use_gemmi);
      if (asc.read_success) {
         coot::protein_geometry geom;
         geom.init_standard();
         int read_number = 55;
         geom.init_refmac_mon_lib("PC1.cif", read_number);

         int imol = 0;
         coot::grid_balls_t gb(imol, asc.mol, &geom, 1.4, 4.5);
         gb.write_cavity_points("cavity-points.table");
         gb.write_subpocket_points("subpocket-points.table");
      }
   }
   return status;

}
