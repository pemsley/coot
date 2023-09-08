
#include "rna-backbone.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/atom-selection-container.hh"

// speed notes 22710 lines (5 lines per log are not base stats) 690 logs
// in 126 seconds. 19260 baes in 125 seconds -> 152 bases/second
// 5.5 structures/second

int main(int argc, char **argv) {

   int status = 0;

   if (argc > 2) {
      std::string pdb_file_name(argv[1]);
      std::string mtz_file_name(argv[2]);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true);
      if (asc.mol) {
         std::string f_col = "DELFWT";
         std::string phi_col = "PHDELWT";
         std::string w_col = "";
         clipper::Xmap<float> xmap;
         coot::util::map_fill_from_mtz(&xmap, mtz_file_name, f_col, phi_col, w_col, 0, 1.8);

         coot::rna_backbone_t rb(asc.mol, xmap);
      }
   }

   return status;
}
