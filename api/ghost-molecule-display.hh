#ifndef GHOST_MOLECULE_DISPLAY_HH
#define GHOST_MOLECULE_DISPLAY_HH

#include <clipper/core/coords.h>
#include "coords/graphical-bonds-container.hh"
#include "ncs.hh"

namespace coot {
   class ghost_molecule_display_t {
   public:
      clipper::RTop_orth rtop;
      int SelectionHandle;
      graphical_bonds_container bonds_box;
      // Mesh mesh; 20221025-PE not in this directory

      std::string name;
      std::string chain_id;
      std::string target_chain_id;  // this operator matches to this chain.
      bool display_it_flag;
      std::vector<int> residue_matches;
      ghost_molecule_display_t() {
         SelectionHandle = -1;
         display_it_flag = false; }
      ghost_molecule_display_t(const clipper::RTop_orth &rtop_in,
                               int SelHnd_in,
                               const std::string &name_in) :
         rtop(rtop_in), SelectionHandle(SelHnd_in), name(name_in) {
         display_it_flag = 1;
      }
      virtual ~ghost_molecule_display_t() = default;
      virtual void update_bonds(mmdb::Manager *mol); // the parent's mol

      bool is_empty() const { return (SelectionHandle == -1); }
      ncs_residue_info_t get_differences(mmdb::Residue *this_residue_p,
                                         mmdb::Residue *master_residue_p,
                                         float main_chain_weight) const;
      friend std::ostream& operator<<(std::ostream &s, const ghost_molecule_display_t &ghost);
   };

}

#endif // GHOST_MOLECULE_DISPLAY_HH
