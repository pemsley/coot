#ifndef COOT_MOLECULE_HH
#define COOT_MOLECULE_HH

#include <utility>

#include <clipper/core/xmap.h>
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/coot-rama.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

namespace coot {

   enum { UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2,
	  COLOUR_BY_CHAIN_BONDS=3,
	  CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
	  BONDS_NO_HYDROGENS=15,
	  CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
	  CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
	  CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS=17,
	  COLOUR_BY_MOLECULE_BONDS=8,
	  COLOUR_BY_RAINBOW_BONDS=9,
	  COLOUR_BY_B_FACTOR_BONDS=10,
	  COLOUR_BY_OCCUPANCY_BONDS=11,
	  COLOUR_BY_USER_DEFINED_COLOURS____BONDS=12,
	  COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS=13 };

   class molecule_t {

      class molecule_save_info_t {
      public:
         std::pair<time_t, unsigned int> last_saved;
         unsigned int modification_index;
         molecule_save_info_t() : last_saved(std::make_pair(0,0)), modification_index(0) {}
         void new_modification() {
            modification_index++;
         }
         void made_a_save() {
            // this is called when the server says that it has saved the file
            last_saved.first = time(nullptr);
            last_saved.second = modification_index;
         }
         bool have_unsaved_changes() const {
            return modification_index > last_saved.second;
         }
      };

      molecule_save_info_t save_info;

      void makebonds(coot::protein_geometry *geom, std::set<int> &no_bonds_to_these_atoms);

#if defined __has_builtin
#if __has_builtin (__builtin_FUNCTION)
      void make_bonds_type_checked(coot::protein_geometry *geom, const char *s = __builtin_FUNCTION());
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, const char *s = __builtin_FUNCTION());
#else
      void make_bonds_type_checked(coot::protein_geometry *geom, const char *s = 0);
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, const char *s =0);
#endif
#else // repeat above
      void make_bonds_type_checked(coot::protein_geometry *geom, const char *s = 0);
      void make_bonds_type_checked(coot::protein_geometry *geom, const std::set<int> &no_bonds_to_these_atom_indices, const char *s =0);
#endif

      int bonds_box_type; // public accessable via get_bonds_box_type(); // wass Bonds_box_type()
      graphical_bonds_container bonds_box;
      int get_bonds_box_type() const { return bonds_box_type; }


   public:
      atom_selection_container_t atom_sel;
      molecule_t() {}
      explicit molecule_t(atom_selection_container_t asc) : atom_sel(asc) { bonds_box_type = UNSET_TYPE; }
      clipper::Xmap<float> xmap; // public because the filling function needs access

      // utils

      bool is_valid_model_molecule() const;
      bool is_valid_map_molecule() const;
      std::pair<bool, coot::residue_spec_t> cid_to_residue_spec(const std::string &cid);

      // model utils

      // returns either the specified atom or null if not found
      mmdb::Atom *get_atom(const coot::atom_spec_t &atom_spec) const;
      // returns either the specified residue or null if not found
      mmdb::Residue *get_residue(const coot::residue_spec_t &residue_spec) const;

      bool have_unsaved_changes() const { return save_info.have_unsaved_changes(); }

      // model analysis functions

      std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > ramachandran_validation() const;

      // model-changing functions

      int flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf);

      // map functions

      int writeMap(const std::string &file_name) const;

   };
}


#endif // COOT_MOLECULE_HH
