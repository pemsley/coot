#ifndef PLI_DOTS_REPRESENTATION_INFO_HH
#define PLI_DOTS_REPRESENTATION_INFO_HH

#include <vector>
#include <clipper/core/coords.h>
#include "geometry/residue-and-atom-specs.hh"
#include "solvent-exposure-difference-helper.hh"
#include "utils/colour-holder.hh"
#include "api/instancing.hh"
#include "api/coot-colour.hh"

namespace pli {

   // Is this for the 3D representation of the ligand? I think so.

   class dots_representation_info_t {

      bool is_closed;
      void pure_points(mmdb::Manager *mol); // don't surface mol, the surface points *are* the
                                            // (synthetic) atoms in mol.

      std::vector<std::pair<coot::colour_t, std::vector<clipper::Coord_orth> > > points;
      double get_radius(const std::string &ele) const;
      coot::colour_t get_colour(const std::string &ele) const;

   public:
      dots_representation_info_t() { is_closed = false; }
      explicit dots_representation_info_t(const std::string &n) : name(n) { is_closed = false; }
      explicit dots_representation_info_t(mmdb::Manager *mol);
      dots_representation_info_t(mmdb::Manager *mol, mmdb::Manager *mol_exclude);

      std::string name;
      std::string get_name() const { return name; }
      coot::instanced_mesh_t imm; // was Instanced_Markup_Mesh when it was in src

      std::vector<std::pair<coot::atom_spec_t, float> >
      solvent_accessibilities(mmdb::Residue *res_ref,
                              const std::vector<mmdb::Residue *> &filtered_residues) const;

      std::vector<std::pair<mmdb::Atom *, float> > solvent_exposure(int SelHnd_in, mmdb::Manager *mol) const;

      std::vector<solvent_exposure_difference_helper_t>
      solvent_exposure_differences(mmdb::Residue *res_ref, const std::vector<mmdb::Residue *> &residues) const;

      bool is_open_p() const {
         int r = 1 - is_closed;
         return r;
      }
      void close_yourself() {
         imm.clear();
         is_closed = true;
      }

      // mol_exclude can be NULL.
      void add_dots(int SelHnd_in, mmdb::Manager *mol, mmdb::Manager *mol_exclude,
                    double dots_density, const coot::colour_t &single_col, bool use_single_col);

   };

}


#endif // DOTS_REPRESENTATION_INFO_HH
