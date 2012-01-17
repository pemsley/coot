
#ifndef HAVE_VECTOR
#include <vector>
#define HAVE_VECTOR
#endif

#include "clipper/core/coords.h"

#include <mmdb/mmdb_manager.h>
#include "atom-quads.hh"

namespace coot {


   class alt_confed_chi_angles {
   public:
      alt_confed_chi_angles(const std::string &alt_conf_in, const std::vector<std::pair<int, float> > &chi_angles_in) {
         alt_conf = alt_conf_in;
         chi_angles = chi_angles_in;
      }
      std::string alt_conf;
      std::vector<std::pair<int, float> > chi_angles;
      void PROs_to_one_chi();
   };

   class alt_confed_atom_index_quad {
   public:
      alt_confed_atom_index_quad(const std::string &alt_conf_in, std::vector<atom_index_quad> &quad_in) {
         quad = quad_in;
         alt_conf = alt_conf_in;
      }
      std::string alt_conf;
      std::vector<atom_index_quad> quad;
   };

   class residue_named_chi_angle_atom_name_quad_set_t {
   public:
      std::string residue_type;
      std::vector<atom_name_quad> name_quad;
      residue_named_chi_angle_atom_name_quad_set_t() {};
      residue_named_chi_angle_atom_name_quad_set_t(const std::string &residue_type_in,
                                                   const std::vector<atom_name_quad> &name_quad_set_in) {
         residue_type = residue_type_in;
         name_quad = name_quad_set_in;
      }
      residue_named_chi_angle_atom_name_quad_set_t(const std::string &residue_type_in,
                                                   const atom_name_quad &name_quad_in) {
         residue_type = residue_type_in;
         name_quad.push_back(name_quad_in);
      }
      void add_torsion_bond_by_name(const std::string &atom_name_1,
                                    const std::string &atom_name_2,
                                    const std::string &atom_name_3,
                                    const std::string &atom_name_4);
   };

   class primitive_chi_angles {
      std::string residue_name;
      CResidue *residue;
      std::vector<atom_name_quad> get_atom_name_quads() const;
      std::vector<coot::atom_index_quad>
      get_quads(const std::vector<coot::atom_name_quad> &atom_name_quads,
                CResidue *residue) const;
      std::vector<coot::alt_confed_atom_index_quad>
      get_quads_using_altconfs(const std::vector<coot::atom_name_quad> &atom_name_quads,
			       CResidue *residue) const;
      clipper::Coord_orth atom_to_co(CAtom *at) const;
      std::vector<residue_named_chi_angle_atom_name_quad_set_t> chi_angle_atoms_for_residue_type;
      std::vector<coot::atom_index_quad>
      get_atom_index_quads(const std::vector<coot::atom_name_quad> &atom_name_pairs,
                           const PPCAtom atoms, int nresatoms) const;
      void setup_chi_atom_quads();
      void add_chi_quad(const std::string &residue_type,
                        const std::string &atom_name_1,
                        const std::string &atom_name_2,
                        const std::string &atom_name_3,
                        const std::string &atom_name_4);

   public:
      primitive_chi_angles(CResidue *res_in) {
	 residue = res_in;
	 if (residue)
	    residue_name = residue->GetResName();
	 setup_chi_atom_quads();
      }
      // This throws a runtime_error exception on ALA gly or other
      // unknown residue type.
      std::vector<alt_confed_chi_angles> get_chi_angles();
   };

}
