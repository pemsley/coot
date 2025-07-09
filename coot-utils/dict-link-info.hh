
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   // trivial helper class to get specs and distance for atoms when a
   // link is made.
   //
   class dict_link_info_t {
      bool check_for_order_switch(mmdb::Residue *residue_ref,
                                  mmdb::Residue *residue_new,
                                  const std::string &link_type,
                                  const protein_geometry &geom) const;
   public:
      // this can throw a std::runtime_error
      dict_link_info_t (mmdb::Residue *residue_ref, mmdb::Residue *residue_new,
                        const std::string &link_type, const protein_geometry &geom);
      atom_spec_t spec_ref;
      atom_spec_t spec_new;
      double dist;
   };
}
