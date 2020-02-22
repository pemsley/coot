
#include <vector>
#include <utility>
#include <string>
#include <unordered_set>

#include <clipper/core/coords.h>
#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"


namespace coot {

   class daca {

   public:

      class box_index_t {
      public:
            int idx_x;
            int idx_y;
            int idx_z;
            box_index_t(const clipper::Coord_orth &pos);
      bool operator<(const box_index_t &other) const;
      clipper::Coord_orth coord_orth() const;
      };

   private:
      // the atoms in the reference fragments have a particular order and (for each fragmemt)
      // are centred on the origin.
      std::map<std::string, std::vector<std::vector<clipper::Coord_orth> > > reference_fragments;
      void fill_reference_fragments();
      std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > > boxes;
      std::vector<mmdb::Residue *> helical_residues; // fill this before calling calculate_daca()
      void fill_helix_flags(mmdb::Model *model_p, mmdb::Manager *mol);
      clipper::RTop_orth
      get_frag_to_reference_rtop(const std::string &res_name,
                                 const unsigned int &frag_idx,
                                 const std::vector<mmdb::Atom *> &fragment_atoms) const;
      void add_to_box(const std::string &res_name,
                      bool residue_is_helical_flag,
                      unsigned int frag_index,
                      const box_index_t &box_index,
                      const std::string &atom_type);
      std::vector<std::vector<std::string> > atom_names_for_fragments(const std::string &res_name) const;
      std::vector<std::vector<mmdb::Atom *> > get_daca_fragments(mmdb::Residue *reference_residue_p) const;
      std::vector<std::pair<mmdb::Atom *, std::string> >
         make_typed_atoms(mmdb::Model *model_p, const protein_geometry &geom) const;
      void calculate_daca(mmdb::Residue *reference_residue_p,
                          const std::vector<std::pair<mmdb::Atom *, std::string> > &typed_atoms);
      bool atom_is_close_to_a_residue_atom(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const;
      bool atom_is_neighbour_mainchain(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const;
      void debug_boxes() const;
      void write_tables() const;

   public:
      daca() { fill_reference_fragments(); }
      void write_tables_using_reference_structures_from_dir(const std::string &dir_name);

   };


}
