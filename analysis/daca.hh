
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
         float box_width;
         box_index_t(int ii, int jj, int kk) : idx_x(ii), idx_y(jj), idx_z(kk) {
            box_width = 1.0;
         }
         box_index_t(const clipper::Coord_orth &pos);
         bool operator<(const box_index_t &other) const;
         clipper::Coord_orth coord_orth() const;
         // friend std::ostream &operator<<(std::ostream &s, const box_index_t &bi);
         float d() const;
         float d_squared() const;
      };


      double get_radius(const std::string &ele) const;
      // std::ostream &operator<<(std::ostream &s, const box_index_t &bi);
      std::vector<std::pair<mmdb::Atom *, float> >
      solvent_exposure_old_version(int selhnd, mmdb::Manager *mol) const;
      std::vector<std::pair<mmdb::Residue *, float> >
      solvent_exposure_old_version_v2(mmdb::Manager *mol, bool side_chain_only=true) const;
      std::vector<std::pair<mmdb::Residue *, float> >
      solvent_exposure(mmdb::Manager *mol, bool side_chain_only=true) const;

   private:
      enum mode_t {REFERENCE, ANALYSIS};
      // the atoms in the reference fragments have a particular order and (for each fragmemt)
      // are centred on the origin.
      std::map<std::string, std::vector<std::vector<clipper::Coord_orth> > > reference_fragments;
      void fill_reference_fragments();
      bool boxes_have_been_resized = false;
      std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > > boxes; // reference
      std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > > boxes_for_testing; // this pdb
      std::vector<mmdb::Residue *> helical_residues; // fill this before calling calculate_daca()
      std::map<float, float> envelope_distance_map; // the key is squared values
      void fill_helix_flags(mmdb::Model *model_p, mmdb::Manager *mol);
      std::pair<bool, clipper::RTop_orth>  // return status also, as this can fail.
      get_frag_to_reference_rtop(const std::string &res_name,
                                 const unsigned int &frag_idx,
                                 const std::vector<mmdb::Atom *> &fragment_atoms) const;
      void add_to_box(mode_t mode,
                      const std::string &res_name,
                      bool residue_is_helical_flag,
                      unsigned int frag_index,
                      const box_index_t &box_index,
                      const std::string &atom_type,
                      unsigned int counts=1);
      int get_reference_counts(const std::string &res_name,
                               bool residue_is_helical_flag,
                               unsigned int frag_index,
                               const box_index_t &box_index,
                               const std::string &atom_type) const;
      std::vector<std::vector<std::string> > atom_names_for_fragments(const std::string &res_name) const;
      std::vector<std::vector<mmdb::Atom *> > get_daca_fragments(mmdb::Residue *reference_residue_p) const;
      std::vector<std::pair<mmdb::Atom *, std::string> >
         make_typed_atoms(mmdb::Model *model_p, const protein_geometry &geom) const;

      // maybe this should be a float
      int calculate_daca(mmdb::Residue *reference_residue_p,
                         const std::vector<std::pair<mmdb::Atom *, std::string> > &typed_atoms,
                         mode_t mode);
      bool atom_is_close_to_a_residue_atom(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const;
      bool atom_is_neighbour_mainchain(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const;
      void debug_boxes(const std::string &debug_prefix="") const;
      float gompertz_scale(const float &dist);
      void compare_boxes() const;
      void presize_boxes(mode_t mode=REFERENCE);
      void normalize();
      void normalize_v2();
      void envelope();
      void smooth();

   public:
      daca() { fill_reference_fragments(); boxes_have_been_resized = false; }
      void write_tables_using_reference_structures_from_dir(const std::string &input_pdb_files_dir_name,
                                                            const std::string &output_tables_dir);
      void read_tables(const std::string &dir);
      void read_many_tables(const std::vector<std::string> &dirs);
      void write_tables(const std::string &dir_name) const;
      void score_molecule(const std::string &pdb_file_name);
      void cook();

   };


}
