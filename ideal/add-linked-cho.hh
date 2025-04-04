
#include <clipper/core/xmap.h>
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   namespace cho {

      // (define *add-linked-residue-tree-correlation-cut-off* 0.50)

      //! checks to see if residue_pp is well fitting and doesn't clash or symmetry clash
      bool is_well_fitting(mmdb::Residue *residue_p,
                           mmdb::Manager *mol,
                           clipper::Xmap<float> &xmap,
                           const protein_geometry &geom);

      //! \brief return 1 if this residue clashes with the symmetry-related
      //!  atoms of the same molecule.
      //!
      //! 0 means that it did not clash,
      //! -1 means that the residue or molecule could not be found or that there
      //!    was no cell and symmetry.
      int clashes_with_symmetry(mmdb::Manager *mol, const residue_spec_t  &res_spec, float clash_dist,
                                const protein_geometry &geom);


      //! do the thing. This function calls the others.
      //! res-pair is new-link-type and new-res-type
      residue_spec_t
      add_linked_residue_add_cho_function(mmdb::Manager *mol,
                                          int imol, // need for dictionaries
                                          const residue_spec_t &parent,
                                          const std::pair<std::string, std::string> &res_pair,
                                          float new_atoms_b_factor,
                                          protein_geometry &geom,
                                          const clipper::Xmap<float> *xmap); // ca be null

      residue_spec_t
      add_linked_residue(mmdb::Manager *mol,
                         int imol, // because dictionaries
                         const residue_spec_t &parent,
                         const std::pair<std::string, std::string> &new_link_types,
                         float new_atoms_b_factor,
                         int mode,
                         protein_geometry &geom,
                         const clipper::Xmap<float> *xmap);

      // make this a coord utils function? (And others here, perhaps).
      // This adds a header - that is all..
      void make_link(mmdb::Manager *mol, const coot::atom_spec_t &spec_1,
                     const coot::atom_spec_t &spec_2,
                     const std::string &link_name, float length,
                     const coot::protein_geometry &geom);

      void asn_hydrogen_position_swap(std::vector<std::pair<bool, mmdb::Residue *> > residues);

      residue_spec_t add_linked_residue_by_atom_torsions(mmdb::Manager *mol,
                                                         const residue_spec_t &parent,
                                                         const std::pair<std::string, std::string> &new_link_types,
                                                         protein_geometry &geom,
                                                         float new_atoms_b_factor);

      std::pair<bool, mmdb::Residue *> add_residue(mmdb::Manager *mol,
                                                   mmdb::Residue *new_res,
                                                   const std::string &chain_id_in);

      mmdb::Residue * copy_and_add_residue_to_chain(mmdb::Manager *mol,
                                                    mmdb::Chain *this_model_chain,
                                                    mmdb::Residue *add_model_residue,
                                                    bool new_resno_by_hundreds_flag);
      // return state, max_resno + 1, or 0, 1 of no residues in chain.
      //
      // new_res_no_by_hundreds is default false
      std::pair<short int, int> next_residue_number_in_chain(mmdb::Chain *w,
                                                             bool new_res_no_by_hundreds);

      bool is_het_residue(mmdb::Residue *residue_p);

      void replace_coords(mmdb::Manager *fragment_mol, mmdb::Manager *mol);

   }
}
