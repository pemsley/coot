#include "coot-map-utils.hh"
#include "coot-coord-utils.hh"
#include "atom-overlaps.hh"

namespace coot {
   namespace cho {
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
   }
}

int
coot::cho::clashes_with_symmetry(mmdb::Manager *mol, const coot::residue_spec_t &res_spec, float clash_dist,
                                 const coot::protein_geometry &geom) {

   int r = -1;
   mmdb::Residue *residue_p = util::get_residue(res_spec, mol);
   if (mol) {
      if (residue_p) {
         std::vector<mmdb::Residue *> dummy; // neighbours
         atom_overlaps_container_t ao(residue_p, dummy, mol, &geom);
         std::vector<coot::atom_overlap_t> v = ao.symmetry_contacts(clash_dist);
         if (v.empty())
            r = 0;
         else
            r = 1;
      }
   }
   return r;
}

bool
coot::cho::is_well_fitting(mmdb::Residue *residue_p,
                           mmdb::Manager *mol,
                           clipper::Xmap<float> &xmap,
                           const coot::protein_geometry &geom) {

   float add_linked_residue_tree_correlation_cut_off = 0.50;
   float clash_dist = 2.0;
   float atom_radius = 1.6;

   bool status = false;
   float radius = 4.0;
   residue_spec_t res_spec(residue_p);
   std::vector<mmdb::Residue *> neighbours = residues_near_residue(residue_p, mol, radius);
   std::vector<residue_spec_t> residues_for_masking;
   for(mmdb::Residue *r : neighbours)
      residues_for_masking.push_back(residue_spec_t(r));
   std::vector<residue_spec_t> residues_for_cc = { res_spec };
   unsigned short int atom_mask_mode = 0; // all atom

   float c = util::map_to_model_correlation(mol, residues_for_cc, residues_for_masking, atom_mask_mode, atom_radius, xmap);
   if (c > add_linked_residue_tree_correlation_cut_off) {
      int symm_clash = clashes_with_symmetry(mol, res_spec, clash_dist, geom);
      if (symm_clash == 0) {
         status = true;
      }
   }
   return status;
}
