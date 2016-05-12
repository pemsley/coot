
#ifndef CHEMICAL_FEATURE_CLUSTERS_HH
#define CHEMICAL_FEATURE_CLUSTERS_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <vector>
#include "coot-utils/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"
// #include "use-rdkit.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

namespace coot {

   class chem_feat_solvated_ligand_spec {

   public:
      residue_spec_t ligand_spec;
      std::vector<residue_spec_t> waters;
      mmdb::Manager *mol;
      chem_feat_solvated_ligand_spec(const residue_spec_t &ligand_spec_in,
				     const std::vector<residue_spec_t> &waters_in,
				     mmdb::Manager *mol_in) {
	 ligand_spec = ligand_spec_in;
	 waters = waters_in;
	 mol = mol_in;
      }
   };

   class chem_feat_solvated_ligand : public chem_feat_solvated_ligand_spec {

      void init_residue();

   public:
      chem_feat_solvated_ligand(const residue_spec_t &ligand_spec_in,
				const std::vector<residue_spec_t> &waters_in,
				mmdb::Manager *mol_in) : chem_feat_solvated_ligand_spec(ligand_spec_in, waters_in, mol_in) { init_residue(); }
      chem_feat_solvated_ligand(const chem_feat_solvated_ligand_spec &ligand_spec_in) :
				chem_feat_solvated_ligand_spec(ligand_spec_in) { init_residue(); }
      
      mmdb::Residue *residue;
   };

   class chem_feat_clust {

   public:
      class water_attribs {
      public:
	 water_attribs(unsigned int ligand_idx_in,
		       unsigned int water_spec_idx_in,
		       mmdb::Atom *water_atom_in,
		       const clipper::Coord_orth &pt) {
	    ligand_idx = ligand_idx_in;
	    water_spec_idx = water_spec_idx_in;
	    atom_p = water_atom_in;
	    pos = pt;
	 }
	 unsigned int ligand_idx;
	 unsigned int water_spec_idx;
	 mmdb::Atom *atom_p;
	 clipper::Coord_orth pos;
      };

   private:

      bool setup_success;
      std::vector<chem_feat_solvated_ligand> ligands;
      const protein_geometry *geometry_p;

      // for the water positions, what do we want to store along with the coord orth?
      // 
      // the idx in the ligands array, so see trivial class above
      // 
      std::vector<water_attribs> water_positions;
      
      bool get_chemical_features(unsigned int idx, residue_spec_t lig_spec,
				 mmdb::Manager *mol);

      bool cluster_waters(const std::vector<chem_feat_solvated_ligand_spec> &ligands);

      // check that the ligand specs point to real residues
      bool fill_ligands(const std::vector<chem_feat_solvated_ligand_spec> &ligands_in);

      bool check_dictionaries(); // check that we have the dictiionaries for the given ligands

      void fill_waters();

      void align();

   public:

      // protein_residues includes waters
      // 
      // align other molecules onto the residues of the first molecule
      // if needed, moving the atoms of the given Manager as needed
      // (therefore, if do_alignment is true, the molecule in Manager
      // will need to be rebonded and redrawn after this constructor
      // is used)
      //
      // Waters (HOHs) are not used for alignment.
      // 
      chem_feat_clust(const std::vector<residue_spec_t> &protein_residues,
		      const std::vector<chem_feat_solvated_ligand_spec> &ligands,
		      const protein_geometry *geometry, // maybe should be non-const pointer
		      bool do_alignment = false);

      void cluster_waters();

      void cluster_ligand_chemical_features();

      std::vector<water_attribs> get_water_positions() const { return water_positions; }
      
   };

}


namespace chemical_features {

   RDKit::MolChemicalFeatureFactory *get_feature_factory();

}

#endif // MAKE_ENHANCED_LIGAND_TOOLS

#endif // CHEMICAL_FEATURE_CLUSTERS_HH

