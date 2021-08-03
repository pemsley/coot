/* lidia-core/chemical-feature-clusters.hh
 * 
 * Copyright 2012 by the University of Oxford
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef CHEMICAL_FEATURE_CLUSTERS_HH
#define CHEMICAL_FEATURE_CLUSTERS_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <vector>
#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"
// #include "use-rdkit.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>



namespace chemical_features {
   RDKit::MolChemicalFeatureFactory *get_feature_factory();
}


namespace coot {

   class simple_chemical_feature_attributes {
   public:
      simple_chemical_feature_attributes() {}
      simple_chemical_feature_attributes(const std::string &type_in,
					 const clipper::Coord_orth &pos_in,
					 int imol_in,
					 const residue_spec_t &rs) {
	 type = type_in;
	 pos = pos_in;
	 imol = imol_in;
	 residue_spec = rs;
      }
      std::string type;
      clipper::Coord_orth pos;
      int imol;
      residue_spec_t residue_spec;
   };

   class chem_feat_solvated_ligand_spec {

   public:
      residue_spec_t ligand_spec;
      std::vector<residue_spec_t> waters;
      mmdb::Manager *mol;
      int imol; // the molecule number that this molecule came from
                // (needed so that we can navigate in the button-press callback)
      chem_feat_solvated_ligand_spec(const residue_spec_t &ligand_spec_in,
				     const std::vector<residue_spec_t> &waters_in,
				     mmdb::Manager *mol_in,
				     int imol_in) {
	 ligand_spec = ligand_spec_in;
	 waters = waters_in;
	 mol = mol_in;
	 imol = imol_in;
      }
   };

   class chem_feat_solvated_ligand : public chem_feat_solvated_ligand_spec {

      void init_residue();

   public:
      chem_feat_solvated_ligand(const residue_spec_t &ligand_spec_in,
				const std::vector<residue_spec_t> &waters_in,
				mmdb::Manager *mol_in,
				int imol_in) : chem_feat_solvated_ligand_spec(ligand_spec_in, waters_in, mol_in, imol_in) { init_residue(); }
      explicit chem_feat_solvated_ligand(const chem_feat_solvated_ligand_spec &ligand_spec_in) :
         chem_feat_solvated_ligand_spec(ligand_spec_in) { init_residue(); }

      mmdb::Residue *residue;
   };

   class chem_feat_clust {

   public:
      class water_attribs {
      public:
	 water_attribs(int imol_in,
                       unsigned int ligand_idx_in,
		       unsigned int water_spec_idx_in,
		       mmdb::Atom *water_atom_in,
		       const clipper::Coord_orth &pt) : imol(imol_in), pos(pt) {
	    ligand_idx = ligand_idx_in;
	    water_spec_idx = water_spec_idx_in;
	    atom_p = water_atom_in;
	 }
         int imol;
	 unsigned int ligand_idx;
	 unsigned int water_spec_idx;
	 mmdb::Atom *atom_p;
	 clipper::Coord_orth pos;
	 residue_spec_t residue_spec() const { return residue_spec_t(atom_p->residue) ; }
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
      
      std::vector<simple_chemical_feature_attributes>
      get_chemical_features(int imol, residue_spec_t lig_spec, mmdb::Manager *mol);

      // check that the ligand specs point to real residues
      bool fill_ligands(const std::vector<chem_feat_solvated_ligand_spec> &ligands_in);

      bool check_dictionaries(); // check that we have the dictiionaries for the given ligands

      void fill_waters();

      void align();

      double water_dist_cutoff;
      std::vector<clipper::Coord_orth> get_ligands_coords() const;
      bool is_close_to_a_ligand_atom(const clipper::Coord_orth &pt_test,
				     const std::vector<clipper::Coord_orth> &ligand_atom_positions) const;
      
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
      // waters in water_positions must be with dist_cutoff_in of any ligand atom
      // (not just its own) - say 5A.
      // 
      chem_feat_clust(const std::vector<residue_spec_t> &protein_residues,
		      const std::vector<chem_feat_solvated_ligand_spec> &ligands,
		      double water_dist_cutoff_in,
		      const protein_geometry *geometry, // maybe should be non-const pointer
		      bool do_alignment = false);

      std::vector<water_attribs> get_water_positions() const { return water_positions; }

      // not const because protein_geometry might change
      std::vector<simple_chemical_feature_attributes> get_chemical_features();

   };

}

#endif // MAKE_ENHANCED_LIGAND_TOOLS

#endif // CHEMICAL_FEATURE_CLUSTERS_HH

