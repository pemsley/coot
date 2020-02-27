

namespace coot {

      /* 
Model composition
   Non-hydrogen atoms
   Protein residues
   Zinc ions (or other ligands) 

RMS deviations
   Bonds
   Angles

Validation
   Number of C-beta deviation outliers
   Atom overlap score [Clash score]
   Good rotamers (%)
   Ramachandran plot
       Favoured
       Outliers
   */

   class model_composition_stats_t {
   public:
      int n_non_hydrogen_atoms;
      int n_atoms;
      int n_protein_residues;
      int n_rna_residues;
      int n_dna_residues;
      int n_hetgroups;
      int n_metal_atoms;
      float rmsd_bonds;
      float rmsd_angles;
      float overall_atom_overlap_volume;
      float fraction_good_rotamers;
      float fraction_ramachandran_favoured;
      float fraction_ramachandran_outliers;
      model_composition_stats_t() {
	 n_non_hydrogen_atoms = -1;
	 n_atoms = -1;
	 n_rna_residues = -1;
	 n_protein_residues = -1;
	 n_dna_residues = -1;
	 n_hetgroups = -1;
	 n_metal_atoms = -1;
	 rmsd_bonds = -1;
	 rmsd_angles = -1;
	 overall_atom_overlap_volume = -1;
	 fraction_good_rotamers = -1;
	 fraction_ramachandran_favoured = -1;
	 fraction_ramachandran_outliers = -1;
      }
   };

}
