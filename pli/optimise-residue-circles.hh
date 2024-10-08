
#ifndef OPTIMIZE_RESIDUE_CIRCLES_HH
#define OPTIMIZE_RESIDUE_CIRCLES_HH

#include <vector>
#include <gsl/gsl_multimin.h>

#include "residue-circle.hh"
#include "lidia-core/svg-molecule.hh"

namespace pli {

   class optimise_residue_circles {
   private:

      class angle {
      public:
	 angle(int i_atom_index_in, int ires_1_in, int ires_2_in) {
	    i_atom_index = i_atom_index_in;
	    ires_1_index = ires_1_in;
	    ires_2_index = ires_2_in;
	 }
	 int i_atom_index;
	 int ires_1_index;
	 int ires_2_index;
      };
      int status; 
      static double f(const gsl_vector *v, void *params);
      static void  df(const gsl_vector *v, void *params, gsl_vector *df); 
      static void fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
      
      std::vector<residue_circle_t> starting_circles;
      std::vector<residue_circle_t>  current_circles;

      svg_molecule_t mol;

      void numerical_gradients(gsl_vector *x, gsl_vector *df, void *params) const;
      std::vector<angle> angles;
      void setup_angles(); // uses current_circles and mol and fills angles

      // these can be unstaticed by passing an lbg_info_t to f(), df()
      // as part of params.  Hmmm... it is already.  Should be easy to
      // fix then.  But not now.
      // 
      static bool score_vs_ligand_atoms;
      static bool score_vs_ring_centres;
      static bool score_vs_other_residues;
      static bool score_vs_other_residues_for_angles;
      static bool score_vs_original_positions;
      static bool score_vs_ligand_atom_bond_length;

      double score_vs_ligand_atoms_rk;
      double score_vs_ligand_atoms_exp_scale;
      double score_vs_other_residues_kk;
      double score_vs_other_residues_exp_scale;
      double score_vs_original_positions_kk;
      double score_vs_ligand_atom_bond_length_kk;

   public:
      // we pass two vectors here because (for trajectory-view) we
      // don't want to restart the minimisation with the current
      // positions - we want to minimise against the (constant)
      // original positions.
      optimise_residue_circles(const std::vector<residue_circle_t> &r, // starting points
			       const std::vector<residue_circle_t> &c, // current points
			       const svg_molecule_t &mol,
			       const std::vector<int> &primary_indices);
      std::pair<int, std::vector<residue_circle_t> > solution() const;
      // return GSL minimisation status;
      int get_gsl_min_status() const { return status; }
      std::vector<int> primary_indices;  // indices in residues_circles
                                         // of the primary residues.
   };

}


#endif // OPTIMIZE_RESIDUE_CIRCLES_HH

