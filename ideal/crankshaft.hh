
#ifndef CRANKSHAFT_HH
#define CRANKSHAFT_HH

#include <vector>
#include "gsl/gsl_multimin.h"

#include <clipper/core/coords.h>

#include <geometry/residue-and-atom-specs.hh>
#include "phi-psi.hh"
#include "zo-rama.hh"

namespace coot {

   class crankshaft_set {
      // do it with atom indices of atoms in an atom selection
      // i.e.
      // atom-idx 0: res-0-C # non-moving
      // atom-idx 1: res-1-N # non-moving
      // atom-idx 2: res-1-C
      // atom-idx 3: res-1-O
      // atom-idx 4: res-2-N
      // atom-idx 5: res-2-H # optional, otherwise null
      // atom-idx 6: res-C-C # non-moving
      // atom-idx 7: res-3-N # non-moving

      // both phi,psi pairs - for the 2 CAs involved
      //
      std::pair<phi_psi_t, phi_psi_t> get_phi_phis(const clipper::Coord_orth &C_pos,
						   const clipper::Coord_orth &N_pos) const;

      phi_psi_t phi_psi(const clipper::Coord_orth &C_pos,
			const clipper::Coord_orth &N_pos) const;

   public:
      crankshaft_set() {
	 ca_1 = 0;
	 ca_2 = 0;
      }
      // cranshaft res_1 and res_2. res_0 and res_3 are for determining the 
      // phi and psi angles
      //
      // throw runtime_error when no all atoms are present
      //
      crankshaft_set(mmdb::Residue *res_0,
		     mmdb::Residue *res_1,
		     mmdb::Residue *res_2,
		     mmdb::Residue *res_3);
      std::vector<mmdb::Atom *> v;
      mmdb::Atom *ca_1;
      mmdb::Atom *ca_2;
      std::pair<phi_psi_t, phi_psi_t> phi_psis(float ang) const;
      phi_psi_t phi_psi(float ang) const;
      void move_the_atoms(float ang);
   };

   // perhaps this can be extended to an arbitrary number of residues in due course
   //
   class triple_crankshaft_set {

      crankshaft_set cs[3];
      std::vector<std::string> residue_types;

   public:
      triple_crankshaft_set(mmdb::Residue *res_0,
			    mmdb::Residue *res_1,
			    mmdb::Residue *res_2,
			    mmdb::Residue *res_3,
			    mmdb::Residue *res_4,
			    mmdb::Residue *res_5,
			    const std::vector<std::string> &rts);
      // maybe we can save zorts in this class with a const ref?
      // above constructor will need changing however.
      triple_crankshaft_set(const residue_spec_t &spec_first_residue,
			    const zo::rama_table_set &zorts,
			    mmdb::Manager *mol);
      const crankshaft_set &operator[](unsigned int i) const { return cs[i]; }
      phi_psi_t phi_psi(unsigned int peptide_idx, float angle) const;
      std::pair<phi_psi_t, phi_psi_t> phi_psis_last(float ang_third) const;
      void move_the_atoms(float angles[]);
      const std::string &residue_type(unsigned int idx) const { return residue_types[idx]; }
      float log_prob(const phi_psi_t &pp, unsigned int peptide_index, const zo::rama_table_set &zorts) const {
	 return zorts.value(pp, residue_types[peptide_index]);
      }
   };

   // this class does not do the right thing when used for residues with alt confs
   //
   class crankshaft {

      class optimize_a_triple {
      public:
	 static double f(const gsl_vector *v, void *params);
	 static void  df(const gsl_vector *v, void *params, gsl_vector *df);
	 static void fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
      };

      class param_holder_t {
      public:
	 param_holder_t(const zo::rama_table_set &rts_in,
			const triple_crankshaft_set &tcs_in)
	    : zorts(rts_in), tcs(tcs_in) { }
	 const zo::rama_table_set &zorts;
	 const triple_crankshaft_set &tcs; // aha! const ref in constructed object
      };

      mmdb::Manager *mol;
      // caller ensures that res_1 and res_2 are valid
      std::vector<mmdb::Atom *> get_mainchain_atoms(mmdb::Residue *res_1, mmdb::Residue *res_2) const;
      mmdb::Atom *get_atom(mmdb::Residue *res_1, const std::string &atom_name) const;
   public:
      crankshaft(mmdb::Manager *mol_in) { mol = mol_in; }

      std::pair<float, float> probability_of_spin_orientation(const std::pair<phi_psi_t, phi_psi_t> &ppp,
							      const std::string &residue_table_type_1,
							      const std::string &residue_table_type_2,
							      const zo::rama_table_set &zorts) const;
      // for just the first (C-N-CA-C-N) phi,psi we don't need the second residue type
      float probability_of_spin_orientation(const phi_psi_t &pp,
					    const std::string &residue_table_type_1,
					    const zo::rama_table_set &zorts) const;

      // pass the spec of the lower residue number (of a pair of residues in
      // the peptide)
      //
      std::vector<std::pair<float, float> > spin_search(const residue_spec_t &spec,
							const zo::rama_table_set &zorts,
							int n_samples=60) const;

      // change the return type at some stage
      // not const because we can change the atoms of mol if apply_best_angles_flag is set
      //
      // Note that a 3D grid-search is not the best way to find the minimum. Also, it finds
      // the grid point close to the minimum, no the actual (non-integer) solution.
      //
      void triple_spin_search(const residue_spec_t &spec_first_residue,
			      const zo::rama_table_set &zorts,
			      bool apply_best_angles_flag,
			      int n_samples=60);

      // so triple_spin_search doesn't find the "correct" solution for tutorial residue A40.
      // Let's try to find a number of solutions.
      // We can refine all of them and then chose perhaps.
      // std::vector<scored_angle_set_t>
      // scored_angle_set_t is a score and a set of peptide rotation angles.
      //
      void find_maxima(const residue_spec_t &spec_first_residue,
		       const zo::rama_table_set &zorts,
		       unsigned int n_samples=60); // there are perhaps 50 maxima

      void run_optimizer(float start_angles[],
			 const coot::triple_crankshaft_set &tcs,
			 const zo::rama_table_set &zorts);

      // spin-search test the individual residues of input mol
      void test() const;

   };

}

#endif // CRANKSHAFT_HH
