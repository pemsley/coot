
#ifndef CRANKSHAFT_HH
#define CRANKSHAFT_HH

#include <vector>
#include <iostream>
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

      bool is_cis() const;

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
      void make_trans_from_non_pro_cis_if_needed(); // move the atoms of of res_1 and res_2 in mol
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
      triple_crankshaft_set() {}
      const crankshaft_set &operator[](unsigned int i) const { return cs[i]; }
      crankshaft_set &operator[](unsigned int i) { return cs[i]; } // for move_the_atoms()
      phi_psi_t phi_psi(unsigned int peptide_idx, float angle) const;
      std::pair<phi_psi_t, phi_psi_t> phi_psis_last(float ang_third) const;
      void move_the_atoms(float angles[]);
      const std::string &residue_type(unsigned int idx) const { return residue_types[idx]; }
      // the peptide_index is used index the residue types (i.e. starts at 1).
      float log_prob(const phi_psi_t &pp, unsigned int peptide_index, const zo::rama_table_set &zorts) const {
	 return zorts.value(pp, residue_types[peptide_index]);
      }
   };


   // now extended to an arbitrary number of residues
   //
   class nmer_crankshaft_set {

      std::vector<crankshaft_set> cs; // changed from array[3]
      std::vector<std::string> residue_types;

   public:
      nmer_crankshaft_set(const std::vector<mmdb::Residue *> residues_in,
			  const std::vector<std::string> &rts);
      // maybe we can save zorts in this class with a const ref?
      // above constructor will need changing however.
      nmer_crankshaft_set(const residue_spec_t &spec_first_residue,
			  const zo::rama_table_set &zorts,
			  mmdb::Manager *mol);
      nmer_crankshaft_set() {}
      const crankshaft_set &operator[](unsigned int i) const { return cs[i]; }
      crankshaft_set &operator[](unsigned int i) { return cs[i]; } // for move_the_atoms()
      phi_psi_t phi_psi(unsigned int peptide_idx, float angle) const;
      std::pair<phi_psi_t, phi_psi_t> phi_psis_last(float ang_third) const;
      void move_the_atoms(const std::vector<float> &angles_in);
      const std::string &residue_type(unsigned int idx) const { return residue_types[idx]; }
      // the peptide_index is used index the residue types (i.e. starts at 1).
      float log_prob(const phi_psi_t &pp, unsigned int peptide_index, const zo::rama_table_set &zorts) const {
	 return zorts.value(pp, residue_types[peptide_index]);
      }
      unsigned int size() const { return cs.size(); }
      unsigned int n_peptides() const { return size(); }
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

      class optimize_an_nmer {
      public:
	 static double f(const gsl_vector *v, void *params);
	 static void  df(const gsl_vector *v, void *params, gsl_vector *df);
	 static void fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
      };

      class triple_set_param_holder_t {
      public:
	 triple_set_param_holder_t(const zo::rama_table_set &rts_in,
				  const triple_crankshaft_set &tcs_in)
	    : zorts(rts_in), tcs(tcs_in) { }
	 const zo::rama_table_set &zorts;
	 const triple_crankshaft_set &tcs; // aha! const ref in constructed object
      };

      class nmer_set_param_holder_t {
      public:
	 nmer_set_param_holder_t(const zo::rama_table_set &rts_in,
				 const nmer_crankshaft_set &cs_in)
	    : zorts(rts_in), cs(cs_in) { }
	 const zo::rama_table_set &zorts;
	 const nmer_crankshaft_set &cs; // aha! const ref in constructed object
      };

      mmdb::Manager *mol;
      // caller ensures that res_1 and res_2 are valid
      std::vector<mmdb::Atom *> get_mainchain_atoms(mmdb::Residue *res_1, mmdb::Residue *res_2) const;
      mmdb::Atom *get_atom(mmdb::Residue *res_1, const std::string &atom_name) const;
   public:
      crankshaft(mmdb::Manager *mol_in) {
	 mol = new mmdb::Manager;
	 mol->Copy(mol_in, mmdb::MMDBFCM_All);
      }

      // a scored_angle_set_t needs to contain the info about the atoms
      // so that a scored_angle_set_t can be used to move the atoms
      //
      class scored_triple_angle_set_t : public triple_crankshaft_set {
      public:
	 scored_triple_angle_set_t() { minus_log_prob = 0; }
	 scored_triple_angle_set_t(const triple_crankshaft_set &tcs_in,
				   const std::vector<float> &angles_in, float lp) : triple_crankshaft_set(tcs_in), angles(angles_in), minus_log_prob(lp) {};
	 std::vector<float> angles;
	 float minus_log_prob;
	 bool is_close(const scored_triple_angle_set_t &sas_in) const {
	    float big_delta = clipper::Util::d2rad(5.0);
	    bool same = true;
	    for (std::size_t i=0; i<angles.size(); i++) {
	       if (std::abs(sas_in.angles[i] - angles[i]) > big_delta) {
		  same = false;
		  break;
	       }
	    }
	    return same;
	 }
	 bool operator<(const scored_triple_angle_set_t &sas_in) const {
	    return (minus_log_prob < sas_in.minus_log_prob);
	 }
	 bool filled() { return (angles.size() > 0); }
	 friend std::ostream &operator<<(std::ostream &s, const scored_triple_angle_set_t &r);
      };

      // a scored_angle_set_t needs to contain the info about the atoms
      // so that a scored_angle_set_t can be used to move the atoms
      //
      class scored_nmer_angle_set_t : public nmer_crankshaft_set {
      public:
	 scored_nmer_angle_set_t() { minus_log_prob = 0; }
	 scored_nmer_angle_set_t(const nmer_crankshaft_set &nmer_in,
				 const std::vector<float> &angles_in, float lp) : nmer_crankshaft_set(nmer_in), angles(angles_in), minus_log_prob(lp) {};
	 std::vector<float> angles;
	 float minus_log_prob;
	 bool is_close(const scored_triple_angle_set_t &sas_in) const {
	    float big_delta = clipper::Util::d2rad(5.0);
	    bool same = true;
	    for (std::size_t i=0; i<angles.size(); i++) {
	       if (std::abs(sas_in.angles[i] - angles[i]) > big_delta) {
		  same = false;
		  break;
	       }
	    }
	    return same;
	 }
	 bool operator<(const scored_triple_angle_set_t &sas_in) const {
	    return (minus_log_prob < sas_in.minus_log_prob);
	 }
	 bool filled() { return (angles.size() > 0); }
	 friend std::ostream &operator<<(std::ostream &s, const scored_triple_angle_set_t &r);
      };

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
      std::vector<scored_triple_angle_set_t>
      find_maxima(const residue_spec_t &spec_first_residue,
		  const zo::rama_table_set &zorts,
		  unsigned int n_samples=60); // there are perhaps 50 maxima

      scored_triple_angle_set_t run_optimizer(float start_angles[],
					      const coot::triple_crankshaft_set &tcs,
					      const zo::rama_table_set &zorts);

      // restores the atom positions in mol after write
      // sas is not const because we move the atoms (non-const of a crankshaft_set).
      void move_the_atoms_write_and_restore(scored_triple_angle_set_t sas, const std::string &pdb_file_name);

      // move the atoms, create a copy of mol, restore the atom positions
      mmdb::Manager *new_mol_with_moved_atoms(scored_triple_angle_set_t sas);

      // spin-search test the individual residues of input mol
      void test() const;

   };

   std::ostream &operator<<(std::ostream &s, const crankshaft::scored_triple_angle_set_t &r);
   
}

#endif // CRANKSHAFT_HH
