
#ifndef CRANKSHAFT_HH
#define CRANKSHAFT_HH

#include <vector>
#include <clipper/core/coords.h>

#include <geometry/residue-and-atom-specs.hh>

namespace coot {

   // consolidate with phi_psi_pair (i.e. make that disappear)
   class phi_psi_t {
      // and now with tau!
   public:
      phi_psi_t(float a, float b, float c) {
	 phi = a;
	 psi = b;
	 tau = c;
      }
      float phi;
      float psi;
      float tau;
   };

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


   public:
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
      std::pair<float, float> probability_of_spin_orientation(float a) const;
      
   };

   // this class does not do the right thing when used for residues with alt confs
   //
   class crankshaft {
      mmdb::Manager *mol;
      // caller ensures that res_1 and res_2 are valid
      std::vector<mmdb::Atom *> get_mainchain_atoms(mmdb::Residue *res_1, mmdb::Residue *res_2) const;
      mmdb::Atom *get_atom(mmdb::Residue *res_1, const std::string &atom_name) const;
   public:
      crankshaft(mmdb::Manager *mol_in) { mol = mol_in; }

      // pass the spec of the lower residue number (of a pair of residues in
      // the peptide)
      //
      std::vector<std::pair<float, float> > spin_search(const residue_spec_t &spec) const;

      // test the residues of input mol
      void test() const;

   };

}

#endif // CRANKSHAFT_HH
