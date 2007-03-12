
#ifndef DUNBRACK_HH
#define DUNBRACK_HH

// #ifndef HAVE_STRING
// #define HAVE_STRING
// #include <string>
// #endif

// #ifndef HAVE_VECTOR
// #define HAVE_VECTOR
// #include <vector>
// #endif

// #include "mmdb_manager.h"

#include "chi-angles.hh"

namespace coot {


   class dunbrack : public chi_angles {

      float probability_limit; 
      const CMMDBManager *stored_mol;

      // void add_all_rotamers();  // an autogen function

      contact_info getcontacts(const atom_selection_container_t &asc) const;

      std::vector<CAtom *> ordered_residue_atoms(CResidue *residue_p) const;
      float d2rad(float degrees) const;
      std::vector<std::vector<std::string> >
      rotamer_atoms(const std::string &residue_name) const; 
      double chi_torsion(const std::vector<int> &chi_angle_atom_indices,
			 PCAtom *residue_atoms);
      double probability_score(double chi_angle, int ichi, const coot::simple_rotamer &rot);

      std::vector<coot::simple_rotamer>
      get_all_rotamers(const std::string &res_type) const;
      std::pair<short int, double> probability_of_this_rotamer(const std::vector<double> &chi_angles,
							       const std::vector<coot::simple_rotamer> &rots) const;
      std::vector<std::vector<int> > rotamer_atom_names_to_indices(const std::vector<std::vector<std::string> > &residue_rotamer_atoms, PCAtom *residue_atoms, int n_residue_atoms) const;

      short int similar_rotamer_chi(double target, double model) const {
	 short int is = 0;
	 double diff = target - model;
	 while (diff > 180.0)
	    diff -= 360.0;
	 while (diff < -180.0)
	    diff += 360.0;

	 if (fabs(diff) < 40.0)
	    is = 1;

	 return is;
      }

      // for penultime rotamer library
      static short int is_a_residue_name(const std::string &line_part);
      static short int end_of_a_rotamer_p(const std::vector<std::string> &parts);
      static std::string convert_residue_name(const std::string &name_in);
      simple_rotamer parse_prl_rotamer_line(const std::string &line, const std::vector<std::string> &line_parts);

   public:
      // We must be passed a deep copy of residue, the constructor
      // only copies the pointer.
      //
      // We only copy the mol pointer.  Why is it even needed you
      // might ask...
      // 
      // Well (sigh) it's needed in the calculation of the bonds which
      // uses mol->SeekContacts, even though we have a perfectly good
      // atom selection. There should be a version of SeekContacts
      // that does not need a CMMDBManager... Grumble grumble...
      //
      dunbrack(CResidue *residue,
	       CMMDBManager *mol,
	       float lowest_probability) :
	 chi_angles(residue, 0) {
	 probability_limit = lowest_probability;
	 stored_mol = mol;
      	 // add_all_rotamers();
	 // setup_chi_atom_pairs(); // in dunbrack at least,
		 		 // setup_chi_atom_pairs should happen
				 // after add_all_rotamers().
      }

      dunbrack(CResidue *residue,
	       CMMDBManager *mol,
	       float lowest_probability,
	       short int add_extra_PHE_and_TYR_rotamers_flag) :
	 chi_angles(residue, add_extra_PHE_and_TYR_rotamers_flag) {
	 probability_limit = lowest_probability;
	 stored_mol = mol;
      	 // add_all_rotamers();
	 // setup_chi_atom_pairs(); // in dunbrack at least,
		 		 // setup_chi_atom_pairs should happen
				 // after add_all_rotamers().
      }

            

      // For use with Z-score (which is analysis only: we don't move anything)
      // 
      dunbrack(CResidue *residue) : chi_angles(residue, 0) {}
      
      // Return NULL if no residues available for this residue type
      // 
      CResidue *GetResidue(int i_rot) const; // rotamer/button number
      std::vector<float> probabilities() const;
      std::vector<coot::simple_rotamer> rotamers(const std::string &res_type, float prob_cut) const; 

      void info() const;
      float Chi1(int i) const; // chi1 for the ith rotamer

      std::pair<short int, double> probability_of_this_rotamer(); // can't const - mmdb CResidue issues...

      //
      // LEU, VAL, THR have "nomenclature" (or real) chiral centres -
      // they are not dealt with here.
      // 
      // We deal with bifurcated symmetric non-chiral side chains (PHE, ASP,
      // GLU, THR)
      // 
      int optimize_rotamer_by_atom_names();

      // maybe this will need to be a static, or a constructor that
      // gets passed to some setup function... not sure yet.  Or maybe
      // it should return some sort of internal data.
      void read_penultimate_library(const std::string &filename);
   };
}

#endif // DUNBRACK_HH
