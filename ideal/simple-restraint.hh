// -*-c++-*-
/* ideal/simple-resetraint.hh
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
 * Copyright 2008 by The University of Oxford
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

#ifndef HAVE_SIMPLE_RESTRAINT_HH
#define HAVE_SIMPLE_RESTRAINT_HH

#include <vector>
#include <string>
#include <stdexcept>

#include <mmdb2/mmdb_manager.h>
#include "coot-utils/bonded-pairs.hh"

#include "parallel-planes.hh"

#include "extra-restraints.hh"

// refinement_results_t is outside of the GSL test because it is
// needed to make the accept_reject_dialog, and that can be compiled
// without the GSL.
// 
namespace coot {

   enum { UNSET_INDEX = -1 };

   class refinement_lights_info_t {
   public:
      std::string name;   // e.g. "Bonds" or "Angles"
      std::string label;  // e.g. "Bonds:  6.543" 
      float value;        // e.g. 6.543
      refinement_lights_info_t(const std::string &name_in, const std::string label_in, float value_in) {
	 name = name_in;
	 label = label_in;
	 value = value_in;
      }
   };


   class rama_triple_t {
   public:
      mmdb::Residue *r_1; 
      mmdb::Residue *r_2; 
      mmdb::Residue *r_3;
      bool fixed_1; 
      bool fixed_2; 
      bool fixed_3; 
      rama_triple_t(mmdb::Residue *r1, mmdb::Residue *r2, mmdb::Residue *r3) {
	 r_1 = r1;
	 r_2 = r2;
	 r_3 = r3;
	 fixed_1 = 0;
	 fixed_2 = 0;
	 fixed_3 = 0;
      }
      rama_triple_t(mmdb::Residue *r1, mmdb::Residue *r2, mmdb::Residue *r3,
		    bool fixed_1_in, bool fixed_2_in, bool fixed_3_in) {
	 r_1 = r1;
	 r_2 = r2;
	 r_3 = r3;
	 fixed_1 = fixed_1_in;
	 fixed_2 = fixed_2_in;
	 fixed_3 = fixed_3_in;
      }
   };


   class distortion_torsion_gradients_t {
   public:
      double theta; // the torsion angle
      // x
      double dD_dxP1;
      double dD_dxP2;
      double dD_dxP3;
      double dD_dxP4;
      
      // y
      double dD_dyP1;
      double dD_dyP2;
      double dD_dyP3;
      double dD_dyP4;
      
      // z
      double dD_dzP1;
      double dD_dzP2;
      double dD_dzP3;
      double dD_dzP4;
   };

   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   //     class refinement_results_t, helper class for sending text
   //     results back to invoking function.  Returned by minimize()
   //     function.
   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   
   class refinement_results_t { 
   public:
      bool found_restraints_flag; // if we found restraints or not.
      int progress; // GSL_ENOPROG, GSL_CONTINUE, GSL_SUCCESS, GSL_ENOPROG (no progress)
      std::string info;
      std::vector<refinement_lights_info_t> lights;      
      refinement_results_t(bool frf, int prog_in,
			   const std::vector<refinement_lights_info_t> &lights_in) {
	 found_restraints_flag = frf;
	 info = ""; // not used
	 progress = prog_in;
	 lights = lights_in;
     }
      refinement_results_t(bool frf, int prog_in, const std::string &info_in) {
	 found_restraints_flag = frf;
	 info = info_in;
	 progress = prog_in;
      }
      refinement_results_t() {
	 info = "";
	 found_restraints_flag = 0;
      }
      refinement_results_t(const std::string &s_in) {
	 info = s_in;
	 found_restraints_flag = 0;
      }
   };
}

// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL  

#include <map>

#include "gsl/gsl_multimin.h"

// #include "Cartesian.h"
#include "coords/mmdb-extras.h" // for atom_selection_container_t, this and
			 // the interface that uses this can be
			 // deleted, I think, when we move to clipper
			 // (so that clipper does not get infected
			 // with
			 // MyMMDBManager/atom_selection_container_t
			 // stuff).

// for map stuff
#include "clipper/core/xmap.h"
#include "clipper/core/map_interp.h"

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh" // for atom_spec_t

// for protein dictionary container:
#include "geometry/protein-geometry.hh"

// For Kevin's (Log) Ramachandran Plot and derivatives
#include "lograma.h"

#ifndef DEGTORAD 
#define DEGTORAD 0.017453293
#endif
#ifndef RADTODEG
#define RADTODEG 57.29577793
#endif

namespace coot {

   
   // restraint types:
   // 
   enum restraint_type_t {BOND_RESTRAINT=1, ANGLE_RESTRAINT=2, TORSION_RESTRAINT=4, PLANE_RESTRAINT=8,
			  NON_BONDED_CONTACT_RESTRAINT=16, CHIRAL_VOLUME_RESTRAINT=32, RAMACHANDRAN_RESTRAINT=64,
			  START_POS_RESTRAINT=128,
			  TARGET_POS_RESTRANT=256, // restraint to make an atom be at a position
			  PARALLEL_PLANES_RESTRAINT=512,
			  GEMAN_MCCLURE_DISTANCE_RESTRAINT=1024};

   enum pseudo_restraint_bond_type {NO_PSEUDO_BONDS, HELIX_PSEUDO_BONDS,
				    STRAND_PSEUDO_BONDS};

   
   // restraints_usage_flag
   // 
   // use with restraints_usage_flag & 1 for BONDS or 
   //          restraints_usage_flag & 2 for ANGLES etc.      
   //
   enum restraint_usage_Flags { NO_GEOMETRY_RESTRAINTS = 0,
				BONDS = 1, BONDS_AND_ANGLES = 3,
				TORSIONS = 4,
				NON_BONDED = 16,
				CHIRAL_VOLUMES = 32,
				RAMA = 64, 
				BONDS_ANGLES_AND_TORSIONS = 7,
				BONDS_ANGLES_TORSIONS_AND_PLANES = 15,
				BONDS_AND_PLANES = 9,
				BONDS_ANGLES_AND_PLANES = 11, // no torsions
				BONDS_AND_NON_BONDED = 17,
				BONDS_ANGLES_AND_NON_BONDED = 19,
				BONDS_ANGLES_TORSIONS_AND_NON_BONDED = 23,
				BONDS_ANGLES_PLANES_AND_NON_BONDED = 27,
				BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS = 59,
				BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED = 31,
				BONDS_ANGLES_TORSIONS_PLANES_AND_CHIRALS = 47,
				BONDS_ANGLES_PLANES_AND_CHIRALS = 43,
				BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS = 63,
				
				BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA = 127,
				
				BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_PARALLEL_PLANES = 187,
				BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_PARALLEL_PLANES = 191,
				BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_RAMA_AND_PARALLEL_PLANES = 255,
				BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_GEMAN_MCCLURE_DISTANCES = 63+1024,
				TYPICAL_RESTRAINTS               = 1+2+  8+16+32+128+512,
				TYPICAL_RESTRAINTS_WITH_TORSIONS = 1+2+4+8+16+32+128+512,
				ALL_RESTRAINTS = 1+2+4+8+16+32+64+128+256+512+1024
				
   };

   enum peptide_restraints_usage_Flags { OMEGA_TORSION = 1,
					 OMEGA_PHI_PSI_TORSION = 3 };
					 
   enum { BONDS_MASK = 1,  ANGLES_MASK = 2, TORSIONS_MASK = 4, PLANES_MASK = 8, 
          NON_BONDED_MASK = 16,
	  CHIRAL_VOLUME_MASK = 32,
	  RAMA_PLOT_MASK = 64,
	  PARALLEL_PLANES_MASK = 256,
	  GEMAN_MCCLURE_DISTANCE_MASK = 512 };


   class ramachandran_restraint_flanking_residues_helper_t {
   public:
      int resno_first;
      int resno_third;
      std::vector<bool> is_fixed;
      ramachandran_restraint_flanking_residues_helper_t() {
	 is_fixed.resize(3,0);
      } 
   };


   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   //     class simple_restraint
   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   
   class simple_restraint {

   public: 
   
      int atom_index_1, atom_index_2, atom_index_3, atom_index_4, atom_index_5, atom_index_6;
      int atom_index_centre;
      // index and weight
      std::vector <std::pair<int, double> > plane_atom_index; // atom_index values can return negative (-1) for planes
      std::vector <std::pair<int, double> > atom_index_other_plane; // for the second plane in parallel planes
      double target_value; 
      double sigma; 
      float observed_value;    
      restraint_type_t restraint_type;
      int periodicity; 
      int chiral_volume_sign;
      double target_chiral_volume;
      int chiral_hydrogen_index; // if exactly one H attached to this chiral
                                 // centre, then the atom index, otherwise -1.
      // otherwise this is -1.
      std::vector<bool> fixed_atom_flags;
      std::vector<bool> fixed_atom_flags_other_plane;
      bool is_user_defined_restraint;
      // for mouse pull on an atom: this is where the user wants the atom to be
      //
      clipper::Coord_orth atom_pull_target_pos;

      // allocator for geometry_distortion_info_t
      simple_restraint() { is_user_defined_restraint = 0; }
      
      // Bond
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2, 
		       const std::vector<bool> &fixed_atom_flags_in,
		       float tar, 
		       float sig, float obs){
	 
	 restraint_type = rest_type; 
	 atom_index_1 = atom_1; 
	 atom_index_2 = atom_2;
	 observed_value = obs; 
	 sigma = sig; 
	 target_value = tar; 
	 fixed_atom_flags = fixed_atom_flags_in;
	 is_user_defined_restraint = false;
	 
	 // This finds a coding error
	 if (rest_type != BOND_RESTRAINT) { 
	    std::cout << "BOND ERROR" << std::endl; 
	    exit(1); 
	 }
      };

      // Geman-McClure distance (no obs)
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2, 
		       const std::vector<bool> &fixed_atom_flags_in,
		       float tar, float sig){
	 
	 restraint_type = rest_type; 
	 atom_index_1 = atom_1; 
	 atom_index_2 = atom_2;
	 observed_value = -1;  // not given or used
	 sigma = sig; 
	 target_value = tar; 
	 fixed_atom_flags = fixed_atom_flags_in;
	 is_user_defined_restraint = true;

	 // This finds a coding error
	 if (rest_type != restraint_type_t(GEMAN_MCCLURE_DISTANCE_MASK)) { 
	    std::cout << "BOND ERROR" << std::endl; 
	    exit(1); 
	 }
      };
      
    
      // Angle
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2, 
		       int atom_3, 
		       const std::vector<bool> &fixed_atom_flags_in,
		       float tar, 
		       float sig, float obs){
	 
	 restraint_type = rest_type; 
	 atom_index_1 = atom_1; 
	 atom_index_2 = atom_2;
	 atom_index_3 = atom_3;
	 observed_value = obs; 
	 sigma = sig; 
	 target_value = tar; 
	 fixed_atom_flags = fixed_atom_flags_in;
	 is_user_defined_restraint = 0;
	 if (rest_type != ANGLE_RESTRAINT) { 
	    std::cout << "ERROR::::: PROGRAM ERROR - ANGLE ERROR" << std::endl; 
	 }
      };

      // Torsion
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2, 
		       int atom_3, int atom_4, 
		       const std::vector<bool> &fixed_atom_flags_in,
		       float tar, 
		       float sig, float obs, int periodicity_in){

	 restraint_type = rest_type; 
	 atom_index_1 = atom_1; 
	 atom_index_2 = atom_2;
	 atom_index_3 = atom_3;
	 atom_index_4 = atom_4;
	 observed_value = obs; 
	 sigma = sig; 
	 target_value = tar; 
	 fixed_atom_flags = fixed_atom_flags_in;
	 periodicity = periodicity_in; 
	 is_user_defined_restraint = 0;
	 if (rest_type != TORSION_RESTRAINT) { 
	    std::cout << "ERROR::::: PROGRAM ERROR - TORSION ERROR" << std::endl; 
	 }
      }

      // Rama
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2, 
		       int atom_3, int atom_4, int atom_5, 
		       const std::vector<bool> &fixed_atom_flags_in) { 

	 restraint_type = rest_type; 
	 atom_index_1 = atom_1; 
	 atom_index_2 = atom_2;
	 atom_index_3 = atom_3;
	 atom_index_4 = atom_4;
	 atom_index_5 = atom_5;
	 fixed_atom_flags = fixed_atom_flags_in;
	 is_user_defined_restraint = 0;
	 if (rest_type != RAMACHANDRAN_RESTRAINT) { 
	    std::cout << "RAMACHANDRAN_RESTRAINT ERROR" << std::endl; 
	    exit(1); 
	 }
      }
      
      // Old Plane
      simple_restraint(restraint_type_t restraint_type_in,
		       const std::vector<int> &atom_index_in,
		       const std::vector<bool> &fixed_atom_flags_in,
		       float sig) {
	 
	 // Check restraint_type?
	 // Currently the only thing with an std::vector <int>
	 // atom_index_in is a plane restraint.  This could well
	 // change in the future.
	 // 
	 restraint_type = restraint_type_in;

	 plane_atom_index.resize(atom_index_in.size());
	 for (unsigned int i=0; i<atom_index_in.size(); i++)
	    plane_atom_index[i] = std::pair<int, double> (atom_index_in[i], sig);
	 
	 target_value = 0.0; // not needed for planes
	 sigma = sig;
	 fixed_atom_flags = fixed_atom_flags_in;
	 is_user_defined_restraint = 0;
      }

      // modern (atoms individually weighted) Plane
      // 
      simple_restraint(restraint_type_t restraint_type_in,
		       const std::vector<std::pair<int, double> > &atom_index_sigma_in,
		       const std::vector<bool> &fixed_atom_flags_in) {

	 //
	 // Check restraint_type?
	 // Currently the only thing with an std::vector <int>
	 // atom_index_in is a plane restraint.  This could well
	 // change in the future.
	 //
	 restraint_type = restraint_type_in; 

	 plane_atom_index = atom_index_sigma_in;
	 
	 target_value = 0.0; // not needed for planes
	 sigma = 0.02; // hack 
	 fixed_atom_flags = fixed_atom_flags_in;
	 is_user_defined_restraint = 0;
	 
      }
      

      // Parallel planes (actually angle-between-planes,typically-zero)
      simple_restraint(restraint_type_t restraint_type_in,
		       const std::vector<int> &atom_index_plane_1_in,
		       const std::vector<int> &atom_index_plane_2_in,
		       const std::vector<bool> &fixed_atom_flags_plane_1_in,
		       const std::vector<bool> &fixed_atom_flags_plane_2_in,
		       double target_angle_in,
		       double sigma_in) {
	 restraint_type = restraint_type_in; 
	 // atom_index = atom_index_plane_1_in;
	 // atom_index_other_plane = atom_index_plane_2_in;
	 plane_atom_index.resize(atom_index_plane_1_in.size());
	 atom_index_other_plane.resize(atom_index_plane_2_in.size());
	 for (unsigned int i=0; i<atom_index_plane_1_in.size(); i++)
	    plane_atom_index[i] = std::pair<int, double> (atom_index_plane_1_in[i], sigma_in);
	 for (unsigned int i=0; i<atom_index_plane_2_in.size(); i++)
	    atom_index_other_plane[i] = std::pair<int, double> (atom_index_plane_2_in[i], sigma_in);
	 
	 fixed_atom_flags             = fixed_atom_flags_plane_1_in;
	 fixed_atom_flags_other_plane = fixed_atom_flags_plane_2_in;
	 target_value = target_angle_in;
	 sigma = sigma_in;
	 is_user_defined_restraint = 1;
      } 

      // Non-bonded
      simple_restraint(restraint_type_t restraint_type_in, 
		       int index_1, 
		       int index_2,
		       const std::string &atom_1_type,
		       const std::string &atom_2_type,
		       const std::vector<bool> &fixed_atom_flags_in,
		       const protein_geometry &geom) { 
	 
	 if (restraint_type_in == NON_BONDED_CONTACT_RESTRAINT) { 
	    restraint_type = restraint_type_in;
	    atom_index_1 = index_1;
	    atom_index_2 = index_2;
	    // Enable shortening if there might be an H-bond between them.
	    std::pair<bool, double> nbc_dist = get_nbc_dist(atom_1_type, atom_2_type, geom); 
	    if (nbc_dist.first) {
	       target_value = nbc_dist.second;
	    } else {
	       // short/standard value
	       target_value = 2.5;
	    } 
	    sigma = 0.02;
	    fixed_atom_flags = fixed_atom_flags_in;
	    is_user_defined_restraint = 0;
	 } else { 
	    std::cout << "ERROR:: bad simple_restraint constructor usage "
		      << "- should be non-bonded\n";
	 } 
      }

      // Non-bonded v2 
      simple_restraint(restraint_type_t restraint_type_in, 
		       int index_1, 
		       int index_2,
		       const std::string &atom_1_type,
		       const std::string &atom_2_type,
		       const std::vector<bool> &fixed_atom_flags_in,
		       double dist_min) { 
	 
	 if (restraint_type_in == NON_BONDED_CONTACT_RESTRAINT) { 
	    restraint_type = restraint_type_in;
	    atom_index_1 = index_1;
	    atom_index_2 = index_2;
	    target_value = dist_min;
	    sigma = 0.02;
	    fixed_atom_flags = fixed_atom_flags_in;
	    is_user_defined_restraint = 0;
	 } else { 
	    std::cout << "ERROR:: bad simple_restraint constructor usage "
		      << "- should be non-bonded\n";
	 } 
      }

      // Chiral
      simple_restraint(restraint_type_t restraint_type_in, 
		       int atom_centre_idx_in,
		       int atom_idx_1_in,
		       int atom_idx_2_in, 
		       int atom_idx_3_in,
		       int volume_sign_in,
		       double target_volume_in,
		       double target_volume_sigma_in,
		       const std::vector<bool> &fixed_atom_flags_in,
		       int chiral_hydrogen_index_in) { 
	 
	 if (restraint_type_in == CHIRAL_VOLUME_RESTRAINT) { 
	    restraint_type = restraint_type_in;
	    atom_index_1 = atom_idx_1_in;
	    atom_index_2 = atom_idx_2_in;
	    atom_index_3 = atom_idx_3_in;
	    atom_index_centre = atom_centre_idx_in;
	    chiral_volume_sign = volume_sign_in;
	    target_chiral_volume = target_volume_in;
	    sigma = target_volume_sigma_in;
	    fixed_atom_flags = fixed_atom_flags_in;
	    chiral_hydrogen_index = chiral_hydrogen_index_in;
	    is_user_defined_restraint = 0;
	 } 
      }
      
      //start pos
      simple_restraint(restraint_type_t rest_type, int atom_1,
		       bool fixed_atom_flag_in,
		       float sig, float obs){
	 
	 restraint_type = rest_type; 
	 atom_index_1 = atom_1; 
	 observed_value = obs; 
	 sigma = sig; 
	 fixed_atom_flags = std::vector<bool> (1,fixed_atom_flag_in);
	 is_user_defined_restraint = 0;
	 
	 if (rest_type != START_POS_RESTRAINT) { 
	    std::cout << "START POS ERROR" << std::endl; 
	 }
      }

      // target_pos
      simple_restraint(restraint_type_t rest_type, int atom_idx, const clipper::Coord_orth &pos) {
	 restraint_type = rest_type;
	 atom_index_1 = atom_idx;
	 atom_pull_target_pos = pos;
	 if (rest_type != TARGET_POS_RESTRANT) {
	    std::cout << "TARGET POS ERROR" << std::endl; 
	 }
      }

      std::pair<bool, double> get_nbc_dist(const std::string &atom_1_type,
					   const std::string &atom_2_type, 
					   const protein_geometry &geom);

      double torsion_distortion(double model_torsion) const; 
      
      friend std::ostream &operator<<(std::ostream &s, const simple_restraint &r);
   };
   std::ostream &operator<<(std::ostream &s, const simple_restraint &r);
   bool target_position_eraser(const simple_restraint &r);

   class turn_off_when_close_target_position_restraint_eraser {
      int n_atoms;
      mmdb::PAtom *atoms;
      double close_dist;
   public:
      turn_off_when_close_target_position_restraint_eraser(mmdb::PAtom *atoms_in, int n_atoms_in) {
	 atoms = atoms_in;
	 n_atoms = n_atoms_in;
	 close_dist = 0.4;
      }
      bool operator() (const simple_restraint &r) const {
	 bool v = false;
	 if (r.restraint_type == restraint_type_t(TARGET_POS_RESTRANT)) { 
	    clipper::Coord_orth p_1 = co(atoms[r.atom_index_1]);
	    double d = sqrt((p_1-r.atom_pull_target_pos).lengthsq());
	    std::cout << "turn_off_when_close_target_position_restraint_eraser: dist " << d
		      << std::endl;
	    if (d < close_dist) v = true;
	 }
	 return v;
      }
   };

   

   // We need something to quickly convert between atom name,
   // sequence number, chain id to index into the atom selection
   // array (get_asc_index)
   //
   // Let's use an associative array to do that (map)
   // 
      


   // so that we can use distortion_score_plane_internal for both the
   // distortion score and in the generation of the plane gradients:
   // 
   class plane_distortion_info_t {
   public:
      std::vector<double> abcd;
      clipper::Coord_orth centre_1;
      clipper::Coord_orth centre_2;
      plane_distortion_info_t() {
	 distortion_score = 0.0;
      }
      double distortion_score;
   };

   class geometry_distortion_info_t {
   public:
      geometry_distortion_info_t(double distortion_in,
				 const simple_restraint &rest_in,
				 residue_spec_t &residue_spec_in) {
	 distortion_score = distortion_in;
	 restraint = rest_in;
	 residue_spec = residue_spec_in;
	 set = 1;
      }
      geometry_distortion_info_t() {
	 set = 0;
      }
      bool set;
      double distortion_score;
      simple_restraint restraint;
      std::vector<int> atom_indices;
      residue_spec_t residue_spec;
      friend std::ostream &operator<<(std::ostream &s, geometry_distortion_info_t);

      // This comparison can throw an exception, (when either of the
      // arguments is not initialised)
      // 
      bool operator<(const geometry_distortion_info_t &gdi) const {
	 if (! gdi.initialised_p())
	    throw std::runtime_error("unitialised passed geometry_distortion_info_t");
	 if (! initialised_p())
	    throw std::runtime_error("unitialised this geometry_distortion_info_t");
	 return (distortion_score < gdi.distortion_score);
      }
      bool operator>(const geometry_distortion_info_t &gdi) const {
	 if (! gdi.initialised_p())
	    throw std::runtime_error("unitialised passed geometry_distortion_info_t");
	 if (! initialised_p())
	    throw std::runtime_error("unitialised this geometry_distortion_info_t");
	 return (distortion_score > gdi.distortion_score);

      }
      bool initialised_p() const { return set; }
   };
   std::ostream &operator<<(std::ostream &s, geometry_distortion_info_t);

   class geometry_distortion_info_container_t {
   public:
      std::string chain_id;
      std::vector<geometry_distortion_info_t> geometry_distortion;
      mmdb::PAtom *atom;
      int n_atoms;
      int min_resno; 
      int max_resno;
      geometry_distortion_info_container_t(const std::vector<geometry_distortion_info_t> &geometry_distortion_in, 
					   mmdb::PAtom *atom_in, 
					   const std::string &chain_id_in) { 
	 geometry_distortion = geometry_distortion_in;
	 atom = atom_in;
	 chain_id = chain_id_in;
      }
      geometry_distortion_info_container_t(mmdb::PAtom *atom_in, int n_atoms_in, 
					   const std::string &chain_id_in) { 
	 atom = atom_in;
	 n_atoms = n_atoms_in;
	 chain_id = chain_id_in;
      }
      void set_min_max(int min_resno_in, int max_resno_in) { 
	 min_resno = min_resno_in;
	 max_resno = max_resno_in;
      }
      int size () const { return geometry_distortion.size(); }
      double print() const;  // return the total distortion
      friend std::ostream &operator<<(std::ostream &s, geometry_distortion_info_container_t);
   };
   std::ostream &operator<<(std::ostream &s, geometry_distortion_info_container_t gdic);

   class omega_distortion_info_t {
   public:
      int resno;
      double distortion;
      std::string info_string;
      omega_distortion_info_t(int resno_in, double distortion_in, const std::string &s) {
	 resno = resno_in;
	 distortion = distortion_in;
	 info_string = s;
      } 
   };
   
   class omega_distortion_info_container_t {
   public:
      std::string chain_id;
      std::vector<omega_distortion_info_t> omega_distortions; // in degrees away from 180
      int min_resno;
      int max_resno;
      omega_distortion_info_container_t(const std::string &chain_id_in, int min_resno_in, int max_resno_in) {
	 chain_id = chain_id_in;
	 min_resno = min_resno_in;
	 max_resno = max_resno_in;
      }
   };


   double distortion_score(const gsl_vector *v, void *params); 
   double distortion_score_bond(const simple_restraint &bond_restraint,
				const gsl_vector *v); 
   double distortion_score_geman_mcclure_distance(const simple_restraint &bond_restraint,
						  const gsl_vector *v,
						  const double &alpha); 
   double distortion_score_angle(const simple_restraint &angle_restraint,
				 const gsl_vector *v);
   // torsion score can throw a std::runtime_error if there is a problem calculating the torsion.
   double distortion_score_torsion(const simple_restraint &torsion_restraint,
				    const gsl_vector *v); 
   double distortion_score_plane(const simple_restraint &plane_restraint,
				  const gsl_vector *v); 
   double distortion_score_chiral_volume(const simple_restraint &chiral_restraint,
					 const gsl_vector *v); 
   double distortion_score_rama(const simple_restraint &chiral_restraint,
				const gsl_vector *v,
				const LogRamachandran &lograma);
   double distortion_score_start_pos(const simple_restraint &start_pos_restraint,
			    void *params,
			    const gsl_vector *v);
   double distortion_score_target_pos(const simple_restraint &start_pos_restraint,
				      void *params,
				      const gsl_vector *v);
   double distortion_score_non_bonded_contact(const simple_restraint &plane_restraint,
					      const gsl_vector *v);
   double distortion_score_parallel_planes(const simple_restraint &plane_restraint,
					   const gsl_vector *v); 
   void fix_chiral_atom_maybe (const simple_restraint &chiral_restraint,
			       gsl_vector *v);
   void fix_chiral_atom_internal (const simple_restraint &chiral_restraint,
				  gsl_vector *v);
   
   plane_distortion_info_t
   distortion_score_plane_internal(const simple_restraint &plane_restraint,
				   const gsl_vector *v);
   plane_distortion_info_t
   distortion_score_2_planes(const std::vector<std::pair<int, double> > &atom_index_set_1,
			     const std::vector<std::pair<int, double> > &atom_index_set_2,
			     const double &restraint_sigma,
			     const gsl_vector *v);


   distortion_torsion_gradients_t
   fill_distortion_torsion_gradients(const clipper::Coord_orth &P1,
				     const clipper::Coord_orth &P2,
				     const clipper::Coord_orth &P3,
				     const clipper::Coord_orth &P4);
   /* The gradients of f, df = (df/dx(k), df/dy(k) .. df/dx(l) .. ). */
   void my_df (const gsl_vector *v, void *params, gsl_vector *df); 
   // just the bond terms: 
   void my_df_bonds(const gsl_vector *v, void *params, gsl_vector *df); 
   // just the bond terms: 
   void my_df_geman_mcclure_distances(const gsl_vector *v, void *params, gsl_vector *df); 
   // just the angle terms: 
   void my_df_angles(const gsl_vector *v, void *params, gsl_vector *df); 
   //  just the torsion terms:
   void my_df_torsions(const gsl_vector *v, void *params, gsl_vector *df);
   void my_df_torsions_internal(const gsl_vector *v, void *params, gsl_vector *df, bool do_rama_torsions);
   //  just the ramachandran plot gradient terms:
   void my_df_rama(const gsl_vector *v, void *params, gsl_vector *df); 
   //  the plane deviation from terms:
   void my_df_planes(const gsl_vector *v, void *params, gsl_vector *df); 
   //  the non-bonded contacts
   void my_df_non_bonded(const gsl_vector *v, void *params, gsl_vector *df); 
   //  the chiral volumes
   void my_df_chiral_vol(const gsl_vector *v, void *params, gsl_vector *df); 
   //  the deviation from starting point terms:
   void my_df_start_pos(const gsl_vector *v, void *params, gsl_vector *df); 
   //  the deviation from atom pull point 
   void my_df_target_pos(const gsl_vector *v, void *params, gsl_vector *df); 
   //  20131012 the parallel plane deviation from terms:
   void my_df_parallel_planes(const gsl_vector *v, void *params, gsl_vector *df); 
   // Compute both f and df together.
   void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df);

   // replace this function, to test if things go faster with
   // alternative implementations?
   // 
   inline double f_inv_fsqrt(double v) {
      //
      return 1.0/sqrt(v);
   } 
   
   // debugging function
   void 
   numerical_gradients(gsl_vector *v, void *params, gsl_vector *df);
   
   // Getting adventurous:
   // 
   double electron_density_score(const gsl_vector *v, void *params);
   // new style Grad_map/Grad_orth method
   void my_df_electron_density (const gsl_vector *v, void *params, gsl_vector *df);
   // old style numerical method
   void my_df_electron_density_old (gsl_vector *v, void *params, gsl_vector *df); 


   // non-refinement function: just checking geometry:
   // 
   // return a vector because there could be many alt confs to this
   // residue, only one of which is wrong.  return a pair: the first
   // part is the good/bad flag the second tells us what the altconf
   // of the (bad) chiral atom was.
   // 
   std::vector<std::pair<short int, atom_spec_t> >
   is_inverted_chiral_atom_p(const dict_chiral_restraint_t &chiral_restraint,
			mmdb::Residue *res);

   std::pair<std::vector<std::string> , std::vector <atom_spec_t> >
   inverted_chiral_volumes(mmdb::Manager *mol, protein_geometry *geom_p,
			   int cif_dictionary_read_number);
      
   

   // -------------------------------------------------------------------------
   // -------------------------------------------------------------------------
   //                      restraints_container_t
   // -------------------------------------------------------------------------
   // -------------------------------------------------------------------------

   class restraints_container_t {

      enum reporting_level_t { QUIET, NORMAL, VERBOSE };

   public:

      class restraint_counts_t {
      public:
	 int n_bond_restraints; 
	 int n_angle_restraints; 
	 int n_plane_restraints; 
	 int n_chiral_restr;
	 int n_torsion_restr;
	 restraint_counts_t() {
	    n_bond_restraints = 0;
	    n_angle_restraints = 0; 
	    n_plane_restraints =0; 
	    n_chiral_restr = 0;
	    n_torsion_restr = 0;
	 }
	 void operator+=(const restraint_counts_t &r) {
	    n_bond_restraints += r.n_bond_restraints;
	    n_angle_restraints += r.n_angle_restraints;
	    n_plane_restraints += r.n_plane_restraints;
	    n_chiral_restr += r.n_chiral_restr;
	    n_torsion_restr += r.n_torsion_restr;
	 }
	 void report(bool do_residue_internal_torsions) {
	    std::cout << "created " << n_bond_restraints   << " bond       restraints " << std::endl;
	    std::cout << "created " << n_angle_restraints  << " angle      restraints " << std::endl;
	    std::cout << "created " << n_plane_restraints  << " plane      restraints " << std::endl;
	    std::cout << "created " << n_chiral_restr << " chiral vol restraints " << std::endl;
	    if (do_residue_internal_torsions)
	       std::cout << "created " << n_torsion_restr << " torsion restraints " << std::endl;
	 } 
      };

   private:

      std::vector<simple_restraint> restraints_vec; 
      mmdb::PPAtom atom;
      bool from_residue_vector;
      int SelHnd_atom; // the selection handle for the atom array.
		       // Note to self: when restraints_container_t
		       // goes out of scope, we should do a
		       // mol->DeleteSelection(SelHnd_atom).
      // atom_selection_container_t asc;
      double *par; 
      int n_atoms;
      gsl_vector *x; // these are the variables, x_k, y_k, z_k, x_l etc.
      bool are_all_one_atom_residues; 
      mmdb::Manager *mol;
      
      // The bool is the "atoms of this residue are fixed" flag.
      std::vector<std::pair<bool,mmdb::Residue *> > residues_vec;
      int udd_bond_angle;  // for is a bond, angle or not (0).
      int udd_atom_index_handle; // for indexing into the atoms array.

      // and the new (7Nov2003) residue selection (used in non-bonded stuff)
      //
      mmdb::PPResidue SelResidue_active;
      int nSelResidues_active;

      // using residue_vector
      bool is_a_moving_residue_p(mmdb::Residue *r) const;

      void filter_non_bonded_by_distance(const std::vector<std::vector<int> > &non_bonded_atom_indices, double dist);
      // don't include bonds which are make to flanking atoms/residues
      void construct_non_bonded_contact_list(const bonded_pair_container_t &bpc,
					     const protein_geometry &geom); // fills filtered_non_bonded_atom_indices;
      // which uses
      void construct_non_bonded_contact_list_by_res_vec(const bonded_pair_container_t &bpc,
							const protein_geometry &geom);
      // and
      void construct_nbc_for_moving_non_moving_bonded(unsigned int i, unsigned int j, const std::string &link_type, const protein_geometry &geom);
      // (which modifies (adds to) filtered_non_bonded_atom_indices).

      // or 
      void construct_non_bonded_contact_list_conventional();
      // to make
      std::vector<std::vector<int> > filtered_non_bonded_atom_indices;
      

      // we are given these at the constructor and they are needed in
      // make_restraints():
      int istart_res;
      int iend_res;  // inclusive.
      
      short int istart_minus_flag; // were there flanking residues in the molecule?
      short int iend_plus_flag;    // (set in init_from_mol);
      // likewise for the chain too
      std::string chain_id_save;

      // Here we store the pointers to residues that flank the
      // moving_atom_selection.  We need them to provide fixed points
      // from which we generate restraints, not all the atoms of which
      // are variables.
      // 
      mmdb::PResidue previous_residue;
      mmdb::PResidue next_residue;
      
      reporting_level_t verbose_geometry_reporting;
   
      // we will store in here the initial positions as parameters
      // 
      std::vector<double> initial_position_params_vec;

      // print chi_squared values (after refinement)
      // return a string that can be added into a dialog.
      std::vector<refinement_lights_info_t>
      chi_squareds(std::string title, const gsl_vector *v) const;

      // all the alt confs should either be the same as each other or ""
      // 
      short int check_altconfs_for_plane_restraint(const std::vector<std::string> &altconfs) const {
	 short int stat = 1;
	 std::string A("");

	 for (unsigned int iat=0; iat<altconfs.size(); iat++) {
	    if (altconfs[iat] != A) {
	       A = altconfs[iat];
	       break;
	    }
	 }
	 if (A == "") // there are no alt confs
	    return 1;  // OK

	 for (unsigned int iat=0; iat<altconfs.size(); iat++) {
	    if ((altconfs[iat] != "") && (altconfs[iat] != A)) {
	       stat = 0;  // mixed altconfs, badness.
	       break;
	    }
	 }
	 return stat;
      }

      // because we may want to refine just the bonds at the beginning
      // and then introduce angles and other stuff after things have
      // settled down
      // 

      // needs to be public for gsl_multimin_fdfminimizer_set()
      // (but we could make that used in a class function...)
      // 
      gsl_multimin_function_fdf multimin_func; 

      short int include_map_terms_flag;

      LogRamachandran lograma;

      // internal function, most of the job of the constructor:
      void init_from_mol(int istart_res_in, int iend_res_in,
			 short int have_flanking_residue_at_start,
			 short int have_flanking_residue_at_end,
			 short int have_disulfide_residues,
			 const std::string &altloc,
			 const char *chain_id,
			 mmdb::Manager *mol_in, 
			 const std::vector<atom_spec_t> &fixed_atom_specs);

      void init_from_residue_vec(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
				 const protein_geometry &geom,
				 mmdb::Manager *mol_in,
				 const std::vector<atom_spec_t> &fixed_atom_specs);

      void init_shared_pre(mmdb::Manager *mol_in);

      void init_shared_post(const std::vector<atom_spec_t> &fixed_atom_specs);
      // neighbour residues already are fixed.
      void add_fixed_atoms_from_flanking_residues(const bonded_pair_container_t &bpc);

   
      // 
      clipper::Xmap<float> map; 
      double map_weight; 

      void add(restraint_type_t rest_type, int atom_1, int atom_2, 
	       const std::vector<bool> &fixed_atom_flags,
	       float tar, 
	       float sig, float obs){

	 if (sig > 0.0) { 
	    simple_restraint r(rest_type, atom_1, atom_2, fixed_atom_flags, tar, sig, obs);
	    restraints_vec.push_back(r);
	 }
      }

      bool add(restraint_type_t rest_type, int atom_1, int atom_2, int atom_3, 
	       const std::vector<bool> &fixed_atom_flags,
	       float tar, 
	       float sig, float obs){
         
	 bool r = 0;
	 if (sig > 0.0) { 
	    restraints_vec.push_back(simple_restraint(rest_type, atom_1, atom_2, 
						      atom_3,
						      fixed_atom_flags,
						      tar, sig, obs));
	    r = 1;
	 }
	 return r;
      }

      bool add(restraint_type_t rest_type, int atom_1, int atom_2, 
	       int atom_3, int atom_4,
	       const std::vector<bool> &fixed_atom_flags,
	       float tar, 
	       float sig, float obs, int periodicty){

	 bool r = 0;
	 if (sig > 0.0) { 
	    restraints_vec.push_back(simple_restraint(rest_type, atom_1, atom_2, 
						      atom_3, atom_4,
						      fixed_atom_flags,
						      tar, sig, obs, periodicty));
	    r = 1;
	 }
	 return r;
      }

      void add_user_defined_torsion_restraint(restraint_type_t rest_type, int atom_1, int atom_2, 
					      int atom_3, int atom_4,
					      const std::vector<bool> &fixed_atom_flags,
					      float tar, 
					      float sig, float obs, int periodicty) {
	 bool r = add(rest_type, atom_1, atom_2, atom_3, atom_4,
		      fixed_atom_flags, tar, sig, obs, periodicty);
	 if (r) {
	    // bleugh.
	    restraints_vec.back().is_user_defined_restraint = 1;
	 }
      }
      
      void add_user_defined_angle_restraint(restraint_type_t rest_type, int atom_1, int atom_2, 
					    int atom_3,
					    const std::vector<bool> &fixed_atom_flags,
					    float tar, 
					    float sig, float obs) {
	 bool r = add(rest_type, atom_1, atom_2, atom_3,
		      fixed_atom_flags, tar, sig, obs);
	 if (r) {
	    // bleugh.
	    restraints_vec.back().is_user_defined_restraint = 1;
	 }
      }

      

      // used for Ramachandran restraint
      void add(restraint_type_t rest_type,
	       int atom_1, int atom_2, int atom_3, 
	       int atom_4, int atom_5, 
	       const std::vector<bool> &fixed_atom_flag){
    
	 restraints_vec.push_back(simple_restraint(rest_type,
						   atom_1, atom_2, atom_3,
						   atom_4, atom_5, 
						   fixed_atom_flag));
      }

      
      void add_non_bonded(int index1, int index2,
			  const std::string &atom_type_1, 
			  const std::string &atom_type_2, 
			  const std::vector<bool> &fixed_atom_flag,
			  const protein_geometry &geom) { 
	 restraints_vec.push_back(simple_restraint(NON_BONDED_CONTACT_RESTRAINT,
						   index1, index2,
						   atom_type_1, atom_type_2,
						   fixed_atom_flag, geom));
      } 

      void add_non_bonded(int index1, int index2,
			  const std::string &atom_type_1, 
			  const std::string &atom_type_2, 
			  const std::vector<bool> &fixed_atom_flag,
			  double dist_min) { 
	 restraints_vec.push_back(simple_restraint(NON_BONDED_CONTACT_RESTRAINT,
						   index1, index2,
						   atom_type_1, atom_type_2,
						   fixed_atom_flag, dist_min));
      }

      void add_target_position_restraint(int idx, clipper::Coord_orth &target_pos) {
	 simple_restraint r(TARGET_POS_RESTRANT, idx, target_pos);
	 restraints_vec.push_back(r);
      }

      // construct a restraint and add it to restraints_vec
      //
      // this assumes the sigmas in atom_index_sigma_in are sensible - so the calling 
      // routines needs to make sure that this is the case.
      void add_plane(const std::vector<std::pair<int, double> > atom_index_sigma_in,
		     const std::vector<bool> &fixed_atom_flags) {
	 restraints_vec.push_back(simple_restraint(PLANE_RESTRAINT,
						   atom_index_sigma_in, 
						   fixed_atom_flags));
      }
      
      //used for start pos restraints
      bool add(restraint_type_t rest_type, int atom_1,
	       bool fixed_atom_flag,
	       float sig, float obs){

	 bool r = 0;
	 if (sig > 0.0) { 
	    restraints_vec.push_back(simple_restraint(rest_type, atom_1, 
						      fixed_atom_flag,
						      sig, obs));
	    r = 1;
	 }
	 return r;
      }
      
      void add_user_defined_start_pos_restraint(restraint_type_t rest_type, int atom_1,
						bool fixed_atom_flag, float sig, float obs) {
	 bool r = add(rest_type, atom_1, fixed_atom_flag, sig, obs);
	 if (r) {
	    restraints_vec.back().is_user_defined_restraint = 1;
	 }
      }
      
      void add_geman_mcclure_distance(restraint_type_t rest_type,
				      int atom_1, int atom_2, 
				      const std::vector<bool> &fixed_atom_flags,
				      float tar, float sig) { 

	 if (sig > 0.0) { 
	    simple_restraint r(rest_type, atom_1, atom_2, fixed_atom_flags, tar, sig);
	    restraints_vec.push_back(r);
	 }
      }


      // atom indexing stuff
      //
      void setup_asc_indexing(); 
      
      int get_asc_index_old(const std::string &at_name,
			int resno,
			const char *chain_id) const; 

      int get_asc_index_new(const char *at_name,
			    const char *alt_loc,
			    int resno,
			    const char *ins_code,
			    const char *chain_id) const; 

      int get_asc_index(const char *at_name,
			const char *alt_loc,
			int resno,
			const char *ins_code,
			const char *chain_id) const;

      int get_asc_index(const atom_spec_t &spec) const;

      int add_bonds(int idr, mmdb::PPAtom res_selection,
		    int i_no_res_atoms,
		    mmdb::PResidue SelRes,
		    const protein_geometry &geom); 

      int add_angles(int idr, mmdb::PPAtom res_selection,
		     int i_no_res_atoms,
		     mmdb::PResidue SelRes,
		     const protein_geometry &geom); 
   
      int add_torsions(int idr, mmdb::PPAtom res_selection,
		       int i_no_res_atoms,
		       mmdb::PResidue SelRes,
		       const protein_geometry &geom);

      int add_chirals(int idr, mmdb::PPAtom res_selection,
		      int i_no_res_atoms,
		      mmdb::PResidue SelRes,
		      const protein_geometry &geom);

      int add_planes  (int idr, mmdb::PPAtom res_selection,
		       int i_no_res_atoms,
		       mmdb::PResidue SelRes,
		       const protein_geometry &geom);

      restraint_counts_t 
      apply_mods(int idr, mmdb::PPAtom res_selection,
		 int i_no_res_atoms,
		 mmdb::PResidue SelRes,
		 const protein_geometry &geom);

      void
      apply_mod(const std::string &mod_name,
		const protein_geometry &geom,
		int idr,
		mmdb::PResidue residue_p);

      void apply_mod_bond(const chem_mod_bond &mod_bond,
			  mmdb::PResidue residue_p);

      void apply_mod_angle(const chem_mod_angle &mod_angle,
			  mmdb::PResidue residue_p);

      void apply_mod_plane(const chem_mod_plane &mod_plane,
			   mmdb::PResidue residue_p);

      void mod_bond_add(const chem_mod_bond &mod_bond,
			mmdb::PResidue residue_p);

      void mod_bond_change(const chem_mod_bond &mod_bond,
			   mmdb::PResidue residue_p);

      void mod_bond_delete(const chem_mod_bond &mod_bond,
			   mmdb::PResidue residue_p);

      void mod_angle_add(const chem_mod_angle &mod_angle,
			mmdb::PResidue residue_p);

      void mod_angle_change(const chem_mod_angle &mod_angle,
			   mmdb::PResidue residue_p);

      void mod_angle_delete(const chem_mod_angle &mod_angle,
			   mmdb::PResidue residue_p);

      void mod_plane_add(const chem_mod_plane &mod_plane,
			 mmdb::PResidue residue_p);

      void mod_plane_delete(const chem_mod_plane &mod_plane,
			    mmdb::PResidue residue_p);
      
      bool dictionary_name_matches_coords_resname(const std::string &comp_id,
						  const std::string &resname) const {
	 std::string r = resname;
	 if (r.length() > 2) 
	    if (r[2] == ' ')
	       r = resname.substr(0,2);
	 return (r == comp_id);
      }

      int make_link_restraints          (const protein_geometry &geom,
					 bool do_rama_plot_retraints);

      // which uses either:
      int make_link_restraints_by_linear(const protein_geometry &geom,
					 bool do_rama_plot_retraints);

      bonded_pair_container_t
      make_link_restraints_from_res_vec(const protein_geometry &geom,
					    bool do_rama_plot_retraints);
      // both of which use:
      int make_link_restraints_by_pairs(const protein_geometry &geom,
					const bonded_pair_container_t &bonded_residue_pairs,
					std::string link_or_flanking_link_string);

      void add_rama_links(int SelHnd, const protein_geometry &geom);
					
      void add_rama_links_from_res_vec(const bonded_pair_container_t &bonded_residue_pairs,
				       const protein_geometry &geom);
					

      int make_monomer_restraints       (const protein_geometry &geom,
					 short int do_residue_internal_torsions);
      //uses:
      int make_monomer_restraints_by_linear (const protein_geometry &geom,
					     bool do_residue_internal_torsions);
      int make_monomer_restraints_from_res_vec (const protein_geometry &geom,
						bool do_residue_internal_torsions);
      // which use:
      restraint_counts_t make_monomer_restraints_by_residue(mmdb::Residue *residue_p,
							    const protein_geometry &geom,
							    bool do_residue_internal_torsions);
      
      bonded_pair_container_t
      make_flanking_atoms_restraints(const protein_geometry &geom,
				     bool do_rama_plot_retraints); // no torsions
      // uses the following
      bonded_pair_container_t
      bonded_flanking_residues(const protein_geometry &geom) const;
   
      // new flanking residue search
      bonded_pair_container_t bonded_flanking_residues_by_residue_vector(const protein_geometry &geom) const;
      // old style linear search (n +/- 1) selection for flanking residues
      bonded_pair_container_t bonded_flanking_residues_by_linear(const protein_geometry &geom) const;
      // find residues in the neighbourhood that are not in the refining set
      // and are not already marked as bonded flankers.
      // 
      std::vector<mmdb::Residue *> non_bonded_neighbour_residues;
      // set by this function:
      void set_non_bonded_neighbour_residues_by_residue_vector(const bonded_pair_container_t &bonded_flanking_pairs, const protein_geometry &geom);

      
      int make_flanking_atoms_rama_restraints(const protein_geometry &geom);

      // return a container of all the bonded residues (as pairs) from
      // the given atom selection
      bonded_pair_container_t bonded_residues_conventional(int SelResHnd,
							    const protein_geometry &geom) const;

      // this is to make the make_link_restraints work as it used to,
      // simply by bonding residues that are next to each other by
      // seqNum or residue index.
      bonded_pair_container_t bonded_residues_by_linear(int SelResHnd,
							 const protein_geometry &geom) const;
      std::vector<rama_triple_t> make_rama_triples(int SelResHnd,
						   const protein_geometry &geom) const;
      std::pair<bool,float> closest_approach(mmdb::Residue *r1, mmdb::Residue *r2) const;

      // find simple (tandem residue) links (based on residue-name and
      // restraints group type).
      // 
      std::string find_link_type(mmdb::Residue *first,
				 mmdb::Residue *second,
				 const protein_geometry &geom) const;

      // a pair, first is if C and N are close and second if and order
      // switch is needed to make it so.
      // 
      std::pair<bool, bool> peptide_C_and_N_are_close_p(mmdb::Residue *r1, mmdb::Residue *r2) const;


      // a pair, first is if C and N are close and second if and order
      // switch is needed to make it so.
      //
      // return "" as first if no close link found.
      // 
      std::pair<std::string, bool> general_link_find_close_link(const std::vector<std::pair<chem_link, bool> > &li,
								mmdb::Residue *r1, mmdb::Residue *r2,
								bool order_switch_flag,
								const protein_geometry &geom) const;

      std::string general_link_find_close_link_inner(const std::vector<std::pair<chem_link, bool> > &li,
						     mmdb::Residue *r1, mmdb::Residue *r2,
						     bool order_switch_flag,
						     const protein_geometry &geom) const; 
      
      
      void make_helix_pseudo_bond_restraints();
      void make_strand_pseudo_bond_restraints();

      bool link_infos_are_glycosidic_p(const std::vector<std::pair<chem_link, bool> > &link_infos) const;

      // return "" on failure to find link
      std::string find_glycosidic_linkage_type(mmdb::Residue *first, mmdb::Residue *second,
					       const protein_geometry &geom,
					       bool use_links_in_molecule) const;
   
      int add_link_bond(std::string link_type,
			mmdb::PResidue first, mmdb::PResidue second,
			short int is_fixed_first_res,
			short int is_fixed_second_res,
			const protein_geometry &geom);

      int add_link_angle(std::string link_type,
			 mmdb::PResidue first, mmdb::PResidue second,
			 short int is_fixed_first_res,
			 short int is_fixed_second_res,
			 const protein_geometry &geom);

      int add_link_torsion(std::string link_type,
			   int phi_psi_restraints_type,
			   mmdb::PResidue first, mmdb::PResidue second,
			   short int is_fixed_first, short int is_fixed_second,
			   const protein_geometry &geom);

      int add_rama(std::string link_type,
		   mmdb::PResidue prev,
		   mmdb::PResidue this_res,
		   mmdb::PResidue post,
		   bool is_fixed_first_res,
		   bool is_fixed_second_res,
		   bool is_fixed_third_res,
		   const protein_geometry &geom);

      int add_link_plane(std::string link_type,
			 mmdb::PResidue first, mmdb::PResidue second,
			 short int is_fixed_first_res,
			 short int is_fixed_second_res,
			 const protein_geometry &geom);

      int add_parallel_planes(mmdb::Residue *first, mmdb::Residue *second,
			      short int is_fixed_first_res,
			      short int is_fixed_second_res,
			      const protein_geometry &geom);

      int add_link_plane_tmp(std::string link_type,
			 mmdb::PResidue first, mmdb::PResidue second,
			 short int is_fixed_first_res,
			 short int is_fixed_second_res,
			 const protein_geometry &geom);
      void symmetry_non_bonded_contacts(bool p);
      std::vector<std::vector<int> > bonded_atom_indices;

      class reduced_angle_info_container_t {
      public:
	 reduced_angle_info_container_t(const std::vector<simple_restraint> &r);
	 std::map<int, std::vector<std::pair<int, int> > > angles;
	 bool is_1_4(int i, int j) const;
	 void write_angles_map(const std::string &file_name) const;
      };
      bool check_for_1_4_relation(int i, int j) const;
      bool check_for_1_4_relation(int i, int j, const reduced_angle_info_container_t &ai) const;
      bool check_for_O_C_1_5_relation(mmdb::Atom *at_1, mmdb::Atom *at_2) const;  // check either way round


      int make_non_bonded_contact_restraints(const bonded_pair_container_t &bpc, const protein_geometry &geom);
      int make_non_bonded_contact_restraints(const bonded_pair_container_t &bpc,
					     const reduced_angle_info_container_t &ai,
					     const protein_geometry &geom);
      
      //! Set a flag that we have an OXT and we need to position it
      //after the refinement.
      void mark_OXT(const protein_geometry &geom);
      void position_OXT(); // called from update atoms if we have an OXT.
      bool have_oxt_flag;
      int oxt_index;
      std::vector<mmdb::Residue *> residues_with_OXTs;
      // these are in the order N, CA, C, O
      std::vector<clipper::Coord_orth> oxt_reference_atom_pos;
      // short int is_nucleotide(mmdb::Residue *res_p);
      bool do_numerical_gradients_flag;

      // validation:
      geometry_distortion_info_container_t
      distortion_vector(const gsl_vector *v) const;

      std::vector<bool>  make_fixed_flags(int index1, int index2) const;
      std::vector<bool>  make_non_bonded_fixed_flags(int index1, int index2) const;
      std::vector<bool>  make_fixed_flags(int index1, int index2, int index3) const;
      std::vector<bool>  make_fixed_flags(int index1, int index2, int index3, int index4) const;
      std::vector<bool>  make_fixed_flags(const std::vector<int> &indices) const;
      bool fixed_check(int i) const; // simple version of above functions.

      // return in milliseconds
      //
      // a coot-util function?
      // double time_diff(const timeval &current, const timeval &start) const;


      // a single H on this chiral centre is on the wrong side of the
      // chiral centre?
      // 
      // (not sure that this is the best name for this function).
      // 
      // Called from the minimize() function
      // 
      bool chiral_hydrogen_needs_pushing(const simple_restraint &chiral_restraint, const gsl_vector *v) const;
      bool check_pushable_chiral_hydrogens(gsl_vector *v); // and push them if needed (non-const *v)
      bool check_through_ring_bonds(gsl_vector *v); // and shorten them if needed (non-const *v)
      void push_chiral_hydrogen(const simple_restraint &chiral_restraint, gsl_vector *v);
      int get_chiral_hydrogen_index(int indexc, int index1, int index2, int index3) const;
      bool has_inverted_chiral_centre(const simple_restraint &chiral_restraint,
				      const gsl_vector *v) const;
      bool has_tiny_chiral_centre_volume(const simple_restraint &chiral_restraint,
					 const gsl_vector *v) const;
      bool bond_is_very_long(const simple_restraint &bond_restraint,
			     const gsl_vector *v) const;

      bool is_hydrogen(mmdb::Atom *at_p) const {
	 std::string ele = at_p->element;
	 if ((ele == "H") || (ele == " H"))
	    return 1;
	 else
	    return 0;
      }

      // return "" on no type found
      std::string get_type_energy(mmdb::Atom *at, const protein_geometry &geom) const {
	 std::string r;
	 if (at) { 
	    std::string atom_name = at->name;
	    const char *rn = at->GetResName();
	    if (rn) {
	       std::string residue_name = rn;
	       r = geom.get_type_energy(atom_name, residue_name);
	    } 
	 }
	 return r;
      }

      bonded_pair_container_t bonded_pairs_container;
      

   public: 

      enum link_torsion_restraints_type { NO_LINK_TORSION = 0, 
					  LINK_TORSION_RAMACHANDRAN_GOODNESS = 1,
					  LINK_TORSION_ALPHA_HELIX = 2,
					  LINK_TORSION_BETA_STRAND = 3 };

      // my_df_electron_density and electron_density_score need access
      // to fixed_atom_indices.
      std::vector<int> fixed_atom_indices; 

      // In all of these constructors the mmdb::PPAtom that is passed, either
      // explicitly or as part of an atom_selection_container_t has the
      // atoms to which is points *changed* by the minimization.
      // 
      // So, given that you may well want to keep your initial atom
      // positions, or reject the results of the idealization, you should
      // pass a copy of the atoms here, not the atoms themselves.
      //

      // This interface is withdrawn currently.  This is because
      // make_link_restraints does some residue selection and thus
      // needs a MMDBManager rather than a mmdb::PPAtom.
      // 
//       restraints_container_t(mmdb::PPAtom atoms_in, int n_at) {
// 	 verbose_geometry_reporting = 0;
// 	 atom = atoms_in; 
// 	 n_atoms = n_at; 
// 	 asc.mol = NULL; 
// 	 include_map_terms_flag = 0;
// 	 initial_position_params_vec.resize(3*n_at); 
// 	 for (int i=0; i<n_at; i++) {
// 	    initial_position_params_vec[3*i  ] = atoms_in[i]->x; 
// 	    initial_position_params_vec[3*i+1] = atoms_in[i]->y; 
// 	    initial_position_params_vec[3*i+2] = atoms_in[i]->z; 
// 	 }
//       }

      restraints_container_t(atom_selection_container_t asc_in) { 
	 verbose_geometry_reporting = NORMAL;
	 n_atoms = asc_in.n_selected_atoms; 
	 mol = asc_in.mol;
	 n_atoms = asc_in.n_selected_atoms;
	 atom = asc_in.atom_selection;
	 include_map_terms_flag = 0;
	 have_oxt_flag = 0;
	 do_numerical_gradients_flag = 0;
	 std::cout << "asc_in.n_selected_atoms " << asc_in.n_selected_atoms
		   << std::endl; 
	 initial_position_params_vec.resize(3*asc_in.n_selected_atoms);
	 lograma.init(LogRamachandran::All, 2.0, true);
	 from_residue_vector = 0;

	 for (int i=0; i<asc_in.n_selected_atoms; i++) {
	    initial_position_params_vec[3*i  ] = asc_in.atom_selection[i]->x; 
	    initial_position_params_vec[3*i+1] = asc_in.atom_selection[i]->y; 
	    initial_position_params_vec[3*i+2] = asc_in.atom_selection[i]->z; 
	 }
      }

      // for omega distortion info:
      restraints_container_t(atom_selection_container_t asc_in,
			     const std::string &chain_id);

      // iend_res is inclusive, so that 17,17 selects just residue 17.
      // 
      // Interface used by Regularize button callback:
      // 
      restraints_container_t(int istart_res, int iend_res,
			     short int have_flanking_residue_at_start,
			     short int have_flanking_residue_at_end,
			     short int have_disulfide_residues,
			     const std::string &altloc,
			     const char *chain_id,
			     mmdb::Manager *mol, // const in an ideal world
			     const std::vector<atom_spec_t> &fixed_atom_specs);

      // Interface used by Refine button callback:
      // 
      restraints_container_t(int istart_res, int iend_res,
			     short int have_flanking_residue_at_start,
			     short int have_flanking_residue_at_end,
			     short int have_disulfide_residues,
			     const std::string &altloc,
			     const char *chain_id,
			     mmdb::Manager *mol, // const in an ideal world
			     const std::vector<atom_spec_t> &fixed_atom_specs,
			     const clipper::Xmap<float> &map_in,
			     float map_weight);

      // For validation.
      //
      // The whole chain is selected (without flanking atoms) and we
      // use geometry_distortion() function.
      // 
      restraints_container_t(mmdb::PResidue *SelResidues, int nSelResidues,
			     const std::string &chain_id,
			     mmdb::Manager *mol);

      // 20081106 construct from a vector of residues, each of which
      // has a flag attached that denotes whether or not it is a fixed
      // residue (it would be set, for example in the case of flanking
      // reisdues).
      //
      // Consider also a regularize version of this (without map and
      // weight).
      // 
      // 20100210 that is what we have now, we use add_map() for the
      // restraints that fit to a map.
      //
      // if you need linkrs too, change links to a coot container
      // class for LINKs and LINKRs
      // const &link_container &links
      // 
      restraints_container_t(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
			     const std::vector<mmdb::Link> &links,
			     const protein_geometry &geom,			     
			     mmdb::Manager *mol,
			     const std::vector<atom_spec_t> &fixed_atom_specs);

      // 
      // geometric_distortions not const because we set restraints_usage_flag:
      //
      // return data useful for making the graphs:
      geometry_distortion_info_container_t
      geometric_distortions(restraint_usage_Flags flags);

      // Here we use the internal flags.  Causes crash currently (no inital atom positions?)
      geometry_distortion_info_container_t
      geometric_distortions() const;

      omega_distortion_info_container_t
      omega_trans_distortions(int mark_cis_peptides_as_bad_flag);
      

      // So, we provide easy(?) access to the atoms of next and
      // previous residues (to those in the atom selection
      // moving_residue_atoms).  It is also possible to select "fixed"
      // atoms in the graphics (so that they don't move).  Let's
      // provide a vector of indices in the the moving_residue_atoms
      // array to define those (lovely mixture of styles - heh).
      // 
      restraints_container_t(mmdb::PPAtom moving_residue_atoms, // e.g. atom of residue 16,17,18
			     int n_moving_residue_atoms, // e.g. 21
			     mmdb::PResidue previous_residue, // e.g. residue 15
			     mmdb::PResidue next_atom,
			     const std::vector<int> &fixed_atom_indices);

      restraints_container_t(){
	 from_residue_vector = 0;
	 include_map_terms_flag = 0;
      };

      ~restraints_container_t() {
	 if (from_residue_vector) {
	    if (atom) { 
	       // this is constructed manually.

	       // Oh we can't do this here because we copy the
	       // restraints in simple_refine_residues() and that
	       // shallow copies the atom pointer - the original
	       // restriants go out of scope and call this destructor.
	       //
	       // We need a new way to get rid of atom - c.f. the
	       // linear/conventional way?
	       
	       // delete [] atom;
	       // atom = NULL;
	    } 
	 } else {
	       // member data item mmdb::PPAtom atom is constructed by an
	       // mmdb SelectAtoms()/GetSelIndex() (which includes
	       // flanking atoms).
	    // 20081207: don't do this here now - because the
	    // memory/selection is deleted again in
	    // clear_up_moving_atoms(). It *should* be done here of
	    // course, but we'll save that for the future.
	    // 
	    //if (atom) { 
	    // mol->DeleteSelection(SelHnd_atom);
	    // atom = NULL;
	    // } 
	 } 
      }

      mmdb::Atom *get_atom(int i) const {
	 if (atom) 
	    return atom[i];
	 else
	    return NULL;
      }

      atom_spec_t get_atom_spec(int atom_index) const;

      void set_verbose_geometry_reporting() { 
	 verbose_geometry_reporting = VERBOSE;
      }

      void set_quiet_reporting() {
	 verbose_geometry_reporting = QUIET;
      } 

      void assign_fixed_atom_indices(const std::vector<atom_spec_t> &fixed_atom_specs);

      double initial_position(int i) {
	 return initial_position_params_vec[i]; 
      }

      int init_positions_size() {
	 return initial_position_params_vec.size(); 
      } 

      // Using this, we can mask out the restraints we don't want to
      // use using e.g. BONDS_MASK, BONDS_ANGLES_AND_PLANES etc...
      // 
      short int restraints_usage_flag;

      double starting_structure_diff_score(const gsl_vector *v, void *params); 


      short int include_map_terms() { 
	 return include_map_terms_flag; 
      }

      double Map_weight() const { 
	 return map_weight; 
      } 

      void setup_multimin_func() { 

	 multimin_func.f   = &distortion_score; 
	 multimin_func.df  = &my_df; 
	 multimin_func.fdf = &my_fdf; 
	 multimin_func.n = n_variables(); 
	 multimin_func.params = (double *) this; 
      } 

      int n_variables() { 
	 // return 3 * the number of atoms
	 return 3*n_atoms; 
      }

      // no longer done in the constructor:
      void add_map(const clipper::Xmap<float> &map_in, float map_weight_in) {
	 map = map_in;
	 map_weight = map_weight_in;
	 include_map_terms_flag = 1;
      }

      int size() const { return restraints_vec.size(); }

      // return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
      // 
      refinement_results_t minimize(restraint_usage_Flags);
      refinement_results_t minimize(restraint_usage_Flags, int nsteps, short int print_chi_sq_flag);
      void fix_chiral_atoms_maybe(gsl_vector *s);

      simple_restraint& operator[] (unsigned int i) {
	 return restraints_vec[i]; 
      } 

      // because chi_squareds is const:
      const simple_restraint& operator[] (const unsigned int &i) const { 
	 return restraints_vec[i]; 
      } 
  
      void setup_gsl_vector_variables();

      double electron_density_score_at_point(const clipper::Coord_orth &ao) const;
      clipper::Grad_orth<double> electron_density_gradient_at_point(const clipper::Coord_orth &ao) const; 
      
      // We need to fill restraints_vec (which is a vector of
      // simple_restraint) using the coordinates () and the dictionary of
      // restraints, protein_geometry geom.
      int make_restraints(const protein_geometry &geom,
			  restraint_usage_Flags flags,
			  bool do_residue_internal_torsions,
			  float rama_plot_target_weight,
			  bool do_rama_plot_retraints, 
			  pseudo_restraint_bond_type sec_struct_pseudo_bonds);

      unsigned int test_function(const protein_geometry &geom);
      unsigned int inline_const_test_function(const protein_geometry &geom) const {
	 std::cout << "----- inline_const_test_function() with geom of size : " << geom.size()
		   << std::endl;
	 std::cout << "    geom ref pointer " << &geom << std::endl;
	 return geom.size();
      } 
      unsigned int const_test_function(const protein_geometry &geom) const;

      mmdb::Atom *add_atom_pull_restraint(atom_spec_t spec, clipper::Coord_orth pos);
      void clear_atom_pull_restraint();
      
      void add_extra_restraints(const extra_restraints_t &extra_restraints,
				const protein_geometry &geom);
      // and that calls:
      void add_extra_bond_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_angle_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_torsion_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_start_pos_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_parallel_plane_restraints(const extra_restraints_t &extra_restraints,
					       const protein_geometry &geom);

      // old code:
      // Read restraints from the refmac .rst file
      int make_restraints(std::string ciffilename);

      void update_atoms(gsl_vector *s);
      // return the WritePDBASCII() status, or -1 if mol was 0.
      int write_new_atoms(std::string pdb_file_name);
      void info() const;

      // The restraints have flags in them to tell us whether to
      // include the gradients for this atom, however, those cannot be
      // used when we are considering the gradients from the map.  So
      // let's create another array of flags (atom indexed (not
      // parameter indexed)) that tells us whether these atoms of the
      // atom selection are OK for atom gradients (i.e. atoms of the
      // flanking residues will have zeroes here).
      //
      std::vector<bool> use_map_gradient_for_atom;
      std::vector<double> atom_z_weight;  // weight e.d. fit by atomic number

      // Make a MMDBManager from the selection and return this new
      // mmdb::PManager.  It is the users responsibility to delete it.
      mmdb::Manager *results() const;
      void adjust_variables(const atom_selection_container_t &asc);

      LogRamachandran LogRama() const { return lograma; };

      // here phi and psi are in clipper units (radians).
      double rama_prob(const double &phi_rads, const double &psi_rads) const {
	 return lograma.interp(phi_rads, psi_rads);
      }
      LogRamachandran::Lgrad rama_grad(const double &phir, const double &psir) const {
	 return lograma.interp_grad(phir, psir);
      } 


      // Allow public access to this - the general method for knowing if 2 residues have a (dictionary) link.
      // 
      // find disulphides, protein-glycan bonds etc.
      // 
      // Return the link type and a residue order switch flag.
      // Return link_type as "" if not found.
      // 
      // and with a flag, using CLinks and SSBond, rather than
      // guessing (carbohydrates and disulfides):
      // 
      std::pair<std::string, bool> find_link_type_complicado(mmdb::Residue *first,
							     mmdb::Residue *second,
							     const protein_geometry &geom) const;

      // Allow public access to this - we need it to find the links
      // between residues when all we have to go on is the refmac
      // dictionary - no LINKRs and no user input.
      bonded_pair_container_t bonded_residues_from_res_vec(const protein_geometry &geom) const;

      // Using bonded pairs internal copy, modify residues as needed
      // by deleting atoms in chem mods.
      //
      // But what about the mod_OXT code?  How does that fit in here?
      //
      void apply_link_chem_mods(const protein_geometry &geom);
      
      double geman_mcclure_alpha; // = 0.02 or something set in init_shared_pre(). // needed for derivative calculation
                                  // (which is not done in this class)

      // return true when turned off
      bool turn_off_when_close_target_position_restraint();

      // more debugging interface:
      //
      void set_do_numerical_gradients() { do_numerical_gradients_flag = 1;}
      bool do_numerical_gradients_status() { return do_numerical_gradients_flag; }

   }; 

   void simple_refine(mmdb::Residue *residue_p,
		      mmdb::Manager *mol,
		      const dictionary_residue_restraints_t &dict_restraints);

} // namespace coot


#endif // HAVE_GSL
#endif // HAVE_SIMPLE_RESTRAINT_HH
