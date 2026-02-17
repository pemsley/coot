// -*-c++-*-
/* ideal/simple-resetraint.hh
 *
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
 * Copyright 2008 by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifndef HAVE_SIMPLE_RESTRAINT_HH
#define HAVE_SIMPLE_RESTRAINT_HH

#include <vector>
#include <list>
#include <string>
#include <stdexcept>
#include <memory>

#ifndef __NVCC__
#include <thread>
#include <atomic>
#endif // __NVCC__

#ifndef __NVCC__
#ifdef HAVE_BOOST

#ifdef HAVE_BOOST_THREAD // part of DEFS in Makefile
#define HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#endif


/*
  if I comment out the #include "compat/coot-getopt.h"
  then this compilation error goes away:

In file included from /Users/pemsley/autobuild/build-refinement-pre-release-gtk2-python/include/boost/config/posix_features.hpp:18:
/usr/include/unistd.h:503:6: error: conflicting types for 'getopt'
int      getopt(int, char * const [], const char *) __DARWIN_ALIAS(getopt);
         ^
../../coot/compat/coot-getopt.h:153:12: note: previous declaration is here
extern int getopt ();
           ^
*/

// #include "compat/coot-getopt.h"
#include "utils/ctpl.h"
#endif // HAVE_BOOST
#endif // __NVCC__

#include "refinement-results-t.hh"

#include <mmdb2/mmdb_manager.h>
#include "coot-utils/bonded-pairs.hh"

#include "parallel-planes.hh"

#include "extra-restraints.hh"

#include "model-bond-deltas.hh"

#include "compat/coot-sysdep.h"

// refinement_results_t is outside of the GSL test because it is
// needed to make the accept_reject_dialog, and that can be compiled
// without the GSL.
//
namespace coot {

   enum { UNSET_INDEX = -1 };

   class rama_triple_t {
   public:
      mmdb::Residue *r_1;
      mmdb::Residue *r_2;
      mmdb::Residue *r_3;
      std::string link_type;
      bool fixed_1;
      bool fixed_2;
      bool fixed_3;
      rama_triple_t(mmdb::Residue *r1, mmdb::Residue *r2, mmdb::Residue *r3,
                    const std::string &link_type_in) : link_type(link_type_in) {
         r_1 = r1;
         r_2 = r2;
         r_3 = r3;
         fixed_1 = 0;
         fixed_2 = 0;
         fixed_3 = 0;
      }
      rama_triple_t(mmdb::Residue *r1, mmdb::Residue *r2, mmdb::Residue *r3,
                    const std::string &link_type_in,
                    bool fixed_1_in, bool fixed_2_in, bool fixed_3_in) : link_type(link_type_in) {
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
      bool zero_gradients;
      double theta; // the torsion angle
      double tan_theta; // store the tan too
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

}

// we don't want to compile anything if we don't have gsl

#include <map>

#include "gsl/gsl_multimin.h"

// #include "Cartesian.h"
#include "coot-utils/atom-selection-container.hh" // for atom_selection_container_t, this and
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

#include "new-linked-residue-t.hh"


// For Kevin's (Log) Ramachandran Plot and derivatives
#include "lograma.h"
// For ZO's Ramachandran Plot and derivatives
#include "zo-rama.hh"

#ifndef DEGTORAD
#define DEGTORAD 0.017453293
#endif
#ifndef RADTODEG
#define RADTODEG 57.29577793
#endif

// for conversion between the old atom index and new atom index for extra bond restraints
#define ATOM_INDEX_MAX 500000 // is that enough?

namespace coot {


   // restraint types:
   //
   enum restraint_type_t {BOND_RESTRAINT=1, ANGLE_RESTRAINT=2, TORSION_RESTRAINT=4, PLANE_RESTRAINT=8,
                          NON_BONDED_CONTACT_RESTRAINT=16, CHIRAL_VOLUME_RESTRAINT=32, RAMACHANDRAN_RESTRAINT=64,
                          START_POS_RESTRAINT=128,
                          TARGET_POS_RESTRAINT=256, // restraint to make an atom be at a position
                          PARALLEL_PLANES_RESTRAINT=512,
                          GEMAN_MCCLURE_DISTANCE_RESTRAINT=1024,
                          TRANS_PEPTIDE_RESTRAINT=2048,
                          IMPROPER_DIHEDRAL_RESTRAINT=4096
   };

   enum pseudo_restraint_bond_type {NO_PSEUDO_BONDS, HELIX_PSEUDO_BONDS,
                                    STRAND_PSEUDO_BONDS};


   // restraints_usage_flag
   //
   // use with restraints_usage_flag & 1 for BONDS or
   //          restraints_usage_flag & 2 for ANGLES etc.
   //
   enum restraint_usage_Flags { NO_GEOMETRY_RESTRAINTS = 0,
                                BONDS = 1,
                                ANGLES = 2,
                                BONDS_AND_ANGLES = 3,
                                TORSIONS = 4,
                                NON_BONDED = 16,
                                CHIRAL_VOLUMES = 32,
                                //PLANES = 8,
                                RAMA = 64,
                                GEMAN_MCCLURE_DISTANCE_RESTRAINTS=1024,
                                TRANS_PEPTIDE_RESTRAINTS=2048,
                                BONDS_ANGLES_AND_TORSIONS = 7,
                                BONDS_ANGLES_TORSIONS_AND_PLANES = 15,
                                BONDS_AND_PLANES = 9,
                                BONDS_ANGLES_AND_PLANES = 11, // no torsions
                                BONDS_ANGLES_AND_CHIRALS = 35, // no torsions
                                BONDS_AND_NON_BONDED = 17,
                                BONDS_ANGLES_AND_NON_BONDED = 19,
                                BONDS_ANGLES_CHIRALS_AND_NON_BONDED = 19 + 32, // pre-sanitize
                                BONDS_ANGLES_TORSIONS_AND_NON_BONDED = 23,
                                BONDS_ANGLES_PLANES_AND_NON_BONDED = 27,
                                BONDS_ANGLES_PLANES_NON_BONDED_AND_TRANS_PEPTIDE_RESTRAINTS = 27 + 2048,
                                BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS = 59,
                                BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_TRANS_PEPTIDE_RESTRAINTS = 59+2048,
                                BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED = 31,
                                BONDS_ANGLES_TORSIONS_PLANES_AND_CHIRALS = 47,
                                BONDS_ANGLES_PLANES_AND_CHIRALS = 43,
                                BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS = 63,
                                JUST_RAMAS = 64,
                                BONDS_ANGLES_TORSIONS_NON_BONDED_AND_CHIRALS = 55,
                                BONDS_ANGLES_TORSIONS_NON_BONDED_CHIRALS_AND_PLANES = 55 + 8,
                                BONDS_ANGLES_TORSIONS_NON_BONDED_CHIRALS_AND_TRANS_PEPTIDE_RESTRAINTS = 55+2048,

                                BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA = 127,

                                BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_PARALLEL_PLANES = 187,
                                BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_PARALLEL_PLANES = 191,
                                BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_RAMA_AND_PARALLEL_PLANES = 255,

                                // These become pushed along one slot to fit in target pos restraints
                                // GEMAN_MCCLURE_DISTANCE_RESTRAINTS = 512,
                                // BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_GEMAN_MCCLURE_DISTANCES = 63+512,
                                // TYPICAL_RESTRAINTS               = 1+2+  8+16+32+128+256+512,
                                // typical restraints add trans-peptide restraints
                                // TYPICAL_RESTRAINTS               = 1+2+  8+16+32+128+256+512+1024,
                                // TYPICAL_RESTRAINTS_WITH_TORSIONS = 1+2+4+8+16+32+128+256+512+1024,
                                // TYPICAL_NO_PLANES = 1+2+4 +16+32+128+256+512+1024,
                                // ALL_RESTRAINTS = 1+2+4+8+16+32+64+128+256+512+1024

                                BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_GEMAN_MCCLURE_DISTANCES = 63+1024,
                                // typical restraints add trans peptide restraints and geman-mcclure
                                TYPICAL_RESTRAINTS                = 1+2+  8+16+32+128+256+512+1024+2048,
                                TYPICAL_RESTRAINTS_NO_PLANES      = 1+2+    16+32+128+256+512+1024+2048,
                                TYPICAL_RESTRAINTS_WITH_IMPROPERS = 1+2+  8+16+32+128+256+512+1024+2048+4096,  // for testing planes -> improper-dihedrals
                                TYPICAL_RESTRAINTS_WITH_TORSIONS = 1+2+4+8+16+32+128+256+512+1024+2048,
                                ALL_RESTRAINTS = 1+2+4+8+16+32+64+128+256+512+1024+2048+4096 // adds torsions and ramas

   };

   enum peptide_restraints_usage_Flags { OMEGA_TORSION = 1,
                                         OMEGA_PHI_PSI_TORSION = 3 };

   enum { BONDS_MASK = 1,  ANGLES_MASK = 2, TORSIONS_MASK = 4, PLANES_MASK = 8,
          NON_BONDED_MASK = 16,
          CHIRAL_VOLUME_MASK = 32,
          RAMA_PLOT_MASK = 64,
          START_POS_RESTRAINT_MASK = 128,
          PARALLEL_PLANES_MASK = 256,
          GEMAN_MCCLURE_DISTANCE_MASK = 1024,
          TRANS_PEPTIDE_MASK = 2048,
          IMPROPER_DIHEDRALS_MASK = 4096
   };


   class ramachandran_restraint_flanking_residues_helper_t {
   public:
      int resno_first;
      int resno_third;
      std::vector<bool> is_fixed;
      ramachandran_restraint_flanking_residues_helper_t() {
         is_fixed.resize(3, false);
         resno_first = -1;
         resno_third = -1;
      }
   };

   bool residue_sorter(const std::pair<bool, mmdb::Residue *> &r1,
                       const std::pair<bool, mmdb::Residue *> &r2);

   // ---------------------------------------------------------------
   // ---------------------------------------------------------------
   //     class simple_restraint
   // ---------------------------------------------------------------
   // ---------------------------------------------------------------

   class simple_restraint {

      void init() {
         restraint_type = BOND_RESTRAINT;
         restraint_index = -1;
         restraints_index = -1; // same thing?
         atom_index_1 = -1;
         atom_index_2 = -1;
         observed_value = 0.0;
         sigma = 0.0;
         target_value = 0.0;
         is_user_defined_restraint = false;
         is_H_non_bonded_contact = false;
         is_hydrogen_bond = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         nbc_function = HARMONIC; // not used
         torsion_restraint_weight = 1.0;
         target_chiral_volume = 0.0;
         chiral_volume_sign = 1;
         chiral_hydrogen_index = -1;
         is_closed = false;
         n_atoms_from_all_restraints = 0;
         periodicity = 0;
         atom_index_centre = -1;
         atom_index_3 = -1;
         atom_index_4 = -1;
         atom_index_5 = -1;
         atom_index_6 = -1;
      }

   public:

      int restraint_index;
      int atom_index_1, atom_index_2, atom_index_3, atom_index_4, atom_index_5, atom_index_6;
      int atom_index_centre;
      bool is_closed; // so that we can "remove" restraints without deleting them
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
                                 // centre, then the atom index,
                                 // otherwise this is -1.
      atom_spec_t atom_spec; // for pull atoms (so that we can on the fly delete this restraints)
      int n_atoms_from_all_restraints; // for debugging GSL/atom index errors
      int restraints_index;            // ditto. Is this the same thing as restraints_index?
      std::vector<bool> fixed_atom_flags;
      std::vector<bool> fixed_atom_flags_other_plane;
      bool is_user_defined_restraint;
      bool is_H_non_bonded_contact;
      bool is_hydrogen_bond;
      bool is_single_Hydrogen_atom_angle_restraint;
      double torsion_restraint_weight;
      //
      // for mouse pull on an atom: this is where the user wants the atom to be
      //
      clipper::Coord_orth atom_pull_target_pos;

      std::string rama_plot_residue_type; // so that we look up the correct residue type
                                          // for this (middle-of-three) residue

      enum nbc_function_t { LENNARD_JONES, HARMONIC};
      nbc_function_t nbc_function;

      // allocator for geometry_distortion_info_t
      simple_restraint() {
         init();
      }

      // Bond
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2,
                       const std::vector<bool> &fixed_atom_flags_in,
                       float tar,
                       float sig, float obs) : fixed_atom_flags(fixed_atom_flags_in) {

         init();
         restraint_type = rest_type;
         restraint_index = -1;
         atom_index_1 = atom_1;
         atom_index_2 = atom_2;
         observed_value = obs;
         sigma = sig;
         target_value = tar;
         is_user_defined_restraint = false;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         nbc_function = HARMONIC; // not used
         torsion_restraint_weight = 1.0;

         // This finds a coding error
         if (rest_type != BOND_RESTRAINT) {
            std::cout << "BOND ERROR in simple_restraint()" << std::endl;
         }
      };

      // Geman-McClure distance (no obs)
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2,
                       const std::vector<bool> &fixed_atom_flags_in,
                       float tar, float sig) : fixed_atom_flags(fixed_atom_flags_in) {

         init();
         restraint_type = rest_type;
         restraint_index = -1;
         atom_index_1 = atom_1;
         atom_index_2 = atom_2;
         observed_value = -1;  // not given or used
         sigma = sig;
         target_value = tar;
         is_user_defined_restraint = true;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         torsion_restraint_weight = 1.0;

         if (rest_type != restraint_type_t(GEMAN_MCCLURE_DISTANCE_MASK)) {
            std::cout << "BOND ERROR (Geman McClure) in simple_restraint()"
                      << std::endl;
         }
      };

      // Angle
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2,
                       int atom_3,
                       const std::vector<bool> &fixed_atom_flags_in,
                       float tar,
                       float sig, bool is_single_Hydrogen_atom_angle_restraint_in) : fixed_atom_flags(fixed_atom_flags_in) {

         init();
         restraint_type = rest_type;
         restraint_index = -1;
         atom_index_1 = atom_1;
         atom_index_2 = atom_2;
         atom_index_3 = atom_3;
         sigma = sig;
         target_value = tar;
         is_user_defined_restraint = 0;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = is_single_Hydrogen_atom_angle_restraint_in;
         torsion_restraint_weight = 1.0;
         if (rest_type != ANGLE_RESTRAINT) {
            std::cout << "ERROR::::: PROGRAM ERROR - ANGLE ERROR" << std::endl;
         }
      };

      // Torsion
      simple_restraint(restraint_type_t rest_type, int atom_1, int atom_2, int atom_3, int atom_4,
                       const std::vector<bool> &fixed_atom_flags_in,
                       float tar,
                       float sig, float weight, int periodicity_in) : fixed_atom_flags(fixed_atom_flags_in) {

         init();
         restraint_type = rest_type;
         restraint_index = -1;
         atom_index_1 = atom_1;
         atom_index_2 = atom_2;
         atom_index_3 = atom_3;
         atom_index_4 = atom_4;
         observed_value = 0.0;
         torsion_restraint_weight = weight;
         sigma = sig;
         target_value = tar;
         periodicity = periodicity_in;
         is_user_defined_restraint = 0;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         if ((rest_type != TORSION_RESTRAINT) && (rest_type != TRANS_PEPTIDE_RESTRAINT)) {
            std::cout << "ERROR::::: PROGRAM ERROR - TORSION/TRANSP-PEP ERROR" << std::endl;
         }
      }

      // Rama
      simple_restraint(restraint_type_t rest_type,
                       const std::string &rama_plot_zo_residue_type,
                       int atom_1, int atom_2, int atom_3, int atom_4, int atom_5,
                       const std::vector<bool> &fixed_atom_flags_in) : fixed_atom_flags(fixed_atom_flags_in)  {

         init();
         restraint_type = rest_type;
         restraint_index = -1;
         rama_plot_residue_type = rama_plot_zo_residue_type;
         atom_index_1 = atom_1;
         atom_index_2 = atom_2;
         atom_index_3 = atom_3;
         atom_index_4 = atom_4;
         atom_index_5 = atom_5;
         is_user_defined_restraint = 0;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         torsion_restraint_weight = 1.0;
         if (rest_type != RAMACHANDRAN_RESTRAINT) {
            std::cout << "ERROR:: RAMACHANDRAN_RESTRAINT ERROR" << std::endl;
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
         restraint_index = -1;
         restraint_type = restraint_type_in;

         plane_atom_index.resize(atom_index_in.size());
         for (unsigned int i=0; i<atom_index_in.size(); i++)
            plane_atom_index[i] = std::pair<int, double> (atom_index_in[i], sig);

         target_value = 0.0; // not needed for planes
         sigma = sig;
         fixed_atom_flags = fixed_atom_flags_in;
         is_user_defined_restraint = 0;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         torsion_restraint_weight = 1.0;
      }

      // modern (atoms individually weighted) Plane
      //
      simple_restraint(restraint_type_t restraint_type_in,
                       const std::vector<std::pair<int, double> > &atom_index_sigma_in,
                       const std::vector<bool> &fixed_atom_flags_in) :
         plane_atom_index(atom_index_sigma_in),
         fixed_atom_flags(fixed_atom_flags_in) {

         init();
         //
         // Check restraint_type?
         // Currently the only thing with an std::vector <int>
         // atom_index_in is a plane restraint.  This could well
         // change in the future.
         //
         restraint_type = restraint_type_in;
         restraint_index = -1;

         target_value = 0.0; // not needed for planes
         sigma = 0.02; // hack
         is_user_defined_restraint = 0;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         torsion_restraint_weight = 1.0;
      }


      // Parallel planes (actually angle-between-planes,typically-zero)
      simple_restraint(restraint_type_t restraint_type_in,
                       const std::vector<int> &atom_index_plane_1_in,
                       const std::vector<int> &atom_index_plane_2_in,
                       const std::vector<bool> &fixed_atom_flags_plane_1_in,
                       const std::vector<bool> &fixed_atom_flags_plane_2_in,
                       double target_angle_in,
                       double sigma_in) {
         init();
         restraint_type = restraint_type_in;
         restraint_index = -1;
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
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         torsion_restraint_weight = 1.0;
      }

     // improper dihedral (new-style 4 atom plane restraints)
     simple_restraint(restraint_type_t restraint_type_in,
                      int index_1, int index_2, int index_3, int index_4,
                      float sigma_in, const std::vector<bool> &fixed_atom_flags_in) : fixed_atom_flags(fixed_atom_flags_in) {

        init();
        restraint_type = restraint_type_in;
        restraint_index = -1;

        atom_index_1 = index_1;
        atom_index_2 = index_2;
        atom_index_3 = index_3;
        atom_index_4 = index_4;

        sigma = sigma_in;
        is_user_defined_restraint = false;
        is_H_non_bonded_contact = false;
        is_single_Hydrogen_atom_angle_restraint = false;
        torsion_restraint_weight = 1.0;
     }

      // Non-bonded
      //
      //- are you sure that this is the constructor that you want?
      //
      simple_restraint(restraint_type_t restraint_type_in,
                       int index_1,
                       int index_2,
                       const std::string &atom_1_type,
                       const std::string &atom_2_type,
                       const std::vector<bool> &fixed_atom_flags_in,
                       const protein_geometry &geom) {

         init();
         restraint_index = -1;
         torsion_restraint_weight = 1.0;
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
            sigma = 0.02; // probably not this constructor....
            fixed_atom_flags = fixed_atom_flags_in;
            is_user_defined_restraint = 0;
            is_H_non_bonded_contact = false;
            is_single_Hydrogen_atom_angle_restraint = false;
         } else {
            std::cout << "ERROR:: bad simple_restraint constructor usage "
                      << "- should be non-bonded\n";
         }
      }

      // Non-bonded v2
      simple_restraint(restraint_type_t restraint_type_in,
                       const nbc_function_t &nbc_func_type,
                       int index_1,
                       int index_2,
                       bool is_H_non_bonded_contact_in,
                       const std::vector<bool> &fixed_atom_flags_in,
                       double dist_min) {

         init();
         restraint_index = -1;
         torsion_restraint_weight = 1.0;
         if (restraint_type_in == NON_BONDED_CONTACT_RESTRAINT) {
            restraint_type = restraint_type_in;
            atom_index_1 = index_1;
            atom_index_2 = index_2;
            target_value = dist_min;
            // nbc_function = HARMONIC;
            nbc_function = nbc_func_type;
            sigma = 0.06;
            fixed_atom_flags = fixed_atom_flags_in;
            is_user_defined_restraint = 0;
            is_single_Hydrogen_atom_angle_restraint = false;
            is_H_non_bonded_contact = is_H_non_bonded_contact_in;
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

         init();
         restraint_index = -1;
         torsion_restraint_weight = 1.0;
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
            is_H_non_bonded_contact = false;
            is_single_Hydrogen_atom_angle_restraint = false;
         }
      }

      // start pos
      simple_restraint(restraint_type_t rest_type, int atom_1,
                       bool fixed_atom_flag_in,
                       float sig, float obs){

         init();
         torsion_restraint_weight = 1.0;
         restraint_index = -1;
         restraint_type = rest_type;
         atom_index_1 = atom_1;
         observed_value = obs;
         sigma = sig;
         fixed_atom_flags = std::vector<bool> (1,fixed_atom_flag_in);
         is_user_defined_restraint = 0;
         is_H_non_bonded_contact = false;
         is_single_Hydrogen_atom_angle_restraint = false;
         is_closed = false;

         if (rest_type != START_POS_RESTRAINT) {
            std::cout << "ERROR:: START POS ERROR" << std::endl;
         }
      }

      // target_position, including pull_atoms
      simple_restraint(restraint_type_t rest_type, int atom_idx,
                       const atom_spec_t &spec_in,
                       const clipper::Coord_orth &pos) :
         atom_spec(spec_in),
         atom_pull_target_pos(pos) {

         init();
         restraint_index = -1;
         restraint_type = rest_type;
         atom_index_1 = atom_idx;
         is_closed = false;
         torsion_restraint_weight = 1.0;
         if (rest_type != TARGET_POS_RESTRAINT) {
            std::cout << "ERROR:: TARGET POS ERROR" << std::endl;
         }
      }

      void close() { is_closed = true; }

      std::pair<bool, double> get_nbc_dist(const std::string &atom_1_type,
                                           const std::string &atom_2_type,
                                           const protein_geometry &geom);

      double torsion_distortion(double model_torsion) const;
      // distortion (penalty score) and bond length delta from target value
      std::pair<double, double> distortion(mmdb::PAtom *atoms, const double &lj_epsilon) const; // i.e. atom
      std::string type() const; // a string representation of the restraint type
      friend std::ostream &operator<<(std::ostream &s, const simple_restraint &r);
      std::string format(mmdb::PAtom *atoms_vec, double distortion) const;
      void set_torsion_restraint_weight(const double &tw) { torsion_restraint_weight = tw; }
   };

   // ------------------------------ end of simple_restraint ----------------------------------------


   std::ostream &operator<<(std::ostream &s, const simple_restraint &r);
   bool target_position_eraser(const simple_restraint &r); // this is static, I guess

   // a good example for erase... remove_if (another is the crankshaft eraser)
   class target_position_for_atom_eraser {
   public:
      explicit target_position_for_atom_eraser(const atom_spec_t &spec_in) : spec(spec_in) {}
      atom_spec_t spec;
      bool operator() (const simple_restraint &r) const {
#ifdef SWIG
         // 20221028-PE it feels like a trap
#else
         if (r.restraint_type == restraint_type_t(TARGET_POS_RESTRAINT)) {
            if (r.atom_spec == spec) {
               return true;
            }
         }
#endif
         return false;
      }
   };

   class turn_off_when_close_target_position_restraint_eraser {
      int n_atoms;
      mmdb::PAtom *atoms;
      double close_dist;
      atom_spec_t exclude_spec;
   public:
      turn_off_when_close_target_position_restraint_eraser(double close_dist_in, mmdb::PAtom *atoms_in, int n_atoms_in,
                                                           const atom_spec_t &exclude_spec_in) : exclude_spec(exclude_spec_in) {
         atoms = atoms_in;
         n_atoms = n_atoms_in;
         close_dist = close_dist_in; // 0.6; // was 0.5; // was 0.4
      }
      bool operator() (const simple_restraint &r) const {
         bool v = false;
         if (r.restraint_type == restraint_type_t(TARGET_POS_RESTRAINT)) {
            clipper::Coord_orth p_1 = co(atoms[r.atom_index_1]);
            double d = sqrt((p_1-r.atom_pull_target_pos).lengthsq());
            if (d < close_dist) {
               atom_spec_t t(atoms[r.atom_index_1]);
#ifdef SWIG
               // 20221028-PE more trap?
#else
               if (t != exclude_spec)
                  v = true;
#endif
            }
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
                                 const residue_spec_t &residue_spec_in) :
         restraint(rest_in), residue_spec(residue_spec_in) {
         distortion_score = distortion_in;
         is_set = true;
      }
      geometry_distortion_info_t() {
         is_set = false;
         distortion_score = 0.0;
      }
      bool is_set;
      double distortion_score;
      simple_restraint restraint;
      std::vector<int> atom_indices;
      std::vector<atom_spec_t> atom_specs;
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
      bool initialised_p() const { return is_set; }
   };
   std::ostream &operator<<(std::ostream &s, geometry_distortion_info_t);

   class geometry_distortion_info_pod_t {
   public:
      geometry_distortion_info_pod_t(double distortion_in,
                                     const simple_restraint &rest_in,
                                     const residue_spec_t &residue_spec_in) :
         restraint(rest_in), residue_spec(residue_spec_in) {
         distortion_score = distortion_in;
         is_set = true;
      }
      geometry_distortion_info_pod_t() {
         is_set = false;
         distortion_score = 0.0;
      }
      bool initialised_p;
      simple_restraint restraint;
      bool is_set;
      std::vector<atom_spec_t> atom_specs;
      residue_spec_t residue_spec;
      double distortion_score;
      double get_distortion() const { return distortion_score; }
   };

   class geometry_distortion_info_container_t {
   public:
      std::vector<geometry_distortion_info_t> geometry_distortion;
      std::string chain_id;
      mmdb::PAtom *atom;
      int n_atoms;
      int min_resno;
      int max_resno;
      geometry_distortion_info_container_t() {
         atom = 0;
         n_atoms = 0;
         max_resno = 0;
         min_resno = 0;
      }
      geometry_distortion_info_container_t(const std::vector<geometry_distortion_info_t> &geometry_distortion_in,
                                           mmdb::PAtom *atom_in,
                                           const std::string &chain_id_in) :
         geometry_distortion(geometry_distortion_in), chain_id(chain_id_in) {
         atom = atom_in;
         n_atoms = 0; // this is worrying - why is this not passed. Who uses this constructor?
         max_resno = 0;
         min_resno = 0;
      }
      geometry_distortion_info_container_t(mmdb::PAtom *atom_in, int n_atoms_in,
                                           const std::string &chain_id_in) : chain_id(chain_id_in) {
         atom = atom_in;
         n_atoms = n_atoms_in;
         max_resno = 0;
         min_resno = 0;
      }
      void set_min_max(int min_resno_in, int max_resno_in) {
         min_resno = min_resno_in;
         max_resno = max_resno_in;
      }
      unsigned int size () const { return geometry_distortion.size(); }
      double print() const;  // return the total distortion
      double print_using_atom_specs() const;  // safe version using stored atom_specs
      double distortion() const;  // return the total distortion
      double distortion_sum() const; // return the sum of the distortions from the restraints - no calculation
      geometry_distortion_info_t get_geometry_distortion_info(unsigned int idx) const;
      friend std::ostream &operator<<(std::ostream &s, geometry_distortion_info_container_t);
   };
   std::ostream &operator<<(std::ostream &s, geometry_distortion_info_container_t gdic);

   class geometry_distortion_info_pod_container_t {
   public:
      std::vector<geometry_distortion_info_pod_t> geometry_distortion;
      std::string chain_id;
      int n_atoms;
      int min_resno;
      int max_resno;
      geometry_distortion_info_pod_container_t() {
         n_atoms = 0;
         max_resno = 0;
         min_resno = 0;
      }
      geometry_distortion_info_pod_container_t(const std::vector<geometry_distortion_info_pod_t> &geometry_distortion_in,
                                               const std::string &chain_id_in) :
         geometry_distortion(geometry_distortion_in), chain_id(chain_id_in) {
         n_atoms = 0; // this is worrying - why is this not passed. Who uses this constructor?
         max_resno = 0;
         min_resno = 0;
      }
      double print() const;  // return the total distortion
      unsigned int size () const { return geometry_distortion.size(); }
      geometry_distortion_info_pod_t get_geometry_distortion_info(unsigned int idx) const;
      double get_distortion() const;
   };

   class omega_distortion_info_t {
   public:
      int resno;
      double distortion;
      std::string info_string;
      omega_distortion_info_t(int resno_in, double distortion_in, const std::string &s) : info_string(s) {
         resno = resno_in;
         distortion = distortion_in;
      }
   };

   class omega_distortion_info_container_t {
   public:
      std::string chain_id;
      std::vector<omega_distortion_info_t> omega_distortions; // in degrees away from 180
      int min_resno;
      int max_resno;
      omega_distortion_info_container_t(const std::string &chain_id_in, int min_resno_in, int max_resno_in) : chain_id(chain_id_in) {
         min_resno = min_resno_in;
         max_resno = max_resno_in;
      }
   };

   // params can't be const because distortion_score function is an argument to
   // GSL multimin and that takes a function that doesn't have const void *params.
   double distortion_score(const gsl_vector *v, void *params);
#ifdef HAVE_CXX_THREAD
#ifndef __NVCC__
   // return value in distortion
   /* restraints_indices version - seems to slow things down!?
   void distortion_score_multithread(int thread_id,
                                     const gsl_vector *v, void *params,
                                     const std::vector<std::size_t> &restraint_indices,
                                     double *distortion,
                                     std::atomic<unsigned int> &done_count); */
   void distortion_score_multithread(int thread_id,
                                     const gsl_vector *v, void *params,
                                     int start, int end,
                                     double *distortion,
                                     std::atomic<unsigned int> &done_count);
#endif // __NVCC__
#endif // HAVE_CXX_THREAD
   void distortion_score_single_thread(const gsl_vector *v, void *params,
                                       int idx_start, int idx_end, double *distortion);
   double distortion_score_bond(const simple_restraint &bond_restraint,
                                const gsl_vector *v);
   double distortion_score_geman_mcclure_distance(const simple_restraint &bond_restraint,
                                                  const gsl_vector *v,
                                                  const double &alpha);
   double distortion_score_angle(const simple_restraint &angle_restraint,
                                 const gsl_vector *v);
   // torsion score can throw a std::runtime_error if there is a problem calculating the torsion.

   double distortion_score_torsion(unsigned int idx_restraint,
                                   const simple_restraint &torsion_restraint,
                                   const gsl_vector *v);
   double distortion_score_torsion_fourier_series(const simple_restraint &torsion_restraint,
                                                  const gsl_vector *v);
   double distortion_score_trans_peptide(const int &restraint_index,
                                         const simple_restraint &trans_peptide_restraint,
                                         const gsl_vector *v);
   double distortion_score_plane(const simple_restraint &plane_restraint,
                                  const gsl_vector *v);
   double distortion_score_improper_dihedral(const simple_restraint &id_restraint,
                                  const gsl_vector *v);
   double distortion_score_chiral_volume(const simple_restraint &chiral_restraint,
                                         const gsl_vector *v);
   double distortion_score_rama(const simple_restraint &rama_restraint,
                                const gsl_vector *v,
                                const LogRamachandran &lograma,
                                double rama_plot_weight);
   double distortion_score_rama(const simple_restraint &rama_restraint,
                                const gsl_vector *v,
                                const zo::rama_table_set &rama,
                                float rama_plot_weight);
   double distortion_score_start_pos(const simple_restraint &start_pos_restraint,
                            void *params,
                            const gsl_vector *v);
   double distortion_score_target_pos(const simple_restraint &start_pos_restraint,
                                      double scale_factor,
                                      const gsl_vector *v);
   double distortion_score_non_bonded_contact(const simple_restraint &plane_restraint,
                                              const double &lennard_jones_epsilon,
                                              const gsl_vector *v);
   double distortion_score_non_bonded_contact_lennard_jones(const simple_restraint &plane_restraint,
                                                            const double &lennard_jones_epsilon,
                                                            const gsl_vector *v);
   double distortion_score_parallel_planes(const simple_restraint &plane_restraint,
                                           const gsl_vector *v);
   void fix_chiral_atom_maybe (const simple_restraint &chiral_restraint,
                               gsl_vector *v);
   void fix_chiral_atom_internal (const simple_restraint &chiral_restraint,
                                  gsl_vector *v);

   plane_distortion_info_t
   distortion_score_plane_internal(const simple_restraint &plane_restraint,
                                   const gsl_vector *v,
                                   bool calculate_distortion_flag);
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
   // GM terms
   void my_df_geman_mcclure_distances_old(const gsl_vector *v, void *params, gsl_vector *df);
   void my_df_geman_mcclure_distances(const gsl_vector *v, void *params, gsl_vector *df); // possible multi-thread
   // just the angle terms:
   void my_df_angles(const gsl_vector *v, void *params, gsl_vector *df);
   //  just the torsion terms:
   void my_df_torsions(const gsl_vector *v, void *params, gsl_vector *df);
   void my_df_torsions_internal(const gsl_vector *v, void *params, gsl_vector *df, bool do_rama_torsions);
   void my_df_trans_peptides(const gsl_vector *v, void *params, gsl_vector *df);
   //  just the ramachandran plot gradient terms:
   void my_df_rama(const gsl_vector *v, void *params, gsl_vector *df);
   //  the plane deviation from terms:
   void my_df_planes(const gsl_vector *v, void *params, gsl_vector *df);
   //  the non-bonded contacts
   void my_df_non_bonded(const gsl_vector *v, void *params, gsl_vector *df);
   //
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
   inline double f_inv_fsqrt(const double &v) {
      //
      return 1.0/sqrt(v);
   }

   // debugging function.
   // v needs to be non-const, because gsl_vector_set().
   // if gradients_file_name is not of length 0, then
   // write the gradients to the given file and not to the screen.
   void
   numerical_gradients(gsl_vector *v, void *params, gsl_vector *df,
                       std::string file_name=std::string());

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
   inverted_chiral_volumes(int imol, mmdb::Manager *mol, protein_geometry *geom_p,
                           int cif_dictionary_read_number);

   // make restraints and get distortions. chiral_volume_limit_for_outlier
   // should/might be about 2.0.
   // the spec for the chiral atom and its distortion
   std::pair<std::vector<std::string> , std::vector<std::pair<atom_spec_t, double> > >
   distorted_chiral_volumes(int imol, mmdb::Manager *mol, protein_geometry *geom_p,
                            int cif_dictionary_read_number,
                            double chiral_volume_limit_for_outlier);


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
         int n_improper_dihedral_restr;
         restraint_counts_t() {
            n_bond_restraints = 0;
            n_angle_restraints = 0;
            n_plane_restraints =0;
            n_chiral_restr = 0;
            n_torsion_restr = 0;
            n_improper_dihedral_restr = 0;
         }
         void operator+=(const restraint_counts_t &r) {
            n_bond_restraints += r.n_bond_restraints;
            n_angle_restraints += r.n_angle_restraints;
            n_plane_restraints += r.n_plane_restraints;
            n_chiral_restr += r.n_chiral_restr;
            n_torsion_restr += r.n_torsion_restr;
            n_improper_dihedral_restr += r.n_improper_dihedral_restr;
         }
         void report(bool do_residue_internal_torsions) const;
      };

   private:

      std::vector<simple_restraint> restraints_vec;
      int n_atoms;
      int n_atoms_limit_for_nbc; // the neighbours in non_bonded_contacts_atom_indices are only useful
                                 // for the moving atoms.
      mmdb::PPAtom atom;
      bool atom_array_needs_to_be_deleted_at_end;
      bool model_has_hydrogen_atoms;
      std::vector<bool> atom_is_metal;
      std::vector<bool> atom_is_hydrogen;
      std::vector<int>  old_atom_index_to_new_atom_index;
      void fill_old_to_new_index_vector(); // for fast extra bond (GM) restraints
      bool from_residue_vector;
      int SelHnd_atom; // the selection handle for the atom array.
                       // Note to self: when restraints_container_t
                       // goes out of scope, we should do a
                       // mol->DeleteSelection(SelHnd_atom).


      gsl_multimin_fdfminimizer *m_s;
      double m_initial_step_size;
      double m_tolerance;
      double *par;
      double m_grad_lim;
      gsl_vector *x; // these are the variables, x_k, y_k, z_k, x_l etc.
      bool are_all_one_atom_residues;
      mmdb::Manager *mol;
      void setup_minimize();
      unsigned int n_refiners_refining;
      bool needs_reset; // needs reset when an atom pull restraint gets *added*.

      // The bool is the "atoms of this residue are fixed" flag.
      std::vector<std::pair<bool,mmdb::Residue *> > residues_vec; // these are sorted in
                                                                  // the constructor
      std::set<mmdb::Residue *> residues_vec_moving_set;
      std::map<mmdb::Residue *, std::set<mmdb::Residue *> > fixed_neighbours_set;
      void debug_sets() const;
      int udd_bond_angle;  // for is a bond, angle or not (0).
      int udd_atom_index_handle; // for indexing into the atoms array.

      // and the new (7Nov2003) residue selection (used in non-bonded stuff)
      //
      mmdb::PPResidue SelResidue_active;
      int nSelResidues_active;
      bool apply_H_non_bonded_contacts;

      void init() {
         verbose_geometry_reporting = NORMAL;
         n_refiners_refining = 0;
         n_atoms = 0;
         n_atoms_limit_for_nbc = 0; // needs to be set in every constructor.
         x = 0;
         mol = 0;
         n_atoms = 0;
         atom = 0;
         atom_array_needs_to_be_deleted_at_end = false; // it's not allocated
         model_has_hydrogen_atoms = true;
         include_map_terms_flag = 0;
         have_oxt_flag = 0;
         do_numerical_gradients_flag = 0;
         n_threads = 0;
         apply_H_non_bonded_contacts = true;
         lograma.init(LogRamachandran::All, 2.0, true);
         // when zo_rama is a static, this is already done
         //          try {
         //             zo_rama.init();
         //          }
         //          catch (const std::runtime_error &rte) {
         //             std::cout << "ERROR:: ZO Rama tables failed. " << rte.what() << std::endl;
         //          }
         from_residue_vector = 0;
         rama_type = RAMA_TYPE_LOGRAMA;
         rama_plot_weight = 40.0;
         do_hydrogen_atom_refinement = false;
         do_neutron_refinement = false;

         refinement_results_add_details = true;
         torsion_restraints_weight = 1.0;

#ifndef __NVCC__
         restraints_lock = false; // not locked
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
         thread_pool_p = 0; // null pointer
#endif
#endif // __NVCC__
      }


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
      //
      // sometimes we want the return value without printing the table
      // (hence print_table_flag) - perhaps that should be its own function.
      std::vector<refinement_lights_info_t>
      chi_squareds(std::string title, const gsl_vector *v, bool print_table_flag=true) const;

      // like above but not quite. "Energy" scoring
      void distortion_score_each_restraint(const gsl_vector *v) const;

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
      static zo::rama_table_set zo_rama;
      double rama_plot_weight; // get_rama_plot_weight() is public
      // rama_type is public

      // internal function, most of the job of the constructor:
      void init_from_mol(int istart_res_in, int iend_res_in,
                         bool have_flanking_residue_at_start,
                         bool have_flanking_residue_at_end,
                         short int have_disulfide_residues,
                         const std::string &altloc,
                         const std::string &chain_id,
                         mmdb::Manager *mol_in,
                         const std::vector<atom_spec_t> &fixed_atom_specs);

      void init_from_residue_vec(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
                                 const protein_geometry &geom,
                                 mmdb::Manager *mol_in,
                                 const std::vector<atom_spec_t> &fixed_atom_specs);

      void set_has_hydrogen_atoms_state(); // unset has hydrogens if there are none in the atoms.
      void init_shared_pre(mmdb::Manager *mol_in);

      void init_shared_post(const std::vector<atom_spec_t> &fixed_atom_specs);
      void set_fixed_during_refinement_udd(); // uses fixed_atom_indices
      // the unset function is public - there is no destructor of simple_restraint (hmm)
      // void unset_fixed_during_refinement_udd(); // unsets atoms that were maked as fixed using FixedDuringRefinement
      // neighbour residues already are fixed.
      void add_fixed_atoms_from_flanking_residues(const bonded_pair_container_t &bpc);
      void add_fixed_atoms_from_non_bonded_neighbours(); // use non_bonded_neighbour_residues

      // 20171012 - to help with debugging gradients, we want to know what the fixed
      // atoms are in the atom list when we have a linear residue selection. So
      // set them with the function - called from init_from_mol().
      void add_fixed_atoms_from_flanking_residues(bool have_flanking_residue_at_start,
                                                  bool have_flanking_residue_at_end,
                                                  int iselection_start_res, int iselection_end_res);

      // man this is tricky
      //
      // So. We get a crash when trying to use the const ref xmap for electron density
      // score - it seems as though it is going out of scope (but actually, can't be
      // of course).
      // Maybe it's because the clipper libraries with which I am linking this are
      // compiler with an older compiler?

      // 20180131:
      // this xmap seems to go out of scope if it's a const reference
      // (AFAICS it *can't* go out of scope), but there is a crash.
      // when we try to use pull restraints
      const clipper::Xmap<float> *xmap_p;

      double map_weight;

      double torsion_restraints_weight;

      void add_h_bond(restraint_type_t rest_type, int atom_1, int atom_2,
                      const std::vector<bool> &fixed_atom_flags,
                      float tar, float sig) {

         if (sig > 0.0) {
            float obs = -1; // dummy value
            simple_restraint r(rest_type, atom_1, atom_2, fixed_atom_flags, tar, sig, obs);
            r.is_hydrogen_bond = true;
            restraints_vec.push_back(r);
         }
      }
      void add(restraint_type_t rest_type, int atom_1, int atom_2,
               const std::vector<bool> &fixed_atom_flags,
               float tar, float sig, float obs) {

         if (sig > 0.0) {
            simple_restraint r(rest_type, atom_1, atom_2, fixed_atom_flags, tar, sig, obs);
            restraints_vec.push_back(r);
         }
      }

      bool add(restraint_type_t rest_type, int atom_1, int atom_2, int atom_3,
               const std::vector<bool> &fixed_atom_flags,
               float tar,
               float sig, bool is_single_Hydrogen_atom_angle_restraint){

         bool r = 0;
         if (sig > 0.0) {
            restraints_vec.push_back(simple_restraint(rest_type, atom_1, atom_2, atom_3,
                                                      fixed_atom_flags, tar, sig,
                                                      is_single_Hydrogen_atom_angle_restraint));
            r = 1;
         }
         return r;
      }

      bool add(restraint_type_t rest_type, int atom_1, int atom_2,
               int atom_3, int atom_4,
               const std::vector<bool> &fixed_atom_flags,
               float tar, float sig, float torsion_restraint_weight, int periodicty) {

         bool r = 0;
         if (sig > 0.0) {

            restraints_vec.push_back(simple_restraint(rest_type,
                                                      atom_1, atom_2, atom_3, atom_4,
                                                      fixed_atom_flags, tar, sig, torsion_restraint_weight, periodicty));
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
               const std::string &rama_plot_zo_residue_type,
               int atom_1, int atom_2, int atom_3,
               int atom_4, int atom_5,
               const std::vector<bool> &fixed_atom_flag){

         restraints_vec.push_back(simple_restraint(rest_type,
                                                   rama_plot_zo_residue_type,
                                                   atom_1, atom_2, atom_3, atom_4, atom_5,
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
                          const simple_restraint::nbc_function_t &nbcf,
                          const std::string &atom_type_1,
                          const std::string &atom_type_2,
                          bool is_H_non_bonded_contact,
                          const std::vector<bool> &fixed_atom_flag,
                          double dist_min) {

         restraints_vec.push_back(simple_restraint(NON_BONDED_CONTACT_RESTRAINT, nbcf,
                                                   index1, index2,
                                                   is_H_non_bonded_contact,
                                                   fixed_atom_flag, dist_min));
      }

      void add_target_position_restraint(int idx, const atom_spec_t &spec, clipper::Coord_orth &target_pos);

      // construct a restraint and add it to restraints_vec
      //
      // this assumes the sigmas in atom_index_sigma_in are sensible - so the calling
      // routines needs to make sure that this is the case.
      void add_plane(const std::vector<std::pair<int, double> > &atom_index_sigma_in,
                     const std::vector<bool> &fixed_atom_flags) {
         if (! convert_plane_restraints_to_improper_dihedral_restraints_flag)
            restraints_vec.push_back(simple_restraint(PLANE_RESTRAINT,
                                                      atom_index_sigma_in,
                                                      fixed_atom_flags));
         else
            convert_plane_restraints_to_improper_dihedral_restraints(atom_index_sigma_in, fixed_atom_flags);
      }

      void convert_plane_restraints_to_improper_dihedral_restraints(const std::vector<std::pair<int, double> > &atom_index_sigma_in,
                     const std::vector<bool> &fixed_atom_flags);

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

      void add_user_defined_target_position_restraint(restraint_type_t rest_type, int atom_idx,
                                                      const atom_spec_t &spec,
                                                      const clipper::Coord_orth &pos, float weight) {
         // weight not used yet
         simple_restraint r(rest_type, atom_idx, spec, pos);
         r.is_user_defined_restraint = 1;
         restraints_vec.push_back(r);
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

      int get_asc_index(mmdb::Atom *at);

      int add_bonds(int idr, mmdb::PPAtom res_selection,
                    int i_no_res_atoms,
                    mmdb::PResidue SelRes,
                    const protein_geometry &geom);

      restraint_counts_t add_N_terminal_residue_bonds_and_angles_to_hydrogens(mmdb::Residue *residue_p);
      int get_N_index(mmdb::Residue *residue_p) const;
      int get_CA_index(mmdb::Residue *residue_p) const;
      int get_atom_index(const std::string &atom_name_in, mmdb::Residue *residue_p) const;

      int add_angles(int idr, mmdb::PPAtom res_selection,
                     int i_no_res_atoms,
                     mmdb::PResidue SelRes,
                     const protein_geometry &geom);

      int add_torsions(int idr, mmdb::PPAtom res_selection,
                       int i_no_res_atoms,
                       mmdb::PResidue SelRes,
                       const protein_geometry &geom,
                       const double &torsion_restraints_weight);

      bool add_torsion_internal(const coot::dict_torsion_restraint_t &torsion_restraint,
                                mmdb::PPAtom res_selection, int i_no_res_atoms,
                                const double &torsion_restraints_weight);

      bool
      replace_torsion_restraint(const dict_torsion_restraint_t &new_torsion_restraint,
                                mmdb::PPAtom res_selection, int i_no_res_atoms,
                                const std::vector<unsigned int> &torsion_restraint_indices);

      std::vector<unsigned int> make_torsion_restraint_indices_vector() const;


      // helper function for that:
      int get_atom_index_for_restraint_using_alt_conf(const std::string &atom_name,
                                                      const std::string &alt_conf,
                                                      mmdb::PPAtom res_selection, int num_res_atoms) const;


      int add_chirals(int idr, mmdb::PPAtom res_selection,
                      int i_no_res_atoms,
                      mmdb::PResidue SelRes,
                      const protein_geometry &geom);

      int add_planes  (int idr, mmdb::PPAtom res_selection,
                       int i_no_res_atoms,
                       mmdb::PResidue SelRes,
                       const protein_geometry &geom);

      // called by above.
      //
      // this was the way planes were done 2004-2019
      int add_planes_multiatom_eigen(int idr, mmdb::PPAtom res_selection,
                                     int i_no_res_atoms,
                                     mmdb::PResidue SelRes,
                                     const protein_geometry &geom);

      int add_planes_as_improper_dihedrals(int idr, mmdb::PPAtom res_selection,
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
                                         bool do_rama_plot_retraints,
                                         bool do_trans_peptide_restraints);

      // which uses either:
      int make_link_restraints_by_linear(const protein_geometry &geom,
                                         bool do_rama_plot_retraints,
                                         bool do_trans_peptide_restraints);

      bonded_pair_container_t
      make_link_restraints_from_res_vec(const protein_geometry &geom,
                                        bool do_rama_plot_retraints,
                                        bool do_trans_peptide_restraints);

      bonded_pair_container_t
      make_link_restraints_by_distance(const protein_geometry &geom,
                                        bool do_rama_plot_retraints,
                                        bool do_trans_peptide_restraints);

      bonded_pair_container_t
      make_link_restraints_from_links(const protein_geometry &geom);

      // both of which use:
      int make_link_restraints_by_pairs(const protein_geometry &geom,
                                        const bonded_pair_container_t &bonded_residue_pairs,
                                        bool do_trans_peptide_restraints,
                                        std::string link_or_flanking_link_string);

      void add_rama_links(int SelHnd, const protein_geometry &geom);

      void add_rama_links_from_res_vec(const bonded_pair_container_t &bonded_residue_pairs,
                                       const protein_geometry &geom);


      int make_monomer_restraints(int imol,
                                  const protein_geometry &geom,
                                  short int do_residue_internal_torsions);
      //uses:
      int make_monomer_restraints_by_linear(int imol,
                                            const protein_geometry &geom,
                                            bool do_residue_internal_torsions);
      int make_monomer_restraints_from_res_vec(int imol,
                                               const protein_geometry &geom,
                                               bool do_residue_internal_torsions);
      // which use:
      restraint_counts_t make_monomer_restraints_by_residue(int imol,
                                                            mmdb::Residue *residue_p,
                                                            const protein_geometry &geom,
                                                            bool do_residue_internal_torsions);

      bonded_pair_container_t
      make_flanking_atoms_restraints(const protein_geometry &geom,
                                     bool do_rama_plot_retraints,
                                     bool do_trans_peptide_restraints); // no torsions
      // uses the following
      bonded_pair_container_t
      bonded_flanking_residues(const protein_geometry &geom) const;

      bonded_pair_container_t bonded_flanking_residues_by_residue_vector(const protein_geometry &geom) const;
      // new flanking residue search
      bonded_pair_container_t bonded_flanking_residues_by_residue_vector(const std::map<mmdb::Residue *, std::set<mmdb::Residue *> > &resm,
                                                                         const protein_geometry &geom) const;
      // old style linear search (n +/- 1) selection for flanking residues
      bonded_pair_container_t bonded_flanking_residues_by_linear(const protein_geometry &geom) const;
      // find residues in the neighbourhood that are not in the refining set
      // and are not already marked as bonded flankers.
      //
      std::vector<mmdb::Residue *> non_bonded_neighbour_residues;
      // set by this function:
      // old version 20180224
      void set_non_bonded_neighbour_residues_by_residue_vector(const bonded_pair_container_t &bonded_flanking_pairs,
                                                               const protein_geometry &geom);
      // new version
      void set_non_bonded_neighbour_residues_by_residue_vector(const std::map<mmdb::Residue *, std::set<mmdb::Residue *> > &resm,
                                                               const bonded_pair_container_t &bonded_flanking_pairs,
                                                               const protein_geometry &geom);

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
      // the first value means the following:
      // 1: is
      // 0: is not
      // -1: can't decide here (if that;s the case, use peptide_C_and_N_are_close_p())
      //
      enum peptide_order_info_t { IS_PEPTIDE=1, IS_NOT_PEPTIDE=0, UNKNOWN=-1 };
      //
      std::pair<peptide_order_info_t, bool> peptide_C_and_N_are_in_order_p(mmdb::Residue *r1, mmdb::Residue *r2) const;


      // a pair, first is if C and N are close and second if and order
      // switch is needed to make it so.
      //
      // return "" as first if no close link found.
      //
      std::pair<std::string, bool> general_link_find_close_link(const std::vector<chem_link> &li,
                                                                mmdb::Residue *r1, mmdb::Residue *r2,
                                                                bool order_switch_flag,
                                                                const protein_geometry &geom) const;

      std::string general_link_find_close_link_inner(const std::vector<chem_link> &li,
                                                     mmdb::Residue *r1, mmdb::Residue *r2,
                                                     bool order_switch_flag,
                                                     const protein_geometry &geom) const;



      // restraint_addition_mode can be AUTO_HELIX - restrain anything that looks like a helix (alpha currently)
      // // or EVERYTHING_HELICAL - add helix restrains to residue with the same chain id and in a residue range
      // // that matches a H-bonded residue pair of a helix.
      //
      enum restraint_addition_mode_t { AUTO_HELIX, EVERYTHING_HELICAL};
      void make_helix_pseudo_bond_restraints();
      void make_strand_pseudo_bond_restraints();
      void make_helix_pseudo_bond_restraints_from_res_vec();
      void make_helix_pseudo_bond_restraints_from_res_vec_auto();
      void make_h_bond_restraints_from_res_vec_auto(const protein_geometry &geom, int imol);

      bool link_infos_are_glycosidic_by_name_p(const std::vector<chem_link> &link_infos) const;

      // return "" on failure to find link
      std::string find_glycosidic_linkage_type(mmdb::Residue *first, mmdb::Residue *second,
                                               const protein_geometry &geom,
                                               bool use_links_in_molecule) const;

      // -------------------------- ng restraint generation ------------------------

      class link_restraints_counts {
         void init() {
            n_link_bond_restr = 0;
            n_link_angle_restr = 0;
            n_link_trans_peptide = 0;
            n_link_torsion_restr = 0;
            n_link_plane_restr = 0;
            n_link_improper_dihedral_restr = 0;
            link_type = "link";
         }
      public:
         link_restraints_counts() {
            init();
         }
         explicit link_restraints_counts(const std::string &s) {
            init();
            link_type = s; // e.g. "flank"
         }
         std::string link_type;
         unsigned int n_link_bond_restr;
         unsigned int n_link_angle_restr;
         unsigned int n_link_plane_restr;
         unsigned int n_link_torsion_restr;
         unsigned int n_link_trans_peptide;
         unsigned int n_link_improper_dihedral_restr;
         void add(const link_restraints_counts &lrc) {
            n_link_bond_restr    += lrc.n_link_bond_restr;
            n_link_angle_restr   += lrc.n_link_angle_restr;
            n_link_plane_restr   += lrc.n_link_plane_restr;
            n_link_trans_peptide += lrc.n_link_trans_peptide;
            n_link_torsion_restr += lrc.n_link_torsion_restr;
            n_link_improper_dihedral_restr += lrc.n_link_improper_dihedral_restr;
         }
         void report() const;
      };

      class reduced_angle_info_container_t {
      public:
         reduced_angle_info_container_t() {}
         reduced_angle_info_container_t(const std::vector<simple_restraint> &r);
         reduced_angle_info_container_t(const std::vector<std::vector<simple_restraint> > &rvv); // needs init() also?
         void init(const std::vector<simple_restraint> &r);
         std::map<int, std::set<int> > bonds;
         std::map<int, std::vector<std::pair<int, int> > > angles;
         bool is_1_4(int i, int j, const std::vector<bool> &fixed_atom_flags) const;
         void write_angles_map(const std::string &file_name) const;
      };

      reduced_angle_info_container_t raic;
      void
      make_link_restraints_ng(const protein_geometry &geom,
                              bool do_rama_plot_retraints,
                              bool do_trans_peptide_restraints,
                              std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > *residue_link_vector_map_p,
                              std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > *residue_pair_link_set_p);
      void make_polymer_links_ng(const protein_geometry &geom,
                                 bool do_rama_plot_restraints,
                                 bool do_trans_peptide_restraints,
                                 std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > *residue_link_count_map_p,
                                 std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > *residue_pair_link_set_p);
      bool N_and_C_are_close_ng(mmdb::Residue *res_1, mmdb::Residue *res_2, float d_crit) const;
      bool O3prime_and_P_are_close_ng(mmdb::Residue *res_1, mmdb::Residue *res_2, float d_crit) const;

      std::pair<bool, link_restraints_counts> try_make_peptide_link_ng(const coot::protein_geometry &geom,
                                                                       std::pair<bool, mmdb::Residue *> res_1,
                                                                       std::pair<bool, mmdb::Residue *> res_2,
                                                                       bool do_trans_peptide_restraints);
      std::pair<bool, link_restraints_counts> try_make_phosphodiester_link_ng(const coot::protein_geometry &geom,
                                                                              std::pair<bool, mmdb::Residue *> res_1,
                                                                              std::pair<bool, mmdb::Residue *> res_2);
      link_restraints_counts make_other_types_of_link(const coot::protein_geometry &geom,
                                                      const std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > &residue_link_count_map,
                                                      const std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > &residue_pair_link_set);
      link_restraints_counts make_link_restraints_for_link_ng(const std::string &link_type,
                                                              mmdb::Residue *res_1,
                                                              mmdb::Residue *res_2,
                                                              bool is_fixed_first_residue,
                                                              bool is_fixed_second_residue,
                                                              bool do_trans_peptide_restraints,
                                                              const protein_geometry &geom);

      link_restraints_counts make_link_restraints_for_link_ng(const new_linked_residue_t &nlr,
                                                              const protein_geometry &geom);
      void make_header_metal_links_ng(const protein_geometry &geom);
      void add_header_metal_link_bond_ng(const atom_spec_t &atom_spec_1,
                                         const atom_spec_t &atom_spec_2,
                                         double  dist);


      std::string find_peptide_link_type_ng(mmdb::Residue *res_1,
                                            mmdb::Residue *res_2,
                                            const coot::protein_geometry &geom) const;

      unsigned int  make_non_bonded_contact_restraints_ng(int imol, const protein_geometry &geom);
      void make_non_bonded_contact_restraints_using_threads_ng(int imol, const protein_geometry &geom);
      std::vector<std::set<int> > non_bonded_contacts_atom_indices; // these can now get updated on the fly.
                                                       // We need to keep a record of what has
                                                       // already been added as a restraint
                                                       // before we add a new one.

      void analyze_for_bad_restraints(restraint_type_t r_type, double interesting_distortion_limit);
#ifndef __NVCC__
      // threaded workpackage
      static
      void make_non_bonded_contact_restraints_workpackage_ng(int ithread,
                                                             int imol,
                                                             const coot::protein_geometry &geom,
                                                             const std::vector<std::set<int> > &bonded_atom_indices,
                                                             const reduced_angle_info_container_t &raic,
                                                             const std::vector<std::set<unsigned int> > &vcontacts,
                                                             std::pair<unsigned int, unsigned int> atom_index_range_pair,
                                                             const std::set<int> &fixed_atom_indices,
                                                             const std::vector<std::string> &energy_type_for_atom,
                                                             bool use_extended_atom_mode,
                                                             mmdb::PPAtom atom,
                                                             const std::vector<bool> &atom_is_metal,
                                                             const std::vector<bool> &atom_is_hydrogen,
                                                             const std::vector<bool> &H_atom_parent_atom_is_donor_vec,
                                                             const std::vector<bool> &atom_is_acceptor_vec,
                                                             std::vector<std::set<int> > *non_bonded_contacts_atom_indices_p,
                                                             std::vector<simple_restraint> *nbc_restraints_fragment_p,
                                                             std::atomic<unsigned int> &done_count);
#endif

      // update residue_link_vector_map_p and residue_pair_link_set if new links are made
      //
      void make_flanking_atoms_restraints_ng(const coot::protein_geometry &geom,
                                             std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > *residue_link_vector_map_p,
                                             std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > *residue_pair_link_set_p,
                                             bool do_rama_plot_restraints,
                                             bool do_trans_peptide_restraints);

      void make_base_pairing_and_stacking_restraints_ng(int imol, const protein_geometry &geom);

      void make_rama_plot_restraints(const std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > &residue_link_vector_map,
                                     const std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > &residue_pair_link_set,
                                     const protein_geometry &geom);

      void make_rama_plot_restraints_ng(const std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > &residue_link_vector_map,
                                        const std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > &residue_pair_link_set,
                                        const protein_geometry &geom); // uses restraints_vec

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
                           mmdb::Residue *first,
                           mmdb::Residue *second,
                           short int is_fixed_first,
                           short int is_fixed_second,
                           const protein_geometry &geom);

      int add_link_torsion_for_phi_psi(std::string link_type,
                                       int phi_psi_restraints_type,
                                       mmdb::Residue *first,
                                       mmdb::Residue *second,
                                       short int is_fixed_first,
                                       short int is_fixed_second,
                                       const protein_geometry &geom);

      // a strong restraint on w to make it 180, period 1.
      //
      int add_link_trans_peptide(mmdb::Residue *first,
                                 mmdb::Residue *second,
                                 bool is_fixed_first,
                                 bool is_fixed_second,
                                 bool add_even_if_cis);

      int add_rama(std::string link_type,
                   mmdb::PResidue prev,
                   mmdb::PResidue this_res,
                   mmdb::PResidue post,
                   bool is_fixed_first_res,
                   bool is_fixed_second_res,
                   bool is_fixed_third_res,
                   const protein_geometry &geom);

      int add_rama(const rama_triple_t &rt, const coot::protein_geometry &geom);

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
      // bonded_atom_indices also contain 1-3 atoms of angles
      std::vector<std::set<int> > bonded_atom_indices;

      bool check_for_1_4_relation(int i, int j) const;
      bool check_for_1_4_relation(int i, int j, const reduced_angle_info_container_t &ai) const;
      static bool check_for_O_C_1_5_relation(mmdb::Atom *at_1, mmdb::Atom *at_2);  // check either way round


      int make_non_bonded_contact_restraints(int imol, const bonded_pair_container_t &bpc, const protein_geometry &geom);
      int make_non_bonded_contact_restraints(int imol,
                                             const bonded_pair_container_t &bpc,
                                             const reduced_angle_info_container_t &ai,
                                             const protein_geometry &geom);
      static bool is_in_same_ring(int imol, mmdb::Residue *residue_p,
                                  std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > &residue_ring_map_cache,
                                  const std::string &atom_name_1,
                                  const std::string &atom_name_2,
                                  const coot::protein_geometry &geom);

      bool is_acceptor(const std::string &energy_type,
                       const coot::protein_geometry &geom) const;


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

      bool do_hydrogen_atom_refinement;
      bool do_neutron_refinement;
      void set_do_hydrogen_atom_refinement(bool state) { do_hydrogen_atom_refinement = state; }
      void set_do_neutron_refinement(bool state) {
         do_neutron_refinement = state;
         if (state) set_z_occ_weights(); // needs neutron weights now
      }

      // validation:
      geometry_distortion_info_container_t
      distortion_vector(const gsl_vector *v, bool keep_distortion_for_hydrogen_atom_restraintsb) const;

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
      int get_chiral_hydrogen_index(int indexc, int index1, int index_2, int index_3, const dict_chiral_restraint_t &dcr) const;
      bool has_inverted_chiral_centre(const simple_restraint &chiral_restraint,
                                      const gsl_vector *v) const;
      bool has_tiny_chiral_centre_volume(const simple_restraint &chiral_restraint,
                                         const gsl_vector *v) const;
      bool bond_is_very_long(const simple_restraint &bond_restraint,
                             const gsl_vector *v) const;

      bool is_hydrogen(mmdb::Atom *at_p) const {
         std::string ele = at_p->element;
         if ((ele == "H") || (ele == " H"))
            return true;
         else
            return ((ele == "D") || (ele == " D"));
         }

      // return "" on no type found
      std::string get_type_energy(int imol, mmdb::Atom *at, const protein_geometry &geom) const {
         std::string r;
         if (at) {
            std::string atom_name = at->name;
            const char *rn = at->GetResName();
            if (rn) {
               std::string residue_name = rn;
               r = geom.get_type_energy(atom_name, residue_name, imol);
            }
         }
         return r;
      }

      bonded_pair_container_t bonded_pairs_container;

      void pre_sanitize_as_needed(std::vector<refinement_lights_info_t> lights);

      model_bond_deltas resolve_bonds(const gsl_vector *v) const;

      void make_restraint_types_index_limits();

      float dist_crit_for_bonded_pairs;

      // regenerate the restraints_indices df_by_thread_results
      //
      void post_add_new_restraint();
      void post_add_new_restraints();

      // return false if any of the atoms are fixed
      bool none_are_fixed_p(const std::vector<bool> &fixed_atom_indices) const;

      unsigned int n_times_called; // so that we can do certain things only the first time
      unsigned int n_small_cycles_accumulator;

      // what is the energy type of the atom to which the Hydrogen atom is bonded?
      //
      std::map<mmdb::Atom *, hb_t> H_atom_parent_energy_type_atom_map;
      bool H_parent_atom_is_donor(mmdb::Atom *at); // adds to the above map potentially

      std::vector<mmdb::Link> links; // worry about deleting these. Are they just shall copies?
      void fill_links(mmdb::Manager *mol); // if they were not passed in the constructor.

   public:

      enum link_torsion_restraints_type { NO_LINK_TORSION = 0,
                                          LINK_TORSION_RAMACHANDRAN_GOODNESS = 1,
                                          LINK_TORSION_ALPHA_HELIX = 2,
                                          LINK_TORSION_BETA_STRAND = 3 };

      // my_df_electron_density and electron_density_score need access
      // to fixed_atom_indices.
      std::set<int> fixed_atom_indices;

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
//          verbose_geometry_reporting = 0;
//          atom = atoms_in;
//          n_atoms = n_at;
//          asc.mol = NULL;
//          include_map_terms_flag = 0;
//          initial_position_params_vec.resize(3*n_at);
//          for (int i=0; i<n_at; i++) {
//             initial_position_params_vec[3*i  ] = atoms_in[i]->x;
//             initial_position_params_vec[3*i+1] = atoms_in[i]->y;
//             initial_position_params_vec[3*i+2] = atoms_in[i]->z;
//          }
//       }

      restraints_container_t(atom_selection_container_t asc_in, const clipper::Xmap<float> *xmap_p_in)
         : xmap_p(xmap_p_in) {

         // xmap = xmap_in;

         init();
         mol = asc_in.mol;
         n_atoms = asc_in.n_selected_atoms;
         atom = asc_in.atom_selection;
         initial_position_params_vec.resize(3*asc_in.n_selected_atoms);
         dist_crit_for_bonded_pairs = 3.0;

         for (int i=0; i<asc_in.n_selected_atoms; i++) {
            initial_position_params_vec[3*i  ] = asc_in.atom_selection[i]->x;
            initial_position_params_vec[3*i+1] = asc_in.atom_selection[i]->y;
            initial_position_params_vec[3*i+2] = asc_in.atom_selection[i]->z;
         }
      }

      // for omega distortion info:
      restraints_container_t(atom_selection_container_t asc_in,
                             const std::string &chain_id,
                             const clipper::Xmap<float> *xmap_in);


      // For validation.
      //
      // The whole chain is selected (without flanking atoms) and we
      // use geometry_distortion() function.
      //
      restraints_container_t(mmdb::PResidue *SelResidues, int nSelResidues,
                             const std::string &chain_id,
                             mmdb::Manager *mol,
                             const clipper::Xmap<float> *xmap_p_in);

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
                             const std::vector<atom_spec_t> &fixed_atom_specs,
                             const clipper::Xmap<float> *xmap_p_in);

      restraints_container_t(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
                             const protein_geometry &geom,
                             mmdb::Manager *mol,
                             const clipper::Xmap<float> *xmap_p_in);

      unsigned int df_by_thread_results_size() const;

      //
      // geometric_distortions not const because we set restraints_usage_flag:
      //
      // return data useful for making the graphs:
      //
      // 20181231 - yeah, I don't think that reseting restraints_usage_flags is a
      //            good idea. Remove this function.
      // geometry_distortion_info_container_t
      // geometric_distortions(restraint_usage_Flags flags);

      // Here we use the internal flags.  Causes crash currently (no inital atom positions?)
      //
      geometry_distortion_info_container_t geometric_distortions(bool keep_distortion_for_hydrogen_atom_restraints=true);

      // Here we use the internal flags.
      //
      geometry_distortion_info_pod_container_t geometric_distortions_pod(bool include_distortion_for_hydrogen_atom_restraints=true);

      omega_distortion_info_container_t
      omega_trans_distortions(const protein_geometry &geom,
                              bool mark_cis_peptides_as_bad_flag);


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
                             const std::set<int> &fixed_atom_indices,
                             clipper::Xmap<float> *map_p_in);

      explicit restraints_container_t(const clipper::Xmap<float> *map_p_in) : xmap_p(map_p_in) {
         init();
         from_residue_vector = 0;
         include_map_terms_flag = 0;

      };

      ~restraints_container_t();

      void unset_fixed_during_refinement_udd(); // unsets atoms that were maked as fixed using FixedDuringRefinement
                                                // now called in the destructor

      mmdb::Atom *get_atom(int i) const {
         if (atom)
            return atom[i];
         else
            return NULL;
      }

      atom_spec_t get_atom_spec(int atom_index) const;

      // rama_type_in is either RAMA_TYPE_ZO or RAMA_TYPE_LOGRAMA
      //
      void set_rama_type(int rama_type_in) {
         rama_type = rama_type_in;
      }

      void set_rama_plot_weight(float w) {
         rama_plot_weight = w;
      }

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
      int restraints_usage_flag;

      double starting_structure_diff_score(const gsl_vector *v, void *params);

      bool apply_H_non_bonded_contacts_state() const { return apply_H_non_bonded_contacts; }

      void set_apply_H_non_bonded_contacts(bool state) { apply_H_non_bonded_contacts = state; }

      short int include_map_terms() {
         return include_map_terms_flag;
      }

      double Map_weight() const {
         return map_weight;
      }

      void set_map_weight(const double &mw) {
         map_weight = mw;
      }

      void set_torsion_restraints_weight(double w);

      double get_torsion_restraints_weight() const { return torsion_restraints_weight; }

      void setup_multimin_func() {

         multimin_func.f   = &distortion_score;
         multimin_func.df  = &my_df;
         multimin_func.fdf = &my_fdf;
         multimin_func.n = n_variables();
         multimin_func.params = (double *) this;
      }

#ifndef __NVCC__
      // we should not update the atom pull restraints while the refinement is running.
      // we shouldn't refine when the atom pull restraints are being updated.
      // we shouldn't clear the gsl_vector x when o
      std::atomic<bool> restraints_lock;
      void get_restraints_lock();
      void release_restraints_lock();
      static std::atomic<bool> print_lock;
      static void get_print_lock();
      static void release_print_lock();
#endif

      unsigned int get_n_atoms() const { return n_atoms; } // access from split_the_gradients_with_threads()
      unsigned int n_variables() const {
         // return 3 * the number of atoms
         return 3*n_atoms;
      }

      // no longer done in the constructor:
      void add_map(float map_weight_in) {
         map_weight = map_weight_in;
         include_map_terms_flag = true;
      }

      int size() const { return restraints_vec.size(); }

      // return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
      //
      // We now have access to n_times_called: we want to do pre-sanitization
      // only when n_times_called is 1.
      //
      // print_chi_sq_flag is not used, so one of these functions is redundant.
      refinement_results_t minimize(restraint_usage_Flags, int n_steps_max = 1000);
      refinement_results_t minimize(restraint_usage_Flags, int nsteps, short int print_chi_sq_flag);
      refinement_results_t minimize(int imol, restraint_usage_Flags usage_flags,
                                    int nsteps_max, short int print_initial_chi_sq_flag,
                                    const protein_geometry &geom);
      refinement_results_t minimize_inner(restraint_usage_Flags, int nsteps);

      refinement_results_t get_refinement_results(); // not const because setup_minimize()

      void simulated_annealing();

      void free_delete_reset();

      bool refinement_results_add_details;
      void add_details_to_refinement_results(refinement_results_t *rr) const;

      void fix_chiral_atoms_maybe(gsl_vector *s);

      refinement_lights_info_t::the_worst_t
      find_the_worst(const std::vector<refinement_lights_info_t> &lights) const;

      simple_restraint& operator[] (unsigned int i) {
         return restraints_vec[i];
      }

      // because chi_squareds is const:
      const simple_restraint& operator[] (const unsigned int &i) const {
         return restraints_vec[i];
      }

      const simple_restraint& at(const unsigned int &i) const{
         return restraints_vec[i];
      }

      void setup_gsl_vector_variables();

      double electron_density_score_at_point(const clipper::Coord_orth &ao) const;
      clipper::Grad_orth<double> electron_density_gradient_at_point(const clipper::Coord_orth &ao) const;

      // sometime we want to anneal bonds for atoms that are "far" apart. Default
      // distance is 3.0.
      void set_dist_crit_for_bonded_pairs(float dist);

      // We need to fill restraints_vec (which is a vector of
      // simple_restraint) using the coordinates () and the dictionary of
      // restraints, protein_geometry geom.
      int make_restraints(int imol,
                          const protein_geometry &geom,
                          restraint_usage_Flags flags,
                          bool do_residue_internal_torsions,
                          bool do_trans_peptide_restraints,
                          float rama_plot_target_weight,
                          bool do_rama_plot_retraints,
                          bool do_auto_helix_restraints,
                          bool do_auto_strand_restraints,
                          bool do_auto_h_bond_restraints,
                          pseudo_restraint_bond_type sec_struct_pseudo_bonds,
                          bool do_link_restraints=true,
                          bool do_flank_restraints=true);

      int make_restraints_ng(int imol,
                          const protein_geometry &geom,
                          restraint_usage_Flags flags,
                          bool do_residue_internal_torsions,
                          bool do_trans_peptide_restraints,
                          float rama_plot_target_weight,
                          bool do_rama_plot_retraints,
                          bool do_auto_helix_restraints,
                          bool do_auto_strand_restraints,
                          bool do_auto_h_bond_restraints,
                          pseudo_restraint_bond_type sec_struct_pseudo_bonds,
                          bool do_link_restraints=true,
                          bool do_flank_restraints=true);

      bool add_or_replace_torsion_restraints_with_closest_rotamer_restraints(const std::vector<std::pair<mmdb::Residue *, std::vector<dict_torsion_restraint_t> > > &rotamer_torsions);

      unsigned int test_function(const protein_geometry &geom);
      unsigned int inline_const_test_function(const protein_geometry &geom) const {
         std::cout << "----- inline_const_test_function() with geom of size : " << geom.size()
                   << std::endl;
         std::cout << "    geom ref pointer " << &geom << std::endl;
         return geom.size();
      }
      unsigned int const_test_function(const protein_geometry &geom) const;

      mmdb::Atom *add_atom_pull_restraint(const atom_spec_t &spec, clipper::Coord_orth pos);
      void clear_atom_pull_restraint(const atom_spec_t &spec); // clear any previous restraint for this atom.
      void clear_all_atom_pull_restraints();
      unsigned int n_atom_pull_restraints() const; // counts closed pull restraints also.

      void add_extra_restraints(int imol,
                                const std::string &description,
                                const extra_restraints_t &extra_restraints,
                                const protein_geometry &geom);
      // and that calls:
      void add_extra_geman_mcclure_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_bond_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_angle_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_torsion_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_start_pos_restraints(const extra_restraints_t &extra_restraints);
      void add_extra_target_position_restraints(const extra_restraints_t &extra_restraints); // not pull-atoms
      void add_extra_parallel_plane_restraints(int imol,
                                               const extra_restraints_t &extra_restraints,
                                               const protein_geometry &geom);
      // can I find the atoms using the atom indices from the original molecule?
      bool try_add_using_old_atom_indices(const extra_restraints_t::extra_bond_restraint_t &ebr);
      // ditto
      bool try_add_using_old_atom_indices(const extra_restraints_t::extra_geman_mcclure_restraint_t &ebr);

      // rama_type is public, maybe instead use get_rama_type()
      enum { RAMA_TYPE_ZO, RAMA_TYPE_LOGRAMA };
      int rama_type;
      float get_rama_plot_weight() const { return rama_plot_weight; }

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
      std::vector<double> atom_z_occ_weight;  // weight e.d. fit by atomic number and occ
      std::map<std::string, double> neutron_occupancy_map;
      void set_z_occ_weights();
      void init_neutron_occupancies();
      double neutron_occupancy(const std::string &element, int formal_charge) const;

      // Make a MMDBManager from the selection and return this new
      // mmdb::PManager.  It is the users responsibility to delete it.
      mmdb::Manager *results() const;
      void adjust_variables(const atom_selection_container_t &asc);

      const LogRamachandran &LogRama() const { return lograma; };

      const zo::rama_table_set &ZO_Rama() const { return zo_rama; };

      // here phi and psi are in clipper units (radians).
      double rama_prob(const double &phi_rads, const double &psi_rads) const {
         return lograma.interp(phi_rads, psi_rads);
      }
      LogRamachandran::Lgrad rama_grad(const double &phir, const double &psir) const {
         return lograma.interp_grad(phir, psir);
      }

      // calling function should also provide the plot type
      // residue type eg "ALL!nP" "ALLnP" "GLY!nP"  "GLYnP" "PRO!nP"
      //
      float zo_rama_prob(const std::string &residue_type, const double &phir, const double &psir) const {
         return zo_rama.value(residue_type, phir, psir);
      }

      // calling function should now also provides the plot/residue type
      //
      std::pair<float, float> zo_rama_grad(const std::string &residue_type,
                                           const double &phir, const double &psir) const {
         return zo_rama.df(residue_type, phir, psir);
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

      // 20221120-PE find_link_type_complicado is too complicated for me to understand.
      // Let's try again
      std::pair<std::string, bool> find_link_type_2022(mmdb::Residue *first,
                                                       mmdb::Residue *second,
                                                       const protein_geometry &geom) const;

      // which calls
      bool have_intermediate_residue_by_seqnum(mmdb::Residue *first,
                                               mmdb::Residue *second) const;

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

      // return a vector of the specs of the restraints  if the restraint was turned off.
      // Never include that atom that the user is dragging.
      //
      std::vector<atom_spec_t> turn_off_atom_pull_restraints_when_close_to_target_position(const atom_spec_t &dragged_atom);

      void pull_restraint_displace_neighbours(mmdb::Atom *at,
                                              const clipper::Coord_orth &new_pull_atom_target_position,
                                              float radius_of_effect);
      bool use_proportional_editing;
      float pull_restraint_neighbour_displacement_max_radius;
      void set_use_proportional_editing(bool state);


      bool cryo_em_mode; // for weighting fit to density of atoms (side-chains and others are down-weighted)

      // more debugging interface:
      //
      void set_do_numerical_gradients() { do_numerical_gradients_flag = 1;}
      bool do_numerical_gradients_status() { return do_numerical_gradients_flag; }
      void debug_atoms() const;

      model_bond_deltas resolve_bonds(); // calls setup_gsl_vector_variables()

      // we set these so that the functions, for particular types, that loop over the restraints
      // need only loop over the relevant range
      // i.e. say we have 2000 restraints, only the top 100 of which might contain
      // bond restraints. make_restraint_types_index_limits() is called at the end of
      // make_restraints()
      //
      // these are public because they are used in the my_df_xxx functions
      //
      std::pair<unsigned int, unsigned int> restraints_limits_bonds;
      std::pair<unsigned int, unsigned int> restraints_limits_angles;
      std::pair<unsigned int, unsigned int> restraints_limits_torsions;
      std::pair<unsigned int, unsigned int> restraints_limits_chirals;
      std::pair<unsigned int, unsigned int> restraints_limits_planes;
      std::pair<unsigned int, unsigned int> restraints_limits_parallel_planes;
      std::pair<unsigned int, unsigned int> restraints_limits_non_bonded_contacts;
      std::pair<unsigned int, unsigned int> restraints_limits_geman_mclure;
      std::pair<unsigned int, unsigned int> restraints_limits_start_pos;
      std::pair<unsigned int, unsigned int> restraints_limits_trans_peptide;
      // std::pair<unsigned int, unsigned int> restraints_limits_target_pos; // atom pull

      // so that each thread gets a more or less similar number of plane restraints.
      // (not needed I think - that happens anyway)
      void disperse_plane_restraints();

      // plane restraint should be sets of 4-atom chiral-like restraints?
      bool convert_plane_restraints_to_improper_dihedral_restraints_flag;
      void set_convert_plane_restraints_to_improper_dihedral_restraints(bool state) {
         convert_plane_restraints_to_improper_dihedral_restraints_flag = state;
      }

      void set_geman_mcclure_alpha(double alpha_in) { geman_mcclure_alpha = alpha_in; }

      double lennard_jones_epsilon; // 0.1 default values

      void set_lennard_jones_epsilon(const double &e) { lennard_jones_epsilon = e; }

      // when we dynamically do a cis-trans conversion on the moving/intermediate atoms
      void add_trans_peptide_restraint(mmdb::Residue *first, mmdb::Residue *second);
      void remove_trans_peptide_restraint(mmdb::Residue *first, mmdb::Residue *second);

      // Is it sane to have threads without a thread pool?
      //
      // I think so - for example crankshaft
      //
      unsigned int n_threads;
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
      // thread pool!
      //
      ctpl::thread_pool *thread_pool_p;
      // std::atomic<unsigned int> &done_count_for_threads;
      void thread_pool(ctpl::thread_pool *tp_in, int n_threads_in) {
         thread_pool_p = tp_in;
         n_threads = n_threads_in;
         // std::cout << "##### thread_pool called with n_thread " << n_threads << std::endl;
      }

      // we can't have a non-pointer thread pool because restraints are copied in
      // update_refinement_atoms() (graphics-info-modelling.cc)
      // and to do that we need a copy operator for thread_pool (and that is
      // deleted in the header).
      //
      // ctpl::thread_pool another_thread_pool;

      void make_df_restraints_indices();
      void make_distortion_electron_density_ranges();
      void clear_df_by_thread_results();

      // generated by make_distortion_electron_density_ranges():
      std::vector<std::pair<unsigned int, unsigned int> > m_atom_index_ranges;

      std::vector<std::vector<double> > df_by_thread_results;
      std::vector<std::vector<std::size_t> > df_by_thread_atom_indices; // for electron density
      // pull restraints are dynamically added to the end of restraints_indices
      std::vector<std::vector<std::size_t> > restraints_indices;

#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

      void clear();

      double log_cosh_target_distance_scale_factor;
      void set_log_cosh_target_distance_scale_factor(double sf) {
         log_cosh_target_distance_scale_factor = sf;
      }

      void copy_from(int i);

      // friend?
      void copy_from(const restraints_container_t &other);

      void set_use_harmonic_approximations_for_nbcs(bool flag);

      // allow the calling function to tell us that the atoms have moved and they
      // need to be updated from the intermediate atoms molecule (moving_atom_asc)
      //
      // e.g. JED refine, cis-trans, pepflip will set this
      //
      void set_needs_reset() { needs_reset = true; }

      double get_distortion_score() const;

      enum analyze_bad_restraints_mode { NORMAL_BAD_RESTRAINTS_ANALYSIS, BAD_RESTRAINT_ANALYSIS_INCLUDE_ANGLES };
      void analyze_for_bad_restraints(analyze_bad_restraints_mode mode = NORMAL_BAD_RESTRAINTS_ANALYSIS);

   };


#ifndef __NVCC__
   void my_df_non_bonded_thread_dispatcher(int thread_idx,
                                           const gsl_vector *v,
                                           gsl_vector *df,
                                           restraints_container_t *restraints_p,
                                           int idx_start,
                                           int idx_end,
                                           std::atomic<unsigned int> &done_count);

   void my_df_geman_mcclure_distances_thread_dispatcher(int thread_index, const gsl_vector *v,
                                                        gsl_vector *df,
                                                        restraints_container_t *restraints_p,
                                                        int idx_start,
                                                        int idx_end,
                                                        std::atomic<unsigned int> &done_count);
   // parallel version of my_df()
   void split_the_gradients_with_threads(const gsl_vector *v,
                                         restraints_container_t *restraints_p,
                                         gsl_vector *df);

   // this should be in process_df_in_range.hh perhaps? Anyway, splitting it up
   // doesn't seem to speed things up.
   void consolidate_derivatives(unsigned int thread_index,
                                unsigned int n_restraints_sets,
                                unsigned int variable_idx_start,
                                unsigned int variable_idx_end,  // stop before this end, e.g. 0, 10
                                const std::vector<std::vector<double> > &df_sets_from,
                                gsl_vector *df,
                                std::atomic<unsigned int> &done_count_for_threads);

#endif


   double electron_density_score(const gsl_vector *v, void *params);
   // interestingly, this needs a different name to the above so that std::async()
   // in distortion_score() refers to the correct function.
   double electron_density_score_from_restraints(const gsl_vector *v, coot::restraints_container_t *restraints_p);
   double electron_density_score_from_restraints_simple(const gsl_vector *v, coot::restraints_container_t *restraints_p);

#ifndef __NVCC__
   // The version of electron_density_score_from_restraints that can be used with a thread pool.
   // The calling function needs to push this onto the queue, one for each thread,
   // where the atom_index_range splits up the atom indices
   // e.g. for 3 threads and 36 atoms, the atom_index_ranges would be:
   // (0, 12)  (12,24) (24,36)
   // but 2001 and 13 threads... errr...?  There is function in split-indices to do this now.
   //
   // atom_index_range works "as expected"
   // so given atom_index_range of 0,10 we start at the first value (0) and check that the
   // current value is less than the atom_index_range.second (10):
   // ie. density values for atom indices 0 to 9 inclusive are added.
   void electron_density_score_from_restraints_using_atom_index_range(int thread_idx,
                                                 const gsl_vector *v,
                                                 const std::pair<unsigned int, unsigned int> &atom_index_range,
                                                 restraints_container_t *restraints_p,
                                                 double *result,
                                                 std::atomic<unsigned int> &done_count);
#endif // __NVCC__

   // new style Grad_map/Grad_orth method
   void my_df_electron_density(const gsl_vector *v, void *params, gsl_vector *df);
   // pre-threaded
   void my_df_electron_density_old_2017(const gsl_vector *v, void *params, gsl_vector *df);
   // old style numerical method
   void my_df_electron_density_old(gsl_vector *v, void *params, gsl_vector *df);


   void my_df_non_bonded_single(const gsl_vector *v,
                                gsl_vector *df,
                                const simple_restraint &this_restraint
                                // const restraints_container_t &restraints // for debugging
                                );

   void my_df_non_bonded_lennard_jones(const gsl_vector *v,
                                       gsl_vector *df,
                                       const simple_restraint &this_restraint,
                                       const double &lj_epsilon);

   void my_df_geman_mcclure_distances_single(const gsl_vector *v,
                                             gsl_vector *df,
                                             const simple_restraint &this_restraint,
                                             const double &alpha);
   void my_df_non_bonded_single(const gsl_vector *v,
                                gsl_vector *df,
                                const simple_restraint &this_restraint
                                // const restraints_container_t &restraints // for debugging
                                );

   void my_df_electron_density_single(const gsl_vector *v,
                                      restraints_container_t *restraints,
                                      gsl_vector *df, int idx_start, int idx_end);

#ifdef HAVE_CXX_THREAD
#ifndef __NVCC__
   // done_count_for_threads is modified
   //
   void my_df_electron_density_threaded_single(int thread_idx, const gsl_vector *v,
                                               restraints_container_t *restraints,
                                               gsl_vector *df,
                                               int atom_idx_start, int atom_idx_end,
                                               std::atomic<unsigned int> &done_count_for_threads);
#endif // __NVCC__
#endif // HAVE_CXX_THREAD

   void simple_refine(mmdb::Residue *residue_p,
                      mmdb::Manager *mol,
                      const dictionary_residue_restraints_t &dict_restraints);

} // namespace coot


#endif // HAVE_SIMPLE_RESTRAINT_HH
