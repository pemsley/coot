/* coot-utils/coot-coord-utils.hh
 *
 * Copyright 2006, 2007, by The University of York
 * Copyright 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2013 by Medical Research Council
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

#ifndef HAVE_COOT_COORD_UTILS_HH
#define HAVE_COOT_COORD_UTILS_HH

#include <stdexcept>
#include <sstream>
#include <set>

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include <map>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

#include "mini-mol/atom-quads.hh"
#include "coot-lsq-types.hh"

// How should I do this better?
#define CXX_UNSET_CHARGE -99.8

#include "geometry/residue-and-atom-specs.hh"

#include "arc-info.hh"

#include "cis-peptide-info.hh"

namespace coot {

   // a generally useful class to be used with std::map where the
   // return type is an index (int).  We use this to determine if the
   // map returns a useful result (if the map used a starndard int, it
   // would return a int (typically 0) when the map failed (no indexer
   // in map).
   //

   // for canonical_base_name();
   enum base_t { RNA, DNA };

   // Perhaps this should be a class function of a class derived from mmdb::Manager?
   int write_coords_pdb(mmdb::Manager *mol, const std::string &file_name);
   int write_coords_cif(mmdb::Manager *mol, const std::string &file_name);

   std::string pad_atom_name(const std::string &atom_name_in,
                             const std::string &element);


   class torsion {
   public:
      // imol torsion pairs
      std::pair<int, atom_spec_t> atom_1;
      std::pair<int, atom_spec_t> atom_2;
      std::pair<int, atom_spec_t> atom_3;
      std::pair<int, atom_spec_t> atom_4;
      torsion(const std::pair<int, atom_spec_t> &atom_1_in,
              const std::pair<int, atom_spec_t> &atom_2_in,
              const std::pair<int, atom_spec_t> &atom_3_in,
              const std::pair<int, atom_spec_t> &atom_4_in) :
         atom_1(atom_1_in),
         atom_2(atom_2_in),
         atom_3(atom_3_in),
         atom_4(atom_4_in) {}
      torsion(int imol,
              const atom_spec_t &atom_1_in,
              const atom_spec_t &atom_2_in,
              const atom_spec_t &atom_3_in,
              const atom_spec_t &atom_4_in) :
         atom_1(std::pair<int, atom_spec_t> (imol, atom_1_in)),
         atom_2(std::pair<int, atom_spec_t> (imol, atom_2_in)),
         atom_3(std::pair<int, atom_spec_t> (imol, atom_3_in)),
         atom_4(std::pair<int, atom_spec_t> (imol, atom_4_in)) {}

      // Find 4 atoms in residue that match the torsion spec.  If the
      // returning vector is not of size 4, then this function has
      // failed.
      std::vector<mmdb::Atom *> matching_atoms(mmdb::Residue *residue);
   };


   // Return dist in Angstroms, can throw an exception if any of the
   // atoms is null.
   //
   double distance(mmdb::Atom *at_1, mmdb::Atom *at_2);

   // Return angle in degrees, can throw an exception if any of the
   // atoms is null.
   //
   double angle(mmdb::Atom *at_1, mmdb::Atom *at_2, mmdb::Atom *at_3);

   // Return torsion angle in degrees, can throw an exception if any of the
   // atoms is null.
   //
   // double torsion() do this by making an atom_quad.

   // must be a valid pointer.
   clipper::Coord_orth co(mmdb::Atom *at);
   //
   void update_position(mmdb::Atom *at, const clipper::Coord_orth &pos);

   class lsq_range_match_info_t {
   public:
      int model_number_reference; // usually 0, meaning unset.
      int model_number_matcher;   // usually 0;
      int to_reference_start_resno;
      int to_reference_end_resno;
      int from_matcher_start_resno;
      int from_matcher_end_resno; // is this used? // for validation?
      std::string reference_chain_id;
      std::string matcher_chain_id;
      int match_type_flag; // CA/Main/All/{N,CA,C}
      bool is_single_atom_match;
      std::string reference_atom_name;
      std::string reference_alt_conf;
      std::string matcher_atom_name;
      std::string matcher_alt_conf;
      lsq_range_match_info_t() {
         is_single_atom_match = false;
         match_type_flag = 0;
         from_matcher_start_resno = -1;
         from_matcher_end_resno = -1;
         to_reference_start_resno = -1;
         to_reference_end_resno = -1;
         model_number_reference = -1;
         model_number_matcher = -1;
      };
      lsq_range_match_info_t(int to_reference_start_resno_in,
                             int to_reference_end_resno_in,
                             const std::string &reference_chain_id_in,
                             int from_matcher_start_resno_in,
                             int from_matcher_end_resno_in,
                             const std::string &matcher_chain_id_in,
                             short int match_type_flag_in) :
         reference_chain_id(reference_chain_id_in), matcher_chain_id(matcher_chain_id_in) {
         is_single_atom_match = 0;
         match_type_flag = match_type_flag_in;
         to_reference_start_resno = to_reference_start_resno_in;
         to_reference_end_resno = to_reference_end_resno_in;
         from_matcher_start_resno = from_matcher_start_resno_in;
         from_matcher_end_resno = from_matcher_end_resno_in;
         model_number_matcher = 0;
         model_number_reference = 0;
      }

      lsq_range_match_info_t(const std::string &reference_chain_id_in,
                             int reference_resno_in,
                             const std::string &reference_insertion_code_in,
                             const std::string &reference_atom_name_in,
                             const std::string &reference_alt_conf_in,
                             const std::string &matcher_chain_id_in,
                             int matcher_resno_in,
                             const std::string &matcher_insertion_code_in,
                             const std::string &matcher_atom_name_in,
                             const std::string &matcher_alt_conf_in) :
         reference_chain_id(reference_chain_id_in),
         matcher_chain_id(matcher_chain_id_in),
         reference_atom_name(reference_atom_name_in),
         reference_alt_conf(reference_alt_conf_in),
         matcher_atom_name(matcher_atom_name_in),
         matcher_alt_conf(matcher_alt_conf_in)
      {
         is_single_atom_match = 1;
         match_type_flag = lsq_t::ALL; // COOT_LSQ_ALL;
         to_reference_start_resno = reference_resno_in;
         to_reference_end_resno   = reference_resno_in;
         from_matcher_start_resno = matcher_resno_in;
         from_matcher_end_resno   = matcher_resno_in;
         model_number_matcher = 0;
         model_number_reference = 0;
         std::string somethinga = reference_insertion_code_in; // make orange lines above be
         std::string somethingb = matcher_insertion_code_in;   // green lines here
      }
      void set_model_number_reference(int ino) {
         model_number_reference = ino;
      }
      void set_model_number_matcher(int ino) {
         model_number_matcher = ino;
      }
      friend std::ostream&  operator<<(std::ostream&  s, const lsq_range_match_info_t &q);
   };
   std::ostream&  operator<<(std::ostream&  s, const lsq_range_match_info_t &q);

   float get_position_hash(mmdb::Manager *mol);

   bool sort_chains_util(const std::pair<mmdb::Chain *, std::string> &a,
                         const std::pair<mmdb::Chain *, std::string> &b);

   // return -1 on badness
   int get_selection_handle(mmdb::Manager *mol, const atom_spec_t &at);

   // return a selection handle.
   //
   // mode:
   // 0: all
   // 1: mainchain
   // 2: not mainchain
   // 3: not mainchain or CB
   // caller deletes the selection.
   int
   specs_to_atom_selection(const std::vector<coot::residue_spec_t> &specs,
                           mmdb::Manager *mol,
                           int atom_mask_mode);

   // Return the lsq deviation of the pt atom and (in second) the rms
   // deviation of the atoms in the plane (not the including pt of
   // course)
   //
   class lsq_plane_info_t {
      std::vector<double> abcd;
      double rms;
      clipper::Coord_orth centre_;
   public:
      lsq_plane_info_t() {}
      // can throw an exception
      explicit lsq_plane_info_t(const std::vector<clipper::Coord_orth> &v);
      // can throw an exception
      double plane_deviation(const clipper::Coord_orth &pt) const {
         if (abcd.size() == 4)
            return abcd[0]*pt.x() + abcd[1]*pt.y() + abcd[2]*pt.z() - abcd[3];
         else
            throw std::runtime_error("no plane defined");
      }
      //! project the point pt onto the lsq plane
      clipper::Coord_orth projected_point(const clipper::Coord_orth &pt) {
         const double &a = abcd[0];
         const double &b = abcd[1];
         const double &c = abcd[2];
         const double &d = abcd[3];
         double top = a * pt.x() + b * pt.y() + c * pt.z() + d;
         double bot = a * a + b * b + c * c;
         double f = top / bot;
         double x_p = pt.x() - a * f;
         double y_p = pt.y() - b * f;
         double z_p = pt.z() - c * f;
         return clipper::Coord_orth(x_p, y_p, z_p);
      }
      //! return the r.m.s. deviation of the point in the plane
      double plane_atoms_rms() const {
         return rms;
      }
      //! return the lsq plane normal
      clipper::Coord_orth normal() const {
         return clipper::Coord_orth(abcd[0], abcd[1], abcd[2]);
      }
      clipper::Coord_orth centre() const {
         return centre_;
      }
      double angle(const clipper::Coord_orth &vect) const {
         clipper::Coord_orth uv(vect.unit());
         clipper::Coord_orth abc(abcd[0], abcd[1], abcd[2]);
         clipper::Coord_orth abc_uv(abc.unit());
         double cos_theta = clipper::Coord_orth::dot(uv, abc_uv);
         double theta = acos(cos_theta) * 2 * M_PI;
         return theta;
      }
      double a() const { return abcd[0]; }
      double b() const { return abcd[1]; }
      double c() const { return abcd[2]; }
      double d() const { return abcd[3]; }
   };


   std::pair<double, double>
   lsq_plane_deviation(const std::vector<clipper::Coord_orth> &v,
                       const clipper::Coord_orth &pt);

   bool is_member_p(const std::vector<mmdb::Residue *> &v, mmdb::Residue *a);

   bool is_hydrogen_atom(mmdb::Atom *at);

   // Throw an exception if there is no consistent seg id for the
   // atoms in the given residue.
   std::string residue_atoms_segid(mmdb::Residue *residue_p);

   // Throw an exception if there is no consistent seg id for the
   // atoms in the given chain.
   std::string chain_atoms_segid(mmdb::Chain *chain_p);

   // Use the above function the get the segid and insert it into all
   // the atoms of receiver.
   bool copy_segid(mmdb::Residue *provider, mmdb::Residue *receiver);

   // residues are in increasing number order (if they are not, this
   // causes the residue selection of mmdb to fail at the misplaced
   // residue).
   //
   bool residues_in_order_p(mmdb::Chain *chain_p);

   //! get the centre of the molecule
   //!
   //! The mass of the atoms is not used.
   //!
   //! @ return success status as first element
   std::pair<bool, clipper::Coord_orth> centre_of_molecule(mmdb::Manager *mol);

   //! get the centre of the molecule, using atom masses
   //!
   //! @ return success status as first element
   std::pair<bool, clipper::Coord_orth> centre_of_molecule_using_masses(mmdb::Manager *mol);

   //! get the radius of gyration - using the centre from above
   std::pair<bool, double> radius_of_gyration(mmdb::Manager *mol);

   std::pair<bool, clipper::Coord_orth> centre_of_residues(const std::vector<mmdb::Residue *> &residues);

   // convert atoms in residue to HETATMs if residue is not of
   // standard PDB type
   int hetify_residue_atoms_as_needed(mmdb::Residue *res);
   int hetify_residues_as_needed(mmdb::Manager *mol);
   void put_amino_acid_residue_atom_in_standard_order(mmdb::Residue *residue_p);

   std::vector<mmdb::Atom *> atoms_with_zero_occupancy(mmdb::Manager *mol);

   std::vector<mmdb::Residue *> residues_with_alt_confs(mmdb::Manager *mol);

   // convert atoms in residue to HETATMs.  Return the number of HET
   // atoms.
   int hetify_residue_atoms(mmdb::Residue *res);

   // return residue specs for residues that have atoms that are
   // closer than radius Angstroems to any atom in the residue
   // specified by res_in.
   //
   std::vector<residue_spec_t> residues_near_residue(const residue_spec_t &res_in,
                                                     mmdb::Manager *mol,
                                                     float radius);

   std::vector<mmdb::Residue *> residues_near_residue(mmdb::Residue *res_ref, mmdb::Manager *mol,
                                                 float radius);

   // calling residues_near_residue for every residue in a chain is slow.
   // Let's make a map to store the results of just one selection
   //
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > residues_near_residues(const std::vector<std::pair<bool,mmdb::Residue *> > &residues_vec,
                                                                                mmdb::Manager *mol,
                                                                                float dist_crit);

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > residues_near_residues(mmdb::Manager *mol, float dist_crit);

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > residues_near_residues_for_residues(const std::map<mmdb::Residue *, std::set<mmdb::Residue *> > &all_molecule_map,
                                                                                             const std::set<mmdb::Residue *> &limit_to_these_residues_vec);

   std::pair<std::set<mmdb::Residue *>, std::set<mmdb::Residue *> > interface_residues(mmdb::Manager *mol,
                                                                                       const std::string &chain_A,
                                                                                       const std::string &chain_B,
                                                                                       float min_dist);

   std::vector<mmdb::Residue *> residues_near_position(const clipper::Coord_orth &pt,
                                                  mmdb::Manager *mol,
                                                  double radius);


   // Don't include residues that are HOH residues that are not bonded to
   // the protein (if bonding to res_ref but not protein, then reject.
   // Reject waters that are not within water_dist_max to any atom in
   // res_ref.
   //
   std::vector<mmdb::Residue *> filter_residues_by_solvent_contact(mmdb::Residue *res_ref,
                                                              mmdb::Manager *mol,
                                                              const std::vector<mmdb::Residue *> &residues,
                                                              const double &water_dist_max);

   // filter out residues that are HOH and are not in contact with
   // res_ref (i.e. have an atom closer than water_dist_max)
   //
   // std::vector<mmdb::Residue *>
   // filter_residues_by_solvent_contact_to_residue(mmdb::Residue *res_ref,
   // const std::vector<mmdb::Residue *> &residues,
   // const double &water_dist_max);
   //
   // we don't do this here - we'll do it in the viewer - were we can
   // toggle on an off the waters (if they are not rejected at an
   // early stage).


   // Return a pair, the bool of which is set if the float is sensible.
   //
   std::pair<bool,float> closest_approach(mmdb::Manager *mol,
                                          mmdb::Residue *r1, mmdb::Residue *r2);

   mmdb::Residue *nearest_residue_by_sequence(mmdb::Manager *mol,
                                         const residue_spec_t &spec);

  // create a new molecule.  rtop_frac comes from the symmetry
  // operator.  If pre_shift_abc is not of size 3 then don't apply it
  // (if it is, do so, of course).
  mmdb::Manager *mol_by_symmetry(mmdb::Manager *mol,
                                clipper::Cell cell,
                                clipper::RTop_frac rtop_frac,
                                std::vector<int> pre_shift_abc);

   std::vector<std::vector<mmdb::Chain *> > ncs_related_chains(mmdb::Manager *mol, int model_number);

   class close_residues_from_different_molecules_t {
      // Interacting Residues: Return all residues of mol1, mol2 that
      // have atoms that are closer that dist to atoms of mol2/mol1.
      //
      mmdb::Manager *combined_mol;
   public:
      close_residues_from_different_molecules_t() { combined_mol = 0; };
      std::pair<std::vector<mmdb::Residue *>, std::vector<mmdb::Residue *> >
         close_residues(mmdb::Manager *mol1, mmdb::Manager *mol2, float dist);
      void clean_up() { delete combined_mol; }
   };


   // Fiddle with mol.
   //
   // sort chains in lexographical order.  Calls PDBCleanup() and FinishStructEdit().
   void sort_chains(mmdb::Manager *mol);

   // sort residues
   // Calls PDBCleanup() and FinishStructEdit().
   void sort_residues(mmdb::Manager *mol);


   // split a molecule into "bricks" - cubic sets of atom indices that don't overlap
   // return vector needs to be multiples of 8 (8 cubes will coverer all space)
   // If the atom max radius is 3A, then the brick should have length 6A.
   // brick-id: 0,8,16,24 (etc) are done at the same time
   // then
   // brick-id: 1,9,17,25 (etc) are done at the same time
   // then
   // brick-id: 2,10,18,26 (etc) are done at the same time.
   //
   // pretty tricky - not implemented yet
   // (ints because using mmdb atom selection)
   std::vector<std::vector<int> > molecule_to_bricks(mmdb::Manager *mol, int SelectionHandle,
                                                     float max_radius);
   int get_brick_id(const clipper::Coord_orth &pt, const clipper::Coord_orth &pt_minimums,
                    int nx_grid, int ny_grid, int nz_grid,
                    float brick_length);

   int get_brick_id_inner(int x_idx, int y_idx, int z_idx,
                          int nx_grid, int ny_grid, int nz_grid);

   // pucker analysis helper class - for mark-up/rendering.
   class pucker_analysis_markup_info_t {
      public:
         clipper::Coord_orth base_ring_centre;
         clipper::Coord_orth base_ring_normal;
         clipper::Coord_orth phosphorus_position;
         clipper::Coord_orth projected_point; // of the phosphorus position onto the the plane
                                              // of the base
         double phosphate_distance_to_base_plane;
         pucker_analysis_markup_info_t() {
            base_ring_centre    = clipper::Coord_orth(0,0,0);
            base_ring_normal    = clipper::Coord_orth(0,0,0);
            phosphorus_position = clipper::Coord_orth(0,0,0);
            projected_point     = clipper::Coord_orth(0,0,0);
            phosphate_distance_to_base_plane = 0.0;
         }
   };
   // Pukka puckers?
   //
   // Throw an exception if it is not possible to generate pucker info
   //
   class pucker_analysis_info_t {
      enum PUCKERED_ATOM_T { NONE, C1_PRIME, C2_PRIME, C3_PRIME, C4_PRIME, O4_PRIME};
      PUCKERED_ATOM_T puckered_atom_;
      std::string altconf;
      void assign_base_atom_coords(mmdb::Residue *residue_p);
      mmdb::Atom *N1_or_9;
      mmdb::Atom *C1_prime;
   public:
      float plane_distortion;
      float out_of_plane_distance;
      std::vector<clipper::Coord_orth> ribose_atoms_coords;
      std::vector<clipper::Coord_orth> base_atoms_coords; // the perpendicular distance of the
                                                          // following phosphate is to the *BASE*
                                                          // atoms (stupid boy).

      // Throw an exception if it is not possible to generate pucker info
      //
      pucker_analysis_info_t(mmdb::Residue *res, std::string altconf_in);

      // Use the 3' phosphate of the following residue to calculate
      // its out of plane distance.  Decide from that if this should
      // have been 3' or 2'.  Check vs the actual puckering.
      //
      float phosphate_distance(mmdb::Residue *following_res);
      float phosphate_distance_to_base_plane(mmdb::Residue *following_res);
      pucker_analysis_markup_info_t markup_info; // fill this
      std::string puckered_atom() const;
      std::string to_json() const;
   };



   class graph_match_info_t {
   public:
      // atom_match_names: ((atom_name_wrk alt_conf_work) (atom_name_ref alt_conf_ref))
      std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > matching_atom_names;
      bool success;
      clipper::RTop_orth rtop;
      int n_match;
      double dist_score;
      graph_match_info_t() {
         n_match = 0;
         success = 0;
         dist_score = 0;
      }
      // Change the names in res_moving_names to match those in
      // res_reference as much as possible.  When there is a name
      // collision (i.e. the name maped from the res_reference is
      // already in the res_moving_names (and that is not to be
      // replace by anything)), invent a new name for the atom.
      // Use the internal matching_atom_names.
      void match_names(mmdb::Residue *res_moving_names);
   private:
      // Find a new name for name_in that is not already in the residue
      //
      std::string invent_new_name(const std::string &name_in,
                                  const std::string &ele,
                                  const std::vector<std::string> &residue_atom_name) const;
   };

   // Match on graph
   //
   // Return the orientation matrix moving res_moving to res_reference
   // and a flag letting us know that the match worked OK.  Also
   // return the vector of matching atom names, which will be used in
   // the ligand analysis (where the rtop is not applied)
   //
   // Sometimes, we don't want to do an automatic rotation/translation
   // of a test ligand onto the reference structure, we want to
   // measure the distances from exactly where they are (e.g. in
   // testing ligand placement).
   //
   // Usually though apply_rtop_flag will be 1.
   //
   graph_match_info_t graph_match(mmdb::Residue *res_moving, mmdb::Residue *res_reference, bool apply_rtop_flag,
                                  bool match_hydrogens_also);


   //
   bool mol_has_symmetry(mmdb::Manager *mol);

   // (if it has any atom that is anisotropic, return true)
   bool mol_is_anisotropic(mmdb::Manager *mol);


   // This class doesn't work with alt confs - it just finds the first
   // atom names that matches the strings in the given quad.
   //
   class position_residue_by_internal_coordinates {
      // return NULL on atom not found (res_1 and res_2 are guarenteed
      // to be non-null)
      mmdb::Atom *get_atom(mmdb::Residue *res_1, mmdb::Residue *res_2, const atom_name_quad &quad, int atom_index);
   public:
      position_residue_by_internal_coordinates(mmdb::Residue *residue_ref, mmdb::Residue *residue_moving,
                                               const atom_name_quad &quad,
                                               const double &bond_length,
                                               const double &bond_angle,
                                               const double &bond_torsion); // degrees
      // return success status (0 = fail)
      bool move_moving_residue();
   };

   // Return the 4th atom attached to at_centre that is not any of the arguments.
   // Checked by sane distance.  Return the closes atom.
   //
   // Can return null if no close atoms.
   //
   mmdb::Atom *chiral_4th_atom(mmdb::Residue *residue_p, mmdb::Atom *at_centre,
                          mmdb::Atom *at_1, mmdb::Atom *at_2, mmdb::Atom *at_3);

   namespace util {

      class stats_data {
      public:
         explicit stats_data(const std::vector<float> &d);
         explicit stats_data(const std::vector<double> &d);
         float mean;
         float sd;
         float iqr;
      };

      class qq_plot_t {
         std::vector<double> data;
         int n_bins;
         double gaussian(const double &mean, const double &sd, const double &v);
      public:
         explicit qq_plot_t(const std::vector<double> &d) : data(d) {
            n_bins = 50;
         }
         std::vector<std::pair<double, double> > qq_norm(); // sorts data.
      };

      float interquartile_range(const std::vector<float> &v);

      class peptide_torsion_angles_info_t {
      public:
         short int status;
         // angles in radians from clipper::Coord_orth::torsion
         double omega;
         double phi;
         double psi;
         std::string altconf;
      };

      class atom_spec_and_button_info_t {
      public:
         atom_spec_t as;
         std::string button_label;
         std::string callback_func;
         atom_spec_and_button_info_t(atom_spec_t as_in,
                                     const std::string &button_label_in,
                                     const std::string &callback_func_in) :
            as(as_in), button_label(button_label_in), callback_func(callback_func_in) {}
      };


      // weighted RTop_orth with deviance
      class w_rtop_orth : public clipper::RTop_orth {
      public:
         w_rtop_orth() : clipper::RTop_orth() {
            deviance = 0;
            weight = 0;
         }
         clipper::RTop_orth rtop;
         double deviance;
         double weight;
      };

      class quaternion {
      public:
         float q0, q1, q2, q3;
         quaternion(const float &q0in, const float &q1in,
                    const float &q2in, const float &q3in) {
            q0 = q0in;
            q1 = q1in;
            q2 = q2in;
            q3 = q3in;
         }
         explicit quaternion(const clipper::Mat33<double> &mat_in);
         clipper::Mat33<double> matrix() const;
         float convert_sign(const float &x, const float &y) const;
         friend std::ostream&  operator<<(std::ostream&  s, const quaternion &q);
         friend std::ofstream& operator<<(std::ofstream& s, const quaternion &q);

         void normalize();
         static bool close_float_p (const float &f1, const float &f2) { //testing func
            float d = fabsf(f1-f2);
            if (d < 0.001)
               return 1;
            else
               return 0;
         }
         // angle in radians
         quaternion rotate(double angle, const clipper::Coord_orth &vec) const;
         quaternion inverse() const;
         bool is_similar_p(const quaternion &q) {
            bool r = 0;
            if (close_float_p(q.q0, q0) &&
                close_float_p(q.q1, q1) &&
                close_float_p(q.q2, q2) &&
                close_float_p(q.q3, q3)) {
               r = 1;
            }
            return r;
         }
         static void test_quaternion(); // test yourself

         clipper::RTop_orth centroid_rtop(const std::vector<std::pair<clipper::RTop_orth,float> > &rtops);
         clipper::RTop_orth centroid_rtop(const std::vector<std::pair<clipper::RTop_orth,float> > &rtops,
                                          bool robust_filter);
         static bool deviance_sorter(const w_rtop_orth &a, const w_rtop_orth &b);
      };
      std::ostream&  operator<<(std::ostream&  s, const quaternion &q);
      std::ofstream& operator<<(std::ofstream& s, const quaternion &q);

      // chain-split the residues, dont just rely on the sequence number
      bool residues_sort_function(mmdb::Residue *first, mmdb::Residue *second);

      class chain_id_residue_vec_helper_t {
      public:
         std::vector<mmdb::Residue *> residues;
         std::string chain_id;
         void sort_residues();
         static bool residues_sort_func(mmdb::Residue *first, mmdb::Residue *second);
         bool operator<(const chain_id_residue_vec_helper_t &c) const;
      };


      double refmac_atom_radius(mmdb::Atom *at);

      std::string single_letter_to_3_letter_code(char code);

      std::string single_letter_to_3_letter_code(const std::string &slc);

      short int is_nucleotide(mmdb::Residue *r); // by refmac restraint naming
                                            // e.g. Td Ur Gr Ad etc.

      bool nucleotide_is_DNA(mmdb::Residue *r);  // test for presence of O2'

      bool residue_has_hydrogens_p(mmdb::Residue *res);

      // return 0 for no, 1 for yes, -1 for NULL residue or 0 atoms;
      int residue_has_hetatms(mmdb::Residue *residue_p);

      // Return NULL on residue not found in this molecule.
      //
      mmdb::Residue *get_residue(const std::string &chain_id, int res_no,
                                 const std::string &insertion_code,
                                 mmdb::Manager *mol);

      mmdb::Residue *get_residue_by_binary_search(const std::string &chain_id, int res_no,
                                                  const std::string &insertion_code,
                                                  mmdb::Manager *mol);

      mmdb::Residue *get_first_residue(mmdb::Manager *mol);

      mmdb::Residue *get_first_residue_in_chain(mmdb::Chain *chain_p);
      mmdb::Residue *get_last_residue_in_chain(mmdb::Chain *chain_p);
      mmdb::Residue *get_second_residue_in_chain(mmdb::Chain *chain_p);
      mmdb::Residue *get_penultimate_residue_in_chain(mmdb::Chain *chain_p);

      std::vector<mmdb::Residue *> get_hetgroups(mmdb::Manager *mol, bool include_waters=false);

      mmdb::Residue *get_biggest_hetgroup(mmdb::Manager *mol);

      // Trivial helper function.
      // Can return NULL.
      //
      // Typically called with nth is 1, 2 or 3.
      // Will return 0 if called with nth is 0.
      // i.e. if you want the first residue in mol, nth is 1.
      mmdb::Residue *get_nth_residue(int nth, mmdb::Manager *mol);

      // convenience interface to above
      mmdb::Residue *get_residue(const residue_spec_t &rs, mmdb::Manager *mol);

      // return first false on failure to find residue
      std::pair<bool, clipper::Coord_orth> get_residue_mid_point(mmdb::Manager *mol, const coot::residue_spec_t &rs);

      // get this and next residue - either can be null - both need testing
      std::pair<mmdb::Residue *, mmdb::Residue *> get_this_and_next_residues(const residue_spec_t &rs, mmdb::Manager *mol);

      // Return NULL on atom not found in this molecule
      //
      mmdb::Atom *get_atom(const atom_spec_t &spec, mmdb::Manager *mol);

      // Return NULL on atom not found in this molecule.
      // Here, if the first search files, then pad the atom name in various ways to try to find a match
      mmdb::Atom *get_atom_using_fuzzy_search(const atom_spec_t &spec, mmdb::Manager *mol);

      // Return NULL on atom not found in this residue
      //
      mmdb::Atom *get_atom(const atom_spec_t &spec, mmdb::Residue *residue_p);

      // can throw an exception if atom_p is NULL
      clipper::Coord_orth get_coords(mmdb::Atom *at);

      // Return NULL on residue not found in this molecule.
      //
      mmdb::Residue *get_following_residue(const residue_spec_t &rs,
                                      mmdb::Manager *mol);

      // Return NULL on residue not found in this molecule.
      //
      mmdb::Residue *get_previous_residue(const residue_spec_t &rs,
                                          mmdb::Manager *mol);


      std::pair<bool, clipper::Coord_orth> get_residue_centre(mmdb::Residue *res);

      std::pair<bool, clipper::Coord_orth> get_CA_position_in_residue(mmdb::Residue *residue_p);

      std::pair<bool, clipper::Coord_orth> get_CB_position_in_residue(mmdb::Residue *residue_p);

      // current view has a particular orientation of the mainchain on the screen -
      // I want to move to the residue_next (which could be the previous residue)
      // for shift-space - what rotation/translation do I need?
      //
      // @return pair.first false if we don't have a valid RTop.
      //
      std::pair<bool, clipper::RTop_orth> get_reorientation_matrix(mmdb::Residue *residue_current,
                                                                   mmdb::Residue *residue_next);


      // now in mol-utils
      // std::vector<std::string> get_residue_alt_confs(mmdb::Residue *res);


      std::vector<std::string> residue_types_in_molecule(mmdb::Manager *mol);
      // non-standard means not one of the standard protein amino acid
      // residues.
      std::vector<std::string> non_standard_residue_types_in_molecule(mmdb::Manager *mol);
      // a utility function
      std::vector<std::string> standard_residue_types();

      std::vector<std::string> PDB_standard_residue_types();

      std::vector<std::string> residue_types_in_chain(mmdb::Chain *chain_p);
      std::vector<std::string> residue_types_in_residue_vec(const std::vector<mmdb::Residue *> &residues);

      std::vector<std::string> chains_in_molecule(mmdb::Manager *mol);
      std::vector<mmdb::Residue *> residues_in_molecule(mmdb::Manager *mol);
      std::vector<mmdb::Residue *> residues_in_molecule_of_type(mmdb::Manager *mol, const std::string &residue_type);
      std::vector<mmdb::Residue *> residues_in_chain(mmdb::Manager *mol, const std::string &chain_id_in);
      std::vector<mmdb::Residue *> residues_in_chain(mmdb::Chain *chain_p);
      int number_of_residues_in_molecule(mmdb::Manager *mol);
      // Return -1 on badness
      int max_number_of_residues_in_chain(mmdb::Manager *mol);

      // get the number of residue in chain, protein first.
      std::pair<unsigned int, unsigned int> get_number_of_protein_or_nucleotides(mmdb::Chain *chain_p);

      std::vector<std::string> alt_confs_in_molecule(mmdb::Manager *mol);

      // Return NULL on no such chain:
      mmdb::Chain *chain_only_of_type(mmdb::Manager *mol, const std::string &residue_type);

      clipper::RTop_orth matrix_convert(mmdb::mat44 mat);

      // return 9999,-9999 on failure
      // (test for failure: second being less than first)
      std::pair<int, int> min_and_max_residues(mmdb::Chain *chain_p);

      // Return -1 on badness.
      //
      // So that we can calculate the length of the graph x axis - there may
      // be gaps, which is why max_number_of_residues_in_chain would fail.
      //
      int max_min_max_residue_range(mmdb::Manager *mol);

      // return first == 0 if residues not found in chain
      std::pair<bool, int> min_resno_in_chain(mmdb::Chain *chain_p);
      std::pair<bool, int> max_resno_in_chain(mmdb::Chain *chain_p);
      std::pair<bool, int> max_resno_in_molecule(mmdb::Manager *mol);
      // like the above, but don't count ligands and waters that might have
      // high residue numbers
      std::pair<bool, std::pair<int, int> > min_max_residues_in_polymer_chain(mmdb::Chain *chain_p);

      std::vector<mmdb::Chain *> chains_in_atom_selection(mmdb::Manager *mol, int model_number, const std::string &atom_selection);

      // Return -1 on badness (actually, number of chains in the last model)
      int number_of_chains(mmdb::Manager *mol);

      // The sum of all occupancies:
      float occupancy_sum(mmdb::PAtom *atoms, int n_atoms);

      float median_temperature_factor(mmdb::PPAtom atom_selection,
                                      int n_atoms,
                                      float low_cutoff,
                                      float high_cutoff,
                                      bool apply_low_cutoff,
                                      bool apply_high_cuttoff);
      float average_temperature_factor(mmdb::PPAtom atom_selection,
                                       int n_atoms,
                                       float low_cutoff,
                                       float high_cutoff,
                                       short int apply_low_cutoff,
                                       short int apply_high_cuttoff);
      float standard_deviation_temperature_factor(mmdb::PPAtom atom_selection,
                                                  int n_atoms,
                                                  float low_cutoff,
                                                  float high_cutoff,
                                                  short int apply_low_cutoff,
                                                  short int apply_high_cuttoff);

      class contact_atoms_info_t {
      public:
         enum ele_index_t { ELE_UNASSIGNED, ELE_NA, ELE_K, ELE_MG2, ELE_LI, ELE_CA2 };
         class contact_atom_t {
         public:
            double dist;
            mmdb::Atom *at;
            mmdb::mat44 mat;
            contact_atom_t(mmdb::Atom *contactor, mmdb::Atom *central_atom);
            contact_atom_t(mmdb::Atom *contactor, mmdb::Atom *central_atom, const mmdb::mat44 &m);
         };
      private:
         std::vector<contact_atom_t> contact_atoms;
         mmdb::Atom *at;
         ele_index_t metal_type;
      public:
         contact_atoms_info_t() {
            at = NULL;
            metal_type = ELE_UNASSIGNED;
         }
         contact_atoms_info_t(mmdb::Atom *at_central_in,
                              const contact_atom_t &con_at) :
            at(at_central_in) {
            contact_atoms.push_back(con_at);
            metal_type = ELE_UNASSIGNED;
         }
         unsigned int size() const {
            return contact_atoms.size();
         }
         contact_atom_t &operator[](int i) {
            return contact_atoms[i];
         }
         bool matches_atom(mmdb::Atom *at_in) {
            return at_in == at;
         }
         void add(const contact_atom_t &con_at) {
            contact_atoms.push_back(con_at);
         }
         bool has_contacts(unsigned int n_contacts, double dist_max) const {
            bool r = 0;
            if (size() >= n_contacts) {
               unsigned int n_local = 0;
               for (unsigned int i=0; i<contact_atoms.size(); i++) {
                  if (contact_atoms[i].dist <= dist_max) {
                     n_local++;
                  }
               }
               if (n_local >= n_contacts)
                  r = 1;
            }
            return r;
         }
         mmdb::Atom *central_atom() const { return at; }
         double smallest_contact_dist() const {
            if (contact_atoms.size() == 0)
               throw std::runtime_error("zero contacts");
            double d = 999999999999.9;
            for (unsigned int i=0; i<contact_atoms.size(); i++) {
               if (contact_atoms[i].dist < d)
                  d = contact_atoms[i].dist;
            }
            return d;
         }
         bool test_for_na() const;
         bool test_for_ele(ele_index_t ele_index) const;
         static std::string ele_to_string(ele_index_t ele) {
            std::string r = "Unknown";
            if (ele == ELE_UNASSIGNED)
               r = "Unassigned";
            if (ele == ELE_NA)
               r = "Sodium (Na+)";
            if (ele == ELE_K)
               r = "Potassium (K+)";
            if (ele == ELE_LI)
               r = "Lithium (Li+)";
            if (ele == ELE_MG2)
               r = "Magnesium (Mg+2)";
            if (ele == ELE_CA2)
               r = "Calcium (Ca+2)";
            return r;
         }
      };

      // This does naive symmetry expansion, if the mol has symmetry
      // then this class should be used with (typically) a copy of the
      // molecule that has been moved as close as possible to the
      // origin.
      //
      class water_coordination_t {
         std::vector<contact_atoms_info_t> atom_contacts;
         void add_contact(mmdb::Atom *atom_contactor,
                          mmdb::Atom *atom_central,
                          const mmdb::mat44 &m);

         void add_contacts(mmdb::Manager *mol,
                           mmdb::PAtom *water_selection, int n_water_atoms,
                           mmdb::PAtom *atom_selection, int n_selected_atoms,
                           mmdb::realtype min_dist, mmdb::realtype max_dist,
                           const mmdb::mat44 &my_mat);
         void sort_contacts(std::vector<contact_atoms_info_t> *v) const;
         static bool sort_contacts_func(const contact_atoms_info_t &first,
                                        const contact_atoms_info_t &second);
         void init_internal(mmdb::Manager *mol,
                            mmdb::realtype radius_limit,
                            bool do_metals_only_flag);
      public:
         water_coordination_t(mmdb::Manager *mol, mmdb::realtype radius_limit, bool do_metals_only_flag);
         water_coordination_t(mmdb::Manager *mol, mmdb::realtype radius_limit);
         water_coordination_t() {};

         // I don't know what the "get" interface to
         // water_coordination_t should be. So I'll make one up:
         //
         std::vector<contact_atoms_info_t>
         get_highly_coordinated_waters(int n_contacts,  // at least n_contacts
                                       double dist_max) const; // within dist_max

         // This checks against built-in values from the literature
         //
         std::vector<std::pair<util::contact_atoms_info_t,
                               util::contact_atoms_info_t::ele_index_t> > metals() const;

         std::vector<contact_atoms_info_t> get_contacts() const { return atom_contacts; }

         void transform_atom(int i, int j);

      };

      // a simple full copy, caller deletes.
      mmdb::Manager *copy_molecule(mmdb::Manager *mol);

      // copy cell, symm, origin and scale cards from m1 to m2 (if possible)
      // return success status.
      bool copy_cell_and_symm_headers(mmdb::Manager *m1, mmdb::Manager *m2);

      // All headers except (optionally) cell, symmetry, origin and scale.
      bool copy_headers(mmdb::Manager *m_from, mmdb::Manager *m_to, bool include_cryst);

      // adjust the atoms of residue_p
      void delete_alt_confs_except(mmdb::Residue *residue_p, const std::string &alt_conf);

      // helper function for below create_mmdbmanager function
      void transfer_links(mmdb::Manager *mol_old, mmdb::Manager *mol_new);

      // The flanking residues (if any) are in the residue selection (SelResidues).
      // The flags are not needed now we have made adjustments in the calling
      // function.
      //
      // create_mmdbmanager_from_res_selection must make adjustments
      //
      // Note: there is now a molecule-class-info version of this - perhaps
      // we should call it?  Next bug fix here: move over to the function call.
      //
      // We pass the original molecule because here we do atom index transfer.
      //
      std::pair<mmdb::Manager *, int>
      create_mmdbmanager_from_res_selection(mmdb::Manager *orig_mol,
                                            mmdb::PResidue *SelResidues,
                                            int nSelResidues,
                                            int have_flanking_residue_at_start,
                                            int have_flanking_residue_at_end,
                                            const std::string &altconf,
                                            const std::string &chain_id_1,
                                            short int residue_from_alt_conf_split_flag);

      // We don't mess with the chain ids (give as we get), but also
      // return the handle for the atom index transfer.
      std::pair<mmdb::Manager *, int> create_mmdbmanager_from_mmdbmanager(mmdb::Manager *mol);


      // we pass the mol_old so that selected header info can be transfered also
      // currently only LINKs.
      //
      // Also add the index of the reference residue (the one in molecules[imol].atom_selection.mol)
      // to the molecule that we are construction here. So that we can properly link
      // the residues in restraints_container (there we rather need to know the references indices,
      // not the indices from the fragment molecule). The label (for lookup later) is
      // "index from reference residue".
      //
      // the first of the retuned pair should be the success status - currently it is set to true
      //
      std::pair<bool, mmdb::Manager *>
      create_mmdbmanager_from_residue_vector(const std::vector<mmdb::Residue *> &res_vec,
                                             mmdb::Manager *mol_old,
                                             const std::pair<bool,std::string> &use_alt_conf = std::pair<bool, std::string>(false, ""));

      // ignore atom index transfer, return NULL on error.
      //
      // create a copy of res and add it to a newly created molecule.
      // Return the new molecule. Return NULL on error.
      //
      mmdb::Manager *create_mmdbmanager_from_residue(mmdb::Residue *res);

      mmdb::Manager *create_mmdbmanager_from_atom_selection(mmdb::Manager *orig_mol,
                                                           int SelectionHandle,
                                                           bool invert_selection_flag=0);
      // uses the following:
      mmdb::Manager *create_mmdbmanager_from_atom_selection_straight(mmdb::Manager *orig_mol,
                                                                    int SelectionHandle);
      // Beware: This destroys (inverts) the atom selection as passed.
      mmdb::Manager *create_mmdbmanager_from_inverted_atom_selection(mmdb::Manager *orig_mol,
                                                                    int SelectionHandle);

      mmdb::Manager *create_mmdbmanager_from_atom(mmdb::Atom *at);

      // a new residue for each point.  Caller deletes.
      //
      mmdb::Manager *create_mmdbmanager_from_points(const std::vector<clipper::Coord_orth> &pts, float b_factor=20.0);

      // calling function deletes
      //
      mmdb::Manager *create_mmdbmanager_from_residue_specs(const std::vector<residue_spec_t> &r1,
                                                          mmdb::Manager *mol);

      // split a NRM model into multiple models all with MODEL 1.
      std::vector<mmdb::Manager *> split_multi_model_molecule(mmdb::Manager *mol);

      void add_copy_of_atom(mmdb::Manager *mol, mmdb::Atom *atom);

      // Important for bonding in refinement.
      // This sets the internal index of the residue in the chain array (it doesn't touch atoms).
      void pdbcleanup_serial_residue_numbers(mmdb::Manager *mol);

      // return success status, 1 is good, 0 is fail.  Use clipper::Coord_orth constructor
      //
      bool add_atom(mmdb::Residue *res,
                    const std::string &atom_name_1,
                    const std::string &atom_name_2,
                    const std::string &atom_name_3,
                    const std::string &alt_conf,
                    double length,
                    double angle, // degrees
                    double torsion,
                    const std::string &new_atom_name,
                    const std::string &new_atom_ele,
                    float new_atom_occ,
                    float new_atom_b_factor); // degrees

      // copy the chain and add it to a new molecule hierarchy
      std::pair<mmdb::Chain *, mmdb::Manager *> copy_chain(mmdb::Chain *chain_p);
      // to be used with copy_atoms_from_chain_to_chain() - this needs the same residue and atoms - the
      // difference between the chains is the positions of the atoms
      void copy_atoms_from_chain_to_chain(mmdb::Chain *from_chain, mmdb::Chain *to_chain);

      // add or delete residues and atoms as needed.
      void replace_chain_contents_with_atoms_from_chain(mmdb::Chain *from_chain,
                                                        mmdb::Manager *from_mol_orig,
                                                        mmdb::Chain *to_chain,
                                                        bool do_finishstructedit);

      // utility function for above:
      mmdb::Residue* deep_copy_this_residue_add_chain(mmdb::Residue *residue,
                                                 const std::string &altconf,
                                                 bool whole_residue_flag,
                                                 bool attach_to_new_chain_flag);

      // This doesn't add the new (returned) residue to the chain of
      // passed residue.  Simple copy of residue and atoms.
      //
      mmdb::Residue *deep_copy_this_residue(mmdb::Residue *residue);

      // As above but use the alt conf flags to filter copied atoms.
      // Can return 0 if there are no atoms copied
      //
      // If use_alt_conf first is true then
      //    if use_alt_conf second is not blank
      // then copy atoms that have blank alt conf and those
      // with altconfs that match use_alt_conf.second.
      mmdb::Residue *deep_copy_this_residue(mmdb::Residue *residue,
                                            const std::pair<bool,std::string> &use_alt_conf);

      mmdb::Residue *copy_and_delete_hydrogens(mmdb::Residue *residue_in);

      // utility function for above:
      mmdb::Residue* deep_copy_this_residue_with_atom_index_and_afix_transfer(mmdb::Manager *std_mol,
                                                                         const mmdb::Residue *residue,
                                                                         const std::string &altconf,
                                                                         short int whole_residue_flag,
                                                                         int new_atom_index_udd_handle,
                                                                         int new_afix_handle);

      // return a vector of residue that are in this fragment.
      // Fragments are marked by consecutively numbered residues.  A
      // gap in the sequence numbers marks the end/beginning of a
      // fragment.
      std::vector<mmdb::PResidue> get_residues_in_fragment(mmdb::Chain *clicked_residue_chain_p,
                                                      residue_spec_t clicked_residue);

      // this can return an empty list if the residues for resno_start or resnoend  are not found
      std::vector<mmdb::Residue *> get_residues_in_range(mmdb::Manager *mol, const std::string &chain_id, int resno_start, int resno_end);

      // deleted by calling process
      std::pair<mmdb::Manager *, std::vector<residue_spec_t> >
      get_fragment_from_atom_spec(const atom_spec_t &atom_spec, mmdb::Manager *mol);

      // return true if something was removed from header info
      bool delete_residue_references_in_header_info(mmdb::Residue *residue_p, mmdb::Manager *mol);

      // transform the atoms in mol that are in moving_chain
      // it seems (for some reason) that atom::Transform(mat) now needs a (non-const)
      // argument or else clang complains
      //
      void transform_chain(mmdb::Manager *mol, mmdb::Chain *moving_chain,
                           int n_atoms, mmdb::PAtom *atoms, mmdb::mat44 &my_matt);

      void transform_chain(mmdb::Chain *moving_chain, const clipper::RTop_orth &rtop);

      // transform atoms in residue
      void transform_atoms(mmdb::Residue *res, const clipper::RTop_orth &rtop);

      // transform all the atom in mol
      void transform_mol(mmdb::Manager *mol, const clipper::RTop_orth &rtop);

      // transform the atom selection (provided by handle SelHnd) in mol
      void transform_selection(mmdb::Manager *mol, int SelHnd, const clipper::RTop_orth &rtop);


      // Rotate position round vector, return a position.
      //
      // Note that angle is in radians.
      //
      clipper::Coord_orth rotate_around_vector(const clipper::Coord_orth &direction,
                                               const clipper::Coord_orth &position,
                                               const clipper::Coord_orth &origin_shift,
                                               double angle);
      // angle in radians
      //
      void rotate_residue(mmdb::Residue *residue_p,
                          const clipper::Coord_orth &direction,
                          const clipper::Coord_orth &origin_shift,
                          double angle);

      // move the coordinates of at:
      // angle in radians
      void rotate_atom_about(const clipper::Coord_orth &direction,
                             const clipper::Coord_orth &origin_shift,
                             double angle, mmdb::Atom *at);

      // This presumes that a_residue_p and b_residue_p are valid.
      std::vector<std::pair<int, int> > pair_residue_atoms(mmdb::Residue *a_residue_p,
                                                           mmdb::Residue *b_residue_p);

      //! CBs in GLY etc
      void delete_anomalous_atoms(mmdb::Manager *mol);

      //! delete all carbohydrate
      bool delete_all_carbohydrate(mmdb::Manager *mol);

      // A useful function that was (is) in molecule_class_info_t
      //
      // If this residue has a CA, return a pointer to that atoms, if
      // not, return the first atom, if no atoms, return NULL.
      //
      mmdb::Atom* intelligent_this_residue_mmdb_atom(mmdb::Residue *res_p);

      // A function to find the previous (next) residue.  Find residue
      // by previous (next) serial number.
      // Return NULL if prevous (next) resiude not found.
      // this_residue must be part of a molecule hierarchy.
      //
      // (c.f. get_following_residue()).
      //
      // (used in pepflip)
      mmdb::Residue *previous_residue(mmdb::Residue *this_residue);
      mmdb::Residue *next_residue(mmdb::Residue *this_residue);

      // normal sequence codes, X for non-protein
      std::string three_letter_to_one_letter(const std::string &resname);
      // as above, but allow specials, currently HOH -> "~"
      std::string three_letter_to_one_letter_with_specials(const std::string &resname);


      // for sequence/sequence alignment
      //
      // Take into account the insertion code too:
      std::vector<std::pair<mmdb::Residue *, int> > sort_residues_by_seqno(mmdb::PResidue *residues,
                                                                           int nResidues);

      // Use the results of the above to give us a sequence string:
      std::string model_sequence(const std::vector<std::pair<mmdb::Residue *, int> > &sa,
                                 bool allow_ligands = true);
      bool compare_residues(const std::pair<mmdb::Residue *, int> &a,
                            const std::pair<mmdb::Residue *, int> &b);

      // extents
      std::pair<clipper::Coord_orth, clipper::Coord_orth> extents(mmdb::Manager *mol);
      std::pair<clipper::Coord_orth, clipper::Coord_orth> extents(mmdb::Manager *mol,
                                                                  int SelectionHandle);
      std::pair<clipper::Coord_orth, clipper::Coord_orth> extents(mmdb::Manager *mol,
                                                                  const std::vector<residue_spec_t> &specs);

      // pair.second = 0 for failure
      // pair.first  = 1 for success
      //
      // How to convert a standard orientation residue to res position.
      //
      // return an orientation for each alt conf that you can find.
      //
      std::map<std::string, clipper::RTop_orth> get_ori_to_this_res(mmdb::Residue *res);


      // residues with insertion codes
      std::vector<mmdb::Residue *> residues_with_insertion_codes(mmdb::Manager *mol);

      // LSQing
      //
      std::pair<short int, clipper::RTop_orth>
      get_lsq_matrix(mmdb::Manager *mol1,
                     mmdb::Manager *mol2,
                     const std::vector<lsq_range_match_info_t> &matches,
                     int every_nth,
                     bool summary_to_screen=1);
      // used by above
      // On useful return, first.length == second.length and first.length > 0.
      //
      std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> >
      get_matching_indices(mmdb::Manager *mol1,
                           mmdb::Manager *mol2,
                           const lsq_range_match_info_t &match,
                           int every_nth);

      // Return the status in the first position
      // Return the angle in radians.
      std::pair<bool, double> omega_torsion(mmdb::Residue *C_residue, mmdb::Residue *N_residue,
                                            const std::string &altconf);
      bool is_cis(const double &omega_torsion); // in radians
      // and along the same lines:
      // Return the angles in radians.
      peptide_torsion_angles_info_t peptide_torsions(mmdb::Residue *C_residue, mmdb::Residue *N_residue,
                                                     const std::string &altconf);

      // return "" on no canonical name found
      std::string canonical_base_name(const std::string &res_name_in, base_t rna_or_dna);

      // Return the RTop that matches moving to reference.  Don't move
      // moving though.
      std::pair<bool,clipper::RTop_orth> base_to_base(mmdb::Residue *reference,
                                                      mmdb::Residue *moving);

      // Return the RTop that matches moving to reference.  Don't move
      // moving.  Lite above, but add phosphate and furanose atoms.
      std::pair<bool,clipper::RTop_orth> nucleotide_to_nucleotide(mmdb::Residue *reference,
                                                                  mmdb::Residue *moving,
                                                                  bool use_old_style_naming); // Cd, Ar
                                                                                // (modern is C and DA).

      // mutate the atoms of the residue with altLoc alt_conf (this
      // should be called multiple times for residues with alt_confs).
      // Return mutated state.
      //
      int mutate(mmdb::Residue *residue_p, mmdb::Residue *std_res_unoriented, const std::string &alt_conf,
                 short int shelx_flag, float b_factor=20.0);

      int mutate(mmdb::Residue *residue_p, const std::string &new_res_name);

      // a deep copy
      mmdb::Residue *get_standard_residue_instance(const std::string &res_name);
      static std::map<std::string, mmdb::Residue *> standard_residues_map_cache;

      // given a std residue oriented over residue, make the mutation
      // to std_residue
      void mutate_internal(mmdb::Residue *residue,
                           mmdb::Residue * std_residue_oriented,
                           const std::string &alt_conf,
                           short int is_from_shelx_ins_flag,
                           float b_factor=20.);

      void mutate_base(mmdb::Residue *residue, mmdb::Residue *std_base,
                       bool use_old_style_naming,
                       bool print_match_stats=false,
                       float b_factor=20.0);
      // which calls
      std::string convert_base_name(const std::string &std_base_name, bool use_old_style_naming);

      // For use with interesting-things-gui, make the list argument
      // from a vector of atom specs.
      //
      // use the user data in the atom spec to give us the molecule number
      // and the button label
      //
      std::string interesting_things_list(const std::vector<atom_spec_t> &v);
      std::string interesting_things_list_py(const std::vector<atom_spec_t> &v);

      // error_type is e.g. "Z score", "Clash gap"
      std::string
      interesting_things_list_with_fix(const std::vector<atom_spec_and_button_info_t> &v,
                                       const std::string &error_type);
      std::string
      interesting_things_list_with_fix_py(const std::vector<atom_spec_and_button_info_t> &v,
                                          const std::string &error_type);


      // return a set of residue specifiers and a string to put on the
      // button.  Atsushi Nakagawa told me about the problem of these
      // flips and how they could be identified by B factor
      // outlierness.
      std::vector<std::pair<atom_spec_t, std::string> >
      gln_asn_b_factor_outliers(mmdb::Manager *mol);

      std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > peptide_C_N_pairs(mmdb::Chain *chain_p);
      void standardize_peptide_C_N_distances(const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &C_N_pairs);

      // this is the more modern interface of the two
      int cis_trans_conversion(mmdb::Residue *res_first, mmdb::Residue *res_second,
                               mmdb::Manager *mol, mmdb::Manager *standard_residues_mol);

      // alter the model to fix the cis peptide of the given atom
      int cis_trans_conversion(mmdb::Atom *at, bool is_N_flag, mmdb::Manager *mol, mmdb::Manager *standard_residues_mol);

      // mol_residues, trans_residues, cis_residues must be at least of length 2.
      int cis_trans_convert(std::pair<mmdb::Residue *, mmdb::Residue *> mol_residues, // for mc conversion
                            mmdb::PResidue *trans_residues,
                            mmdb::PResidue *cis_residues);

      // return the number of cis peptides in mol:
      int count_cis_peptides(mmdb::Manager *mol);

      // more info on the real cis peptides derived from atom positions:
      std::vector<cis_peptide_info_t>
      cis_peptides_info_from_coords(mmdb::Manager *mol);

      // remove wrong cis_peptides from the header records
      void remove_wrong_cis_peptides(mmdb::Manager *mol);

      // remove 1555 links between atoms that are more than dist_min
      void remove_long_links(mmdb::Manager *mol, mmdb::realtype dist_min);

      // LINKs now have distances - let's make sure that they are correct
      //
      void correct_link_distances(mmdb::Manager *mol);

      // return the number of changed links
      unsigned int change_chain_in_links(mmdb::Model *model_p,
                                         const std::string &from_chain_id,
                                         const std::string &to_chain_id);

      // move hetgroups round protein.  Find the centres of each
      // hetgroup and move it to the protein.  Waters are handled individually.
      // Fiddle with mol.
      //
      void move_hetgroups_around_protein(mmdb::Manager *mol);

      // move waters round protein, fiddle with mol.
      // return the number of moved waters.
      int move_waters_around_protein(mmdb::Manager *mol);
      // above uses this more primitive interface
      //
      std::vector<std::pair<mmdb::Atom *, clipper::Coord_orth> >
      symmetry_move_atoms(const std::vector<clipper::Coord_orth> &protein_coords,
                          const std::vector<std::pair<mmdb::Atom*, clipper::Coord_orth> > &water_atoms,
                          clipper::Cell cell,
                          clipper::Spacegroup spacegroup);


      // Throw an std::runtime_error exception on
      // not-able-to-extract-cell/symm-info.  (In such a case, we convert a
      // clipper::Message_base to a std::runtime_error).
      //
      std::pair<clipper::Cell, clipper::Spacegroup> get_cell_symm(mmdb::Manager *mol);

      // shove a cell from a clipper cell into the passed mol.
      bool set_mol_cell(mmdb::Manager *mol, clipper::Cell cell);


      // caller must check that others has some points in it.
      //
      double min_dist_to_points(const clipper::Coord_orth &pt,
                                const std::vector<clipper::Coord_orth> &others);

      // Return the fractional shift needed to translate the protein
      // as close as possible to the origin (do not apply the shift).
      //
      // Throw an exception (e.g. no cell or symmetry).
      //
      clipper::Coord_frac shift_to_origin(mmdb::Manager *mol);
      clipper::Coord_frac shift_to_origin(const std::vector<clipper::Coord_orth> &protein_coords,
                                          clipper::Cell cell,
                                          clipper::Spacegroup spacegroup);
      bool is_000_shift(const clipper::Coord_frac &cf_shift);

      //
      clipper::Mat33<double> residue_orientation(mmdb::Residue *residue_p,
                                                 const clipper::Mat33<double> &orientation_in);

      clipper::Coord_orth average_position(std::vector<clipper::Coord_orth> &pts);

      clipper::Coord_orth average_position(mmdb::Residue *r);

      // Return the median position.  Throw an exception on failure
      // (e.g no atoms).
      //
      clipper::Coord_orth median_position(mmdb::Manager *mol);
      // also throws an exception
      clipper::Coord_orth median_position(const std::vector<clipper::Coord_orth> &pts);

      void shift(mmdb::Manager *mol, clipper::Coord_orth pt);
      //
      clipper::Coord_orth
      translate_close_to_origin(const clipper::Coord_orth water_pos_pre,
                                const clipper::Cell &cell);

      void translate_close_to_origin(mmdb::Manager *mol);

      // Print secondary structure info:
      void print_secondary_structure_info(mmdb::Model *model_p);

      // return a string description of MMDB SSE values
      std::string sse_to_string(int sse);

      // Can we find CO orientation oddities like cablam does?
      std::vector<std::pair<mmdb::Residue *, double> > CO_orientations(mmdb::Manager *mol);

      void multi_parse_prosmart_log_and_gen_CO_plot();

      void
      parse_prosmart_log_and_gen_CO_plot(const std::string &prosmart_log_file_helix,
                                         const std::string &prosmart_log_file_strand,
                                         const std::string &data_points_file_name,
                                         const std::string &pdb_file_name,
                                         const std::string &chain_id);

   } // namespace util
   std::ostream&  operator<<(std::ostream&  s, const util::quaternion &q);
   std::ofstream& operator<<(std::ofstream& s, const util::quaternion &q);
   std::ostream& operator<<(std::ostream& s, const atom_spec_t &spec);
   std::ostream& operator<<(std::ostream& s, const residue_spec_t &spec);

   std::vector<clipper::RTop_orth> mtrix_info(const std::string &file_name);


} // namespace coot

#endif // HAVE_COOT_COORD_UTILS_HH
