/*
 * src/molecule-class-info.h
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2008, 2009, 2010, 2011, 2012 by the University of Oxford
 * Copyright 2013, 2014, 2015 by Medical Research Council
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


#ifndef MOLECULE_CLASS_INFO_T
#define MOLECULE_CLASS_INFO_T

#include "coords/phenix-geo.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "stereo-eye.hh"
#include <ctime>
#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include <deque>
#include <iomanip>

#include "compat/coot-sysdep.h"

enum {CONTOUR_UP, CONTOUR_DOWN};

// needs:

//#include "mmdb_manager.h"
//#include "mmdb-extras.h"
//#include "mmdb.h"

#include <epoxy/gl.h>

#include <glm/glm.hpp>

#ifndef __NVCC__
#ifdef HAVE_BOOST_THREAD // now consistent with ideal/simple-restraint.hh
#define HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#include "utils/ctpl.h"
#endif // HAVE_CXX_THREAD
#endif // __NVCC__


#ifdef USE_MOLECULES_TO_TRIANGLES
// #include <MoleculesToTriangles/CXXClasses/RendererGL.h>
#include <MoleculesToTriangles/CXXClasses/Light.h>
#include <MoleculesToTriangles/CXXClasses/Camera.h>
// #include <CXXClasses/CameraPort.h>
#include <MoleculesToTriangles/CXXClasses/SceneSetup.h>
#include <MoleculesToTriangles/CXXClasses/ColorScheme.h>
#include <MoleculesToTriangles/CXXClasses/MyMolecule.h>
#include <MoleculesToTriangles/CXXClasses/RepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.h>
#endif

#include <clipper/ccp4/ccp4_map_io.h>

#include "coords/Cartesian.hh"
#include "coords/mmdb-extras.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Bond_lines.hh"

#include "gtk-manual.h"

#include "mini-mol/mini-mol.hh"
#include "build/CalphaBuild.hh"
#include "coot-render.hh" // 20220723-PE no graphics for WebAssembly build

// #include "coot-surface/coot-surface.hh" dead now
#include "coot-utils/coot-align.hh"
#include "utils/coot-fasta.hh"
#include "coot-utils/coot-shelx.hh"
#include "utils/coot-utils.hh"
#include "utils/pir-alignment.hh"
#include "api/ghost-molecule-display.hh"
#include "drawn-ghost-molecule-display.hh"

#include "protein_db/protein_db_utils.h"

#include "select-atom-info.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-coord-extras.hh"
#include "coot-utils/xmap-stats.hh"
#include "coot-utils/atom-tree.hh"
#include "crunch-model.hh"

#include "geometry/protein-geometry.hh"

#include "ideal/simple-restraint.hh" // for extra restraints.

#include "ligand/rotamer.hh" // in ligand, for rotamer probabilty tables


#include "validation-graphs/validation-graphs.hh"  // GTK things, now part of
				 // molecule_class_info_t, they used
				 // to be part of graphics_info_t before the
                                 // array->vector change-over.

#include "ligand/dipole.hh"
#include "density-contour/density-contour-triangles.hh"

#include "coot-utils/sfcalc-genmap.hh"

#include "gl-bits.hh"

#include "pli/flev-annotations.hh" // animated ligand interactions

#include "pick.hh"

// #include "dots-representation.hh"
#include "pli/dots-representation-info.hh"
#include "named-rotamer-score.hh"

#include "c-interface-sequence.hh"
#include "map-statistics.hh"
#include "animated-ligand.hh"

#include "array-2d.hh"

namespace molecule_map_type {
   enum { TYPE_SIGMAA=0, TYPE_2FO_FC=1, TYPE_FO_FC=2, TYPE_FO_ALPHA_CALC=3,
	  TYPE_DIFF_SIGMAA=4 };
}

#include "new-centre.hh"
#include "model-view.hh"
#include "ncs.hh"
#include "atom-selection.hh"
#include "atom-attribute.hh"
#include "extra-restraints-representation.hh"
#include "additional-representation.hh"
#include "fragment-info.hh"
#include "atom-name-bits.hh"
#include "rama-rota-score.hh"
#include "api/merge-molecule-results-info-t.hh"
#include "density-results-container-t.hh"

#include "Shader.hh"

#include "updating-map-params.hh"

#include "updating-coordinates-molecule-parameters.hh"
#include "cmtz-interface.hh" // for udating molecules

#include "clipper-ccp4-map-file-wrapper.hh"
#include "model-composition-statistics.hh"

#include "coot-utils/g_triangle.hh"
#include "model-molecule-meshes.hh"

#include "fresnel-settings.hh"

glm::vec3 cartesian_to_glm(const coot::Cartesian &c);

namespace coot {

   // c.f. Bond_lines.h coords_bond_colour_t and api/bond-colour.hh - what a mess.
   //
   enum { UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2,
	  COLOUR_BY_CHAIN_BONDS=3,
	  CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
	  BONDS_NO_HYDROGENS=15,
	  CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
	  CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
	  CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS=17,
	  COLOUR_BY_MOLECULE_BONDS=8,
	  COLOUR_BY_RAINBOW_BONDS=9,
	  COLOUR_BY_B_FACTOR_BONDS=10,
	  COLOUR_BY_OCCUPANCY_BONDS=11,
	  COLOUR_BY_USER_DEFINED_COLOURS____BONDS=12,
	  COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS=13,
          COLOUR_BY_CHAIN_GOODSELL_BONDS=21
   };

   // this is now in api/coot_molecule.hh
   // enum { RESIDUE_NUMBER_UNSET = -1111};


   // a helper class - provide filenames and status for dialog widget
   //
   class backup_file_info_t {
   public:
      bool valid_status;
      int imol;
      std::string name;
      std::string description;
      std::string backup_file_name;
      timespec ctime;
      std::string get_timespec_string() const;
      backup_file_info_t() {
	 valid_status = false; // initially no backup reported
         imol = -1;
      }
      backup_file_info_t(const std::string &file_name,
                         const std::string &descr) {
         valid_status = false;
         imol = -1;
         backup_file_name = file_name;
         description = descr;
      }
   };

   std::ostream& operator<<(std::ostream &s, const ghost_molecule_display_t &ghost);


   class display_list_object_info {
   public:
      bool is_closed;
      GLuint tag_1;
      GLuint tag_2;
      int type;
      std::string atom_selection;
      int atom_selection_handle;
      bool display_it;
      bool operator==(const display_list_object_info &dloi) const {
	 return (dloi.tag_1 == tag_1);
      }
      display_list_object_info() : tag_1(0), tag_2(0) {
         atom_selection_handle = -1;
         type = 0;
	 display_it = 1;
	 is_closed = 0;
      }
      void close_yourself() {
	 is_closed = 1;
      }
   };

   class at_dist_info_t {

   public:
      int imol;
      mmdb::Atom *atom;
      float dist;
      at_dist_info_t(int imol_in, mmdb::Atom *atom_in, float dist_in) {
	 imol = imol_in;
	 atom = atom_in;
	 dist = dist_in;
      }
   };


   class goto_residue_string_info_t {
   public:
      bool res_no_is_set;
      bool chain_id_is_set;
      int res_no;
      std::string chain_id;
      goto_residue_string_info_t(const std::string &goto_residue_string, mmdb::Manager *mol);
   };
} // namespace coot


#include "generic-vertex.hh"

#include "coot-utils/cylinder.hh"

#include "Mesh.hh"
#include "LinesMesh.hh"
#include "Instanced-Markup-Mesh.hh"

bool trial_results_comparer(const std::pair<clipper::RTop_orth, float> &a,
			    const std::pair<clipper::RTop_orth, float> &b);


// Forward declaration
class graphics_info_t;

// should be arrays so that we can store lots of molecule
// informations.
//
class molecule_class_info_t {

   // This is needed because we use side by side stereo with display
   // lists.  When the maps get updated we need to generate display
   // lists - and the display lists are specific to the glarea
   // (GLContext).  So, when we compile_density_map_display_list() we
   // need to know where to store the returned diplay list index in
   // theMapContours (in first or second).
   //
   // Note SIDE_BY_SIDE_MAIN refers to the normal mono glarea.
   //
   enum { SIDE_BY_SIDE_MAIN = 0, SIDE_BY_SIDE_SECONDARY = 1};

   // we will use a pointer to this int so that we can (potentially
   // amongst other things) use it to construct gtk widget labels.
   // We need it to not go away before the whole
   // (molecule_class_info_t) object goes away.
   //
   int imol_no;

   float data_resolution_;
   float map_sigma_;
   float map_mean_;
   float map_max_;
   float map_min_;
   float sharpen_b_factor_;
   float sharpen_b_factor_kurtosis_optimised_;
   short int is_dynamically_transformed_map_flag;
   coot::ghost_molecule_display_t map_ghost_info;

   const unsigned int VAO_NOT_SET = 999999;

   // display flags:
   //
   int bonds_box_type; // public accessable via Bonds_box_type();
   short int manual_bond_colour;  // bond colour was set by hand and
				  // shouldn't be overridden by
				  // color-by-molecule-number.

   float bond_width;
   std::vector<float> bond_colour_internal;
   //
   int draw_hydrogens_flag;
   //
   short int molecule_is_all_c_alphas() const;

   std::string make_symm_atom_label_string(mmdb::PAtom atom,
					   const std::pair <symm_trans_t, Cell_Translation> &symm_trans) const;
   std::string make_atom_label_string(mmdb::PAtom atom, int brief_atom_labels_flag,
				      short int seg_ids_in_atom_labels_flag) const;

   // rebuild/save state command
   std::vector<std::string> save_state_command_strings_;
   std::string single_quote(const std::string &s) const;
   short int have_unsaved_changes_flag;
   int coot_save_index; // initially 0, incremeneted on save.

   // saving temporary files (undo)
   //
   std::string get_save_molecule_filename(const std::string &dir);
   int make_backup(const std::string &descr); // changes history details
   int make_maybe_backup_dir(const std::string &filename) const;
   bool backup_this_molecule;

   // history_index and max_history_index tell us where we are in the
   // history of modifications of this molecule.  When we "Undo"
   // history_index is decreased (by 1) (and we restore from backup).
   //
   int history_index;
   int max_history_index;
   void save_history_file_name(const std::string &file, const std::string &description);
   std::vector<coot::backup_file_info_t> history_filename_vec;
   std::string save_time_string;
   // return success status.
   bool restore_from_backup(int history_offset, const std::string &cwd);

   public: // FIXME later

   /*! \brief Make a backup for a model molecule
    *
    * @param imol the model molecule index
    * @description a description that goes along with this back point
    */
   int make_backup_checkpoint(const std::string &description);

   /*! \brief Restore molecule from backup
    * 
    * restore model @p imol to checkpoint backup @p backup_index
    *
    * @param imol the model molecule index
    * @param backup_index the backup index to restore to
    */
   int restore_to_backup_checkpoint(int backup_index);

   /*! \brief Compare current model to backup
    * 
    * @param imol the model molecule index
    * @param backup_index the backup index to restore to
    * @return a list of residue specs for residues that have
    *         at least one atom in a different place.
    *   the first says is the backup_index was valid.
    */
   std::pair<bool, std::vector<coot::residue_spec_t> > compare_current_model_to_backup(int backup_index);

   /*! \brief Get backup info
    * 
    * @param imol the model molecule index
    * @param backup_index the backup index to restore to
    * @return a Python list of the given description (str)
    *         and a timestamp (str).
    */
   coot::backup_file_info_t get_backup_info(int backup_index);

   void print_backup_history_info() const;

   private:

   // map tools
   //
   void set_initial_contour_level(); // tinker with the class data.
				     // Must be called after sigma_
				     // and is_diff_map has been set

   // fourier
   std::string fourier_f_label;
   std::string fourier_phi_label;
   std::string fourier_weight_label;   // tested vs "" (unset).

   // refmac
   short int have_sensible_refmac_params; // has public interface;
   std::string refmac_mtz_filename;
   std::string refmac_file_mtz_filename;
   std::string refmac_fobs_col;
   std::string refmac_sigfobs_col;
   std::string refmac_r_free_col;
   int refmac_r_free_flag_sensible;
   int refmac_count;

   short int have_refmac_phase_params; // has public interface;
   std::string refmac_phi_col;
   std::string refmac_fom_col;
   std::string refmac_hla_col;
   std::string refmac_hlb_col;
   std::string refmac_hlc_col;
   std::string refmac_hld_col;

   // generic

   // change asc.
   atom_selection_container_t filter_atom_selection_container_CA_sidechain_only(atom_selection_container_t asc) const;

   // skeleton colouring flag
   short int colour_skeleton_by_random;

   // Display list map contours id.  For left and right in
   // side-by-side stereo.  For mono (or hardware stereo) use just theMapContours.first()
   //
   std::pair<GLuint, GLuint> theMapContours;

   // Noble surface display list id
   GLuint theSurface;
   // If we want the surface transparent, we have to do it in
   // immediate mode. transparent_molecular_surface_flag is 0 by
   // default.
   void draw_transparent_molecular_surface(); // the function to draw surface transparently

   // difference map negative level colour relative to positive level:
   float rotate_colour_map_for_difference_map; // 240.0 default colour_map_rotation

   std::pair<float, std::vector<mmdb::Atom *> > get_clash_score(const coot::minimol::molecule &a_rotamer, bool score_hydrogen_atoms_flag, int water_interaction_modeb) const;
   std::pair<double, clipper::Coord_orth> get_minimol_pos(const coot::minimol::molecule &a_rotamer) const;

   // backup is done in the wrappers.  (This factorization needed for
   // add_alt_conf/INSERT_CHANGE_ALTCONF supprt).
   void insert_coords_internal(const atom_selection_container_t &asc);

   // In this instance, we don't want to install a whole residue, we want
   // to install atoms in this residue (alt conf B) into a a atom_sel mol
   // residue that contains (say) "" and "A".
   //
   // pass shelx_occ_fvar_number -1 if we don't have it.
   void insert_coords_atoms_into_residue_internal(const atom_selection_container_t &asc,
						  int shelx_occ_fvar_number);

   bool input_molecule_was_in_mmcif;
   // public acces to this is below

   void unalt_conf_residue_atoms(mmdb::Residue *residue_p);

   // return status and a chain id [status = 0 when there are 26 chains...]
   //
   std::pair<bool, std::string> unused_chain_id() const;
   coot::atom_name_bits_t atom_name_to_atom_name_plus_res_name() const;

   int set_coot_save_index(const std::string &s);

   // move this where they belong (from globjects)
   void set_symm_bond_colour_mol(int i);
   void set_symm_bond_colour_mol_and_symop(int icol, int isymop);
   void set_symm_bond_colour_mol_rotate_colour_map(int icol, int isymop);
   void rotate_rgb_in_place(float *rgb, const float &amount) const;
   void convert_rgb_to_hsv_in_place(const float *rgb, float *hsv) const;
   void convert_hsv_to_rgb_in_place(const float* hsv, float *rgb) const;

   float combine_colour (float v, int col_part_index);

   std::vector<coot::dipole> dipoles;

   // make fphidata lie within the resolution limits of reso.  Do we
   // need a cell to do this?
   //
   void filter_by_resolution(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata,
			     const float &reso_low,
			     const float &reso_high) const;
   // Retard the phases for use with anomalous data.
   void shift_90_anomalous_phases(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata) const;


   // merge molecules helper function
   std::vector<std::string> map_chains_to_new_chains(const std::vector<std::string> &adding_model_chains,
						     const std::vector<std::string> &this_model_chains) const;
   // that's too complicated for try_add_by_consolidation(), we just want this:
   std::string suggest_new_chain_id(const std::string &current_chain_id) const;

   // returned the copied residue (possibly can return NULL on failure).
   mmdb::Residue* copy_and_add_residue_to_chain(mmdb::Chain *this_model_chain,
						mmdb::Residue *add_model_residue,
						bool new_res_no_by_hundreds=false);
   void copy_and_add_chain_residues_to_chain(mmdb::Chain *new_chain, mmdb::Chain *this_molecule_chain);


   short int ligand_flip_number;

   // NCS ghost molecules:
   //
   std::vector<drawn_ghost_molecule_display_t> ncs_ghosts;
   void update_ghosts();
   bool show_ghosts_flag; // i.e. draw_it_for_ncs_ghosts
   float ghost_bond_width;
   bool ncs_ghost_chain_is_a_target_chain_p(const std::string &chain_id) const;
   // throw an exception when the matrix is not defined (e.g. no atoms).
   coot::ncs_matrix_info_t find_ncs_matrix(int SelHandle1, int SelHandle2) const;
   short int ncs_ghosts_have_rtops_flag;
   // have to take into account the potential built/non-built offsets:
   bool ncs_chains_match_p(const std::vector<std::pair<std::string, int> > &v1,
			   const std::vector<std::pair<std::string, int> > &v2,
			   float exact_homology_level,
			   bool allow_offset_flag) const;
   bool ncs_chains_match_with_offset_p(const std::vector<std::pair<std::string, int> > &v1,
				       const std::vector<std::pair<std::string, int> > &v2,
				       float exact_homology_level) const;
   bool last_ghost_matching_target_chain_id_p(int i_match,
					      const std::vector<drawn_ghost_molecule_display_t> &ncs_ghosts) const;
   void delete_ghost_selections();
   // void debug_ghosts() const; public

   std::vector<coot::ghost_molecule_display_t> strict_ncs_info;
   std::vector<coot::coot_mat44> strict_ncs_matrices;

   std::vector <coot::atom_spec_t>
   find_water_baddies_AND(float b_factor_lim, const clipper::Xmap<float> &xmap_in,
			  float map_sigma,
			  float outlier_sigma_level,
			  float min_dist, float max_dist,
			  short int part_occ_contact_flag,
			  short int zero_occ_flag);

   std::vector <coot::atom_spec_t>
   find_water_baddies_OR(float b_factor_lim, const clipper::Xmap<float> &xmap_in,
			 float map_sigma,
			 float outlier_sigma_level,
			 float min_dist, float max_dist,
			 short int part_occ_contact_flag,
			 short int zero_occ_flag);

   // alignment
   coot::chain_mutation_info_container_t
   align_on_chain(const std::string &chain_id,
		  mmdb::PResidue *SelResidues, int nSelResidues,
		  const std::string &target,
		  mmdb::realtype wgap,
		  mmdb::realtype wspace,
		  bool is_nucleic_acid_flag = false,
		  bool console_output = true) const;

   std::string
   make_model_string_for_alignment(mmdb::PResidue *SelResidues,
				   int nSelResidues) const;
   // renumber_reidues starting at 1 and removing insertion codes
   // (no backup).  For use in alignment (maybe other places)
   void simplify_numbering_internal(mmdb::Chain *chain_p);

   std::string output_alignment_in_blocks(const std::string &aligned,
					  const std::string &target,
					  const std::string &matches) const;

   // String munging helper function (for reading mtz files).
   // Return a pair.first string of length 0 on error to construct dataname(s).
   std::pair<std::string, std::string> make_import_datanames(const std::string &fcol,
							     const std::string &phi_col,
							     const std::string &weight_col,
							     int use_weights) const;

   // change chain id internal function
   std::pair<int, std::string>
   change_chain_id_with_residue_range(const std::string &from_chain_id,
				      const std::string &to_chain_id,
				      int start_resno,
				      int end_resno);

   void change_chain_id_with_residue_range_helper_insert_or_add(mmdb::Chain *to_chain_p, mmdb::Residue *new_residue);

   // private nomenclature fixing (maybe) function.
   int test_and_fix_PHE_TYR_nomenclature_errors(mmdb::Residue *res);

   // shelx stuff
   bool is_from_shelx_ins_flag;
   coot::ShelxIns shelxins;

   // for NCS copying of A chain onto the other chains.
   int copy_chain(mmdb::Chain *from_chain, mmdb::Chain *to_chain,
		  clipper::RTop_orth a_to_b_transform);

   std::vector<pli::dots_representation_info_t> dots;
   coot::colour_t dots_colour;
   bool dots_colour_set;

   // save the data used for the fourier, so that we can use it to
   // sharpen the map:
   // uncommenting the following line causes a crash in the multi-molecule
   // (expand molecule space) test.
   bool original_fphis_filled;
   bool original_fobs_sigfobs_filled;
   bool original_fobs_sigfobs_fill_tried_and_failed;
   clipper::HKL_data< clipper::datatypes::F_phi<float> >  *original_fphis_p;
   clipper::HKL_data< clipper::datatypes::F_sigF<float> > *original_fobs_sigfobs_p;
   clipper::HKL_data< clipper::data32::Flag> *original_r_free_flags_p;


   // is the CCP4 map a EM map? (this is so that we can fill the
   // NXmap, not the xmap)
   //
   // bool is_em_map(const clipper::CCP4MAPfile &file) const;
   bool set_is_em_map(const clipper_map_file_wrapper &file, const std::string &file_name);

   // for quads/triangle strip for the bond representation (rather
   // than gl_lines).
   coot::Cartesian get_vector_pependicular_to_screen_z(const coot::Cartesian &front,
						       const coot::Cartesian &back,
						       const coot::Cartesian &bond_dir,
						       float zoom,
						       float p_bond_width) const;

   // remove TER record from residue
   bool residue_has_TER_atom(mmdb::Residue *res_p) const;
   void remove_TER_internal(mmdb::Residue *res_p);
   void remove_TER_on_last_residue(mmdb::Chain *chain_p);
   void remove_TER_on_residue_if_last_residue(mmdb::Chain *chain_p, mmdb::Residue *residue_p);

   // Return a new orienation, to be used to set the view orientation/quaternion.
   //
   std::pair<bool, clipper::RTop_orth>
   apply_ncs_to_view_orientation_forward(const clipper::Mat33<double> &current_view_mat,
					 const clipper::Coord_orth &current_position,
					 const std::string &current_chain,
					 const std::string &next_ncs_chain) const;
   std::pair<bool, clipper::RTop_orth>
   apply_ncs_to_view_orientation_backward(const clipper::Mat33<double> &current_view_mat,
					  const clipper::Coord_orth &current_position,
					  const std::string &current_chain,
					  const std::string &next_ncs_chain) const;

   // single model view
   int single_model_view_current_model_number;

   // return a non-null string on a problem
   std::string jed_flip_internal(coot::atom_tree_t &tree,
				 const std::vector<coot::dict_torsion_restraint_t> &interesting_torsions,
				 const std::string &atom_name,
				 int clicked_atom_idx,
				 bool invert_selection);

   // return a non-null string on a problem
   std::string jed_flip_internal(coot::atom_tree_t &tree,
				 const coot::dict_torsion_restraint_t &torsion,
				 const std::string &atom_name,
				 int clicked_atom_idx,
				 bool invert_selection);

   void asn_hydrogen_position_swap(std::vector<std::pair<bool, mmdb::Residue *> > residues);

   // ----------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------
public:        //                      public
   // ----------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------

   // we should dump this constructor in the skip then.
   // it does not set imol_no;
   molecule_class_info_t() :
      map_as_mesh(Mesh("map-as-mesh")),
      map_as_mesh_gl_lines_version(Mesh("map-as-mesh-gl-lines")),
      // mesh_for_symmetry_atoms(Mesh("mesh-for-symmetry-atoms")),
      molecule_as_mesh_rama_balls(Mesh("molecule_as_mesh_rama_balls")),
      molecule_as_mesh_rota_dodecs(Mesh("molecule_as_mesh_rota_dodecs"))
   {
      setup_internal();
      // give it a pointer to something at least *vagely* sensible
      // (better than pointing at -129345453).
      imol_no = -1;
      //
   }

   // See graphics_info_t::initialize_graphics_molecules()
   //
   molecule_class_info_t(int i) :
      map_as_mesh(Mesh("map-as-mesh")),
      map_as_mesh_gl_lines_version(Mesh("map-as-mesh-gl-lines")),
      // mesh_for_symmetry_atoms(Mesh("mesh-for-symmetry-atoms")),
      molecule_as_mesh_rama_balls(Mesh("molecule_as_mesh_rama_balls")),
      molecule_as_mesh_rota_dodecs(Mesh("molecule_as_mesh_rota_dodecs"))
   {
      setup_internal();
      imol_no = i;
   }

   // Why did I add this?  It causes a bug in "expanding molecule space"
//    ~molecule_class_info_t() {
//       delete imol_no_ptr;
//       bonds_box.clear_up();
//       symmetry_bonds_box.clear_up();

//       // also the atom labels

//    }

   ~molecule_class_info_t() {

      // using this destructor causes a redraw, it seems.
      //
      // Note that the constructor used in expand_molecule_space()
      // new_molecules[i] = molecules[i];
      // does a shallow copy of the pointers.  So we can't delete them here.
      //

      // give back the memory from the map, so that we don't get
      // clipper leak message?
      draw_it = 0;
      draw_it_for_extra_restraints = false;
      draw_it_for_parallel_plane_restraints = false;
      draw_it_for_map = 0;  // don't display this thing on a redraw!
      draw_it_for_map_standard_lines = 0;
   }

   void setup_internal();

   int handle_read_draw_molecule(int imol_no_in,
				 std::string filename,
				 std::string cwd,
				 coot::protein_geometry *geom_p,
				 short int recentre_rotation_centre,
				 short int is_undo_or_redo,
				 bool allow_duplseqnum,
				 bool convert_to_v2_atom_names_flag,
				 float bond_width_in,
				 int bonds_box_type,
				 bool warn_about_missing_symmetry_flag);

   int update_molecule(std::string file_name, std::string cwd);

   void label_symmetry_atom(int i);

   // used for raster3d (where we need to know the position of the label)
   std::pair<std::string, clipper::Coord_orth>
   make_atom_label_string(unsigned int ith_labelled_atom,
			  int brief_atom_labels_flag,
			  short int seg_ids_in_atom_labels_flag) const;

   void draw_atom_label(int atom_index,
                        int brief_atom_labels_flag,
                        short int seg_ids_in_atom_labels_flag,
                        const glm::vec4 &atom_label_colour,
                        stereo_eye_t eye,
                        const glm::mat4 &mvp,
                        const glm::mat4 &view_rotation);

   void draw_symm_atom_label(int atom_index, const std::pair <symm_trans_t, Cell_Translation> &st,
                             const glm::vec4 &atom_label_colour,
                             const glm::mat4 &mvp,
                             const glm::mat4 &view_rotation);

   // don't count the mainchain of the peptide-linked residues
   //
   std::vector<mmdb::Atom *> closest_atoms_in_neighbour_residues(mmdb::Residue *residue_p, float radius) const;

   void debug_selection() const;
   void debug(bool debug_atoms_also_flag=false) const;

   void set_bond_colour_by_mol_no(int icolour,
				  bool against_a_dark_background);  // not const because
					   	                    // we also set
						                    // bond_colour_internal.
   void set_bond_colour_for_goodsell_mode(int icol, bool against_a_dark_background);

   // return the colour, don't call glColor3f();
   coot::colour_t get_bond_colour_basic(int colour_index, bool against_a_dark_background) const;
   // return the colour, don't call glColor3f();
   // make this static? so that get_glm_colour_func() works from Mesh::make_from_graphical_bonds()?
   coot::colour_t get_bond_colour_by_mol_no(int colour_index, bool against_a_dark_background) const;

   // 20220214-PE modern colour
   glm::vec4 get_bond_colour_by_colour_wheel_position(int i, int bonds_box_type) const;
   void set_bond_colour_by_colour_wheel_position(int i, int bonds_box_type);
   bool use_bespoke_grey_colour_for_carbon_atoms;
   coot::colour_t bespoke_carbon_atoms_colour;
   void set_use_bespoke_carbon_atom_colour(bool state) {
      use_bespoke_grey_colour_for_carbon_atoms = state;
      make_bonds_type_checked("set_use_bespoke_carbon_atom_colour");
   }
   void set_bespoke_carbon_atom_colour(const coot::colour_t &col) {
      bespoke_carbon_atoms_colour = col;
      make_bonds_type_checked("set_bespoke_carbon_atom_colour");
   }

   std::string name_;
   std::string get_name() const { return name_; }

   int MoleculeNumber() const { return imol_no; }

   std::string save_mtz_file_name;
   std::string save_f_col;
   std::string save_phi_col;
   std::string save_weight_col;
   int save_use_weights;
   int save_is_anomalous_map_flag;
   int save_is_diff_map_flag;
   float save_high_reso_limit;
   float save_low_reso_limit;
   int save_use_reso_limits;

   void map_fill_from_mtz(std::string mtz_file_name,
			  std::string cwd,
			  std::string f_col,
			  std::string phi_col,
			  std::string weight_col,
			  int use_weights,
			  int is_diff_map,
			  float map_sampling_rate,
                          bool updating_existing_map_flag=false);

   void map_fill_from_mtz(const coot::mtz_to_map_info_t &mmi, const std::string &wcd, float sampling_rate);

   void map_fill_from_mtz_with_reso_limits(std::string mtz_file_name,
					   std::string cwd,
					   std::string f_col,
					   std::string phi_col,
					   std::string weight_col,
					   int use_weights,
					   short int is_anomalous_flag,
					   int is_diff_map,
					   short int use_reso_flag,
					   float low_reso_limit,
					   float high_reso_limit,
					   float map_sampling_rate,
                                           bool updating_existing_map_flag=false);

   // return succes status, if mtz file is broken or empty, or
   // non-existant, return 0.
   //
   bool map_fill_from_cns_hkl(std::string cns_file_name,
			      std::string f_col,
			      int is_diff_map,
			      float map_sampling_rate);

   bool make_patterson(std::string mtz_file_name,
		       std::string f_col,
		       std::string sigf_col,
		       float map_sampling_rate);

   bool make_patterson_using_intensities(std::string mtz_file_name,
					 std::string i_col,
					 std::string sigi_col,
					 float map_sampling_rate);

   // return -1 on not found
   // 20231025-PE this is now public because it is used in
   //             on_model_toolbar_edit_chi_angles_button_clicked().
   //             That can be reworked if needed to use an api
   //             that doesn't take an atom pointer.
   int get_atom_index(mmdb::Atom *atom) {
     int idx = -1;
     if (has_model()) {
       int ic = -1;
       if (atom->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
	 idx = ic;
       }
     }
     return idx;
   }

   atom_selection_container_t atom_sel;

   // Shall we draw anything for this molecule?
   //
   int draw_it; // used by Molecule Display control, toggled using toggle fuctions.
   bool draw_model_molecule_as_lines; // default false
   void set_draw_model_molecule_as_lines(bool state); // redo the bonding if state is different
   bool draw_it_for_map;
   bool draw_it_for_map_standard_lines; // was draw_it_for_map
   int pickable_atom_selection;  // ditto (toggling).

   std::string show_spacegroup() const;

   void set_mol_triangles_is_displayed(int state);

   void set_mol_is_displayed(int state) {
      if (atom_sel.n_selected_atoms > 0) {
	 draw_it = state;
      }
      set_mol_triangles_is_displayed(state);
   }

   int get_mol_is_displayed() const { return draw_it; }

   void set_mol_is_active(int state) {
      if (atom_sel.n_selected_atoms > 0)
	 pickable_atom_selection = state;
      else
	 pickable_atom_selection = 0;
   }

   void set_map_is_displayed(int state); // 20250216-PE moved out of header, to handle
                                         // expired map contours
   bool get_map_is_displayed() const { return draw_it_for_map; }

   void set_map_is_displayed_as_standard_lines(short int state) {
      draw_it_for_map_standard_lines = state;
   }

   void do_solid_surface_for_density(short int on_off_flag);

   int atom_selection_is_pickable() const {
      return pickable_atom_selection && (atom_sel.n_selected_atoms > 0);
   }

   bool show_symmetry; // now we have individual control of the
	               // symmetry display (as well as a master
		       // control).

   void set_show_symmetry(short int istate) {
      show_symmetry = istate;
   }

   // pass 0 or 1
   void set_display_extra_restraints(int state) {
      draw_it_for_extra_restraints = state;
      if (state)
	 update_extra_restraints_representation();
   }

   // pass 0 or 1
   void set_display_parallel_plane_restraints(int state) {
      draw_it_for_parallel_plane_restraints = state;
   }

   void delete_extra_restraints_for_residue(const coot::residue_spec_t &rs);
   void delete_extra_restraints_worse_than(const double &n_sigma);

   // Unit Cell (should be one for each molecule)
   //
   void set_show_unit_cell(bool state);
   bool show_unit_cell_flag;
   bool have_unit_cell;
   void set_have_unit_cell_flag_maybe(bool warn_about_missing_symmetry_flag);

   void draw_molecule(short int do_zero_occ_spots,
		      bool against_a_dark_background,
		      bool show_cis_peptide_markups);
   void zero_occupancy_spots() const;
   void deuterium_spots() const;
   void set_occupancy_residue_range(const std::string &chain_id, int ires1, int ires2, float occ_val);
   void draw_cis_peptide_markups() const;
   void draw_bad_CA_CA_dist_spots() const;


   void set_b_factor_residue_range(const std::string &chain_id, int ires1, int ires2, float b_val);
   void set_b_factor_atom_selection(const atom_selection_container_t &asc, float b_val, bool moving_atoms);
   void set_b_factor_residues(const std::vector<std::pair<coot::residue_spec_t, double> > &rbs); // all atoms of specified
   void set_b_factor_residue(coot::residue_spec_t spec, float bf);
   void change_b_factors_of_residue_by(coot::residue_spec_t spec, float bf);


   std::vector<coot::atom_spec_t> fixed_atom_specs;
   std::vector<coot::Cartesian>   fixed_atom_positions; // updated on make_bonds_type_checked()
   void update_fixed_atom_positions();
   void update_additional_representations(const gl_context_info_t &gl_info, const coot::protein_geometry *geom);
   void update_mols_in_additional_representations(); //uses atom_sel.mol
   void draw_fixed_atom_positions() const;
   void clear_all_fixed_atoms();
   std::vector<coot::atom_spec_t> get_fixed_atoms() const;

   // we could combine these 2, because a bonds_box contains symmetry
   // data members, but we do not, because update_symmetry is more
   // elegant without it.
   //
   // Note that display_bonds now used symmetry_bonds_box.
   //
   // Tripping up here?  Need Bond_lines.h
   graphical_bonds_container bonds_box;
   std::vector<coot::additional_representations_t> add_reps;
   void remove_display_list_object_with_handle(int handle);

   std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > > symmetry_bonds_box;
   std::vector<std::pair<graphical_bonds_container, symm_trans_t> > strict_ncs_bonds_box;
   void clear_up_all_symmetry() {
      for (unsigned int i=0; i<symmetry_bonds_box.size(); i++)
	 symmetry_bonds_box[i].first.clear_up();
      for (unsigned int i=0; i<strict_ncs_bonds_box.size(); i++)
	 strict_ncs_bonds_box[i].first.clear_up();
   }


   // Return a copy of the pointer to the chain (only).  Return NULL
   // on chain with given chain ID not found.
   //
   mmdb::Chain *get_chain(const std::string &chain_id) const;

   // Return a copy of the pointer (only).  Return NULL on residue not
   // found.
   //
   mmdb::Residue *get_residue(const std::string &chain_id,
                              int reso,
                              const std::string &insertion_code) const;

   // Return a copy of the pointer (only).  Return NULL on residue not
   // found.
   //
   mmdb::Residue *get_residue(const coot::residue_spec_t &rs) const;

   std::string get_residue_name(const coot::residue_spec_t &rs) const;

   // Return a copy of the pointer (only) of the residue following
   // that of the given spec.  Return NULL on residue not found.
   //
   mmdb::Residue *get_following_residue(const coot::residue_spec_t &rs) const;

   // Useful when we know that the molecule is just one residue
   mmdb::Residue *get_first_residue();

   // A function of molecule-class that returns the position
   // (if posible) of an atom given a string e.g.
   // "A" -> nearest atom in the "A" chain
   // "a" -> nearest atom in the "a" chain, or failing that the "A" chain.
   // "43" residue 43 of the chain we are looking at.
   //"21b" -> residue 21 of the "b" chain, otherwise residue 21 of the "B"
   //         chain
   //
   // For keyboarding go-to atom.
   //
   // Return NULL on not-able-to-find-atom.
   //
   mmdb::Atom *get_atom(const std::string &go_to_residue_string,
		   const coot::atom_spec_t &active_atom_spec,
		   const coot::Cartesian &pt) const;

   mmdb::Atom *get_atom(const coot::atom_spec_t &atom_spec) const;

   mmdb::Atom *get_atom(int idx) const;
   mmdb::Atom *get_atom(const pick_info &pi) const;

   bool have_atom_close_to_position(const coot::Cartesian &pos) const;

   // return the maximum residue number in the chain. first of false means failure to do so.
   //
   std::pair<bool,int> max_res_no_in_chain(mmdb::Chain *chain_p) const;

   std::pair<bool,int> max_res_no_in_chain(const std::string &chain_id) const;
   std::pair<bool,int> min_res_no_in_chain(const std::string &chain_id) const;

   void set_draw_hydrogens_state(int i) {
      if (draw_hydrogens_flag != i) {
	      draw_hydrogens_flag = i;
	      make_bonds_type_checked();
	      update_symmetry();
      }
   }

   int draw_hydrogens() const {
      return draw_hydrogens_flag;
   }

   void makebonds(const coot::protein_geometry *geom_p, const std::set<int> &no_bonds_to_these_atom_indices);
   void makebonds(float max_dist, const coot::protein_geometry *geom_p); // maximum distance for bond (search)
   void makebonds(float min_dist, float max_dist, const coot::protein_geometry *geom_p);
   void make_ca_bonds(float min_dist, float max_dist);
   void make_ca_bonds(float min_dist, float max_dist, const std::set<int> &no_bonds_to_these_atom_indices);
   void make_ca_bonds();
   void make_ca_plus_ligands_bonds(coot::protein_geometry *pg);
   void make_ca_plus_ligands_and_sidechains_bonds(coot::protein_geometry *pg);
   void make_colour_by_chain_bonds(bool rebonding_is_needed); // simple/usual interfce to below function
   void make_colour_by_chain_bonds(const std::set<int> &no_bonds_to_these_atoms, bool c_only_flag, bool goodsell_mode, bool rebonding_is_needed);
   void make_colour_by_ncs_related_chains(bool goodsell_mode); // presume rebonding *is* needed.
   void make_colour_by_molecule_bonds(bool rebonding_is_needed);
   void bonds_no_waters_representation();
   void bonds_sec_struct_representation();
   void ca_plus_ligands_sec_struct_representation(coot::protein_geometry *pg);
   void ca_plus_ligands_rainbow_representation(coot::protein_geometry *pg);
   void ca_representation(bool rebonding_is_needed);
   void ca_plus_ligands_representation(coot::protein_geometry *pg, bool rebonding_is_needed);
   void ca_plus_ligands_and_sidechains_representation(coot::protein_geometry *pg);
   void b_factor_representation();
   void b_factor_representation_as_cas();
   void occupancy_representation();
   void user_defined_colours_representation(coot::protein_geometry *geom_p,
                                            bool all_atoms_mode,
                                            bool draw_missing_loops_flag);

   void alt_conf_view_next_alt_conf(const std::string &current_alt_conf);


   // This doesn't catch the case when__builtin_FUNCTION exists but __has_builtin does not
   // (as seems to be the case in g++ 9.2.1)

#if defined __has_builtin
#if __has_builtin (__builtin_FUNCTION)
   void make_bonds_type_checked(const char *s = __builtin_FUNCTION());
   void make_bonds_type_checked(const std::set<int> &no_bonds_to_these_atom_indices, const char *s = __builtin_FUNCTION());
   void make_glsl_bonds_type_checked(const char *s = __builtin_FUNCTION());
#else
   void make_bonds_type_checked(const char *s = 0);
   void make_bonds_type_checked(const std::set<int> &no_bonds_to_these_atom_indices, const char *s =0);
   void make_glsl_bonds_type_checked(const char *s = 0);
#endif
#else // repeat above
   void make_bonds_type_checked(const char *s = 0);
   void make_bonds_type_checked(const std::set<int> &no_bonds_to_these_atom_indices, const char *s =0);
   void make_glsl_bonds_type_checked(const char *s = 0);
#endif

   void add_to_non_drawn_bonds(const std::string &cid);
   // this clears the old no_bonds_to_these_atom_indices set and replaces it with a new one - and regens bonds.
   void set_new_non_drawn_bonds(const std::string &cid);
   std::set<int> no_bonds_to_these_atom_indices;
   void clear_non_drawn_bonds(bool regen_bonds);

   float atom_radius_scale_factor; // 3 is quite nice, 1 by default.
   void set_atom_radius_scale_factor(float sf); // regenerate

   // atom labels and symmetry atom labels, that is.
   //
   void draw_atom_labels(int brief_atom_labels_flag,
                         short int seg_ids_in_atom_labels_flag,
                         const glm::vec4 &atom_label_colour,
                         stereo_eye_t eye,
                         const glm::mat4 &mvp,
                         const glm::mat4 &view_rotation);

   //
   void update_molecule_after_additions(); // cleanup, new
					   // atom_selection, sets
					   // unsaved changes flag,
					   // makes bonds.

   void update_symmetry();
   void make_glsl_symmetry_bonds();
   void update_strict_ncs_symmetry(const coot::Cartesian &centre_point,
				   const molecule_extents_t &extents); // in m-c-i-ncs.cc
   void old_draw_anisotropic_atoms(); // old OpenGL function
   bool show_atoms_as_aniso_flag;
   bool show_aniso_atoms_as_ortep_flag;
   void set_show_atoms_as_aniso(bool state) {
      if (state != show_atoms_as_aniso_flag) {
         show_atoms_as_aniso_flag = state;
         make_bonds_type_checked("set_show_atoms_as_aniso()");
      }
   }
   void set_show_aniso_atoms_as_ortep(bool state) {
      if (state)
         show_atoms_as_aniso_flag = true;
      if (state != show_aniso_atoms_as_ortep_flag) {
         show_aniso_atoms_as_ortep_flag = state;
         make_bonds_type_checked("set_show_aniso_atoms_as_ortep()");
      }
   }

   // void draw_coord_unit_cell(const coot::colour_holder &cell_colour);
   // void draw_map_unit_cell(const coot::colour_holder &cell_colour);
   // void draw_unit_cell_internal(float rsc[8][3]);

   LinesMesh lines_mesh_for_cell;
   void setup_unit_cell();
   void draw_unit_cell(Shader *shader_p, const glm::mat4 &mvp);

   void draw_dots(Shader *shader_p,
                  const glm::mat4 &mvp,
                  const glm::mat4 &view_rotation_matrix,
                  const std::map<unsigned int, lights_info_t> &lights,
                  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                  const glm::vec4 &background_colour,
                  bool do_depth_fog);

   // return the status of whether or not the dots were cleared.
   bool clear_dots(int dots_handle);
   // clear the first open dots object with the given name.
   // return the status of whether or not the dots were cleared.
   bool clear_dots(const std::string &dots_object_name);
   int  make_dots(const std::string &atom_selection_str,
		  const std::string &dots_name,
		  float dot_density, float atom_radius_scale);
   int n_dots_sets() const {  // return the number of sets of dots.
      return dots.size();
   }
   void unset_dots_colour() {
      dots_colour_set = false; // back to default
   }

   void initialize_map_things_on_read_molecule(std::string name,
					       bool is_diff_map,
					       bool is_anomalous_map,
					       bool swap_difference_map_colours);
   void initialize_coordinate_things_on_read_molecule(std::string name);
   void initialize_coordinate_things_on_read_molecule_internal(std::string name,
							       short int is_undo_or_redo);
   void install_model(int imol_no_in,
		      atom_selection_container_t asc,
		      const coot::protein_geometry *geom_p,
		      const std::string &mol_name,
		      short int display_in_display_control_widget_status,
		      bool is_from_shelx_ins=false,
		      bool warn_about_missing_symmetry_flag=true);

   void install_model(int imol_no_in,
		      mmdb::Manager *mol,
		      const coot::protein_geometry *geom_p,
		      const std::string &mol_name,
		      short int display_in_display_control_widget_status,
		      bool is_from_shelx_ins=false,
		      bool warn_about_missing_symmetry_flag=true);

   void install_model_with_ghosts(int imol_no_in,
				  atom_selection_container_t asc,
				  const coot::protein_geometry *geom_p,
				  const std::string &mol_name,
				  short int display_in_display_control_widget_status,
				  const std::vector<coot::ghost_molecule_display_t> &ncs_ghosts_in,
				  bool is_from_shelx_ins=false,
				  bool warn_about_missing_symmetry_flag=true,
				  bool generate_ghost_info=true);

   int copy_residue_range(mmdb::Chain *from_chain, mmdb::Chain *to_chain,
			  int residue_range_1,
			  int residue_range_2,
			  clipper::RTop_orth a_to_b_transform);

   // const coot::CartesianPair* draw_vectors;
   // int n_draw_vectors;
   //
   // now coot makes many draw_vectors by sending off a "set" - sets of planes - calculated in threads.
   // no need for consolidation before draw time.
   // just lines: std::vector<coot::CartesianPairInfo> draw_vector_sets;
   std::vector<coot::density_contour_triangles_container_t> draw_vector_sets;
   std::vector<std::pair<int, TRIANGLE> > map_triangle_centres; // with associated mid-points and indices
   void sort_map_triangles(const clipper::Coord_orth &eye_position);
   static void depth_sort();
   // we only need to sort the triangles if the eye position has moved. So store the eye position
   // of the previous time the triangles were sorted
   clipper::Coord_orth previous_eye_position;

   static std::atomic<bool> draw_vector_sets_lock; // not here because implicitly deleted copy constructor(?)
   // const coot::CartesianPair* diff_map_draw_vectors;
   // int n_diff_map_draw_vectors;
   std::vector<coot::density_contour_triangles_container_t> draw_diff_map_vector_sets;

   void set_use_vertex_gradients_for_map_normals(bool state);
   bool use_vertex_gradients_for_map_normals_flag;

   coot::Cartesian  centre_of_molecule() const;
   float size_of_molecule() const; // return the standard deviation of
				   // the length of the atoms from the
				   // centre of the molecule



   // return -1 if atom not found.
   int atom_index(const char *chain_id, int iresno, const char *atom_id);

   // atom labels
   //

   void gtk3_draw();

   // These are not const because set_bond_colour_by_mol_no() gets called.
   // maybe needs fixing.  Similarly  set_symm_bond_colour_mol_and_symop().
   bool display_stick_mode_atoms_flag;
   void display_bonds(bool against_a_dark_background);
   void display_symmetry_bonds();
   void display_bonds(const graphical_bonds_container &bonds_box, float bond_width_in, bool against_a_dark_background);
   void display_bonds_stick_mode_atoms(const graphical_bonds_container &bonds_box,
				       const coot::Cartesian &front,
				       const coot::Cartesian &back,
				       bool against_a_dark_background);

   void draw_ghost_bonds(int ighost);
   void set_display_stick_mode_atoms(bool f) {
      display_stick_mode_atoms_flag = f;
   }

   void debug_ghosts() const;

   std::vector<int> labelled_atom_index_list;
   // a functor to remove them
   class labelled_atom_remover {
   public:
      int max_idx;
      labelled_atom_remover(int max_idx_in) { max_idx = max_idx_in; }
      bool operator()(int idx) const { return (idx >= max_idx); }
   };

   //
   // Symmetery atom labels.
   //
   //    int* labelled_symm_atom_index_list;
   std::vector<int> labelled_symm_atom_index_list;
   //    int  n_labelled_symm_atoms;
   std::vector<std::pair<symm_trans_t, Cell_Translation> > labelled_symm_atom_symm_trans_;

   // Atom Labelling Interface Functions
   //
   void add_to_labelled_atom_list(int atom_index);

   // Return the atom index of the i'th atom to be labelled
   // (0 indexed).
   //
   int labelled_atom(int i);
   //
   // how many atoms do we want to index?
   // int max_labelled_atom();  // old pointer stuff
   //
   void unlabel_atom(int i);
   // is the i'th atom in the list of atoms to be labelled?
   //
   bool is_in_labelled_list(int i);
   //
   void unlabel_last_atom(); // remove the last atom from the list (if
			     // there *are* atoms in the list, else do
			     // nothing)/

   void trim_atom_label_table(); // when we delete a residue and have
				 // new atoms, we don't want to label
				 // atoms that are over the end of the
				 // (atom label) list.


   //
   void add_atom_to_labelled_symm_atom_list(int atom_index,
					    const symm_trans_t &symm_trans,
					    const Cell_Translation &pre_shift_cell_trans);

   // Return the atom index of the i'th symmetry atom to be labelled.
   //
   int labelled_symm_atom(int i);
   //
   // unlabell the i'th symm atom:
   //
   void unlabel_symm_atom(int i);
   //
   // int max_labelled_symm_atom(); // old pointer stuff
   bool is_in_labelled_symm_list(int i);
   std::pair<symm_trans_t, Cell_Translation> labelled_symm_atom_symm_trans(int i);
   //
   // Added from c-interface (the guile interface) simply named
   // functions, add a label to the atom with the characteristics
   // (using atom_index).
   //
   int    add_atom_label(const char *chain_id, int iresno, const char *atom_id);
   int remove_atom_label(const char *chain_id, int iresno, const char *atom_id);
   void remove_atom_labels(); // and symm labels
   int add_atom_labels_for_residue(mmdb::Residue *residue_p);

   void add_labels_for_all_CAs();

   void local_b_factor_display(bool state, const coot::Cartesian &screen_centre);

   // xmap information
   //
   // We have problems using static vectors, so to Kevin's irritation,
   // we will use a pointer to xmaps that gets renewed, copied and
   // deleted (if necessary).
   //
   // We use xmap_is_filled[0] to see if this molecule is a map.
   //
   clipper::Xmap<float> xmap;
   bool xmap_is_diff_map;
   clipper::NXmap<float> nxmap;
   bool is_patterson;  // for (at least) contour level protection
   std::string map_name;


   float contour_level;
   short int contour_by_sigma_flag;
   float contour_sigma_step;
   //
   GdkRGBA map_colour;
   GdkRGBA map_colour_negative_level;
   GdkRGBA previous_map_colour;
   void save_previous_map_colour();
   void restore_previous_map_colour();
   GdkRGBA radius_to_colour(float radius, float min_radius, float max_radius);
   GdkRGBA fraction_to_colour(float fraction);

   float other_map_for_colouring_min_value;
   float other_map_for_colouring_max_value;
   std::vector<coot::colour_t> other_map_for_colouring_colour_table;
   // use the above values to generate a colour given a value (typically, a correlation)
   GdkRGBA value_to_colour_using_colour_table(float value);

   std::vector<coot::display_list_object_info> display_list_tags;
   void update_map_internal();
   void update_map(bool auto_recontour_map_flag);
   void compile_density_map_display_list(short int first_or_second);

   void draw_surface();
   void draw_dipoles() const;
   bool has_display_list_objects();
   int draw_display_list_objects(int GL_context); // return number of display list objects drawn
   // return the display list object index
   int make_ball_and_stick(const std::string &atom_selection_str,
 			   float bond_thickness, float sphere_size,
 			   bool do_spheres_flag, gl_context_info_t gl_info,
			   const coot::protein_geometry *geom);
   // return the display list object info
   coot::display_list_object_info
   make_ball_and_stick(const std::string &atom_selection_str,
		       float bond_thickness, float sphere_size,
		       bool do_spheres_flag, bool is_second_context,
		       coot::display_list_object_info dloi,
		       const coot::protein_geometry *geom);
   void clear_display_list_object(GLuint tag);

   // the charges for the surface come from the dictionary.
   void make_surface(int on_off_flag, const coot::protein_geometry &geom, float colour_scale);
   // the interface function, converting between a residue specs set
   // and a selection handle
   void make_surface(const std::vector<coot::residue_spec_t> &res_specs_vec,
		     const coot::protein_geometry &geom,
		     float col_scale);
   // SelHnd_selection is the selection of the environment (residues
   // of the active site, say).  SelHnd_all is all tha atoms
   // contributing to the charge (typically all the atoms of the
   // chain).
   void make_surface(int SelHnd_selection, int SelHnd_all, const coot::protein_geometry &geom,
		     float col_scale);

   // a generic function to convert from a residue_spec_vec to a
   // selection handle. Caller creates the SelHnd_selection so that it
   // is clearer where the SelHnd_selection should be deleted.
   //
   void fill_residue_selection(int SelHnd_selection,
			       const std::vector<coot::residue_spec_t> &res_specs_vec,
			       bool allow_waters_flag);


   // void dynamically_transform(coot::CartesianPairInfo v);
   void dynamically_transform(coot::density_contour_triangles_container_t *dctc);

   void clear_draw_vecs();
   void clear_diff_map_draw_vecs();

   bool map_contours_outdated;

   // void add_draw_vecs_to_set(const coot::CartesianPairInfo &cpi);

   // for negative the other map.
   //
   void set_diff_map_draw_vecs(const coot::CartesianPair* c, int n);

   void update_map_triangles(float radius, coot::Cartesian centre);


   // skeleton
   //
   int greer_skeleton_draw_on;
   int fc_skeleton_draw_on;
   void draw_skeleton(bool is_dark_background);
   // void update_skeleton();  bye. use update_clipper_skeleton instead
   graphical_bonds_container greer_skel_box;
   graphical_bonds_container fc_skel_box;
   void draw_fc_skeleton();
   void unskeletonize_map();
   void set_skeleton_bond_colour(float f);
   void set_colour_skeleton_by_segment(); // use random colouring
   void set_colour_skeleton_by_level(); // use given colour (by shading it)
   void update_fc_skeleton_old();
   void update_clipper_skeleton();
   clipper::Xmap<int> xskel_cowtan;
   short int xskel_is_filled;

   void reverse_map();

   //
   void update_map_colour_menu_maybe(int imol);

   void handle_map_colour_change(GdkRGBA map_col,
                                 bool swap_difference_map_colours_flag,
                                 bool main_or_secondary,
                                 clipper::Coord_orth centre,
                                 float radius);

   void handle_map_colour_change_rotate_difference_map(bool swap_difference_map_colours_flag);

   int next_free_map();

   //
   void check_static_vecs_extents();

   //
   int read_ccp4_map(std::string f, int is_diff_map_flag,
		     const std::vector<std::string> &map_glob_extensions); // return -1 on error

   //
   int make_map_from_phs(std::string pdb_filename,
			 std::string phs_filename);

   int make_map_from_phs_using_reso(std::string phs_filename,
				    const clipper::Spacegroup &sg,
				    const clipper::Cell &cell,
				    float reso_limit_low,
				    float reso_limit_high,
				    float map_sampling_rate);

   int make_map_from_phs(const clipper::Spacegroup &sg,
			 const clipper::Cell &cell,
			 std::string phs_filename);

   int make_map_from_cns_data(const clipper::Spacegroup &sg,
			      const clipper::Cell &cell,
			      std::string cns_data_filename);

   int make_map_from_cif_sigmaa(int imol_no_in,
				std::string cif_filename, int sigma_map_type); // has phases
   int make_map_from_cif(int imol_no_in,
			 std::string cif_file_name); // uses above with type SIGMAA
   int make_map_from_cif_diff_sigmaa(int imol_no_in,
				     std::string cif_file_name); // uses above with TYPE_DIFF_SIGMAA

   int make_map_from_cif(int imol_no_in,
			 std::string cif_filename,  // generate phases from mol
			  int imol_coords);
   int make_map_from_cif(int imol_no_in,
			 std::string cif_file_name,
			 atom_selection_container_t SelAtom);
   int make_map_from_cif_2fofc(int imol_no_in,
			       std::string cif_file_name,
			       atom_selection_container_t SelAtom);
   int make_map_from_cif_2fofc(int imol_no_in,
			       std::string cif_file_name,
			       int imol_coords);
   int make_map_from_cif_fofc(int imol_no_in,
			      std::string cif_file_name,
			      int imol_coords);
   int make_map_from_cif_generic(int imol_no_in,
				 std::string cif_file_name,
				 atom_selection_container_t SelAtom,
				 short int is_2fofc_type);
   int make_map_from_cif_nfofc(int imol_map_in,
			       std::string cif_file_name,  // Virgina request
			       int map_type,
			       short int swap_difference_map_colours);
   int calculate_sfs_and_make_map(int imol_no_in,
				  const std::string &mol_name,
				  const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &myfsigf,
				  atom_selection_container_t SelAtom,
				  short int is_2fofc_type);
   int make_map_from_mtz_by_calc_phases(int imol_no_in,
					const std::string &mtz_file_name,
					const std::string &f_col,
					const std::string &sigf_col,
					atom_selection_container_t SelAtom,
					short int is_2fofc_type);
   // use this molecules mol and the passed data to make a map for some other
   // molecule
   int sfcalc_genmap(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                     const clipper::HKL_data<clipper::data32::Flag> &free,
                     clipper::Xmap<float> *xmap_p);
   coot::util::sfcalc_genmap_stats_t
   sfcalc_genmaps_using_bulk_solvent(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                     const clipper::HKL_data<clipper::data32::Flag> &free,
                                     clipper::Xmap<float> *xmap_2fofc_p,
                                     clipper::Xmap<float> *xmap_fofc_p);
   void fill_fobs_sigfobs(); // re-reads MTZ file (currently 20210816-PE)
   bool sanity_check_atoms(mmdb::Manager *mol); // sfcalc_genmap crashes after merge of ligand.
                                                // Why? Something wrong with the atoms after merge?
                                                // Let's diagnose.... Return false on non-sane.


   void update_map_in_display_control_widget() const;
   void new_coords_mol_in_display_control_widget() const;  // for a new molecule.
   void update_mol_in_display_control_widget() const; // changing the name of
                                                      // this mol (e.g on save).

   void update_mol_in_simple_display_control_menu(GtkWidget *model_menu, int map_coords_mol_flag);

   //
   float map_mean()  const { return map_mean_;  }
   float map_sigma() const { return map_sigma_; } // cached
   float get_map_sigma_current(); // regen stats and update map_sigma_

   map_statistics_t map_statistics() const;

   //
   float sharpen_b_factor() const { return sharpen_b_factor_; }
   float sharpen_b_factor_kurtosis_optimised() const { return sharpen_b_factor_kurtosis_optimised_; }
   void set_sharpen_b_factor_kurtosis_optimised(float b_factor_in) {
      sharpen_b_factor_kurtosis_optimised_ = b_factor_in;
   }
   void save_original_fphis_from_map();

   // for debugging.
   int test_function();

   // output
   // return the mmdb exit status (as write_pdb_file)
   int export_coordinates(std::string filename) const;

   // for labelling from guile
   int atom_spec_to_atom_index(std::string chain, int reso,
			       std::string atom_name) const;

   // return -1 if atom not found
   //
   int full_atom_spec_to_atom_index(const std::string &chain,
				    int reso,
				    const std::string &insersion_code,
				    const std::string &atom_name,
				    const std::string &alt_conf) const;

   // return -1 if atom not found
   //
   int full_atom_spec_to_atom_index(const coot::atom_spec_t &spec) const;

   int atom_to_atom_index(mmdb::Atom *at) const;

   // Does atom at from moving atoms match atom_sel.atom_selection[this_mol_index_maybe]?
   // or has atom_sel changed in the mean time?
   bool moving_atom_matches(mmdb::Atom *at, int this_mol_index_maybe) const;

   int atom_index_first_atom_in_residue(const std::string &chain_id,
					int iresno,
					const std::string &ins_code) const;

   int atom_index_first_atom_in_residue(const std::string &chain_id,
					int iresno,
					const std::string &ins_code,
					const std::string &altconf) const;

   int atom_index_first_atom_in_residue_internal(const std::string &chain_id,
						 int iresno,
						 const std::string &ins_code,
						 const std::string &altconf,
						 bool test_alt_conf_flag) const;

   void install_ghost_map(const clipper::Xmap<float> &mapin, std::string name,
			  const coot::ghost_molecule_display_t &ghost_info,
			  int is_diff_map_flag,
			  int swap_difference_map_colours_flag,
			  float contour_level_in);

   void install_new_map(const clipper::Xmap<float> &mapin, std::string name, bool is_em_map_in);

   void install_new_map_with_contour_level(const clipper::Xmap<float> &mapin, std::string name, float contour_level, bool is_em_map_in);

   void set_name(std::string name); // you are encouraged not to use
				    // this (only for use after having
				    // imported an xmap).


   // regularization results:
   void replace_coords(const atom_selection_container_t &asc,
		       bool change_altconf_occs_flag,
		       bool replace_coords_with_zero_occ_flag);
   // helper function for above function
   bool movable_atom(mmdb::Atom *mol_atom, bool replace_coords_with_zero_occ_flag) const;

   int add_terminal_residue_using_phi_psi(const std::string &chain_id,
					  int res_no,
					  const std::string &residue_type,
					  float phi, float psi);

   // either rama-search for a protein residue or simple add for nucleic acid.
   void add_terminal_residue_wrapper(const coot::residue_spec_t &res_spec,
                                     const std::string &residue_type);

   // When a new residue is added to the C-terminus of a chain/fragment, we will need to move
   // the O of this one to make a proper peptide plane (the position of the next residue
   // was not dependent on the position of the O of this one).
   // (note: read as 'added-to' residue)
   void move_O_atom_of_added_to_residue(mmdb::Residue *res_p, const std::string &chain_id);

   // extra modelling results, e.g. waters, terminal residues, etc
   void add_coords(const atom_selection_container_t &asc);
   //
   // Add a residue "in the middle" of some other residues
   // (uses mmdb::Chain::InsResidue()).
   //
   // This relies on the mol of the asc being different to the mol of the
   // atom_sel.
   //
   // If it is the same, then error and do nothing.
   //
   // Typically, this is called by fit terminal residue, which has its
   // mol created from a pcmmdbmanager() from the molecule of the
   // residue, so this is fine in this case.
   //
   void insert_coords(const atom_selection_container_t &asc);
   void insert_coords_change_altconf(const atom_selection_container_t &asc);
   // a utility function for insert_coords, return the residue (that
   // is currently at res_index) too.  Return a index of -1 (and a
   // mmdb::Residue of NULL) when no residue found.

   std::pair<int, mmdb::Residue *> find_serial_number_for_insert(int seqnum_for_new,
								 const std::string &ins_code_for_new,
								 const std::string &chain_id) const;

   void update_molecule_to(std::vector<coot::scored_skel_coord> &pos_position);

   //
   // Add this molecule (typically of waters to this
   // molecule... somehow).  All the atoms of water_mol need to be in
   // a chain that has a different chain id to all the chains in this
   // molecule.  Else fail (return status 0).
   //
   // But we should try to put the waters into (add/append to) a chain
   // of waters in this molecule, if it has one.
   //
   int insert_waters_into_molecule(const coot::minimol::molecule &water_mol, const std::string &res_name);
   int append_to_molecule(const coot::minimol::molecule &water_mol);
   mmdb::Residue *residue_from_external(int reso, const std::string &insertion_code,
					const std::string &chain_id) const;

   // used to be a const ref. Now return the whole thing!. Caller must call
   // fill_fobs_sigfobs() directly before using this function - meh, not a good API.
   // Return a *pointer* to the data so that we don't get this hideous non-reproducable
   // crash when we access this data item after the moelcule vector has been resized
   // 20210816-PE.
   clipper::HKL_data<clipper::data32::F_sigF> *get_original_fobs_sigfobs() const {
      if (!original_fobs_sigfobs_filled) {
         std::string m("Original Fobs/sigFobs is not filled");
         throw(std::runtime_error(m));
      }
      return original_fobs_sigfobs_p;
   }

   clipper::HKL_data<clipper::data32::Flag> *get_original_rfree_flags() const {
      if (!original_fobs_sigfobs_filled) {
         std::string m("Original Fobs/sigFobs is not filled - so no RFree flags");
         throw(std::runtime_error(m));
      }
      return original_r_free_flags_p;
   }


   // for the "Render As: " menu items:
   //
   void bond_representation(const coot::protein_geometry *geom_p, bool rebonding_is_needed);
   //

   float bonds_colour_map_rotation; // OpenGL1
   void update_bonds_colour_using_map_rotation(float f); // modern OpenGL
   short int bonds_rotate_colour_map_flag;
   int Bonds_box_type() const { return bonds_box_type; }

   //
   void pepflip(int atom_index);
   int pepflip_residue(const std::string &chain_id,
		       int ires_seqno,
		       const std::string &ins_code,
		       const std::string &alt_conf);
   void pepflip(const coot::atom_spec_t &atom_spec_t);

   int do_180_degree_side_chain_flip(const std::string &chain_id,
				     int resno,
				     const std::string &inscode,
				     const std::string &altconf,
				     coot::protein_geometry *geom_p);

   int do_180_degree_side_chain_flip_protein(const std::string &chain_id,
					     int resno,
					     const std::string &inscode,
					     const std::string &altconf,
					     coot::protein_geometry *geom_p);

   int do_180_degree_side_chain_flip_nucleic_acid(const std::string &chain_id,
						  int resno,
						  const std::string &inscode,
						  const std::string &altconf,
						  coot::protein_geometry *geom_p);

   // return "N', "C" or "not-terminal-residue"
   std::string get_term_type_old(int atom_index);
   std::string get_term_type(int atom_index) const;
   std::string get_term_type(mmdb::Atom *atom) const;
   // by alignment (against asigned pir seq file) return, "HIS", "ALA" etc, if we can.
   std::pair<bool, std::string> find_terminal_residue_type(const std::string &chain_id, int resno,
							   mmdb::realtype alignment_wgap,
							   mmdb::realtype alignment_wspace,
							   bool is_nucleic_acid_flag = false) const;

   std::vector<std::string> get_types_in_molecule() const;


   //
   graphical_bonds_container make_environment_bonds_box(int atom_index,
							coot::protein_geometry *protein_geom_p) const;
   graphical_bonds_container make_symmetry_environment_bonds_box(int atom_index,
								 coot::protein_geometry *protein_geom_p) const;
   graphical_bonds_container make_environment_bonds_box(const coot::residue_spec_t &residue_spec,
							coot::protein_geometry *protein_geom_p) const;

   bool has_xmap() const { return ! xmap.is_null(); }

   bool has_nxmap() const { return ! nxmap.is_null(); }

   bool has_model() const { return (atom_sel.n_selected_atoms > 0); }

   // for use when we want to tell the difference between a molecule
   // with no atoms that is open and one that is closed. I guess that
   // most has_model() usaged should be changed to open_molecule_p()
   // usage...
   bool open_molecule_p() const {
      if (atom_sel.mol)
	 return 1;
      else
	 return 0;
   }

   bool is_displayed_p() const {
      bool i;
      if (has_model()) {
	 if (draw_it) {
	    i = 1;
	 } else {
	    i = 0;
	 }
      } else {
	 if (has_xmap()) { // NXMAP-FIXME
	    if (draw_it_for_map) {
	       i = 1;
	    } else {
	       i = 0;
	    }
	 } else {
	    i = 0;
	 }
      }
      return i;
   }

   bool map_is_too_blue_p() const;


   // delete residue, typically for waters
   //
   // return the success status (0: failure to delete, 1 is deleted)
   //
   // if you want to delete a/the residue without having to specify
   // the model number (which is typically the case) pass mmdb::MinInt4 as
   // the model_number.
   //
   short int delete_residue(int model_number,
			    const std::string &chain_id, int resno,
                            const std::string &inscode);

   // delete from all models
   // wraps above
   short int delete_residue(const coot::residue_spec_t &spec);

   // Delete only the atoms of the residue that have the same altconf (as
   // the selected atom).  If the selected atom has altconf "", you
   // should call simply delete_residue().
   //
   // Return 1 if at least one atom was deleted, else 0.
   //
   short int delete_residue_with_full_spec(int imodel,
					   const std::string &chain_id,
					   int resno,
					   const std::string &inscode,
					   const std::string &altconf);

   // Return 1 if at least one atom was deleted, else 0.
   //
   short int delete_residues(const std::vector<coot::residue_spec_t> &specs);

   short int delete_residue_sidechain(const std::string &chain_id,
				      int resno,
				      const std::string &inscode);
   short int delete_residue_sidechain(const coot::residue_spec_t &rs);

   bool delete_atom(const std::string &chain_id,
		    int resno,
		    const std::string &ins_code,
		    const std::string &atname,
		    const std::string &altconf);

   bool delete_atom(const coot::atom_spec_t &atom_spec);

   short int delete_residue_hydrogens(const std::string &chain_id, int resno,
				      const std::string &ins_code,
				      const std::string &altloc);

   int delete_atoms(const std::vector<coot::atom_spec_t> &atom_specs);

   int delete_hydrogens(); // return status of atoms deleted (0 -> none deleted).

   int delete_waters(); // return status of atoms deleted (0 -> none deleted).

   int delete_chain(const std::string &chain_id);

   bool delete_sidechain(mmdb::Residue *residue_p);

   int delete_sidechains_for_chain(const std::string &chain_id);

   int delete_sidechain_range(const coot::residue_spec_t &res_1,
			      const coot::residue_spec_t &res_2);

   int delete_water(const coot::atom_spec_t &atom_spec); // residue type is checked here

   // closing molecules, delete maps and atom sels as appropriate
   // and unset "filled" variables.  Set name_ to "".
   void close_yourself();

   // "Interactive" function that does make_backup().
   int mutate(int atom_index, const std::string &residue_type, short int do_stub_flag);
   // an alternate interface to the above:
   int mutate(int resno, const std::string &insertion_code,
	       const std::string &chain_id, const std::string &residue_type);
   // and another:
   int mutate(mmdb::Residue *res, const std::string &residue_type, bool verbose=true);

   // Here is something that does DNA/RNA
   int mutate_base(const coot::residue_spec_t &res_spec, std::string type,
		   bool use_old_style_naming);

   int mutate_by_overlap(const std::string &chain_id, int res_no, const std::string &new_type);

   // and the biggie: lots of mutations/deletions/insertions from an
   // alignment:
   //
   // We return a mutation container so that the calling function can
   // do a autofit rotamer on the mutated residues.
   //
   coot::chain_mutation_info_container_t align_and_mutate(const std::string chain_id,
							  const coot::fasta &seq,
							  bool renumber_residues_flag,
							  mmdb::realtype wgap,
							  mmdb::realtype wspace);

   void mutate_chain(const std::string &chain_id,
		     const coot::chain_mutation_info_container_t &mut_cont_info,
		     mmdb::PResidue *SelResidues,
		     int nSelResidues,
		     bool renumber_residues_flag);

   std::pair<bool, std::string>
   residue_type_next_residue_by_alignment(const coot::residue_spec_t &clicked_residue,
					  mmdb::Chain *clicked_residue_chain_p,
					  short int is_n_term_addition,
					  mmdb::realtype alignment_wgap,
					  mmdb::realtype alignment_wspace) const;

   // return a status flag (alignments done)
   std::pair<bool, std::vector<coot::chain_mutation_info_container_t> >
   residue_mismatches(mmdb::realtype alignment_wgap, mmdb::realtype aligment_wspace) const;

   // 20180302 add PIR alignment parsing
   void associate_pir_alignment(const std::string &chain_id, const std::string &alignment);
   // apply the alignment
   void apply_pir_alignment(const std::string &chain_id);
   void apply_pir_renumber(const coot::pir_alignment_t &a, mmdb::Chain *chain_p);
   // this is where the PIR alignments are stored, the key is the chain-id
   std::map<std::string, coot::pir_alignment_t> pir_alignments;


   // Try to align on all chains - pick the best one and return it in
   // the second.  If there is no chain that matches within match_frag
   // (e.g. 0.95) then return 0 as first and a blank in second.
   //
   std::pair<bool, std::pair<std::string, coot::chain_mutation_info_container_t> >
   try_align_on_all_chains(const std::string &target, float match_fragment,
			   mmdb::realtype wgap, mmdb::realtype wspace) const;

   void make_backup_from_outside(); // when we have a multi mutate, we
				    // want the wrapper to make a
				    // backup when we start and set
				    // changes when when finish.
				    // Rather crap that this needs to
				    // be done externally, I think.
   bool backups_state() const { return backup_this_molecule; }
   void set_have_unsaved_changes_from_outside();

   void mutate_internal(mmdb::Residue *residue, mmdb::Residue *std_residue, const std::string &alt_conf);
   bool progressive_residues_in_chain_check_by_chain(const char *chain_id) const;
   void mutate_base_internal(mmdb::Residue *residue, mmdb::Residue *std_base, bool use_old_style_naming);

   // multiple mutate: we don't want to put the results into an asc,
   // we simply want to make the replacement withouth user
   // intervention.  This is the externally visible function (via
   // c-interface imol selector wrapper).
   //
   // Importantly, this function is called by some sort of wrapper
   // that is doing multiple mutations and therefore doesn't do a
   // backup.  However, backup should be done in the wrapping function
   //
   int mutate_single_multipart(int ires_serial, const std::string &chain_id,
			       const std::string &target_res_type);

   // mutate and autofit the residues
   //
   int nudge_residue_sequence(const std::string &chain_id,
			      int resno_range_start,
			      int resno_range_end,
			      int offset,
			      short int nudge_residue_numbers_also);

   //
   // and the functions that mutate functions uses:
   // (which returns success status - up from get_ori_to_this_res):
   //
   short int move_std_residue(mmdb::Residue* moving_residue, const mmdb::Residue *reference_residue) const;
   //
   // Get a deep copy:
// return NULL on failure
   mmdb::Residue *get_standard_residue_instance(const std::string &res_type);

   // Return the atom index of the "next" atom
   // -1 on failure.
   int intelligent_next_atom(const std::string &chain,
			     int resno,
			     const std::string &atom_name,
			     const std::string &ins_code,
			     const coot::Cartesian &rc);
   int intelligent_previous_atom(const std::string &chain,
				 int resno,
				 const std::string &atom_name,
				 const std::string &ins_code,
				 const coot::Cartesian &rc);
   mmdb::Residue *next_residue_missing_residue(const coot::residue_spec_t &spec) const;
   mmdb::Atom *atom_intelligent(const std::string &chain_id, int resno,
			   const std::string &ins_code) const;

   // is point close (< 1A) to any atom in the given residue?
   bool close_to_residue(mmdb::Residue *residue_p, coot::Cartesian point) const;

   // If there is a CA in this residue then return the index of that
   // atom, if not, then return the index of the first atom in the
   // residue.
   //
   // Return -1 on no atoms in residue.
   //
   int intelligent_this_residue_atom(mmdb::Residue *res_p) const;
   coot::atom_spec_t intelligent_this_residue_atom(const coot::residue_spec_t &rs) const;
   mmdb::Atom *intelligent_this_residue_mmdb_atom(mmdb::Residue *res_p) const; // has null res_p protection.

   // pointer atoms:
   void add_pointer_atom(coot::Cartesian pos);
   // if there is another atom very close, then don't allow the addition of a new atom
   // return status false denotes failure.
   std::pair<bool,std::string> add_typed_pointer_atom(coot::Cartesian pos, const std::string &type);

   // dummy atoms (not bonded)
   void add_dummy_atom(coot::Cartesian pos);
   void add_pointer_multiatom(mmdb::Residue *res_p, const coot::Cartesian &pos, const std::string &type);
   void add_multiple_dummies(mmdb::Chain *chain_p,
			     const std::vector<coot::scored_skel_coord> &pos_position);
   void add_multiple_dummies(const std::vector<coot::scored_skel_coord> &pos_position);
   void add_multiple_dummies(const std::vector<coot::Cartesian> &positions);




   // Return a vector of upto 3 positions of the most latestly added
   // atoms with the most lastest atom addition (that is the passed atom)
   // in the back() slot of the vector.
   //
   // Are we looking back along the chain (i.e. we are building forward (direction = 1))
   // or building backward (direction = 0)?
   //
   std::vector<clipper::Coord_orth> previous_baton_atom(mmdb::Atom* latest_atom_addition,
							short int direction) const;

   clipper::Xmap<coot::SkeletonTreeNode> skeleton_treenodemap;
   short int skeleton_treenodemap_is_filled;
   void fill_skeleton_treenodemap();

   // baton atoms:
   //
   // 20100125 now the chain-id is passed
   //
   mmdb::Atom *add_baton_atom(coot::Cartesian pos,
			 int i_chain_start_resno,
			 const std::string &chain_id,
			 short int i_start_resno_active, // dont ignore it this time?
			 short int direction_flag); // return a pointer to
                                                    // the just added atom.

   // Return the position of the previous atom.
   //
   std::pair<short int, mmdb::Atom *> baton_build_delete_last_residue();

   std::vector<coot::scored_skel_coord>
   next_ca_by_skel(const std::vector<clipper::Coord_orth> &previous_ca_positions,
		   const clipper::Coord_grid &coord_grid_start,
		   short int use_coord_grid_start_flag,
		   float ca_ca_bond_length,
		   float map_cut_off,
		   int max_skeleton_search_depth) const;
   std::pair<short int, clipper::Coord_grid> search_for_skeleton_near(const coot::Cartesian &pos) const;
   float density_at_point(const clipper::Coord_orth &pt) const; // return the values of map_list[0]
                                                                // at that point.


   //
   std::vector<std::string> save_state_command_strings() const {
      return save_state_command_strings_;
   }


   void set_map_colour(GdkRGBA col) { map_colour = col; update_map(true); /* for now */ }
   std::pair<GdkRGBA, GdkRGBA> get_map_colours() const;

   std::vector<std::string> set_map_colour_strings() const;
   void colour_map_using_map(const clipper::Xmap<float> &xmap);
   void colour_map_using_map(const clipper::Xmap<float> &xmap, float table_bin_start, float table_bin_size,
                             const std::vector<coot::colour_t> &colours);
   const clipper::Xmap<float> *other_map_for_colouring_p;
   void turn_off_other_map_for_colouring() {
      other_map_for_colouring_p = NULL;
      colour_map_using_other_map_flag = false;
   }
   fresnel_settings_t fresnel_settings;
   void set_fresnel_colour(const glm::vec4 &col_in);

   // save yourself and update have_unsaved_changes_flag status
   //
   int save_coordinates(const std::string &filename,
			bool save_hydrogens=1,
			bool save_aniso_records=1,
			bool save_conect_records=0);
   int quick_save(); // save to default file name if has unsaved changes.  Return non-zero on problem.
   std::string stripped_save_name_suggestion(); // sets coot_save_index maybe
   int Have_unsaved_changes_p() const;
   bool Have_modifications_p() const { return history_index > 0 ? 1 : 0;}
   bool Have_redoable_modifications_p() const ;
   int get_history_index() const;
   void turn_off_backup() { backup_this_molecule = false; }
   void turn_on_backup()  { backup_this_molecule = true;  }
   int apply_undo(const std::string &cwd);
   int apply_redo(const std::string &cwd);
   // Called from outside, if there is a backup file more recent
   // than the file name_ then restore from it.
   //
   // Return 1 if the restore happened, 0 if not.
   //
   short int execute_restore_from_recent_backup(std::string backup_file_name,
						std::string cwd);
   coot::backup_file_info_t recent_backup_file_info() const;

   // For model view (go to atom)
   //
   std::vector<coot::model_view_residue_button_info_t>  // old style
   model_view_residue_button_labels() const;

   std::vector<coot::model_view_atom_tree_chain_t>
   model_view_residue_tree_labels(bool include_water_residue_flag, bool ligands_ony_flag) const;

   std::vector<coot::model_view_atom_button_info_t>
   model_view_atom_button_labels(const std::string &chain_id,
				 int seqno,
				 const std::string &ins_code) const;

   // return the number of residues in chain with chain_id, return -1 on error
   //
   int chain_n_residues(const char *chain_id) const;

   // return the number of residues in the molecule. return -1 on error.
   int n_residues() const;

   // return the number of atoms in the molecule. return -1 on error.
   int n_atoms() const;

   // Fourier stuff
   std::string Fourier_f_label()      const { return fourier_f_label; }
   std::string Fourier_phi_label()    const { return fourier_phi_label; }
   std::string Fourier_weight_label() const { return fourier_weight_label; }

   // refmac
   //
   // sets private have_sensible_refmac_params
   //
   void store_refmac_params(const std::string &mtz_filename,
			    const std::string &fobs_col,
			    const std::string &sigfobs_col,
			    const std::string &r_free_col,
			    int r_free_flag_sensible);
   std::vector<coot::atom_attribute_setting_help_t> get_refmac_params() const;

   void store_refmac_mtz_filename(const std::string &mtz_filename);

   void store_refmac_phase_params(const std::string &phi,
				  const std::string &fom,
				  const std::string &hla,
				  const std::string &hlb,
				  const std::string &hlc,
				  const std::string &hld);

   // more refmac stuff
   //
   int write_pdb_file(const std::string &filename); // not const because of shelx/name manip
   int write_cif_file(const std::string &filename); // not const because of shelx/name manip
   short int Have_sensible_refmac_params() const { return has_xmap() && have_sensible_refmac_params; }
   short int Have_refmac_phase_params()    const { return have_refmac_phase_params; }
   void increment_refmac_count() { refmac_count++; }
   int Refmac_count() const { return refmac_count; }
   std::string Refmac_mtz_filename() const { return refmac_mtz_filename; }
   std::string Refmac_file_mtz_filename() const { return refmac_file_mtz_filename; }
   std::string Refmac_fobs_col() const { return refmac_fobs_col; }
   std::string Refmac_sigfobs_col() const { return refmac_sigfobs_col; }
   std::string Refmac_r_free_col() const { return refmac_r_free_col; }
   std::string Refmac_phi_col() const { return refmac_phi_col; }
   std::string Refmac_fom_col() const { return refmac_fom_col; }
   std::string Refmac_hla_col() const { return refmac_hla_col; }
   std::string Refmac_hlb_col() const { return refmac_hlb_col; }
   std::string Refmac_hlc_col() const { return refmac_hlc_col; }
   std::string Refmac_hld_col() const { return refmac_hld_col; }
   short int Refmac_r_free_sensible() const { return refmac_r_free_flag_sensible; }
   void set_refmac_counter(int i) { refmac_count = i; }
   void set_refmac_save_state_commands(std::string mtz_file_name,
				       std::string f_col,
				       std::string phi_col,
				       std::string weight_col,
				       bool use_weights,
				       bool is_diff_map,
				       std::string refmac_fobs_col,
				       std::string refmac_sigfobs_col,
				       std::string refmac_r_free_col,
				       bool refmac_r_free_flag_sensible);
   std::string Refmac_name_stub() const;
   std::string Refmac_in_name() const;
   std::string Refmac_out_name() const;
   std::string Refmac_mtz_out_name() const;
   std::string name_sans_extension(short int include_path_flag) const;
   mmdb::Manager *get_residue_range_as_mol(const std::string &chain_id,
					  int resno_start,
					  int resno_end) const;


   // (best fit) rotamer stuff:
   //
   // This is the generic interface.  Inside auto_fit_best_rotamer, we
   // look at the rotamer_seach_mode and decide there if we want to
   // call the function that runs the backrub_rotamer mode.
   //
   float auto_fit_best_rotamer(int rotamer_seach_mode,
			       int resno,
			       const std::string &altloc,
			       const std::string &insertion_code,
			       const std::string &chain_id, int imol_map, int clash_flag,
			       float lowest_probability,
			       const coot::protein_geometry &pg);
   // (best fit) rotamer stuff:
   //
   float auto_fit_best_rotamer(int resno,
			       const std::string &altloc,
			       const std::string &insertion_code,
			       const std::string &chain_id, int imol_map, int clash_flag,
			       float lowest_probability,
			       const coot::protein_geometry &pg);

   // interface from atom picking (which simply gets the resno, altloc
   // etc form the atom_index and calls the above function.
   float auto_fit_best_rotamer(int rotamer_search_mode,
			       int atom_index, int imol_map, int clash_flag,
			       float lowest_probability,
			       const coot::protein_geometry &pg);

   // Internal.  Return succes status and score (we need status
   // because on failure, we should fall back to conventional rotamer
   // search).
   std::pair<bool,float> backrub_rotamer(const std::string &chain_id,
					 int res_no,
					 const std::string &ins_code,
					 const std::string &alt_conf,
					 const coot::protein_geometry &pg);

   // calls above
   void backrub_rotamer_residue_range(const std::string &chain_id, int resno_start, int resno_end, const coot::protein_geometry &pg);

   // a chain-base version of the above would be useful (currently a scripting function)


   int set_residue_to_rotamer_number(coot::residue_spec_t res_spec,
				     const std::string &alt_conf,
				     int rotamer_number,
				     const coot::protein_geometry &pg);

   int set_residue_to_rotamer_name(coot::residue_spec_t res_spec,
				   const std::string &alt_conf,
				   const std::string &rotamer_name,
				   const coot::protein_geometry &pg);

   // internal function for above two functions
   int set_residue_to_rotamer_move_atoms(mmdb::Residue *res, mmdb::Residue *moving_res);

   std::vector<coot::named_rotamer_score> score_rotamers(const std::string &chain_id,
							 int res_no,
							 const std::string &ins_code,
							 const std::string &alt_conf,
							 int clash_flag,
							 float lowest_probability,
							 const clipper::Xmap<float> &xmap,
							 const coot::protein_geometry &pg);


   // Add OXT atom:  Return status, 0 = fail, 1 = worked.
   // (use get_residue() to get the residue for this);
   //
   short int add_OXT_to_residue(mmdb::Residue *residue, coot::protein_geometry *geom_p);
   short int add_OXT_to_residue(int reso, const std::string &insertion_code,
				const std::string &chain_id,
				coot::protein_geometry *geom_p); // external usage
   bool residue_has_oxt_p(mmdb::Residue *residue) const; // used by above.  Dont add if returns true.

   std::pair<bool, int>  last_residue_in_chain(const std::string &chain_id) const;
   std::pair<bool, int> first_residue_in_chain(const std::string &chain_id) const;
   std::pair<bool, int>  last_protein_residue_in_chain(const std::string &chain_id) const;

   // return NULL on no last residue.
   mmdb::Residue *last_residue_in_chain(mmdb::Chain *chain_p) const;

   //! @return -1 on failure to find residue
   int residue_serial_number(const std::string &chain_id, int reso, const std::string &insertion_code) const;

   std::string res_name_from_serial_number(std::string chain_id, unsigned int serial_number) const;

   // cell and symmetry swapping
   //
   // return an empty vector on failure, a vector of size 6 on success:
   std::pair<std::vector<float>, std::string> get_cell_and_symm() const;
   //
   void set_mmdb_cell_and_symm(std::pair<std::vector<float>, std::string>);
   // return the success status of the set
   bool set_mmdb_symm(const std::string &s);

   //
   void set_bond_thickness(float t) { bond_width = t; }
   int get_bond_thickness() const { return int(bond_width); }

   //
   void apply_atom_edit(const coot::select_atom_info &sai);
   void apply_atom_edits(const std::vector <coot::select_atom_info> &saiv);


   mmdb::Chain *water_chain() const; // which chain contains just waters
				// (or is empty)?  return 0 if none.

   mmdb::Chain *water_chain_from_shelx_ins() const; // the single chain

   std::pair<bool, std::string> chain_id_for_shelxl_residue_number(int resno) const;

   // return state, max_resno + 1, or 0, 1 of no residues in chain.
   //
   std::pair<short int, int> next_residue_number_in_chain(mmdb::Chain *w,
							  bool new_res_no_by_hundreds=false) const;

   // For environment distance application, we need to find the atom
   // nearest the centre of rotation
   std::pair<float, int> nearest_atom(const coot::Cartesian &pos) const;

   // molecule-class-info-other functions:
   //
   // For edit phi psi, tinker with the passed (moving atoms) asc
   // Return status, -1 on error (unable to change).
   short int residue_edit_phi_psi(atom_selection_container_t residue_asc,
				  int atom_index, double phi, double psi);

   // Return the residue thats for moving_atoms_asc as a molecule.
   //

   // whole_res_flag is for specifying if we want the whole residue
   // (1) or just the atoms of this residue with alt conf (and atoms
   // with altconf "") (0).
   atom_selection_container_t edit_residue_pull_residue(int atom_index,
							short int whole_res_flag);
   clipper::Coord_orth to_coord_orth(mmdb::Atom *atom) const;

   std::pair<double, double> get_phi_psi(int atom_index) const;

   // edit chi angles:
   //
   // Return the number of chi angles for this residue:
   // Do not be mislead this is only a flag, *not* the number of chi
   // angles for this residue, i.e. does this residue have chi angles?
   //
   int N_chis(int atom_index);

   // return the resulting torsion value
   double set_torsion(const std::string &chain_id,
		      int res_no,
		      const std::string &insertion_code,
		      const std::string &alt_conf,
		      const std::string &atom_name_1,
		      const std::string &atom_name_2,
		      const std::string &atom_name_3,
		      const std::string &atom_name_4,
		      double tors,
		      const coot::protein_geometry &geom);

   // manipulate torsion angles of first residue in the molecule to
   // match those of the passed (reference residue (from a different
   // molecule, typically).
   //
   int match_torsions(mmdb::Residue *res_ref,
		      const std::vector <coot::dict_torsion_restraint_t> &tr_ligand,
		      const coot::protein_geometry &geom);

   // So that we can move around all the atoms of a ligand (typically)
   void translate_by(float x, float y, float z);
   void translate_by_internal(const clipper::Coord_orth &co, mmdb::Residue *residue_p);
   void transform_by(mmdb::mat44 mat); // can't make this const: mmdb probs.
   void transform_by(const clipper::RTop_orth &rtop);
   void transform_by(const clipper::RTop_orth &rtop, mmdb::Residue *res);
   // called by above (no backup or bonds update here in internal function).
   void transform_by_internal(const clipper::RTop_orth &rtop, mmdb::Residue *res);
   // we add the arg make_backup_flag so that "splinter mode" doesn't do backups.
   void transform_zone_by(const std::string &chain_id, int resno_start, int resno_end,
			  const std::string &ins_code,
			  const clipper::RTop_orth &rtop,
			  bool make_backup_flag);

   //
   std::string name_for_display_manager() const; // stripped of path maybe
   std::string dotted_chopped_name() const;

   float get_contour_level() const {
      if (! has_xmap()) {
	 return 0;
      } else {
	 // NXMAP-FIXME
	 return contour_level;
      }
   }


   float get_contour_level_by_sigma() const;
   // external contour control (from saving parameters):
   void set_contour_level(float f);
   void set_contour_level_by_sigma(float f);
   short int change_contour(int direction); // return status 0 if it didn't happen.
   // for the state file:
   std::vector <std::string> get_map_contour_strings() const;
   // for the state fle
   std::vector <std::string> get_map_contour_sigma_step_strings() const;
   void set_contour_by_sigma_step(float v, short int state);
   // we ask of this molecule a question: contour_by_sigma?
   short int contoured_by_sigma_p() const;

   // a -other function (return the number of trimmed atoms):
   int trim_by_map(const clipper::Xmap<float> &xmap_in,
		   float map_level, short int delete_or_zero_occ_flag);

   int trim_molecule_by_b_factor(float limit, bool keep_higher_flag);

   void pLDDT_to_b_factor();

   // logical_operator_and_or_flag 0 for AND 1 for OR.
   //
   std::vector <coot::atom_spec_t>
   find_water_baddies(float b_factor_lim, const clipper::Xmap<float> &xmap_in,
		      float map_sigma,
		      float outlier_sigma_level,
		      float min_dist, float max_dist,
		      short int part_occ_contact_flag,
		      short int zero_occ_flag,
		      short int logical_operator_and_or_flag);

   // return a vector of outliers
   std::vector <coot::atom_spec_t>
   check_waters_by_difference_map(const clipper::Xmap<float> &xmap_in,
				  float outlier_sigma_level) const;


   //
   void set_map_is_difference_map(bool flag);
   bool is_difference_map_p() const;


   // Scripting Refinement:
   int does_residue_exist_p(const std::string &chain_id, int resno, const std::string &inscode) const;


   // LMB surface
   // coot::surface *cootsurface; // dead for now

   // for widget label:
   std::string cell_text_with_embeded_newline() const;

   // for pointer distances
   std::vector<clipper::Coord_orth> distances_to_point(const clipper::Coord_orth &pt,
						       double min_dist,
						       double max_dist);

   // we may decide in future that sequence can be more sophisticated
   // than a simple string.
   //
   // pair: chain_id sequence.  We need public access to these now
   // that they are going into the state script.
   std::vector<std::pair<std::string, std::string> > input_sequence;
   bool is_fasta_aa(const std::string &a) const;
   bool is_pir_aa  (const std::string &a) const;

   // add the sequence the file (read depending on file name) to input_sequence vector (chain-id is blank
   // as it could apply to any chain)
   void associate_sequence_from_file(const std::string &seq_file_name);

   // sequence [a -other function]
   void assign_fasta_sequence(const std::string &chain_id, const std::string &seq); // add to input_sequence vector

   // this is not assigning the sequence! This is adding a PIR file for a particular chain id!
   void assign_pir_sequence(const std::string &chain_id, const std::string &seq);

   void assign_sequence(const clipper::Xmap<float> &xmap, const std::string &chain_id);
   std::vector<std::pair<std::string, std::string> > sequence_info() const { return input_sequence; };

   // this does an alignment! How confusing
   void assign_sequence_from_file(const std::string &filename);

   // Apply to NCS-related chains too, if present
   void assign_sequence_to_NCS_related_chains_from_string(const std::string &chain_id, const std::string &seq);
   void assign_sequence_from_string_simple(const std::string &chain_id, const std::string &seq);

   void delete_all_sequences_from_molecule();

   void delete_sequence_by_chain_id(const std::string &chain_id);

   // render option (other functions)
   coot::ray_trace_molecule_info fill_raster_model_info(bool against_a_dark_background); // messes with bond_colour_internal
   coot::ray_trace_molecule_info fill_raster_map_info(short int lev) const;
   coot::ray_trace_molecule_info fill_raster_additional_info() const;

   // return a list of bad chiral volumes for this molecule
   // (first is a vector of bad chiral volume types (residues for which we don't have
   //  a dictionary).
   //
   std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> > bad_chiral_volumes() const;
   std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> > inverted_chiral_volumes() const;
   std::pair<std::vector<std::string>, std::vector<std::pair<coot::atom_spec_t, double> > > distorted_chiral_volumes(double chiral_volume_distortion_limit) const;

   // a other function
   float score_residue_range_fit_to_map(int res1, int res2, std::string altloc,
					std::string chain_id, int imol_for_map);
   // likewise:
   void fit_residue_range_to_map_by_simplex(int res1, int res2, std::string altloc,
					    std::string chain_id, int imol_for_map);

   // residue splitting (add alt conf)
   std::pair<bool,std::string> split_residue(int atom_index, int alt_conf_split_type);

   // internal (private) functions:
   std::pair<bool,std::string>
   split_residue_internal(mmdb::Residue *residue,
			  const std::string &altconf,
			  const std::vector<std::string> &all_altconfs,
			  atom_selection_container_t residue_mol,
			  short int use_residue_mol_flag);
   void split_residue_then_rotamer(mmdb::Residue *residue, const std::string &altconf,
				   const std::vector<std::string> &all_altconfs,
				   atom_selection_container_t residue_mol,
				   short int use_residue_mol_flag);
   void split_residue_internal(mmdb::Residue *residue);
   // We just added a new atom to a residue, now we need to adjust the
   // occupancy of the other atoms (so that we don't get residues with
   // atoms whose occupancy is greater than 1.0 (Care for SHELX molecule?)).
   void adjust_occupancy_other_residue_atoms(mmdb::Atom *at, mmdb::Residue *residue,
					     short int force_one_sum_flag);
   std::vector<std::string> get_residue_alt_confs(mmdb::Residue *res) const;
   std::string make_new_alt_conf(const std::vector<std::string> &residue_alt_confs,
				 int alt_conf_split_type_in) const;
   short int have_atoms_for_rotamer(mmdb::Residue *res) const;
   mmdb::Manager * create_mmdbmanager_from_res_selection(mmdb::PResidue *SelResidues,
							int nSelResidues,
							int have_flanking_residue_at_start,
							int have_flanking_residue_at_end,
							const std::string &altconf,
							const std::string &chain_id_1,
							short int residue_from_alt_conf_split_flag);

   // merge molecules

   std::pair<int, std::vector<merge_molecule_results_info_t> > merge_molecules(const std::vector<atom_selection_container_t> &add_molecules);
   std::pair<bool, std::vector<std::string> > try_add_by_consolidation(mmdb::Manager *adding_mol);
   bool merge_molecules_just_one_residue_homogeneous(atom_selection_container_t molecule_to_add);
   // try to add the ligand at the given spec, if not (say the spec was not filled or there
   // was already a ligand at the given spec) then return false.
   bool merge_molecules_just_one_residue_at_given_spec(atom_selection_container_t molecule_to_add,
						       coot::residue_spec_t target_spec);
   std::pair<bool, coot::residue_spec_t> merge_ligand_to_near_chain(mmdb::Manager *mol); // return success status and spec if new residue if possible.

   // merge change/fragments of this molecule
   // return 1 if a merge was done;
   int merge_fragments();

   bool is_intermediate_atoms_molecule;

   int renumber_residue_range(const std::string &chain_id,
			      int start_resno, int last_resno, int offset);

   int renumber_residue_range_old(const std::string &chain_id,
			      int start_resno, int last_resno, int offset);

   int change_residue_number(const std::string &chain_id_str,
			     int current_resno,
			     const std::string &current_inscode_str,
			     int new_resno,
			     const std::string &new_inscode_str);

   int renumber_waters(); // renumber all solvent changes so that
			  // their waters start at residue 1 and
			  // continue monotonically.

   void split_water(std::string chain_id, int res_no, std::string ins_code,
		    const clipper::Xmap<float> &xmap,
		    float sigma);

   void mark_atom_as_fixed(const coot::atom_spec_t &atom_spec, bool state);

   // validation
   void find_deviant_geometry(float strictness);


   // ===================== NCS ghosts ===============================

   void set_show_ghosts(short int istate);
   void set_ghost_bond_thickness(float f);
   int update_ncs_ghosts();
   float ghost_bond_thickness() const {return int(ghost_bond_width);}
   int draw_ncs_ghosts_p() const { // needed for setting the Bond Parameters checkbutton
      return show_ghosts_flag;
   }

   void draw_ncs_ghosts(Shader *shader_for_meshes,
                        stereo_eye_t eye,
                        const glm::mat4 &mvp,
                        const glm::mat4 &model_rotation_matrix,
                        const std::map<unsigned int, lights_info_t> &lights,
                        const glm::vec3 &eye_position,
                        const glm::vec4 &background_colour);

   std::vector<drawn_ghost_molecule_display_t> NCS_ghosts() const;

   std::vector<std::vector<std::string> > ncs_ghost_chains() const;

   // Not const because we may modify ncs_ghosts by adding their ncs
   // operators: (and recall that this function is a coordinates
   // molecule function, so it uses its ncs operators on the given
   // map) to make new maps):
   // BL note:: made extra arg a string, so that we cannot interfere
   // with imol_map
   std::vector<std::pair<clipper::Xmap<float>, std::string> >
     ncs_averaged_maps(const clipper::Xmap<float> &xmap_in, float homology_lev, std::string &imol_map_name);
   short int has_ncs_p() { return (ncs_ghosts.size() > 0) ? 1 : 0; }
   short int ncs_ghosts_have_rtops_p() const { return ncs_ghosts_have_rtops_flag;}
   int fill_ghost_info(short int do_rtops_flag,
		       float homology_lev);   // button callback fills ghosts
					      // if requested to display them.
   void add_ncs_ghost(const std::string &chain_id,
		      const std::string &target_chain_id,
		      const clipper::RTop_orth &rtop);

   void clear_ncs_ghost_matrices();

   // and the other way for CNS NCS users
   void add_strict_ncs_matrix(const std::string &chain_id,
			      const std::string &target_chain_id,
			      const coot::coot_mat44 &m);

   void add_molecular_symmetry(const clipper::Mat33<double> &mol_symm,
                               const clipper::Coord_orth &molecular_origin);

   // and that add to this:
   // (consider using a class)
   std::vector<std::pair<clipper::Mat33<double>, clipper::Coord_orth> > molecular_symmetry_matrices;

   // trivial helper class for add_molecular_symmetry_matrices()
   class quad_d_t {
   public:
      double x, y, z, t;
      quad_d_t(const double &x_in, const double &y_in, const double &z_in, const double &t_in) {
	 x = x_in; y = y_in; z = z_in; t = t_in;
      }
      quad_d_t() {}
   };
   void add_molecular_symmetry_matrices(); // process REMARK 350s

   // the first value is if we should apply the matrix or not (we may not have ghosts)
   std::pair<bool, clipper::RTop_orth>
     apply_ncs_to_view_orientation(const clipper::Mat33<double> &current_view_mat,
				   const clipper::Coord_orth &current_position,
				   const std::string &current_chain,
				   const std::string &next_ncs_chain,
				   bool backward_flag) const;

   short int show_strict_ncs_flag;

   void add_strict_ncs_from_mtrix_from_file(const std::string &file_name);
   void add_strict_ncs_from_mtrix_from_self_file();

   // New style EM molecular symmetry
   void add_molecular_symmetry_from_mtrix_from_self_file();
   void add_molecular_symmetry_from_mtrix_from_file(const std::string &file_name);

   // Not 'const' because we can do a fill_ghost_info if the NCS ghosts
   // do not have rtops.
   coot::ncs_differences_t ncs_chain_differences(std::string master_chain_id,
						 float main_chain_weight);

   // if we are are the centre of a given chain_id, how big a radius
   // do we need to encompass all atoms of that chain?
   std::pair<clipper::Coord_orth, double> chain_centre_and_radius(const std::string &chain_id) const;




   // ====================== SHELX stuff ======================================
   std::pair<int, std::string> write_shelx_ins_file(const std::string &filename);
   int read_shelx_ins_file(const std::string &filename);
   // return the success status, 0 for fail, 1 for good.
   int add_shelx_string_to_molecule(const std::string &str);
   bool is_from_shelx_ins() const { return is_from_shelx_ins_flag; }

   // data resolution, in A (or a negative number on error)
   float data_resolution() const { return data_resolution_; }

   // Change chain id
   // return -1 on a conflict
   //
   std::pair<int, std::string> change_chain_id(const std::string &from_chain_id,
					       const std::string &to_chain_id,
					       bool use_resno_range,
					       int start_resno,
					       int end_resno);

   // return 1 on successful deletion, 0 on fail
   int delete_zone(const coot::residue_spec_t &res1, const coot::residue_spec_t &res2);

   //
   void spin_search(clipper::Xmap<float> &xmap,
		    const std::string &chain_id,
		    int resno,
		    const std::string &ins_code,
		    const std::pair<std::string, std::string> &direction_atoms,
		    const std::vector<std::string> &moving_atoms_list);


    density_results_container_t
    spin_atom(const clipper::Xmap<float> &xmap,
              const coot::residue_spec_t &spec,
              const std::string &direction_atoms_ref,  //e.g. N
              const std::string &direction_atoms_base, //e.g. CA
              const std::string &direction_atoms_tip,  //e.g. CB, where moving atom is CG
              const std::vector<std::string> &moving_atoms_list) const;


   std::vector<std::pair<coot::residue_spec_t, float> >
   em_ringer(const clipper::Xmap<float> &xmap) const;

   // nomenclature errors
   // return a vector of the changed residues (used for updating the rotamer graph)
   std::vector<mmdb::Residue *> fix_nomenclature_errors(coot::protein_geometry *geom_p);
                                                      // by looking for bad rotamers in
				                      // some residue types and alter ing
                            			      // the atom names to see if they get
				                      // more likely rotamers
   // the residue type and the spec
   std::vector<std::pair<std::string, coot::residue_spec_t> > list_nomenclature_errors(const coot::protein_geometry *geom_p);

   // ---- cis <-> trans conversion
   int cis_trans_conversion(const std::string &chain_id, int resno, const std::string &inscode,
			    mmdb::Manager *standard_residues_mol);
   int cis_trans_conversion(mmdb::Atom *at, short int is_N_flag, mmdb::Manager *standard_residues_mol);
   int cis_trans_convert(mmdb::PResidue *mol_residues,   // internal function, make private
			 mmdb::PResidue *trans_residues, // or move into utils?
			 mmdb::PResidue *cis_residues);


   // ---- baton build redirection ----
   short int reverse_direction_of_fragment(const std::string &chain_id,
					   int resno);

   // ---- missing atoms ----
   //
   // Return a vector of residues that have missing atoms by dictionary
   // search.  missing_hydrogens_flag reflects if we want to count
   // residues that have missing hydrogens as residues with missing
   // atoms that should be part of the returned vector. Most of the
   // time, we don't care about hydrogens and the flag is 0.
   //
   // Pass and potentially add to geom_p
   coot::util::missing_atom_info
   missing_atoms(short int missing_hydrogens_flag, coot::protein_geometry *geom_p) const;

   // The function that uses missing atom info:
   coot::util::missing_atom_info
   fill_partial_residues(coot::protein_geometry *geom_p, int imol_refinement_map);
   // return 1 if the residue was filled, 0 if the residue was not found
   int fill_partial_residue(const coot::residue_spec_t &residue_spec,
                            const coot::protein_geometry *geom_p, int imol_refinement_map);

   std::vector<std::string> get_chain_ids() const;

   // Ribosome People:
   int exchange_chain_ids_for_seg_ids();

   // public interface to chain copying
   void copy_chain(const std::string &from_chain, const std::string &to_chain);
   int copy_from_ncs_master_to_others(const std::string &master_chain_id);
   int copy_from_ncs_master_to_specific_other_chains(const std::string &master_chain_id,
                                                     const std::vector<std::string> &other_chain_ids);
   int copy_residue_range_from_ncs_master_to_other_using_ghost(std::string from_chain_id,
							       std::string to_chain_id,
							       int residue_range_1,
							       int residue_range_2);
   int copy_residue_range_from_ncs_master_to_others(const std::string &master_chain_id,
						    int resno_start, int resno_end);
   int copy_residue_range_from_ncs_master_to_chains(const std::string &master_chain_id,
						    int resno_start, int resno_end,
						    const std::vector<std::string> &chain_ids);
   int copy_from_ncs_master_to_chains(const std::string &master_chain_id,
				      const std::vector<std::string> &chain_ids);
   int set_ncs_master_chain(const std::string &new_master_chain_id, float homology_lev);
   // add ghosts
   void add_ncs_ghosts_no_explicit_master(const std::vector<std::string> &chain_ids,
					  const std::vector<std::vector<std::pair<std::string, int> > > &residue_types,
					  std::vector<short int> first_chain_of_this_type,
					  const std::vector<int> &chain_atom_selection_handles,
					  short int do_rtops_flag,
					  float homology_lev,
					  bool allow_offset_flag);
   void add_ncs_ghosts_using_ncs_master(const std::string &master_chain_id,
					const std::vector<std::string> &chain_ids,
					const std::vector<std::vector<std::pair<std::string, int> > > &residue_types,
					const std::vector<int> &chain_atom_selection_handles,
					float homology_lev);

   // symmetry control
   //
   // We fill the frame that's passed.  It is used to fill the
   // symmetry control widget (requested by Frank von Delft)

   void fill_symmetry_control_frame(GtkWidget *dialog) const;

   int   symmetry_whole_chain_flag;
   int   symmetry_as_calphas;
   int   symmetry_colour_by_symop_flag;
   short int symmetry_rotate_colour_map_flag; // do we want symmetry of other
						     // molecules to have a different
						     // colour [MOL]?

   // ncs control

   void move_reference_chain_to_symm_chain_position(coot::Symm_Atom_Pick_Info_t naii);
   void fill_ncs_control_frame(GtkWidget *dialog) const; // called for every coords mol
   void fill_ncs_control_frame_internal(GtkWidget *dialog) const; // called if needed.
   void old_fill_ncs_control_frame_internal(GtkWidget *dialog) const; // delete one day
   void ncs_control_change_ncs_master_to_chain_update_widget(GtkWidget *w, int ichain) const;

   void set_display_ncs_ghost_chain(int ichain, int state);
   // return status 0 if ncs master chain was not set.
   std::pair<bool, std::string> first_ncs_master_chain_id() const; // for ncs graphs use
   std::vector<std::string> ncs_master_chains() const;


   std::vector<std::string> get_symop_strings() const;


   // Replace the atoms in this molecule by those in the given atom selection.
   int replace_fragment(atom_selection_container_t asc);

   int swap_atom_alt_conf(std::string chain_id, int res_no, std::string ins_code, std::string atom_name,
                          std::string alt_conf);

   std::vector<std::string> alt_confs_in_molecule() const;

   int swap_residue_alt_confs(const std::string &chain_id, int res_no, const std::string &ins_code);

   int set_atom_attribute(std::string chain_id, int resno, std::string ins_code,
			  std::string atom_name, std::string alt_conf,
			  std::string attribute_name, float val);

   int set_atom_string_attribute(std::string chain_id, int resno, std::string ins_code,
				 std::string atom_name, std::string alt_conf,
				 std::string attribute_name, std::string val_str);

   int set_atom_attributes(const std::vector<coot::atom_attribute_setting_t> &v);

   void set_residue_name(std::string chain_id, int res_no, std::string ins_code, std::string new_name);


   coot::at_dist_info_t closest_atom(const coot::Cartesian &pt,
				     bool ca_check_flag) const;
   coot::at_dist_info_t closest_atom(const coot::Cartesian &pt,
				     bool ca_check_flag,
				     const std::string &chain_id,
				     bool use_this_chain_id) const;
   coot::at_dist_info_t closest_atom(const coot::Cartesian &pt) const;


   // cleaner interface to molecule's attributes:
   std::pair<bool, clipper::Spacegroup> space_group() const;
   std::pair<bool, clipper::Cell> cell() const;


   //
   clipper::Coord_orth find_peak_along_line(const clipper::Coord_orth &p1,
					    const clipper::Coord_orth &p2) const;
   //  Throw an exception if peak not found (about the contour level for this map)
   clipper::Coord_orth find_peak_along_line_favour_front(const clipper::Coord_orth &p1,
							 const clipper::Coord_orth &p2) const;

   coot::minimol::molecule eigen_flip_residue(const std::string &chain_id, int resno);

   // return value is an error string that we can put in the status bar
   std::string jed_flip(coot::residue_spec_t &spec, const std::string &atom_name, const std::string &alt_conf,
			bool invert_selection,
			coot::protein_geometry *geom);

   // replace molecule
   int replace_molecule(mmdb::Manager *mol);
   int replace_models(std::deque<mmdb::Model *> model_list);


   // add a factor to scale the colours in b factor representation:.
   // It goes into the atom_sel.mol
   void set_b_factor_bonds_scale_factor(float f);

   int
   apply_sequence(int imol_map, mmdb::Manager *mol,
		  std::vector<coot::residue_spec_t> mmdb_residues,
		  std::string best_seq, std::string chain_id,
		  int resno_offset,
		  const coot::protein_geometry &pg);

   int delete_all_except_res(mmdb::Residue *res);

   // EM map function
   int scale_cell(float fac_u, float fac_v, float fac_w);

   // Shall we try to use gompertz correction factor (based on f/sigf
   // (should be I/sigI, ideally).  Try to use gompertz if it is
   // available.  Gompertz factor is by default set so that an f/sigf
   // of 3 will result in 50% correction.
   //
   void sharpen(float b_factor, bool gompertz, float gompertz_factor);

   // number of chains. Return -1 on failure
   int number_of_chains() const;

   // reorder the chains in the models
   void sort_chains();

   // reorder the residues in the models
   void sort_residues();

   int add_additional_representation(int representation_type,
				     const int &bonds_box_type_in,
				     float bonds_width,
				     bool draw_hydrogens_flag,
				     const coot::atom_selection_info_t &info,
				     GtkWidget *display_control_window,
				     const gl_context_info_t &glci,
				     const coot::protein_geometry *geom);

   int adjust_additional_representation(int representation_number,
					const int &bonds_box_type_in,
					float bonds_width,
					bool draw_hydrogens_flag,
					const coot::atom_selection_info_t &info,
					bool show_it_flag_in);

   void clear_additional_representation(int representation_number);
   void set_show_additional_representation(int representation_number, bool on_off_flag);
   void set_show_all_additional_representations(bool on_off_flag);
   void all_additional_representations_off_except(int rep_no,
						  bool ball_and_sticks_off_too_flag);
   graphical_bonds_container get_bonds_representation() { make_bonds_type_checked(); return bonds_box; }
   //
   std::vector<coot::residue_spec_t> residues_near_residue(const coot::residue_spec_t &rspec, float radius) const;
   void label_closest_atoms_in_neighbour_atoms(coot::residue_spec_t residue_spec, float radius);

   void remove_ter_atoms(const coot::residue_spec_t &spec); // from all models

   // c.f. progressive_residues_in_chain_check_by_chain()
   // bool residues_in_order_p(std::string &chain_id) const;


   // Only apply charges if the molecule contains lots of hydrogens or
   // there were few (< 100) atoms in the molecule.
   //
   // so return a flag, whether or not the charges were applied.
   //
   // More than 15% of the atoms have to be hydrogens for use to set
   // the charges on all the atoms (to something other than
   // CXX_UNSET_CHARGE).
   //
   bool apply_charges(const coot::protein_geometry &geom);

   // return the dipole (a copy) and its number
   //
   std::pair<coot::dipole, int>
   add_dipole(const std::vector<coot::residue_spec_t> &res_specs,
	      const coot::protein_geometry &geom);

   void delete_dipole(int dipole_number);

   // return the number of new hetatoms
   int assign_hetatms();
   bool is_het_residue(mmdb::Residue *residue_p) const;

   // move waters so that they are around H-bonders (non-C) in protein.
   // return the number of moved atoms
   int move_waters_to_around_protein();

   void move_hetgroups_to_around_protein();


   // Return the maximum minimum distance of waters to protein atoms.
   // return something negative when we can't do above (no protein
   // atoms or no water atoms).
   float max_water_distance();

   float fit_chain_to_map_by_random_jiggle(const std::string &chain_id, const clipper::Xmap<float> &xmap,
                                           float map_sigma,
                                           int n_trials, float jiggle_scale_factor);
   float fit_molecule_to_map_by_random_jiggle(const clipper::Xmap<float> &xmap,
                                              float map_sigma, int n_trias, float jiggle_scale_factor);

   // jiggle residue (a specific, useful/typical interface to jiggling).
   float fit_to_map_by_random_jiggle(coot::residue_spec_t &spec,
				     const clipper::Xmap<float> &xmap,
				     float map_sigma,
				     int n_trials,
				     float jiggle_scale_factor);

   // Random rotation and translation (translations scaled by
   // jiggle_scale_factor).
   //
   // return the z-weighted fit to density score of the atom
   // selection.
   //
   // called by above
   //
   // if chain_for_moving is not empty, apply the transformation
   // the the atoms of chain_for_moving rather than to the atom of atom_selection
   //
   float fit_to_map_by_random_jiggle(mmdb::PPAtom atom_selection,
				     int n_atoms,
				     const clipper::Xmap<float> &xmap,
				     float map_sigma,
				     int n_trials,
				     float jiggle_scale_factor,
				     bool use_biased_density_scoring,
				     std::vector<mmdb::Chain *> chains_for_moving);

#ifdef HAVE_CXX_THREAD
   static void test_jiggle_fit_func(unsigned int thread_index,
				    unsigned int i_trial,
				    unsigned int n_trials,
				    mmdb::PPAtom atom_selection,
				    int n_atoms,
				    const std::vector<mmdb::Atom *> &initial_atoms,
				    const clipper::Coord_orth &centre_pt,
				    const std::vector<std::pair<std::string, int> > &atom_numbers,
				    const clipper::Xmap<float> *xmap_masked,
				    float jiggle_scale_factor);
   static void jiggle_fit_multi_thread_func_1(int thread_index,
					      unsigned int i_trial,
					      unsigned int n_trials,
					      mmdb::PPAtom atom_selection,
					      int n_atoms,
					      const std::vector<mmdb::Atom *> &initial_atoms,
					      const clipper::Coord_orth &centre_pt,
					      float jiggle_scale_factor,
					      const std::vector<std::pair<std::string, int> > &atom_numbers,
					      const clipper::Xmap<float> *xmap_masked_p,
					      float (*density_scoring_function)(const coot::minimol::molecule &mol,
										const std::vector<std::pair<std::string, int> > &atom_number_list,
										const clipper::Xmap<float> &map),
					      std::pair<clipper::RTop_orth, float> *trail_results_p);
   static void jiggle_fit_multi_thread_func_2(int thread_index,
					      const coot::minimol::molecule &direct_mol,
					      const clipper::Xmap<float> &xmap_masked,
					      float map_sigma,
					      const clipper::Coord_orth &centre_pt,
					      const std::vector<std::pair<std::string, int> > &atom_numbers,
					      float trial_results_pre_fit_score_for_trial,
					      float (*density_scoring_function)(const coot::minimol::molecule &mol,
										const std::vector<std::pair<std::string, int> > &atom_number_list,
										const clipper::Xmap<float> &map),
					      std::pair<clipper::RTop_orth, float> *post_fix_scores_p);
#endif

   // return a fitted molecule
   coot::minimol::molecule rigid_body_fit(const coot::minimol::molecule &mol_in,
					  const clipper::Xmap<float> &xmap,
					  float map_sigma) const;

   // ---- utility function --- (so that we know to delete hydrogens
   // from HETATM molecule before merging with this one
   //
   bool molecule_has_hydrogens() const;

   bool get_input_molecule_was_in_mmcif_state() const {
      return input_molecule_was_in_mmcif;
   }


   // -------- simply print it (at the moment) --------------
   void print_secondary_structure_info();

   // --------- pisa surface make dots --------------------
   // here we add the dots to the dots vector of this molecule
   void add_dots(const pli::dots_representation_info_t &dots_in) {
      dots.push_back(dots_in);
   }

   // ---- colour dots surface -----------
   //
   void set_dots_colour(float r, float g, float b) {
      dots_colour.set(r,g,b);
      dots_colour_set = true;
   }

   // ---- extra restraints (currently only bonds) -----------
   //
   bool draw_it_for_extra_restraints;
   bool draw_it_for_parallel_plane_restraints;
   coot::extra_restraints_t extra_restraints;
   bool extra_restraints_representation_for_bonds_go_to_CA;
   void set_extra_restraints_representation_for_bonds_go_to_CA(bool val) {
      if (val != extra_restraints_representation_for_bonds_go_to_CA) {
	      extra_restraints_representation_for_bonds_go_to_CA = val;
	      update_extra_restraints_representation();
      }
   }
   coot::extra_restraints_representation_t extra_restraints_representation;
   void draw_extra_restraints_representation();
   void draw_parallel_plane_restraints_representation();
   void set_extra_restraints_prosmart_sigma_limits(double limit_low, double limit_high);

   // return an index of the new restraint
   int add_extra_bond_restraint(coot::atom_spec_t atom_1,
				coot::atom_spec_t atom_2,
				double bond_dist, double esd);
   int add_extra_geman_mcclure_restraint(coot::atom_spec_t atom_1,
                                         coot::atom_spec_t atom_2,
                                         double bond_dist, double esd);
   int add_extra_bond_restraints(const std::vector<coot::extra_restraints_t::extra_bond_restraint_t> &bond_specs);
   int add_extra_geman_mcclure_restraints(const std::vector<coot::extra_restraints_t::extra_geman_mcclure_restraint_t> &bond_specs);
   int add_extra_angle_restraint(coot::atom_spec_t atom_1,
				 coot::atom_spec_t atom_2,
				 coot::atom_spec_t atom_3,
				 double angle, double esd);
   int add_extra_torsion_restraint(coot::atom_spec_t atom_1,
				   coot::atom_spec_t atom_2,
				   coot::atom_spec_t atom_3,
				   coot::atom_spec_t atom_4,
				   double torsion_angle, double esd, int period);
   int add_extra_start_pos_restraint(coot::atom_spec_t atom_1,
				     double esd);

   // extra target position restraints are like pull atom restraints
   int add_extra_target_position_restraint(coot::atom_spec_t &spec,
					   const clipper::Coord_orth &pos,
					   float weight);
   int add_extra_target_position_restraints(const std::vector<std::tuple<coot::atom_spec_t, const clipper::Coord_orth , float > > &etprs);

   // the atom specs do not need to be in order for bond restraints only
   void remove_extra_bond_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2);
   void remove_extra_geman_mcclure_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2);
   void remove_extra_start_pos_restraint(coot::atom_spec_t atom_1);
   void remove_extra_angle_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2,
                                      coot::atom_spec_t atom_3);
   void remove_extra_torsion_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2,
                                      coot::atom_spec_t atom_3, coot::atom_spec_t atom_4);
   void update_extra_restraints_representation(); // called from make_bonds_type_checked()
   void update_extra_restraints_representation_bonds();
   void update_extra_restraints_representation_geman_mcclure();
   void update_extra_restraints_representation_bonds_internal(const coot::extra_restraints_t::extra_bond_restraint_t &res);
   void update_extra_restraints_representation_parallel_planes();
   void add_refmac_extra_restraints(const std::string &file_name);
   void remove_extra_target_position_restraints(coot::atom_spec_t &spec);

   // make them yourself - easy as pie.
   void generate_self_restraints(float local_dist_max,
				 const coot::protein_geometry &geom);
   void generate_local_self_restraints(float local_dist_max,
				       const std::string &chain_id,
				       const coot::protein_geometry &geom);
   void generate_local_self_restraints(float local_dist_max,
				       const std::vector<coot::residue_spec_t> &residue_specs,
				       const coot::protein_geometry &geom);
   void generate_local_self_restraints(int selHnd, float local_dist_max,
				       const coot::protein_geometry &geom);

   void add_parallel_plane_restraint(coot::residue_spec_t spec_1,
				     coot::residue_spec_t spec_2);
   // which uses:
   std::vector<std::string> nucelotide_residue_name_to_base_atom_names(const std::string &rn) const;
   // for non-bases, normal amino acids (simple-minded, currently).
   std::vector<std::string> residue_name_to_plane_atom_names(const std::string &rn) const;

   void clear_extra_restraints();

   // --------- (transparent) solid rendering of density ------------------
   // bool draw_it_for_solid_density_surface;  // everything is "solid" now (calculate normals)

   // 20200312-PE Now we use draw_vector_sets.
   // coot::density_contour_triangles_container_t tri_con;


   coot::density_contour_triangles_container_t tri_con_diff_map_neg; // negative contour

   // these functions are old and need some thought to see if they need to exist.
   void display_solid_surface_triangles(const coot::density_contour_triangles_container_t &tri_con, bool do_flat_shading) const;
   void draw_solid_density_surface(bool do_flat_shading);
   void set_draw_solid_density_surface(bool state);

   // new
   void post_process_map_triangles();
   void setup_glsl_map_rendering(const clipper::Coord_orth &centre, float radius);
   std::pair<std::vector<coot::api::vertex_with_rotation_translation>, std::vector<g_triangle> >
   make_generic_vertices_for_atoms(const std::vector<glm::vec4> &index_to_colour, float atom_radius_scale_factor=1.0) const;
   std::pair<std::vector<coot::api::vertex_with_rotation_translation>, std::vector<g_triangle> >
   make_generic_vertices_for_rama_balls(float ball_scale_factor, const glm::vec3 &screen_up_dir) const;
   std::pair<std::vector<coot::api::vertex_with_rotation_translation>, std::vector<g_triangle> >
   make_generic_vertices_for_bad_CA_CA_distances() const;
   std::pair<std::vector<coot::api::vertex_with_rotation_translation>, std::vector<g_triangle> > make_end_cap(float z);
   std::pair<std::vector<coot::api::vertex_with_rotation_translation>, std::vector<g_triangle> > fun(float radius_scale) const;

   void setup_glsl_bonds_buffers(const std::vector<coot::api::vertex_with_rotation_translation> &vertices,
                                 const std::vector<g_triangle> &triangles);

   GLuint m_VertexArrayID_for_map;
   GLuint m_VertexArrayID_for_map_cap;

   GLuint n_vertices_for_map_VertexArray;

   GLuint n_indices_for_triangles;
   GLuint n_indices_for_lines;
   GLuint m_VertexBufferID;
   GLuint m_IndexBuffer_for_map_lines_ID;
   GLuint m_IndexBuffer_for_map_triangles_ID; // solid and transparent surfaces

   // 20220211-PE pre map-as-mesh rewrite.
   // GLuint m_NormalBufferID; // is this map or model - or something else? Be clear!
   // GLuint m_ColourBufferID; // Likewise.

   Mesh map_as_mesh;
   Mesh map_as_mesh_gl_lines_version;

   GLuint m_VertexArray_for_model_ID;
   GLuint n_vertices_for_model_VertexArray;
   GLuint n_indices_for_model_triangles;
   GLuint m_VertexBuffer_for_model_ID;
   GLuint m_IndexBuffer_for_model_ID;
   GLuint m_NormalBuffer_for_model_ID;
   GLuint m_ColourBuffer_for_modelID;
   GLuint m_ModelMatrix_for_model_ID;
   GLuint m_VertexBuffer_for_map_cap_ID;  // more on map cap below
   GLuint m_IndexBuffer_for_map_cap_ID;
   GLuint n_vertices_for_map_cap;

   bool map_mesh_first_time;
   bool model_mesh_first_time;

   float density_surface_opacity;
   bool is_an_opaque_map() const { return density_surface_opacity == 1.0; } // needs explicit assignment to 1.0
                                                                            // elsewhere in the code, e.g. in
                                                                            // the adjustment handler.

   Material material_for_maps;
   Material material_for_models;

void draw_map_molecule(stereo_eye_t eye,
                          bool draw_transparent_maps,
                          Shader &shader, // unusual reference.. .change to pointer for consistency?
                          const glm::mat4 &mvp,
                          const glm::mat4 &view_rotation,
                          const glm::vec3 &eye_position,
                          const glm::vec4 &ep,
                          const std::map<unsigned int, lights_info_t> &lights,
                          const glm::vec3 &background_colour,
                          bool perspective_projection_flag);
   // A map is not a Mesh at the moment, so this needs a new function
   void draw_map_molecule_for_ssao(Shader *shader_p, const glm::mat4 &model_matrix, const glm::mat4 &view_matrix, const glm::mat4 &proj_matrix);

   // using current contour level,
   // return world coordinates and normals
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
   make_map_cap(const clipper::Coord_orth &base_point,
                const clipper::Coord_orth &x_axis_uv,  // unit vectors
                const clipper::Coord_orth &y_axis_uv,
                double x_axis_step_size,
                double y_axis_step_size,
                unsigned int n_axis_points,
                unsigned int y_axis_points) const;

   void setup_map_cap(Shader *shader_p,
                      const clipper::Coord_orth &base_pt, // Bring it into this class.
                      const clipper::Coord_orth &x_axis_uv, // Of the cap plane, of course.
                      const clipper::Coord_orth &y_axis_uv,
                      double x_axis_step_size,
                      double y_axis_step_size,
                      unsigned int n_axis_points,
                      unsigned int y_axis_points);
   void draw_map_cap(Shader *shader_p,
                     const glm::mat4 &mvp,
                     const glm::mat4 &world_rotation_matrix,
                     const glm::mat4 &world_rotation_translation_matrix,
                     const std::map<unsigned int, lights_info_t> &lights,
                     const glm::vec3 &eye_position);

   // non molecular mesh
   Shader shader_for_draw_map_normals;
   void draw_map_normals(const glm::mat4 &mvp);

   // uses molecular meshes/graphical molecules
   void draw_normals(const glm::mat4 &mvp); // defunct now I think
   void mesh_draw_normals(const glm::mat4 &mvp);

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > make_map_mesh();

   float shader_shininess;
   float shader_specular_strength;

   int get_square_type(const unsigned int &ii, const unsigned int &jj, // contouring
                       const coord_array_2d &arr, const float &contour_level) const;

   void setup_density_surface_material(bool solid_mode, float opacity,
				       bool is_negative_level = 0); // shininess, material colour etc.
   bool transparent_molecular_surface_flag; // 0 by default.

   // return a error status and an error message.
   std::pair<bool, std::string>
   sprout_hydrogens(const std::string &chain_id,
		    int res_no,
		    const std::string &ins_code,
		    const coot::protein_geometry &geom);
   // which uses:
   // (and this returns a status: 0 if we did not correct a chiral, 1 if we did)
   bool sprout_hydrogens_correct_chirals_maybe(mmdb::Residue *residue_cp_p,
					       const std::string &alt_conf,
					       const coot::dictionary_residue_restraints_t &rp);
   void sprout_hydrogens_transfer_hydrogen_positions(mmdb::Residue *from_res, mmdb::Residue *to_res,
						     const std::string &alt_conf);



   std::vector<std::string> no_dictionary_for_residue_type_as_yet(const coot::protein_geometry &geom) const;

   // ------------- helper function to orient_view() --------------------
   // the vector is from the central residue (atom) to the neighbour_residue (atom);
   //
   // can throw an std::runtime  exception;
   //
   clipper::Coord_orth get_vector(const coot::residue_spec_t &central_residue_spec, // ligand typically
				  const coot::residue_spec_t &neighbour_residue_spec) const;

   // --------- match ligand atom names ------------------
   void match_ligand_atom_names(const std::string &chain_id, int res_no, const std::string &ins_code,
				mmdb::Residue *res_ref);

   // --------- molecule probability scoring ------------
   coot::rama_score_t get_all_molecule_rama_score() const;
   coot::rama_score_t get_all_molecule_rama_score_old() const; // 20230611-PE delete this one day
   coot::rotamer_score_t get_all_molecule_rotamer_score(const coot::rotamer_probability_tables &rpt) const;

   // --------- lsq-improve ------------
   void lsq_improve(mmdb::Manager *mol_ref, const std::string &ref_selection_str,
		    const std::string &moving_selection_str, int n_res, float dist_crit);

   std::vector<ProteinDB::Chain> protein_db_loops(const std::vector<coot::residue_spec_t> &residue_specs,
						  int nfrags, const clipper::Xmap<float> &xmap);
   ProteinDB::Chain make_fragment_chain(const std::vector<coot::residue_spec_t> &residue_specs) const;

   // ---------- LSQ fit by atom pairs -------------
   // return success (0 is fail)
   int lsq_fit_by_atom_pairs(const std::vector<coot::atom_spec_t> &mov_atom_specs,
			     const std::vector<clipper::Coord_orth> &ref_coords);


   // --------- HETATMs ------------
   int residue_has_hetatms(const std::string &chain_id, int resno, const std::string &ins_code) const;

   // change them if they are not standard PDB 'ATOM' residues
   int hetify_residue_atoms(const std::string &chain_id, int resno, const std::string &ins_code);

   // so that we only draw things that need redrawing on modifying ligands
   bool has_residue_with_name(const std::string comp_id) const; // unimplemented

   // --------- LINKs ---------------
   void make_link(const coot::atom_spec_t &spec_1, const coot::atom_spec_t &spec_2,
		  const std::string &link_name, float length,
		  const coot::protein_geometry &geom);
   void delete_any_link_containing_residue(const coot::residue_spec_t &res_spec);
   // this will not do a update of bonds, caller should do that.
   void update_any_link_containing_residue(const coot::residue_spec_t &old_spec,
					   const coot::residue_spec_t &new_spec);
   void delete_link(mmdb::Link *link, mmdb::Model *model_p);


   // ------------ watson crick pair additions  ---------
   int watson_crick_pair_for_residue_range(const std::string &chain_id,
					   int resno_start, int resno_end,
					   mmdb::Manager *standard_residues_mol);

   // --------- Pretty (hopefully) animated ligand interactions -----------

   std::vector<coot::animated_ligand_interactions_t> animated_ligand_interactions_vec;
   void add_animated_ligand_interaction(const  pli::fle_ligand_bond_t &lb);

   void draw_animated_ligand_interactions(const gl_context_info_t &gl,
					  const long &start_time) const;

   bool draw_animated_ligand_interactions_flag; // tweaked by outside function
   void add_hydrogens_from_file(const std::string &reduce_pdb_out);

   void add_hydrogen_atoms_to_residue(const coot::residue_spec_t &rs);

   std::pair<bool, clipper::Coord_orth>
   residue_centre(const coot::residue_spec_t &spec) const;
   std::pair<bool, clipper::Coord_orth>
   residue_centre(const std::string &chain_id, int resno, const std::string &ins_code) const;
   // which calls:
   std::pair<bool, clipper::Coord_orth> residue_centre(mmdb::Residue *residue_p) const;
   // related (distance between residue_centres), return negative number if not valid:
   float distance_between_residues(mmdb::Residue *r1, mmdb::Residue *r2) const;

   // ------------------- ligand centre ---------------------
   // we want a button that goes to the ligand when we click it.
   // (if we are already on a ligand, go to the next one).
   //
   // the first value flags if this contains a useful return value.
   //
   // the first value flags if this contains a useful return value.
   // 1: normal case, go ahead and use the coord
   // 0:  No ligands found
   // -1: No movement because we are at the (single) ligand already.
   //
   coot::new_centre_info_t
   new_ligand_centre(const clipper::Coord_orth &current_centre, int n_atoms_min) const;

   // Add a LINK record if link_type is not blank (link_type is for
   // example "NAG-ASN")
   // return a pair, success-status and the added residue (it is a deep copy of the res_new)
   std::pair<bool, mmdb::Residue *> add_residue(mmdb::Residue *res_new, const std::string &chain_id);

   // return the number of added atoms
   int add_residue_with_atoms(const coot::residue_spec_t &residue_spec, const std::string &res_name, const std::vector<coot::minimol::atom> &list_of_atoms);

   // Add a LINK record if link_type is not blank (link_type is for
   // example "NAG-ASN")
   //
   //  (the protein_geometry is passes so that we // can look up the
   //  bonded atoms in it for the link).
   //
   // return the spec of the new residue (possibly unset).
   //
   coot::residue_spec_t
   add_linked_residue_by_beam_in(const coot::residue_spec_t &spec_in,
				 const std::string &new_residue_comp_id,
				 const std::string &link_type,
				 coot::protein_geometry *geom_p);

   // 20140429 As above, but using the atom-by-torsion template method
   coot::residue_spec_t
   add_linked_residue_by_atom_torsions(const coot::residue_spec_t &spec_in,
				       const std::string &new_residue_comp_id,
				       const std::string &link_type,
				       coot::protein_geometry *geom_p,
				       float b_factor_new_atoms);

   // n-models
   int n_models() const;

   // single model view
   void single_model_view_model_number(int imodel);
   int single_model_view_this_model_number() const;
   int single_model_view_next_model_number(); // changes the representation
   int single_model_view_prev_model_number(); //    ditto.

   // multi-residue torsion map fitting interface
   void multi_residue_torsion_fit(const std::vector<coot::residue_spec_t> &residue_specs,
				  const clipper::Xmap<float> &xmap,
				  int n_trials,
				  coot::protein_geometry *geom_p);


   // export map fragment (.ext)
   //
   void export_map_fragment(float radius,
			    clipper::Coord_orth centre,
			    const std::string &file_name) const;

   // shift "bottom left" to the origin and make sure that it's on a grid that is
   // "nice" (acceptable?) for molrep
   //
   int export_map_fragment_with_origin_shift(float radius,
					     clipper::Coord_orth centre,
					     const std::string &file_name) const;

   coot::residue_spec_t get_residue_by_type(const std::string &residue_type) const;

   std::vector<coot::residue_spec_t> get_residues_by_type(const std::string &residue_type) const;

   std::vector<coot::residue_spec_t> all_residues() const;

   std::vector<coot::residue_spec_t> het_groups() const;

   std::vector<mmdb::Residue *> get_all_protein_residues() const;

   // return null on failure.  seq_trip is something like "ACE".
   mmdb::Atom *get_centre_atom_from_sequence_triplet(const std::string &seq_trip) const;
   // which uses (like align's make_model_string).  Ignores waters.
   // The length of the string is guaranteed to the the length of the vector.
   std::pair<std::string, std::vector<mmdb::Residue *> > sequence_from_chain(mmdb::Chain *chain_p) const;

   std::string get_sequence_as_block(const std::string &chain_id) const;

   std::vector<coot::chain_mutation_info_container_t>
   sequence_comparison_to_chains(const std::string &sequence) const;

   void rotate_residue(const coot::residue_spec_t &rs,
		       const clipper::Coord_orth &around_vec,
		       const clipper::Coord_orth &origin_offset,
		       double angle);


   // for morphing
   class morph_rtop_triple {
   public:
      bool valid;
      clipper::Coord_orth co;
      clipper::RTop_orth rtop;
   public:
      morph_rtop_triple() { valid = false; }
      morph_rtop_triple(bool valid_in,
		  const clipper::Coord_orth &co_in,
		  const clipper::RTop_orth &rtop_in) {
	 valid = valid_in;
	 co = co_in;
	 rtop = rtop_in;
      }
      morph_rtop_triple(const clipper::Coord_orth &co_in,
		  const std::pair<bool, clipper::RTop_orth> &rtop_in) {
	 valid = rtop_in.first;
	 co = co_in;
	 rtop = rtop_in.second;
      }
   };


   // model morphing (average the atom shift by using shifts of the
   // atoms within shift_average_radius A of the central residue)
   //
   int morph_fit_all(const clipper::Xmap<float> &xmap_in, float shift_average_radius);
   int morph_fit_residues(std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > > moving_residues,
			  const clipper::Xmap<float> &xmap_in, float transformation_average_radius);
   int morph_fit_residues(const std::vector<coot::residue_spec_t> &residue_specs,
			  const clipper::Xmap<float> &xmap_in, float transformation_average_radius);
   int morph_fit_chain(const std::string &chain_id,
		       const clipper::Xmap<float> &xmap_in, float transformation_average_radius);
   void morph_show_shifts(const std::map<mmdb::Residue *, morph_rtop_triple> &simple_shifts,
			  const std::map<mmdb::Residue *, morph_rtop_triple> &smooth_shifts) const;

   crunch_model_t morph_fit_crunch_analysis(const std::map<mmdb::Residue *, morph_rtop_triple> &smooth_shifts) const;
   void morph_fit_uncrunch(std::map<mmdb::Residue *, morph_rtop_triple> *shifts, // modify the shifts
			   crunch_model_t crunch_model);

   // I fail to make a function that does a good "average" of RTops,
   // so do it long-hand by generating sets of coordinates by applying
   // each rtop to each atom - weights are transfered in the second part of the pair
   void morph_residue_atoms_by_average_rtops(mmdb::Residue *this_residue,
					     const std::vector<std::pair<clipper::RTop_orth, float> > &rtops);

   int morph_fit_by_secondary_structure_elements(const std::string &chain_id,
						 const clipper::Xmap<float> &xmap_in,
						 float map_rmsd);

   // Return a map of RTops for the residues in the fragment (and the
   // local fragment centre by which we need to move atoms when
   // applying the RTops (and move the atoms back from the origin
   // after of course).
   //
   std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> >
   morph_fit_by_secondary_structure_fragment(mmdb::Chain *chain_p, const std::string &chain_id,
					     int initSeqNum, int endSeqNum,
					     const clipper::Xmap<float> &xmap_in,
					     float map_rmsd,
					     bool simple_move);


   // fragment info for Alan:
   std::vector<coot::fragment_info_t> get_fragment_info(bool screen_output_also) const;

   std::pair<mmdb::Residue *, coot::dictionary_residue_restraints_t>
   invert_chiral_centre(const std::string &chain_id, int res_no,
			const std::string &ins_code,
			const std::string &atom_name,
			const coot::protein_geometry &geom);

   void update_bonds_using_phenix_geo(const coot::phenix_geo::phenix_geometry &b);

   void export_map_fragment_to_plain_file(float radius,
					  clipper::Coord_orth centre,
					  const std::string &filename) const;

   void globularize();

   // is_EM_map is used (1) to know not to contour outside the box
   //                   (2) to know which entry to use in Map Parameters dialog
   //                   so, if this is an-origin-centred map, then
   //                   those answers don't correspond. Needs some thought.
   bool is_EM_map() const;
   short int is_em_map_cached_flag; // -1 mean unset (so set it, 0 means no, 1 means yes)
   short int is_em_map_cached_state(); // set is_em_map_cached_flag if not set

   // user-setting over-ride internal rules for P1&909090 means EM
   void set_map_has_symmetry(bool is_em_map);

   void residue_partial_alt_locs_split_residue(coot::residue_spec_t spec,
					       int i_bond,
					       double theta,  // degrees
					       bool wag_the_dog,
					       coot::protein_geometry *geom);

   void set_user_defined_colour_indices_by_selections(const std::vector<std::pair<std::string, unsigned int> > &cis);
   void set_user_defined_colour_indices(const std::vector<std::pair<coot::atom_spec_t, int> > &cis);
   void clear_user_defined_atom_colours();

   void switch_HIS_protonation(coot::residue_spec_t res_spec);
   void reduce(coot::protein_geometry *geom_p);

   std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> >
   get_contours(float contour_level, float radius, const coot::Cartesian &centre) const;
   std::string map_units() const { std::string u = "e/A^3";
                                   if (is_EM_map()) u = "V";
                                   return u; }


#ifdef USE_MOLECULES_TO_TRIANGLES
   std::vector<std::shared_ptr<MolecularRepresentationInstance> > molrepinsts;
#endif // USE_MOLECULES_TO_TRIANGLES

   std::vector<std::pair<std::string, float> > M2T_float_params;
   std::vector<std::pair<std::string, int> >   M2T_int_params;
   //! Update float parameter for MoleculesToTriangles molecular mesh
   void M2T_updateFloatParameter(const std::string &param_name, float value);

   //! Update int parameter for MoleculesToTriangles molecular mesh
   void M2T_updateIntParameter(const std::string &param_name, int value);

   // return the index in the molrepinsts vector (can be negative for failure)
   int make_molecularrepresentationinstance(const std::string &atom_selection,
					    const std::string &colour_scheme,
					    const std::string &style);
   int add_molecular_representation(const std::string &atom_selection,
				    const std::string &colour_scheme,
				    const std::string &style,
                                    int secondary_structure_usage_flag);

   // for AlphaFold pLDDT colouring
   void add_ribbon_representation_with_user_defined_residue_colours(const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours,
                                                                    const std::string &mesh_name);
   void remove_molecular_representation(int idx);

   void delete_all_carbohydrate();

   // carbohydrate validation tools
   void glyco_tree_internal_distances_fn(const coot::residue_spec_t &base_residue_spec,
					 coot::protein_geometry *geom_p,
					 const std::string &file_name);

   // carbohydrate building - WTA for the moment
   void add_named_glyco_tree(const std::string &glycosylation_type,
                             coot::protein_geometry *geom_p,
                             const coot::residue_spec_t &res_spec,
                             const clipper::Xmap<float> &xmap);

   // hacky function to retrive the atom based on the position
   // (silly thing to do)
   mmdb::Atom *get_atom_at_pos(const coot::Cartesian &pt) const;

   void add_secondary_structure_header_records(bool overwrite=false);

   // angle in degrees.
   void spin_N(const coot::residue_spec_t &residue_spec, float angle);

   // place the O (because we have added a new residue)
   bool move_atom(const std::string &atom_name, mmdb::Residue *res_p, const clipper::Coord_orth &new_O_pos);

   int pending_contour_level_change_count;

   void crankshaft_peptide_rotation_optimization(const coot::residue_spec_t &rs,
						 unsigned int n_peptides,
						 const clipper::Xmap<float> &xmap,
						 float map_weight,
						 int n_samples,
						 ctpl::thread_pool *thread_pool_p, int n_threads);

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > peptide_C_N_pairs(const std::vector<mmdb::Residue *> &residues) const;

   mean_and_variance<float> map_histogram_values;
   mean_and_variance<float> set_and_get_histogram_values(unsigned int n_bins, bool ignore_pseudo_zeroes); // fill above

   void resolve_clashing_sidechains_by_deletion(const coot::protein_geometry *geom_p);
   void resolve_clashing_sidechains_by_rebuilding(const coot::protein_geometry *geom_p,
                                                  int imol_refinement_map);

   static int watch_mtz(gpointer data); // return 0 to stop watching
   bool continue_watching_mtz;
   updating_map_params_t updating_map_previous;
   clipper::Xmap<float> updating_map_previous_difference_map;
   int update_map_from_mtz_if_changed(const updating_map_params_t &rump);
   void update_self_from_file(const std::string &file_name);
   void update_self(const coot::mtz_to_map_info_t &mmi);

   static int watch_coordinates_file(gpointer data);
   bool continue_watching_coordinates_file;
   updating_coordinates_molecule_parameters_t updating_coordinates_molecule_previous;
   int update_coordinates_molecule_if_changed(const updating_coordinates_molecule_parameters_t &p);

   // watch for changes in the model of a different molecule (denoted by change of backup index)
   // and act on it (by updating this difference map). This needs to be static because its
   // called by g_timeout.
   static int watch_coordinates_updates(gpointer);  // for just the difference map

   static int updating_coordinates_updates_genmaps(gpointer); // oh dear, the triggers for this work the other way.
                                                              // i.e. a change in the coordinates forces
                                                              // a change in the maps, not (as above) where a
                                                              // map molecule looks for a change in the model.
                                                              // In this case, both maps are calculated together.

   int previous_backup_index;
   int other_molecule_backup_index;
   int get_other_molecule_backup_index() const { return other_molecule_backup_index; }

   // allow this to be called from the outside, when this map gets updated (by sfcalc_genmap)
   void set_mean_and_sigma(bool show_terminal=true, bool ignore_pseudo_zeroes=false);

   std::string pdb_string() const;

   coot::model_composition_stats_t get_model_composition_statistics() const;

   void shiftfield_b_factor_refinement(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                       const clipper::HKL_data<clipper::data32::Flag> &free);

   void shiftfield_xyz_factor_refinement(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                         const clipper::HKL_data<clipper::data32::Flag> &free);

   // radial colouring
   class radial_colouring_info_t {
   public:
      float radius; // from min_radius (0) to max_radius (1)
      glm::vec4 colour;
      radial_colouring_info_t(const float &r, const glm::vec4 &c) : radius(r), colour(c) {}
      float fraction_of_range(float min_radius, float max_radius) const {
         float delta = max_radius - min_radius;
         float x = radius - min_radius;
         return x/delta;
      }
   };

   class radial_colouring_info_container_t {
   public:
      std::vector<radial_colouring_info_t> colour_stops;
      void sort_colour_stops(); // smallest at the top
   };

   void set_radial_map_colouring_centre(float x, float y, float z);
   void set_radial_map_colouring_min_radius(float r);
   void set_radial_map_colouring_max_radius(float r);
   void set_radial_map_colouring_invert(bool invert_state);
   void set_radial_map_colouring_saturation(float saturation);
   void set_radial_map_colouring_do_radial_colouring(bool state) {
      if (state != radial_map_colouring_do_radial_colouring) {
         radial_map_colouring_do_radial_colouring = state;
         update_map(true);
      }
   }
   bool radial_map_colouring_do_radial_colouring;
   clipper::Coord_orth radial_map_colour_centre;
   double radial_map_colour_radius_min;
   double radial_map_colour_radius_max;
   double radial_map_colour_saturation;
   bool   radial_map_colour_invert_flag;

   // colour by other map (e.g. correlation)
   bool colour_map_using_other_map_flag;
   void set_colour_map_using_other_map(bool state) {
      colour_map_using_other_map_flag = state;
   }

   GdkRGBA position_to_colour_using_other_map(const clipper::Coord_orth &position);


   coot::density_contour_triangles_container_t export_molecule_as_x3d() const;
   bool export_molecule_as_obj(const std::string &file_name);
   bool export_map_molecule_as_obj(const std::string &file_name) const;
   bool export_model_molecule_as_obj(const std::string &file_name);
   bool export_molecule_as_gltf(const std::string &file_name) const;
   bool export_map_molecule_as_gltf(const std::string &file_name) const;
   bool export_model_molecule_as_gltf(const std::string &file_name) const;

   void export_these_as_3d_object(const std::vector<coot::api::vertex_with_rotation_translation> &vertices,
                                  const std::vector<g_triangle> &triangles);

   bool write_model_vertices_and_triangles_to_file_mode;
   bool export_vertices_and_triangles_func(const std::vector<coot::api::vertex_with_rotation_translation> &vertices,
                                           const std::vector<g_triangle> &triangles);
   std::string export_vertices_and_triangles_file_name_for_func;

   // These meshes are not the way coot 0.9 organized generic display objects.
   //
   // meshes are drawn with draw_meshed_generic_display_object_meshes()
   // and instanced_meshes are drawn with draw_instanced_meshes().
   //
   // these are for specific molecule-based objects using regular Mesh

   std::vector<Mesh> meshes;
   // these are for specific molecule-based objects using instancing Mesh
   std::vector<Instanced_Markup_Mesh> instanced_meshes;
   Instanced_Markup_Mesh &find_or_make_new(const std::string &mesh_name);
   // Mesh mesh_for_symmetry_atoms;
   model_molecule_meshes_t meshes_for_symmetry_atoms;

   // And now symmetry atoms are displayed as a Mesh
   bool this_molecule_has_crystallographic_symmetry;

   // either we have licorice/ball-and-stick (licorice is a form of ball-and-stick) or big-ball-no-bonds
   Mesh::representation_mode_t model_representation_mode;
   void set_model_molecule_representation_style(unsigned int mode);

   // These meshes are the molecule, replacing the inital way of representing the molecule. Uses
   // instances of cylinders and spheres and hemispheres. Put them in a Model at some stage.
   std::vector<glm::vec4> make_colour_table() const;

   static glm::vec4 get_glm_colour_func(int idx_col, int bonds_box_type);

   // 20230813-PE user-defined colours that come together with the instanced bonds.
   // adding colours using the functions below add into user_defined_colours
   std::vector<coot::colour_holder> user_defined_bond_colours;

   //! user-defined colour-index to colour
   //! (internallly, this converts the `colour_map` to the above vector of colour holders, so it's probably a good idea
   //! if the colour (index) keys are less than 200 or so.
   void set_user_defined_bond_colours(const std::map<unsigned int, std::array<float, 3> > &colour_map);

   //! user-defined atom selection to colour index
   void set_user_defined_atom_colour_by_selection(const std::vector<std::pair<std::string, unsigned int> > &indexed_residues_cids,
                                                  bool apply_to_non_carbon_atoms);

   void make_mesh_from_bonds_box();
   void make_meshes_from_bonds_box_instanced_version(); // fills the below meshes (for instancing)
   // 20230826-PE here we go with the new "instancing.hh" wrapper class
   model_molecule_meshes_t model_molecule_meshes;
   void set_material_in_molecules_as_mesh(const Material &material) {
      model_molecule_meshes.set_material(material);
   }
   Mesh molecule_as_mesh_rama_balls;
   Mesh molecule_as_mesh_rota_dodecs;
   // pass this function to the Mesh so that we can determine the atom and bond colours
   void draw_molecule_as_meshes(Shader *shader_p,
                                stereo_eye_t eye,
                                const glm::mat4 &mvp,
                                const glm::mat4 &view_rotation_matrix,
                                const std::map<unsigned int, lights_info_t> &lights,
                                const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                const glm::vec4 &background_colour,
                                bool do_depth_fog);

   // instanced models
   void draw_molecule_as_meshes_for_ssao(Shader *shader_for_meshes_for_ssao,
                                         Shader *shader_for_instanced_meshes_for_ssao,
                                         const glm::mat4 &model_matrix,
                                         const glm::mat4 &view_matrix,
                                         const glm::mat4 &projection_matrix);
   // instanced models
   void draw_molecule_as_meshes_with_shadows(Shader *shader,
                                             const glm::mat4 &mvp,
                                             const glm::mat4 &model_rotation_matrix,
                                             const std::map<unsigned int, lights_info_t> &lights,
                                             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                             float opacity,
                                             const glm::vec4 &background_colour,
                                             bool do_depth_fog,
                                             const glm::mat4 &light_view_mvp,
                                             unsigned int shadow_depthMap,
                                             float shadow_strength,
                                             unsigned int shadow_softness, // 1, 2 or 3.
                                             bool show_just_shadows);

   void draw_symmetry(Shader *shader_p,
                      const glm::mat4 &mvp,
                      const glm::mat4 &mouse_based_rotation_matrix,
                      const std::map<unsigned int, lights_info_t> &lights,
                      const glm::vec3 &eye_position,
                      const glm::vec4 &background_colour,
                      bool do_depth_fog);


   // float scale_factor 4 , float offset 3
   void recolour_ribbon_by_map(const clipper::Xmap<float> &xmap, float scale_factor, float offset);

   void fill_chiral_volume_outlier_marker_positions(int state);
   void set_show_non_bonded_contact_baddies_markers(int state);
   // 20230829-PE
   // draw_chiral_volume_outlier_markers_flag is per molecule
   // draw_bad_nbc_atom_pair_markers is global (only one). Maybe this is a mistake
   bool draw_chiral_volume_outlier_markers_flag;
   std::vector<glm::vec3> chiral_volume_outlier_marker_positions;

   std::vector<glm::vec3> unhappy_atom_marker_positions;

   bool read_nef(const std::string &file_name);


};

#endif // MOLECULE_CLASS_INFO_T
