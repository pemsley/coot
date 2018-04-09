/* ligand/ligand.hh
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Copyright 2009 by The University of Oxford.
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
 * 02110-1301, USA.
 */

#ifndef COOT_LIGAND_HH
#define COOT_LIGAND_HH

#include <string>
#include <vector>

#include <mmdb2/mmdb_manager.h>

#include "mini-mol/mini-mol.hh"
#include "geometry/protein-geometry.hh" // for wiggly ligands

#include "clipper/core/xmap.h"
#include "clipper/contrib/skeleton.h" // neighbs is in the recursive function call
#include "clipper/core/map_utils.h" // map stats (returned value)

#include "scored-coord.hh"
#include "ligand-score-card.hh"
#include "scored-molecule.hh"
#include "map-point-cluster.hh"

namespace coot {

   // return a molecule of atoms (first) used to mask the map
   // (i.e. those not in the selection) and (second) a molecule of
   // atoms in the selection [strange but true]
   // 
   std::pair<minimol::molecule, minimol::molecule>
      make_mols_from_atom_selection_string(mmdb::Manager *mol,
					   std::string atom_selection_string,
					   bool fill_masking_molecule_flag);

   class ligand {

   private:
      enum { OK_GOLDILOCKS, TOO_FAR, TOO_CLOSE, WATER_STATUS_UNKNOWN};
      clipper::Xmap<float> xmap_pristine;
      clipper::Xmap<float> xmap_cluster; // the map that gets masked/scribbled on
      clipper::Xmap<float> xmap_masked;  // the post-masking, pre-clustering map 
      std::vector<clipper::Coord_orth> points;
      int n_clusters;
      std::vector<map_point_cluster> cluster;
      short int do_cluster_size_check_flag;
      short int do_chemically_sensible_test_flag;
      short int do_sphericity_test_flag;
      float water_molecule_volume;

      // There is no operator= for clipper::Map_stats
      // values: set, mean, std_dev:
      std::pair<short int, std::pair<float, float> > xmap_masked_stats;

      short int verbose_reporting;
      coot::minimol::molecule protein_atoms; // keep this for water
				             // testing.
      
      clipper::Coord_orth protein_centre;
      std::vector<clipper::Coord_orth> initial_ligand_model_centre;
      std::vector<clipper::Mat33<double> > initial_ligand_eigenvectors;
      std::vector<std::vector<double> >    initial_ligand_eigenvalues;
      
      // For resolving the contributions to \alpha_x, \alpha_y,
      // \alpha_z of the gradients perpendicular to the vector between
      // the atom position (x,y,z) and the centre of rotation:
      clipper::RTop_orth rotation_component[3];

      // high resolution maps have much bigger gradients than low
      // resolution maps.  So we need to scale the gradients.  Here's
      // where we keep that scale:
      double gradient_scale;
      void calculate_gradient_scale(); 
      
      float cut_off; // For efficiency purposes, normally I wouldn't
		     // approve (I would pass it as an argument).
		     // This is the level above which grid points are
		     // considered part of a cluster


      // Here we put the 4 orientations (one of which is the identity)
      // which we rotate the ligand to match the various ways in which
      // the eigenvectors of the ligand molecule and match the
      // eigenvectors of the density cluster.  They are set at
      // construction time.
      std::vector<clipper::Mat33<double> > origin_rotations; 
      std::vector<clipper::Mat33<double> > origin_unrotations; 
      std::string ligand_filename(int n_count, int ior) const;
      
      void trace_along(const clipper::Coord_grid &cg_start,
		       const clipper::Skeleton_basic::Neighbours &neighb,
		       int n_clusters);
      void calculate_cluster_centres_and_eigens();
      void move_ligand_sites_close_to_protein(int iclust);
      void move_ligand_site_close_to_protein_using_shape(int iclust,
							 const std::vector<clipper::Coord_orth> &sampled_protein_coords);
      // and a function that is used there:
      double min_dist_to_protein(const clipper::Coord_orth &point,
				 const std::vector<clipper::Coord_orth> &sampled_protein_coords) const;
      // make the sampled coordinates from the protein mask:
      std::vector <clipper::Coord_orth> make_sample_protein_coords() const;

      
      void move_ligand_centres_close_to_protein(const std::vector<clipper::Coord_orth> &sampled_protein_coords);
      short int cluster_ligand_size_match(int iclust, int ilig);
      // recommended interface:
      void mask_around_coord(const clipper::Coord_orth &co, float atom_radius);
      void mask_around_coord(const clipper::Coord_orth &co, float atom_radius,
			     clipper::Xmap<float> *xmap_p);  // alter xmap_p
      void mask_around_coord(const clipper::Coord_orth &co, float atom_radius,
			     clipper::Xmap<int> *xmap_p);  // alter xmap_p, template?
      int make_selected_atoms(mmdb::PPAtom *atoms_p, mmdb::Manager *mol); 

      // for graphics "refine one residue" fitting, we don't want to
      // test orientations other than the first (identity) one.
      short int dont_test_rotations;

      // Map that is close to protein has this value:
      //
      float masked_map_val;

      // Map that is as close as this or closer is masked:
      // default: 
      float map_atom_mask_radius; 

      // For graphics, we want to pass back the solution
      // (get_solution), not write them out.
      // 
      short int write_solutions;
      short int write_orientation_solutions;
      short int write_raw_waters;

      // solvent stuff
      clipper::Coord_orth move_atom_to_peak(const clipper::Coord_orth &a,
					    const clipper::Xmap<float> &search_map) const;
      short int cluster_is_possible_water(int i) const; 
      short int cluster_is_possible_water(const coot::map_point_cluster &mpc) const;
      short int has_sphericalish_density(const clipper::Coord_orth &a,
					 const clipper::Xmap<float> &search_map) const;

      // is close to protein, H-bonding distance to an O or N atom.
      short int water_pos_is_chemically_sensible(clipper::Coord_orth new_centre) const;
      short int water_pos_is_chemically_sensible(const clipper::Coord_orth &water_centre,
						 const std::vector<clipper::Coord_orth> &extra_sites) const;
      float density_at_point(const clipper::Coord_orth &a,
			     const clipper::Xmap<float> &search_map) const;


      // make this private?
      void fit_ligands_to_cluster(int ilig);

      // clipper::RTop_orth ligand_transformation(int i_cluster) const; // old
      clipper::Coord_orth transform_ligand_atom(clipper::Coord_orth a_in,
						int ilig, int iclust, int ior) const;
      clipper::Coord_orth transform_ligand_atom(clipper::Coord_orth a_in,
						int ilig, int iclust, int ior,
						const clipper::RTop_orth &eigen_ori) const;
      clipper::Coord_orth transform_ligand_atom(clipper::Coord_orth a_in,
						int ilig,
						const clipper::RTop_orth &cluster_rtop,
						int ior) const;

      std::vector <std::vector<minimol::molecule> > fitted_ligand_vec;

      // return  \alpha_x, \alpha_y, \alpha_z in radians:
      clipper::Vec3<double>
      get_rigid_body_angle_components(const std::vector<minimol::atom *> &atoms,
				      const clipper::Coord_orth &mean_pos,
				      const std::vector<clipper::Grad_orth<float> > &grad_vec) const;
      void apply_angles_to_ligand(const clipper::Vec3<double> &angles,
				  const std::vector<minimol::atom *> *atoms_p,
				  const clipper::Coord_orth &mean_pos);

      clipper::Coord_orth 
	mean_ligand_position(const std::vector<minimol::atom *> &atoms) const;

      clipper::Mat33<double> mat33(const clipper::Matrix<double> &mat) const;

      void make_ligand_properties(int ilig);
      short int similar_eigen_values(int iclust, int ilig) const;

      // water stuff
      coot::minimol::molecule water_molecule;  // a molecule of (many) waters
      minimol::molecule new_ligand_with_centre(const clipper::Coord_orth &pos);
      void write_waters(const std::vector<clipper::Coord_orth> &water_list,
			const std::string &filename) const;

      // called by the public interface(s)
      void mask_map(bool mask_waters_flag); 

      // initially 0.75
      float fit_fraction;

      // return the score (of rigid body fitting - you may want to
      // expand this to each cluster in future):
      //
      std::vector<ligand_score_card> save_ligand_score; // a vector over n-clusters,
                                                        // indexed with iclust

      // tinker with mmmol (minimol molecule)
      void set_cell_and_symm(coot::minimol::molecule *mmmol) const; 

      int n_grid_limit_for_water_cluster() const;

      // blobs that are too big to be waters
      // 
      std::vector<std::pair<clipper::Coord_orth, double> > keep_blobs;

   protected:

      float default_b_factor;
      const clipper::Xmap<float> & Xmap() const { return xmap_cluster;}

      // keep a record of the map rms so that we can kludge the
      // gradient scale factor.

      float map_rms;

      // So that we only sample the coords once, not for every round
      // of water picking.
      void find_clusters_internal(float z_cutoff,
				  const std::vector<clipper::Coord_orth> &sampled_protein_coords);

      // we use these in residue_by_phi_psi:
      //

      ligand_score_card
      fit_ligand_copy(int iclust, int ilig, int ior); // fit a copy of the
	   				      // initial_ligand to the
					      // cluster i.
      ligand_score_card
      fit_ligand_copy(int iclust, int ilig, int ior,
		      const clipper::RTop_orth &eigen_ori); // fit a copy of the
	   				      // initial_ligand to the
					      // cluster iclust.

      // we want these available to a derived class
      std::vector<minimol::molecule> initial_ligand;

      // We used to have a final_ligand for every iclust.  Now we have
      // a vector of solutions for every iclust (so that we can pass
      // back the top 10 (say) solutions for a given position and
      // allow the user to examine).
      // std::vector<coot::minimol::molecule>   final_ligand;
      std::vector<std::vector<std::pair<coot::minimol::molecule, ligand_score_card> > > final_ligand;
      void sort_final_ligand(unsigned int iclust);
      // which uses comparison function
      // return "is first less than second?" status
      static bool compare_scored_ligands(const std::pair<minimol::molecule, ligand_score_card> &sl_1,
					 const std::pair<minimol::molecule, ligand_score_card> &sl_2);
      static bool
      compare_scored_ligands_using_correlation(const std::pair<minimol::molecule, ligand_score_card> &sl_1,
					       const std::pair<minimol::molecule, ligand_score_card> &sl_2);

      double get_correl(const minimol::molecule &mol) const;
      void write_orientation_solution(unsigned int iclust,
				      unsigned int ilig,
				      unsigned int i_eigen_ori,
				      unsigned int ior,
				      const minimol::molecule &mol) const;
				      


      // Consider doing away with this sort of manipulation.  It is
      // not very elegant, and instead we could (should) make the
      // transformations by altering the Coord_orth (from the origin
      // to the density) only.  When the refinement is complete - then
      // we move the (saved in fitted_ligand_vec) atoms. 
      void rigid_body_refine_ligand(std::vector<minimol::atom *> *atoms, 
				    const clipper::Xmap<float> &xmap_fitting);
      ligand_score_card
      score_orientation(const std::vector<minimol::atom*> &atoms, 
			const clipper::Xmap<float> &xmap_fitting,
			bool use_linear_interpolation=false) const;
      float
      score_position(const clipper::Coord_orth &pt,
		     const clipper::Xmap<float> &xmap_fitting) const;

      // Here's a function that someone else might want to use:
      //
      std::string get_first_residue_name(const coot::minimol::molecule &mol) const;

      std::vector<clipper::Coord_orth> water_fit_internal(float sigma_cutoff, int n_cycle);

   float var_limit;
      float water_to_protein_distance_lim_max;
      float water_to_protein_distance_lim_min;

      // single entities for rigid body refinement do not need size
      // matching.  It only causes problems with grid sampling...
      short int do_size_match_test;

      // save the z cut off for waters
      float z_cut_off_in_save;

      bool close_to_another(const clipper::Coord_orth &p1,
			    const std::vector<clipper::Coord_orth> &ref,
			    const double &d_crit) const;

      
   public:
      ligand();
      // ~ligand(); // causes core dump currently (badness).
      // void mask_map(const clipper::DBAtom_selection &atoms);

      void mask_map(mmdb::Manager *mol, short int mask_waters_flag); // used in c-interface-ligand.cc

      // This is the prefered interface
      void mask_map(mmdb::Manager *mol, 
		    int SelectionHandle,
		    short int invert_flag); // if invert_flag is 0 put the map to 0
                                            // where the atoms are.  If 1, put map to
                                            // 0 where the atoms are not (e.g. ligand
                                            // density figure)

      void mask_map(const minimol::molecule &mol,
		    short int mask_waters_flag); // coot usage (coot
						   // has prepared a
						   // mol that does
						   // not have in it
						   // the residue of
						   // interest).

      void import_map_from(const clipper::Xmap<float> &map_in);
      // Create a new form of the function that sets the rms of the map,
      // which then gets used in a kludgey way to set the gradient scale.
      void import_map_from(const clipper::Xmap<float> &map_in, float rms_map_in);
      // export map. Ouch, heavy!
      // Xmap_plus_bits masked_map_plus_bits() const; not today.
      const clipper::Xmap<float> &masked_map() const;
      std::string masked_map_name() const {
	 return std::string("Masked (by protein)"); }

      // This is a console/testing function.  Should not be used in a
      // real graphics program.  Use instead import_map_from() with a
      // precalculated map.
      // 
      // return 0 on error, 1 on success
      // 
      short int map_fill_from_mtz(std::string mtz_file_name,
				  std::string f_col,
				  std::string phi_col,
				  std::string weight_col,
				  short int use_weights,
				  short int is_diff_map,
				  float map_sampling_rate=1.5);

      // This is a console/testing function.  Use instead mask_map()
      // in a real graphics program (mask_by_atoms runs mask_map).
      //
      // Return status of the map read and if the masking succeeded.
      // 1 is success.
      // 
      int mask_by_atoms(std::string pdb_filename);
      int mask_by_atoms_mapview_way(std::string pdb_filename);
      void cluster_test();
      void find_clusters(float z_cut_off); // number of standard deviations
      void find_clusters_old(float z_cut_off);
      void find_clusters_int(float z_cut_off);

      // and now we want to give a cluster centre (actually, a point
      // close to a cluster centre).  We need to create a
      // map_point_cluster and calulate its centre and eigens, like we
      // do when given a whole map.
      // 
      void cluster_from_point(clipper::Coord_orth pt, float cut_off);

      // cluster size maxes out, starts new.
      void find_clusters_water_flood(float cut_off,
				     const std::vector <clipper::Coord_orth> &sampled_protein_coords);

      // Call this when rigid body refining (rather than ligand search):
      // 
      void find_centre_by_ligand(short int do_size_match_flag);
                                    // use the installed ligand to
				    // provide a centre.  This is for
				    // use in "rigid body refinement",
				    // where the "ligand" (range of
				    // residues (typically 1)) is at
				    // the place from where we start.
                                    //
                                    // And in that case, we don't want
                                    // to test for the size of the
                                    // ligand (model) vs the size of
                                    // the density (cluster).
                                    // So the flag should be set to 0.

      // Is this function of any use really?  Each new ligand generates a
      // new cluster against which the ligands (all of them) are tested.
      // 
      void find_centres_by_ligands(); 
      void find_centre_by_ligand_internal(int ilig); 
      
      void output_centres(); 
      void print_cluster_details(bool print_grid_points_flag=false) const; 
      clipper::Map_stats map_statistics();
      
      // a check for writability is made.
      void output_map(std::string filename) const; // write_map, perhaps (xmap_cluster).
      // must be writable or crashes, debugging only.
      void output_map(const clipper::Xmap<float> &xmap, const std::string &filename) const; 
      
      // Turn the cluster centres into anisotropic atoms, so that we can see
      // the direction eigenvalues as if it were anisotropy of an atom.
      // This involves some tricky conversion of the eigenvectors to
      // anisotropic Us (the "diagonalized form"?).
      void make_pseudo_atoms();

      void install_ligand(std::string filename);
      void install_ligand(mmdb::Manager *mol);
      void install_ligand(const minimol::molecule &ligand);

      // return whether we successfully installed wiggly ligands or
      // that the ligand couldn't wiggle for some reason
      // 
//       short int install_wiggly_ligands(const coot::protein_geometry &pg,
// 				       const minimol::molecule &ligand,
// 				       int n_samples);

      // we use 1 for max_placements for rigid body refinement
      // 
      // Go down the cluster list, starting at biggest cluster fitting ligands
      // until the cluter number is max_placements
      // 
      void fit_ligands_to_clusters(int max_placements); // for the top n_clusters

      // Yet more tinkering parameters now accessable from the outside
      void set_map_atom_mask_radius(float f) { map_atom_mask_radius = f; }
      
      // Eleanor wants the peak search (unfiltered) results:
      void set_write_raw_waters() {
	 write_raw_waters = 1; }
      
      // When we are rigid body refining residues, we don't want to
      // twist them arround.  Used typically in conjuction with
      // find_centre_by_ligand. 
      // 
      void set_dont_test_rotations() { dont_test_rotations = 1; }

      //
      void set_masked_map_value(float val); 

      // For graphics we use get_solution, not write out pdb files.
      // 
      void set_dont_write_solutions() { write_solutions = 0;}
      //
      // this is the number of clusters, not the number of solutions
      // for a given cluster. 
      // 
      unsigned int n_clusters_final() const { return final_ligand.size(); }

      // return the number of scored solutions with score of more than
      // frac_limit_of_peak_score (say, 0.5)
      // 
      unsigned int n_ligands_for_cluster(unsigned int iclust,
					 float frac_limit_of_peak_score) const;

      // generate correlation scores for the top n_sol solutions and
      // re-sort.
      //
      // clusters should already be in "score" order by now (biggest
      // score at the top (0th position)).
      //
      void score_and_resort_using_correlation(unsigned int iclust, unsigned int n_sol);

      // get the ith solution for the iclust cluster
      coot::minimol::molecule get_solution(unsigned int isolution,
					   unsigned int iclust) const;
      //
      ligand_score_card get_solution_score(int iclust, int isolution) const {
	 return final_ligand[iclust][isolution].second;
      }

      // People like to tweak the water finding parameters:
      void set_variance_limit(float f) { var_limit = f; }
      void set_water_to_protein_distance_limits(float f1, float f2) {
	 water_to_protein_distance_lim_max = f1;
	 water_to_protein_distance_lim_min = f2;
	 if (water_to_protein_distance_lim_max < water_to_protein_distance_lim_min) {
	    float tmp = water_to_protein_distance_lim_max;
	    water_to_protein_distance_lim_max = water_to_protein_distance_lim_min;
	    water_to_protein_distance_lim_min = tmp;
	 }
      }

      double possible_water_volume_limit() const;

      // 
      void water_fit_old(int n_cycles);
      void water_fit(float sigma_cutoff, int n_cycles);
      void flood();
      void flood2(float n_sigma); // the return of flood(), this time
				  // by peak picking. // find peaks
				  // that are in density greater than
				  // n_sigma
      // which uses:
      std::vector<clipper::Coord_orth>
      move_waters_close_to_protein(const std::vector<clipper::Coord_orth> &water_list,
				   const std::vector<clipper::Coord_orth> &sampled_protein_coords) const;
      clipper::Coord_orth move_water_close_to_protein(const clipper::Coord_orth &water_pos,
						      const std::vector<clipper::Coord_orth> &sampled_protein_coords) const;


      coot::minimol::molecule water_mol() { return water_molecule; }

      // a number between 0 (none) and 1.0 (all) atoms need to be
      // positive (initally set to 0.75)
      // 
      void set_acceptable_fit_fraction(float f) { fit_fraction = f; }
      //
      
      //
      void set_verbose_reporting() { verbose_reporting = 1; }

      void set_default_b_factor(float f);

      std::vector<std::pair<clipper::Coord_orth, double> > big_blobs() const { 
	 return keep_blobs;
      }

      // this should only be run post-sort (post-correlation sort)
      void limit_solutions(unsigned int iclust,
			   float frac_max_correl_lim,
			   int max_n_solutions,
			   float tolerance=0.0,
			   bool filter_by_torsion_match=false);

      // This list of peaks above a particular level (and one day
      // optionally below the same but negative level):
      // 
      std::vector<std::pair<clipper::Coord_orth, float> > cluster_centres() const; 

      // stop cluster size checking.  We want to do this (i.e. no
      // check) when we are trying to flood a map with coordinates,
      // rather than just the water sites.
      // 
      void set_cluster_size_check_off() { do_cluster_size_check_flag = 0; }; 
      void set_chemically_sensible_check_off() { do_chemically_sensible_test_flag = 0; }; 
      void set_sphericity_test_off() { do_sphericity_test_flag = 0;}

      // sets protein atoms for moving ASU peaks to around the protein
      void import_reference_coords(mmdb::Manager *mol);
      
      // external error checking
      // 
      short int masking_molecule_has_atoms() const;


      /* flip the ligand (usually active residue) around its eigen
	 vectors to the next flip number.  Immediate replacement (like
	 flip peptide). We need to undo the current flip number first
	 though (if flip_number is not 1). */
      coot::minimol::molecule flip_ligand(short int flip_number) const;
   }; 


} // namespace coot

#endif // COOT_LIGAND_HH

