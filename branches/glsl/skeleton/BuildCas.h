
// mode: -*-c++-*-

#include <string>
#include <mmdb/mmdb_manager.h>

#include "mmdb-extras.h"
#include "clipper/core/xmap.h"
#include "clipper/contrib/skeleton.h"

#include "Cartesian.h"

#include "AngleInfo.h"

#ifndef RADTODEG
#define RADTODEG 57.2957795147
#endif

#ifndef PI_BY_2
#define PI_BY_2   1.570796327
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

class cv_list_t { 

 public:

   coot::Cartesian geometrical_vector; 
   int atom_index;

   cv_list_t(coot::Cartesian vec, int a_index) { 
      geometrical_vector = vec; 
      atom_index = a_index; 
   } 

}; 


// Waaah... I want multiple return values... grumble....
// 
class score_and_cart { 
    
public: 

   coot::Cartesian pos;
   float score; 
   clipper::Coord_grid near_grid_point; 

   score_and_cart() { 

      // set to some obviously bogus initial values: 
      // 
      coot::Cartesian a(-0.999, -0.999, -0.999); 

      pos = a;
      score = -888.8; 
   } 
    
}; 

// 
class scores { 
public: 
   float score; 
   float density_score_val; 
   float deviation_from_ideal_length_score_val; 
   float branch_point_proximity_score_val; 

}; 

// These are atom positions, used in recursive Ca model building. 
// 
class TreeNode { 

public:
   const TreeNode *parent; 
   coot::Cartesian pos; 
   clipper::Coord_grid near_grid_point; 

   void setup (const TreeNode *par, coot::Cartesian inpos) { 
      parent = par; 
      pos = inpos; 
   }

   void setup (const TreeNode *par, score_and_cart sc) { 
      parent = par; 
      pos = sc.pos;
      near_grid_point = sc.near_grid_point; 
   }

   void setup (const TreeNode *par, coot::Cartesian inpos, clipper::Coord_grid cg) { 
      parent = par; 
      pos = inpos;
      near_grid_point = cg; 
   }
   
   TreeNode(TreeNode *par, coot::Cartesian inpos) { 
      setup(par, inpos); 
   }
   
   TreeNode() {
      parent = NULL;
   } 

   TreeNode(coot::Cartesian pos_in) { 
      parent = NULL;
      pos = pos_in; 
   }

   TreeNode(score_and_cart sc) { 

      parent = NULL; 
      pos = sc.pos; 
      near_grid_point = sc.near_grid_point; 
   }
};


// I want multiple return values again:
// 
class small_triangle_thing { 

public: 
   coot::Cartesian pos; 
   short int is_small_triangle;
   short int is_new_small_triangle; 

}; 

// These are points (clipper::Coord_grids) in a skeleton map and are
// used to see if target point is "close to" the current point (by
// being less than some (gridding dependent) critical value number of
// nodes away).
// 
class SkeletonTreeNode { 

public: 

   // I get stuck with de-referencing the ix of a neighbour (I mean I
   // don't know how to set it - we have a Coord_grid, not a
   // Map_reference_index).
   // 
   // clipper::Xmap_base::Map_reference_index ix;

   // clipper::Coord_grid cg; 

   // please be fast enough...
   // 
   // Do not include yourself.
   vector<clipper::Coord_grid> neighbs; 

   clipper::Coord_grid near_grid_point; 

}; 

namespace coot { 
   class Cartesian_and_Grid { 

   public:
      coot::Cartesian pos; 
      clipper::Coord_grid near_grid_point; 

      Cartesian_and_Grid(coot::Cartesian p, const clipper::Coord_grid &cg) { 
	 pos = p; 
	 near_grid_point = cg; 
      } 
   };
}


class asc_and_grids { 

public:
   atom_selection_container_t asc; 
   vector<clipper::Coord_grid> grid_points; 
   
   asc_and_grids(atom_selection_container_t a, vector<clipper::Coord_grid> gp) {
      
      asc = a; 
      grid_points = gp; 
   }

}; 



class BuildCas { 

   // a vector of vectors of Ca builds:
   // 
   vector<vector <score_and_cart> > build; 
   // an index into the build array for the "string" we are building now.
   // 
   int i_current_build, i_max_build;  // i_max_build is not the maximum
				      // number of builds, but the
				      // number of the top build in the
          		              // build array.

   int segment_map_filled; 

   float map_cut_off; 

   clipper::Xmap<int> segment_map; 
  
  // 
  int n_fitted_in_current_segment; 

  float grid_dependent_distance_param; 
  short int branch_point_have_been_expanded_flag; 

  vector<coot::Cartesian> branch_points;
  vector<coot::Cartesian> branch_points_symm_expanded; 

   // introduce AngleInfo
   AngleInfo angle_info;

   // OK, I give in.  Usually I don't like doing this sort of things,
   // I prefer to pass the parameter.  But this map gets everywhere.
   // So lets keep a pointer to it.
   // 
   // This gets set when we: 
   // set_density_map

   clipper::Xmap<float> *d_map_p; 

   // New fangled (experimental) construction 
   // Lets see if it is of any use.  7/10/2002 - PE.
   // 
   clipper::Xmap<SkeletonTreeNode> treenodemap; 
   short int treenodemap_is_filled; 

   // 
   coot::Cartesian expansion_centre; // the "centre" of the molecule
   short int expansion_centre_is_set; 

 public:

   vector <coot::Cartesian> big_ball;
   vector <clipper::Coord_grid> big_ball_grid; 

   BuildCas() {
      setup_internal();
   };


   // the prefered constructor:
   // 
   BuildCas(clipper::Xmap<float> &map_in, float map_cut) { 
      setup_internal(); 
      set_density_map_and_cut(map_in, map_cut); 
   }

   short int 
   old_is_same_gridpoint(const clipper::Coord_grid &g1, const clipper::Coord_grid &g2) const { 

      cout << "is_same_gridpoint: " 
	   << g1.unit(d_map_p->grid_sampling()).format() << " and "
	   << g2.unit(d_map_p->grid_sampling()).format() << endl; 
      return (g1.unit(d_map_p->grid_sampling()) == g2.unit(d_map_p->grid_sampling())); 

   } 

   short int 
   is_same_gridpoint(const clipper::Coord_grid &g1, const clipper::Coord_grid &g2) const {

//       cout << "is_same_gridpoint: " 
// 	   << g1.unit(d_map_p->grid_sampling()).format() << " and "
// 	   << g2.unit(d_map_p->grid_sampling()).format() << endl; 
      return (treenodemap.get_data(g1).near_grid_point == treenodemap.get_data(g2).near_grid_point); 
   } 
   

   void setup_internal(); 

   void set_density_map_and_cut(clipper::Xmap<float> &map_in, float map_cut); 
  
  // a wrapper, depending on what we already have in this segment:
  // 
  score_and_cart old_fit(const clipper::Xmap<float> &map);

  // We have no points in the current segment
  // 
  // There is no prior information.
  // 
  // Pick any old branch point? Or something like that.
  // 
  score_and_cart fit_first_in_segment (const clipper::Xmap<float> &map);

  // We have 1 point in the current segment.
  // Prior probability density is distance only.
  // 
  score_and_cart fit_second_in_segment(const clipper::Xmap<float> &map);
  
  // We have 2 ponts in the current segment
  //
  // Prior probability density is in the form of selection is distance and angle.
  // 
  score_and_cart fit_third_in_segment (const clipper::Xmap<float> &map);

  // We have 3 other points in the current segment.
  // 
  // Prior probability density is in the form of torsion, angle and distance
  // 
  score_and_cart fit_forth_in_segment (const clipper::Xmap<float> &map);

   // The "normal" version (case) where we have the previous 2 torsions
   // available. Prior probability density is in the form of previous 2
   // torsions and angles and distance. 
   // 
   score_and_cart fit_next_in_segment  (const clipper::Xmap<float> &map);
   
   // debugging 
   int n_fitted_in_current_seg() { return n_fitted_in_current_segment; }; 
   
   vector<coot::Cartesian_and_Grid>
   cluster_centres(vector<vector<coot::Cartesian_and_Grid> > cluster_vec) const;
   
   vector<vector<coot::Cartesian_and_Grid> > 
   cluster_bones_points(vector<coot::Cartesian_and_Grid> big_ball, coot::Cartesian centre_point) const; 
   
   vector <coot::Cartesian>
   point_list_by_symmetry(atom_selection_container_t AtomSel,
			  const vector<clipper::Coord_grid> &grids,
			  coot::Cartesian current_point, 
			  float radius, short int use_grids = 1); 
   
   atom_selection_container_t
   convert_to_atoms(const clipper::Xmap<int> &l1,
		    const vector<coot::Cartesian> &c, 
		    std::string molecule_name) const; 
   
   atom_selection_container_t
   convert_to_atoms(const clipper::Xmap<float> &l1,
		    const vector<coot::Cartesian> &c, 
		    std::string molecule_name) const;
   
   atom_selection_container_t
   convert_to_atoms(const vector<coot::Cartesian> &c,
		    std::string molecule_name) const; 
   
   
   atom_selection_container_t
   convert_to_atoms_internal(clipper::Spacegroup spg, clipper::Cell cell, 
			     const vector<coot::Cartesian> &c,
			     short int diff_residue_flag,
			     std::string molecule_name) const; 


  vector<coot::Cartesian_and_Grid>
     select_by_distance(coot::Cartesian start_point, float near, float far) const; // 3.7A +/- a bit

   // atom_selection_container_t 
   asc_and_grids
     all_skel_pts_in_asu(const clipper::Xmap<float> &map,
			 const clipper::Xmap<int>   &l1,
			 float cut_off) const; 

   asc_and_grids toplevel_skel_pts_in_asu() const; 

  coot::Cartesian
    move_by_symmetry(coot::Cartesian target_point, 
		     coot::Cartesian middle_mol, 
		     CMMDBCryst *cryst_p) const ; 

  atom_selection_container_t 
   build_big_ball(const clipper::Xmap<float> &map, 
		  atom_selection_container_t asc,
		  const vector<clipper::Coord_grid> &grids);
 
  int 
    count_and_mark_segments(const clipper::Xmap<int>   &skel,
			    const clipper::Xmap<float> &map,
			    float cut_off); 

  void
    trace_along(const clipper::Coord_grid &c_g_start, 
		const clipper::Skeleton_basic::Neighbours &neighb,
		int i_segment_number,
		int i_max_level, 
		float cut_off); 

  coot::Cartesian
    position_by_torsion(float theta_2, float torsion, float dist) const; 

  coot::Cartesian
    position_by_torsion(coot::Cartesian Atom_1, coot::Cartesian Atom_2, coot::Cartesian Atom_3,
			float theta_2, float torsion, float dist) const; 

   
   score_and_cart peak_search_simple() const; 

  score_and_cart
    peak_search_distance(coot::Cartesian previous_atom, 
			 coot::Cartesian point) const;


   score_and_cart
   peak_search_distance_theta_2(coot::Cartesian ith_point,
				coot::Cartesian ith_plus_one_point, 
				coot::Cartesian point //trial_centre_point
				) const; 
   score_and_cart
   peak_search_distance_theta_2(const TreeNode *node) const;

   score_and_cart
   peak_search_distance_angle_torsion(const TreeNode *node) const; 


   score_and_cart 
   old_peak_search_wrapper(coot::Cartesian point, int ith_res) const; 

  float deviation_from_ideal_length_score(float length) const; 

  vector<coot::Cartesian> branch_pts() const { return branch_points; } 

  vector<coot::Cartesian> find_branch_points(const clipper::Xmap<float> &map, 
				       const clipper::Xmap<int> &l1,
				       float cut_off); 
  
  short isSmallTriangle(const clipper::Xmap<int> &l1,
			const clipper::Xmap<float> &map,
			float cut_off, 
			const clipper::Skeleton_basic::Neighbours &fd_neighb,
			const clipper::Skeleton_basic::Neighbours &edge_neighb,
			const clipper::Coord_grid &pos) const; 
  
//   small_triangle_thing
//   isSmallTriangle_new(const clipper::Xmap<int> &l1,
// 		      const clipper::Xmap<float> &map,
// 		      float cut_off, 
// 		      const clipper::Skeleton_basic::Neighbours &fd_neighb,
// 		      const clipper::Skeleton_basic::Neighbours &edge_neighb,
// 		      const clipper::Coord_grid &pos,
// 		      const vector<coot::Cartesian> &previous,
// 		      vector<clipper::Coord_grid> *vcg_running) const; 
   
  small_triangle_thing
  isSmallTriangle_new(const clipper::Xmap<int> &l1,
		      const clipper::Xmap<float> &map,
		      float cut_off, 
		      const clipper::Skeleton_basic::Neighbours &fd_neighb,
		      const clipper::Skeleton_basic::Neighbours &edge_neighb,
		      const clipper::Coord_grid &pos,
		      const vector<coot::Cartesian> &previous) const; 
   
   
   coot::Cartesian SmallTriangle_to_branch_point(const clipper::Xmap<int> &l1,
					   const clipper::Skeleton_basic::Neighbours &fd_neighb,
					   const clipper::Coord_grid &pos) const; 
   
   float density_at_point(coot::Cartesian trial_point) const; 
   
   float mid_point_density_score(coot::Cartesian prev, 
				 coot::Cartesian trial) const; 
      
   float mid_points_density_score(coot::Cartesian prev, 
				  coot::Cartesian trial) const; 
      
   void ca_grow_recursive(); 

   score_and_cart build_first_recursive(); 

   score_and_cart build_first_cheat(); 

   atom_selection_container_t grown_Cas() const; 
   
   float branch_point_proximity_score(coot::Cartesian trial_point) const;
   
   float prebuilt_exclusion_score(coot::Cartesian trial_point) const; 
   
   float segment_score(const clipper::Coord_grid &c_g_point, 
		       const clipper::Coord_grid &c_g_previous_atom) const; 
   
   void symmetry_expand_branch_points(); 
   
   atom_selection_container_t 
   symmetry_expanded_branch_points(const clipper::Xmap<float> &map) const
   {return convert_to_atoms(map, 
			    branch_points_symm_expanded, 
			    "symmetry expanded branch_points");}; 
   
   scores non_angle_micro_point_score(coot::Cartesian previous_atom, 
				      coot::Cartesian trial_point) const; 
   
   float theta_2_score(coot::Cartesian ith_plus_1_point,
		       coot::Cartesian ith_plus_2_point, 
		       coot::Cartesian ith_plus_3_point) const; 

   float angle_torsion_score(const TreeNode *node) const; 

   void check_angle_torsion(atom_selection_container_t asc) const; 


   
   
   float interconnectedness(int ntips) const; 

   float maximum_gridding(const clipper::Xmap<float> &map) const; 
   float maximum_gridding() const; 

   score_and_cart recursive_build(int ith_res, int depth); 

   score_and_cart recursive_build(const TreeNode *node, int ith_res, int depth); 

   // vector<coot::Cartesian_and_Grid> fitting_targets(int ires, float gridding) const; 

   vector<coot::Cartesian_and_Grid> fitting_targets(const TreeNode *node, 
					      float max_gridding) const; 


   score_and_cart
   peak_search_wrapper(const TreeNode *node, int ith_res, int depth); 

   void make_tree_node_map(); 



   short int 
   depth_search_skeleton_internal(const clipper::Coord_grid &current,
				  const clipper::Coord_grid &previous,
				  const clipper::Coord_grid &prev_prev,
				  const clipper::Coord_grid &target, 
				  int depth, int length) const; 

   short int
   depth_search_skeleton(const clipper::Coord_grid &start, 
			 const clipper::Coord_grid &target) const; 
      
   void depth_search_skeleton_testing();

   void depth_search_skeleton_testing_2();

   void show_segment_map(); 

   void transfer_segment_map(clipper::Xmap<int> *skel_p) const; 

   void export_coordinates(atom_selection_container_t asc, 
			   std::string filename) const; 

};



