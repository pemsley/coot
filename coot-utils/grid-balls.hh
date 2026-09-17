
#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

#include "geometry/protein-geometry.hh"


namespace coot {

   class grid_balls_t {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> get_extents(mmdb::Manager *mol) const;
   public:
      // 20260722-PE we pass the imol so that the correct ligand geometry can be looked up.
      grid_balls_t(int imol, mmdb::Manager *mol, const protein_geometry *geom_p,
                   float small_ball_radius, float big_ball_radius);
      class triple_index_t {
      public:
         int ix, iy, iz;
         triple_index_t(const int &xx, const int &yy, const int &zz) { ix = xx; iy = yy; iz = zz; }
         triple_index_t() { ix = 0; iy = 0; iz = 0; }
      };
      class point_3d_t { // float class (unlike coord_orth)
      public:
         point_3d_t(const float &xx, const float &yy, const float &zz) { x = xx; y = yy; z = zz; }
         point_3d_t() { x = 0.0f; y = 0.0f; z = 0.0f; }
         float x, y, z;
      };
      class grid_point_t {
      public:
         float v;
         grid_point_t(const float &vv) { v = vv; }
         grid_point_t() { v = 0.0; }
         std::vector<mmdb::Atom *> atoms;
      };
      class cavity_t { // a distinct cavity (or subpocket): a cluster of cavity points
      public:
         int id;                        // 1-based, ranked by volume (largest first)
         int parent_id;                 // 0 for a top-level cavity; else the parent cavity id
         float volume;                  // A^3
         float surface_area;            // A^2 (exposed voxel faces)
         point_3d_t centre;             // centroid, mol space
         std::vector<int> grid_indices; // the cavity's grid points
         std::vector<mmdb::Residue *> lining_residues; // residues lining the cavity
         cavity_t() { id = 0; parent_id = 0; volume = 0.0f; surface_area = 0.0f; }
      };
      int imol;
      mmdb::Manager *mol; // not owned
      const protein_geometry *geom_p;
      int nx, ny, nz;
      std::vector<grid_point_t> grid;
      float n_grids_per_angstrom;
      float probe_in_radius;  // small ball (KVFinder Probe In,  ~1.4A)
      float probe_out_radius; // big ball   (KVFinder Probe Out, ~4A)
      std::vector<unsigned char> cavity_grid; // 1 = cavity point; filled by find_cavities()
      float grid_min_x;
      float grid_min_y;
      float grid_min_z;
      float grid_max_x;
      float grid_max_y;
      float grid_max_z;
      void test_grid() const; // test that position transformation correct
      void test_index_deindex() const; //  test that I have the indexing and deindexing correct
      triple_index_t mol_space_to_grid_point(const point_3d_t &p) const;
      point_3d_t grid_point_to_mol_space(const triple_index_t &t) const;
      int grid_index(const triple_index_t &t) const;
      triple_index_t deindex(int idx) const;
      void brick_the_model(mmdb::Manager *mol);
      // two-probe (KVFinder) cavity detection:
      std::vector<unsigned char> make_blocked_grid(float probe_radius) const;
      std::vector<unsigned char> flood_fill_exterior(const std::vector<unsigned char> &blocked) const;
      std::vector<unsigned char> make_accessible_grid(float probe_radius) const;
      void find_cavities();
      std::vector<cavity_t> cavities; // filled by cluster_cavities()
      // cluster the cavity points into distinct cavities, with per-cavity volume,
      // dropping any below volume_min (A^3) and ranking the rest by volume.
      void cluster_cavities(float volume_min = 0.0);
      // write a plain table of cavity point positions (cavity_id x y z), one point
      // per line, with a per-cavity comment header. Intended for a script that adds
      // the points to a generic display object via to_generic_object_add_points().
      void write_cavity_points(const std::string &file_name) const;
      // for each cavity, find the residues with an atom within contact_dist (A) of
      // any of its points, and store them (sorted) in cavity_t::lining_residues.
      void compute_lining_residues(float contact_dist = 5.0);
      void compute_lining_residues_for(std::vector<cavity_t> &cavs, float contact_dist);

      // subpocket segmentation (distance-transform watershed):
      std::vector<cavity_t> subpockets; // filled by segment_cavities()
      // split one cavity into subpockets; returns them (empty if it does not split).
      // min_prominence (A): basins whose separating neck is shallower than this (a
      // low saddle relative to the peak) are merged - controls over-segmentation.
      std::vector<cavity_t> segment_one_cavity(const cavity_t &c, int min_seed_depth,
                                               float merge_dist, float min_prominence) const;
      // segment every cavity above min_volume_to_segment (A^3) into subpockets.
      void segment_cavities(float min_volume_to_segment = 100.0, int min_seed_depth = 1,
                            float merge_dist = 4.0, float min_prominence = 0.9);
      void write_subpocket_points(const std::string &file_name) const;
   };
}

