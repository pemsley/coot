#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

namespace coot {

   class grid_balls_t {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> get_extents(mmdb::Manager *mol) const;
   public:
      grid_balls_t(mmdb::Manager *mol, float small_ball_radius, float big_ball_radius);
      class triple_index_t {
      public:
         // int ix,iy, iz;
         float ix,iy, iz;
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
      int nx, ny, nz;
      std::vector<grid_point_t> grid;
      float n_grids_per_angstrom;
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
   };
}

