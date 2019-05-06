
#include "coot-coord-utils.hh"

int
coot::get_brick_id_inner(int x_idx, int y_idx, int z_idx,
			 int nx_grid, int ny_grid, int nz_grid) {

   // This is a non-intuitive function

   // for a 6x4x..
   //
   // 26   27    34  35   42 43
   // 24   25    32  33   40 41
   //  2    3    10  11   18 19
   //  0    1     8   9   16 17

   // and the layer (in z) above that:
   //
   // 30   31    38  39   46 47
   // 28   29    36  37   44 44
   //  6    7    14  15   22 23
   //  4    5    12  13   20 21

   int even_x = (x_idx%2==0);
   int even_y = (y_idx%2==0);
   int even_z = (z_idx%2==0);

//    // convert 2d layout to 1D
//    //
//    // act as if even, then correct
//    //
//    int id_1 = 2 * x_idx;
//    if (!even_x) id_1 = 2*(x_idx-1)+1;

//    int id_2 = y_idx * nx_grid; // ! confusing(?)
//    if (!even_y) id_2 = (y_idx-1) * nx_grid + 2;

   // act as if even, then correct
   //
   int id_1 = 4 * x_idx;
   if (!even_x) id_1 = 4*(x_idx-1)+1;
   int id_2 = y_idx * nx_grid * 2;
   if (!even_y) id_2 = (y_idx - 1) * nx_grid * 2 + 2;
   int id_3 = z_idx * nx_grid * ny_grid;
   if (!even_z) id_3 = (z_idx-1) * nx_grid * ny_grid + 4;

   int id_sum = id_1 + id_2 + id_3;

   // convert 3d layout to 1D

   if (false)
      std::cout << "debug:: get_brick_id_inner() "
		<< x_idx << " " << y_idx << " " << z_idx << " "
		<< " -> even-x even-y "
		<< even_x << " " << even_y << " id_1, id_2: " << id_1 << " " << id_2
		<< " -> id_sum " << id_sum << std::endl;

   return id_sum;

}

int
coot::get_brick_id(const clipper::Coord_orth &pt, const clipper::Coord_orth &pt_minimums,
		   int nx_grid, int ny_grid, int nz_grid, float brick_length) {

   clipper::Coord_orth delta = pt - pt_minimums;
   int x_idx(delta.x()/brick_length);
   int y_idx(delta.y()/brick_length);
   int z_idx(delta.z()/brick_length);

   return get_brick_id_inner(x_idx, y_idx, z_idx, nx_grid, ny_grid, nz_grid);
}

// split a molecule into "bricks" - cubic sets of atom indices that don't overlap
// return vector needs to be multiples of 8 (8 cubes will coverer all space)
// If the atom max radius is 3A, then the brick should have length 6A.
// brick-id: 0,8,16,24 (etc) are done at the same time
// then
// brick-id: 1,9,17,25 (etc) are done at the same time
// then
// brick-id: 2,10,18,26 (etc) are done at the same time.
//
std::vector<std::vector<int> >
coot::molecule_to_bricks(mmdb::Manager *mol, int SelectionHandle,
			       float atom_max_radius) {

   std::vector<std::vector<int> > v;
   float brick_length = atom_max_radius * 2.0;

   std::pair<clipper::Coord_orth, clipper::Coord_orth> e = util::extents(mol, SelectionHandle);
   clipper::Coord_orth pt_minimums = e.first;

   clipper::Coord_orth delta = e.second - e.first;
   int n_x(delta.x()/brick_length);
   int n_y(delta.x()/brick_length);
   int n_z(delta.x()/brick_length);

   // make them all even
   if (n_x%2 != 0) n_x++;
   if (n_y%2 != 0) n_y++;
   if (n_z%2 != 0) n_z++;

   std::cout << "----------- here with brick dimension " << n_x << " " << n_y << " " << n_z
	     << std::endl;

   try {
      v.reserve(n_x*n_y*n_z);
      mmdb::Atom **selected_atoms = 0;
      int n_selected_atoms = 0;
      mol->GetSelIndex(SelectionHandle, selected_atoms, n_selected_atoms);
      for (int iat=0; iat<n_selected_atoms; iat++) {
	 mmdb::Atom *at = selected_atoms[iat];
	 clipper::Coord_orth pt = co(at);
	 unsigned int brick_idx = get_brick_id(pt, pt_minimums, n_x, n_y, n_z, brick_length);
	 try {
	    if (brick_idx >= v.size()) {
	       v.resize(brick_idx+1);
	       v[brick_idx].reserve(20);
	    }
	    v[brick_idx].push_back(iat);
	 }
	 catch (...) {
	    std::cout << "caught an default exception for brick idx " << brick_idx << std::endl;
	 }
      }
   }
   catch (...) {
      std::cout << "caught an default exception on reserve " << std::endl;
   }
   return v;
}
