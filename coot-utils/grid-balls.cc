
#include "grid-balls.hh"

std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::grid_balls_t::get_extents(mmdb::Manager *mol) const {

   float mol_x_min =  1.0e30;
   float mol_y_min =  1.0e30;
   float mol_z_min =  1.0e30;
   float mol_x_max = -1.0e30;
   float mol_y_max = -1.0e30;
   float mol_z_max = -1.0e30;

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        if (at->x < mol_x_min) mol_x_min = at->x;
                        if (at->y < mol_y_min) mol_y_min = at->y;
                        if (at->z < mol_z_min) mol_z_min = at->z;
                        if (at->x > mol_x_max) mol_x_max = at->x;
                        if (at->y > mol_y_max) mol_y_max = at->y;
                        if (at->z > mol_z_max) mol_z_max = at->z;
                     }
                  }
               }
            }
         }
      }
   }

   // test here for sane min and max values

   return std::make_pair(clipper::Coord_orth(mol_x_min, mol_y_min, mol_z_min),
                         clipper::Coord_orth(mol_x_max, mol_y_max, mol_z_max));

}

int
coot::grid_balls_t::grid_index(const triple_index_t &t) const {
   return t.ix + t.iy * nx + t.iz * nx * ny;
}

coot::grid_balls_t::triple_index_t
coot::grid_balls_t::deindex(int i) const {
   triple_index_t t;
   t.ix = i % nx;
   t.iy = (i / nx) % ny;
   t.iz = i / (nx * ny);
   return t;
}

coot::grid_balls_t::grid_balls_t(mmdb::Manager *mol, float small_ball_radius, float big_ball_radius) {

   n_grids_per_angstrom = 1.0; //testing

   // get extents of molecule
   std::pair<clipper::Coord_orth, clipper::Coord_orth> ee = get_extents(mol);

   // calculate grid size
   float mol_x_min = ee.first.x();
   float mol_y_min = ee.first.y();
   float mol_z_min = ee.first.z();
   float mol_x_max = ee.second.x();
   float mol_y_max = ee.second.y();
   float mol_z_max = ee.second.z();

   float extra_extents = 5.0; // angstroms

   grid_min_x = mol_x_min - extra_extents;
   grid_min_y = mol_y_min - extra_extents;
   grid_min_z = mol_z_min - extra_extents;
   grid_max_x = mol_x_max + extra_extents;
   grid_max_y = mol_y_max + extra_extents;
   grid_max_z = mol_z_max + extra_extents;

   nx = int((grid_max_x - grid_min_x + extra_extents) * n_grids_per_angstrom) + 1;
   ny = int((grid_max_y - grid_min_y + extra_extents) * n_grids_per_angstrom) + 1;
   nz = int((grid_max_z - grid_min_z + extra_extents) * n_grids_per_angstrom) + 1;

   grid.resize(nx * ny * nz);

   test_grid();
   test_index_deindex();

   brick_the_model(mol);

}

coot::grid_balls_t::triple_index_t
coot::grid_balls_t::mol_space_to_grid_point(const point_3d_t &p) const {

   triple_index_t t;
   float x = p.x - grid_min_x;
   float y = p.y - grid_min_y;
   float z = p.z - grid_min_z;

   t.ix = float(x * n_grids_per_angstrom - 0);
   t.iy = float(y * n_grids_per_angstrom - 0);
   t.iz = float(z * n_grids_per_angstrom - 0);

   if (false) {
      // std::cout << "grid_min_x " << grid_min_x << " grid_min_y " << grid_min_y << " grid_min_z " << grid_min_z << std::endl;
      std::cout << "p.x " << p.x << " p.y " << p.y << " p.z " << p.z << std::endl;
      std::cout << "x " << x << " y " << y << " z " << z << std::endl;
      std::cout << "t.ix " << t.ix << " t.iy " << t.iy << " t.iz " << t.iz << std::endl;
   }

   return t;

}

coot::grid_balls_t::point_3d_t
coot::grid_balls_t::grid_point_to_mol_space(const triple_index_t &t) const {

   point_3d_t p;
   p.x = float(t.ix+0) / n_grids_per_angstrom + grid_min_x;
   p.y = float(t.iy+0) / n_grids_per_angstrom + grid_min_y;
   p.z = float(t.iz+0) / n_grids_per_angstrom + grid_min_z;

   return p;
}


void
coot::grid_balls_t::test_grid() const {

   std::cout << "testing grid to space..." << std::endl;
   int n_correct = 0;
   int n_wrong = 0;

   for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
         for (int iz=0; iz<nx; iz++) {
            point_3d_t p = grid_point_to_mol_space(triple_index_t(ix, iy, iz));
            triple_index_t t = mol_space_to_grid_point(p);
            int tix_as_int =  int(round(t.ix));
            int tiy_as_int =  int(round(t.iy));
            int tiz_as_int =  int(round(t.iz));
            if (ix != tix_as_int) {
               std::cout << "Error in grid indexing X: input: " << ix << " " << iy << " " << iz
                         << " as_int: " << tix_as_int << " " << tiy_as_int << " " << tiz_as_int
                         << " result: " << t.ix << " " << t.iy << " " << t.iz
                         << "\n";
               n_wrong++;
            }
            if (iy != tiy_as_int) {
               std::cout << "Error in grid indexing Y: input: " << ix << " " << iy << " " << iz
                         << " as_int: " << tix_as_int << " " << tiy_as_int << " " << tiz_as_int
                         << " result: " << t.ix << " " << t.iy << " " << t.iz
                         << "\n";
               n_wrong++;
            }
            if (iz != tiz_as_int) {
               std::cout << "Error in grid indexing Z: input " << ix << " " << iy << " " << iz
                         << " as_int: " << tix_as_int << " " << tiy_as_int << " " << tiz_as_int
                         << " result: " << t.ix << " " << t.iy << " " << t.iz
                         << "\n";
               n_wrong++;
            }


            if (t.ix != ix || t.iy != iy || t.iz != iz) {
            } else {
               n_correct++;
            }
         }
      }
   }
   int n_total = n_correct + n_wrong;
   std::cout << "testing done. n_correct: " << n_correct << " n_wrong " << n_wrong
             << "  " << 100.0 * static_cast<float>(n_wrong)/static_cast<float>(n_total) << " %" << std::endl;

}


void
coot::grid_balls_t::test_index_deindex() const {

   std::cout << "testing index/deindex..." << std::endl;
   int n_correct = 0;
   int n_wrong = 0;

   for (int i=0; i<nx*ny*nz; i++) {
      triple_index_t t = deindex(i);
      int i_as_int = grid_index(t);
      if (i != i_as_int) {
         std::cout << "Error in index/deindex: input: " << i << " as_int: " << i_as_int
                   << " result: " << t.ix << " " << t.iy << " " << t.iz
                   << "\n";
         n_wrong++;
      } else {
         n_correct++;
      }
   }
   int n_total = n_correct + n_wrong;
   std::cout << "testing for index/deindex done: n_correct: " << n_correct << " n_wrong " << n_wrong
             << "  " << 100.0 * static_cast<float>(n_wrong)/static_cast<float>(n_total) << " %" << std::endl;

}


void
coot::grid_balls_t::brick_the_model(mmdb::Manager *mol) {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        point_3d_t p(at->x, at->y, at->z);

                     }
                  }
               }
            }
         }
      }
   }

}
