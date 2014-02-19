

#include "molecule-class-info.h"

// shift "bottom left" to the origin and make sure that it's on a grid that is
// "nice" (acceptable?) for molrep.
//
// Now the resulting map has gompertz scaling of density (drops off
// towards the edge) and border added (8% of radius).
// 
int
molecule_class_info_t::export_map_fragment_with_origin_shift(float radius,
							     clipper::Coord_orth centre,
							     const std::string &file_name) const {

   // centre is the centre of the map of this molecule (that we are extracting from)

   int r = 0;
   if (has_map()) {
      clipper::Xmap<float>  &xmap = xmap_list[0];
      clipper::Cell          xmap_cell = xmap.cell();
      clipper::Grid_sampling xmap_grid_sampling = xmap.grid_sampling();
      clipper::Coord_orth centre_moved = centre_moved - clipper::Coord_orth(0.1, 0.1, 0.1);

      clipper::Grid_range gr0(xmap_cell, xmap_grid_sampling, radius);
      clipper::Grid_range gr1(gr0.min() + centre.coord_frac(xmap_cell).coord_grid(xmap_grid_sampling),
			      gr0.max() + centre.coord_frac(xmap_cell).coord_grid(xmap_grid_sampling));
      
      int nx_u = 2*gr0.max().u();
      int nx_v = 2*gr0.max().v();
      int nx_w = 2*gr0.max().w();
      
      clipper::Coord_grid nxmap_grid(nx_u, nx_v, nx_w);
      clipper::Coord_grid nxmap_grid_max(nx_u-1, nx_v-1, nx_w-1);
      clipper::Coord_grid nxmap_origin(0,0,0);
      clipper::Grid_range nxmap_grid_range(nxmap_origin, nxmap_grid_max);
      clipper::Grid_sampling nxmap_grid_sampling(nx_u, nx_v, nx_w);
      std::cout << "--------------- nxmap grid_sampling init with " << nx_u << " " << nx_v << " " << nx_w
		<< std::endl;
      
      std::cout << "--------------- gr0.min() " << gr0.min().format() << std::endl;
      std::cout << "--------------- gr0.max() " << gr0.max().format() << std::endl;
      std::cout << "--------------- nxmap_grid_sampling "   << nxmap_grid_sampling.format() << std::endl;
      std::cout << "--------------- nxmap grid range min: " << nxmap_grid_range.min().format() << std::endl;
      std::cout << "--------------- nxmap grid range max: " << nxmap_grid_range.max().format() << std::endl;
      

      double nx_alpha = xmap.cell().descr().alpha();
      double nx_beta  = xmap.cell().descr().beta();
      double nx_gamma = xmap.cell().descr().gamma();
      double nx_a = xmap.cell().descr().a() * double(nxmap_grid_sampling.nu())/double(xmap_grid_sampling.nu());
      double nx_b = xmap.cell().descr().b() * double(nxmap_grid_sampling.nv())/double(xmap_grid_sampling.nv());
      double nx_c = xmap.cell().descr().c() * double(nxmap_grid_sampling.nw())/double(xmap_grid_sampling.nw());
      clipper::Cell_descr nxmap_cell_descr(nx_a, nx_b, nx_c, nx_alpha, nx_beta, nx_beta);
      clipper::Cell nxmap_cell(nxmap_cell_descr);

      // init nxmap
      clipper::NXmap<float> nxmap(nxmap_cell, nxmap_grid_sampling, nxmap_grid_range);
      
      clipper::Xmap<float>::Map_reference_coord ix(xmap);
      clipper::Coord_orth centre_radius(centre.x() - radius,
					centre.y() - radius,
					centre.z() - radius);
      clipper::Coord_orth nxmap_centre(radius, radius, radius);
      clipper::Coord_grid offset = xmap.coord_map(centre_radius).coord_grid();
      
      std::cout << "--------------- xmap offset to centre        " << centre_radius.format() << std::endl;
      std::cout << "--------------- xmap offset to centre (grid) " << offset.format() << std::endl;

      typedef clipper::NXmap<float>::Map_reference_index NRI;
      double limited_radius = radius * 0.92;
      for (NRI inx = nxmap.first(); !inx.last(); inx.next()) {
	 clipper::Coord_orth p = inx.coord().coord_frac(nxmap_grid_sampling).coord_orth(nxmap_cell);
	 double d_to_c_sq = clipper::Coord_orth(p-nxmap_centre).lengthsq();
	 if (d_to_c_sq > limited_radius*limited_radius) {
	    nxmap[inx] = 0.0;
	    if (0) 
	       std::cout << " inx " << inx.coord().format() << " " << d_to_c_sq << "  " << p.format() << " "
			 << centre.format() << " vs " << limited_radius*limited_radius << " is outside "
			 << std::endl;
	 } else {
	    // ix indexes the xmap
	    clipper::Coord_grid gp = p.coord_frac(xmap_cell).coord_grid(xmap_grid_sampling);
	    ix.set_coord(gp + offset);

	    // make a function that is y=1 around x=0 and y=1 around x=1 and
	    // falls of to 0.5 around x=0.8 or so.
	    //
	    double x = sqrt(d_to_c_sq)/limited_radius;
	    double gompertz_a = 0.14;
	    double gompertz_b = 0.1; 
	    double gompertz_c = 3;
	    double gompertz_scale = 1 - (-gompertz_a*1.1 + gompertz_a * exp (gompertz_b * exp(gompertz_c * x)));
	    nxmap[inx] = xmap[ix] * gompertz_scale;
	    if (0)
	       std::cout << " inx " << inx.coord().format() << " " << d_to_c_sq << "  " << p.format() << " " << centre.format()
			 << " vs " << limited_radius*limited_radius << " " << gompertz_scale << std::endl;
	 } 
	    
      }
      clipper::CCP4MAPfile mapout;
      mapout.open_write(file_name);
      mapout.set_cell(nxmap_cell);
      mapout.export_nxmap(nxmap);
      mapout.close_write();
      std::cout << "Exported map " << file_name << std::endl;
   }

   return r;
}
