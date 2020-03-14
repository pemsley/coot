// Header file for shift field refinement
/* Copyright 2018 Kevin Cowtan & University of York all rights reserved */


#include "shiftfield.h"

bool Shift_field_refine::shift_field_coord(
  const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
  clipper::Xmap<float>& x1map,      clipper::Xmap<float>& x2map,      clipper::Xmap<float>& x3map,
  float rad, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;
  const clipper::Xmap<float>&   mmap = mask;

  { // being map preparation
    clipper::FFTmap_p1 cfftx(grid), cffty(grid), cfftz(grid);
    clipper::Xmap<float>::Map_reference_coord i0( cmap, clipper::Coord_grid(0,0,0) );
    clipper::Xmap<float>::Map_reference_coord iu, iv, iw;

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	cfftx.real_data( iw.coord() ) = cmap[iw];
        }

    // calculate map coefficients
    cfftx.fft_x_to_h( cmap.cell().volume() );

    // calculate gradient map coefficients
    const clipper::Grid_sampling& g = cfftx.grid_real();
    const clipper::Grid&         gr = cfftx.grid_reci();
    clipper::Coord_grid ch( g.nu()/2, g.nv()/2, g.nw()/2 );
    clipper::Coord_grid c;
    std::complex<float> i(0.0,1.0);
    for ( c.u() = 0; c.u() < gr.nu(); c.u()++ )
      for ( c.v() = 0; c.v() < gr.nv(); c.v()++ )
        for ( c.w() = 0; c.w() < gr.nw(); c.w()++ ) {
          const clipper::HKL hkl = clipper::HKL( ( c + ch ).unit( g ) - ch );
          const std::complex<float> cdata = i * cfftx.cplx_data(c);
          cfftx.cplx_data(c) = float(clipper::Util::twopi()*hkl.h()) * cdata;
          cffty.cplx_data(c) = float(clipper::Util::twopi()*hkl.k()) * cdata;
          cfftz.cplx_data(c) = float(clipper::Util::twopi()*hkl.l()) * cdata;
        }

    // calculate gradient maps
    cfftx.fft_h_to_x( 1.0 / cmap.cell().volume() );
    cffty.fft_h_to_x( 1.0 / cmap.cell().volume() );
    cfftz.fft_h_to_x( 1.0 / cmap.cell().volume() );

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	x1map[iw] = cfftx.real_data( iw.coord() );
        	x2map[iw] = cffty.real_data( iw.coord() );
        	x3map[iw] = cfftz.real_data( iw.coord() );
        }

    /*
    // output intermediate maps
    clipper::CCP4MAPfile mapout;
    mapout.open_write( "grad1.map" );
    mapout.export_xmap( x1map );
    mapout.close_write();
    mapout.open_write( "grad2.map" );
    mapout.export_xmap( x2map );
    mapout.close_write();
    mapout.open_write( "grad3.map" );
    mapout.export_xmap( x3map );
    mapout.close_write();
    */

  } // end map preparation

  // make xmap
  clipper::Xmap<float> y1map( spgr, cell, grid ); 
  clipper::Xmap<float> y2map( spgr, cell, grid ); 
  clipper::Xmap<float> y3map( spgr, cell, grid );  
  clipper::Xmap<float> x11map( spgr, cell, grid );
  clipper::Xmap<float> x12map( spgr, cell, grid );
  clipper::Xmap<float> x13map( spgr, cell, grid );
  clipper::Xmap<float> x22map( spgr, cell, grid );
  clipper::Xmap<float> x23map( spgr, cell, grid );
  clipper::Xmap<float> x33map( spgr, cell, grid );

  // calculate the term XTY and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y2map[ix]  = x2map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y3map[ix]  = x3map[ix]*ymap[ix]*mmap[ix];

  // calculate the term XTX and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x12map[ix] = x1map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x13map[ix] = x1map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x22map[ix] = x2map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x23map[ix] = x2map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x33map[ix] = x3map[ix]*x3map[ix]*mmap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // apply the convolution filter g to each map
  fltr( y1map, y1map );
  fltr( y2map, y2map );
  fltr( y3map, y3map );
  fltr( x11map, x11map );
  fltr( x12map, x12map );
  fltr( x13map, x13map );
  fltr( x22map, x22map );
  fltr( x23map, x23map );
  fltr( x33map, x33map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) {

    std::vector<float> v(3);
    v[0] = y1map[ix];
    v[1] = y2map[ix];
    v[2] = y3map[ix];

    clipper::Matrix<float> m(3,3);
    m(0,0) = x11map[ix];
    m(0,1) = m(1,0) = x12map[ix];
    m(0,2) = m(2,0) = x13map[ix];
    m(1,1) = x22map[ix];
    m(1,2) = m(2,1) = x23map[ix];
    m(2,2) = x33map[ix];

    v = m.solve(v);

    x1map[ix] = v[0];
    x2map[ix] = v[1];
    x3map[ix] = v[2];
  }

  return true;
}



bool Shift_field_refine::shift_field_coord_const( 
  const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
  clipper::Xmap<float>& x1map,      clipper::Xmap<float>& x2map,      clipper::Xmap<float>& x3map,
  float rad, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;
  const clipper::Xmap<float>&   mmap = mask;

  { // being map preparation
    clipper::FFTmap_p1 cfftx(grid), cffty(grid), cfftz(grid);
    clipper::Xmap<float>::Map_reference_coord i0( cmap, clipper::Coord_grid(0,0,0) );
    clipper::Xmap<float>::Map_reference_coord iu, iv, iw;

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	cfftx.real_data( iw.coord() ) = cmap[iw];
        }

    // calculate map coefficients
    cfftx.fft_x_to_h( cmap.cell().volume() );

    // calculate gradient map coefficients
    const clipper::Grid_sampling& g = cfftx.grid_real();
    const clipper::Grid&         gr = cfftx.grid_reci();
    clipper::Coord_grid ch( g.nu()/2, g.nv()/2, g.nw()/2 );
    clipper::Coord_grid c;
    std::complex<float> i(0.0,1.0);
    for ( c.u() = 0; c.u() < gr.nu(); c.u()++ )
      for ( c.v() = 0; c.v() < gr.nv(); c.v()++ )
        for ( c.w() = 0; c.w() < gr.nw(); c.w()++ ) {
          const clipper::HKL hkl = clipper::HKL( ( c + ch ).unit( g ) - ch );
          const std::complex<float> cdata = i * cfftx.cplx_data(c);
          cfftx.cplx_data(c) = float(clipper::Util::twopi()*hkl.h()) * cdata;
          cffty.cplx_data(c) = float(clipper::Util::twopi()*hkl.k()) * cdata;
          cfftz.cplx_data(c) = float(clipper::Util::twopi()*hkl.l()) * cdata;
        }

    // calculate gradient maps
    cfftx.fft_h_to_x( 1.0 / cmap.cell().volume() );
    cffty.fft_h_to_x( 1.0 / cmap.cell().volume() );
    cfftz.fft_h_to_x( 1.0 / cmap.cell().volume() );

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	x1map[iw] = cfftx.real_data( iw.coord() );
        	x2map[iw] = cffty.real_data( iw.coord() );
        	x3map[iw] = cfftz.real_data( iw.coord() );
        }

    /*
    // output intermediate maps
    clipper::CCP4MAPfile mapout;
    mapout.open_write( "grad1.map" );
    mapout.export_xmap( x1map );
    mapout.close_write();
    mapout.open_write( "grad2.map" );
    mapout.export_xmap( x2map );
    mapout.close_write();
    mapout.open_write( "grad3.map" );
    mapout.export_xmap( x3map );
    mapout.close_write();
    */

  } // end map preparation

  // make xmap
  clipper::Xmap<float> y0map( spgr, cell, grid );
  clipper::Xmap<float> y1map( spgr, cell, grid ); 
  clipper::Xmap<float> y2map( spgr, cell, grid ); 
  clipper::Xmap<float> y3map( spgr, cell, grid );  
  clipper::Xmap<float> x00map( spgr, cell, grid );
  clipper::Xmap<float> x01map( spgr, cell, grid );
  clipper::Xmap<float> x02map( spgr, cell, grid );
  clipper::Xmap<float> x03map( spgr, cell, grid );
  clipper::Xmap<float> x11map( spgr, cell, grid );
  clipper::Xmap<float> x12map( spgr, cell, grid );
  clipper::Xmap<float> x13map( spgr, cell, grid );
  clipper::Xmap<float> x22map( spgr, cell, grid );
  clipper::Xmap<float> x23map( spgr, cell, grid );
  clipper::Xmap<float> x33map( spgr, cell, grid );

  // calculate the term XTY and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y0map[ix]  = ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y2map[ix]  = x2map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y3map[ix]  = x3map[ix]*ymap[ix]*mmap[ix];

  // calculate the term XTX and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x00map[ix] = mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x01map[ix] = x1map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x02map[ix] = x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x03map[ix] = x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x12map[ix] = x1map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x13map[ix] = x1map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x22map[ix] = x2map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x23map[ix] = x2map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x33map[ix] = x3map[ix]*x3map[ix]*mmap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // apply the convolution filter g to each map
  fltr( y0map, y0map );
  fltr( y1map, y1map );
  fltr( y2map, y2map );
  fltr( y3map, y3map );
  fltr( x00map, x00map );
  fltr( x01map, x01map );
  fltr( x02map, x02map );
  fltr( x03map, x03map );
  fltr( x11map, x11map );
  fltr( x12map, x12map );
  fltr( x13map, x13map );
  fltr( x22map, x22map );
  fltr( x23map, x23map );
  fltr( x33map, x33map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) {

    std::vector<float> v(4);
    v[0] = y0map[ix];
    v[1] = y1map[ix];
    v[2] = y2map[ix];
    v[3] = y3map[ix];

    clipper::Matrix<float> m(4,4);
    m(0,0) = x00map[ix];
    m(0,1) = m(1,0) = x01map[ix];
    m(0,2) = m(2,0) = x02map[ix];
    m(0,3) = m(3,0) = x03map[ix];
    m(1,1) = x11map[ix];
    m(1,2) = m(2,1) = x12map[ix];
    m(1,3) = m(3,1) = x13map[ix];
    m(2,2) = x22map[ix];
    m(2,3) = m(3,2) = x23map[ix];
    m(3,3) = x33map[ix];

    v = m.solve(v);

    x1map[ix] = v[1];
    x2map[ix] = v[2];
    x3map[ix] = v[3];
  }

  return true;
}



bool Shift_field_refine::shift_field_u_iso(
  const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
  clipper::Xmap<float>& x1map,
  float rad, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;
  const clipper::Xmap<float>&   mmap = mask;

  { // being map preparation
    clipper::FFTmap_p1 cfftx(grid);
    clipper::Xmap<float>::Map_reference_coord i0( cmap, clipper::Coord_grid(0,0,0) );
    clipper::Xmap<float>::Map_reference_coord iu, iv, iw;

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	cfftx.real_data( iw.coord() ) = cmap[iw];
        }

    // calculate map coefficients
    cfftx.fft_x_to_h( cmap.cell().volume() );

    // calculate gradient map coefficients
    const clipper::Grid_sampling& g = cfftx.grid_real();
    const clipper::Grid&         gr = cfftx.grid_reci();
    clipper::Coord_grid ch( g.nu()/2, g.nv()/2, g.nw()/2 );
    clipper::Coord_grid c;
    for ( c.u() = 0; c.u() < gr.nu(); c.u()++ )
      for ( c.v() = 0; c.v() < gr.nv(); c.v()++ )
        for ( c.w() = 0; c.w() < gr.nw(); c.w()++ ) {
          const clipper::HKL hkl = clipper::HKL( ( c + ch ).unit( g ) - ch );
          const float scl = clipper::Util::twopi2() * hkl.invresolsq(cell);
          cfftx.cplx_data(c) = scl * cfftx.cplx_data(c);
        }

    // calculate gradient maps
    cfftx.fft_h_to_x( 1.0 / cmap.cell().volume() );

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	x1map[iw] = cfftx.real_data( iw.coord() );
        }

    /*
    // output intermediate maps
    clipper::CCP4MAPfile mapout;
    mapout.open_write( "grad1.map" );
    mapout.export_xmap( x1map );
    mapout.close_write();
    */

  } // end map preparation

  // make xmap
  clipper::Xmap<float> y1map( spgr, cell, grid ); 
  clipper::Xmap<float> x11map( spgr, cell, grid );

  // calculate the term XTY and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix]*mmap[ix];

  // calculate the term XTX and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix]*mmap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // apply the convolution filter g to each map
  fltr( y1map, y1map );
  fltr( x11map, x11map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x1map[ix] = y1map[ix] / x11map[ix];

  return true;
}



bool Shift_field_refine::shift_field_u_iso_const(
  const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
  clipper::Xmap<float>& x1map,
  float rad, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;
  const clipper::Xmap<float>&   mmap = mask;

  { // being map preparation
    clipper::FFTmap_p1 cfftx(grid);
    clipper::Xmap<float>::Map_reference_coord i0( cmap, clipper::Coord_grid(0,0,0) );
    clipper::Xmap<float>::Map_reference_coord iu, iv, iw;

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	cfftx.real_data( iw.coord() ) = cmap[iw];
        }

    // calculate map coefficients
    cfftx.fft_x_to_h( cmap.cell().volume() );

    // calculate gradient map coefficients
    const clipper::Grid_sampling& g = cfftx.grid_real();
    const clipper::Grid&         gr = cfftx.grid_reci();
    clipper::Coord_grid ch( g.nu()/2, g.nv()/2, g.nw()/2 );
    clipper::Coord_grid c;
    for ( c.u() = 0; c.u() < gr.nu(); c.u()++ )
      for ( c.v() = 0; c.v() < gr.nv(); c.v()++ )
        for ( c.w() = 0; c.w() < gr.nw(); c.w()++ ) {
          const clipper::HKL hkl = clipper::HKL( ( c + ch ).unit( g ) - ch );
          const float scl = clipper::Util::twopi2() * hkl.invresolsq(cell);
          cfftx.cplx_data(c) = scl * cfftx.cplx_data(c);
        }

    // calculate gradient maps
    cfftx.fft_h_to_x( 1.0 / cmap.cell().volume() );

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	x1map[iw] = cfftx.real_data( iw.coord() );
        }

    /*
    // output intermediate maps
    clipper::CCP4MAPfile mapout;
    mapout.open_write( "grad1.map" );
    mapout.export_xmap( x1map );
    mapout.close_write();
    */

  } // end map preparation

  // make xmap
  clipper::Xmap<float> y0map( spgr, cell, grid );
  clipper::Xmap<float> y1map( spgr, cell, grid ); 
  clipper::Xmap<float> x00map( spgr, cell, grid );
  clipper::Xmap<float> x01map( spgr, cell, grid );
  clipper::Xmap<float> x11map( spgr, cell, grid );

  // calculate the term XTY and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y0map[ix]  = ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix]*mmap[ix];

  // calculate the term XTX and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x00map[ix] = mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x01map[ix] = x1map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix]*mmap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // apply the convolution filter g to each map
  fltr( y0map, y0map );
  fltr( y1map, y1map );
  fltr( x00map, x00map );
  fltr( x01map, x01map );
  fltr( x11map, x11map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) {

    std::vector<float> v(2);
    v[0] = y0map[ix];
    v[1] = y1map[ix];

    clipper::Matrix<float> m(2,2);
    m(0,0) = x00map[ix];
    m(0,1) = m(1,0) = x01map[ix];
    m(1,1) = x11map[ix];

    v = m.solve(v);

    x1map[ix] = v[1];
  }

  return true;
}



bool Shift_field_refine::shift_field_u_aniso(
  const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
	clipper::Xmap<float>& x1map, clipper::Xmap<float>& x2map, clipper::Xmap<float>& x3map,
  clipper::Xmap<float>& x4map, clipper::Xmap<float>& x5map, clipper::Xmap<float>& x6map,
  float rad, clipper::Resolution& reso, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;
  const clipper::Xmap<float>&   mmap = mask;

  // calculate map coefficients
  clipper::HKL_info hkls (spgr, cell, reso, true);
  clipper::HKL_data<clipper::data32::F_phi> cphi( hkls );
  cmap.fft_to(cphi);

  // calculate the gradient map coefficients (by forming a difference quotient in each anisotropic direction)
  clipper::HKL_data<clipper::data32::F_phi> uphi0( hkls ), uphi1( hkls ), uphi2( hkls ), uphi3( hkls ), uphi4( hkls ), uphi5( hkls );

  clipper::U_aniso_orth u0(0.001,0.00,0.00,0.00,0.00,0.00);
  clipper::U_aniso_orth u1(0.00,0.001,0.00,0.00,0.00,0.00);
  clipper::U_aniso_orth u2(0.00,0.00,0.001,0.00,0.00,0.00);
  clipper::U_aniso_orth u3(0.00,0.00,0.00,0.001,0.00,0.00);
  clipper::U_aniso_orth u4(0.00,0.00,0.00,0.00,0.001,0.00);
  clipper::U_aniso_orth u5(0.00,0.00,0.00,0.00,0.00,0.001);

  uphi0.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u0) );
  uphi1.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u1) );
  uphi2.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u2) );
  uphi3.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u3) );
  uphi4.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u4) );
  uphi5.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u5) );

  uphi0 = uphi0 - cphi;
  uphi1 = uphi1 - cphi;
  uphi2 = uphi2 - cphi;
  uphi3 = uphi3 - cphi;
  uphi4 = uphi4 - cphi;
  uphi5 = uphi5 - cphi;

  uphi0.compute( uphi0, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi1.compute( uphi1, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi2.compute( uphi2, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi3.compute( uphi3, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi4.compute( uphi4, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi5.compute( uphi5, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );

  // calculate gradient maps
  // x0map = 1 so we don't need to calculate it
  x1map.fft_from( uphi0 );
  x2map.fft_from( uphi1 );
  x3map.fft_from( uphi2 );
  x4map.fft_from( uphi3 );
  x5map.fft_from( uphi4 );
  x6map.fft_from( uphi5 );

  /*
  // output intermediate maps
  clipper::CCP4MAPfile mapout;
  mapout.open_write( "grad1.map" );
  mapout.export_xmap( x1map );
  mapout.close_write();
  mapout.open_write( "grad2.map" );
  mapout.export_xmap( x2map );
  mapout.close_write();
  mapout.open_write( "grad3.map" );
  mapout.export_xmap( x3map );
  mapout.close_write();
  mapout.open_write( "grad4.map" );
  mapout.export_xmap( x4map );
  mapout.close_write();
  mapout.open_write( "grad5.map" );
  mapout.export_xmap( x5map );
  mapout.close_write();
  mapout.open_write( "grad6.map" );
  mapout.export_xmap( x6map );
  mapout.close_write();
  */

 // make xmap - this allocates the maps for future use
  clipper::Xmap<float> y1map( spgr, cell, grid );
  clipper::Xmap<float> y2map( spgr, cell, grid );
  clipper::Xmap<float> y3map( spgr, cell, grid );
  clipper::Xmap<float> y4map( spgr, cell, grid );
  clipper::Xmap<float> y5map( spgr, cell, grid );
  clipper::Xmap<float> y6map( spgr, cell, grid );
  clipper::Xmap<float> x11map( spgr, cell, grid );
  clipper::Xmap<float> x12map( spgr, cell, grid );
  clipper::Xmap<float> x13map( spgr, cell, grid );
  clipper::Xmap<float> x14map( spgr, cell, grid );
  clipper::Xmap<float> x15map( spgr, cell, grid );
  clipper::Xmap<float> x16map( spgr, cell, grid );
  clipper::Xmap<float> x22map( spgr, cell, grid );
  clipper::Xmap<float> x23map( spgr, cell, grid );
  clipper::Xmap<float> x24map( spgr, cell, grid );
  clipper::Xmap<float> x25map( spgr, cell, grid );
  clipper::Xmap<float> x26map( spgr, cell, grid );
  clipper::Xmap<float> x33map( spgr, cell, grid );
  clipper::Xmap<float> x34map( spgr, cell, grid );
  clipper::Xmap<float> x35map( spgr, cell, grid );
  clipper::Xmap<float> x36map( spgr, cell, grid );
  clipper::Xmap<float> x44map( spgr, cell, grid );
  clipper::Xmap<float> x45map( spgr, cell, grid );
  clipper::Xmap<float> x46map( spgr, cell, grid );
  clipper::Xmap<float> x55map( spgr, cell, grid );
  clipper::Xmap<float> x56map( spgr, cell, grid );
  clipper::Xmap<float> x66map( spgr, cell, grid );
  
  // calculate the term XTY and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y2map[ix]  = x2map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y3map[ix]  = x3map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y4map[ix]  = x4map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y5map[ix]  = x5map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y6map[ix]  = x6map[ix]*ymap[ix]*mmap[ix];

  // calculate the term XTX and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x12map[ix] = x1map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x13map[ix] = x1map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x14map[ix] = x1map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x15map[ix] = x1map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x16map[ix] = x1map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x22map[ix] = x2map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x23map[ix] = x2map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x24map[ix] = x2map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x25map[ix] = x2map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x26map[ix] = x2map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x33map[ix] = x3map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x34map[ix] = x3map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x35map[ix] = x3map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x36map[ix] = x3map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x44map[ix] = x4map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x45map[ix] = x4map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x46map[ix] = x4map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x55map[ix] = x5map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x56map[ix] = x5map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x66map[ix] = x6map[ix]*x6map[ix]*mmap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // apply the convolution filter g to each map
  fltr( y1map, y1map );
  fltr( y2map, y2map );
  fltr( y3map, y3map );
  fltr( y4map, y4map );
  fltr( y5map, y5map );
  fltr( y6map, y6map );
  fltr( x11map, x11map );
  fltr( x12map, x12map );
  fltr( x13map, x13map );
  fltr( x14map, x14map );
  fltr( x15map, x15map );
  fltr( x16map, x16map );
  fltr( x22map, x22map );
  fltr( x23map, x23map );
  fltr( x24map, x24map );
  fltr( x25map, x25map );
  fltr( x26map, x26map );
  fltr( x33map, x33map );
  fltr( x34map, x34map );
  fltr( x35map, x35map );
  fltr( x36map, x36map );
  fltr( x44map, x44map );
  fltr( x45map, x45map );
  fltr( x46map, x46map );
  fltr( x55map, x55map );
  fltr( x56map, x56map );
  fltr( x66map, x66map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) {

    std::vector<float> v(6);
    v[0] = y1map[ix];
    v[1] = y2map[ix];
    v[2] = y3map[ix];
    v[3] = y4map[ix];
    v[4] = y5map[ix];
    v[5] = y6map[ix];

    clipper::Matrix<float> m(6,6);
    m(0,0) = x11map[ix];
    m(0,1) = m(1,0) = x12map[ix];
    m(0,2) = m(2,0) = x13map[ix];
    m(0,3) = m(3,0) = x14map[ix];
    m(0,4) = m(4,0) = x15map[ix];
    m(0,5) = m(5,0) = x16map[ix];
    m(1,1) = x22map[ix];
    m(1,2) = m(2,1) = x23map[ix];
    m(1,3) = m(3,1) = x24map[ix];
    m(1,4) = m(4,1) = x25map[ix];
    m(1,5) = m(5,1) = x26map[ix];
    m(2,2) = x33map[ix];
    m(2,3) = m(3,2) = x34map[ix];
    m(2,4) = m(4,2) = x35map[ix];
    m(2,5) = m(5,2) = x36map[ix];
    m(3,3) = x44map[ix];
    m(3,4) = m(4,3) = x45map[ix];
    m(3,5) = m(5,3) = x46map[ix];
    m(4,4) = x55map[ix];
    m(4,5) = m(5,4) = x56map[ix];
    m(5,5) = x66map[ix];

    v = m.solve(v);

    x1map[ix] = v[0];
    x2map[ix] = v[1];
    x3map[ix] = v[2];
    x4map[ix] = v[3];
    x5map[ix] = v[4];
    x6map[ix] = v[5];
  }

  return true;
}



bool Shift_field_refine::shift_field_u_aniso_const(
  const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
	clipper::Xmap<float>& x1map, clipper::Xmap<float>& x2map, clipper::Xmap<float>& x3map,
  clipper::Xmap<float>& x4map, clipper::Xmap<float>& x5map, clipper::Xmap<float>& x6map,
  float rad, clipper::Resolution& reso, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;
  const clipper::Xmap<float>&   mmap = mask;

  // calculate map coefficients
  clipper::HKL_info hkls (spgr, cell, reso, true);
  clipper::HKL_data<clipper::data32::F_phi> cphi( hkls );
  cmap.fft_to(cphi);

  // calculate the gradient map coefficients (by forming a difference quotient in each anisotropic direction)
  clipper::HKL_data<clipper::data32::F_phi> uphi0( hkls ), uphi1( hkls ), uphi2( hkls ), uphi3( hkls ), uphi4( hkls ), uphi5( hkls );

  clipper::U_aniso_orth u0(0.001,0.00,0.00,0.00,0.00,0.00);
  clipper::U_aniso_orth u1(0.00,0.001,0.00,0.00,0.00,0.00);
  clipper::U_aniso_orth u2(0.00,0.00,0.001,0.00,0.00,0.00);
  clipper::U_aniso_orth u3(0.00,0.00,0.00,0.001,0.00,0.00);
  clipper::U_aniso_orth u4(0.00,0.00,0.00,0.00,0.001,0.00);
  clipper::U_aniso_orth u5(0.00,0.00,0.00,0.00,0.00,0.001);

  uphi0.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u0) );
  uphi1.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u1) );
  uphi2.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u2) );
  uphi3.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u3) );
  uphi4.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u4) );
  uphi5.compute( cphi, clipper::data32::Compute_scale_u_aniso_fphi(1.0,u5) );

  uphi0 = uphi0 - cphi;
  uphi1 = uphi1 - cphi;
  uphi2 = uphi2 - cphi;
  uphi3 = uphi3 - cphi;
  uphi4 = uphi4 - cphi;
  uphi5 = uphi5 - cphi;

  uphi0.compute( uphi0, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi1.compute( uphi1, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi2.compute( uphi2, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi3.compute( uphi3, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi4.compute( uphi4, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );
  uphi5.compute( uphi5, clipper::data32::Compute_scale_u_iso_fphi(1000.0,0.0) );

  // calculate gradient maps
  // x0map = 1 so we don't need to calculate it
  x1map.fft_from( uphi0 );
  x2map.fft_from( uphi1 );
  x3map.fft_from( uphi2 );
  x4map.fft_from( uphi3 );
  x5map.fft_from( uphi4 );
  x6map.fft_from( uphi5 );

  /*
  // output intermediate maps
  clipper::CCP4MAPfile mapout;
  mapout.open_write( "grad1.map" );
  mapout.export_xmap( x1map );
  mapout.close_write();
  mapout.open_write( "grad2.map" );
  mapout.export_xmap( x2map );
  mapout.close_write();
  mapout.open_write( "grad3.map" );
  mapout.export_xmap( x3map );
  mapout.close_write();
  mapout.open_write( "grad4.map" );
  mapout.export_xmap( x4map );
  mapout.close_write();
  mapout.open_write( "grad5.map" );
  mapout.export_xmap( x5map );
  mapout.close_write();
  mapout.open_write( "grad6.map" );
  mapout.export_xmap( x6map );
  mapout.close_write();
  */

 // make xmap - this allocates the maps for future use
  clipper::Xmap<float> y0map( spgr, cell, grid );
  clipper::Xmap<float> y1map( spgr, cell, grid );
  clipper::Xmap<float> y2map( spgr, cell, grid );
  clipper::Xmap<float> y3map( spgr, cell, grid );
  clipper::Xmap<float> y4map( spgr, cell, grid );
  clipper::Xmap<float> y5map( spgr, cell, grid );
  clipper::Xmap<float> y6map( spgr, cell, grid );
  clipper::Xmap<float> x00map( spgr, cell, grid );
  clipper::Xmap<float> x01map( spgr, cell, grid );
  clipper::Xmap<float> x02map( spgr, cell, grid );
  clipper::Xmap<float> x03map( spgr, cell, grid );
  clipper::Xmap<float> x04map( spgr, cell, grid );
  clipper::Xmap<float> x05map( spgr, cell, grid );
  clipper::Xmap<float> x06map( spgr, cell, grid );
  clipper::Xmap<float> x11map( spgr, cell, grid );
  clipper::Xmap<float> x12map( spgr, cell, grid );
  clipper::Xmap<float> x13map( spgr, cell, grid );
  clipper::Xmap<float> x14map( spgr, cell, grid );
  clipper::Xmap<float> x15map( spgr, cell, grid );
  clipper::Xmap<float> x16map( spgr, cell, grid );
  clipper::Xmap<float> x22map( spgr, cell, grid );
  clipper::Xmap<float> x23map( spgr, cell, grid );
  clipper::Xmap<float> x24map( spgr, cell, grid );
  clipper::Xmap<float> x25map( spgr, cell, grid );
  clipper::Xmap<float> x26map( spgr, cell, grid );
  clipper::Xmap<float> x33map( spgr, cell, grid );
  clipper::Xmap<float> x34map( spgr, cell, grid );
  clipper::Xmap<float> x35map( spgr, cell, grid );
  clipper::Xmap<float> x36map( spgr, cell, grid );
  clipper::Xmap<float> x44map( spgr, cell, grid );
  clipper::Xmap<float> x45map( spgr, cell, grid );
  clipper::Xmap<float> x46map( spgr, cell, grid );
  clipper::Xmap<float> x55map( spgr, cell, grid );
  clipper::Xmap<float> x56map( spgr, cell, grid );
  clipper::Xmap<float> x66map( spgr, cell, grid );
  
  // calculate the term XTY and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y0map[ix]  = ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y2map[ix]  = x2map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y3map[ix]  = x3map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y4map[ix]  = x4map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y5map[ix]  = x5map[ix]*ymap[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y6map[ix]  = x6map[ix]*ymap[ix]*mmap[ix];

  // calculate the term XTX and apply the mask
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x00map[ix] = mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x01map[ix] = x1map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x02map[ix] = x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x03map[ix] = x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x04map[ix] = x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x05map[ix] = x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x06map[ix] = x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x12map[ix] = x1map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x13map[ix] = x1map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x14map[ix] = x1map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x15map[ix] = x1map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x16map[ix] = x1map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x22map[ix] = x2map[ix]*x2map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x23map[ix] = x2map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x24map[ix] = x2map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x25map[ix] = x2map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x26map[ix] = x2map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x33map[ix] = x3map[ix]*x3map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x34map[ix] = x3map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x35map[ix] = x3map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x36map[ix] = x3map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x44map[ix] = x4map[ix]*x4map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x45map[ix] = x4map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x46map[ix] = x4map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x55map[ix] = x5map[ix]*x5map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x56map[ix] = x5map[ix]*x6map[ix]*mmap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x66map[ix] = x6map[ix]*x6map[ix]*mmap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // apply the convolution filter g to each map
  fltr( y0map, y0map );
  fltr( y1map, y1map );
  fltr( y2map, y2map );
  fltr( y3map, y3map );
  fltr( y4map, y4map );
  fltr( y5map, y5map );
  fltr( y6map, y6map );
  fltr( x00map, x00map );
  fltr( x01map, x01map );
  fltr( x02map, x02map );
  fltr( x03map, x03map );
  fltr( x04map, x04map );
  fltr( x05map, x05map );
  fltr( x06map, x06map );
  fltr( x11map, x11map );
  fltr( x12map, x12map );
  fltr( x13map, x13map );
  fltr( x14map, x14map );
  fltr( x15map, x15map );
  fltr( x16map, x16map );
  fltr( x22map, x22map );
  fltr( x23map, x23map );
  fltr( x24map, x24map );
  fltr( x25map, x25map );
  fltr( x26map, x26map );
  fltr( x33map, x33map );
  fltr( x34map, x34map );
  fltr( x35map, x35map );
  fltr( x36map, x36map );
  fltr( x44map, x44map );
  fltr( x45map, x45map );
  fltr( x46map, x46map );
  fltr( x55map, x55map );
  fltr( x56map, x56map );
  fltr( x66map, x66map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) {

    std::vector<float> v(7);
    v[0] = y0map[ix];
    v[1] = y1map[ix];
    v[2] = y2map[ix];
    v[3] = y3map[ix];
    v[4] = y4map[ix];
    v[5] = y5map[ix];
    v[6] = y6map[ix];

    clipper::Matrix<float> m(7,7);
    m(0,0) = x00map[ix];
    m(0,1) = m(1,0) = x01map[ix];
    m(0,2) = m(2,0) = x02map[ix];
    m(0,3) = m(3,0) = x03map[ix];
    m(0,4) = m(4,0) = x04map[ix];
    m(0,5) = m(5,0) = x05map[ix];
    m(0,6) = m(0,6) = x06map[ix];
    m(1,1) = x11map[ix];
    m(1,2) = m(2,1) = x12map[ix];
    m(1,3) = m(3,1) = x13map[ix];
    m(1,4) = m(4,1) = x14map[ix];
    m(1,5) = m(5,1) = x15map[ix];
    m(1,6) = m(6,1) = x16map[ix];
    m(2,2) = x22map[ix];
    m(2,3) = m(3,2) = x23map[ix];
    m(2,4) = m(4,2) = x24map[ix];
    m(2,5) = m(5,2) = x25map[ix];
    m(2,6) = m(6,2) = x26map[ix];
    m(3,3) = x33map[ix];
    m(3,4) = m(4,3) = x34map[ix];
    m(3,5) = m(5,3) = x35map[ix];
    m(3,6) = m(6,3) = x36map[ix];
    m(4,4) = x44map[ix];
    m(4,5) = m(5,4) = x45map[ix];
    m(4,6) = m(6,4) = x46map[ix];
    m(5,5) = x55map[ix];
    m(5,6) = m(6,5) = x56map[ix];
    m(6,6) = x66map[ix];

    v = m.solve(v);

    x1map[ix] = v[1];
    x2map[ix] = v[2];
    x3map[ix] = v[3];
    x4map[ix] = v[4];
    x5map[ix] = v[5];
    x6map[ix] = v[6];
  }

  return true;
}

