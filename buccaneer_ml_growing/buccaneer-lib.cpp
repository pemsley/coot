/*! \file buccaneer-lib.cpp buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-lib.h"

#include <clipper/clipper-contrib.h>

#include <algorithm>


/*! The target is constructed and initialised to zero. It can then be
  filled by accumulation or by loading the target and weight.
  \param rad Radius in which LLK function will be used, in Angstroms.
  \param sampling Sampling spacing for the LLK function, in Angstroms. */
void LLK_map_target::init( const clipper::ftype& rad, const clipper::ftype& sampling, TYPE type )
{
  radius = rad;

  clipper::ftype extent = radius + 2.0 * sampling;  // model box size
  
  clipper::ftype ng = rint( extent / sampling );
  clipper::Grid grid( 2*int(ng)+1, 2*int(ng)+1, 2*int(ng)+1 );
  clipper::RTop<> rtop( clipper::Mat33<>( 1.0/sampling, 0.0, 0.0,
					  0.0, 1.0/sampling, 0.0,
					  0.0, 0.0, 1.0/sampling ),
			clipper::Vec3<>( ng, ng, ng ) );

  target.init( grid, rtop );
  weight.init( grid, rtop );

  type_ = type;
  fasttgt.set_type( type );
  slowtgt.set_type( type );

  tgt_scl = 1.0;
  tgt_off = 0.0;

  naccum = 0;  // number of samples
}

/*! This most be called after loading and before saving or using the
  LLK targets in any way. i.e. if you accumulate the log likelihood
  target in one program and use it in another, then you must call
  prep_llk() both before saving the targets in the first program and
  after loading in the second program. */
void LLK_map_target::prep_llk()
{
  clipper::NXmap_base::Map_reference_index ix;

  // first make a llk target, if necessary
  if ( naccum != 0 ) {
    // aliases for maps
    clipper::NXmap<float>& mrho = target;
    clipper::NXmap<float>& mrho2 = weight;
    // calculate density moments for whole map
    clipper::ftype sn = 0.0, sr = 0.0, sr2 = 0.0;
    for ( ix = mrho.first(); !ix.last(); ix.next() )
      if ( mrho2[ix] > 0.0 ) {
	sn  += naccum;        // overall desnity stats
	sr  += mrho[ix];
	sr2 += mrho2[ix];
      }
    float rmap = sr/sn;
    float smap = sqrt( sn*sr2 - sr*sr )/sn;
    // and for individual positions within map
    for ( ix = mrho.first(); !ix.last(); ix.next() )
      if ( mrho2[ix] > 0.0 ) {
	mrho[ix]  /= naccum;  // local density stats
	mrho2[ix] /= naccum;  // limit std to no less than 3% of map std
	mrho2[ix] = sqrt( std::max( mrho2[ix] - mrho[ix]*mrho[ix],
					      0.001f*smap*smap ) );
      }
    // convert density moments to llk target
    float rho, std, v1, v2, w1, w2;
    for ( ix = mrho.first(); !ix.last(); ix.next() )
      if ( mrho2[ix] > 0.0 ) {
	rho = mrho[ix];
	std = mrho2[ix];
	v1 = std*std;
	v2 = smap*smap;
	w1 = std::max( v2/v1-1.0, 0.001 );  // weight must be +ve
	w2 = std::min( 1.0/w1, 2.0 );  // limit density upweighting
	target[ix] = rho + w2*(rho-rmap);
	weight[ix] = w1 * 0.5/v2;
      }
    // done
    naccum = 0;
  }

  // now truncate the functions at the limit radii
  for ( ix = target.first(); !ix.last(); ix.next() )
    if ( ix.coord_orth().lengthsq() > radius * radius )
      weight[ix] = target[ix] = 0.0;

  // now assemble optimised llk lists
  clipper::Coord_grid cg, dg;
  clipper::Coord_orth co;
  clipper::Coord_map cm;
  const clipper::ftype r0 = radius * 3.0 / 8.0;
  clipper::Coord_grid cg0( clipper::Coord_map( target.operator_orth_grid() * clipper::Coord_map( 0.0, 0.0, 0.0 ) ).coord_grid() );
  clipper::Coord_grid cg1( clipper::Coord_map( target.operator_orth_grid() * clipper::Coord_map( radius, 0.0, 0.0 ) ).coord_grid() );
  int irad = cg1.u() - cg0.u();

  // Make list of sites for the approximate LLK fn
  for ( cg.u() = -1; cg.u() <= 1; cg.u()++ )
    for ( cg.v() = -1; cg.v() <= 1; cg.v()++ )
      for ( cg.w() = -1; cg.w() <= 1; cg.w()++ )
	if ( (cg.u()+cg.v()+cg.w())%2 == 0 ) {
	  co = clipper::Coord_orth( r0 * cg.coord_map() );
	  cm = target.coord_map( co );
	  cm = clipper::Coord_map( target.operator_orth_grid() * co );
	  fasttgt.insert( co, target.interp<clipper::Interp_cubic>( cm ),
                              weight.interp<clipper::Interp_cubic>( cm ) );
	}
  // Make list of sites for the accurate LLK fn
  for ( ix = target.first(); !ix.last(); ix.next() ) {
    cg = ix.coord();
    if ( (cg.u()+cg.v()+cg.w())%2 == 0 ) {  // only use alternate grids
      dg = cg - cg0;
      if ( dg*dg <= irad*irad ) {
	co = target.coord_orth( cg.coord_map() );
	slowtgt.insert( co, target[ix], weight[ix] );
      }
    }
  }
}

/*! Accumulate the statistics for a log-likelihood target using a known
  map and operator maping into that map. After accumulation is
  completed, you must call prep_llk().
  \param xmap The known map from which to accumulate density statistics.
  \param rtop The operator from the target map at the origin into the xmap. */
void LLK_map_target::accumulate( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop )
{
  // aliases for maps
  clipper::NXmap<float>& mrho = target;
  clipper::NXmap<float>& mrho2 = weight;

  // zero maps if necessary
  if ( naccum == 0.0 ) target = weight = float(0.0);
  naccum++;

  // accumulate stats
  clipper::ftype extentsq =
    pow( mrho.operator_grid_orth().rot()(0,0) * (mrho.grid().nu()-1)/2, 2 );
  clipper::NXmap_base::Map_reference_index ix;
  float rho;
  for ( ix = mrho.first(); !ix.last(); ix.next() )
    if ( ix.coord_orth().lengthsq() <= extentsq ) {
      rho = xmap.interp<clipper::Interp_cubic>( (rtop*ix.coord_orth()).coord_frac(xmap.cell()) );
      mrho[ix]  += rho;
      mrho2[ix] += rho*rho;
    }
}


/*! Construct a list of rotations for a given spacegroup and sampling.
  \param spgr The spacegroup.
  \param step The search andle step (in degrees, NOT radians). */
std::vector<clipper::RTop_orth> LLK_map_target::rtop_list( const clipper::Spacegroup& spgr, const clipper::ftype& step ) const
{
  std::vector<clipper::RTop_orth>  ops;
  // make a list of rotation ops to try
  float glim = 360.0;  // gamma
  float blim = 180.0;  // beta
  float alim = 360.0;  // alpha
  // reduce search angles by symmetry rotations
  alim /= float( spgr.order_of_symmetry_about_axis( clipper::Spacegroup::C ) );
  if ( spgr.order_of_symmetry_about_axis( clipper::Spacegroup::A ) % 2 == 0 ||
       spgr.order_of_symmetry_about_axis( clipper::Spacegroup::B ) % 2 == 0 )
    blim /= 2.0;
  // do a uniformly sampled search of orientation space
  float anglim = std::min( alim, glim );
  for ( float bdeg=step/2; bdeg < 180.0; bdeg += step ) {
    float beta = clipper::Util::d2rad(bdeg);
    float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
    float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
    for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
      for ( float thmi=smi/2; thmi < 360.0; thmi += smi ) {
	float adeg = clipper::Util::mod(0.5*(thpl+thmi),360.0);
	float gdeg = clipper::Util::mod(0.5*(thpl-thmi),360.0);
	if ( adeg <= alim && bdeg <= blim && gdeg <= glim ) {
	  float alpha = clipper::Util::d2rad(adeg);
	  float gamma = clipper::Util::d2rad(gdeg);
	  const clipper::Euler_ccp4 euler( alpha, beta, gamma );
          const clipper::Rotation   rotn(euler);
          const clipper::RTop_orth  op(rotn.matrix().inverse());
	  ops.push_back(op);  // Note: inverse above for back compatibility
	}
      }
  }
  // return the rotations
  return ops;
}


/*! A log-likelihood FFFear search is performed for the target in the given map.
  \param resultscr The best scores.
  \param resultrot The best rotations.
  \param resulttrn The best translations.
  \param xmap The map to search.
  \param rtops The oprientations to search. */
void LLK_map_target::search( clipper::Xmap<float>& resultscr, clipper::Xmap<int>& resultrot, clipper::Xmap<int>& resulttrn, const clipper::Xmap<float>& xmap, const std::vector<clipper::RTop_orth>& rtops ) const
{
  // set up results
  const clipper::Spacegroup&    spgr = xmap.spacegroup();
  const clipper::Cell&          cell = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();
  resultscr.init( spgr, cell, grid );
  resultrot.init( spgr, cell, grid );
  resulttrn.init( spgr, cell, grid );
  resultscr = 1.0e20;

  // now search for ML target in each orientation in turn
  clipper::Xmap<float> resultp1( clipper::Spacegroup::p1(), cell, grid );
  clipper::Xmap<float>::Map_reference_index i1(resultp1);
  clipper::Xmap<float>::Map_reference_coord ix(resultscr);

  // set up z scoring
  clipper::FFFear_fft<float> srch( xmap );
  clipper::NX_operator nxop( xmap, target, rtops[0] );
  srch( resultp1, target, weight, nxop );
  clipper::Map_stats zstats( resultp1 );

  // loop over orientations
  for ( int op = 0; op < rtops.size(); op++ ) {
    // do the fffear search
    clipper::NX_operator nxop( xmap, target, rtops[op].inverse() );
    srch( resultp1, target, weight, nxop );

    // store best scores
    for ( i1 = resultp1.first(); !i1.last(); i1.next() ) {
      ix.set_coord( i1.coord() );
      float score = ( resultp1[i1] - zstats.mean() ) / zstats.std_dev();
      if ( score < resultscr[ix] ) {
	resultscr[ix] = score;
	resultrot[ix] = op;
	resulttrn[ix] = grid.index( i1.coord() );
      }
    }
  }
}


clipper::String LLK_map_target::format() const
{
  clipper::String result = "";
  clipper::NXmap_base::Map_reference_index ix;
  clipper::ftype y = 0.0;
  for ( ix = target.first(); !ix.last(); ix.next() ) y += pow(target[ix],2);
  y = sqrt( y / double( target.grid().nu() * target.grid().nv() * 
		        target.grid().nw() / 2 ) );
  for ( int w = 0; w < target.grid().nw(); w++ ) {
    result += "\n w=" + clipper::String(w,3);
    for ( int v = 0; v < target.grid().nv(); v++ ) {
      result += "\n  v=" + clipper::String(v,3) + " ";
      for ( int u = 0; u < target.grid().nu(); u++ ) {
	clipper::Coord_grid g( u, v, w );
	float x = target.get_data(g) / y;
	if ( x < 0.0 ) result += "-";
	else if ( x > 1.0 ) result += "#";
	else if ( x > 0.3 ) result += "+";
	else if ( x > 0.1 ) result += ".";
	else result += " ";
      }
    }
  }
  result += "\n";
  return result;
}


void LLK_map_target::prep_llk_distribution( const clipper::Xmap<float>& xmap )
{
  const int ndist = 5000;
  llkdist.resize( ndist );
  // Magic ints, m3 = 2*m1-1, diffs are 13, 8 from Fibonacci, m1^3 ~ 10K 
  const double m1(22.0), m2(35.0), m3(43.0); // useful up to 10K samples
  const double twopi = clipper::Util::twopi();
  for ( int i = 0; i < ndist; i++ ) {
    const double x = double(i)/double(ndist);
    const double x1 = clipper::Util::mod( m1*x, 1.0 );
    const double x2 = clipper::Util::mod( m2*x, 1.0 );
    const double x3 = clipper::Util::mod( m3*x, 1.0 );
    const clipper::Euler_ccp4 rot( twopi*x1, twopi*x2, twopi*x3 );
    const clipper::Rotation   ro( rot );
    const clipper::Coord_frac cf( x1, x2, x3 );
    const clipper::Coord_orth trn( cf.coord_orth( xmap.cell() ) );
    const clipper::RTop_orth rtop( ro.matrix(), trn );
    llkdist[i] = llk( xmap, rtop );
  }
  std::sort( llkdist.begin(), llkdist.end() );
}


clipper::ftype LLK_map_target::llk_distribution( const clipper::ftype& ordinal ) const
{
  if ( llkdist.size() == 0 ) return clipper::Util::nan();
  const int i = clipper::Util::intf( ordinal*double(llkdist.size()) + 0.001 );
  const int j = clipper::Util::bound( 0, i, int(llkdist.size()-1) );
  return llkdist[j];
}


void LLK_map_target::Sampled::insert( clipper::Coord_orth coord, clipper::ftype tgt, clipper::ftype wgt ) {
  repxyz.push_back( coord );
  reptgt.push_back( tgt );
  repwgt.push_back( wgt );
}

/* \return The log likelihood */
clipper::ftype LLK_map_target::Sampled::llk( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const {
  clipper::ftype r( 0.0 ), s( 0.0 );
  for ( int i = 0; i < repxyz.size(); i++ ) {
    r += repwgt[i] * pow( xmap.interp<clipper::Interp_linear>( (rtop*repxyz[i]).coord_frac(xmap.cell()) ) - reptgt[i], 2 );
    s += repwgt[i];
  }
  return r/s;
}

/* \return The negative of the correlation */
clipper::ftype LLK_map_target::Sampled::correl( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const {
  clipper::ftype x, y, w, sw, swx, swy, swxx, swyy, swxy;
  sw = swx = swy = swxx = swyy = swxy = 0.0;
  for ( int i = 0; i < repxyz.size(); i++ ) {
    w = repwgt[i];
    x = reptgt[i];
    y = xmap.interp<clipper::Interp_linear>( (rtop*repxyz[i]).coord_frac(xmap.cell()) );
    sw   += w;
    swx  += w * x;
    swy  += w * y;
    swxx += w * x * x;
    swyy += w * y * y;
    swxy += w * x * y;
  }
  return -( sw*swxy - swx*swy ) / sqrt( std::max(
      ( sw*swxx - swx*swx ) * ( sw*swyy - swy*swy ), 1.0e-20 ) );
}

