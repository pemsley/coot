/*! Nautilus target main chain database */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#include "nautilus-target.h"
#include "nautilus-tools.h"
#include "nautilus-ss-find.h"

#include <algorithm>


void NucleicAcidTarget::init( const float c_hi[][3], const float c_lo[][3], const float c_repr[3][3], const int ncoord )
{
  target_.resize( ncoord );
  standard_.resize( 3 );
  for ( int i = 0; i < ncoord; i++ ) {
    target_[i].first  = clipper::Coord_orth(c_hi[i][0],c_hi[i][1],c_hi[i][2]);
    target_[i].second = clipper::Coord_orth(c_lo[i][0],c_lo[i][1],c_lo[i][2]);
  }
  for ( int i = 0; i < 3; i++ )
    standard_[i] = clipper::Coord_orth(c_repr[i][0],c_repr[i][1],c_repr[i][2]);
}


void NucleicAcidTarget::init_stats( const clipper::Xmap<float>& xmap, const int nmax )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  double s = 1.0/double(nmax);
  smin.resize(nmax*nmax*nmax);
  ssum.resize(nmax*nmax*nmax);
  for ( int ix = 0; ix < nmax; ix++ )
    for ( int iy = 0; iy < nmax; iy++ )
      for ( int iz = 0; iz < nmax; iz++ ) {
	const clipper::Coord_frac cf(s*double(ix),s*double(iy),s*double(iz));
	const clipper::Coord_orth co = cf.coord_orth( xmap.cell() );
	const clipper::RTop_orth rtop( clipper::Mat33<>::identity(), co );
	smin[iz+nmax*(iy+nmax*ix)] = score_min( xmap, rtop );
	ssum[iz+nmax*(iy+nmax*ix)] = score_sum( xmap, rtop );
      }
  std::sort( smin.begin(), smin.end() );
  std::sort( ssum.begin(), ssum.end() );
}


double NucleicAcidTarget::radius() const
{
  std::vector<clipper::Coord_orth> coords;
  for ( int i = 0; i < target_.size(); i++ ) {
    coords.push_back( target_[i].first  );
    coords.push_back( target_[i].second );
  }
  double r2 = 0.0;
  for ( int a = 0; a < coords.size(); a++ ) {
    double d2 = coords[a].lengthsq();
    if ( d2 > r2 ) r2 = d2;
  }
  return sqrt( r2 ) + 1.0;
}


float NucleicAcidTarget::score_min( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const
{
  typedef clipper::Interp_cubic I;
  float mn(0.0), mx(0.0);
  for ( int i = 0; i < target_.size(); i++ ) {
    mx = std::min(mx,xmap.interp<I>(xmap.coord_map(rtop*target_[i].first )));
    mn = std::min(mn,xmap.interp<I>(xmap.coord_map(rtop*target_[i].second)));
  }
  return ( mx - mn );
}


float NucleicAcidTarget::score_sum( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const
{
  typedef clipper::Interp_cubic I;
  float mn(0.0), mx(0.0);
  for ( int i = 0; i < target_.size(); i++ ) {
    mx += xmap.interp<I>(xmap.coord_map(rtop*target_[i].first ));
    mn += xmap.interp<I>(xmap.coord_map(rtop*target_[i].second));
  }
  return ( mx - mn );
}


float NucleicAcidTarget::cutoff_min( double p ) const
{
  return smin[ int( ( 1.0 - p ) * double( smin.size() ) ) ];
}


float NucleicAcidTarget::cutoff_sum( double p ) const
{
  return ssum[ int( ( 1.0 - p ) * double( ssum.size() ) ) ];
}


NucleicAcidTargets::NucleicAcidTargets()
{
  const float s_hi[][3] = {
    {  2.227, -0.009, -0.765},
    {  2.235,  0.008,  0.754},
    {  0.874,  0.626,  1.088},
    {  0.009,  0.008, -0.011},
    {  0.806, -0.063, -1.139},
    {  2.864, -0.873, -1.591},
    {  3.015,  0.682,  1.121},
    { -0.692, -1.342,  0.468},
    { -1.647, -2.580,  0.645},
  };
  const float s_lo[][3] = {
    {  4.750, -1.251,  0.004},
    {  3.329,  1.100, -1.718},
    { -0.349, -3.939, -1.101},
    {  1.256,  1.985, -1.337},
    { -0.564,  1.006,  2.432},
    { -1.256,  1.492, -0.346},
    {  2.104, -2.081,  1.712},
    {  1.708, -3.639,  0.375},
    {  0.171, -0.858,  2.479},
  };
  const float s_repr[][3] = {
    {  2.223,  0.000,  0.764}, //  C3'
    {  0.000,  0.000,  0.000}, //  C1'
    {  2.214,  0.000, -0.761}, //  C4'
  };
  const float p_hi[][3] = {
    { -0.017, -0.003,  0.006},
    {  0.742,  0.004,  1.279},
    {  0.745, -0.008, -1.287},
    { -1.250, -1.021,  0.189},
    { -0.974,  1.225, -0.128},
    { -2.115,  1.543, -0.851},
    { -2.186, -2.913,  0.169},
    { -3.088,  2.435, -1.325},
  };
  const float p_lo[][3] = {
    {  1.300, -2.336,  0.987},
    { -0.488, -2.412, -1.727},
    { -3.883, -0.391,  1.565},
    { -2.280,  3.656,  0.753},
    {  0.332,  3.223,  0.764},
    {  2.734,  1.014, -0.556},
    { -1.997, -0.173, -2.950},
    { -0.815,  0.128,  3.003},
  };
  const float p_repr[][3] = {
    { -0.952,  1.277, -0.017}, //  O3'
    { -0.000, -0.000,  0.000}, //  P  
    { -1.027, -1.214,  0.023}, //  O5'
  };
  target_s.init( s_hi, s_lo, s_repr, sizeof(s_hi)/sizeof(s_hi[0]) );
  target_p.init( p_hi, p_lo, p_repr, sizeof(p_hi)/sizeof(p_hi[0]) );
  rad = -1.0;
}


void NucleicAcidTargets::add_pdb( const clipper::String& file )
{
  nadb.add_pdb( file );

  // initialise representative coordinate
  clipper::Coord_orth cc2(0.0,0.0,0.0), co3(0.0,0.0,0.0),
    co5(0.0,0.0,0.0), cc5(0.0,0.0,0.0), cp (0.0,0.0,0.0);
  NucleicAcidDB::Chain nas;
  for ( int r = 0; r < nadb.size(); r++ ) {
    NucleicAcidDB::NucleicAcid na = nadb[r];
    if ( !na.coord_c2().is_null() && !na.coord_o3().is_null() &&
	 !na.coord_o5().is_null() && !na.coord_c5().is_null() &&
	 !na.coord_p().is_null() ) {
      std::vector<clipper::Coord_orth>vf(3);
      vf[0] = na.coord_c3();
      vf[1] = na.coord_c1();
      vf[2] = na.coord_c4();
      const clipper::RTop_orth rtop( vf, target_s.standard() );
      na.transform( rtop );
      nas.add_monomer( na );
    }
  }
  for ( int r = 0; r < nas.size(); r++ ) {
    cc2 += nas[r].coord_c2();
    co3 += nas[r].coord_o3();
    co5 += nas[r].coord_o5();
    cc5 += nas[r].coord_c5();
    cp  += nas[r].coord_p();
  }
  double s = 1.0/double(nas.size());
  cc2 = s * cc2;
  co3 = s * co3;
  co5 = s * co5;
  cc5 = s * cc5;
  cp  = s * cp;
  double d2min = 1.0e20;
  for ( int r = 0; r < nas.size(); r++ ) {
    double d2 = ( ( nas[r].coord_c2() - cc2 ).lengthsq() +
		  ( nas[r].coord_o3() - co3 ).lengthsq() +
		  ( nas[r].coord_o5() - co5 ).lengthsq() +
		  ( nas[r].coord_c5() - cc5 ).lengthsq() +
		  ( nas[r].coord_p()  - cp  ).lengthsq() );
    if ( d2 < d2min ) {
      d2min = d2;
      narepr = nas[r];
    }
  }
}


void NucleicAcidTargets::init_stats( const clipper::Xmap<float>& xmap, const int nmax )
{
  // initialise map stats
  target_s.init_stats( xmap, nmax );
  target_p.init_stats( xmap, nmax );
}


void NucleicAcidTargets::set_search_sphere( const clipper::Coord_orth cent, double radi )
{
  cen = cent;
  rad = radi;
}


void NucleicAcidTargets::superpose_sugar( NucleicAcidDB::Chain& frag, int posn, const NucleicAcidDB::NucleicAcid& na )
{
  std::vector<clipper::Coord_orth> v1(3), v2(3);
  v1[0] = na.coord_c3();
  v1[1] = na.coord_c1();
  v1[2] = na.coord_c4();  
  v2[0] = frag[posn].coord_c3();
  v2[1] = frag[posn].coord_c1();
  v2[2] = frag[posn].coord_c4();
  clipper::RTop_orth rtop( v2, v1 );
  frag.transform( rtop );
}


float NucleicAcidTargets::score_sugar( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const
{
  std::vector<clipper::Coord_orth>vf(3);
  vf[0] = na.coord_c3();
  vf[1] = na.coord_c1();
  vf[2] = na.coord_c4();
  const clipper::RTop_orth rtop( target_s.standard(), vf );
  return target_s.score_sum( xmap, rtop );
}


float NucleicAcidTargets::score_phosphate( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na1, const NucleicAcidDB::NucleicAcid& na2 ) const
{
  std::vector<clipper::Coord_orth> vf(3);
  vf[0] = na1.coord_o3();
  vf[1] = na2.coord_p();
  vf[2] = na2.coord_o5();
  const clipper::RTop_orth rtop(  target_p.standard(), vf );
  return target_p.score_sum( xmap, rtop );
}


NucleicAcidDB::NucleicAcid NucleicAcidTargets::next_na_group( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const
{
  NucleicAcidDB::NucleicAcid nmax;
  float smax = -1.0e20;
  for ( int p = 0; p < nadb.size()-1; p++ ) {
    NucleicAcidDB::Chain frag = nadb.extract( p, 2 );
    if ( frag.is_continuous() ) {
      superpose_sugar( frag, 0, na );
      float score = ( score_phosphate( xmap, frag[0], frag[1] ) +
		      score_sugar    ( xmap, frag[1] ) );
      if ( score > smax ) {
	smax = score;
	nmax = frag[1];
      }
    }
  }
  return nmax;
}


NucleicAcidDB::NucleicAcid NucleicAcidTargets::prev_na_group( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const
{
  NucleicAcidDB::NucleicAcid nmax;
  float smax = -1.0e20;
  for ( int p = 0; p < nadb.size()-1; p++ ) {
    NucleicAcidDB::Chain frag = nadb.extract( p, 2 );
    if ( frag.is_continuous() ) {
      superpose_sugar( frag, 1, na );
      float score = ( score_phosphate( xmap, frag[0], frag[1] ) +
		      score_sugar    ( xmap, frag[0] ) );
      if ( score > smax ) {
	smax = score;
	nmax = frag[0];
      }
    }
  }
  return nmax;
}


const NucleicAcidDB::Chain NucleicAcidTargets::join_sugars( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na1, const NucleicAcidDB::NucleicAcid& na2, int len, double rmsdlim ) const
{
  typedef clipper::Interp_cubic I;

  const int l = len - 1;
  const clipper::Cell&       cell = xmap.cell();

  double rhobest = -1.0e20;
  NucleicAcidDB::Chain chnbest;

  std::vector<clipper::Coord_orth> v1(6), v2(6);
  v1[0] = na1.coord_c4();
  v1[1] = na1.coord_c1();
  v1[2] = na1.coord_c3();
  v1[3] = na2.coord_c4();
  v1[4] = na2.coord_c1();
  v1[5] = na2.coord_c3();
  // check for missing atoms
  for ( int i = 0; i < 6; i++ ) if ( v1[i].is_null() ) return chnbest;
  // search DB
  for ( int p = 0; p < nadb.size()-l; p++ ) {
    NucleicAcidDB::Chain frag = nadb.extract( p, len );
    if ( frag.is_continuous() ) {
      // map fragment onto endpoints
      NucleicAcidDB::NucleicAcid nf1 = frag[0];
      NucleicAcidDB::NucleicAcid nf2 = frag[l];
      v2[0] = nf1.coord_c4();
      v2[1] = nf1.coord_c1();
      v2[2] = nf1.coord_c3();
      v2[3] = nf2.coord_c4();
      v2[4] = nf2.coord_c1();
      v2[5] = nf2.coord_c3();
      // check for sufficiently good fit
      clipper::RTop_orth rtdb( v2, v1 );
      double rmsd = 0.0;
      for ( int i = 0; i < 6; i++ )
	rmsd += ( rtdb*v2[i] - v1[i] ).lengthsq();
      rmsd = sqrt( rmsd / 6.0 );
      if ( rmsd < rmsdlim ) {
	frag.transform( rtdb );
	// now adjust coordinates
	typedef NucleicAcidTools NAT;
	const clipper::Coord_orth cc3 = na1.coord_c3();
	const clipper::Coord_orth cf3 = frag[0].coord_c3();
	const clipper::Coord_orth cc4 = na2.coord_c4();
	const clipper::Coord_orth cf4 = frag[l].coord_c4();
	frag[0] = NucleicAcidDB::NucleicAcid
	  ( na1.coord_p() , na1.coord_o5(), na1.coord_c5(),
	    na1.coord_c4(), na1.coord_o4(), na1.coord_c3(),
	    NAT::coord_adjust( frag[0].coord_o3(), cc3, cf3, cc4, cf4, 4.0 ),
	    na1.coord_c2(), na1.coord_c1(), na1.coord_n() ,
	    std::string(1,na1.type()) );
	for ( int i = 1; i < l; i++ )
	  frag[i] = NucleicAcidDB::NucleicAcid
	    ( NAT::coord_adjust( frag[i].coord_p() , cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_o5(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_c5(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_c4(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_o4(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_c3(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_o3(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_c2(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_c1(), cc3, cf3, cc4, cf4, 4.0 ),
	      NAT::coord_adjust( frag[i].coord_n() , cc3, cf3, cc4, cf4, 4.0 ),
	      "?" );
	frag[l] = NucleicAcidDB::NucleicAcid
	  ( NAT::coord_adjust( frag[l].coord_p() , cc3, cf3, cc4, cf4, 4.0 ),
	    NAT::coord_adjust( frag[l].coord_o5(), cc3, cf3, cc4, cf4, 4.0 ),
	    NAT::coord_adjust( frag[l].coord_c5(), cc3, cf3, cc4, cf4, 4.0 ),
	    na2.coord_c4(),
	    na2.coord_o4(), na2.coord_c3(), na2.coord_o3(),
	    na2.coord_c2(), na2.coord_c1(), na2.coord_n() ,
	    std::string(1,na2.type()) );
	// now score the fragment
	double rho = 0.0;
	for ( int i = 0; i < frag.size(); i++ ) {
	  rho += xmap.interp<I>( frag[i].coord_p().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_o5().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_c5().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_c4().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_o4().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_c3().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_o3().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_c2().coord_frac(cell) );
	  rho += xmap.interp<I>( frag[i].coord_c1().coord_frac(cell) );
	}
	if ( rho > rhobest ) {
	  rhobest = rho;
	  chnbest = frag;
	}
      }
    }
  }
  return chnbest;
}


const clipper::MiniMol NucleicAcidTargets::find( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, int nsugar, int nphosp, double step )
{
  clipper::MiniMol mol_new = mol;

  // make a list of rotations
  std::vector<clipper::RTop_orth> rots;
  // make a list of rotation ops to try
  float glim = 360.0;  // gamma
  float blim = 180.0;  // beta
  float alim = 360.0;  // alpha
  // do a uniformly sampled search of orientation space
  float anglim = clipper::Util::min( alim, glim );
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
	  clipper::Euler_ccp4 euler( alpha, beta, gamma );
	  rots.push_back(clipper::RTop_orth(clipper::Rotation(euler).matrix()));
	}
      }
  }

  // get cutoff (for optimisation)
  clipper::Map_stats stats( xmap );
  double sigcut = stats.mean() + 1.0*stats.std_dev();

  // feature search
  if ( found_s.size() == 0 || found_p.size() == 0 ) {
    SSfind ssfind;
    ssfind.prep_xmap( xmap, std::max( target_sugar().radius(),
				    target_phosphate().radius() ) + 1.0 );
    if ( rad <= 0.0 ) ssfind.prep_search( xmap );
    else              ssfind.prep_search( xmap, sigcut, rad, cen );
    found_s = ssfind.search( target_sugar().target(), rots, sigcut, 0.0 );
    found_p = ssfind.search( target_phosphate().target(), rots, sigcut, 0.0 );

    std::sort( found_s.begin(), found_s.end() );
    std::reverse( found_s.begin(), found_s.end() );
    std::sort( found_p.begin(), found_p.end() );
    std::reverse( found_p.begin(), found_p.end() );
    //std::cout << found_s.size() << "\t" << found_p.size() << std::endl;
    //for ( int i = 0; i < clipper::Util::min( int(found_p.size()), 100 ); i++ )
    //  std::cout << i << ":\tSgr: " << found_s[i].score << "\t" << found_s[i].rot << " " << found_s[i].trn << "\tPho: " << found_p[i].score << "\t" << found_p[i].rot << " " << found_p[i].trn << std::endl;
  }

  const clipper::Grid_sampling& grid = xmap.grid_sampling();
  clipper::MAtomNonBond nb( mol_new, 4.0 );

  // filter lists on translation
  std::vector<SearchResult> filter_s, filter_p;
  for ( int i = 0; i < found_s.size(); i++ ) {
    int it = found_s[i].trn;
    clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
    std::vector<clipper::MAtomIndexSymmetry> atoms = nb( trn, 4.0 );
    if ( atoms.size() == 0 ) filter_s.push_back( found_s[i] );
  }
  for ( int i = 0; i < found_p.size(); i++ ) {
    int it = found_p[i].trn;
    clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
    std::vector<clipper::MAtomIndexSymmetry> atoms = nb( trn, 4.0 );
    if ( atoms.size() == 0 ) filter_p.push_back( found_p[i] );
  }
  //std::cout << "Filter: " << mol.atom_list().size() << std::endl;
  //std::cout << found_s.size() << " " << filter_s.size() << std::endl;
  //std::cout << found_p.size() << " " << filter_p.size() << std::endl;
  //std::cout << mol_new.size() << std::endl;

  // build mono-units on sugars from db fragments
  for ( int i = 0; i < std::min(int(filter_s.size()),nsugar); i++ ) {
    std::vector<clipper::Coord_orth> v1(3), v2(3);
    int ir = filter_s[i].rot;
    int it = filter_s[i].trn;
    clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
    clipper::RTop_orth rtop( rots[ir].rot(), trn );
    NucleicAcidDB::NucleicAcid na = narepr;
    v1[0] = rtop * target_s.standard()[0]; // C3'
    v1[1] = rtop * target_s.standard()[1]; // C1'
    v1[2] = rtop * target_s.standard()[2]; // C4'
    v2[0] = na.coord_c3();
    v2[1] = na.coord_c1();
    v2[2] = na.coord_c4();
    clipper::RTop_orth rtdb( v2, v1 );
    na.transform( rtdb );
    clipper::MPolymer mp;
    na.set_type( '?' );
    mp.insert( na.mmonomer() );
    mol_new.insert( mp );
  }
  //std::cout << mol_new.size() << std::endl;

  // build bi-units on phosphates from db fragments
  for ( int i = 0; i < std::min(int(filter_p.size()),nphosp); i++ ) {
    std::vector<clipper::Coord_orth> v1(3), v2(3);
    int ir = filter_p[i].rot;
    int it = filter_p[i].trn;
    clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
    clipper::RTop_orth rtop( rots[ir].rot(), trn );
    v1[0] = rtop * target_p.standard()[0]; // O3'
    v1[1] = rtop * target_p.standard()[1]; // P
    v1[2] = rtop * target_p.standard()[2]; // O5'
    float smax = -1.0e20;
    clipper::MPolymer mpmax;
    for ( int j = 0; j < nadb.size()-1; j++ ) {
      NucleicAcidDB::Chain frag = nadb.extract( j, 2 );
      if ( frag.is_continuous() ) {
	v2[0] = frag[0].coord_o3();
	v2[1] = frag[1].coord_p();
	v2[2] = frag[1].coord_o5();
	clipper::RTop_orth rtdb( v2, v1 );
	frag.transform( rtdb );
	float score = ( score_sugar( xmap, frag[0] ) +
			score_sugar( xmap, frag[1] ) );
	if ( score > smax ) {
	  clipper::MPolymer mpx;
	  frag[0].set_type( '?' );
	  frag[1].set_type( '?' );
	  mpx.insert( frag[0].mmonomer() );
	  mpx.insert( frag[1].mmonomer() );
	  smax = score;
	  mpmax = mpx;
	}
      }
    }
    mol_new.insert( mpmax );
  }
  //std::cout << mol_new.size() << std::endl;

  return mol_new;
}


const clipper::MiniMol NucleicAcidTargets::grow( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, int ngrow, double fcut ) const
{
  clipper::MiniMol mol_new = mol;
  float scut = target_sugar().cutoff_sum( fcut );
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      NucleicAcidDB::NucleicAcid na;
      for ( int i = 0; i < 25; i++ ) {
	na = next_na_group( xmap, mol_new[c][mol_new[c].size()-1] );
	if ( score_sugar( xmap, na ) < scut ) break;
        na.set_type( '?' );
	mol_new[c].insert( na.mmonomer() );
      }
      for ( int i = 0; i < 25; i++ ) {
	na = prev_na_group( xmap, mol_new[c][0] );
	if ( score_sugar( xmap, na ) < scut ) break;
        na.set_type( '?' );
	mol_new[c].insert( na.mmonomer(), 0 );
      }
    }
  }
  return mol_new;
}


const clipper::MiniMol NucleicAcidTargets::prune( const clipper::MiniMol& mol ) const
{
  // set up clash model (with bases removed)
  clipper::MiniMol mol_nb = mol;
  for ( int c = 0; c < mol_nb.size(); c++ ) {
    for ( int r = 0; r < mol_nb[c].size(); r++ ) {
      NucleicAcidDB::NucleicAcid na( mol_nb[c][r] );
      if ( na.flag() != NucleicAcidDB::NucleicAcid::NONE )
	mol_nb[c][r] = na.mmonomer();
    }
  }
  // find clashes
  clipper::MAtomNonBond nb( mol_nb, 4.0 );
  for ( int c1 = 0; c1 < mol_nb.size(); c1++ ) {
    if ( !mol_nb[c1].exists_property( "NON-NA" ) ) {
      for ( int r1 = 0; r1 < mol_nb[c1].size(); r1++ ) {
	for ( int a1 = 0; a1 < mol_nb[c1][r1].size(); a1++ ) {
	  const clipper::Coord_orth co = mol_nb[c1][r1][a1].coord_orth();
	  std::vector<clipper::MAtomIndexSymmetry> atoms = nb( co, 2.0 );
	  for ( int i = 0; i < atoms.size(); i++ ) {
	    int c2 = atoms[i].polymer();
	    int r2 = atoms[i].monomer();
	    if ( !mol_nb[c2].exists_property( "NON-NA" ) ) {  // clash with NA
	      if ( mol_nb[c1][r1].type() != "~~~" &&
		   mol_nb[c2][r2].type() != "~~~" ) {
		if ( c1 != c2 || r1 < r2-1 || r1 > r2+1 ) {
		  if ( mol_nb[c1].size() < mol_nb[c2].size() )
		    mol_nb[c1][r1].set_type( "~~~" );
		  else
		    mol_nb[c2][r2].set_type( "~~~" );
		}
	      }
	    } else {  // clash with non-NA - delete the NA
	      mol_nb[c1][r1].set_type( "~~~" );
	    }
	  }
	}
      }
    }
  }
  clipper::MiniMol mol_new( mol.spacegroup(), mol.cell() );
  for ( int c = 0; c < mol.size(); c++ ) {
    if ( !mol[c].exists_property( "NON-NA" ) ) {
      clipper::MPolymer mp;
      for ( int r = 0; r < mol[c].size(); r++ ) {
	if ( mol_nb[c][r].type() != "~~~" ) {
	  mp.insert( mol[c][r] );
	} else {
	  if ( mp.size() >= 3 ) mol_new.insert( mp );
	  mp = clipper::MPolymer();
	}
      }
      if ( mp.size() >= 3 ) mol_new.insert( mp );
    } else {
      mol_new.insert( mol[c] );
    }
  }

  return NucleicAcidTools::chain_sort( mol_new );
}


const clipper::MiniMol NucleicAcidTargets::rebuild_chain( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const
{
  const clipper::Spacegroup& spgr = xmap.spacegroup();
  const clipper::Cell&       cell = xmap.cell();
  clipper::MiniMol mol_new = mol;
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      // rebuild chain from fragments
      for ( int r = 0; r < mol_new[c].size()-1; r++ ) {
	NucleicAcidDB::NucleicAcid na1( mol_new[c][r  ] );
	NucleicAcidDB::NucleicAcid na2( mol_new[c][r+1] );
	if ( na1.flag() == NucleicAcidDB::NucleicAcid::COMPLETE &&
	     na2.flag() == NucleicAcidDB::NucleicAcid::COMPLETE ) {
	  clipper::Coord_orth   cref = na1.coord_o3();
	  if ( cref.is_null() ) cref = na1.coord_c3();
	  std::vector<clipper::Coord_orth> cwrk;
	  if ( !na2.coord_p().is_null()  ) cwrk.push_back( na2.coord_p()  );
	  if ( !na2.coord_o5().is_null() ) cwrk.push_back( na2.coord_o5() );
	  if ( !na2.coord_c4().is_null() ) cwrk.push_back( na2.coord_c4() );
	  na2.transform(NucleicAcidTools::symmetry_rtop(cwrk,cref,spgr,cell));
	  NucleicAcidDB::Chain chn = join_sugars( xmap, na1, na2, 2, 0.5 );
	  if ( chn.size() == 2 ) {
	    mol_new[c][r  ] = chn[0].mmonomer();
	    mol_new[c][r+1] = chn[1].mmonomer();
	  }
	}
      }
    }
  }
  return mol_new;
}


const clipper::MiniMol NucleicAcidTargets::rebuild_bases( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const
{
  clipper::MiniMol mol_new = mol;

  // torsions
  const double tors[] = {-3.1, -2.9, -2.8, -2.7, -2.5, -2.2, -2.1, -2.0};
  const int    ntor   = sizeof(tors)/sizeof(tors[0]);

  typedef clipper::Interp_cubic I;

  // build bases in default position
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      for ( int r = 0; r < mol_new[c].size(); r++ ) {
	NucleicAcidDB::NucleicAcid na( mol_new[c][r] );
	if ( na.flag() == NucleicAcidDB::NucleicAcid::COMPLETE ) {
	  clipper::MMonomer mmbest = mol_new[c][r];
	  double scbest = -1.0e6;
	  for ( int itor = 0; itor < ntor; itor++ ) {
	    clipper::Coord_orth cc2( na.coord_o4(), na.coord_c1(), na.coord_n(), 1.37, 2.03, tors[itor] );
	    clipper::Coord_orth cn3( na.coord_c1(), na.coord_n(), cc2, 1.35, 2.09, 3.142 );
	    clipper::Coord_orth cc4( na.coord_n(), cc2, cn3, 1.35, 2.09, 0.0 );
	    clipper::Coord_orth cc5( cc2, cn3, cc4, 1.35, 2.09, 0.0 );
	    clipper::Coord_orth cc6( cn3, cc4, cc5, 1.35, 2.09, 0.0 );
	    clipper::MMonomer mm = na.mmonomer();
	    const int a0 = mm.size();
	    clipper::MAtom ma = clipper::MAtom::null();
	    ma.set_occupancy( 1.0 ); ma.set_u_iso( 0.25 );
	    ma.set_element( "C" );
	    ma.set_coord_orth( cc2 ); ma.set_id( " C2 " ); mm.insert( ma );
	    ma.set_element( "N" );
	    ma.set_coord_orth( cn3 ); ma.set_id( " N3 " ); mm.insert( ma );
	    ma.set_element( "C" );
	    ma.set_coord_orth( cc4 ); ma.set_id( " C4 " ); mm.insert( ma );
	    ma.set_coord_orth( cc5 ); ma.set_id( " C5 " ); mm.insert( ma );
	    ma.set_coord_orth( cc6 ); ma.set_id( " C6 " ); mm.insert( ma );
	    double sc = 0.0;
	    for ( int a = a0; a < mm.size(); a++ )
	      sc += xmap.interp<I>(xmap.coord_map(mm[a].coord_orth()));
	    if ( sc > scbest ) {
	      scbest = sc;
	      mmbest = mm;
	    }
	  }
	  mol_new[c][r] = mmbest;
	}
      }
    }
  }
  // prune any clashing bases
  const char clashatoms[5][5] = {" C2 "," N3 "," C4 "," C5 "," C6 "};
  const int nclashatoms = sizeof(clashatoms)/sizeof(clashatoms[0]);
  clipper::MAtomNonBond nb( mol_new, 4.0 );
  for ( int c1 = 0; c1 < mol_new.size(); c1++ ) {
    if ( !mol_new[c1].exists_property( "NON-NA" ) ) {
      for ( int r1 = 0; r1 < mol_new[c1].size(); r1++ ) {
	bool clash = false;
	for ( int a1 = 0; a1 < mol_new[c1][r1].size(); a1++ ) {
	  bool test = false;
	  for ( int i = 0; i < nclashatoms; i++ )
	    if ( mol_new[c1][r1][a1].id() == clashatoms[i] )
	      test = true;
	  if ( test ) {
	    const clipper::Coord_orth co = mol_new[c1][r1][a1].coord_orth();
	    std::vector<clipper::MAtomIndexSymmetry> atoms = nb( co, 2.0 );
	    for ( int i = 0; i < atoms.size(); i++ ) {
	      int c2 = atoms[i].polymer();
	      int r2 = atoms[i].monomer();
	      if ( c1 != c2 || r1 != r2 ) clash = true;
	    }
	  }
	}
	if ( clash )
	  mol_new[c1][r1] =
	    NucleicAcidDB::NucleicAcid( mol_new[c1][r1] ).mmonomer();
      }
    }
  }

  // now add the remaining atoms
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      // insert OP1, OP2
      clipper::MAtom ma = clipper::MAtom::null();
      ma.set_occupancy( 1.0 ); ma.set_u_iso( 0.25 ); ma.set_element( "O" );
      for ( int r = 1; r < mol_new[c].size(); r++ ) {
	const int i3 = mol_new[c][r-1].lookup(" O3'",clipper::MM::ANY);
	const int ip = mol_new[c][r].lookup(" P  ",clipper::MM::ANY);
	const int i5 = mol_new[c][r].lookup(" O5'",clipper::MM::ANY);
	if ( i3 >= 0 && ip >= 0 && i5 >= 0 ) {
	  const clipper::Coord_orth c3 = mol_new[c][r-1][i3].coord_orth();
	  const clipper::Coord_orth cp = mol_new[c][r][ip].coord_orth();
	  const clipper::Coord_orth c5 = mol_new[c][r][i5].coord_orth();
	  const clipper::Coord_orth u = c3 - cp;
	  const clipper::Coord_orth v = c5 - cp;
	  if ( u.lengthsq() < 3.0 && v.lengthsq() < 3.0 ) {
	    const clipper::Coord_orth x( clipper::Vec3<>::cross(u,v).unit() );
	    const clipper::Coord_orth y( (u.unit()+v.unit()).unit() );
	    const clipper::Coord_orth co1 = cp - 0.80*y + 1.25*x;
	    const clipper::Coord_orth co2 = cp - 0.80*y - 1.25*x;
	    ma.set_coord_orth( co2 ); ma.set_id( " OP2" );
	    mol_new[c][r].insert( ma, i5 );
	    ma.set_coord_orth( co1 ); ma.set_id( " OP1" );
	    mol_new[c][r].insert( ma, i5 );
	  }
	}
      }
      // insert O2'
      for ( int r = 0; r < mol_new[c].size(); r++ ) {
	const int i1 = mol_new[c][r].lookup(" C1'",clipper::MM::ANY);
	const int i2 = mol_new[c][r].lookup(" C2'",clipper::MM::ANY);
	const int i3 = mol_new[c][r].lookup(" C3'",clipper::MM::ANY);
	if ( i1 >= 0 && i2 >= 0 && i3 >= 0 ) {
	  const clipper::Coord_orth c1 = mol_new[c][r][i1].coord_orth();
	  const clipper::Coord_orth c2 = mol_new[c][r][i2].coord_orth();
	  const clipper::Coord_orth c3 = mol_new[c][r][i3].coord_orth();
	  const clipper::Coord_orth u = c1 - c2;
	  const clipper::Coord_orth v = c3 - c2;
	  if ( u.lengthsq() < 3.0 && v.lengthsq() < 3.0 ) {
	    const clipper::Coord_orth x( clipper::Vec3<>::cross(u,v).unit() );
	    const clipper::Coord_orth y( (u.unit()+v.unit()).unit() );
	    const clipper::Coord_orth co = c2 - 0.80*y - 1.20*x;
	    ma.set_coord_orth( co ); ma.set_id( " O2'" );
	    mol_new[c][r].insert( ma, i1 );
	  }
	}
      }
    }
  }

  return mol_new;
}
