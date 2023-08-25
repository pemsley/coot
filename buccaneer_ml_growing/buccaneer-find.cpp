/*! \file buccaneer-find.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-find.h"
#include "buccaneer-known.h"


int Ca_find::ncpu = 0;


// methods for refinement target


Target_fn_refine_llk_map_target::Target_fn_refine_llk_map_target( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const double& rot_step, const double& trn_step )
{
  xmap_ = &xmap;
  llktarget_ = &llktarget;
  rot_step_ = rot_step;
  trn_step_ = trn_step;
}

double Target_fn_refine_llk_map_target::operator() ( const clipper::RTop_orth& rtop ) const 
{
  return (*llktarget_).llk( *xmap_, rtop );
}

double Target_fn_refine_llk_map_target::operator() ( const std::vector<double>& args ) const
{
  return (*this)( rtop_orth( args ) );
}

clipper::RTop_orth Target_fn_refine_llk_map_target::rtop_orth( const std::vector<double>& args ) const
{
  return clipper::RTop_orth( clipper::Euler<clipper::Rotation::EulerXYZs>(args[0],args[1],args[2]).rotation().matrix() * rtop_.rot(), clipper::Coord_orth(args[3],args[4],args[5]) + rtop_.trn() );
}

clipper::RTop_orth Target_fn_refine_llk_map_target::refine( const clipper::RTop_orth& rtop )
{
  // store initial rtop
  rtop_ = rtop;

  // calculate initial params
  std::vector<std::vector<double> > args_init;
  std::vector<double> arg(6,0.0);
  // identity
  clipper::Euler<clipper::Rotation::EulerXYZs> euler( 0.0, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  // rotation steps
  double step = 0.5 * rot_step_;
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( step, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, step, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, 0.0, step );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  // translation steps
  step = 0.5 * trn_step_;
  arg = args_init[0];
  arg[3] = step;
  arg[4] = 0.0;
  arg[5] = 0.0;
  args_init.push_back( arg );
  arg[3] = 0.0;
  arg[4] = step;
  arg[5] = 0.0;
  args_init.push_back( arg );
  arg[3] = 0.0;
  arg[4] = 0.0;
  arg[5] = step;
  args_init.push_back( arg );
  // simple refinement
  double tol = 0.005 * (*this)( args_init[0] );
  Optimiser_simplex os( tol, 50, Optimiser_simplex::GRADIENT );
  return rtop_orth( os( *this, args_init ) );
}


SSfind::Target::Target( SSfind::SSTYPE type, int num_residues )
{
  const int sslen = num_residues;

  // set up targets
  float ta2[][2][3] = { { { 0.00, 0.00, 0.00}, { 2.50, 0.25, 0.00} },
                        { { 0.87, 0.00, 1.23}, { 1.00,-1.75,-0.25} },
                        { { 0.83, 0.00,-1.18}, { 0.25, 1.75, 0.50} } };
  float ta3[][2][3] = { { { 0.00, 0.00, 0.00}, { 3.00, 0.00, 0.00} },
                        { { 0.87, 0.00, 1.23}, {-1.00, 2.00, 0.25} } };
  float ta4[][2][3] = { { { 0.00, 0.00, 0.00}, { 0.75,-2.75, 0.00} },
                        { { 0.87, 0.00, 1.23}, { 1.25,-2.75, 0.00} } };
  float tb2[][2][3] = { { { 0.00, 0.00, 0.00}, {-1.00, 0.00,-1.75} },
                        { { 0.87, 0.00, 1.23}, { 2.00, 1.25, 0.50} },
                        { { 0.83, 0.00,-1.18}, { 1.75,-1.50,-0.25} } };
  float tb3[][2][3] = { { { 0.00, 0.00, 0.00}, { 2.00, 2.00,-0.25} },
                        { { 0.87, 0.00, 1.23}, { 3.00, 0.50,-0.25} } };
  float tb4[][2][3] = { { { 0.00, 0.00, 0.00}, { 3.25, 0.50,-0.25} },
                        { { 0.87, 0.00, 1.23}, { 3.25, 0.25, 0.00} } };
  typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;
  std::vector<Pair_coord> rep_co, all_co;
  double phi0, psi0;
  if ( type == ALPHA2 ) {
    for ( int i = 0; i < sizeof(ta2)/sizeof(ta2[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta2[i][0][0],ta2[i][0][1],ta2[i][0][2]),
        clipper::Coord_orth(ta2[i][1][0],ta2[i][1][1],ta2[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  } else if ( type == ALPHA3 ) {
    for ( int i = 0; i < sizeof(ta3)/sizeof(ta3[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta3[i][0][0],ta3[i][0][1],ta3[i][0][2]),
        clipper::Coord_orth(ta3[i][1][0],ta3[i][1][1],ta3[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  } else if ( type == ALPHA4 ) {
    for ( int i = 0; i < sizeof(ta4)/sizeof(ta4[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta4[i][0][0],ta4[i][0][1],ta4[i][0][2]),
        clipper::Coord_orth(ta4[i][1][0],ta4[i][1][1],ta4[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  } else if ( type == BETA2 ) {
    for ( int i = 0; i < sizeof(tb2)/sizeof(tb2[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(tb2[i][0][0],tb2[i][0][1],tb2[i][0][2]),
        clipper::Coord_orth(tb2[i][1][0],tb2[i][1][1],tb2[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-120.0); psi0 = clipper::Util::d2rad(120.0);
  } else if ( type == BETA3 ) {
    for ( int i = 0; i < sizeof(tb3)/sizeof(tb3[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(tb3[i][0][0],tb3[i][0][1],tb3[i][0][2]),
        clipper::Coord_orth(tb3[i][1][0],tb3[i][1][1],tb3[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-120.0); psi0 = clipper::Util::d2rad(120.0);
  } else if ( type == BETA4 ) {
    for ( int i = 0; i < sizeof(tb4)/sizeof(tb4[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(tb4[i][0][0],tb4[i][0][1],tb4[i][0][2]),
        clipper::Coord_orth(tb4[i][1][0],tb4[i][1][1],tb4[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-120.0); psi0 = clipper::Util::d2rad(120.0);
  } else {
    for ( int i = 0; i < sizeof(ta3)/sizeof(ta3[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta3[i][0][0],ta3[i][0][1],ta3[i][0][2]),
        clipper::Coord_orth(ta3[i][1][0],ta3[i][1][1],ta3[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  }

  // build residue
  clipper::Coord_orth coa( 0.00, 0.00, 0.00 );  //!< std C-a
  clipper::Coord_orth coc( 0.87, 0.00, 1.23 );  //!< std C
  clipper::Coord_orth con( 0.83, 0.00,-1.18 );  //!< std N
  std::vector<clipper::Coord_orth> mm;
  mm.push_back( con );
  mm.push_back( coa );
  mm.push_back( coc );

  // build secondary structure
  const double pi = clipper::Util::pi();
  std::vector<std::vector<clipper::Coord_orth> > mp;
  for ( int i = 0; i < sslen; i++ ) {
    mm[0] = clipper::Coord_orth( mm[0], mm[1], mm[2], 1.32, 1.99, psi0 );
    mm[1] = clipper::Coord_orth( mm[1], mm[2], mm[0], 1.47, 2.15, pi   );
    mm[2] = clipper::Coord_orth( mm[2], mm[0], mm[1], 1.53, 1.92, phi0 );
    mp.push_back( mm );
  }

  // get RTops
  const int ssmid = (sslen-1)/2;
  std::vector<clipper::RTop_orth> ssops( sslen );
  for ( int m = 0; m < mp.size(); m++ )
    ssops[m] = clipper::RTop_orth( mp[ssmid], mp[m] );

  // build whole ss repr coords
  target_cs.clear();
  for ( int i = 0; i < rep_co.size(); i++ )
    for ( int m = 0; m < ssops.size(); m++ )
      target_cs.push_back( Pair_coord( ssops[m] * rep_co[i].first,
                                       ssops[m] * rep_co[i].second ) );

  // build ca coords
  calpha_cs.clear();
  for ( int m = 0; m < ssops.size(); m++ )
    calpha_cs.push_back( clipper::Coord_orth( ssops[m].trn() ) );
}


void SSfind::prep_xmap( const clipper::Xmap<float>& xmap, const double radius )
{
  // make a 1-d array of gridded density values covering ASU+border
  grid = xmap.grid_sampling();
  grrot = xmap.operator_orth_grid().rot();
  clipper::Grid_range gr0 = xmap.grid_asu();
  clipper::Grid_range gr1( xmap.cell(), xmap.grid_sampling(), radius );
  mxgr = clipper::Grid_range( gr0.min()+gr1.min(), gr0.max()+gr1.max() );
  mapbox = std::vector<float>( mxgr.size(), 0.0 );

  // make 1d list of densities
  clipper::Xmap<float>::Map_reference_index ix( xmap );
  for ( int i = 0; i < mapbox.size(); i++ ) {
    ix.set_coord( mxgr.deindex( i ) );
    mapbox[i] = xmap[ix];
  }
}


void SSfind::prep_search( const clipper::Xmap<float>& xmap )
{
  // make list of results
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  srctrn.clear();
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    srctrn.push_back( grid.index( ix.coord() ) );
}


void SSfind::prep_search( const clipper::Xmap<float>& xmap, const double rhocut, const double radcut, const clipper::Coord_orth centre )
{
  // make list of results
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  srctrn.clear();
  double r2cut = ( radcut > 0.0 ) ? radcut*radcut : 1.0e20;
  clipper::Coord_frac cf = centre.coord_frac( xmap.cell() );
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    if ( xmap[ix] > rhocut ) {
      clipper::Coord_frac df = ix.coord().coord_frac( xmap.grid_sampling() );
      df = df.symmetry_copy_near( xmap.spacegroup(), xmap.cell(), cf ) - cf;
      double r2 = df.lengthsq( xmap.cell() );
      if ( r2 < r2cut )
        srctrn.push_back( grid.index( ix.coord() ) );
    }
}


std::vector<SearchResult> SSfind::search( const std::vector<Pair_coord>& target_cs, const std::vector<clipper::RTop_orth>& ops, const double rhocut, const double frccut ) const
{
  // make a list of indexed, intergerized, rotated lists
  std::vector<std::vector<std::pair<int,int> > > index_lists;
  int i0 = mxgr.index( clipper::Coord_grid(0,0,0) );
  for ( int r = 0; r < ops.size(); r++ ) {
    clipper::RTop_orth op = ops[r];
    std::vector<std::pair<int,int> > tmp;
    for ( int i = 0; i < target_cs.size(); i++ ) {
      const clipper::Coord_map c1( grrot*(op*target_cs[i].first  ) );
      const clipper::Coord_map c2( grrot*(op*target_cs[i].second ) );
      tmp.push_back( std::pair<int,int>( mxgr.index(c1.coord_grid()) - i0,
                                         mxgr.index(c2.coord_grid()) - i0 ) );
    }
    index_lists.push_back( tmp );
  }

  // make list of results
  SearchResult rsltnull = { 0.0, -1, -1 };
  std::vector<SearchResult> rslts( srctrn.size(), rsltnull );
  for ( int i = 0; i < rslts.size(); i++ ) rslts[i].trn = srctrn[i];

  // find ss elements
  float bestcut = 0.0;  // optimisation: abandon searches where score < bestcut
  const float bestscl( frccut ); 
  for ( int i = 0; i < rslts.size(); i++ ) {  // loop over map
    float bestscr = rslts[i].score;
    int   bestrot = rslts[i].rot;
    float bestlim = ( bestscr > bestcut ) ? bestscr : bestcut;
    clipper::Coord_grid cg = grid.deindex( rslts[i].trn );  // coord in grid
    const int index0 = mxgr.index( cg );                    // index in list
    if ( mapbox[index0] > rhocut ) {
      for ( int r = 0; r < index_lists.size(); r++ ) {      // loop over rotns
        const std::vector<std::pair<int,int> >& index_list( index_lists[r] );
        float hi = mapbox[index0+index_list[0].first ];
        float lo = mapbox[index0+index_list[0].second];
        int i = 1;
        while ( hi - lo > bestlim ) {                     // loop over points
          hi = std::min( hi, mapbox[index0+index_list[i].first ] );
          lo = std::max( lo, mapbox[index0+index_list[i].second] );
          i++;
          if ( !( i < index_list.size() ) ) break;
        }
        if ( hi - lo > bestlim ) {
          bestlim = bestscr = hi - lo;
          bestrot = r;
        }
      }
    }
    rslts[i].score = bestscr;  // store
    rslts[i].rot   = bestrot;
    bestcut = std::max( bestscl*bestscr, bestcut );  // optimisation
  }

  // eliminate any results which would have been eliminated by the cutoff
  for ( int i = 0; i < rslts.size(); i++ )
    if ( rslts[i].score < bestcut ) rslts[i] = rsltnull;

  return rslts;
}


// find Ca groups


void Ca_find::prep_prior( clipper::Xmap<float>& prior, const clipper::MiniMol& mol, const double radius ) const
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = prior.spacegroup();
  const clipper::Cell&          cell = prior.cell();
  const clipper::Grid_sampling& grid = prior.grid_sampling();

  // fast path for empty molecule
  if ( mol.size() == 0 ) { prior = 1.0; return; }

  // non-bond search
  const double dmin = radius;
  clipper::MAtomNonBond nb( mol, 0.5*dmin );
  std::vector<clipper::MAtomIndexSymmetry> atoms;
  clipper::Coord_frac f1, f2;
  for ( MRI ix = prior.first(); !ix.last(); ix.next() ) {
    double d2min = pow( dmin, 2.0 );
    f1 = ix.coord().coord_frac(grid);
    atoms = nb.atoms_near( f1.coord_orth(cell), dmin );
    for ( int i = 0; i < atoms.size(); i++ ) {
      const clipper::MAtom& atom = mol.atom(atoms[i]);
      f2 = atom.coord_orth().coord_frac(cell);
      f2 = spgr.symop(atoms[i].symmetry()) * f2;
      f2 = f2.lattice_copy_near( f1 );
      double d2 = ( f2 - f1 ).lengthsq( cell );
      if ( d2 < d2min ) d2min = d2;
    }
    prior[ix] = prob_dist( sqrt(d2min) );
  }
}


void Ca_find::filter_prior( clipper::Xmap<float>& prior, const int modelindex ) const
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;

  // for zero index, use all results
  if ( modelindex <= 0 ) return;

  // otherwise downweight 50% of results on the basis of position in the ASU
  const int div = (modelindex+1)/2;
  const int rem = (modelindex+1)%2;

  for ( MRI ix = prior.first(); !ix.last(); ix.next() ) {
    clipper::Coord_grid cg = ix.coord();
    const int u = cg.u()/div;
    const int v = cg.v()/div;
    const int w = cg.w()/div;
    if ( ( u + v + w ) % 2 != rem ) prior[ix] += 6.0;
  }
}


bool Ca_find::operator() ( clipper::MiniMol& mol, const KnownStructure& knownstruc, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const TYPE type, const int modelindex )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = xmap.spacegroup();
  const clipper::Cell&          cell = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();
  const clipper::MiniMol mold = mol;

  // make prior map
  clipper::MiniMol mol_prior = clipper::MiniMol( xmap.spacegroup(), xmap.cell() );
  for ( int p = 0; p < mold.size(); p++ ) mol_prior.insert( mold[p] );
  knownstruc.copy_to( mol_prior, false );
  clipper::Xmap<float> prior( spgr, cell, grid );
  if ( type == LIKELIHOOD ) prep_prior( prior, mol_prior, 9.0 );
  else                      prep_prior( prior, mol_prior, 6.0 );

  // turn the prior into a z-score
  for ( MRI ix = prior.first(); !ix.last(); ix.next() )
    prior[ix] = -log( prior[ix] );

  // filter prior based on multi-model index parameter
  filter_prior( prior, modelindex );

  // do the fffear search:
  const double stepff = 24.0;
  const double stepss = 24.0;
  const double dres = 5.0;

  // NOW DO THE SEARCH (FIRST CYCLE ONLY)
  if ( results.size() == 0 ) {
    if ( type == LIKELIHOOD ) {
      // make a list of rotation ops to try
      ops = llktarget.rtop_list( spgr, stepff );
      // do the search
      results = search_llk( xmap, llktarget );
    } else {
      // make a list of rotation ops to try
      ops = llktarget.rtop_list( clipper::Spacegroup::p1(), stepss );
      // do the search
      results = search_sec( xmap, llktarget );
    }
  }

  // now create a long scores list from the maps of results
  double zwt = 2.0;  // EXPECTED Z-DIFF BETWEEN CORRECT AND RANDOM SCORES
  Score_list<clipper::RTop_orth> score_long( 20*nfind );
  for ( int i = 0; i < results.size(); i++ ) {
    clipper::RTop_orth  rtop = ops[ results[i].rot ];
    clipper::Coord_grid cg = grid.deindex( results[i].trn );
    rtop.trn() = xmap.coord_orth( cg.coord_map() );
    double score = zwt*results[i].score + prior.get_data(cg);
    score_long.add( score, rtop );
  }

  // create a pruned scores list omitting near-clashes
  Score_list<clipper::RTop_orth> score_trim( nfind );
  clipper::Coord_frac cf1, cf2;
  for ( int i = 0; i < score_long.size(); i++ ) {
    bool clash = false;
    cf1 = clipper::Coord_orth(score_long[i].trn()).coord_frac(cell);
    for ( int j = 0; j < score_trim.size(); j++ ) {
      cf2 = clipper::Coord_orth(score_trim[j].trn()).coord_frac(cell);
      cf2 = cf2.symmetry_copy_near( spgr, cell, cf1 );
      if ( (cf2-cf1).lengthsq(cell) < dres*dres ) { clash = true; break; }
    }
    if ( !clash ) score_trim.add( score_long.score(i), score_long[i] );
  }

  // now refine the best matches
  Score_list<clipper::RTop_orth> score_list( score_trim.size() );
  for ( int i = 0; i < score_trim.size(); i++ ) {
    Target_fn_refine_llk_map_target tgt( xmap, llktarget, 0.2, 0.2 );
    clipper::RTop_orth rtop = tgt.refine( score_trim[i] );
    double score = llktarget.llk( xmap, rtop );
    score_list.add( score, rtop );
  }

  // Now we build a model for output
  mol = clipper::MiniMol( xmap.spacegroup(), xmap.cell() );
  for ( int p = 0; p < mold.size(); p++ ) mol.insert( mold[p] );

  clipper::MPolymer chain;
  chain.set_id(" ");
  int ires = 0;
  for ( int i = 0; i < score_list.size(); i++ ) {
    clipper::MMonomer residue;
    residue.set_type("UNK");
    clipper::MAtom atom = clipper::Atom::null();
    atom.set_occupancy(1.0);
    atom.set_u_iso( exp(-10.0*(score_list.score(i)-score_list.score(0))) );

    atom.set_element( "N" );
    atom.set_id( "N" );
    atom.set_coord_orth( score_list[i] * Ca_group::std_coord_n() );
    residue.insert( atom );

    atom.set_element( "C" );
    atom.set_id( "CA" );
    atom.set_coord_orth( score_list[i] * Ca_group::std_coord_ca() );
    residue.insert( atom );

    atom.set_id( "C" );
    atom.set_coord_orth( score_list[i] * Ca_group::std_coord_c() );
    residue.insert( atom );

    residue.set_seqnum( ires += 2 );
    chain.insert( residue );
  }
  mol.insert( chain );

  return true;
}


void Ca_find::search_op( std::vector<SearchResult>& results, clipper::Xmap<float> xmap1, const clipper::Xmap<int>& xlookp1, const clipper::FFFear_fft<float>& srch, const LLK_map_target& llktarget, const std::vector<clipper::RTop_orth>& ops, int op )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Grid_sampling grid = xlookp1.grid_sampling();
  srch( xmap1, llktarget.llk_target(), llktarget.llk_weight(), ops[op] );
  for ( MRI ix = xmap1.first(); !ix.last(); ix.next() ) {
    const int i = xlookp1[ix];
    if ( xmap1[ix] < results[i].score ) {
      results[i].score = xmap1[ix];
      results[i].rot = op;
      results[i].trn = grid.index( ix.coord() );
    }
  }
}


std::vector<SearchResult> Ca_find::search_llk( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const
{
  const clipper::Spacegroup&    spgr = xmap.spacegroup();
  const clipper::Cell&          cell = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();

  // set up maps
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  clipper::Xmap<int> xlookp1( clipper::Spacegroup::p1(), cell, grid );
  int lresult = 0;
  {
    clipper::Xmap<int> xlooksg( spgr, cell, grid );
    for ( MRI ix = xlooksg.first(); !ix.last(); ix.next() )
      xlooksg[ix] = lresult++;
    for ( MRI ix = xlookp1.first(); !ix.last(); ix.next() )
      xlookp1[ix] = xlooksg.get_data(ix.coord());
  }

  // set up target
  clipper::FFFear_fft<float> srch( xmap );
  srch.set_fft_type( clipper::FFFear_fft<float>::Sparse );
  srch.set_resolution( clipper::Resolution( resol_ ) );

  // set up z-scoring
  clipper::Xmap<float> xmap1( clipper::Spacegroup::p1(), cell, grid );
  srch( xmap1, llktarget.llk_target(), llktarget.llk_weight(), ops[0] );
  clipper::Map_stats zstats( xmap1 );

  // do the search
  /*
    SearchResult result = { 1.0e20, 0, 0 };
    std::vector<SearchResult> results( lresult, result );
    for ( int op = 0; op < ops.size(); op++ )
    search_op( results, xmap1, xlookp1, srch, llktarget, ops, op );
  */
  Search_threaded sthr( xlookp1, srch, llktarget, ops, lresult );
  sthr( ncpu );
  std::vector<SearchResult> result = sthr.results();
    
  // convert to Z scores
  for ( int i = 0; i < result.size(); i++ )
    result[i].score = ( result[i].score - zstats.mean() ) / zstats.std_dev();

  return result;
}


std::vector<SearchResult> Ca_find::search_sec( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const
{
  const clipper::Grid_sampling& grid = xmap.grid_sampling();

  // control params
  const int sslen = 3;
  const double rad = 1.9 * sslen + 1.0;

  // get cutoff (for optimisation)
  clipper::Map_stats stats( xmap );
  double sigcut = stats.mean() + 1.0 * stats.std_dev();

  // prepare target map
  SSfind ssfind;
  ssfind.prep_xmap( xmap, rad );
  ssfind.prep_search( xmap );

  /*
  // do initial search with itentity op to get stats
  clipper::RTop_orth rtid( clipper::RTop_orth::identity() );
  std::vector<clipper::RTop_orth> opsid( 1, rtid );
  ssfind.prep_results( xmap, -1.0e20 );
  ssfind.prep_target( SSfind::ALPHA3, sslen );
  ssfind.search( opsid, -1.0e20 );
  ssfind.prep_target( SSfind::BETA3, sslen );
  ssfind.search( opsid, -1.0e20 );
  const std::vector<SearchResult> resultz = ssfind.results();
  double s0(0.0), s1(0.0), s2(0.0);
  for ( int i = 0; i < resultz.size(); i++ ) {
    s0 += 1.0;
    s1 += resultz[i].score;
    s2 += resultz[i].score * resultz[i].score;
  }
  s1 /= s0; s2 /= s0;
  s2 = sqrt( s2 - s1*s1 );
  */

  // do full search
  SSfind::Target targeta( SSfind::ALPHA3, sslen );
  SSfind::Target targetb( SSfind::BETA3 , sslen );
  std::vector<SearchResult> resulta, resultb;
  resulta = ssfind.search( targeta.target_coords(), ops, sigcut, 0.0 );
  resultb = ssfind.search( targetb.target_coords(), ops, sigcut, 0.0 );
  std::vector<SearchResult> result( resulta.size() );
  for ( int i = 0; i < result.size(); i++ ) {
    if ( resulta[i].score > resultb[i].score ) result[i] = resulta[i];
    else                                       result[i] = resultb[i];
  }

  // rescore
  double s0(0.0), s1(0.0), s2(0.0);
  int zstep = std::max( int(result.size()/1000), 1 );
  for ( int i = 0; i < result.size(); i++ ) {
    clipper::Coord_orth co =
      xmap.coord_orth( grid.deindex( result[i].trn ).coord_map() );
    // update score for target
    if ( result[i].rot >= 0 ) {
      clipper::RTop_orth rtop( ops[ result[i].rot ].rot(), co );
      result[i].score = llktarget.llk_approx( xmap, rtop );
    } else {
      result[i].score = 1.0e20;
    }
    // and accumulate z-score stats
    if ( i % zstep == 0 ) {
      clipper::RTop_orth rtid( clipper::Mat33<>::identity(), co );
      double score = llktarget.llk_approx( xmap, rtid );
      s0 += 1.0;
      s1 += score;
      s2 += score * score;
    }
  }
  s1 /= s0; s2 /= s0;
  s2 = sqrt( s2 - s1*s1 );

  // convert to Z scores
  for ( int i = 0; i < result.size(); i++ )
    result[i].score = ( result[i].score - s1 )/ s2;

  /*
  // build the result
  const clipper::Cell& cell = xmap.cell();
  std::sort( result.begin(), result.end() );
  clipper::MiniMol mol( xmap.spacegroup(), xmap.cell() );
  const clipper::String chainids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  std::vector<clipper::Coord_orth> calphas = ssfind.calpha_coords();
  clipper::Coord_orth co, cr;
  for ( int i = 0; i < result.size(); i++ ) {
    if ( i >= chainids.size() ) break;
    std::cout << i << " " << result[i].score << " " << result[i].rot << " " << result[i].trn << std::endl;
    clipper::MPolymer mpca;
    co = grid.deindex(result[i].trn).coord_frac(grid).coord_orth(cell);
    for ( int r = 0; r < calphas.size(); r++ ) {
      cr = ops[result[i].rot].inverse() * calphas[r];
      clipper::MAtom maca = clipper::MAtom::null();
      maca.set_id( " CA " );  maca.set_element( "C" );
      maca.set_u_iso( 0.25 ); maca.set_occupancy( 1.00 ); 
      maca.set_coord_orth( co + cr );
      clipper::MMonomer mmca;
      mmca.set_type( "ALA" );
      mmca.set_seqnum( r+1 );
      mmca.insert( maca );
      mpca.insert( mmca );
    }
    mpca.set_id( chainids.substr( i, 1 ) );
    mol.insert( mpca );
  }
  clipper::MMDBfile mfile;
  mfile.export_minimol( mol );
  mfile.write_file( "ssfind.pdb" );
  */

  // remove results where no match was found
  std::vector<SearchResult> result_trim;
  for ( int i = 0; i < result.size(); i++ )
    if ( result[i].rot >= 0 ) result_trim.push_back( result[i] );

  return result_trim;
}


// Threadable search class methods

Search_threaded::Search_threaded( const clipper::Xmap<int>& xlookp1, const clipper::FFFear_fft<float>& srch, const LLK_map_target& llktarget, const std::vector<clipper::RTop_orth>& ops, const int lresult ) : xlookp1_(&xlookp1), srch_(&srch), llktarget_(&llktarget), ops_(ops)
{
  // set up work map
  const clipper::Cell&          cell = xlookp1.cell();
  const clipper::Grid_sampling& grid = xlookp1.grid_sampling();
  xmap1_.init( clipper::Spacegroup::p1(), cell, grid );
  SearchResult result = { 1.0e20, 0, 0 };
  results_ = std::vector<SearchResult>( lresult, result );

  n1_ = n2_ = 0;
  done = false;
}

void Search_threaded::search( const int& op )
{
  Ca_find::search_op( results_, xmap1_, *xlookp1_, *srch_, *llktarget_, ops_,
                      op );
}

bool Search_threaded::operator() ( int nthread )
{
  bool thread = ( nthread > 0 );
  // try running multi-threaded
  const int n = ops_.size();
  if ( thread ) {
    std::vector<Search_threaded> threads( nthread-1, (*this) );
    set_range( 0, n/nthread );
    for ( int i = 0; i < threads.size(); i++ )
      threads[i].set_range( ((i+1)*n)/nthread, ((i+2)*n)/nthread );
    done = true;
    run();  for ( int i = 0; i < threads.size(); i++ ) threads[i].run();
    join(); for ( int i = 0; i < threads.size(); i++ ) threads[i].join();
    for ( int i = 0; i < threads.size(); i++ ) merge( threads[i] );
    // check that it finished
    if ( !done ) thread = false;
  }
  // else run in main thread
  if ( !thread ) {
    int lresult = results_.size();
    SearchResult result = { 1.0e20, 0, 0 };
    results_ = std::vector<SearchResult>( lresult, result );
    for ( int op = 0; op < n; op++ ) search( op );
  }
  return true;
}

void Search_threaded::merge( const Search_threaded& other )
{
  for ( int i = 0; i < results_.size(); i++ )
    if ( other.results_[i].score < results_[i].score )
      results_[i] = other.results_[i];
  done = done && other.done;
}

void Search_threaded::Run()
{
  for ( int n = n1_; n < n2_; n++ ) search( n );
  done = true;
}
