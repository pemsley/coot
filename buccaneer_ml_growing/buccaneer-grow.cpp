/*! \file buccaneer-grow.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-grow.h"
#include "ml-grow.h"
#include <clipper/clipper-contrib.h>


const clipper::ftype PI = clipper::Util::pi();

const int Ca_grow::max_conf1 = 50;
const int Ca_grow::max_conf2 = 30;
int Ca_grow::ncpu = 0;


Ca_grow::Ca_grow( int n_grow )
{
  ngrow = n_grow;
}


bool Ca_grow::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const bool& rama_grow ) const
{
  const clipper::MiniMol mold = mol;

  // Find starting chains to expand
  std::vector<Ca_chain> chains = ProteinTools::ca_chains( mold );

  // establish map statistics to determine stopping value for building
  double cutoff = llktarget.llk_distribution( 0.01 );

  // grow the chains
  /*
  clipper::Ramachandran rama1( clipper::Ramachandran::All );
  clipper::Ramachandran rama2( clipper::Ramachandran::NonGly );
  for ( int chn = 0; chn < chains.size(); chn++ )
    grow( chains[chn], xmap, llktarget,        rama1, rama2, cutoff, ngrow );
  */
  Grow_threaded grow( chains, xmap, llktarget, cutoff, ngrow, rama_grow );
  grow( ncpu );
  chains = grow.result();

  // make a new mmdb
  mol = clipper::MiniMol( mold.spacegroup(), mold.cell() );
  ProteinTools::insert_ca_chains( mol, chains );

  // restore the residue types, if any
  ProteinTools::copy_residue_types( mol, mold );
  return true;
}


void Ca_grow::grow( Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& cutoff, const int& ngrow, const bool& rama_grow )
{
  Ca_group ca;
  for ( int i = 0; i < ngrow; i++ ) {
    ca = Ca_grow::next_ca_group( chain, xmap, llktarget, rama1, rama2, rama_grow );
    if ( ca.is_null() ) break;
    if ( llktarget.llk( xmap, ca.rtop_from_std_ori() ) > cutoff ) break;
    chain.push_back( ca );
  }
  for ( int i = 0; i < ngrow; i++ ) {
    ca = Ca_grow::prev_ca_group( chain, xmap, llktarget, rama1, rama2, rama_grow );
    if ( ca.is_null() ) break;
    if ( llktarget.llk( xmap, ca.rtop_from_std_ori() ) > cutoff ) break;
    chain.push_front( ca );
  }
}


Ca_group Ca_grow::next_ca_group( const Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const bool& rama_grow )
{
  Ca_group ca0, ca1, ca2;
  ca0 = chain.back();  // start residue
  if (rama_grow) {
    Rama_ang1 conf1; Rama_ang2 conf2;
    Score_list<Rama_ang1> scores_l1( max_conf1 );
    Score_list<Rama_ang2> scores_l2( max_conf2 );
    double r1, r2;
    const double deg360 = clipper::Util::d2rad(359.0);
    const double deg20  = clipper::Util::d2rad( 20.0);
    const double deg30  = clipper::Util::d2rad( 30.0);
    double phi0 = -9.999; // prev Ramachandran angle
    if ( chain.size() > 1 ) phi0 = chain.ramachandran_phi( chain.size()-1 );
    // search all conformations of first residue
    for ( conf1.r1.psi = 0.0; conf1.r1.psi < deg360; conf1.r1.psi += deg20 )
      for ( conf1.r1.phi = 0.0; conf1.r1.phi < deg360; conf1.r1.phi += deg20 )
        if ( phi0 < -6.283 || rama1.allowed( phi0, conf1.r1.psi ) ) {
          ca1 = ca0.next_ca_group( conf1.r1.psi, PI, conf1.r1.phi );
          r1 = llktarget.llk_approx( xmap, ca1.rtop_from_std_ori() );
          scores_l1.add( r1, conf1 );
        }
    // seach all conformations of second residue using best confirmations of first
    for ( int l1 = 0; l1 < scores_l1.size(); l1++ ) {
      r1 = scores_l1.score(l1);
      conf2.r1 = scores_l1[l1].r1;
      ca1 = ca0.next_ca_group( conf2.r1.psi, PI, conf2.r1.phi );
      for ( conf2.r2.psi = 0.0; conf2.r2.psi < deg360; conf2.r2.psi += deg20 )
        for ( conf2.r2.phi = 0.0; conf2.r2.phi < deg360; conf2.r2.phi += deg30 )
          if ( rama2.favored( conf2.r1.phi, conf2.r2.psi ) ) {
            ca2 = ca1.next_ca_group( conf2.r2.psi, PI, conf2.r2.phi );
            r2 = llktarget.llk_approx( xmap, ca2.rtop_from_std_ori() );
            scores_l2.add( r1+r2, conf2 );
          }
    }
    //return ca0.next_ca_group( scores_l2[0].r1.psi, scores_l2[0].r1.phi );
    if ( scores_l2.size() == 0 ) return Ca_group::null();
    // now calculate full likelihood scores and pick best
    double ll_best = 1.0e6;
    Rama_ang2 ra_best = scores_l2[0];
    for ( int l2 = 0; l2 < scores_l2.size(); l2++ ) {
      ca1 = ca0.next_ca_group( scores_l2[l2].r1.psi, PI, scores_l2[l2].r1.phi );
      ca2 = ca1.next_ca_group( scores_l2[l2].r2.psi, PI, scores_l2[l2].r2.phi );
      r1 = ( llktarget.llk( xmap, ca1.rtop_from_std_ori() ) +
            llktarget.llk( xmap, ca2.rtop_from_std_ori() ) );
      if ( r1 < ll_best ) {
        ll_best = r1;
        ra_best = scores_l2[l2];
      }
    }
  } else {
    ca1 = ml_grow::next_ca_group( ca0, xmap );
    ca2 = ml_grow::next_ca_group( ca1, xmap );
  }
  // refine the result
  Target_fn_refine_c_terminal_build tgt( xmap, llktarget, rama1, rama2, 0.1 );
  return tgt.refine( chain, ca1, ca2 );
}


Ca_group Ca_grow::prev_ca_group( const Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const bool& rama_grow  )
{
  Ca_group ca0, ca1, ca2;
  ca0 = chain.front();  // start residue
  if (rama_grow) {
    Rama_ang1 conf1; Rama_ang2 conf2;
    Score_list<Rama_ang1> scores_l1( max_conf1 );
    Score_list<Rama_ang2> scores_l2( max_conf2 );
    double r1, r2;
    const double deg360 = clipper::Util::d2rad(359.0);
    const double deg20  = clipper::Util::d2rad( 20.0);
    const double deg30  = clipper::Util::d2rad( 30.0);
    double psi0 = -9.999; // prev Ramachandran angle
    if ( chain.size() > 1 ) psi0 = chain.ramachandran_psi( 0 );
    // search all conformations of first residue
    for ( conf1.r1.phi = 0.0; conf1.r1.phi < deg360; conf1.r1.phi += deg20 ) 
      for ( conf1.r1.psi = 0.0; conf1.r1.psi < deg360; conf1.r1.psi += deg20 )
        if ( psi0 < -6.283 || rama1.allowed( conf1.r1.phi, psi0 ) ) {
          ca1 = ca0.prev_ca_group( conf1.r1.phi, PI, conf1.r1.psi );
          r1 = llktarget.llk_approx( xmap, ca1.rtop_from_std_ori() );
        scores_l1.add( r1, conf1 );
      }
    // seach all conformations of second residue using best confirmations of first
    for ( int l1 = 0; l1 < scores_l1.size(); l1++ ) {
      r1 = scores_l1.score(l1);
      conf2.r1 = scores_l1[l1].r1;
      ca1 = ca0.prev_ca_group( conf2.r1.phi, PI, conf2.r1.psi );
      for ( conf2.r2.phi = 0.0; conf2.r2.phi < deg360; conf2.r2.phi += deg20 )
        for ( conf2.r2.psi = 0.0; conf2.r2.psi < deg360; conf2.r2.psi += deg30 )
          if ( rama2.favored( conf2.r2.phi, conf2.r1.psi) ) {
            ca2 = ca1.prev_ca_group( conf2.r2.phi, PI, conf2.r2.psi );
            r2 = llktarget.llk_approx( xmap, ca2.rtop_from_std_ori() );
            scores_l2.add( r1+r2, conf2 );
          }
    }
    //return ca0.prev_ca_group( scores_l2[0].r1.phi, scores_l2[0].r1.psi );
    if ( scores_l2.size() == 0 ) return Ca_group::null();
    // now calculate full likelihood scores and pick best
    double ll_best = 1.0e6;
    Rama_ang2 ra_best = scores_l2[0];
    for ( int l2 = 0; l2 < scores_l2.size(); l2++ ) {
      ca1 = ca0.prev_ca_group( scores_l2[l2].r1.phi, PI, scores_l2[l2].r1.psi );
      ca2 = ca1.prev_ca_group( scores_l2[l2].r2.phi, PI, scores_l2[l2].r2.psi );
      r1 = ( llktarget.llk( xmap, ca1.rtop_from_std_ori() ) +
            llktarget.llk( xmap, ca2.rtop_from_std_ori() ) );
      if ( r1 < ll_best ) {
        ll_best = r1;
        ra_best = scores_l2[l2];
      }
    }
  } else {
    ca1 = ml_grow::prev_ca_group( ca0, xmap );
    ca2 = ml_grow::prev_ca_group( ca1, xmap );
  }
  // refine the result
  Target_fn_refine_n_terminal_build tgt( xmap, llktarget, rama1, rama2, 0.1 );
  return tgt.refine( chain, ca1, ca2 );
}


// thread methods

int Grow_threaded::count = 0;

Grow_threaded::Grow_threaded( const std::vector<Ca_chain>& chains, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const double& cutoff, const int& n_grow, const bool& rama_grow ) : chains_(chains), xmap_(&xmap), llktarget_(&llktarget), cutoff_(cutoff), ngrow(n_grow), rama_grow_(rama_grow)
{
  rama1 = clipper::Ramachandran( clipper::Ramachandran::All );
  rama2 = clipper::Ramachandran( clipper::Ramachandran::NonGly );

  // flag which chains were grown
  done = std::vector<bool>( chains_.size(), false );

  // init thread count
  count = 0;
}

void Grow_threaded::grow( const int& chn )
{
  Ca_grow::grow( chains_[chn], *xmap_, *llktarget_,
                 rama1, rama2, cutoff_, ngrow, rama_grow_ );
  done[chn] = true;
}

bool Grow_threaded::operator() ( int nthread )
{
  bool thread = ( nthread > 0 );
  // try running multi-threaded
  if ( thread ) {
    std::vector<Grow_threaded> threads( nthread-1, (*this) );
    run();  for ( int i = 0; i < threads.size(); i++ ) threads[i].run();
    join(); for ( int i = 0; i < threads.size(); i++ ) threads[i].join();
    // check that it finished
    if ( count >= chains_.size() ) {
      for ( int i = 0; i < threads.size(); i++ ) merge( threads[i] );
    } else {
      thread = false;
    }
  }
  // else run in main thread
  if ( !thread ) {
    for ( int chn = 0; chn < chains_.size(); chn++ ) grow( chn );
  }
  return true;
}

void Grow_threaded::merge( const Grow_threaded& other )
{
  for ( int chn = 0; chn < chains_.size(); chn++ )
    if ( other.done[chn] )
      chains_[chn] = other.chains_[chn];
}

void Grow_threaded::Run()
{
  while (1) {
    lock();
    int chn = count++;
    unlock();
    if ( chn >= chains_.size() ) break;
    grow( chn );
  }
}


// target function for refinement of terminals

Target_fn_refine_n_terminal_build::Target_fn_refine_n_terminal_build( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& rot_step )
{
  xmap_ = &xmap;
  llktarget_ = &llktarget;
  rama1_ = &rama1;
  rama2_ = &rama2;
  rot_step_ = rot_step;
}

double Target_fn_refine_n_terminal_build::operator() ( const std::vector<double>& args ) const
{
  const Ca_chain& chain = *chain_;
  const clipper::Ramachandran& rama1 = *rama1_;
  const clipper::Ramachandran& rama2 = *rama2_;
  double psi0 = -9.999; // prev Ramachandran angle
  if ( chain.size() > 1 ) psi0 = chain.ramachandran_psi( 0 );
  const Ca_group ca0 = chain.front();
  const Ca_group ca1 = ca0.prev_ca_group( args[0], omega1_, args[1] );
  const Ca_group ca2 = ca1.prev_ca_group( args[2], omega2_, args[3] );
  double r = ( (*llktarget_).llk( *xmap_, ca1.rtop_from_std_ori() ) +
               (*llktarget_).llk( *xmap_, ca2.rtop_from_std_ori() ) );
  if ( !rama1.allowed( args[0], psi0 ) ) 
    if ( psi0 > -6.283 ) r += 10.0;
  if ( !rama2.favored( args[2], args[1] ) )
    r += 10.0;
  return r;
}

Ca_group Target_fn_refine_n_terminal_build::refine( const Ca_chain& chain, Ca_group& ca1, Ca_group& ca2 )
{
  // store initial chain
  chain_ = &chain;
  // calculate initial params
  const Ca_group ca0 = chain.front();
  const clipper::Coord_orth c0 = ca0.coord_c();
  const clipper::Coord_orth a0 = ca0.coord_ca();
  const clipper::Coord_orth n0 = ca0.coord_n();
  const clipper::Coord_orth c1 = ca1.coord_c();
  const clipper::Coord_orth a1 = ca1.coord_ca();
  const clipper::Coord_orth n1 = ca1.coord_n();
  const clipper::Coord_orth c2 = ca2.coord_c();
  const clipper::Coord_orth a2 = ca2.coord_ca();
  const clipper::Coord_orth n2 = ca2.coord_n();
  const double phi1 = clipper::Coord_orth::torsion( c1, n0, a0, c0 );
  omega1_           = clipper::Coord_orth::torsion( a1, c1, n0, a0 );
  const double psi1 = clipper::Coord_orth::torsion( n1, a1, c1, n0 );
  const double phi2 = clipper::Coord_orth::torsion( c2, n1, a1, c1 );
  omega2_           = clipper::Coord_orth::torsion( a2, c2, n1, a1 );
  const double psi2 = clipper::Coord_orth::torsion( n2, a2, c2, n1 );
  std::vector<double> args = { phi1, psi1, phi2, psi2 };
  std::vector<double> arg_init;
  std::vector<std::vector<double> > args_init;
  args_init.push_back( args );
  for ( int i = 0; i < num_params(); i++ ) {
    arg_init = args;
    arg_init[i] += rot_step_;
    args_init.push_back( arg_init );
  }
  // simple refinement
  double tol = 0.005 * (*this)( args_init[0] );
  Optimiser_simplex os( tol, 50, Optimiser_simplex::GRADIENT );
  std::vector<double> args_refined = os( *this, args_init );
  return ca0.prev_ca_group( args_refined[0], omega1_, args_refined[1] );
}


Target_fn_refine_c_terminal_build::Target_fn_refine_c_terminal_build( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& rot_step )
{
  xmap_ = &xmap;
  llktarget_ = &llktarget;
  rama1_ = &rama1;
  rama2_ = &rama2;
  rot_step_ = rot_step;
}

double Target_fn_refine_c_terminal_build::operator() ( const std::vector<double>& args ) const
{
  const Ca_chain& chain = *chain_;
  const clipper::Ramachandran& rama1 = *rama1_;
  const clipper::Ramachandran& rama2 = *rama2_;
  double phi0 = -9.999; // prev Ramachandran angle
  if ( chain.size() > 1 ) phi0 = chain.ramachandran_phi( chain.size()-1 );
  const Ca_group ca0 = chain.back();
  const Ca_group ca1 = ca0.next_ca_group( args[0], omega1_, args[1] );
  const Ca_group ca2 = ca1.next_ca_group( args[2], omega2_, args[3] );
  double r = ( (*llktarget_).llk( *xmap_, ca1.rtop_from_std_ori() ) +
               (*llktarget_).llk( *xmap_, ca2.rtop_from_std_ori() ) );
  if ( !rama1.allowed( phi0,    args[0] ) ) 
    if ( phi0 > -6.283 ) r += 10.0;
  if ( !rama2.favored( args[1], args[2] ) )
    r += 10.0;
  return r;
}

Ca_group Target_fn_refine_c_terminal_build::refine( const Ca_chain& chain, Ca_group& ca1, Ca_group& ca2 )
{
  // store initial chain
  chain_ = &chain;
  // calculate initial params
  const Ca_group ca0 = chain.back();
  const clipper::Coord_orth n0 = ca0.coord_n();
  const clipper::Coord_orth a0 = ca0.coord_ca();
  const clipper::Coord_orth c0 = ca0.coord_c();
  const clipper::Coord_orth n1 = ca1.coord_n();
  const clipper::Coord_orth a1 = ca1.coord_ca();
  const clipper::Coord_orth c1 = ca1.coord_c();
  const clipper::Coord_orth n2 = ca2.coord_n();
  const clipper::Coord_orth a2 = ca2.coord_ca();
  const clipper::Coord_orth c2 = ca2.coord_c();
  const double psi1 = clipper::Coord_orth::torsion( n0, a0, c0, n1 );
  omega1_           = clipper::Coord_orth::torsion( a0, c0, n1, a1 );
  const double phi1 = clipper::Coord_orth::torsion( c0, n1, a1, c1 );
  const double psi2 = clipper::Coord_orth::torsion( n1, a1, c1, n2 );
  omega2_           = clipper::Coord_orth::torsion( a1, c1, n2, a2 );
  const double phi2 = clipper::Coord_orth::torsion( c1, n2, a2, c2 );
  std::vector<double> args = { psi1, phi1, psi2, phi2 };
  std::vector<double> arg_init;
  std::vector<std::vector<double> > args_init;
  args_init.push_back( args );
  for ( int i = 0; i < num_params(); i++ ) {
    arg_init = args;
    arg_init[i] += rot_step_;
    args_init.push_back( arg_init );
  }
  // simple refinement
  double tol = 0.005 * (*this)( args_init[0] );
  Optimiser_simplex os( tol, 50, Optimiser_simplex::GRADIENT );
  std::vector<double> args_refined = os( *this, args_init );
  return ca0.next_ca_group( args_refined[0], omega1_, args_refined[1] );
}
