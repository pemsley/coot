/*! \file simplex-lib.cpp optimiser library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include "simplex-lib.h"
#include <clipper/clipper.h>


Optimiser_simplex::Optimiser_simplex( double tolerance, int max_cycles, TYPE type )
{
  tolerance_ = tolerance;
  max_cycles_ = max_cycles;
  type_ = type;
  debug_mode = false;
}


std::vector<double> Optimiser_simplex::operator() ( const Target_fn_order_zero&
target_fn, const std::vector<std::vector<double> >& args ) const
{
  enum STEP { UNKN, EXTN, NRML, CTRN, CTRX };

  // check input
  int size = target_fn.num_params();
  if ( args.size() != size+1 )
    clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );
  for ( int i = 0; i < args.size(); i++ )
    if ( args[i].size() != size )
      clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );

  // make arrays
  std::vector<std::vector<double> > params( size + 1 );
  std::vector<double> fn( size + 1 );
  // calc target fn
  for ( int i = 0; i < args.size(); i++ ) {
    params[i] = args[i];
    fn[i] = target_fn( args[i] );
  }

  // simplex loop
  int iw, ib;
  iw = ib = 0;
  double f0(0.0), f1(0.0), f2(0.0);
  for ( int cyc = 0; cyc < max_cycles_; cyc++ ) {
    if ( debug_mode )  // DEBUG OUTPUT
      for ( int i = 0; i < params.size(); i++ ) {
        std::cout << i << "  " << clipper::String( fn[i] ) << "  ";
        for ( int j = 0; j < size; j++ )
          std::cout << clipper::String( params[i][j] ) << "\t";
        std::cout << "\n";
      }
    // find worst point
    iw = ib = 0;
    for ( int i = 0; i < params.size(); i++ ) {
      if ( fn[i] > fn[iw] ) iw = i;
      if ( fn[i] < fn[ib] ) ib = i;
    }
    // termination condition
    if ( fn[iw] - fn[ib] < tolerance_ ) break;
    // find centroid of the rest
    std::vector<double> centr( size, 0.0 ),
      shift( size ), t0( size ), t1( size ), t2( size );
    double sumweight = 0.0;
    for ( int i = 0; i < params.size(); i++ )
      if ( i != iw ) {
        double weight = 1.0;
        if ( type_ == GRADIENT ) {
          double r2 = 0.0;
          for ( int j = 0; j < size; j++ )
            r2 += clipper::Util::sqr( params[i][j] - params[iw][j] );
          weight = ( fn[iw] - fn[i] ) / sqrt( r2 );
        }
        for ( int j = 0; j < size; j++ )
          centr[j] += weight * params[i][j];
        sumweight += weight;
      }
    for ( int j = 0; j < size; j++ ) centr[j] = centr[j] / sumweight;
    for ( int j = 0; j < size; j++ ) shift[j] = centr[j] - params[iw][j];
    // calculate first trial point
    for ( int j = 0; j < size; j++ ) {
      t0[j] = centr[j] - 0.5 * shift[j];
      t1[j] = centr[j] + 1.0 * shift[j];
      t2[j] = centr[j] + 2.0 * shift[j];
    }
    f1 = target_fn( t1 );
    // simplex conditions
    STEP step = UNKN;
    if ( !clipper::Util::is_nan( f1 ) ) { // new point is valid
      if ( f1 < fn[iw] ) {    // new point is better than worst
        if ( f1 < fn[ib] ) {  // new point is better than best
          f2 = target_fn( t2 );
          if ( !clipper::Util::is_nan( f2 ) ) {
            if ( f2 < f1 ) {  // extended step is best
              step = EXTN;
            } else {  // normal step better than extended step
              step = NRML;
            }
          } else {  // extended step is invalid
            step = NRML;
          }
        } else {  // normal step is neither best nor worst (default)
          step = NRML;
        }
      }  // else normal step is worse
    }  // else normal step is invalid
    if ( step == UNKN ) {    // normal step is worse or invalid
      f0 = target_fn( t0 );  // contraction step (always valid)
      if ( f0 < fn[iw] ) {   // contraction step is better than worst
	step = CTRN;
      } else {               // contraction step is worse
	step = CTRX;
      }
    }
    if      ( step == EXTN ) { params[iw] = t2; fn[iw] = f2; }  // extension
    else if ( step == NRML ) { params[iw] = t1; fn[iw] = f1; }  // normal
    else if ( step == CTRN ) { params[iw] = t0; fn[iw] = f0; }  // contraction
    else {  // otherwise contract all towards best (slow)
      for ( int i = 0; i < params.size(); i++ )
	if ( i != ib ) {
	  for ( int j = 0; j < size; j++ )
	    params[i][j] = 0.5 * ( params[i][j] + params[ib][j] );
	  fn[i] = target_fn( params[i] );
	}
    }
    if ( debug_mode ) {  // DEBUG OUTPUT
      if      ( step == EXTN ) std::cout << "Extn-step\n";
      else if ( step == NRML ) std::cout << "Nrml-step\n";
      else if ( step == CTRN ) std::cout << "Ctrn-step\n";
      else                     std::cout << "Ctrx-step\n";
    }
  }
  return params[ib];
}
