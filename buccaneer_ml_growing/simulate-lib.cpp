/*! \file simulate-lib.cpp map simulation library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include "simulate-lib.h"

#include <algorithm>


MapSimulate::MapSimulate( int nresbins, int binmin )
{
  n_res_bins = nresbins;
  bin_min = binmin;
}


bool MapSimulate::operator()
  ( clipper::HKL_data<clipper::data32::F_sigF>& sim_f,
    clipper::HKL_data<clipper::data32::ABCD>& sim_hl,
    const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
    const clipper::HKL_data<clipper::data32::ABCD>& ref_hl,
    const clipper::HKL_data<clipper::data32::F_sigF>& wrk_f,
    const clipper::HKL_data<clipper::data32::ABCD>& wrk_hl ) const
{
  // get hkls
  const clipper::HKL_info& ref_hkl = ref_f.base_hkl_info();
  const clipper::HKL_info& wrk_hkl = wrk_f.base_hkl_info();
  clipper::HKL_info::HKL_reference_index ih;

  // get phases
  clipper::HKL_data<clipper::data32::Phi_fom> ref_phi( ref_hkl );
  clipper::HKL_data<clipper::data32::Phi_fom> wrk_phi( wrk_hkl );
  ref_phi.compute( ref_hl, clipper::data32::Compute_phifom_from_abcd() );
  wrk_phi.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );

  // make E's
  const int n_scl_param = 10;
  clipper::HKL_data<clipper::data32::E_sigE> ref_e( ref_hkl );
  clipper::HKL_data<clipper::data32::E_sigE> wrk_e( wrk_hkl );
  ref_e.compute( ref_f, clipper::data32::Compute_EsigE_from_FsigF() );
  wrk_e.compute( wrk_f, clipper::data32::Compute_EsigE_from_FsigF() );
  std::vector<double> params( n_scl_param, 1.0 );
  // scale reference data
  std::vector<double> ref_list_res;
  clipper::BasisFn_spline ref_basis( ref_e, n_scl_param, 2.0 );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> ref_target( ref_e );
  clipper::ResolutionFn ref_scale( ref_hkl, ref_basis, ref_target, params );
  for ( ih = ref_hkl.first(); !ih.last(); ih.next() )
    if ( !ref_e[ih].missing() ) {
      ref_e[ih].scale( sqrt( ref_scale.f(ih) ) );  // scale
      ref_list_res.push_back( ih.invresolsq() );   // store resol
    }
  // scale work data
  std::vector<double> wrk_list_res;
  clipper::BasisFn_spline wrk_basis( wrk_e, n_scl_param, 2.0 );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> wrk_target( wrk_e );
  clipper::ResolutionFn wrk_scale( wrk_hkl, wrk_basis, wrk_target, params );
  for ( ih = wrk_hkl.first(); !ih.last(); ih.next() )
    if ( !wrk_e[ih].missing() ) {
      wrk_e[ih].scale( sqrt( wrk_scale.f(ih) ) );  // scale
      wrk_list_res.push_back( ih.invresolsq() );   // store resol
    }

  // Now scale reference data to match work
  double cellscale = ref_hkl.cell().volume() / wrk_hkl.cell().volume();
  for ( ih = ref_hkl.first(); !ih.last(); ih.next() )
    if ( !ref_f[ih].missing() ) {
      double s = ih.hkl().invresolsq( ref_hkl.cell() );
      sim_f[ih] = ref_f[ih];
      sim_f[ih].scale( sqrt( cellscale
			     * ref_basis.f_s( s, ref_scale.params() )
			     / wrk_basis.f_s( s, wrk_scale.params() ) ) );
    }

  // make ordinals
  clipper::Generic_ordinal ref_ord_res;
  clipper::Generic_ordinal wrk_ord_res;
  ref_ord_res.init( ref_list_res );
  wrk_ord_res.init( wrk_list_res );

  std::vector<std::vector<int> > ref_acen(n_res_bins), ref_ccen(n_res_bins);
  std::vector<std::vector<int> > wrk_acen(n_res_bins), wrk_ccen(n_res_bins);

  // now bin the work structure foms by resolution
  for ( ih = ref_hkl.first(); !ih.last(); ih.next() )
    if ( !ref_e[ih].missing() ) {
      int ires = clipper::Util::intf( double( n_res_bins ) *
				      ref_ord_res.ordinal( ih.invresolsq() ) );
      if ( ih.hkl_class().centric() )
	ref_ccen[ires].push_back( ih.index() );
      else
	ref_acen[ires].push_back( ih.index() );
    }
  for ( ih = wrk_hkl.first(); !ih.last(); ih.next() )
    if ( !wrk_e[ih].missing() ) {
      int ires = clipper::Util::intf( double( n_res_bins ) *
				      wrk_ord_res.ordinal( ih.invresolsq() ) );
      if ( ih.hkl_class().centric() )
	wrk_ccen[ires].push_back( ih.index() );
      else
	wrk_acen[ires].push_back( ih.index() );
    }

  // pad the centric bins if too sparse
  std::vector<std::vector<int> > tmp_ccen = wrk_ccen;
  for ( int i = 0; i < n_res_bins; i++ ) {
    int j = 1;
    while ( tmp_ccen[i].size() < bin_min && j < n_res_bins ) {
      if ( i+j >= 0 && i+j < n_res_bins )
	tmp_ccen[i].insert( tmp_ccen[i].end(),
			    wrk_ccen[i+j].begin(), wrk_ccen[i+j].end() );
      j = ( j > 0 ) ? ( -j ) : ( -j + 1 );
    }
  }
  wrk_ccen = tmp_ccen;
  for ( int i = 0; i < n_res_bins; i++ ) {
    if ( wrk_ccen[i].size() == 0 ) wrk_ccen[i] = wrk_acen[i];
    if ( ref_ccen[i].size() == 0 ) ref_ccen[i] = ref_acen[i];
  }

  // sort fom's by magnitude
  for ( int i = 0; i < n_res_bins; i++ ) {
    std::sort( ref_acen[i].begin(), ref_acen[i].end(), EMagCompare( ref_e ) );
    std::sort( ref_ccen[i].begin(), ref_ccen[i].end(), EMagCompare( ref_e ) );
    std::sort( wrk_acen[i].begin(), wrk_acen[i].end(), EMagCompare( wrk_e ) );
    std::sort( wrk_ccen[i].begin(), wrk_ccen[i].end(), EMagCompare( wrk_e ) );
  }

  // Transfer foms to reference structure
  for ( int i = 0; i < n_res_bins; i++ ) {
    for ( int j = 0; j < ref_acen[i].size(); j++ ) {
      if ( wrk_acen[i].size() != 0 && ref_acen[i].size() != 0 ) {
	int k = clipper::Util::intf(double(j)*double(wrk_acen[i].size())
				             /double(ref_acen[i].size()));
	ref_phi[ref_acen[i][j]].fom() = wrk_phi[wrk_acen[i][k]].fom();
      } else {
	ref_phi[ref_acen[i][j]].set_null();
      }
    }
    for ( int j = 0; j < ref_ccen[i].size(); j++ ) {
      if ( wrk_ccen[i].size() != 0 && ref_ccen[i].size() != 0 ) {
	int k = clipper::Util::intf(double(j)*double(wrk_ccen[i].size())
				             /double(ref_ccen[i].size()));
	ref_phi[ref_ccen[i][j]].fom() = wrk_phi[wrk_ccen[i][k]].fom();
      } else {
	ref_phi[ref_ccen[i][j]].set_null();
      }
    }
  }

  // Make phase errors to go with FOMs
  for ( ih = ref_hkl.first(); !ih.last(); ih.next() ) {
    double dphi;
    double fom = 0.0;
    if ( !ref_phi[ih].missing() ) fom = ref_phi[ih].fom();
    if ( ih.hkl_class().centric() ) {  // deal with centrics
      double prob = 0.5 * ( 1.0 + fom );
      if ( (0.001*double(rand()%1000)) < prob )
	dphi = 0.0;
      else
	dphi = clipper::Util::pi();
    } else {                           // deal with acentrics
      int i;
      double x = clipper::Util::invsim( fom );
      for ( i = 0; i < 100; i++ ) {
	dphi = clipper::Util::twopi()*(0.001*double(rand()%1000));
	double prob = exp(x*(cos(dphi)-1.0));
	if ( (0.001*double(rand()%1000)) < prob ) break;
      }
      if ( i == 100 ) dphi = 0.0;
    }
    if ( !ref_phi[ih].missing() ) ref_phi[ih].phi() += dphi;
  }

  // calculate sim hl coeffs
  sim_hl.compute( ref_phi, clipper::data32::Compute_abcd_from_phifom() );

  return true;
}
