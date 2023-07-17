/*! \file simulate-lib.h map simulation library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */


#ifndef SIMULATE_LIB
#define SIMULATE_LIB

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>


//! Map simulation class
/* This class simulates a set of HL coeffs for a reference structure
   (A & B only) matching the properties of the coefficients of a work
   structure. */
class MapSimulate {
 public:
  MapSimulate( int nresbins = 100, int binmin = 20 );
  bool operator() ( clipper::HKL_data<clipper::data32::F_sigF>& sim_f,
		    clipper::HKL_data<clipper::data32::ABCD>& sim_hl,
		    const clipper::HKL_data<clipper::data32::F_sigF>& ref_f,
		    const clipper::HKL_data<clipper::data32::ABCD>& ref_hl,
		    const clipper::HKL_data<clipper::data32::F_sigF>& wrk_f,
		    const clipper::HKL_data<clipper::data32::ABCD>& wrk_hl ) const;
 private:
  class EMagCompare {
  public:
    EMagCompare(const clipper::HKL_data<clipper::data32::E_sigE>& m ) { p = &m; }
    bool operator() ( const int& i1, const int& i2 ) const { return (*p)[i1].E() < (*p)[i2].E(); }
    const clipper::HKL_data<clipper::data32::E_sigE>* p;
  };
  int n_res_bins, bin_min;
};


#endif
