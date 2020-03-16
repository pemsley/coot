// Clipper app to perform shift field refinement
/* Copyright 2018 Kevin Cowtan & University of York all rights reserved */


#ifndef COOT_SHIFTFIELD
#define COOT_SHIFTFIELD


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


namespace coot {

   // update the temperature-factors of the atoms in mol
   void
   shift_field_b_factor_refinement(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &fo0,
                                   const clipper::HKL_data<clipper::data32::Flag> &free,
                                   mmdb::Manager *mol, int ncycles);

   // update the coordinates of the atoms in mol
   void shift_field_xyz_refinement(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &fo0,
                                   const clipper::HKL_data<clipper::data32::Flag> &free,
                                   mmdb::Manager *mol,
                                   float resolution);

}

#endif
