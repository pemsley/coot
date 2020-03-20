// Header file for shift field refinement
/* Copyright 2018 Kevin Cowtan & University of York all rights reserved */


#ifndef SHIFTFIELD
#define SHIFTFIELD

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

class Shift_field_refine{

 public:

  static bool shift_field_coord(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,      clipper::Xmap<float>& x2map,      clipper::Xmap<float>& x3map,
    float rad, int filter );
  static bool shift_field_coord_const(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,      clipper::Xmap<float>& x2map,      clipper::Xmap<float>& x3map,
    float rad, int filter );

  static bool shift_field_u_iso(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,
    float rad, int filter );
  static bool shift_field_u_iso_const(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,
    float rad, int filter );

  static bool shift_field_u_aniso(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
		clipper::Xmap<float>& x1map, clipper::Xmap<float>& x2map, clipper::Xmap<float>& x3map,
    clipper::Xmap<float>& x4map, clipper::Xmap<float>& x5map, clipper::Xmap<float>& x6map, 
    float rad, clipper::Resolution& reso, int filter );
  static bool shift_field_u_aniso_const(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
		clipper::Xmap<float>& x1map, clipper::Xmap<float>& x2map, clipper::Xmap<float>& x3map,
    clipper::Xmap<float>& x4map, clipper::Xmap<float>& x5map, clipper::Xmap<float>& x6map, 
    float rad, clipper::Resolution& reso, int filter );

};

#endif

