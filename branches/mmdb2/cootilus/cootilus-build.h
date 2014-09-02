/*! \file cootilus-build.h nautilus library */
/* (C) 2012 Kevin Cowtan & University of York all rights reserved */


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


class Coot_nucleic_acid_build {
 public:
  //! Constructor: takes filename for nucleic acid library file
  Coot_nucleic_acid_build( std::string filename );
  //! Build or extend a model
  bool build( CMMDBManager* mmdb, const clipper::Xmap<float>& xmap, const clipper::Coord_orth& centre, double radius ) const;
 private:
  std::string filename_;
};

