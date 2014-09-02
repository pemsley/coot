/*! ProteinDB Top 500 main chain database */
/* (C) 2008-2009 Kevin Cowtan & University of York all rights reserved */


#ifndef PROTEIN_DB_UTILS
#define PROTEIN_DB_UTILS


#include "protein_db.h"


namespace ProteinDB {


/*!
  Utility class for evaluating ProteinDB search results against density
*/
  class ScoreDensity {
  public:
    //! Constructor: from Xmap
    ScoreDensity( const clipper::Xmap<float>& xmap, double sig1=0.0, double sig2=2.0 );
    //! Score fragment against xmap density
    double score( const Chain& frag ) const;
  private:
    static double phi_approx( double z );
    const clipper::Xmap<float>* xmap_;
    double s1, s2;
  };


/*!
  Utility class for evaluating ProteinDB search results for clashes
*/
  class ScoreClashes {
  public:
    //! Constructor: from model
    ScoreClashes( const std::vector<clipper::Coord_orth>& coords, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double rad=2.0 );
    //! Set a region to exclude from clash checking
    void set_exclude( const std::vector<clipper::Coord_orth>& coords, double rad=0.0 );
    //! Set a region to exclude from clash checking
    void set_exclude( const Chain& frag, double rad=0.0 );
    //! Score fragment for clashes against model
    double score( const Chain& frag ) const;
  private:
    clipper::MiniMol mol;
    clipper::MAtomNonBond nnb;
    double rad_;
  };


/*!
  Utility class for producing a fully scored list of fragments
*/
  class ProteinDBSearch {
  public:
    //! Constructor
    ProteinDBSearch( const clipper::String file ) { chaindb.init( file ); }
    //! Search, score and sort fragments
    std::vector<Chain> search( const Chain& frag, const int nfrag, ScoreDensity& score_rho, ScoreClashes& score_cls, const double wdense=1.0, double wclash=1.0 );
    //! Search, score and sort fragments
    std::vector<Chain> search( const Chain& frag, const int nfrag, const clipper::Xmap<float>& xmap, const std::vector<clipper::Coord_orth>& coords, double wdense=1.0, double wclash=1.0, double sig1=0.0, double sig2=2.0, double clashrad=2.0 );
    //! Return scores to match last set of fragments
    std::vector<double> scores() const { return score; }
  private:
    ChainDB chaindb;
    std::vector<double> score;
  };


} // namespace ProteinDB


#endif
