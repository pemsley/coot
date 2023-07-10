/*! \file buccaneer-lib.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_LIB
#define BUCCANEER_LIB

#include <clipper/clipper.h>


//! Score list class.
/*! The class keeps a list of the best fits to some target function,
  scoring both the score and an object of a user defined type
  representing the parameters whuich gave that score. On construction,
  a maximum size is specified. Only the best N matches will be
  stored.

  The implementation could be more efficient. */
template<class T> class Score_list {
 public:
  Score_list() {}  //!< null constructor
  //! constructor: takes maximum size of list
  Score_list( const int& n ) { init( n ); }
  //! initialiser: takes maximum size of list
  void init( const int& n ) { max = n; list.clear(); list.reserve( max ); }
  //! is a given score good enough to be added to the list?
  bool addable( const clipper::ftype& score ) const
    { return ( list.size() < max || score < list.back().first ); }
  //! add a score to the list, if it is good enough
  void add( const clipper::ftype& scr, const T& data ) {
    if ( addable( scr ) ) {
      if ( list.size() >= max )	list.pop_back();
      int i; for ( i = list.size()-1; i >= 0; i-- ) if ( scr > score(i) ) break;
      list.insert( list.begin()+(i+1), std::pair<clipper::ftype,T>(scr,data) );
    }
  }
  //! delete a score from the list
  void del( const int& i ) { list.erase( list.begin()+i ); }
  //! access score
  const clipper::ftype& score( const int& i ) const { return( list[i].first ); }
  //! access list
  const T& operator[] ( const int& i ) const { return( list[i].second ); }
  //! list size
  int size() const { return list.size(); }
 private:
  int max;
  std::vector<std::pair<clipper::ftype,T> > list;
};


//! Log-likelihood map matching target
/*! This class is used in determining the log-likelihood fit of some
  desired density features from some region of a noisy electron
  density map. It contains methods to accumulate the log likelihood
  target from a number of sample density regions from a map with
  similar noise levels to the target map; methods for a FFFear
  6-dimensional (rotation/orientation) search of the target map; and
  methods for testing indivdual sample positions and orientations.

  Note the results from the 6-d search and the fast and full LLK
  calculations are on different scales and so cannot be compared
  directly. */
class LLK_map_target {
 public:
  enum TYPE { NORMAL, CORREL };
  //! null constructor
  LLK_map_target() {}
  //! constructor: provide radius and sampling in A for LLK target
  LLK_map_target( const clipper::ftype& rad, const clipper::ftype& sampling, TYPE type = NORMAL ) { init( rad, sampling, type ); }
  //! initialiser: provide radius and sampling in A for LLK target
  void init( const clipper::ftype& rad, const clipper::ftype& sampling, TYPE type = NORMAL );
  //! prepare LLK target after accumulating density or loading map
  void prep_llk();
  //! accumulate density statistics from a sample orientation in a sample map
  void accumulate( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop );

  //! set the scaling for llk and llk_approx
  void set_scale( const double& scale, const double& offset )
    { tgt_scl = scale; tgt_off = offset; }

  //! return a list of RTops
  std::vector<clipper::RTop_orth> rtop_list( const clipper::Spacegroup& spgr, const clipper::ftype& step ) const;
  //! perform 6-d (RT) search for LLK target in given map with prior
  void search( clipper::Xmap<float>& resultscr, clipper::Xmap<int>& resultrot, clipper::Xmap<int>& resulttrn, const clipper::Xmap<float>& xmap, const std::vector<clipper::RTop_orth>& rtops ) const;

  //! calculate fast approx to LLK for given orientation
  clipper::ftype llk_approx( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const { return fasttgt.target( xmap, rtop )*tgt_scl+tgt_off; }
  //! calculate full LLK for given orientation
  clipper::ftype llk       ( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const { return slowtgt.target( xmap, rtop )*tgt_scl+tgt_off; }
  //! output formatted representation
  clipper::String format() const;

  /*! Class to hold sampled values from LLK map target */
  class Sampled {
  public:
    //! null constructor
    Sampled() { type_ = LLK_map_target::NORMAL; }
    //! insert values
    void insert( clipper::Coord_orth coord,
		 clipper::ftype tgt, clipper::ftype wgt );
    //! evaluate llk
    clipper::ftype llk( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const;
    //! evaluate correlation
    clipper::ftype correl( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const;
    //! evaluate target function
    clipper::ftype target( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const {
      if ( type_ == LLK_map_target::NORMAL ) return llk(xmap,rtop);
      else                                   return correl(xmap,rtop);
    }
    int size() const { return repxyz.size(); }  //!< get number of samples
    void set_type( LLK_map_target::TYPE t ) { type_ = t; }
    clipper::Coord_orth coord_orth( int i ) const { return repxyz[i]; }
    clipper::ftype      target( int i )     const { return reptgt[i]; }
    clipper::ftype      weight( int i )     const { return repwgt[i]; }
  private:
    std::vector<clipper::Coord_orth> repxyz;  //!< fast target lists
    std::vector<clipper::ftype>      reptgt;
    std::vector<clipper::ftype>      repwgt;
    LLK_map_target::TYPE type_;
  };

  //! generate a distribution of LLK values for a given map
  void prep_llk_distribution( const clipper::Xmap<float>& xmap );
  //! return llk value by position in cumulative distribution
  clipper::ftype llk_distribution( const clipper::ftype& ordinal ) const;

  //! access to fine sampled target function
  const Sampled& sampled() const { return slowtgt; }
  //! access to LLK target for load/save
  const clipper::NXmap<float>& llk_target() const { return target; }
  //! access to LLK weight for load/save
  const clipper::NXmap<float>& llk_weight() const { return weight; }
  //! access to LLK target for load/save
  clipper::NXmap<float>& llk_target() { return target; }
  //! access to LLK weight for load/save
  clipper::NXmap<float>& llk_weight() { return weight; }
  //! access to number of samples for load/save
  int& num_samples() { return naccum; }
 private:
  clipper::ftype radius;  //!< density sphere radius
  int naccum;             //!< number of maps accumulated
  clipper::NXmap<float> target;  //!< target map
  clipper::NXmap<float> weight;  //!< weight map
  Sampled slowtgt, fasttgt;
  LLK_map_target::TYPE type_;
  double tgt_off, tgt_scl;
  std::vector<clipper::ftype> llkdist;  //!< distrn of llk vals for given xmap
};


#endif
