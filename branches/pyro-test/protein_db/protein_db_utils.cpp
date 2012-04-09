/*! ProteinDB Top 500 main chain database */
/* (C) 2008-2009 Kevin Cowtan & University of York all rights reserved */


#include "protein_db_utils.h"

#include <algorithm>


namespace ProteinDB {


/*! Construct a density scoring object from the given Xmap.

  Fragments will then be scored based on the position of the densitied
  in a cumulative density distribution based on a Gaussian
  distribution derived from the the mean and variance of the map.

  \param xmap The electron density map.
*/
ScoreDensity::ScoreDensity( const clipper::Xmap<float>& xmap,
			    double sig1, double sig2 )
{
  xmap_ = &xmap;
  clipper::Map_stats map_stats( xmap );
  s1 = map_stats.mean() + sig1 * map_stats.std_dev();
  s2 =                    sig2 * map_stats.std_dev();
}


/*! Score a fragment based on the electron density map.

  The fragment will be scored based on the electron density at the
  nearest grid point to each mainchain (N,CA,C) atom in the
  fragment. The result is in the form of a log-probability.

  \param frag The fragment to score.
  \return The fragment score.
*/
double ScoreDensity::score( const Chain& frag ) const
{
  const clipper::Xmap<float>& xmap   = *xmap_;
  const clipper::Cell& cell          = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();
  double s = 0.0;
  for ( int r = 1; r < frag.size()-1; r++ ) {
    const float r1 = xmap.get_data( frag[r].coord_n ().coord_frac(cell).coord_grid(grid) );
    const float r2 = xmap.get_data( frag[r].coord_ca().coord_frac(cell).coord_grid(grid) );
    const float r3 = xmap.get_data( frag[r].coord_c ().coord_frac(cell).coord_grid(grid) );
    s += ( log( phi_approx( ( r1 - s1 ) / s2 ) ) +
	   log( phi_approx( ( r2 - s1 ) / s2 ) ) +
	   log( phi_approx( ( r3 - s1 ) / s2 ) ) );
  }
  return s;
}


/*! Cumulative distribution function for a Gaussian
  \param z
  \return The probability
*/
double ScoreDensity::phi_approx( double z )
{
  if ( z <= 0.0 )
    return     (exp(-0.5*z*z)/(1.2533141373*(-z+sqrt(z*z+2.546479089470))));
  else
    return 1.0-(exp(-0.5*z*z)/(1.2533141373*( z+sqrt(z*z+2.546479089470))));
}


/*! Construct a clash scoring object from an atomic model.

  Fragments will then be scored based on the number of atoms within a
  given radius of atoms in a reference molecule.

  \param coords The coordinates to check for clashes
  \param rad The clash radius
*/
ScoreClashes::ScoreClashes( const std::vector<clipper::Coord_orth>& coords, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double rad )
{
  rad_ = rad;

  // create non-neighbour bonding object
  mol.init( spgr, cell );
  clipper::MPolymer mp;
  clipper::MMonomer mm;
  clipper::MAtom    ma = clipper::MAtom::null();
  ma.set_id( " C  " );
  ma.set_u_iso( 1.0 );
  ma.set_occupancy( 1.0 );
  for ( int i = 0; i < coords.size(); i++ ) {
    ma.set_coord_orth( coords[i] );
    mm.insert( ma );
  }
  mp.insert( mm );
  mol.insert( mp );
  nnb = clipper::MAtomNonBond( mol, 4.0 );
}


/*! Exclude a set of atoms from the clash detection score.

  When fitting a fragment or loop, the new loop will always score
  clashes against the atoms that it is replacing or those it is bonded
  too. These atoms may be excluded from the clash score by supplying
  the original fragment which is being replaced to this function.

  The exclusion radius should be at least as large as the clash
  radius. If zero, or omitted, the default is clash radius + 0.5.

  \param coords The atoms to exclude.
  \param rad    The radius to exclude.
*/
void ScoreClashes::set_exclude( const std::vector<clipper::Coord_orth>& coords, double rad )
{
  if ( mol.size()    != 1 ) clipper::Message::message(
      clipper::Message_fatal( "ScoreClashes: No model to exclude" ) );
  if ( mol[0].size() != 1 ) clipper::Message::message(
      clipper::Message_fatal( "ScoreClashes: Internal error" ) );
  // get radius
  if ( rad == 0.0 ) rad = rad_ + 0.5;
  // flag all atoms as included
  for ( int a = 0; a < mol[0][0].size(); a++ ) 
    mol[0][0][a].set_occupancy( 1.0 );
  // do a clash search and zero clashing occupancies
  for ( int a = 0; a < coords.size(); a++ ) {
    std::vector<clipper::MAtomIndexSymmetry> atoms = nnb( coords[a], rad );
    for ( int i = 0; i < atoms.size(); i++ )
      mol[atoms[i].polymer()]
	 [atoms[i].monomer()]
	 [atoms[i].atom()].set_occupancy( 0.0 );
  }
}


/*! Exclude a set of atoms from the clash detection score.

  When fitting a fragment or loop, the new loop will always score
  clashes against the atoms that it is replacing or those it is bonded
  too. These atoms may be excluded from the clash score by supplying
  the original fragment which is being replaced to this function.

  \param frag The atoms to exclude.
  \param rad  The radius to exclude.
*/
void ScoreClashes::set_exclude( const Chain& frag, double rad )
{
  std::vector<clipper::Coord_orth> coords;
  for ( int r = 0; r < frag.size(); r++ ) {
    if ( frag[r].flag() == Residue::NORMAL ) {
      coords.push_back( frag[r].coord_n() );
      coords.push_back( frag[r].coord_ca() );
      coords.push_back( frag[r].coord_c() );
    } else if ( frag[r].flag() == Residue::CALPHA ) {
      coords.push_back( frag[r].coord_ca() );
    }
  }
  set_exclude( coords, rad );
}


/*! Score a fragment based on clashes

  The fragment will be scored based on the number of clashes with
  atoms in the base model.

  \param frag The fragment to score.
  \return The fragment score.
*/
double ScoreClashes::score( const Chain& frag ) const
{
  double s = 0.0;
  for ( int r = 0; r < frag.size(); r++ ) {
    std::vector<clipper::MAtomIndexSymmetry>
      atoms = nnb( frag[r].coord_ca(), rad_ );
    for ( int i = 0; i < atoms.size(); i++ )
      if ( mol[atoms[i].polymer()]
	      [atoms[i].monomer()]
	      [atoms[i].atom()].occupancy() > 0.5 ) s += -1.0;
  }
  return s;
}


/*! Search using predefined density and clash scoring classes.

  Combined search, scoring and sorting of fragments. This version
  assumes you have already initialised Score classes, for performing
  multiple searches against the same map or model.

  \param frag The fragment to search for.
  \param nfrag The maximum  number of fragments to return
  \param score_rho The density scoring class.
  \param score_cls The clash scoring class.
  \param wdense (optional) weight for the density score
  \param wclash (optional) weight for the clash score
  \return A vector of chains, with the first chain representing the best score.
*/
  std::vector<Chain> ProteinDBSearch::search( const Chain& frag, const int nfrag, ScoreDensity& score_rho, ScoreClashes& score_cls, const double wdense, double wclash )
{
  std::vector<Chain> result;

  ChainDB fragdb( frag );
  std::vector<Chain> frags = chaindb.match_fragment( fragdb, nfrag );

  if ( frags.size() > 0 ) {

    // exclude the initial fragment environment from clash scoring
    score_cls.set_exclude( frag );

    // score vs density and clashes (omitting first and last residues)
    std::vector<std::pair<double,int> > fragscore( frags.size() );
    for ( int f = 0; f < frags.size(); f++ ) {
      // density score
      double scr_rho = score_rho.score( frags[f] );
      // clash score
      double scr_cls = score_cls.score( frags[f] ); 
      // store
      double s = wdense * scr_rho + wclash * scr_cls;
      fragscore[f] = std::pair<double,int>( -s, f );
    }

    // sort
    std::sort( fragscore.begin(), fragscore.end() );

    result.resize( frags.size() );
    score.resize( frags.size() );
    for ( int f = 0; f < frags.size(); f++ ) {
      result[f] = frags[fragscore[f].second];
      score[f]  = -fragscore[f].first;
    }
  }
  return result;
}


/*! Search from scratch.

  Combined search, scoring and sorting of fragments. This version is
  for general use.

  \param frag The fragment to search for.
  \param nfrag The maximum  number of fragments to return
  \param xmap The electron density map for scoring models.
  \param coords The coordinates of the model for use in clash scoring.
  \param wdense (optional) weight for the density score
  \param wclash (optional) weight for the clash score
  \param sig1 (optional) the sigma offset for density scores
  \param sig2 (optional) the sigma weight for density scores
  \param clashrad (optional) the radius for clash penalties
  \return A vector of chains, with the first chain representing the best score.
*/
  std::vector<Chain> ProteinDBSearch::search( const Chain& frag, const int nfrag, const clipper::Xmap<float>& xmap, const std::vector<clipper::Coord_orth>& coords, double wdense, double wclash, double sig1, double sig2, double clashrad )
{
  ScoreDensity score_rho( xmap, sig1, sig2 );
  ScoreClashes score_cls( coords, xmap.spacegroup(), xmap.cell(), clashrad );
  return search( frag, nfrag, score_rho, score_cls, wdense, wclash );
}


} // namespace ProteinDB
