/*! \file buccaneer-util.cpp buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-util.h"
#include "buccaneer-prot.h"
#include "buccaneer-tidy.h"

#include <fstream>
extern "C" {
#include <stdlib.h>
#include <stdio.h>
}


void BuccaneerUtil::set_reference( clipper::String& mtz, clipper::String& pdb )
{
  const char* clibdptr = getenv( "CLIBD" );
  if ( clibdptr != NULL ) {
    clipper::String clibd( clibdptr );
    clipper::String path;
    std::ifstream file;
    if ( pdb == "NONE" ) {
      path = clibd+"/reference_structures/reference-1tqw.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = clibd+"\\reference_structures\\reference-1tqw.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( mtz == "NONE" ) {
      path = clibd+"/reference_structures/reference-1tqw.mtz";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) mtz = path;
    }
    if ( mtz == "NONE" ) {
      path = clibd+"\\reference_structures\\reference-1tqw.mtz";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) mtz = path;
    }
    if ( pdb == "NONE" || mtz == "NONE" ) 
      clipper::Message::message( clipper::Message_fatal( "No reference data specified and not in $CLIBD" ) );
  } else {
    clipper::Message::message( clipper::Message_fatal( "No reference data specified and $CLIBD not found" ) );
  }
}


void BuccaneerUtil::read_model( clipper::MiniMol& mol, clipper::String file, bool verbose )
{
  const int mmdbflags = ( ::mmdb::MMDBF_IgnoreBlankLines |
                          ::mmdb::MMDBF_IgnoreDuplSeqNum |
                          ::mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                          ::mmdb::MMDBF_IgnoreRemarks );
  clipper::MMDBfile mmdb;
  mmdb.SetFlag( mmdbflags );
  if ( file != "NONE" ) {
    std::vector<clipper::String> files = file.split(",");
    for ( int f = 0; f < files.size(); f++ ) {
      try {
        clipper::MiniMol moltmp;
        mmdb.read_file( files[f] );
        mmdb.import_minimol( moltmp );
        clipper::Atom_list atoms = moltmp.atom_list();
        std::cout << "PDB file: " << files[f] << std::endl;
        std::cout << "  Number of atoms read: " << atoms.size() << std::endl;
        for ( int c = 0; c < moltmp.size(); c++ )
          if ( moltmp[c].id() != "!" ) mol.insert( moltmp[c] );        
        if ( verbose ) for ( int i = 0; i < atoms.size(); i += atoms.size()-1 ) printf("%i6  %4s  %8.3f %8.3f %8.3f\n", i, atoms[i].element().c_str(), atoms[i].coord_orth().x(), atoms[i].coord_orth().y(), atoms[i].coord_orth().z() );
      } catch ( clipper::Message_fatal ) {
        std::cout << "FAILED TO READ PDB FILE: " << file << std::endl;
      }
    }
  }
}


#ifdef BUCCANEER_PROFILE
#include <sys/times.h>
#include <unistd.h>
void BuccaneerLog::log( const clipper::String& id )
{
  int i;
  tms tmst;
  times( &tmst );
  long ticks = sysconf(_SC_CLK_TCK);
  double cpu = double( tmst.tms_utime ) / double( ticks );
  double elapsed = cpu - currentcpu;
  if ( id != "" ) {
    for ( i = 0; i < prof.size(); i++ )
      if ( id == prof[i].first ) break;
    if ( i < prof.size() )
      prof[i].second += elapsed;
    else
      prof.push_back( std::pair<std::string,double>( id, elapsed ) );
  }
  currentcpu = cpu;
}
#else
void BuccaneerLog::log( const clipper::String& id ) {}
#endif


void BuccaneerLog::log( const clipper::String& id, const clipper::MiniMol& mol, bool view )
{
  if ( view ) {
    for ( int c = 0; c < mol.size(); c++ )
      for ( int r = 0; r < mol[c].size(); r++ ) {
        clipper::Coord_orth co(0.0,0.0,0.0);
        for ( int a = 0; a < mol[c][r].size(); a++ )
          co += mol[c][r][a].coord_orth();
        co = (1.0/mol[c][r].size()) * co;
        std::cout << id << " " << c << "\t" << r << "\t" << co.format() << "\n";
        std::cout << id << " " << mol[c][r].type() << " ";
        for ( int a = 0; a < mol[c][r].size(); a++ )
          std::cout << mol[c][r][a].id() << " ";
        std::cout << std::endl;
        int cn = mol[c][r].lookup( " N  ", clipper::MM::ANY );
        int ca = mol[c][r].lookup( " CA ", clipper::MM::ANY );
        int cc = mol[c][r].lookup( " C  ", clipper::MM::ANY );
        if ( ca >= 0 && cn >= 0 ) {
          double d2 = ( mol[c][r][ca].coord_orth() -
                        mol[c][r][cn].coord_orth() ).lengthsq();
          if ( d2 > 6.25 )
            std::cout << "BOND N-CA " << d2 << std::endl;
        }
        if ( ca >= 0 && cc >= 0 ) {
          double d2 = ( mol[c][r][ca].coord_orth() -
                        mol[c][r][cc].coord_orth() ).lengthsq();
          if ( d2 > 6.25 )
            std::cout << "BOND CA-C " << d2 << std::endl;
        }
      }
  }
  log( id );
}

std::vector<double>  BuccaneerLog::evaluate(const clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq){

	log(mol,  mol_mr, seq);
std::vector<double>  nunqnseqnfrgs;
nunqnseqnfrgs.push_back(data[data.size()-1].nunq);
nunqnseqnfrgs.push_back(data[data.size()-1].nseq);
nunqnseqnfrgs.push_back(data[data.size()-1].nfrgs);
nunqnseqnfrgs.push_back(data[data.size()-1].nres);
nunqnseqnfrgs.push_back(data[data.size()-1].nmax);
nunqnseqnfrgs.push_back(data[data.size()-1].cres);
nunqnseqnfrgs.push_back(data[data.size()-1].cchn);
	return nunqnseqnfrgs;

}

clipper::String BuccaneerLog::log( const clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq )
{
  clipper::MiniMol mol_wrk = mol;
  std::vector<int> seqnums = ModelTidy::chain_renumber( mol_wrk, seq );
  std::vector<int> chnnums = ModelTidy::chain_assign( mol_wrk, mol_mr, seqnums, 1.0, 12 );
  // model stats
  int nres, nseq, nmax, nunq, nrch, nfrgs, nchns, nseqs;
  nfrgs = mol_wrk.size();
  nseqs = nchns = nres = nseq = nmax = nunq = nrch = 0;
  // find out how many sequences and chains are present
  for ( int c = 0; c < mol_wrk.size(); c++ ) {
    nseqs = std::max( nseqs, seqnums[c]+1 );
    nchns = std::max( nchns, chnnums[c]+1 );
  }
  // count residues
  for ( int c = 0; c < mol_wrk.size(); c++ ) {
    if ( mol_wrk[c].size() > nmax ) nmax = mol_wrk[c].size();
    for ( int r = 0; r < mol_wrk[c].size(); r++ ) {
      if ( mol_wrk[c][r].lookup( " CA ", clipper::MM::ANY ) >= 0 ) nres++;
      if ( ProteinTools::residue_index_3( mol_wrk[c][r].type() ) >= 0 ) nseq++;
    }
  }
  // count chain completeness
  std::vector<std::vector<int> > seqcounts( nchns );
  for ( int c = 0; c < mol_wrk.size(); c++ )
    if ( chnnums[c] >= 0 && seqnums[c] >= 0 )
      seqcounts[chnnums[c]].resize( seq[seqnums[c]].sequence().length(), 0 );
  for ( int c = 0; c < mol_wrk.size(); c++ ) {
    int nc = chnnums[c];
    if ( nc >= 0 ) {
      for ( int r = 0; r < mol_wrk[c].size(); r++ )
        if ( mol_wrk[c][r].lookup( " CA ", clipper::MM::ANY ) >= 0 ) {
          int nr = mol_wrk[c][r].seqnum() - 1;
          if ( nr >= 0 && nr < seqcounts[nc].size() ) seqcounts[nc][nr]++;
        }
    }
  }
  for ( int nc = 0; nc < seqcounts.size(); nc++ )
    for ( int nr = 0; nr < seqcounts[nc].size(); nr++ ) {
      nrch++;
      if ( seqcounts[nc][nr] == 1 ) nunq++;
    }
  double cres = double( nunq ) / double( std::max(nres,1) );
  double cchn = double( nunq ) / double( std::max(nrch,1) );
  /*
  for ( int nc = 0; nc < seqcounts.size(); nc++ ) {
    std::cout << nc << " ";
    for ( int nr = 0; nr < seqcounts[nc].size(); nr++ )
      std::cout << std::min( seqcounts[nc][nr],9 );
    std::cout << std::endl;
  }
  */
  // store
  cycdat dat;
  dat.nfrgs = nfrgs; dat.nchns = nchns;
  dat.nseq = nseq; dat.nres = nres; dat.nmax = nmax; dat.nunq = nunq;
  dat.cres = cres;
  dat.cchn = cchn;
  data.push_back( dat );

  // standard output
  char s[1000];
  sprintf( s, " %6i residues were built in %3i fragments, the longest having %4i residues.\n %6i residues were sequenced, after pruning.\n %6i residues were uniquely allocated to %3i chains.\n  Completeness by residues built:   %5.1f%%\n  Completeness of chains (number):  %5.1f%%    (%i)\n", nres, nfrgs, nmax, nseq, nunq, nchns, 100.0*cres, 100.0*cchn, nchns );
  return clipper::String(s);
}


void BuccaneerLog::xml( const clipper::String& xml ) const
{
  // xml output
  clipper::String xmltmp = xml+".tmp";
  FILE *f = fopen( xmltmp.c_str(), "w" );
  if ( f == NULL ) clipper::Message::message( clipper::Message_fatal( "Error: Could not open xml temporary file: "+xmltmp ) );
  fprintf( f, "<BuccaneerResult>\n" );
  fprintf( f, " <Title>%s</Title>\n", title_.c_str() );
  fprintf( f, " <Cycles>\n" );
  for ( int c = 0; c < data.size(); c++ ) {
    fprintf( f, "  <Cycle>\n" );
    fprintf( f, "   <Number>%i</Number>\n", c+1 );
    fprintf( f, "   <CompletenessByResiduesBuilt>%f</CompletenessByResiduesBuilt>\n", data[c].cres );
    fprintf( f, "   <CompletenessByChainsBuilt>%f</CompletenessByChainsBuilt>\n", data[c].cchn );
    fprintf( f, "   <ChainsBuilt>%i</ChainsBuilt>\n", data[c].nchns );
    fprintf( f, "   <FragmentsBuilt>%i</FragmentsBuilt>\n", data[c].nfrgs );
    fprintf( f, "   <ResiduesUnique>%i</ResiduesUnique>\n", data[c].nunq );
    fprintf( f, "   <ResiduesBuilt>%i</ResiduesBuilt>\n", data[c].nres );
    fprintf( f, "   <ResiduesSequenced>%i</ResiduesSequenced>\n", data[c].nseq );
    fprintf( f, "   <ResiduesLongestFragment>%i</ResiduesLongestFragment>\n", data[c].nmax );
    fprintf( f, "  </Cycle>\n" );
  }
  fprintf( f, " </Cycles>\n" );
  fprintf( f, " <Final>\n" );
  int c = data.size()-1;
  fprintf( f, "  <CompletenessByResiduesBuilt>%f</CompletenessByResiduesBuilt>\n", data[c].cres );
  fprintf( f, "  <CompletenessByChainsBuilt>%f</CompletenessByChainsBuilt>\n", data[c].cchn );
  fprintf( f, "  <ChainsBuilt>%i</ChainsBuilt>\n", data[c].nchns );
  fprintf( f, "  <FragmentsBuilt>%i</FragmentsBuilt>\n", data[c].nfrgs );
  fprintf( f, "  <ResiduesUnique>%i</ResiduesUnique>\n", data[c].nunq );
  fprintf( f, "  <ResiduesBuilt>%i</ResiduesBuilt>\n", data[c].nres );
  fprintf( f, "  <ResiduesSequenced>%i</ResiduesSequenced>\n", data[c].nseq );
  fprintf( f, "  <ResiduesLongestFragment>%i</ResiduesLongestFragment>\n", data[c].nmax );
  fprintf( f, " </Final>\n" );
  fprintf( f, "</BuccaneerResult>\n" );
  fclose(f);
  rename( xmltmp.c_str(), xml.c_str() );
}


void BuccaneerLog::profile()
{
  if ( prof.size() > 0 ) {
    std::cout << std::endl << "Profile:" << std::endl;
    for ( int i = 0; i < prof.size(); i++ )
      std::cout << prof[i].first << ": " << clipper::String( prof[i].second, 8 ) << " s" << std::endl;
  }
}
