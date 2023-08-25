/*! ProteinDB Top 500 main chain database */
/* (C) 2008-2009 Kevin Cowtan & University of York all rights reserved */


#ifndef PROTEIN_DB
#define PROTEIN_DB


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


namespace ProteinDB {


/*!
  Class for storing compact amino acid main chain information for use
  with Top500 DB. The same class is used for storing residues in the
  DB, and also for storing residues for searching against the DB.
  This class stores all the information concerning a single residue.
*/
class Residue {
 public:
  enum FLAG {NONE,NORMAL,CALPHA};
  //! Null constructor
  Residue() : flg(char(NONE)) {}
  //! Constructor: from C-alpha coord and residue type
  Residue( clipper::Coord_orth& ca, const clipper::String& type );
  //! Constructor: from main chain coords and residue type
  Residue( clipper::Coord_orth& cn, clipper::Coord_orth& ca, clipper::Coord_orth& cc, const clipper::String& type );
  //! Constructor: from MMonomer
  Residue( const clipper::MMonomer& mm );
  clipper::Coord_orth coord_n () const;  //!< N atom coordinate
  clipper::Coord_orth coord_ca() const;  //!< C-alpha atom coordinate
  clipper::Coord_orth coord_c () const;  //!< C atom coordinate
  clipper::MMonomer mmonomer() const;    //!< Build MMonomer 
  //! Transform by rotation-translation operator
  void transform( const clipper::RTop_orth& rtop );
  //! Merge other residue coords with this one using given weights
  bool merge( const Residue& other, const double wn, const double wa, const double wc );
  char type() const { return typ; }           //!< return residue type (1-char)
  FLAG flag() const { return FLAG(flg); }     //!< return flag (null,Ca,all)
  void set_type( char t ) { typ = t; }        //!< set residue type (1-char)
  void set_flag( FLAG t ) { flg = char(t); }  //!< set flag (null,Ca,all)
  void data_import( const char* d );          //!< import from char array
  void data_export( char* d ) const;          //!< export to char array
  bool is_null() const { return flg == char(NONE); }  //!< test for null
  //! return 1-char residue type from 1 or 3 char string. Space on error.
  static char residue_type( const clipper::String& type );

  /*! Class for describing a residue type mask, used to describe a
      list of allowed residue types. */
  class TypeMask {
  public:
    //! Null constructor
    TypeMask() : msk(0) {}
    //! Initialise from a residue code, '?' for all
    explicit TypeMask( const char t ) { msk = msks[t&0x1f]; }
    friend TypeMask operator| ( TypeMask t1, TypeMask t2 ) { return TypeMask( t1.mask() | t2.mask() ); }
    friend TypeMask operator& ( TypeMask t1, TypeMask t2 ) { return TypeMask( t1.mask() & t2.mask() ); }
    friend TypeMask operator! ( TypeMask t1 ) { return TypeMask( !t1.mask() ); }
    //! Return mask. Use '&' to test if a mask matches a type.
    const int& mask() const { return msk; }
  private:
    explicit TypeMask( int i ) : msk(i) {}
    static const int msks[32];
    int msk;
  };

 private:
  static const int ntype = 22;
  static const char rtype1[ntype], rtype3[ntype][4];
  static void unpack_float( const char* d, float& f ) { const short s = ((short(d[0])<<8)&0xFF00) | ((short(d[1]))&0x00FF); f = float(s)/100.0; }
  static void pack_float( char* d, const float& f ) { const short s = rint(100.0*f); d[0] = char((s>>8)&0x00FF); d[1] = char((s)&0x00FF); }
  float nnx, nny, nnz, cax, cay, caz, ccx, ccy, ccz;
  char typ, flg;
};


/*!
  Class for storing compact amino acid main chain information for use
  with Top500 DB. The same class is used for storing residues in the
  DB, and also for storing residues for searching against the DB.
  This class stores a list of residues representing either a complete
  chain or chains (in the DB) or a complete search fragment to be
  searched against the DB.
*/
class Chain {
 public:
  //! Null constructor
  Chain() {}
  //! Add a pdb file to this DBchain (NOTE: CHAINS ARE SHIFTED TO ORIGIN)
  bool add_pdb( const clipper::String file );
  //! Add a residue to this DBchain
  void add_residue( const Residue& r ) { dbresidues.push_back( r ); }
  //! Export DB to binary file
  bool save_db( const clipper::String file ) const;
  //! Import DB from binary file
  bool load_db( const clipper::String file );
  //! Import DB from binary file
  //! Merge other residue coords with this one using given weights
  bool merge( const Chain& other, const std::vector<double>& wgt );
  //! extract fragment of given length
  Chain extract( int offset, int len ) const;
  //! check if fragment is continuous
  bool is_continuous() const;
  //! Transform by rotation-translation operator
  void transform( const clipper::RTop_orth& rtop );
  //! get RTop to fit DB fragment to given fragment
  void lsq_superpose( const Chain& frag );
  //! get RTop to fit DB fragment to given fragment with weights
  void lsq_superpose( const Chain& frag, const std::vector<double>& wgts );
  //! Get rmsd versus other fragment
  double rmsd( const Chain& other ) const;
  //! Get rmsd versus other fragment
  double rmsd( const Chain& other, const std::vector<double>& wgts ) const;  
  //! Get residue by position in list
  const Residue& operator[] ( const int& i ) const { return dbresidues[i]; }
  //! Set residue by position in list
  Residue& operator[] ( const int& i ) { return dbresidues[i]; }
  //! Get number of residues in list
  int size() const { return dbresidues.size(); }
  // output some debug info
  void debug() const;
 protected:
  std::vector<Residue> dbresidues;
};


/*!
  Class for storing compact amino acid main chain information for
  use with Top500 DB. The same class is used for storing residues in
  the DB, and also for storing residues for searching against the DB.
  This class is an extension of Chain, which also adds the fast
  distance matrix entries. These are filled out by calling the
  calc_distances() method after the chain has been completely
  assembled.
*/
class ChainDB : public Chain {
 public:
  static const int ndist = 20;
  struct DistVec { float data[ChainDB::ndist]; };
  //! Null constructor
  ChainDB() {}
  //! Constructor: from Chain
  ChainDB( const Chain& chain ) : Chain( chain ) { calc_distances(); }
  //! Constructor: from binary DB file name
  ChainDB( const clipper::String file ) { init( file ); }
  //! Initialiser: from binary DB file name
  void init( const clipper::String file );
  //! calculate the distance matrix elements
  void calc_distances();
  //! score a fragment against the nth fragment in the DB
  double score_distance( const ChainDB& frag, int offset ) const;
  //! score a fragment against the nth fragment in the DB with cutoff
  double score_distance( const ChainDB& frag, int offset, double scut ) const;
  //! score a fragment against the nth fragment in the DB with residue masks
  double score_distance( const ChainDB& frag, const std::vector<Residue::TypeMask>& types, int offset, double scut ) const;
  //! return list of tentative fragment offsets matching a given fragment
  std::vector<int> match_fragment_preliminary( const ChainDB& fragdb, int nhit ) const;
  //! return list of tentative fragment offsets with residue masks
  std::vector<int> match_fragment_preliminary( const ChainDB& fragdb, const std::vector<Residue::TypeMask>& types, int nhit ) const;
  //! return the best DB fragments matching a given fragment
  std::vector<Chain> match_fragment( const ChainDB& fragdb, int nlsq, int nhit=0 ) const;
  //! return the best DB fragments matching a given fragment with residue masks
  std::vector<Chain> match_fragment( const ChainDB& fragdb, const std::vector<Residue::TypeMask>& types, int nlsq, int nhit=0 ) const;
 protected:
  std::vector<DistVec> dbdistvecs;
};


} // namespace ProteinDB


#endif
