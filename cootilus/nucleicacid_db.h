/*! NADB Top 250 main chain database */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NA_DB
#define NA_DB


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


namespace NucleicAcidDB {


/*!
  Class for storing compact nucleic acid main chain information for use
  with Top500 DB. The same class is used for storing monomers in the
  DB, and also for storing monomers for searching against the DB.
  This class stores all the information concerning a single monomer.
*/
class NucleicAcid {
 public:
  enum FLAG {NONE,COMPLETE,INCOMPLETE};
  //! Null constructor
  NucleicAcid() : flg(char(NONE)) {}
  //! Constructor: from main chain coords and monomer type
  NucleicAcid( const clipper::Coord_orth& cp, const clipper::Coord_orth& co5, const clipper::Coord_orth& cc5, const clipper::Coord_orth& cc4, const clipper::Coord_orth& co4, const clipper::Coord_orth& cc3, const clipper::Coord_orth& co3, const clipper::Coord_orth& cc2, const clipper::Coord_orth& cc1, const clipper::Coord_orth& cn, const clipper::String& type );
  //! Constructor: from MMonomer
  NucleicAcid( const clipper::MMonomer& mm );
  clipper::Coord_orth coord_p () const;  //!< P atom coordinate
  clipper::Coord_orth coord_o5() const;  //!< O5 atom coordinate
  clipper::Coord_orth coord_c5() const;  //!< C5 atom coordinate
  clipper::Coord_orth coord_c4() const;  //!< C4 atom coordinate
  clipper::Coord_orth coord_o4() const;  //!< O4 atom coordinate
  clipper::Coord_orth coord_c3() const;  //!< C3 atom coordinate
  clipper::Coord_orth coord_o3() const;  //!< O3 atom coordinate
  clipper::Coord_orth coord_c2() const;  //!< C2 atom coordinate
  clipper::Coord_orth coord_c1() const;  //!< C1 atom coordinate
  clipper::Coord_orth coord_n () const;  //!< N atom coordinate
  clipper::MMonomer mmonomer() const;    //!< Build MMonomer 
  //! Transform by rotation-translation operator
  void transform( const clipper::RTop_orth& rtop );
  //! Merge other monomer coords with this one using given weights
  bool merge( const NucleicAcid& other, const double wp, const double wo5, const double wc5, const double wc4, const double wc3, const double wo3, const double wc2, const double wc1, const double wo4, const double wn );
  char type() const { return typ; }           //!< return monomer type (1-char)
  FLAG flag() const { return FLAG(flg); }     //!< return flag (null,Ca,all)
  void set_type( char t ) { typ = t; }        //!< set monomer type (1-char)
  void set_flag( FLAG t ) { flg = char(t); }  //!< set flag (null,Ca,all)
  void data_import( const char* d );          //!< import from char array
  void data_export( char* d ) const;          //!< export to char array
  bool is_null() const { return flg == char(NONE); }  //!< test for null

  void set_flag();  //!< set flag on the basis of atoms present
 private:
  static void unpack_float( const char* d, float& f ) { const short s = ((short(d[0])<<8)&0xFF00) | ((short(d[1]))&0x00FF); f = float(s)/100.0; }
  static void pack_float( char* d, const float& f ) { const short s = rint(100.0*f); d[0] = char((s>>8)&0x00FF); d[1] = char((s)&0x00FF); }
  float p_x, p_y, p_z;  // main chain
  float o5x, o5y, o5z;
  float c5x, c5y, c5z;
  float c4x, c4y, c4z;
  float o4x, o4y, o4z;
  float c3x, c3y, c3z;
  float o3x, o3y, o3z;
  float c2x, c2y, c2z;
  float c1x, c1y, c1z;
  float n_x, n_y, n_z;  // base
  char typ, flg;
};


/*!
  Class for storing compact nucleic acid main chain information for use
  with Top500 DB. The same class is used for storing monomers in the
  DB, and also for storing monomers for searching against the DB.
  This class stores a list of monomers representing either a complete
  chain or chains (in the DB) or a complete search fragment to be
  searched against the DB.
*/
class Chain {
 public:
  //! Null constructor
  Chain() {}
  //! Add a pdb file to this DBchain (NOTE: CHAINS ARE SHIFTED TO ORIGIN)
  bool add_pdb( const clipper::String file );
  //! Add a monomer to this DBchain
  void add_monomer( const NucleicAcid& r ) { dbmonomers.push_back( r ); }
  //! Export DB to binary file
  bool save_db( const clipper::String file ) const;
  //! Import DB from binary file
  bool load_db( const clipper::String file );
  //! Import DB from binary file
  //! Merge other monomer coords with this one using given weights
  bool merge( const Chain& other, const std::vector<double>& wgt );
  //! extract fragment of given length
  Chain extract( int offset, int len ) const;
  //! check if fragment is continuous
  bool is_continuous() const;
  //! Transform by rotation-translation operator
  void transform( const clipper::RTop_orth& rtop );
  /*
  //! get RTop to fit DB fragment to given fragment
  void lsq_superpose( const Chain& frag );
  //! get RTop to fit DB fragment to given fragment with weights
  void lsq_superpose( const Chain& frag, const std::vector<double>& wgts );
  //! Get rmsd versus other fragment
  double rmsd( const Chain& other ) const;
  //! Get rmsd versus other fragment
  double rmsd( const Chain& other, const std::vector<double>& wgts ) const;  
  */
  //! Get monomer by position in list
  const NucleicAcid& operator[] ( const int& i ) const { return dbmonomers[i]; }
  //! Set monomer by position in list
  NucleicAcid& operator[] ( const int& i ) { return dbmonomers[i]; }
  //! Get number of monomers in list
  int size() const { return dbmonomers.size(); }
  // output some debug info
  void debug() const;
 protected:
  std::vector<NucleicAcid> dbmonomers;
};


} // namespace ProteinDB


#endif
