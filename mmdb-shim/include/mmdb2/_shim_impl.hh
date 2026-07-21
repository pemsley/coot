// mmdb-shim — architecture B implementation header (single unit; sub-headers
// below are thin macro-guarded wrappers around this). See MMDB_SHIM_Recon_and_Plan.md.
//
// The mmdb:: hierarchy classes ARE the stable wrapper nodes from the hardened
// core (mmdb-recon/spike -> mmdb-shim/core): each holds Manager* + parent* + a
// cached sibling index, resolves to live gemmi via g(), and exposes children as
// a canonical vector<T*> (= the identity cache AND the PPAtom/PPResidue table).
//
// Field access is via accessors (pure B): numeric -> reference-returning x()/...;
// float occ/b_iso -> value get + set_*; char[] -> pstr getter + set_*.
#pragma once

#include <gemmi/model.hpp>
#include <gemmi/resinfo.hpp>
#include <gemmi/symmetry.hpp>   // space-group / symmetry operators

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace mmdb {

// ---- basic MMDB scalar/typedefs (real MMDB: mmdb_mattype.h / mmdb_defs.h) ----
typedef double   realtype;
typedef char    *pstr;
typedef const char *cpstr;
typedef unsigned short word;
typedef char     AtomName[20];
typedef char     ResName[20];
typedef char     InsCode[10];
typedef char     ChainID[10];
typedef char     Element[10];
typedef char     AltLoc[20];
typedef char     SegID[10];
typedef char     LinkRID[20];   // Refmac link ID
typedef unsigned char byte;     // mmdb_mattype.h
typedef int     *ivector;       // mmdb_mattype.h 1-based vectors/matrices
typedef realtype *rvector;
typedef ivector *imatrix;
typedef rvector *rmatrix;
typedef char     maxMMDBName[40];

// WhatIsSet mask flags (mmdb_atom.h ASET_FLAG)
enum ASET_FLAG {
  ASET_Coordinates  = 0x00000001, ASET_Occupancy    = 0x00000002,
  ASET_tempFactor   = 0x00000004, ASET_CoordSigma   = 0x00000010,
  ASET_OccSigma     = 0x00000020, ASET_tFacSigma    = 0x00000040,
  ASET_Charge       = 0x00000080, ASET_Anis_tFac    = 0x00000100,
  ASET_Anis_tFSigma = 0x00001000, ASET_All          = 0x000FFFFF
};

// vector/matrix types (mmdb_defs.h) — plain fixed-size arrays of realtype
typedef realtype vect3[3];
typedef realtype vect4[4];
typedef vect3    mat33[3];      // realtype[3][3]
typedef vect4    mat44[4];      // realtype[4][4]
typedef mat44   *pmat44;
typedef mat44   &rmat44;

enum ERROR_CODE {
  Error_NoError     = 0,
  Error_CantOpenFile = 12,   // matches real MMDB's value
  Error_GeneralError1 = 1
};

// ---- UDData (user-defined data) — real MMDB values (mmdb_uddata.h) ----
enum UDR_TYPE { UDR_ATOM = 0, UDR_RESIDUE = 1, UDR_CHAIN = 2, UDR_MODEL = 3,
                UDR_HIERARCHY = 4 };
enum UDDATA_CODE { UDDATA_Ok = 0, UDDATA_WrongHandle = -1,
                   UDDATA_WrongUDRType = -2, UDDATA_NoData = -3 };

// ---- Selection (real MMDB values: mmdb_selmngr.h) ----
enum SELECTION_TYPE { STYPE_INVALID = -1, STYPE_UNDEFINED = 0, STYPE_ATOM = 1,
                      STYPE_RESIDUE = 2, STYPE_CHAIN = 3, STYPE_MODEL = 4 };
enum SELECTION_KEY  { SKEY_NEW = 0, SKEY_OR = 1, SKEY_AND = 2, SKEY_XOR = 3,
                      SKEY_CLR = 4, SKEY_XAND = 100 };
inline const long int MinInt4 = -2147483647;
inline const long int MaxInt4 =  2147483647;
inline const int      ANY_RES = -2147483647;   // real MMDB: extern const == MinInt4
inline const double   Pi = 3.14159265358979323846;

// PDB/CIF read flags (mmdb_io_file.h). Values are arbitrary distinct bits — the
// shim's SetFlag is a no-op, so only distinctness matters for Coot's bit ops.
enum MMDB_READ_FLAG {
  MMDBF_AutoSerials            = 0x00000001,
  MMDBF_IgnoreDuplSeqNum       = 0x00000002,
  MMDBF_IgnoreBlankLines       = 0x00000004,
  MMDBF_IgnoreRemarks          = 0x00000008,
  MMDBF_IgnoreHash             = 0x00000010,
  MMDBF_IgnoreNonCoorPDBErrors = 0x00000020,
  MMDBF_PrintCIFWarnings       = 0x00000040,
  MMDBF_All                    = 0x0000FFFF
};
enum MMDB_FCM { MMDBFCM_None = 0, MMDBFCM_All = 1, MMDBFCM_Coord = 2,
                MMDBFCM_Cryst = 4, MMDBFCM_SC = 8 };
typedef int COPY_MASK;   // Coot uses `COPY_MASK cm = MMDBFCM_All` + bit arithmetic

// Per-object UDData slots + selection membership bits. Each registered UDData
// handle maps to a (type,kind,slot); the object stores contiguous vectors
// indexed by slot. `_inSel[selHnd-1]` = is this object in selection selHnd
// (maintained by Manager::Select/SelectSphere/DeleteSelection).
struct UDStore {
  std::vector<int>         _udi;
  std::vector<double>      _udr;
  std::vector<std::string> _uds;
  std::vector<bool>        _inSel;
  bool isInSelection(int selHnd) const {
    return selHnd >= 1 && selHnd <= (int)_inSel.size() && _inSel[selHnd - 1];
  }
  void _setInSel(int selHnd, bool v) {
    if ((int)_inSel.size() < selHnd) _inSel.resize(selHnd, false);
    _inSel[selHnd - 1] = v;
  }
};

class Atom; class Residue; class Chain; class Model; class Manager;
typedef Atom    *PAtom;    typedef Atom    **PPAtom;
typedef Residue *PResidue; typedef Residue **PPResidue;
typedef Chain   *PChain;   typedef Chain   **PPChain;
typedef Model   *PModel;   typedef Model   **PPModel;
typedef Manager *PManager; typedef Manager **PPManager;

struct Contact { int id1, id2; long group; realtype dist; };
typedef Contact  *PContact;

// base for records held in MMDB containers (Title compound/author, LINK, …)
class ContainerClass { public: virtual ~ContainerClass() {} };
typedef ContainerClass *PContainerClass;

// LINK record. Public data members mirror real MMDB (Coot reads them directly).
// Not gemmi-backed yet — Model::GetNumberOfLinks currently returns 0 (TODO: map
// gemmi Structure connections), so these are declared for compilation.
class Link : public ContainerClass {
public:
  AtomName atName1{}, atName2{};
  AltLoc   aloc1{}, aloc2{};
  ResName  resName1{}, resName2{};
  ChainID  chainID1{}, chainID2{};
  InsCode  insCode1{}, insCode2{};
  int      seqNum1 = 0, seqNum2 = 0;
  int      s1 = 1, i1 = 0, j1 = 0, k1 = 0;   // symmetry id of 1st atom
  int      s2 = 1, i2 = 0, j2 = 0, k2 = 0;   // symmetry id of 2nd atom
  realtype dist = 0;
  void Copy(PContainerClass o) { if (auto *l = dynamic_cast<Link *>(o)) *this = *l; }
};
typedef Link *PLink; typedef Link **PPLink;

// Refmac LINK record (mmdb_model.h LinkR). Public members mirror real MMDB;
// Model::GetNumberOfLinkRs returns 0 for now (TODO: map gemmi connections).
class LinkR {
public:
  LinkRID  linkRID{};
  AtomName atName1{}, atName2{};
  AltLoc   aloc1{}, aloc2{};
  ResName  resName1{}, resName2{};
  ChainID  chainID1{}, chainID2{};
  int      seqNum1 = 0, seqNum2 = 0;
  InsCode  insCode1{}, insCode2{};
  realtype dist = 0;
};
typedef LinkR *PLinkR; typedef LinkR **PPLinkR;

// CIS-peptide record (mmdb_model.h CisPep). Public members mirror real MMDB;
// Model::GetNumberOfCisPeps returns 0 for now (TODO: map gemmi cispeps).
class CisPep {
public:
  int      serNum = 0;
  ResName  pep1{};
  ChainID  chainID1{};
  int      seqNum1 = 0;
  InsCode  icode1{};
  ResName  pep2{};
  ChainID  chainID2{};
  int      seqNum2 = 0;
  InsCode  icode2{};
  int      modNum = 0;
  realtype measure = 0;
};
typedef CisPep *PCisPep;

// Container of LINK records (mmdb_model.h LinkContainer). Minimal: Coot only
// declares `empty_links_container()` returning one by value; never dereferenced.
class LinkContainer {
public:
  std::vector<ContainerClass *> data;
  int Length() { return (int)data.size(); }
  PContainerClass GetContainerClass(int i) { return (i >= 0 && i < (int)data.size()) ? data[i] : nullptr; }
};
typedef LinkContainer *PLinkContainer;

// PDB title records (mmdb_title.h). Coot subclasses Manager & Title to reach the
// COMPND/AUTHOR line containers. Not gemmi-backed yet (title/header records are
// dropped on round-trip) — just enough surface to compile & run. TODO: map to
// gemmi Structure meta (raw_remarks / metadata).
class Compound : public ContainerClass { public: char Line[256] = {0}; };
typedef Compound *PCompound;
class Author   : public ContainerClass { public: char Line[256] = {0}; };
typedef Author   *PAuthor;
class Journal  : public ContainerClass { public: char Line[256] = {0}; };
typedef Journal  *PJournal;
class TitleContainer {
public:
  std::vector<ContainerClass *> data;
  int Length() { return (int)data.size(); }
  PContainerClass GetContainerClass(int i) {
    return (i >= 0 && i < (int)data.size()) ? data[i] : nullptr;
  }
};
class Title {
public:
  TitleContainer compound, author, journal;   // public so Coot's access_title can reach them
  TitleContainer *GetCompound() { return &compound; }   // real Title exposes these
  TitleContainer *GetAuthor()   { return &author; }     // publicly; access_title
  TitleContainer *GetJournal()  { return &journal; }    // inherits GetJournal()
};

// gzip mode flag (mmdb_io_file.h). Minimal mmdb::io — the shim does I/O via gemmi,
// so only this compression-mode enum is provided (Coot passes it to write calls).
namespace io { enum GZ_MODE { GZM_NONE = 0, GZM_CHECK = 1, GZM_ENFORCE = 2 }; }

// Crystal/symmetry record (mmdb_cryst.h). Minimal — used by Coot as a pointer
// type; symmetry math goes through Manager::GetTMatrix (gemmi TODO).
class Cryst { public:
  virtual ~Cryst() {}
  // symmetry not carried on the bare Cryst (Manager owns gemmi cell/SG) — identity/0.
  int  GetTMatrix(mat44 &T, int Nop, int a, int b, int c) {
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) T[i][j]=(i==j)?1.0:0.0;
    return (Nop==0 && a==0 && b==0 && c==0) ? 0 : 1;
  }
  int  GetNumberOfSymOps() { return 0; }
  pstr GetSymOp(int) { return nullptr; }
};
typedef Cryst *PCryst;

// initialise a 4x4 matrix to identity (mmdb_mattype.h Mat4Init)
inline void Mat4Init(mat44 &A) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) A[i][j] = (i == j) ? 1.0 : 0.0;
}

// mmdb::math graph-matching subsystem — full classes defined in _graph_impl.hh
// (included at end of this file, after Atom/Residue are complete). Only the
// Alignment class (unused by the cootapi build) stays a forward decl.
namespace math { class Alignment; }

struct AtomBond { PAtom atom = nullptr; int order = 0; };
typedef AtomBond *PAtomBond; typedef AtomBond **PPAtomBond;

struct AtomStat {   // selection coordinate statistics (mmdb_atom.h)
  int nAtoms = 0;
  realtype xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
  realtype xm = 0, ym = 0, zm = 0;   // coordinate means (centroid)
  realtype GetMaxSize() {
    realtype dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
    return dx > dy ? (dx > dz ? dx : dz) : (dy > dz ? dy : dz);
  }
};
typedef AtomStat &RAtomStat;

// secondary-structure element codes (mmdb_tables.h)
enum SSE_CODE { SSE_None = 0, SSE_Strand = 1, SSE_Bulge = 2, SSE_3Turn = 3,
                SSE_4Turn = 4, SSE_5Turn = 5, SSE_Helix = 6 };

// PDBCleanup flags (mmdb_root.h) — bit flags OR'd into PDBCleanup(word)
// misc return-code / sort-key enums (mmdb_cryst.h / mmdb_selmngr.h / mmdb_tables.h)
enum { SYMOP_Ok = 0, SYMOP_NoLibFile = -1, SYMOP_UnknownSpaceGroup = -2 };
enum { SSERC_Ok = 0, SSERC_noResidues = 1 };
enum { SORT_CHAIN_ChainID_Asc = 0, SORT_CHAIN_ChainID_Desc = 1 };
enum { CNSORT_OFF = 0, CNSORT_1INC = 1, CNSORT_1DEC = 2, CNSORT_2INC = 3, CNSORT_2DEC = 4 };

enum PDB_CLEAN_FLAG {
  PDBCLEAN_ATNAME       = 0x00000001,
  PDBCLEAN_TER          = 0x00000002,
  PDBCLEAN_CHAIN        = 0x00000004,
  PDBCLEAN_CHAIN_STRONG = 0x00000008,
  PDBCLEAN_ALTCODE      = 0x00000010,
  PDBCLEAN_ALTCODE_STRONG = 0x00000020,
  PDBCLEAN_SERIAL       = 0x00000040,
  PDBCLEAN_SEQNUM       = 0x00000080,
  PDBCLEAN_INDEX        = 0x00000800,
  PDBCLEAN_ELEMENT      = 0x00001000,
  PDBCLEAN_ELEMENT_STRONG = 0x00002000
};

// SS records — minimal public-member structs. Model::GetNumberOf{Helices,Sheets}
// return 0 for now (TODO: map gemmi Structure helices/sheets), so these aren't
// dereferenced; fields present for compilation.
class Helix  { public:
  ChainID initChainID{}, endChainID{}; int initSeqNum = 0, endSeqNum = 0, serNum = 0, helixClass = 0, length = 0;
  ResName initResName{}, endResName{}; InsCode initICode{}, endICode{}; char helixID[20]{}, comment[80]{};
};
class Strand { public:
  ChainID initChainID{}, endChainID{}; int initSeqNum = 0, endSeqNum = 0, strandNo = 0, sense = 0;
  ResName initResName{}, endResName{}; InsCode initICode{}, endICode{}; char sheetID[20]{};
};
class Sheet  { public: int nStrands = 0; Strand **strand = nullptr; char sheetID[20]{}; };
class Sheets { public: int nSheets = 0; Sheet **sheet = nullptr; };  // container (SS TODO)
typedef Helix *PHelix; typedef Strand *PStrand; typedef Sheet *PSheet; typedef Sheets *PSheets;
// container of helices (Model.helices); Coot's access_model subclass fills it.
class Helices { public: std::vector<Helix *> data; void AddData(PHelix h) { if (h) data.push_back(h); } int nHelices = 0; };

// container of symmetry operators (mmdb_symop.h SymOps). Coot fills it from a
// space group; ops are xyz-triplet strings.
class SymOps {
  std::vector<std::string> ops;
  std::deque<std::string>  buf;
public:
  int  AddSymOp(cpstr xyz) { ops.push_back(xyz ? xyz : ""); return 0; }
  int  GetNofSymOps()      { return (int)ops.size(); }
  pstr GetSymOp(int n) {
    if (n < 0 || n >= (int)ops.size()) return nullptr;
    buf.push_back(ops[n]); return (pstr) buf.back().c_str();
  }
  void FreeMemory() { ops.clear(); }
};

[[noreturn]] inline void unimpl(const char *w) {
  throw std::logic_error(std::string("mmdb-shim: unimplemented: ") + w);
}

// ---- free functions (mmdb_tables.h / mmdb_mattype.h) ----
inline void InitMatType() {}   // real MMDB inits static matrix-type tables; no-op here
inline cpstr GetErrorDescription(ERROR_CODE ec) {
  switch (ec) {
    case Error_NoError:      return "no error";
    case Error_CantOpenFile: return "cannot open file";
    default:                 return "MMDB error";
  }
}
inline realtype getVdWaalsRadius(cpstr element) {
  return gemmi::Element(element ? element : "X").vdw_r();
}

// UDData helpers (defined after Manager); each class forwards with its UDR type.
int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, int v);
int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, realtype v);
int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, cpstr v);
int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, int &v);
int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, realtype &v);
int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, pstr &v);

// ===========================================================================
class Atom : public UDStore {
public:
  Manager *mgr = nullptr;
  Residue *res = nullptr;   // parent; null => detached (use _local)
  int ai = 0;               // cached index within parent residue's atoms
  bool alive = true;
  gemmi::Atom _local;       // backing store while detached (see g() resolvers)
  int Het = 0;              // heteroatom flag (MMDB public field; Coot sets it)
  int Ter = 0;              // chain-terminator flag (gemmi has none -> always 0)
  word WhatIsSet = 0;       // ASET_* mask; ASET_Anis_tFac set on load if aniso present
  AtomName label_atom_id{}; // mmcif label_atom_id (shim-owned; Coot sets on build)

  Atom() = default;
  explicit Atom(Residue *r);  // construct + add to residue (out-of-line)

  gemmi::Atom &g() const;   // resolve to live gemmi (defined after Manager)

  // --- rewritten field accessors (pure B) ---
  // Scalar fields -> reference-returning accessors, so a uniform `->field`->
  // `->field()` rewrite covers both reads and writes. (occ/b_iso/charge are
  // narrower than realtype in gemmi, so those refs are float/schar-typed — the
  // rare take-address-of-realtype sites surface at Coot build time.)
  // non-const (writable ref) + const (by value) overloads, so reads work on a
  // `const mmdb::Atom` and writes work through `->x() = v` on a non-const one.
  realtype &x() { return g().pos.x; }   realtype x() const { return g().pos.x; }
  realtype &y() { return g().pos.y; }   realtype y() const { return g().pos.y; }
  realtype &z() { return g().pos.z; }   realtype z() const { return g().pos.z; }
  float &occupancy() { return g().occ; }   float occupancy() const { return g().occ; }
  float &tempFactor() { return g().b_iso; } float tempFactor() const { return g().b_iso; }
  signed char &charge() { return g().charge; } signed char charge() const { return g().charge; }
  int &serNum() { return g().serial; }  int serNum() const { return g().serial; }
  // altLoc is a char[] (C-string) in MMDB; gemmi stores a single char. Return a
  // buffer-backed C-string ("" when unset) so strcmp/strcpy-style code works.
  // The non-const overload returns a WRITABLE buffer so `strncpy(at->altLoc(),..)`
  // compiles; the buffer's first char is pushed back into gemmi by Residue::AddAtom
  // (the buffer is refreshed from gemmi on entry, so reads stay correct).
  pstr        altLoc()       { _altloc_buf[0] = g().altloc; _altloc_buf[1] = '\0'; return _altloc_buf; }
  const char *altLoc() const { _altloc_buf[0] = g().altloc; _altloc_buf[1] = '\0'; return _altloc_buf; }
  void set_occupancy(realtype v) { g().occ = (float)v; }
  void set_tempFactor(realtype v) { g().b_iso = (float)v; }
  void set_altLoc(char c) { g().altloc = c; }
  void SetCharge(realtype ch) { g().charge = (signed char) ch; }
  // coordinate/occupancy/B ESDs (MMDB public fields) — gemmi has none, so shim-
  // owned; reference-returning so the rewritten `->sigX` covers reads and writes.
  float &sigX()    { return _sigx; }   float &sigY()    { return _sigy; }
  float &sigZ()    { return _sigz; }   float &sigOcc()  { return _sigocc; }
  float &sigTemp() { return _sigtemp; }
  bool isMetal() const { return gemmi::Element(g().element).is_metal(); }
  // anisotropic B tensor — gemmi's SMat33<float> aniso. Reference-returning so the
  // rewritten `->u11` covers both reads (bonds display) and writes (SHELX import).
  // NOTE: writing here does not set ASET_Anis_tFac in WhatIsSet (TODO if needed).
  float &u11() { return g().aniso.u11; }
  float &u22() { return g().aniso.u22; }
  float &u33() { return g().aniso.u33; }
  float &u12() { return g().aniso.u12; }
  float &u13() { return g().aniso.u13; }
  float &u23() { return g().aniso.u23; }
  // bonds — not modelled yet (gemmi connections); report none.
  int  GetNBonds() { return 0; }
  void GetBonds(PAtomBond &atomBond, int &n) { atomBond = nullptr; n = 0; }
  int  AddBond(PAtom /*a*/, int /*order*/, int /*nAdd*/ = 1) { return 0; }
  SegID segID{};   // shim-owned (gemmi has no segID); MMDB public char[] field

  // --- method surface (hot subset; rest stubbed) ---
  pstr GetAtomName() const;           // aligned name, MMDB semantics (const: called on const Atom)
  void SetAtomName(const AtomName aName);
  pstr GetElementName();
  void SetElementName(const Element elName);
  pstr  GetChainID();
  int   GetSeqNum();
  pstr  GetInsCode();
  pstr  GetResName();
  Residue *&GetResidue() { return res; }  // ref: rewritten `->residue` is assignable
  void  SetResidue(Residue *r) { res = r; }
  Chain   *GetChain();          // out-of-line (needs complete Residue/Chain)
  Model   *GetModel();          // out-of-line
  int   GetModelNum();
  // residue-delegating accessors (bound by the Python API); out-of-line.
  pstr GetLabelCompID();  pstr GetLabelAsymID();
  int  GetLabelSeqID();   int  GetLabelEntityID();
  int  GetResidueNo();    int  GetSSEType();
  bool isSolvent();       bool isNTerminus();   bool isCTerminus();
  bool  isTer() const { return false; } // gemmi has no TER atoms; see notes
  void  SetCoordinates(realtype xx, realtype yy, realtype zz,
                       realtype occ, realtype tF);
  int   GetIndex();
  void  MakeTer() { Ter = 1; }        // mark as chain terminator
  pstr  GetAtomID(pstr S);            // "/mdl/chain/seq(res).ins/name[elem]:alt" (out-of-line)
  int GetUDData(int h, pstr &v)     { return ud_get(mgr, UDR_ATOM, *this, h, v); }
  // copy another atom's data into this one (mmdb Atom::Copy — no hierarchy refs)
  void  Copy(PAtom a) {
    g() = a->g();
    Het = a->Het; WhatIsSet = a->WhatIsSet;
    std::memcpy(segID, a->segID, sizeof segID);
  }
  // apply a 4x4 (rot+trans) or 3x3+vec to the coordinates (mmdb Atom::Transform)
  void  Transform(const mat44 &tm) {
    gemmi::Position &p = g().pos; double x = p.x, y = p.y, z = p.z;
    p.x = tm[0][0]*x + tm[0][1]*y + tm[0][2]*z + tm[0][3];
    p.y = tm[1][0]*x + tm[1][1]*y + tm[1][2]*z + tm[1][3];
    p.z = tm[2][0]*x + tm[2][1]*y + tm[2][2]*z + tm[2][3];
  }
  void  Transform(const mat33 &tm, vect3 &v) {
    gemmi::Position &p = g().pos; double x = p.x, y = p.y, z = p.z;
    p.x = tm[0][0]*x + tm[0][1]*y + tm[0][2]*z + v[0];
    p.y = tm[1][0]*x + tm[1][1]*y + tm[1][2]*z + v[1];
    p.z = tm[2][0]*x + tm[2][1]*y + tm[2][2]*z + v[2];
  }
  // UDData
  int PutUDData(int h, int v)      { return ud_put(mgr, UDR_ATOM, *this, h, v); }
  int PutUDData(int h, realtype v) { return ud_put(mgr, UDR_ATOM, *this, h, v); }
  int PutUDData(int h, cpstr v)    { return ud_put(mgr, UDR_ATOM, *this, h, v); }
  int GetUDData(int h, int &v)      { return ud_get(mgr, UDR_ATOM, *this, h, v); }
  int GetUDData(int h, realtype &v) { return ud_get(mgr, UDR_ATOM, *this, h, v); }

private:
  friend class Residue;   // AddAtom pushes the strncpy'd altLoc buffer to gemmi
  mutable AtomName _name_buf{}; Element _elem_buf{}; mutable char _altloc_buf[4]{};
  float _sigx = 0, _sigy = 0, _sigz = 0, _sigocc = 0, _sigtemp = 0;
};

// ===========================================================================
class Residue : public UDStore {
public:
  Manager *mgr = nullptr;
  Chain *chain = nullptr;   // parent; null => detached (use _local)
  int ri = 0;
  bool alive = true;
  gemmi::Residue _local;    // backing store while detached
  std::vector<Atom *> atoms;  // canonical child wrappers == PPAtom table
  PPAtom atom = nullptr;    // MMDB public atom-table field; kept = atoms.data()
  int    nAtoms = 0;        // MMDB public field; kept = atoms.size()
  void _sync_atom() { atom = atoms.data(); nAtoms = (int)atoms.size(); }
  // mmcif label_* (shim-owned; Coot sets when building dictionary residues)
  ResName label_comp_id{}; ChainID label_asym_id{}; int label_seq_id = 0, label_entity_id = 0;
  pstr GetLabelCompID()   { return label_comp_id; }
  pstr GetLabelAsymID()   { return label_asym_id; }
  int  GetLabelSeqID()    { return label_seq_id; }
  int  GetLabelEntityID() { return label_entity_id; }
  int  GetResidueNo()     { return ri; }         // 0-based index within its chain
  int  GetNofAltLocations() {                    // distinct non-blank altLocs
    std::set<char> a; for (Atom *at : atoms) { char c = at->g().altloc; if (c && c != ' ') a.insert(c); }
    return a.empty() ? 1 : (int) a.size();
  }
  bool isSugar()  { return false; }              // TODO: gemmi residue classification
  bool isModRes() { return false; }              // TODO

  Residue() = default;
  explicit Residue(Chain *c);   // construct + add to chain (out-of-line)

  // MMDB public char-array fields. Coot reads `residue->name` and writes
  // `strncpy(residue->insCode,..)`. Kept as the interface: synced gemmi->buffer on
  // load (_load_id, in build_from_gemmi) and buffer->gemmi at the adopt point
  // (_store_id, in Chain::Add/InsResidue). SetResName/SetResID keep both in step.
  ResName name{};
  InsCode insCode{};
  void _load_id() {
    std::snprintf(name, sizeof name, "%s", g().name.c_str());
    insCode[0] = g().seqid.icode && g().seqid.icode != ' ' ? g().seqid.icode : '\0';
    insCode[1] = '\0';
  }
  void _store_id() {
    g().name = name;
    g().seqid.icode = insCode[0] ? insCode[0] : ' ';
  }

  gemmi::Residue &g() const;

  int   GetNumberOfAtoms() { return (int)atoms.size(); }
  int   GetNumberOfAtoms(bool /*countTers*/) { return (int)atoms.size(); }
  PAtom GetAtom(int atomNo) {
    return (atomNo >= 0 && atomNo < (int)atoms.size()) ? atoms[atomNo] : nullptr;
  }
  PAtom GetAtom(const AtomName aname, const Element elname = nullptr,
                const AltLoc aloc = nullptr);
  void  GetAtomTable(PPAtom &atomTable, int &n) { atomTable = atoms.data(); n = (int)atoms.size(); }
  PAtom AddAtom(Manager &m, gemmi::Atom a);   // append: O(1)
  // Adopt a detached atom (Coot's `new mmdb::Atom` idiom). Copies the atom's
  // local gemmi into this residue's gemmi (detached or bound, via g()) and
  // rebinds the wrapper. Pushes the strncpy'd altLoc buffer back into gemmi.
  int   AddAtom(PAtom atm) {
    g().atoms.push_back(atm->_local);
    atm->res = this; atm->mgr = mgr; atm->ai = (int)atoms.size();
    if (atm->_altloc_buf[0]) g().atoms[atm->ai].altloc = atm->_altloc_buf[0];
    atoms.push_back(atm);
    _sync_atom();
    return 0;
  }
  void  DeleteAtom(int pos);
  void  TrimAtomTable() {}   // compact after deletions — shim keeps them in sync

  pstr  GetResName();
  void  SetResName(const ResName n) {
    g().name = n ? n : "";
    std::snprintf(name, sizeof name, "%s", n ? n : "");
  }
  void  SetResID(const ResName resName, int seqNo, const InsCode ic) {
    g().name = resName ? resName : "";
    g().seqid.num.value = seqNo;
    g().seqid.icode = (ic && ic[0]) ? ic[0] : ' ';
    std::snprintf(name, sizeof name, "%s", resName ? resName : "");
    insCode[0] = (ic && ic[0]) ? ic[0] : '\0'; insCode[1] = '\0';
  }
  int  &GetSeqNum();          // writable (rewrite maps `->seqNum` reads and writes)
  pstr  GetInsCode();
  pstr  GetChainID();
  int   GetModelNum();
  int  &GetIndex() { return ri; }     // ref: rewritten `->index` is assignable
  Chain   *GetChain() { return chain; }
  Model   *GetModel();                // out-of-line (Chain incomplete here)
  // terminus tests — positional within the chain (approximates MMDB's peptide-bond
  // check; good enough for Coot's terminal-residue handling). TODO: bond-aware.
  bool  isNTerminus() { return chain && ri == 0; }
  bool  isCTerminus();               // last in chain (out-of-line: needs Chain)
  pstr  GetResidueID(pstr S) {       // "seqnum(name):inscode"
    if (S) std::snprintf(S, 100, "%d(%s):%s", GetSeqNum(), name, insCode);
    return S;
  }
  Residue *next = nullptr; // MMDB has this; wired lazily if needed
  int   SSE = SSE_None;    // secondary-structure element (shim-owned public field)
  bool  isAminoacid()  { return gemmi::find_tabulated_residue(g().name).is_amino_acid(); }
  bool  isNucleotide() { return gemmi::find_tabulated_residue(g().name).is_nucleic_acid(); }
  bool  isDNARNA()     { return isNucleotide(); }
  bool  isSolvent()    { return gemmi::find_tabulated_residue(g().name).is_water(); }
  // UDData
  int PutUDData(int h, int v)      { return ud_put(mgr, UDR_RESIDUE, *this, h, v); }
  int PutUDData(int h, realtype v) { return ud_put(mgr, UDR_RESIDUE, *this, h, v); }
  int PutUDData(int h, cpstr v)    { return ud_put(mgr, UDR_RESIDUE, *this, h, v); }
  int GetUDData(int h, int &v)      { return ud_get(mgr, UDR_RESIDUE, *this, h, v); }
  int GetUDData(int h, realtype &v) { return ud_get(mgr, UDR_RESIDUE, *this, h, v); }

private:
  ResName _resname_buf{}; InsCode _inscode_buf{};
};

// ===========================================================================
class Chain : public UDStore {
public:
  Manager *mgr = nullptr;
  Model *model = nullptr;   // parent; null => detached (use _local)
  int ci = 0;
  bool alive = true;
  gemmi::Chain _local;      // backing store while detached
  std::vector<Residue *> residues;

  gemmi::Chain &g() const;

  int      GetNumberOfResidues() { return (int)residues.size(); }
  PResidue GetResidue(int resNo) {
    return (resNo >= 0 && resNo < (int)residues.size()) ? residues[resNo] : nullptr;
  }
  // find by (seqNum, insCode) — MMDB's 2-arg overload
  PResidue GetResidue(int seqNum, const InsCode insCode) {
    char ic = (insCode && insCode[0]) ? insCode[0] : ' ';
    for (Residue *r : residues) {
      gemmi::Residue &gr = r->g();
      char ric = gr.seqid.icode ? gr.seqid.icode : ' ';
      if (gr.seqid.num.value == seqNum && ric == ic) return r;
    }
    return nullptr;
  }
  void     GetResidueTable(PPResidue &t, int &n) { t = residues.data(); n = (int)residues.size(); }
  // delete residue at index: erase gemmi + wrapper, reindex the tail
  void     DeleteResidue(int resNo) {
    if (resNo < 0 || resNo >= (int)residues.size()) return;
    g().residues.erase(g().residues.begin() + resNo);
    residues.erase(residues.begin() + resNo);
    for (int k = resNo; k < (int)residues.size(); ++k) residues[k]->ri = k;
  }
  void     TrimResidueTable() {}   // compact after deletions — shim stays in sync
  void     DeleteResidue(int seqNum, const InsCode ic) {   // by (seqNum, insCode)
    PResidue r = GetResidue(seqNum, ic);
    if (r) DeleteResidue(r->ri);
  }
  pstr     GetChainID();
  pstr     GetChainID(pstr buf) { if (buf) std::snprintf(buf, sizeof(ChainID), "%s", g().name.c_str()); return buf; }
  Manager *GetCoordHierarchy() { return mgr; }   // parent manager
  void     SetChainID(const ChainID id) { g().name = id ? id : ""; }
  Chain() = default;
  Chain(Model *m, const ChainID id);   // construct + add to model (out-of-line)
  void     Copy(PChain src);           // deep-copy subtree (out-of-line: needs Manager)
  void     SortResidues(int /*sortKey*/ = 0) {}  // gemmi keeps file order; no-op
  bool isAminoacidChain();   // defined out-of-line (needs Residue predicates)
  bool isNucleotideChain();
  bool isSolventChain();
  PResidue AddResidue(Manager &m, gemmi::Residue r);       // append
  PResidue InsResidue(Manager &m, int pos, gemmi::Residue r);
  // Adopt a detached residue (its atom wrappers already point at it, so they
  // ride along once its gemmi is copied in and the wrapper is rebound).
  int AddResidue(PResidue res) {
    res->_store_id();                   // push name/insCode buffers into gemmi
    g().residues.push_back(res->g());   // res detached -> its _local (with atoms)
    res->chain = this; res->mgr = mgr; res->ri = (int)residues.size();
    residues.push_back(res);
    return 0;
  }
  int InsResidue(PResidue res, int pos) {
    if (pos < 0) pos = 0;
    if (pos > (int)residues.size()) pos = (int)residues.size();
    res->_store_id();
    g().residues.insert(g().residues.begin() + pos, res->g());
    res->chain = this; res->mgr = mgr; res->ri = pos;
    residues.insert(residues.begin() + pos, res);
    for (int k = pos + 1; k < (int)residues.size(); ++k) residues[k]->ri = k;
    return 0;
  }

private:
  ChainID _chainid_buf{};
};

// ===========================================================================
class Model : public UDStore {
public:
  Manager *mgr = nullptr;   // null => detached (use _local)
  int mi = 0;               // 0-based internal; GetModel is 1-based externally
  gemmi::Model _local{1};   // backing store while detached (gemmi Model num is int)
  std::vector<Chain *> chains;

  gemmi::Model &g() const;

  int    GetNumberOfChains() { return (int)chains.size(); }
  PChain GetChain(int chainNo) {
    return (chainNo >= 0 && chainNo < (int)chains.size()) ? chains[chainNo] : nullptr;
  }
  PChain GetChain(const ChainID chID);
  // Adopt a detached chain (Coot's `new mmdb::Chain` idiom): copy its local
  // gemmi (with any residues/atoms) into this model and rebind, cascading mgr
  // to the sub-tree that was built while detached (mgr was null).
  int    AddChain(PChain chn) {
    g().chains.push_back(chn->g());
    chn->model = this; chn->mgr = mgr; chn->ci = (int)chains.size();
    chains.push_back(chn);
    for (Residue *r : chn->residues) {
      r->mgr = mgr;
      for (Atom *a : r->atoms) a->mgr = mgr;
    }
    return 0;
  }
  int    GetSerNum() { return mi + 1; }
  // delete chain at index: erase gemmi + wrapper, reindex the tail
  void   DeleteChain(int chainNo) {
    if (chainNo < 0 || chainNo >= (int)chains.size()) return;
    g().chains.erase(g().chains.begin() + chainNo);
    chains.erase(chains.begin() + chainNo);
    for (int k = chainNo; k < (int)chains.size(); ++k) chains[k]->ci = k;
  }
  void   DeleteChain(const ChainID chainID) {
    for (int i = 0; i < (int)chains.size(); ++i)
      if (chains[i]->g().name == (chainID ? chainID : "")) { DeleteChain(i); return; }
  }
  void   GetChainTable(PPChain &t, int &n) { t = chains.data(); n = (int)chains.size(); }
  std::vector<Atom *> all_atoms;   // flat, filled by build_from_gemmi
  PPAtom GetAllAtoms() { return all_atoms.data(); }
  int    GetNumberOfAtoms() { return (int)all_atoms.size(); }
  int    GetNumberOfAtoms(bool /*countTers*/) { return (int)all_atoms.size(); }
  int    CalcSecStructure(bool /*flag*/) { return 0; }  // TODO: gemmi SS assignment
  // LINK records — Coot-owned Link* objects stored here (AddLink); GetLink is
  // 1-based like MMDB. (Reading from gemmi connections is a separate TODO.)
  std::vector<Link *> _links;
  int    GetNumberOfLinks() { return (int)_links.size(); }
  PLink  GetLink(int i) { return (i >= 1 && i <= (int)_links.size()) ? _links[i - 1] : nullptr; }
  void   AddLink(PLink link) { if (link) _links.push_back(link); }
  int    GetNumberOfLinkRs() { return 0; }
  PLinkR GetLinkR(int /*i*/) { return nullptr; }
  std::vector<CisPep *> _cispeps;
  int    GetNumberOfCisPeps() { return (int)_cispeps.size(); }
  PCisPep GetCisPep(int i) { return (i >= 1 && i <= (int)_cispeps.size()) ? _cispeps[i - 1] : nullptr; }
  void   AddCisPep(PCisPep cp) { if (cp) _cispeps.push_back(cp); }
  void   RemoveCisPeps() { _cispeps.clear(); }
  // secondary structure — TODO: map from gemmi helices/sheets.
  int     GetNumberOfHelices() { return 0; }
  PHelix  GetHelix(int /*i*/) { return nullptr; }
  int     GetNumberOfSheets() { return 0; }
  PSheet  GetSheet(int /*i*/) { return nullptr; }
  Sheets  sheets;                          // SS records (not gemmi-backed; access_model fills)
  Helices helices;                         // "     "     "
  PSheets GetSheets() { return &sheets; }
  int     GetModelID() { return mi + 1; }
  pstr    GetModelID(pstr buf) { if (buf) std::snprintf(buf, 16, "%d", mi + 1); return buf; }
  int     CalcSecStructure(int /*flag*/, int /*selHnd*/) { return SSERC_Ok; }  // TODO gemmi SS
  void    Copy(PModel src);                       // deep-copy subtree (out-of-line)
  Manager *GetCoordHierarchy() { return mgr; }   // parent manager
  int     GetNumberOfResidues() {
    int n = 0; for (Chain *c : chains) n += c->GetNumberOfResidues(); return n;
  }
  LinkContainer _linkc;
  PLinkContainer GetLinks() {
    _linkc.data.assign(_links.begin(), _links.end()); return &_linkc;
  }
  void    RemoveLinks() { _links.clear(); }
  void    SortChains(int /*sortKey*/ = 0) {}   // gemmi keeps file order; no-op
  PChain  CreateChain(const ChainID id);       // add empty chain (out-of-line: needs Manager)
  int     GetNumberOfStrands(int /*sheetNo*/) { return 0; }
  PStrand GetStrand(int /*sheetNo*/, int /*strandNo*/) { return nullptr; }
};

// ===========================================================================
class Manager {
public:
  gemmi::Structure st;
  // stable-address pools
  std::deque<Atom> atom_pool;
  std::deque<Residue> res_pool;
  std::deque<Chain> chain_pool;
  std::deque<Model> model_pool;
  std::vector<Model *> models;

  Atom    *newAtom() { atom_pool.emplace_back(); return &atom_pool.back(); }
  Residue *newRes()  { res_pool.emplace_back();  return &res_pool.back();  }
  Chain   *newChain(){ chain_pool.emplace_back();return &chain_pool.back();}
  Model   *newModel(){ model_pool.emplace_back();return &model_pool.back();}

  int    GetNumberOfModels() { return (int)models.size(); }
  PModel GetModel(int modelNo) {   // MMDB: 1 <= modelNo <= nModels
    int i = modelNo - 1;
    return (i >= 0 && i < (int)models.size()) ? models[i] : nullptr;
  }
  // per-model chain access (modelNo is 1-based, chainNo 0-based) — mmdb_coormngr.h
  int    GetNumberOfChains(int modelNo) {
    PModel m = GetModel(modelNo); return m ? m->GetNumberOfChains() : 0;
  }
  PChain GetChain(int modelNo, int chainNo) {
    PModel m = GetModel(modelNo); return m ? m->GetChain(chainNo) : nullptr;
  }
  // Re-index/renumber after edits. The shim keeps sibling indices in sync as it
  // mutates, so this is a no-op re-validation for now (TODO: serial renumbering).
  word   PDBCleanup(word /*CleanKey*/) { return 0; }

  // PDB title records — Coot reaches `title` via an access_mol subclass.
  Title  title;
  pstr   GetStructureTitle(pstr T) { if (T) T[0] = '\0'; return T; }

  // symmetry transformation matrix. Real symmetry needs gemmi spacegroup/cell;
  // for now return identity for the no-op (Nop==0, no cell shift) and signal
  // "no symmetry" (nonzero) otherwise so Coot skips symmetry expansion. TODO.
  int    GetTMatrix(mat44 &TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c) {
    Mat4Init(TMatrix);
    return (Nop == 0 && cellshift_a == 0 && cellshift_b == 0 && cellshift_c == 0) ? 0 : 1;
  }

  void build_from_gemmi();

  // adopt a detached model (Coot: `new mmdb::Model` -> AddChain… -> AddModel).
  // Copy its local gemmi into st, rebind, cascade mgr through the sub-tree.
  int AddModel(PModel mw) {
    st.models.push_back(mw->g());
    mw->mgr = this; mw->mi = (int)models.size();
    models.push_back(mw);
    for (Chain *cw : mw->chains) {
      cw->mgr = this;
      for (Residue *rw : cw->residues) {
        rw->mgr = this;
        for (Atom *aw : rw->atoms) { aw->mgr = this; all_atoms.push_back(aw); mw->all_atoms.push_back(aw); }
      }
    }
    return 0;
  }

  // clone another manager's structure (mmdb Manager::Copy(PManager, COPY_MASK)).
  // Copies the whole gemmi Structure and rebuilds all wrappers — clean & correct.
  void Copy(PManager m, int /*CopyMask*/) { if (m) { st = m->st; build_from_gemmi(); } }

  // ---- crystal cell & symmetry (gemmi UnitCell / SpaceGroup) ----
  std::string _sg_buf, _symop_buf;
  void GetCell(realtype &a, realtype &b, realtype &c, realtype &al, realtype &be,
               realtype &ga, realtype &vol, int &orthcode) {
    const gemmi::UnitCell &u = st.cell;
    a = u.a; b = u.b; c = u.c; al = u.alpha; be = u.beta; ga = u.gamma;
    vol = u.volume; orthcode = 1;
  }
  void GetCell(realtype &a, realtype &b, realtype &c, realtype &al, realtype &be,
               realtype &ga, realtype &vol) { int oc; GetCell(a,b,c,al,be,ga,vol,oc); }
  void SetCell(realtype a, realtype b, realtype c, realtype al, realtype be,
               realtype ga, int /*OrthCode*/ = 1) { st.cell.set(a, b, c, al, be, ga); }
  void Orth2Frac(realtype x, realtype y, realtype z, realtype &u, realtype &v, realtype &w) {
    gemmi::Fractional f = st.cell.fractionalize(gemmi::Position(x, y, z));
    u = f.x; v = f.y; w = f.z;
  }
  void Frac2Orth(realtype u, realtype v, realtype w, realtype &x, realtype &y, realtype &z) {
    gemmi::Position p = st.cell.orthogonalize(gemmi::Fractional(u, v, w));
    x = p.x; y = p.y; z = p.z;
  }
  pstr GetSpaceGroup()    { _sg_buf = st.spacegroup_hm; return (pstr) _sg_buf.c_str(); }
  pstr GetSpaceGroupFix() { return GetSpaceGroup(); }
  int  SetSpaceGroup(cpstr sg) { st.spacegroup_hm = sg ? sg : ""; return 0; }
  int  GetNumberOfSymOps() {
    const gemmi::SpaceGroup *sg = gemmi::find_spacegroup_by_name(st.spacegroup_hm);
    return sg ? (int) sg->operations().order() : 0;
  }
  pstr GetSymOp(int Nop) {
    const gemmi::SpaceGroup *sg = gemmi::find_spacegroup_by_name(st.spacegroup_hm);
    if (!sg) return nullptr;
    int i = 0;
    for (gemmi::Op op : sg->operations()) {
      if (i++ == Nop) { _symop_buf = op.triplet(); return (pstr) _symop_buf.c_str(); }
    }
    return nullptr;
  }

  // ---- selection ----
  struct Selection {
    SELECTION_TYPE type = STYPE_UNDEFINED;
    std::vector<Atom *>    atoms;
    std::vector<Residue *> residues;
  };
  std::vector<Selection> selections;  // handle is 1-based index

  int  NewSelection() { selections.emplace_back(); return (int)selections.size(); }
  void DeleteSelection(int selHnd) {
    if (selHnd < 1 || selHnd > (int)selections.size()) return;
    Selection &s = selections[selHnd - 1];
    for (Atom *a : s.atoms) a->_setInSel(selHnd, false);
    for (Residue *r : s.residues) r->_setInSel(selHnd, false);
    s = Selection();
  }
  void GetSelIndex(int selHnd, PPAtom &SelAtom, int &n) {
    Selection &s = selections[selHnd - 1]; SelAtom = s.atoms.data(); n = (int)s.atoms.size();
  }
  void GetSelIndex(int selHnd, PPResidue &SelRes, int &n) {
    Selection &s = selections[selHnd - 1]; SelRes = s.residues.data(); n = (int)s.residues.size();
  }
  // chain selections aren't modelled by the engine yet — return empty. TODO.
  void GetSelIndex(int /*selHnd*/, PPChain &SelChain, int &n) { SelChain = nullptr; n = 0; }
  // select atoms by serial-number range (iSer1..iSer2; 0,0 => all).
  void SelectAtoms(int selHnd, int iSer1, int iSer2, SELECTION_KEY key) {
    if (selHnd < 1 || selHnd > (int)selections.size()) return;
    Selection &s = selections[selHnd - 1];
    std::vector<Atom *> pick;
    for (Atom *a : all_atoms) {
      int sn = a->g().serial;
      if ((iSer1 == 0 && iSer2 == 0) || (sn >= iSer1 && sn <= iSer2)) pick.push_back(a);
    }
    if (key == SKEY_OR) { for (Atom *a : pick) if (!a->isInSelection(selHnd)) s.atoms.push_back(a); }
    else { for (Atom *a : s.atoms) a->_setInSel(selHnd, false); s.atoms = pick; }
    s.type = STYPE_ATOM;
    for (Atom *a : s.atoms) a->_setInSel(selHnd, true);
  }

  // full spatial+CID atom selection (mmdb_selmngr.h) — sphere around (x,y,z) with
  // chain/resname/atomname/element filters ("!X" = exclusion, "*" = any).
  void SelectAtoms(int selHnd, int /*iModel*/, cpstr Chains, int ResNo1, cpstr /*Ins1*/,
                   int ResNo2, cpstr /*Ins2*/, cpstr RNames, cpstr ANames, cpstr Elements,
                   cpstr /*altLocs*/, cpstr /*segIDs*/, cpstr /*charges*/,
                   realtype /*occ1*/, realtype /*occ2*/, realtype x, realtype y, realtype z,
                   realtype radius, SELECTION_KEY key) {
    if (selHnd < 1 || selHnd > (int)selections.size()) return;
    Selection &s = selections[selHnd - 1];
    // self-contained comma-list matcher ("*"=any, "!X"=exclude); `detail::` is
    // declared after Manager, so don't depend on it in this inline body.
    auto inlist = [](cpstr list, const std::string &v) -> bool {
      if (!list || !*list || std::strcmp(list, "*") == 0) return true;
      for (const char *p = list; *p; ) {
        const char *c = std::strchr(p, ',');
        std::string tok(p, c ? (size_t)(c - p) : std::strlen(p));
        size_t a = tok.find_first_not_of(' '), b = tok.find_last_not_of(' ');
        tok = (a == std::string::npos) ? std::string() : tok.substr(a, b - a + 1);
        if (tok == v) return true;
        if (!c) break; p = c + 1;
      }
      return false;
    };
    auto match = [&](cpstr list, const std::string &v) -> bool {
      if (!list || !*list || std::strcmp(list, "*") == 0) return true;
      if (list[0] == '!') return !inlist(list + 1, v);
      return inlist(list, v);
    };
    gemmi::Position pt(x, y, z); double r2 = radius * radius;
    std::vector<Atom *> pick;
    for (Atom *a : all_atoms) {
      if (radius > 0 && a->g().pos.dist_sq(pt) > r2) continue;
      Residue *r = a->res;
      int sn = r->GetSeqNum();
      if (ResNo1 != ANY_RES && sn < ResNo1) continue;
      if (ResNo2 != ANY_RES && sn > ResNo2) continue;
      if (!match(Chains, r->chain->g().name)) continue;
      if (!match(RNames, std::string(r->GetResName()))) continue;
      if (!match(ANames, std::string(a->GetAtomName()))) continue;
      if (!match(Elements, gemmi::Element(a->g().element).name())) continue;
      pick.push_back(a);
    }
    if (key == SKEY_OR) { for (Atom *a : pick) if (!a->isInSelection(selHnd)) s.atoms.push_back(a); }
    else { for (Atom *a : s.atoms) a->_setInSel(selHnd, false); s.atoms = pick; }
    s.type = STYPE_ATOM;
    for (Atom *a : s.atoms) a->_setInSel(selHnd, true);
  }

  // --- misc hierarchy/bond/UDData ops used by Coot ---
  void RemoveBonds() {}                       // gemmi has no persistent bond table
  void Delete(int /*DelKey*/) {}              // partial-hierarchy delete — no-op (TODO)
  void DeleteAllModels() { st.models.clear(); build_from_gemmi(); }   // clears the hierarchy
  void DeleteModel(int modelNo) {   // 1-based; erase model + rebuild wrappers
    int i = modelNo - 1;
    if (i >= 0 && i < (int)st.models.size()) { st.models.erase(st.models.begin() + i); build_from_gemmi(); }
  }
  pstr GetInputBuffer(pstr buf, int &count) { count = 0; if (buf) buf[0] = '\0'; return buf; }
  // place an atom into the flat table (mmdb Manager::PutAtom) — the shim builds
  // hierarchy via Add*/gemmi, so this is a stub returning the index. TODO if a
  // PutAtom-built molecule is needed.
  int  PutAtom(int index, PAtom /*atom*/, int /*serNum*/ = 0) { return index; }
  // hierarchy-level UDData (UDR_HIERARCHY) — Manager owns its own UDStore.
  UDStore _ud;
  int PutUDData(int h, int v)      { return ud_put(this, UDR_HIERARCHY, _ud, h, v); }
  int PutUDData(int h, realtype v) { return ud_put(this, UDR_HIERARCHY, _ud, h, v); }
  int PutUDData(int h, cpstr v)    { return ud_put(this, UDR_HIERARCHY, _ud, h, v); }
  int GetUDData(int h, int &v)      { return ud_get(this, UDR_HIERARCHY, _ud, h, v); }
  int GetUDData(int h, realtype &v) { return ud_get(this, UDR_HIERARCHY, _ud, h, v); }
  int GetUDData(int h, pstr &v)     { return ud_get(this, UDR_HIERARCHY, _ud, h, v); }
  // primary CID-range selection (STYPE via Select; SelectAtoms forwards as STYPE_ATOM)
  void Select(int selHnd, SELECTION_TYPE sType, int iModel, cpstr Chains,
              int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2, cpstr RNames,
              cpstr ANames, cpstr Elements, cpstr altLocs, SELECTION_KEY selKey = SKEY_OR);
  void SelectAtoms(int selHnd, int iModel, cpstr Chains, int ResNo1, cpstr Ins1,
                   int ResNo2, cpstr Ins2, cpstr RNames, cpstr ANames,
                   cpstr Elements, cpstr altLocs, SELECTION_KEY selKey = SKEY_OR) {
    Select(selHnd, STYPE_ATOM, iModel, Chains, ResNo1, Ins1, ResNo2, Ins2,
           RNames, ANames, Elements, altLocs, selKey);
  }
  void SelectSphere(int selHnd, SELECTION_TYPE sType, realtype x, realtype y,
                    realtype z, realtype r, SELECTION_KEY sKey = SKEY_OR);
  // select-from-selection: combine selHnd2's contents into selHnd1 per sKey
  void Select(int selHnd1, SELECTION_TYPE sType, int selHnd2, SELECTION_KEY sKey);
  // atoms within [d1,d2] of any atom in the given set (defined in contacts.cc)
  void SelectNeighbours(int selHnd, SELECTION_TYPE sType, PPAtom atoms, int nAtoms,
                        realtype d1, realtype d2, SELECTION_KEY sKey = SKEY_OR);
  void SetFlag(int /*flags*/) {}          // no-op: read/write behaviour is fixed
  void SetFlag(cpstr /*flags*/) {}
  int  PutPDBString(cpstr /*card*/) { return Error_NoError; } // no-op
  int  MakeBonds(bool /*calc*/) { return 0; }                // TODO: gemmi bonds

  // flat atom access (across the whole hierarchy)
  std::vector<Atom *> all_atoms;
  int   GetNumberOfAtoms() { return (int)all_atoms.size(); }
  int   GetNumberOfAtoms(bool /*countTers*/) { return (int)all_atoms.size(); }
  int   GetNumberOfAtoms(cpstr CID);      // count atoms matching CID (defined below)
  PAtom GetAtomI(int i) { return (i >= 0 && i < (int)all_atoms.size()) ? all_atoms[i] : nullptr; }
  void  GetAtomTable(PPAtom &t, int &n) { t = all_atoms.data(); n = (int)all_atoms.size(); }
  void  GetModelTable(PPModel &t, int &n) { t = models.data(); n = (int)models.size(); }
  void  GetAtomStatistics(int selHnd, RAtomStat AS);   // defined below
  int   MakeSelIndex(int selHnd) {
    return (selHnd >= 1 && selHnd <= (int)selections.size())
               ? (int)selections[selHnd - 1].atoms.size() : 0;
  }
  void  SelectAtom(int selHnd, PAtom atom, SELECTION_KEY sKey, bool makeIndex = true);
  // CID-string selection, e.g. "/1/A/10-20/CA"
  void  Select(int selHnd, SELECTION_TYPE sType, cpstr CID, SELECTION_KEY sKey);

  // ---- contacts (gemmi-free uniform-grid search) ----
  // TMatrix is MMDB's optional symmetry transform applied to the 2nd set; the
  // shim does no symmetry (see contacts.cc image_idx!=0 exclusion), so it's
  // accepted and ignored. TODO: gemmi symmetry-aware contacts.
  void SeekContacts(PPAtom A1, int n1, PPAtom A2, int n2, realtype d1,
                    realtype d2, int seqDist, PContact &contact, int &ncontacts,
                    int maxlen = 0, pmat44 TMatrix = nullptr, long group = 0);
  void SeekContacts(PPAtom A, int n, realtype d1, realtype d2, int seqDist,
                    PContact &contact, int &ncontacts, int maxlen = 0,
                    pmat44 TMatrix = nullptr, long group = 0);
  // single-atom vs selection (forwards to the array overload with a 1-elem array)
  void SeekContacts(PAtom a, PPAtom A2, int n2, realtype d1, realtype d2, int seqDist,
                    PContact &contact, int &ncontacts, int maxlen = 0,
                    pmat44 TMatrix = nullptr, long group = 0) {
    PAtom a1[1] = { a };
    SeekContacts(a1, 1, A2, n2, d1, d2, seqDist, contact, ncontacts, maxlen, TMatrix, group);
  }

  int  FinishStructEdit() { return 0; }  // no-op: wrappers stay in sync eagerly

  // ---- UDData registry ----
  struct UDReg { UDR_TYPE type; int kind; std::string name; int slot; }; // kind:0=int,1=real,2=str
  std::vector<UDReg> ud_regs;
  int ud_counts[5][3] = {{0}};   // [UDR_TYPE][kind] -> next slot

  int RegisterUDInteger(UDR_TYPE t, cpstr name) { return _regUD(t, 0, name); }
  int RegisterUDReal   (UDR_TYPE t, cpstr name) { return _regUD(t, 1, name); }
  int RegisterUDString (UDR_TYPE t, cpstr name) { return _regUD(t, 2, name); }
  int GetUDDHandle(UDR_TYPE t, cpstr name) {
    for (int i = 0; i < (int)ud_regs.size(); ++i)
      if (ud_regs[i].type == t && ud_regs[i].name == name) return i;
    return -1;
  }
private:
  int _regUD(UDR_TYPE t, int kind, cpstr name) {
    ud_regs.push_back({t, kind, name ? name : "", ud_counts[t][kind]++});
    return (int)ud_regs.size() - 1;
  }
public:

  // ---- I/O (defined in mmdb-shim/src/io.cc; keeps heavy gemmi write/read
  //      headers out of the ~229 Coot TUs that include mmdb_manager.h) ----
  ERROR_CODE ReadPDBASCII(cpstr fname);
  ERROR_CODE ReadCoorFile(cpstr fname);    // auto-detects PDB / mmCIF
  ERROR_CODE WritePDBASCII(cpstr fname);
  ERROR_CODE WriteCIFASCII(cpstr fname);
};

// ---- g() resolvers ----
// A wrapper with no parent is "detached" (Coot's `new mmdb::Atom` idiom: build
// standalone, set fields, then Add*() into a parent). While detached, g()
// resolves to a wrapper-owned local gemmi object; Add*() copies that local into
// the parent's gemmi vector and rebinds (sets parent + index). Index-based
// resolution makes the vector push/reallocation harmless for siblings.
inline gemmi::Model   &Model::g()   const { return mgr ? mgr->st.models[mi] : const_cast<Model*>(this)->_local; }
inline gemmi::Chain   &Chain::g()   const { return model ? model->g().chains[ci]     : const_cast<Chain*>(this)->_local; }
inline gemmi::Residue &Residue::g() const { return chain ? chain->g().residues[ri]   : const_cast<Residue*>(this)->_local; }
inline gemmi::Atom    &Atom::g()    const { return res   ? res->g().atoms[ai]         : const_cast<Atom*>(this)->_local; }

// ---- UDData helpers ----
inline Manager::UDReg *_ud_desc(Manager *mgr, UDR_TYPE myType, int handle, int kind,
                                int &err) {
  if (!mgr || handle < 0 || handle >= (int)mgr->ud_regs.size()) { err = UDDATA_WrongHandle; return nullptr; }
  Manager::UDReg &d = mgr->ud_regs[handle];
  if (d.type != myType || d.kind != kind) { err = UDDATA_WrongUDRType; return nullptr; }
  err = UDDATA_Ok; return &d;
}
inline int ud_put(Manager *mgr, UDR_TYPE t, UDStore &s, int h, int v) {
  int e; auto *d = _ud_desc(mgr, t, h, 0, e); if (!d) return e;
  if ((int)s._udi.size() <= d->slot) s._udi.resize(d->slot + 1, 0);
  s._udi[d->slot] = v; return UDDATA_Ok;
}
inline int ud_put(Manager *mgr, UDR_TYPE t, UDStore &s, int h, realtype v) {
  int e; auto *d = _ud_desc(mgr, t, h, 1, e); if (!d) return e;
  if ((int)s._udr.size() <= d->slot) s._udr.resize(d->slot + 1, 0.0);
  s._udr[d->slot] = v; return UDDATA_Ok;
}
inline int ud_put(Manager *mgr, UDR_TYPE t, UDStore &s, int h, cpstr v) {
  int e; auto *d = _ud_desc(mgr, t, h, 2, e); if (!d) return e;
  if ((int)s._uds.size() <= d->slot) s._uds.resize(d->slot + 1);
  s._uds[d->slot] = v ? v : ""; return UDDATA_Ok;
}
inline int ud_get(Manager *mgr, UDR_TYPE t, UDStore &s, int h, int &v) {
  int e; auto *d = _ud_desc(mgr, t, h, 0, e); if (!d) return e;
  if ((int)s._udi.size() <= d->slot) return UDDATA_NoData;
  v = s._udi[d->slot]; return UDDATA_Ok;
}
inline int ud_get(Manager *mgr, UDR_TYPE t, UDStore &s, int h, realtype &v) {
  int e; auto *d = _ud_desc(mgr, t, h, 1, e); if (!d) return e;
  if ((int)s._udr.size() <= d->slot) return UDDATA_NoData;
  v = s._udr[d->slot]; return UDDATA_Ok;
}
inline int ud_get(Manager *mgr, UDR_TYPE t, UDStore &s, int h, pstr &v) {
  int e; auto *d = _ud_desc(mgr, t, h, 2, e); if (!d) return e;
  if ((int)s._uds.size() <= d->slot) return UDDATA_NoData;
  v = (pstr) s._uds[d->slot].c_str(); return UDDATA_Ok;   // borrowed
}

// ---- Atom out-of-line ----
inline pstr Atom::GetAtomName() const {
  std::snprintf(_name_buf, sizeof(_name_buf), "%s", g().name.c_str());
  return _name_buf;
}
inline void Atom::SetAtomName(const AtomName aName) { g().name = aName; }
inline pstr Atom::GetElementName() {
  std::snprintf(_elem_buf, sizeof(_elem_buf), "%s", g().element.name());
  return _elem_buf;
}
inline void Atom::SetElementName(const Element elName) { g().element = gemmi::Element(elName); }
inline pstr Atom::GetChainID()  { return res->GetChainID(); }
inline int   Atom::GetSeqNum()   { return res->GetSeqNum(); }
inline Chain *Atom::GetChain()   { return res ? res->GetChain() : nullptr; }
inline Model *Atom::GetModel()   { return res ? res->chain->model : nullptr; }
inline pstr Atom::GetLabelCompID()   { return res ? res->GetLabelCompID() : nullptr; }
inline pstr Atom::GetLabelAsymID()   { return res ? res->GetLabelAsymID() : nullptr; }
inline int  Atom::GetLabelSeqID()    { return res ? res->GetLabelSeqID() : 0; }
inline int  Atom::GetLabelEntityID() { return res ? res->GetLabelEntityID() : 0; }
inline int  Atom::GetResidueNo()     { return res ? res->GetResidueNo() : 0; }
inline int  Atom::GetSSEType()       { return res ? res->SSE : SSE_None; }
inline bool Atom::isSolvent()        { return res ? res->isSolvent() : false; }
inline bool Atom::isNTerminus()      { return res ? res->isNTerminus() : false; }
inline bool Atom::isCTerminus()      { return res ? res->isCTerminus() : false; }
inline pstr  Atom::GetInsCode()  { return res->GetInsCode(); }
inline pstr  Atom::GetResName()  { return res->GetResName(); }
inline int   Atom::GetModelNum() { return res->GetModelNum(); }
inline int   Atom::GetIndex()    { return ai; }
inline void  Atom::SetCoordinates(realtype xx, realtype yy, realtype zz,
                                  realtype occ, realtype tF) {
  auto &a = g(); a.pos = gemmi::Position(xx, yy, zz); a.occ = (float)occ; a.b_iso = (float)tF;
}

// ---- Residue out-of-line ----
inline pstr Residue::GetResName() {
  std::snprintf(_resname_buf, sizeof(_resname_buf), "%s", g().name.c_str());
  return _resname_buf;
}
inline int &Residue::GetSeqNum()  { return g().seqid.num.value; }
inline pstr Residue::GetInsCode() {
  _inscode_buf[0] = g().seqid.icode == ' ' ? '\0' : g().seqid.icode; _inscode_buf[1] = '\0';
  return _inscode_buf;
}
inline pstr Residue::GetChainID() { return chain->GetChainID(); }
inline int   Residue::GetModelNum() { return chain->model->GetSerNum(); }
inline PAtom Residue::GetAtom(const AtomName aname, const Element elname, const AltLoc aloc) {
  for (Atom *a : atoms) {
    if (a->g().name != aname) continue;
    if (elname && *elname && a->g().element.name() != std::string(elname)) continue;
    if (aloc && *aloc && a->g().altloc != aloc[0]) continue;
    return a;
  }
  return nullptr;
}
inline PAtom Residue::AddAtom(Manager &m, gemmi::Atom a) {
  g().atoms.push_back(std::move(a));
  Atom *aw = m.newAtom(); aw->mgr = &m; aw->res = this; aw->ai = (int)atoms.size();
  atoms.push_back(aw);
  return aw;
}
inline void Residue::DeleteAtom(int pos) {
  if (pos < 0 || pos >= (int)atoms.size()) return;
  g().atoms.erase(g().atoms.begin() + pos);
  atoms[pos]->alive = false; atoms[pos]->ai = -1;
  atoms.erase(atoms.begin() + pos);
  for (int k = pos; k < (int)atoms.size(); ++k) atoms[k]->ai = k;
}

// ---- Chain out-of-line ----
inline bool Chain::isAminoacidChain() {
  for (Residue *r : residues) if (r->isAminoacid()) return true;
  return false;
}
inline bool Chain::isNucleotideChain() {
  for (Residue *r : residues) if (r->isNucleotide()) return true;
  return false;
}
inline bool Chain::isSolventChain() {
  if (residues.empty()) return false;
  for (Residue *r : residues) if (!r->isSolvent()) return false;
  return true;
}
inline pstr Chain::GetChainID() {
  std::snprintf(_chainid_buf, sizeof(_chainid_buf), "%s", g().name.c_str());
  return _chainid_buf;
}
inline PResidue Chain::AddResidue(Manager &m, gemmi::Residue r) {
  g().residues.push_back(std::move(r));
  Residue *rw = m.newRes(); rw->mgr = &m; rw->chain = this; rw->ri = (int)residues.size();
  for (int ai = 0; ai < (int)rw->g().atoms.size(); ++ai) {
    Atom *aw = m.newAtom(); aw->mgr = &m; aw->res = rw; aw->ai = ai; rw->atoms.push_back(aw);
  }
  residues.push_back(rw);
  return rw;
}
inline PResidue Chain::InsResidue(Manager &m, int pos, gemmi::Residue r) {
  g().residues.insert(g().residues.begin() + pos, std::move(r));
  Residue *rw = m.newRes(); rw->mgr = &m; rw->chain = this; rw->ri = pos;
  residues.insert(residues.begin() + pos, rw);
  for (int k = pos + 1; k < (int)residues.size(); ++k) residues[k]->ri = k;
  for (int ai = 0; ai < (int)rw->g().atoms.size(); ++ai) {
    Atom *aw = m.newAtom(); aw->mgr = &m; aw->res = rw; aw->ai = ai; rw->atoms.push_back(aw);
  }
  return rw;
}

// ---- Model out-of-line ----
inline PChain Model::GetChain(const ChainID chID) {
  for (Chain *c : chains) if (c->g().name == chID) return c;
  return nullptr;
}

// ---- Manager out-of-line ----
inline void Manager::build_from_gemmi() {
  models.clear();
  all_atoms.clear();
  for (int mi = 0; mi < (int)st.models.size(); ++mi) {
    Model *mw = newModel(); mw->mgr = this; mw->mi = mi;
    auto &gm = st.models[mi];
    for (int ci = 0; ci < (int)gm.chains.size(); ++ci) {
      Chain *cw = newChain(); cw->mgr = this; cw->model = mw; cw->ci = ci;
      auto &gc = gm.chains[ci];
      for (int ri = 0; ri < (int)gc.residues.size(); ++ri) {
        Residue *rw = newRes(); rw->mgr = this; rw->chain = cw; rw->ri = ri;
        auto &gr = gc.residues[ri];
        for (int ai = 0; ai < (int)gr.atoms.size(); ++ai) {
          Atom *aw = newAtom(); aw->mgr = this; aw->res = rw; aw->ai = ai;
          aw->WhatIsSet = ASET_Coordinates | ASET_Occupancy | ASET_tempFactor;
          const gemmi::SMat33<float> &an = gr.atoms[ai].aniso;
          if (an.u11 != 0.f || an.u22 != 0.f || an.u33 != 0.f) aw->WhatIsSet |= ASET_Anis_tFac;
          rw->atoms.push_back(aw); all_atoms.push_back(aw); mw->all_atoms.push_back(aw);
        }
        rw->_sync_atom();
        rw->_load_id();
        cw->residues.push_back(rw);
      }
      mw->chains.push_back(cw);
    }
    models.push_back(mw);
  }
}

// ---- selection matching ----
namespace detail {
inline bool inList(cpstr list, const std::string &v) {
  if (!list || !*list || std::strcmp(list, "*") == 0) return true;
  const char *p = list;
  while (*p) {
    const char *c = std::strchr(p, ',');
    std::string tok(p, c ? (size_t)(c - p) : std::strlen(p));
    size_t a = tok.find_first_not_of(' '), b = tok.find_last_not_of(' ');
    tok = (a == std::string::npos) ? std::string() : tok.substr(a, b - a + 1);
    if (tok == v) return true;
    if (!c) break; p = c + 1;
  }
  return false;
}
inline bool altMatch(cpstr list, char alt) {
  if (!list || std::strcmp(list, "*") == 0) return true;
  std::string a = alt ? std::string(1, alt) : std::string();
  if (!*list) return a.empty();          // "" -> only blank altLoc
  return inList(list, a);
}
} // namespace detail

inline void Manager::Select(int selHnd, SELECTION_TYPE sType, int iModel,
    cpstr Chains, int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2, cpstr RNames,
    cpstr ANames, cpstr Elements, cpstr altLocs, SELECTION_KEY selKey) {
  (void)Ins1; (void)Ins2;  // insertion-code range filtering: TODO (rare in Coot)
  Selection &sel = selections[selHnd - 1];
  if (sel.type == STYPE_UNDEFINED) sel.type = sType;
  std::vector<Atom *> oldA = sel.atoms; std::vector<Residue *> oldR = sel.residues;

  std::vector<Atom *> mAtoms; std::vector<Residue *> mResidues;
  for (Model *mw : models) {
    if (iModel > 0 && mw->GetSerNum() != iModel) continue;
    for (Chain *cw : mw->chains) {
      if (!detail::inList(Chains, cw->g().name)) continue;
      for (Residue *rw : cw->residues) {
        int sn = rw->g().seqid.num.value;
        if (ResNo1 != ANY_RES && sn < ResNo1) continue;
        if (ResNo2 != ANY_RES && sn > ResNo2) continue;
        if (!detail::inList(RNames, rw->g().name)) continue;
        bool anyAtom = false;
        for (Atom *aw : rw->atoms) {
          if (!detail::inList(ANames, aw->g().name)) continue;
          if (!detail::inList(Elements, aw->g().element.name())) continue;
          if (!detail::altMatch(altLocs, aw->g().altloc)) continue;
          anyAtom = true;
          if (sType == STYPE_ATOM) mAtoms.push_back(aw);
        }
        if (anyAtom && sType == STYPE_RESIDUE) mResidues.push_back(rw);
      }
    }
  }
  auto combine = [&](auto &cur, auto &matched) {
    using Vec = typename std::decay<decltype(cur)>::type;
    std::set<typename Vec::value_type> curset(cur.begin(), cur.end());
    std::set<typename Vec::value_type> mset(matched.begin(), matched.end());
    if (selKey == SKEY_NEW) { cur = matched; }
    else if (selKey == SKEY_OR) { for (auto *x : matched) if (!curset.count(x)) cur.push_back(x); }
    else if (selKey == SKEY_AND) { Vec o; for (auto *x : cur) if (mset.count(x)) o.push_back(x); cur = o; }
    else if (selKey == SKEY_XOR) { Vec o; for (auto *x : cur) if (!mset.count(x)) o.push_back(x);
                                   for (auto *x : matched) if (!curset.count(x)) o.push_back(x); cur = o; }
    else if (selKey == SKEY_CLR) { Vec o; for (auto *x : cur) if (!mset.count(x)) o.push_back(x); cur = o; }
  };
  if (sType == STYPE_ATOM) combine(sel.atoms, mAtoms);
  else if (sType == STYPE_RESIDUE) combine(sel.residues, mResidues);
  for (Atom *a : oldA) a->_setInSel(selHnd, false);
  for (Atom *a : sel.atoms) a->_setInSel(selHnd, true);
  for (Residue *r : oldR) r->_setInSel(selHnd, false);
  for (Residue *r : sel.residues) r->_setInSel(selHnd, true);
}

// select-from-selection: combine selHnd2's contents into selHnd1
inline void Manager::Select(int selHnd1, SELECTION_TYPE sType, int selHnd2,
                            SELECTION_KEY sKey) {
  Selection &s1 = selections[selHnd1 - 1];
  Selection &s2 = selections[selHnd2 - 1];
  if (s1.type == STYPE_UNDEFINED) s1.type = sType;
  std::vector<Atom *> oldA = s1.atoms; std::vector<Residue *> oldR = s1.residues;
  auto combine = [&](auto &cur, auto &m) {
    using Vec = typename std::decay<decltype(cur)>::type;
    std::set<typename Vec::value_type> curset(cur.begin(), cur.end());
    std::set<typename Vec::value_type> mset(m.begin(), m.end());
    if (sKey == SKEY_NEW) cur = m;
    else if (sKey == SKEY_OR) { for (auto *x : m) if (!curset.count(x)) cur.push_back(x); }
    else if (sKey == SKEY_AND) { Vec o; for (auto *x : cur) if (mset.count(x)) o.push_back(x); cur = o; }
    else if (sKey == SKEY_XOR) { Vec o; for (auto *x : cur) if (!mset.count(x)) o.push_back(x);
                                 for (auto *x : m) if (!curset.count(x)) o.push_back(x); cur = o; }
    else if (sKey == SKEY_CLR) { Vec o; for (auto *x : cur) if (!mset.count(x)) o.push_back(x); cur = o; }
  };
  if (sType == STYPE_ATOM) combine(s1.atoms, s2.atoms);
  else if (sType == STYPE_RESIDUE) combine(s1.residues, s2.residues);
  for (Atom *a : oldA) a->_setInSel(selHnd1, false);
  for (Atom *a : s1.atoms) a->_setInSel(selHnd1, true);
  for (Residue *r : oldR) r->_setInSel(selHnd1, false);
  for (Residue *r : s1.residues) r->_setInSel(selHnd1, true);
}

inline void Manager::SelectAtom(int selHnd, PAtom atom, SELECTION_KEY sKey, bool) {
  Selection &sel = selections[selHnd - 1];
  if (sel.type == STYPE_UNDEFINED) sel.type = STYPE_ATOM;
  if (sKey == SKEY_NEW) {
    for (Atom *a : sel.atoms) a->_setInSel(selHnd, false);
    sel.atoms.clear();
  }
  if (atom && !atom->isInSelection(selHnd)) {
    sel.atoms.push_back(atom); atom->_setInSel(selHnd, true);
  }
}

// Pragmatic CID parser: "/model/chain/seqNum1(-seqNum2)/atom" (best-effort;
// strips (resname)/[element]/:altloc suffixes). TODO: full MMDB CID grammar.
inline void Manager::Select(int selHnd, SELECTION_TYPE sType, cpstr CID,
                            SELECTION_KEY sKey) {
  std::string s = CID ? CID : "";
  std::vector<std::string> t;
  size_t p = (!s.empty() && s[0] == '/') ? 1 : 0;
  while (p <= s.size()) {
    size_t q = s.find('/', p);
    t.push_back(s.substr(p, q == std::string::npos ? std::string::npos : q - p));
    if (q == std::string::npos) break;
    p = q + 1;
  }
  auto tok = [&](size_t i) { return i < t.size() ? t[i] : std::string(); };
  auto strip = [](std::string v, const char *seps) {
    size_t c = v.find_first_of(seps); return c == std::string::npos ? v : v.substr(0, c);
  };
  int iModel = 0; std::string m = tok(0);
  if (!m.empty() && m != "*" && m != "0") iModel = atoi(m.c_str());
  std::string chains = tok(1).empty() ? "*" : tok(1);
  int r1 = ANY_RES, r2 = ANY_RES;
  std::string rr = strip(tok(2), "(");            // drop (resname)
  if (!rr.empty() && rr != "*") {
    size_t dash = rr.find('-', rr[0] == '-' ? 1 : 0);
    if (dash == std::string::npos) { r1 = r2 = atoi(rr.c_str()); }
    else { r1 = atoi(rr.substr(0, dash).c_str()); r2 = atoi(rr.substr(dash + 1).c_str()); }
  }
  std::string anames = strip(strip(tok(3), "["), ":");   // drop [element]/:altloc
  if (anames.empty()) anames = "*";
  Select(selHnd, sType, iModel, chains.c_str(), r1, "*", r2, "*", "*",
         anames.c_str(), "*", "*", sKey);
}

inline int Manager::GetNumberOfAtoms(cpstr CID) {
  int h = NewSelection();
  Select(h, STYPE_ATOM, CID, SKEY_NEW);
  int n = (int)selections[h - 1].atoms.size();
  DeleteSelection(h);
  return n;
}

inline void Manager::GetAtomStatistics(int selHnd, RAtomStat AS) {
  AS = AtomStat();
  std::vector<Atom *> &atoms = selections[selHnd - 1].atoms;
  AS.nAtoms = (int)atoms.size();
  if (atoms.empty()) return;
  double sx = 0, sy = 0, sz = 0;
  AS.xmin = AS.xmax = atoms[0]->x();
  AS.ymin = AS.ymax = atoms[0]->y();
  AS.zmin = AS.zmax = atoms[0]->z();
  for (Atom *a : atoms) {
    double X = a->x(), Y = a->y(), Z = a->z();
    sx += X; sy += Y; sz += Z;
    AS.xmin = X < AS.xmin ? X : AS.xmin; AS.xmax = X > AS.xmax ? X : AS.xmax;
    AS.ymin = Y < AS.ymin ? Y : AS.ymin; AS.ymax = Y > AS.ymax ? Y : AS.ymax;
    AS.zmin = Z < AS.zmin ? Z : AS.zmin; AS.zmax = Z > AS.zmax ? Z : AS.zmax;
  }
  AS.xm = sx / atoms.size(); AS.ym = sy / atoms.size(); AS.zm = sz / atoms.size();
}

inline void Manager::SelectSphere(int selHnd, SELECTION_TYPE sType, realtype x,
    realtype y, realtype z, realtype r, SELECTION_KEY sKey) {
  Selection &sel = selections[selHnd - 1];
  if (sel.type == STYPE_UNDEFINED) sel.type = sType;
  std::vector<Atom *> oldA = sel.atoms; std::vector<Residue *> oldR = sel.residues;
  gemmi::Position c(x, y, z); double r2 = r * r;
  std::vector<Atom *> mAtoms; std::vector<Residue *> mResidues;
  for (Model *mw : models)
    for (Chain *cw : mw->chains)
      for (Residue *rw : cw->residues) {
        bool any = false;
        for (Atom *aw : rw->atoms)
          if (aw->g().pos.dist_sq(c) <= r2) { any = true; if (sType == STYPE_ATOM) mAtoms.push_back(aw); }
        if (any && sType == STYPE_RESIDUE) mResidues.push_back(rw);
      }
  auto combine = [&](auto &cur, auto &m) {
    std::set<typename std::decay<decltype(cur)>::type::value_type> cs(cur.begin(), cur.end());
    if (sKey == SKEY_NEW) cur = m;
    else if (sKey == SKEY_OR) { for (auto *p : m) if (!cs.count(p)) cur.push_back(p); }
  };
  if (sType == STYPE_ATOM) combine(sel.atoms, mAtoms);
  else if (sType == STYPE_RESIDUE) combine(sel.residues, mResidues);
  for (Atom *a : oldA) a->_setInSel(selHnd, false);
  for (Atom *a : sel.atoms) a->_setInSel(selHnd, true);
  for (Residue *r : oldR) r->_setInSel(selHnd, false);
  for (Residue *r : sel.residues) r->_setInSel(selHnd, true);
}

// SeekContacts (both overloads) is defined in mmdb-shim/src/contacts.cc using
// gemmi::NeighborSearch — keeps the heavy neighbor.hpp out of Coot's many TUs.

// ---- detached-construction constructors + subtree ops (need complete types) ----
inline Atom::Atom(Residue *r) { if (r) r->AddAtom(this); }
inline Residue::Residue(Chain *c) { if (c) c->AddResidue(this); }
inline Chain::Chain(Model *m, const ChainID id) { if (m) m->AddChain(this); SetChainID(id); }

inline bool Residue::isCTerminus() {
  return chain && ri == (int)chain->residues.size() - 1;
}
inline Model *Residue::GetModel() { return chain ? chain->model : nullptr; }

inline void Chain::Copy(PChain src) {
  Manager *pool = mgr ? mgr : src->mgr;
  g() = src->g();                    // deep gemmi copy (residues + atoms)
  residues.clear();
  if (!pool) return;
  gemmi::Chain &gc = g();
  for (int r = 0; r < (int)gc.residues.size(); ++r) {
    Residue *rw = pool->newRes(); rw->mgr = mgr; rw->chain = this; rw->ri = r;
    for (int a = 0; a < (int)gc.residues[r].atoms.size(); ++a) {
      Atom *aw = pool->newAtom(); aw->mgr = mgr; aw->res = rw; aw->ai = a;
      rw->atoms.push_back(aw);
    }
    rw->_sync_atom(); rw->_load_id();
    residues.push_back(rw);
  }
}

inline void Model::Copy(PModel src) {
  Manager *pool = mgr ? mgr : src->mgr;
  g() = src->g();
  chains.clear();
  if (!pool) return;
  gemmi::Model &gm = g();
  for (int c = 0; c < (int)gm.chains.size(); ++c) {
    Chain *cw = pool->newChain(); cw->mgr = mgr; cw->model = this; cw->ci = c;
    for (int r = 0; r < (int)gm.chains[c].residues.size(); ++r) {
      Residue *rw = pool->newRes(); rw->mgr = mgr; rw->chain = cw; rw->ri = r;
      for (int a = 0; a < (int)gm.chains[c].residues[r].atoms.size(); ++a) {
        Atom *aw = pool->newAtom(); aw->mgr = mgr; aw->res = rw; aw->ai = a;
        rw->atoms.push_back(aw);
      }
      rw->_sync_atom(); rw->_load_id();
      cw->residues.push_back(rw);
    }
    chains.push_back(cw);
  }
}

inline PChain Model::CreateChain(const ChainID id) {
  Chain *c = mgr ? mgr->newChain() : new Chain();
  c->mgr = mgr; c->model = this; c->ci = (int)chains.size();
  g().chains.emplace_back(id ? id : "");
  chains.push_back(c);
  return c;
}

inline pstr Atom::GetAtomID(pstr S) {
  if (S) std::snprintf(S, 100, "/%d/%s/%d(%s)/%s", GetModelNum(), GetChainID(),
                       res ? res->GetSeqNum() : 0, GetResName(), GetAtomName());
  return S;
}

// one-letter residue code (mmdb_tables.h) via gemmi's tabulated residues
inline void Get1LetterCode(cpstr res3, pstr res1) {
  if (!res1) return;
  char c = gemmi::find_tabulated_residue(res3 ? res3 : "").one_letter_code;
  res1[0] = c ? (char) std::toupper((unsigned char) c) : 'X'; res1[1] = '\0';
}
inline void Get1LetterCode(cpstr res3, char &res1) { char b[2]; Get1LetterCode(res3, b); res1 = b[0]; }

// sort a contact array by distance (mmdb_coormngr.h SortContacts) — sortkey ignored
inline void SortContacts(PContact contacts, int nContacts, int /*sortkey*/) {
  if (contacts && nContacts > 1)
    std::sort(contacts, contacts + nContacts,
              [](const Contact &a, const Contact &b) { return a.dist < b.dist; });
}

// centroid of an atom array (mmdb_coormngr.h GetMassCenter)
inline void GetMassCenter(PPAtom A, int nA, realtype &xc, realtype &yc, realtype &zc) {
  double sx = 0, sy = 0, sz = 0; int n = 0;
  for (int i = 0; i < nA; ++i) if (A[i]) { sx += A[i]->x(); sy += A[i]->y(); sz += A[i]->z(); ++n; }
  if (n) { xc = sx / n; yc = sy / n; zc = sz / n; } else { xc = yc = zc = 0; }
}

} // namespace mmdb

// mmdb::mmcif::* (thin veneer over gemmi::cif) — re-opens mmdb{mmcif{...}}.
// pstr/cpstr/realtype are already in scope from the headers above.
#include "_mmcif_impl.hh"

// mmdb::math::{Vertex,Edge,Graph,GraphMatch} — molecular graph + subgraph match.
// Included after the mmdb namespace close so Atom/Residue are complete (MakeGraph).
#include "_graph_impl.hh"
