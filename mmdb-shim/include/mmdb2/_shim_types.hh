// mmdb-shim — layer 1 of 4: scalar types, enums, and small record classes.
//
// This is the base layer of the gemmi-backed MMDB API shim (architecture B; see
// MMDB_SHIM_Recon_and_Plan.md). It provides MMDB's typedefs/enums, the per-object
// UDData+selection base (UDStore), forward declarations of the hierarchy classes,
// and the leaf "record" classes that carry no hierarchy (LINK / CISPEP / Cryst /
// Helix / Sheet / SymOps / Title …). Included first by the layers below.
#pragma once

#include <gemmi/model.hpp>
#include <gemmi/util.hpp>  // trim_str — normalise MMDB-padded names to gemmi's trimmed form
#include <gemmi/resinfo.hpp>
#include <gemmi/symmetry.hpp>  // space-group / symmetry operators

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
   typedef double realtype;
   typedef char *pstr;
   typedef const char *cpstr;
   typedef unsigned short word;
   typedef char AtomName[20];
   typedef char ResName[20];
   typedef char InsCode[10];
   typedef char ChainID[10];
   typedef char Element[10];
   typedef char AltLoc[20];
   typedef char SegID[10];
   typedef char LinkRID[20];    // Refmac link ID
   typedef unsigned char byte;  // mmdb_mattype.h
   typedef int *ivector;        // mmdb_mattype.h 1-based vectors/matrices
   typedef realtype *rvector;
   typedef ivector *imatrix;
   typedef rvector *rmatrix;
   typedef char maxMMDBName[40];

   // WhatIsSet mask flags (mmdb_atom.h ASET_FLAG)
   enum ASET_FLAG {
      ASET_Coordinates = 0x00000001,
      ASET_Occupancy = 0x00000002,
      ASET_tempFactor = 0x00000004,
      ASET_CoordSigma = 0x00000010,
      ASET_OccSigma = 0x00000020,
      ASET_tFacSigma = 0x00000040,
      ASET_Charge = 0x00000080,
      ASET_Anis_tFac = 0x00000100,
      ASET_Anis_tFSigma = 0x00001000,
      ASET_All = 0x000FFFFF
   };

   // vector/matrix types (mmdb_defs.h) — plain fixed-size arrays of realtype
   typedef realtype vect3[3];
   typedef realtype vect4[4];
   typedef vect3 mat33[3];  // realtype[3][3]
   typedef vect4 mat44[4];  // realtype[4][4]
   typedef mat44 *pmat44;
   typedef mat44 &rmat44;

   enum ERROR_CODE {
      Error_NoError = 0,
      Error_CantOpenFile = 12,  // matches real MMDB's value
      Error_GeneralError1 = 1
   };

   // ---- UDData (user-defined data) — real MMDB values (mmdb_uddata.h) ----
   enum UDR_TYPE { UDR_ATOM = 0,
                   UDR_RESIDUE = 1,
                   UDR_CHAIN = 2,
                   UDR_MODEL = 3,
                   UDR_HIERARCHY = 4 };
   enum UDDATA_CODE { UDDATA_Ok = 0,
                      UDDATA_WrongHandle = -1,
                      UDDATA_WrongUDRType = -2,
                      UDDATA_NoData = -3 };

   // ---- Selection (real MMDB values: mmdb_selmngr.h) ----
   enum SELECTION_TYPE { STYPE_INVALID = -1,
                         STYPE_UNDEFINED = 0,
                         STYPE_ATOM = 1,
                         STYPE_RESIDUE = 2,
                         STYPE_CHAIN = 3,
                         STYPE_MODEL = 4 };
   enum SELECTION_KEY { SKEY_NEW = 0,
                        SKEY_OR = 1,
                        SKEY_AND = 2,
                        SKEY_XOR = 3,
                        SKEY_CLR = 4,
                        SKEY_XAND = 100 };
   inline const long int MinInt4 = -2147483647;
   inline const long int MaxInt4 = 2147483647;
   inline const int ANY_RES = -2147483647;  // real MMDB: extern const == MinInt4
   inline const double Pi = 3.14159265358979323846;

   // PDB/CIF read flags (mmdb_io_file.h). Values are arbitrary distinct bits — the
   // shim's SetFlag is a no-op, so only distinctness matters for Coot's bit ops.
   enum MMDB_READ_FLAG {
      MMDBF_AutoSerials = 0x00000001,
      MMDBF_IgnoreDuplSeqNum = 0x00000002,
      MMDBF_IgnoreBlankLines = 0x00000004,
      MMDBF_IgnoreRemarks = 0x00000008,
      MMDBF_IgnoreHash = 0x00000010,
      MMDBF_IgnoreNonCoorPDBErrors = 0x00000020,
      MMDBF_PrintCIFWarnings = 0x00000040,
      MMDBF_All = 0x0000FFFF
   };
   enum MMDB_FCM { MMDBFCM_None = 0,
                   MMDBFCM_All = 1,
                   MMDBFCM_Coord = 2,
                   MMDBFCM_Cryst = 4,
                   MMDBFCM_SC = 8 };
   typedef int COPY_MASK;  // Coot uses `COPY_MASK cm = MMDBFCM_All` + bit arithmetic

   // Per-object UDData slots + selection membership bits. Each registered UDData
   // handle maps to a (type,kind,slot); the object stores contiguous vectors
   // indexed by slot. `_inSel[selHnd-1]` = is this object in selection selHnd
   // (maintained by Manager::Select/SelectSphere/DeleteSelection).
   struct UDStore {
      std::vector<int> _udi;
      std::vector<double> _udr;
      std::vector<std::string> _uds;
      std::vector<bool> _inSel;
      bool isInSelection(int selHnd) const {
         return selHnd >= 1 && selHnd <= (int)_inSel.size() && _inSel[selHnd - 1];
      }
      void _setInSel(int selHnd, bool v) {
         if ((int)_inSel.size() < selHnd) _inSel.resize(selHnd, false);
         _inSel[selHnd - 1] = v;
      }
   };

   class Atom;
   class Residue;
   class Chain;
   class Model;
   class Manager;
   typedef Atom *PAtom;
   typedef Atom **PPAtom;
   typedef Residue *PResidue;
   typedef Residue **PPResidue;
   typedef Chain *PChain;
   typedef Chain **PPChain;
   typedef Model *PModel;
   typedef Model **PPModel;
   typedef Manager *PManager;
   typedef Manager **PPManager;

   struct Contact {
      int id1, id2;
      long group;
      realtype dist;
   };
   typedef Contact *PContact;

   // base for records held in MMDB containers (Title compound/author, LINK, …)
   class ContainerClass {
   public:
      virtual ~ContainerClass() {}
   };
   typedef ContainerClass *PContainerClass;

   // LINK record. Public data members mirror real MMDB (Coot reads them directly).
   // Populated from gemmi Structure::connections on load (Manager::_load_metadata);
   // Coot-created links are appended via Model::AddLink.
   class Link : public ContainerClass {
   public:
      AtomName atName1{}, atName2{};
      AltLoc aloc1{}, aloc2{};
      ResName resName1{}, resName2{};
      ChainID chainID1{}, chainID2{};
      InsCode insCode1{}, insCode2{};
      int seqNum1 = 0, seqNum2 = 0;
      int s1 = 1, i1 = 0, j1 = 0, k1 = 0;  // symmetry id of 1st atom
      int s2 = 1, i2 = 0, j2 = 0, k2 = 0;  // symmetry id of 2nd atom
      realtype dist = 0;
      void Copy(PContainerClass o) {
         if (auto *l = dynamic_cast<Link *>(o)) *this = *l;
      }
   };
   typedef Link *PLink;
   typedef Link **PPLink;

   // Refmac LINK record (mmdb_model.h LinkR). Public members mirror real MMDB;
   // populated from gemmi Connections carrying a link_id (Manager::_load_metadata).
   class LinkR {
   public:
      LinkRID linkRID{};
      AtomName atName1{}, atName2{};
      AltLoc aloc1{}, aloc2{};
      ResName resName1{}, resName2{};
      ChainID chainID1{}, chainID2{};
      int seqNum1 = 0, seqNum2 = 0;
      InsCode insCode1{}, insCode2{};
      realtype dist = 0;
   };
   typedef LinkR *PLinkR;
   typedef LinkR **PPLinkR;

   // CIS-peptide record (mmdb_model.h CisPep). Public members mirror real MMDB;
   // populated from gemmi Structure::cispeps on load (Manager::_load_metadata).
   class CisPep {
   public:
      int serNum = 0;
      ResName pep1{};
      ChainID chainID1{};
      int seqNum1 = 0;
      InsCode icode1{};
      ResName pep2{};
      ChainID chainID2{};
      int seqNum2 = 0;
      InsCode icode2{};
      int modNum = 0;
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
   // COMPND/AUTHOR line containers. The AUTHOR container is filled from gemmi
   // meta.authors on load and the TITLE string comes from Structure::get_info
   // ("_struct.title"); COMPND/JRNL have no structured gemmi home, so those
   // containers stay empty.
   class Compound : public ContainerClass {
   public:
      char Line[256] = {0};
   };
   typedef Compound *PCompound;
   class Author : public ContainerClass {
   public:
      char Line[256] = {0};
   };
   typedef Author *PAuthor;
   class Journal : public ContainerClass {
   public:
      char Line[256] = {0};
   };
   typedef Journal *PJournal;
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
      TitleContainer compound, author, journal;            // public so Coot's access_title can reach them
      TitleContainer *GetCompound() { return &compound; }  // real Title exposes these
      TitleContainer *GetAuthor() { return &author; }      // publicly; access_title
      TitleContainer *GetJournal() { return &journal; }    // inherits GetJournal()
   };

   // gzip mode flag (mmdb_io_file.h). Minimal mmdb::io — the shim does I/O via gemmi,
   // so only this compression-mode enum is provided (Coot passes it to write calls).
   namespace io {
      enum GZ_MODE { GZM_NONE = 0,
                     GZM_CHECK = 1,
                     GZM_ENFORCE = 2 };
   }

   // initialise a 4x4 matrix to identity (mmdb_mattype.h Mat4Init)
   inline void Mat4Init(mat44 &A) {
      for (int i = 0; i < 4; ++i)
         for (int j = 0; j < 4; ++j) A[i][j] = (i == j) ? 1.0 : 0.0;
   }

   // Orthogonal symmetry transformation for operator Nop (0-based) + integer cell
   // shifts, from a gemmi cell + space group. The op acts in fractional space; we
   // conjugate it with the cell frac<->orth transforms so TMatrix maps orthogonal
   // coordinates directly (MMDB semantics). Returns 0 on success, 1 if there is no
   // usable space group / the operator is out of range. Shared by Manager and Cryst.
   inline int gemmi_sym_tmatrix(const gemmi::UnitCell &cell, const std::string &sg_name,
                                mat44 &TMatrix, int Nop, int a, int b, int c) {
      Mat4Init(TMatrix);
      const gemmi::SpaceGroup *sg = gemmi::find_spacegroup_by_name(sg_name);
      if (!sg || !cell.is_crystal()) return 1;
      gemmi::GroupOps gops = sg->operations();
      if (Nop < 0 || Nop >= (int)gops.order()) return 1;
      int i = 0;
      gemmi::Op op;
      for (gemmi::Op o : gops) {
         if (i++ == Nop) {
            op = o;
            break;
         }
      }
      gemmi::Transform sym{gemmi::rot_as_mat33(op),
                           gemmi::tran_as_vec3(op) + gemmi::Vec3(a, b, c)};
      gemmi::Transform t = cell.orth.combine(sym).combine(cell.frac);
      for (int r = 0; r < 3; ++r) {
         for (int cc = 0; cc < 3; ++cc) TMatrix[r][cc] = t.mat.a[r][cc];
         TMatrix[r][3] = t.vec.at(r);
      }
      return 0;
   }

   // Crystal/symmetry record (mmdb_cryst.h). Holds a gemmi cell + space-group name
   // and computes symmetry through the shared helper — same result as Manager for a
   // populated Cryst (Manager is the usual live symmetry path).
   class Cryst {
   public:
      gemmi::UnitCell cell;
      std::string spaceGroup;
      virtual ~Cryst() {}
      int GetTMatrix(mat44 &T, int Nop, int a, int b, int c) {
         return gemmi_sym_tmatrix(cell, spaceGroup, T, Nop, a, b, c);
      }
      int GetNumberOfSymOps() {
         const gemmi::SpaceGroup *sg = gemmi::find_spacegroup_by_name(spaceGroup);
         return sg ? (int)sg->operations().order() : 0;
      }
      pstr GetSymOp(int Nop) {
         const gemmi::SpaceGroup *sg = gemmi::find_spacegroup_by_name(spaceGroup);
         if (!sg) return nullptr;
         int i = 0;
         for (gemmi::Op op : sg->operations())
            if (i++ == Nop) {
               _symop_buf = op.triplet();
               return (pstr)_symop_buf.c_str();
            }
         return nullptr;
      }

   private:
      std::string _symop_buf;
   };
   typedef Cryst *PCryst;

   // mmdb::math graph-matching subsystem — full classes defined in _graph_impl.hh
   // (included at end of this file, after Atom/Residue are complete). Only the
   // Alignment class (unused by the cootapi build) stays a forward decl.
   namespace math {
      class Alignment;
   }

   struct AtomBond {
      PAtom atom = nullptr;
      int order = 0;
   };
   typedef AtomBond *PAtomBond;
   typedef AtomBond **PPAtomBond;

   struct AtomStat {  // selection coordinate statistics (mmdb_atom.h)
      int nAtoms = 0;
      realtype xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
      realtype xm = 0, ym = 0, zm = 0;  // coordinate means (centroid)
      realtype GetMaxSize() {
         realtype dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
         return dx > dy ? (dx > dz ? dx : dz) : (dy > dz ? dy : dz);
      }
   };
   typedef AtomStat &RAtomStat;

   // secondary-structure element codes (mmdb_tables.h)
   enum SSE_CODE { SSE_None = 0,
                   SSE_Strand = 1,
                   SSE_Bulge = 2,
                   SSE_3Turn = 3,
                   SSE_4Turn = 4,
                   SSE_5Turn = 5,
                   SSE_Helix = 6 };

   // PDBCleanup flags (mmdb_root.h) — bit flags OR'd into PDBCleanup(word)
   // misc return-code / sort-key enums (mmdb_cryst.h / mmdb_selmngr.h / mmdb_tables.h)
   enum { SYMOP_Ok = 0,
          SYMOP_NoLibFile = -1,
          SYMOP_UnknownSpaceGroup = -2 };
   enum { SSERC_Ok = 0,
          SSERC_noResidues = 1 };
   enum { SORT_CHAIN_ChainID_Asc = 0,
          SORT_CHAIN_ChainID_Desc = 1 };
   enum { CNSORT_OFF = 0,
          CNSORT_1INC = 1,
          CNSORT_1DEC = 2,
          CNSORT_2INC = 3,
          CNSORT_2DEC = 4 };

   enum PDB_CLEAN_FLAG {
      PDBCLEAN_ATNAME = 0x00000001,
      PDBCLEAN_TER = 0x00000002,
      PDBCLEAN_CHAIN = 0x00000004,
      PDBCLEAN_CHAIN_STRONG = 0x00000008,
      PDBCLEAN_ALTCODE = 0x00000010,
      PDBCLEAN_ALTCODE_STRONG = 0x00000020,
      PDBCLEAN_SERIAL = 0x00000040,
      PDBCLEAN_SEQNUM = 0x00000080,
      PDBCLEAN_INDEX = 0x00000800,
      PDBCLEAN_ELEMENT = 0x00001000,
      PDBCLEAN_ELEMENT_STRONG = 0x00002000
   };

   // SS records — public-member structs. Model::GetNumberOf{Helices,Sheets} are
   // populated from gemmi Structure::{helices,sheets} on load (_load_metadata) and
   // also fillable by Coot's own SS computation via the access_model subclass.
   class Helix {
   public:
      ChainID initChainID{}, endChainID{};
      int initSeqNum = 0, endSeqNum = 0, serNum = 0, helixClass = 0, length = 0;
      ResName initResName{}, endResName{};
      InsCode initICode{}, endICode{};
      char helixID[20]{}, comment[80]{};
   };
   class Strand {
   public:
      ChainID initChainID{}, endChainID{};
      int initSeqNum = 0, endSeqNum = 0, strandNo = 0, sense = 0;
      ResName initResName{}, endResName{};
      InsCode initICode{}, endICode{};
      char sheetID[20]{};
   };
   class Sheet {
   public:
      int nStrands = 0;
      Strand **strand = nullptr;
      char sheetID[20]{};
   };
   class Sheets {
   public:
      int nSheets = 0;
      Sheet **sheet = nullptr;
   };  // filled from gemmi in _load_metadata
   typedef Helix *PHelix;
   typedef Strand *PStrand;
   typedef Sheet *PSheet;
   typedef Sheets *PSheets;
   // container of helices (Model.helices); Coot's access_model subclass fills it.
   class Helices {
   public:
      std::vector<Helix *> data;
      void AddData(PHelix h) {
         if (h) data.push_back(h);
      }
      int nHelices = 0;
   };

   // container of symmetry operators (mmdb_symop.h SymOps). Coot fills it from a
   // space group; ops are xyz-triplet strings.
   class SymOps {
      std::vector<std::string> ops;
      std::deque<std::string> buf;

   public:
      int AddSymOp(cpstr xyz) {
         ops.push_back(xyz ? xyz : "");
         return 0;
      }
      int GetNofSymOps() { return (int)ops.size(); }
      pstr GetSymOp(int n) {
         if (n < 0 || n >= (int)ops.size()) return nullptr;
         buf.push_back(ops[n]);
         return (pstr)buf.back().c_str();
      }
      void FreeMemory() { ops.clear(); }
   };

   [[noreturn]] inline void unimpl(const char *w) {
      throw std::logic_error(std::string("mmdb-shim: unimplemented: ") + w);
   }

   // ---- free functions (mmdb_tables.h / mmdb_mattype.h) ----
   inline void InitMatType() {}  // real MMDB inits static matrix-type tables; no-op here
   inline cpstr GetErrorDescription(ERROR_CODE ec) {
      switch (ec) {
         case Error_NoError:
            return "no error";
         case Error_CantOpenFile:
            return "cannot open file";
         default:
            return "MMDB error";
      }
   }
   inline realtype getVdWaalsRadius(cpstr element) {
      return gemmi::Element(element ? element : "X").vdw_r();
   }

   // Borrowed empty C-string, returned by the delegating accessors when an object
   // is detached (no parent) — real MMDB yields safe defaults, not a crash.
   inline pstr mmdb_empty_pstr() {
      static char e[1] = {0};
      return e;
   }

   // UDData helpers (defined after Manager); each class forwards with its UDR type.
   int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, int v);
   int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, realtype v);
   int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, cpstr v);
   int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, int &v);
   int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, realtype &v);
   int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, pstr &v);

}  // namespace mmdb
