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

#include <cstdio>
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
inline const int      ANY_RES = -2147483647;   // real MMDB: extern const == MinInt4

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

struct Contact { int id1, id2; long group; realtype dist; };
typedef Contact  *PContact;

// LINK record. Public data members mirror real MMDB (Coot reads them directly).
// Not gemmi-backed yet — Model::GetNumberOfLinks currently returns 0 (TODO: map
// gemmi Structure connections), so these are declared for compilation.
class Link {
public:
  AtomName atName1{}, atName2{};
  AltLoc   aloc1{}, aloc2{};
  ResName  resName1{}, resName2{};
  ChainID  chainID1{}, chainID2{};
  InsCode  insCode1{}, insCode2{};
  int      seqNum1 = 0, seqNum2 = 0;
  realtype dist = 0;
};
typedef Link *PLink; typedef Link **PPLink;

[[noreturn]] inline void unimpl(const char *w) {
  throw std::logic_error(std::string("mmdb-shim: unimplemented: ") + w);
}

// UDData helpers (defined after Manager); each class forwards with its UDR type.
int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, int v);
int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, realtype v);
int ud_put(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, cpstr v);
int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, int &v);
int ud_get(Manager *mgr, UDR_TYPE myType, UDStore &s, int handle, realtype &v);

// ===========================================================================
class Atom : public UDStore {
public:
  Manager *mgr = nullptr;
  Residue *res = nullptr;   // parent
  int ai = 0;               // cached index within parent residue's atoms
  bool alive = true;

  gemmi::Atom &g() const;   // resolve to live gemmi (defined after Manager)

  // --- rewritten field accessors (pure B) ---
  // Scalar fields -> reference-returning accessors, so a uniform `->field`->
  // `->field()` rewrite covers both reads and writes. (occ/b_iso/charge are
  // narrower than realtype in gemmi, so those refs are float/schar-typed — the
  // rare take-address-of-realtype sites surface at Coot build time.)
  realtype &x() { return g().pos.x; }
  realtype &y() { return g().pos.y; }
  realtype &z() { return g().pos.z; }
  float &occupancy() { return g().occ; }
  float &tempFactor() { return g().b_iso; }
  char  &altLoc()    { return g().altloc; }
  signed char &charge() { return g().charge; }
  int   &serNum()    { return g().serial; }
  void set_occupancy(realtype v) { g().occ = (float)v; }
  void set_tempFactor(realtype v) { g().b_iso = (float)v; }
  void set_altLoc(char c) { g().altloc = c; }

  // --- method surface (hot subset; rest stubbed) ---
  pstr GetAtomName();                 // aligned name, MMDB semantics
  void SetAtomName(const AtomName aName);
  pstr GetElementName();
  void SetElementName(const Element elName);
  cpstr GetChainID();
  int   GetSeqNum();
  pstr  GetInsCode();
  pstr  GetResName();
  Residue *GetResidue() { return res; }
  int   GetModelNum();
  bool  isTer() const { return false; } // gemmi has no TER atoms; see notes
  void  SetCoordinates(realtype xx, realtype yy, realtype zz,
                       realtype occ, realtype tF);
  int   GetIndex();
  // UDData
  int PutUDData(int h, int v)      { return ud_put(mgr, UDR_ATOM, *this, h, v); }
  int PutUDData(int h, realtype v) { return ud_put(mgr, UDR_ATOM, *this, h, v); }
  int PutUDData(int h, cpstr v)    { return ud_put(mgr, UDR_ATOM, *this, h, v); }
  int GetUDData(int h, int &v)      { return ud_get(mgr, UDR_ATOM, *this, h, v); }
  int GetUDData(int h, realtype &v) { return ud_get(mgr, UDR_ATOM, *this, h, v); }

private:
  AtomName _name_buf{}; Element _elem_buf{};
};

// ===========================================================================
class Residue : public UDStore {
public:
  Manager *mgr = nullptr;
  Chain *chain = nullptr;   // parent
  int ri = 0;
  bool alive = true;
  std::vector<Atom *> atoms;  // canonical child wrappers == PPAtom table

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
  int   AddAtom(PAtom /*atm*/) { unimpl("Residue::AddAtom(PAtom)"); }
  void  DeleteAtom(int pos);

  pstr  GetResName();
  void  SetResName(const ResName n) { g().name = n; }
  int   GetSeqNum();
  pstr  GetInsCode();
  cpstr GetChainID();
  int   GetModelNum();
  int   GetIndex() { return ri; }
  Chain   *GetChain() { return chain; }
  Residue *next = nullptr; // MMDB has this; wired lazily if needed
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
  Model *model = nullptr;   // parent
  int ci = 0;
  bool alive = true;
  std::vector<Residue *> residues;

  gemmi::Chain &g() const;

  int      GetNumberOfResidues() { return (int)residues.size(); }
  PResidue GetResidue(int resNo) {
    return (resNo >= 0 && resNo < (int)residues.size()) ? residues[resNo] : nullptr;
  }
  void     GetResidueTable(PPResidue &t, int &n) { t = residues.data(); n = (int)residues.size(); }
  cpstr    GetChainID();
  PResidue AddResidue(Manager &m, gemmi::Residue r);       // append
  PResidue InsResidue(Manager &m, int pos, gemmi::Residue r);

private:
  ChainID _chainid_buf{};
};

// ===========================================================================
class Model : public UDStore {
public:
  Manager *mgr = nullptr;
  int mi = 0;               // 0-based internal; GetModel is 1-based externally
  std::vector<Chain *> chains;

  gemmi::Model &g() const;

  int    GetNumberOfChains() { return (int)chains.size(); }
  PChain GetChain(int chainNo) {
    return (chainNo >= 0 && chainNo < (int)chains.size()) ? chains[chainNo] : nullptr;
  }
  PChain GetChain(const ChainID chID);
  int    GetSerNum() { return mi + 1; }
  // LINK records — TODO: map from gemmi Structure connections.
  int    GetNumberOfLinks() { return 0; }
  PLink  GetLink(int /*i*/) { return nullptr; }
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

  void build_from_gemmi();

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

  // ---- contacts (gemmi-free uniform-grid search) ----
  void SeekContacts(PPAtom A1, int n1, PPAtom A2, int n2, realtype d1,
                    realtype d2, int seqDist, PContact &contact, int &ncontacts,
                    int maxlen = 0, long group = 0);
  void SeekContacts(PPAtom A, int n, realtype d1, realtype d2, int seqDist,
                    PContact &contact, int &ncontacts, int maxlen = 0, long group = 0);

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
inline gemmi::Model   &Model::g()   const { return mgr->st.models[mi]; }
inline gemmi::Chain   &Chain::g()   const { return model->g().chains[ci]; }
inline gemmi::Residue &Residue::g() const { return chain->g().residues[ri]; }
inline gemmi::Atom    &Atom::g()    const { return res->g().atoms[ai]; }

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

// ---- Atom out-of-line ----
inline pstr Atom::GetAtomName() {
  std::snprintf(_name_buf, sizeof(_name_buf), "%s", g().name.c_str());
  return _name_buf;
}
inline void Atom::SetAtomName(const AtomName aName) { g().name = aName; }
inline pstr Atom::GetElementName() {
  std::snprintf(_elem_buf, sizeof(_elem_buf), "%s", g().element.name());
  return _elem_buf;
}
inline void Atom::SetElementName(const Element elName) { g().element = gemmi::Element(elName); }
inline cpstr Atom::GetChainID()  { return res->GetChainID(); }
inline int   Atom::GetSeqNum()   { return res->GetSeqNum(); }
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
inline int  Residue::GetSeqNum()  { return g().seqid.num.value; }
inline pstr Residue::GetInsCode() {
  _inscode_buf[0] = g().seqid.icode == ' ' ? '\0' : g().seqid.icode; _inscode_buf[1] = '\0';
  return _inscode_buf;
}
inline cpstr Residue::GetChainID() { return chain->GetChainID(); }
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
inline cpstr Chain::GetChainID() {
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
          rw->atoms.push_back(aw);
        }
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

} // namespace mmdb
