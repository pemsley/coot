// -*- mode: c++; -*-
//
// mmdb-shim: mmdb::math::{Vertex,Edge,Graph,GraphMatch} — molecular graph +
// subgraph matching. gemmi has NO subgraph-isomorphism facility (see
// mmdb-graph-matching-for-gemmi.md §3), so per coot-shim-prefer-gemmi this is one
// of the few pieces genuinely hand-written: a branch-and-bound maximum-common-
// (induced-)subgraph matcher. Element typing uses gemmi::Element.
//
// Copyright 2026 by Medical Research Council Laboratory of Molecular Biology
//
// Included at the end of _shim_impl.hh (Residue/Atom complete; pstr/cpstr/
// realtype/ivector/word already in namespace mmdb).
#ifndef COOT_MMDB_SHIM_GRAPH_IMPL_HH
#define COOT_MMDB_SHIM_GRAPH_IMPL_HH

#include <gemmi/elem.hpp>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

namespace mmdb {
namespace math {

// ---- constants (mmdb_math_graph.h) --------------------------------------
enum GRAPH_BOND { BOND_SINGLE = 1, BOND_DOUBLE = 2, BOND_AROMATIC = 3, BOND_TRIPLE = 4 };
enum GRAPH_RC   { MKGRAPH_Ok = 0, MKGRAPH_NoAtoms = -1,
                  MKGRAPH_ChangedAltLoc = 1, MKGRAPH_MaxOccupancy = 2 };
enum GRAPH_MATCH_FLAG { GMF_UniqueMatch = 0x00000001, GMF_NoCombinations = 0x00000002 };
enum VERTEX_EXT_TYPE  { EXTTYPE_Ignore = 0, EXTTYPE_Equal = 1, EXTTYPE_AND = 2,
                        EXTTYPE_OR = 3, EXTTYPE_XOR = 4, EXTTYPE_NotEqual = 5,
                        EXTTYPE_NotAND = 6, EXTTYPE_NotOR = 7 };

namespace gdetail {
  inline std::string trim(cpstr s) {
    std::string t(s ? s : "");
    size_t a = t.find_first_not_of(" \t"), b = t.find_last_not_of(" \t");
    return a == std::string::npos ? std::string() : t.substr(a, b - a + 1);
  }
  inline int bond_from_string(cpstr s) {
    std::string t = trim(s);
    for (auto &c : t) c = (char) std::tolower((unsigned char) c);
    if (t == "single" || t == "sing" || t == "1") return BOND_SINGLE;
    if (t == "double" || t == "doub" || t == "2") return BOND_DOUBLE;
    if (t == "aromatic" || t == "arom" || t == "ar") return BOND_AROMATIC;
    if (t == "triple" || t == "trip" || t == "3") return BOND_TRIPLE;
    return BOND_SINGLE;
  }
}

// =========================================================================
//  Vertex — a graph node (atom). `type` encodes element (atomic number) so
//  GetType() equality means "same element".
// =========================================================================
class Vertex {
 public:
  int type = 0, type_ext = 0, property = 0, nBonds = 0, id = 0, user_id = 0;
  std::string name;

  Vertex() {}
  Vertex(cpstr chem_elem) { SetVertex(chem_elem); }
  Vertex(cpstr chem_elem, cpstr vname) { SetVertex(chem_elem); name = vname ? vname : ""; }
  Vertex(int vtype, cpstr vname) { type = vtype; name = vname ? vname : ""; }
  explicit Vertex(int vtype) { type = vtype; }

  void  SetVertex(cpstr chem_elem) {
    type = (int) gemmi::Element(gdetail::trim(chem_elem).c_str()).atomic_number();
  }
  void  SetVertex(int vtype, cpstr vname) { type = vtype; name = vname ? vname : ""; }
  void  SetVertex(int vtype) { type = vtype; }
  void  SetName(cpstr vname) { name = vname ? vname : ""; }
  void  SetType(int t) { type = t; }
  void  SetTypeExt(int t) { type_ext = t; }
  void  SaveType() { property = type; }
  void  RestoreType() { type = property; }
  void  SetUserID(int u) { user_id = u; }
  cpstr GetName()    { return name.c_str(); }
  int   GetType()    { return type; }
  int   GetTypeExt() { return type_ext; }
  int   GetNBonds()  { return nBonds; }
  int   GetUserID()  { return user_id; }     // 0-based atom index (see MakeVertexIDs/MakeGraph)
  void  Print(int /*PKey*/ = 0) {}
};
typedef Vertex *PVertex; typedef Vertex **PPVertex;

// =========================================================================
//  Edge — a graph connection (bond). v1/v2 are 1-indexed vertex numbers.
// =========================================================================
class Edge {
 public:
  int v1 = 0, v2 = 0, type = 0, property = 0;

  Edge() {}
  Edge(int vx1, int vx2, int btype) { v1 = vx1; v2 = vx2; type = btype; }
  Edge(int vx1, int vx2, cpstr btype) { v1 = vx1; v2 = vx2; type = gdetail::bond_from_string(btype); }

  void SetEdge(int vx1, int vx2, int btype) { v1 = vx1; v2 = vx2; type = btype; }
  void SetEdge(int vx1, int vx2, cpstr btype) { v1 = vx1; v2 = vx2; type = gdetail::bond_from_string(btype); }
  void SetType(int t) { type = t; }
  int  GetVertex1() { return v1; }
  int  GetVertex2() { return v2; }
  int  GetType()    { return type; }
  void Print(int /*PKey*/ = 0) {}
};
typedef Edge *PEdge; typedef Edge **PPEdge;

// =========================================================================
//  Graph — vertices + edges + adjacency matrix (built by Build()).
//  Owns the Vertex/Edge objects handed to AddVertex/AddEdge (mmdb semantics).
// =========================================================================
class Graph {
 public:
  std::string gname;
  std::vector<Vertex *> V;                 // owned; 1-indexed via GetVertex
  std::vector<Edge *>   E;                 // owned
  std::vector<std::vector<int>> adj;       // (n+1)x(n+1), 1-indexed; adj[i][j]=bond type or 0

  Graph() {}
  ~Graph() { for (auto *p : V) delete p; for (auto *e : E) delete e; }
  Graph(const Graph &) = delete;
  Graph &operator=(const Graph &) = delete;

  void  SetName(cpstr n) { gname = n ? n : ""; }
  pstr  GetName() { return (pstr) gname.c_str(); }
  void  AddVertex(PVertex v) { if (v) V.push_back(v); }
  void  AddEdge(PEdge e)     { if (e) E.push_back(e); }
  int   GetNofVertices() { return (int) V.size(); }
  int   GetNofEdges()    { return (int) E.size(); }
  PVertex GetVertex(int i) { return (i >= 1 && i <= (int) V.size()) ? V[i - 1] : nullptr; }
  PEdge   GetEdge(int i)   { return (i >= 1 && i <= (int) E.size()) ? E[i - 1] : nullptr; }
  void  GetVertices(PPVertex &v, int &n) { v = V.data(); n = (int) V.size(); }
  void  GetEdges(PPEdge &e, int &n)      { e = E.data(); n = (int) E.size(); }
  // number vertices; user_id = 0-based position so `residue->atom[V->GetUserID()]`
  // (Coot's manual make_graph path) indexes the residue's 0-based atom table.
  void  MakeVertexIDs() { for (int i = 0; i < (int) V.size(); ++i) { V[i]->id = i + 1; V[i]->user_id = i; } }
  void  Print() {}
  void  Print1() {}
  void  MakeSymmetryRelief(bool /*noCO2*/) {}   // type_ext modifiers — low priority, no-op
  void  IdentifyRings() {}                        // "     "
  void  IdentifyConnectedComponents() {}

  // adjacency matrix; bondOrder=false collapses all bonds to 1 (connectivity only)
  int   Build(bool bondOrder) {
    int n = (int) V.size();
    adj.assign(n + 1, std::vector<int>(n + 1, 0));
    for (Edge *e : E)
      if (e->v1 >= 1 && e->v1 <= n && e->v2 >= 1 && e->v2 <= n) {
        int t = bondOrder ? (e->type > 0 ? e->type : 1) : 1;
        adj[e->v1][e->v2] = t; adj[e->v2][e->v1] = t;
      }
    for (int i = 1; i <= n; ++i) {
      int b = 0; for (int j = 1; j <= n; ++j) if (adj[i][j]) ++b;
      V[i - 1]->nBonds = b;
    }
    return 0;
  }

  int MakeGraph(PPAtom atom, int nAtoms);          // build from atoms (distance bonds)
  int MakeGraph(PResidue R, cpstr altLoc = nullptr);
};
typedef Graph *PGraph;

// =========================================================================
//  GraphMatch — maximum common (induced) subgraph via branch-and-bound.
// =========================================================================
class GraphMatch {
 public:
  struct Match { std::vector<int> f1, f2; };   // 1-indexed vertex numbers
  std::vector<Match> matches;
  int  maxMatch = 0;
  int  maxNofMatches = 1000000;
  bool stopOnMax = false;
  int  timeLimit = 0;                           // seconds; 0 = no limit
  word flags = 0;
  bool Stop = false;
  // stable 1-indexed ivector storage for GetMatch (freed with this object)
  std::vector<std::vector<int>> fv1, fv2;

  void SetFlag(word f)   { flags |= f; }
  void RemoveFlag(word f) { flags &= ~f; }
  void SetMaxNofMatches(int m, bool stopOnMaxN) { maxNofMatches = m > 0 ? m : 1; stopOnMax = stopOnMaxN; }
  void SetTimeLimit(int t = 0) { timeLimit = t; }
  int  GetNofMatches()  { return (int) matches.size(); }
  int  GetMaxMatchSize() { return maxMatch; }
  bool GetStopSignal()  { return Stop; }
  void Reset() { matches.clear(); fv1.clear(); fv2.clear(); maxMatch = 0; }
  void PrintMatches() {}

  void MatchGraphs(PGraph Gh1, PGraph Gh2, int minMatch, bool vertexType = true,
                   VERTEX_EXT_TYPE vertexExt = EXTTYPE_Ignore);
  void GetMatch(int MatchNo, ivector &FV1, ivector &FV2, int &nv, realtype &p1, realtype &p2);
};
typedef GraphMatch *PGraphMatch;

// ---- MatchGraphs: branch-and-bound maximum common induced subgraph -------
inline void GraphMatch::MatchGraphs(PGraph g1, PGraph g2, int minMatch,
                                    bool vertexType, VERTEX_EXT_TYPE vertexExt) {
  Reset(); Stop = false;
  int n1 = g1->GetNofVertices(), n2 = g2->GetNofVertices();
  if (n1 == 0 || n2 == 0) return;
  if ((int) g1->adj.size() != n1 + 1) g1->Build(false);
  if ((int) g2->adj.size() != n2 + 1) g2->Build(false);
  const std::vector<std::vector<int>> &A1 = g1->adj, &A2 = g2->adj;

  auto t0 = std::chrono::steady_clock::now();
  auto timed_out = [&]() -> bool {
    if (timeLimit <= 0) return false;
    if (std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - t0).count() >= timeLimit) { Stop = true; return true; }
    return false;
  };

  std::vector<std::pair<int, int>> cur;         // matched (g1,g2) 1-indexed pairs
  std::vector<char> used2(n2 + 1, 0);
  std::vector<std::vector<std::pair<int, int>>> best;

  auto compatible = [&](int i1, int i2) -> bool {
    if (vertexType && g1->V[i1 - 1]->type != g2->V[i2 - 1]->type) return false;
    if (vertexExt == EXTTYPE_Equal && g1->V[i1 - 1]->type_ext != g2->V[i2 - 1]->type_ext) return false;
    for (const auto &pr : cur)                  // preserve edges (incl. absence) + bond type
      if (A1[i1][pr.first] != A2[i2][pr.second]) return false;
    return true;
  };

  std::function<void(int)> rec = [&](int idx) {
    if (Stop || timed_out()) return;
    // bound: best achievable from here can't exceed cur + remaining g1 vertices
    if ((int) cur.size() + (n1 - idx + 1) < maxMatch) return;
    if (idx > n1) {
      int sz = (int) cur.size();
      if (sz >= minMatch && sz >= maxMatch) {
        if (sz > maxMatch) { maxMatch = sz; best.clear(); }
        if ((int) best.size() < maxNofMatches) best.push_back(cur);
        else if (stopOnMax) Stop = true;
      }
      return;
    }
    for (int i2 = 1; i2 <= n2 && !Stop; ++i2) {   // Option A: map g1[idx]
      if (used2[i2] || !compatible(idx, i2)) continue;
      cur.push_back({idx, i2}); used2[i2] = 1;
      rec(idx + 1);
      used2[i2] = 0; cur.pop_back();
    }
    if (!Stop) rec(idx + 1);                      // Option B: skip g1[idx]
  };
  rec(1);

  for (auto &m : best) {
    Match mm;
    for (auto &pr : m) { mm.f1.push_back(pr.first); mm.f2.push_back(pr.second); }
    matches.push_back(std::move(mm));
  }
}

inline void GraphMatch::GetMatch(int MatchNo, ivector &FV1, ivector &FV2,
                                 int &nv, realtype &p1, realtype &p2) {
  if (MatchNo < 0 || MatchNo >= (int) matches.size()) { FV1 = FV2 = nullptr; nv = 0; p1 = p2 = 0; return; }
  if ((int) fv1.size() != (int) matches.size()) { fv1.resize(matches.size()); fv2.resize(matches.size()); }
  Match &m = matches[MatchNo];
  nv = (int) m.f1.size();
  std::vector<int> &a = fv1[MatchNo]; std::vector<int> &b = fv2[MatchNo];
  a.assign(nv + 1, 0); b.assign(nv + 1, 0);      // 1-indexed: [1..nv] valid
  for (int i = 0; i < nv; ++i) { a[i + 1] = m.f1[i]; b[i + 1] = m.f2[i]; }
  FV1 = a.data(); FV2 = b.data();
  p1 = p2 = (realtype) nv;
}

// ---- MakeGraph from atoms (distance-based bonds) -------------------------
inline int Graph::MakeGraph(PPAtom atom, int nAtoms) {
  for (auto *p : V) delete p; for (auto *e : E) delete e;
  V.clear(); E.clear(); adj.clear();
  if (nAtoms <= 0) return MKGRAPH_NoAtoms;
  for (int i = 0; i < nAtoms; ++i) {
    Vertex *v = new Vertex(atom[i]->GetElementName(), atom[i]->GetAtomName());
    v->user_id = i;                        // 0-based atom index for atom-table lookup
    V.push_back(v);
  }
  // covalent-ish bonds by distance (< 1.9 A) — gemmi has no residue bond table
  for (int i = 0; i < nAtoms; ++i)
    for (int j = i + 1; j < nAtoms; ++j) {
      double dx = atom[i]->x() - atom[j]->x();
      double dy = atom[i]->y() - atom[j]->y();
      double dz = atom[i]->z() - atom[j]->z();
      if (dx * dx + dy * dy + dz * dz < 1.9 * 1.9)
        E.push_back(new Edge(i + 1, j + 1, BOND_SINGLE));
    }
  return MKGRAPH_Ok;
}
inline int Graph::MakeGraph(PResidue R, cpstr /*altLoc*/) {
  if (!R) return MKGRAPH_NoAtoms;
  PPAtom a = nullptr; int n = 0;
  R->GetAtomTable(a, n);
  return MakeGraph(a, n);
}

// =========================================================================
//  Alignment — global (Needleman-Wunsch) sequence alignment (mmdb_math_align.h).
//  Coot uses it to align a model sequence to a target and read back the gapped
//  strings + score (for mutation/indel detection, api/coot-molecule.cc).
//
//  NOTE on prefer-gemmi: gemmi DOES provide sequence alignment (align.hpp /
//  seqalign.hpp), but only over a substitution-scoring matrix. MMDB's aligner —
//  which Coot's mutation detection was tuned against — uses simple identity
//  scoring (match=1 / mismatch=0 / linear gap). Swapping in gemmi's scorer would
//  change which residues are called mutations vs. indels, i.e. it would NOT
//  reproduce the MMDB baseline this shim exists to preserve. So this stays a
//  faithful re-creation of MMDB's *specific* simple-scoring global aligner —
//  the exact form gemmi does not offer — rather than a gratuitous reimplementation.
// =========================================================================
class Alignment {
 public:
  std::string _s, _t;              // aligned (gapped) sequences
  realtype VAchieved = 0;          // alignment score

  void Align(cpstr S, cpstr T, realtype /*VGap*/ = 0.0, realtype VSpace = -1.0) {
    std::string a(S ? S : ""), b(T ? T : "");
    int n = (int) a.size(), m = (int) b.size();
    const realtype MATCH = 1.0, MIS = 0.0, GAP = VSpace != 0.0 ? VSpace : -1.0;
    std::vector<std::vector<realtype>> D(n + 1, std::vector<realtype>(m + 1, 0.0));
    for (int i = 1; i <= n; ++i) D[i][0] = D[i - 1][0] + GAP;
    for (int j = 1; j <= m; ++j) D[0][j] = D[0][j - 1] + GAP;
    for (int i = 1; i <= n; ++i)
      for (int j = 1; j <= m; ++j) {
        realtype diag = D[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? MATCH : MIS);
        realtype up = D[i - 1][j] + GAP, left = D[i][j - 1] + GAP;
        D[i][j] = std::max(diag, std::max(up, left));
      }
    VAchieved = D[n][m];
    std::string as, bs;                                   // traceback
    int i = n, j = m;
    while (i > 0 || j > 0) {
      if (i > 0 && j > 0 &&
          D[i][j] == D[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? MATCH : MIS)) {
        as += a[i - 1]; bs += b[j - 1]; --i; --j;
      } else if (i > 0 && D[i][j] == D[i - 1][j] + GAP) {
        as += a[i - 1]; bs += '-'; --i;
      } else {
        as += '-'; bs += b[j - 1]; --j;
      }
    }
    std::reverse(as.begin(), as.end()); std::reverse(bs.begin(), bs.end());
    _s = as; _t = bs;
  }
  pstr     GetAlignedS() { return (pstr) _s.c_str(); }
  pstr     GetAlignedT() { return (pstr) _t.c_str(); }
  realtype GetScore()    { return VAchieved; }
};

}  // namespace math
}  // namespace mmdb

#endif  // COOT_MMDB_SHIM_GRAPH_IMPL_HH
