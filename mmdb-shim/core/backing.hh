// mmdb-shim core — architecture B foundation.
// See ../../MMDB_SHIM_Recon_and_Plan.md.
//
// A LIVE gemmi::Structure holds the data; a PARALLEL wrapper tree provides the
// MMDB pointer semantics Coot depends on:
//   * wrapper addresses are stable (pool-allocated) -> valid `mmdb::Atom*` and
//     usable as std::set/map keys (identity by address);
//   * GetAtom(i)/GetResidue(i) return the SAME canonical wrapper pointer every
//     call (identity cache = the parent's child-pointer vector);
//   * g() resolves a wrapper to its live gemmi object via parent* + a cached
//     sibling index;
//   * structural edits patch ONLY the shifted siblings in one container
//     (localized) — not the whole pool, unlike the first spike.
//
// This is the hardened successor to mmdb-recon/spike/.

#pragma once
#include <gemmi/model.hpp>

#include <deque>
#include <string>
#include <vector>

namespace shim {

struct Backing;
struct ModelW;
struct ChainW;
struct ResidueW;

// ---------------------------------------------------------------------------
struct AtomW {
  ResidueW *parent = nullptr;
  int ai = 0;               // index within parent residue's atoms (cached)
  bool alive = true;

  gemmi::Atom &g() const;   // resolve to live gemmi::Atom (defined below)

  // numeric fields -> reference-returning accessors (rewrite targets)
  double &x() { return g().pos.x; }
  double &y() { return g().pos.y; }
  double &z() { return g().pos.z; }
  // occ/b_iso are float in gemmi -> value get + setter (can't bind double&)
  double occupancy() const { return g().occ; }
  void set_occupancy(double v) { g().occ = (float)v; }
  double tempFactor() const { return g().b_iso; }
  void set_tempFactor(double v) { g().b_iso = (float)v; }
  // char-array fields -> const char* getter + setter (strcpy sites rewrite here)
  const char *name() const { return g().name.c_str(); }
  void set_name(const char *s) { g().name = s; }
};

// ---------------------------------------------------------------------------
struct ResidueW {
  ChainW *parent = nullptr;
  int ri = 0;                       // index within parent chain's residues
  bool alive = true;
  std::vector<AtomW *> atoms_w;     // canonical child wrappers (identity cache)

  gemmi::Residue &g() const;

  int GetNumberOfAtoms() const { return (int)atoms_w.size(); }
  AtomW *GetAtom(int i) { return (i >= 0 && i < (int)atoms_w.size()) ? atoms_w[i] : nullptr; }

  AtomW *AddAtom(Backing &b, gemmi::Atom a);      // append: O(1), no sibling shift
  void DeleteAtom(int pos);                        // O(atoms in this residue)
};

// ---------------------------------------------------------------------------
struct ChainW {
  ModelW *parent = nullptr;
  int ci = 0;
  bool alive = true;
  std::vector<ResidueW *> residues_w;

  gemmi::Chain &g() const;

  int GetNumberOfResidues() const { return (int)residues_w.size(); }
  ResidueW *GetResidue(int i) { return (i >= 0 && i < (int)residues_w.size()) ? residues_w[i] : nullptr; }

  ResidueW *AddResidue(Backing &b, gemmi::Residue r);     // append: O(1)
  ResidueW *InsResidue(Backing &b, int pos, gemmi::Residue r); // O(residues in chain)
};

// ---------------------------------------------------------------------------
struct ModelW {
  Backing *b = nullptr;
  int mi = 0;
  std::vector<ChainW *> chains_w;

  gemmi::Model &g() const;

  int GetNumberOfChains() const { return (int)chains_w.size(); }
  ChainW *GetChain(int i) { return (i >= 0 && i < (int)chains_w.size()) ? chains_w[i] : nullptr; }
};

// ---------------------------------------------------------------------------
struct Backing {
  gemmi::Structure st;
  // Pools: std::deque keeps element addresses stable across growth.
  std::deque<AtomW> atom_pool;
  std::deque<ResidueW> res_pool;
  std::deque<ChainW> chain_pool;
  std::deque<ModelW> model_pool;
  std::vector<ModelW *> models_w;

  AtomW *newAtom() { atom_pool.emplace_back(); return &atom_pool.back(); }
  ResidueW *newRes() { res_pool.emplace_back(); return &res_pool.back(); }
  ChainW *newChain() { chain_pool.emplace_back(); return &chain_pool.back(); }
  ModelW *newModel() { model_pool.emplace_back(); return &model_pool.back(); }

  ModelW *GetModel(int i) { return (i >= 0 && i < (int)models_w.size()) ? models_w[i] : nullptr; }
  int GetNumberOfModels() const { return (int)models_w.size(); }

  // Build the parallel wrapper tree from the current gemmi::Structure.
  void build_from_gemmi() {
    models_w.clear();
    for (int mi = 0; mi < (int)st.models.size(); ++mi) {
      ModelW *mw = newModel(); mw->b = this; mw->mi = mi;
      auto &gm = st.models[mi];
      for (int ci = 0; ci < (int)gm.chains.size(); ++ci) {
        ChainW *cw = newChain(); cw->parent = mw; cw->ci = ci;
        auto &gc = gm.chains[ci];
        for (int ri = 0; ri < (int)gc.residues.size(); ++ri) {
          ResidueW *rw = newRes(); rw->parent = cw; rw->ri = ri;
          auto &gr = gc.residues[ri];
          for (int ai = 0; ai < (int)gr.atoms.size(); ++ai) {
            AtomW *aw = newAtom(); aw->parent = rw; aw->ai = ai;
            rw->atoms_w.push_back(aw);
          }
          cw->residues_w.push_back(rw);
        }
        mw->chains_w.push_back(cw);
      }
      models_w.push_back(mw);
    }
  }
};

// ---- g() resolvers (walk parent + cached index into the live gemmi tree) ----
inline gemmi::Model &ModelW::g() const { return b->st.models[mi]; }
inline gemmi::Chain &ChainW::g() const { return parent->g().chains[ci]; }
inline gemmi::Residue &ResidueW::g() const { return parent->g().residues[ri]; }
inline gemmi::Atom &AtomW::g() const { return parent->g().atoms[ai]; }

// ---- edits (localized patching) ----
inline AtomW *ResidueW::AddAtom(Backing &b, gemmi::Atom a) {
  g().atoms.push_back(std::move(a));               // append -> existing ai valid
  AtomW *aw = b.newAtom();
  aw->parent = this; aw->ai = (int)atoms_w.size();
  atoms_w.push_back(aw);
  return aw;
}

inline void ResidueW::DeleteAtom(int pos) {
  if (pos < 0 || pos >= (int)atoms_w.size()) return;
  g().atoms.erase(g().atoms.begin() + pos);
  atoms_w[pos]->alive = false; atoms_w[pos]->ai = -1;   // tombstone (reclaim later)
  atoms_w.erase(atoms_w.begin() + pos);
  for (int k = pos; k < (int)atoms_w.size(); ++k) atoms_w[k]->ai = k;  // shift -1
}

inline ResidueW *ChainW::AddResidue(Backing &b, gemmi::Residue r) {
  g().residues.push_back(std::move(r));
  ResidueW *rw = b.newRes();
  rw->parent = this; rw->ri = (int)residues_w.size();
  for (int ai = 0; ai < (int)rw->g().atoms.size(); ++ai) {
    AtomW *aw = b.newAtom(); aw->parent = rw; aw->ai = ai; rw->atoms_w.push_back(aw);
  }
  residues_w.push_back(rw);
  return rw;
}

inline ResidueW *ChainW::InsResidue(Backing &b, int pos, gemmi::Residue r) {
  g().residues.insert(g().residues.begin() + pos, std::move(r));
  ResidueW *rw = b.newRes();
  rw->parent = this; rw->ri = pos;
  residues_w.insert(residues_w.begin() + pos, rw);
  for (int k = pos + 1; k < (int)residues_w.size(); ++k) residues_w[k]->ri = k;  // shift +1
  for (int ai = 0; ai < (int)rw->g().atoms.size(); ++ai) {
    AtomW *aw = b.newAtom(); aw->parent = rw; aw->ai = ai; rw->atoms_w.push_back(aw);
  }
  return rw;  // atoms of OTHER residues are untouched (parent ri updated, ai unchanged)
}

} // namespace shim
