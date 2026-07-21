// mmdb-shim contact search — Manager::SeekContacts via gemmi::NeighborSearch
// (grid-accelerated; handles crystallographic PBC and non-crystal bounding box).
// Compiled separately so gemmi/neighbor.hpp stays out of Coot's many TUs.
#define COOT_USE_MMDB_SHIM 1
#include <mmdb2/mmdb_manager.h>

#include <gemmi/neighbor.hpp>

#include <cmath>
#include <unordered_map>

namespace mmdb {
namespace {

bool seqNeglect(Atom *a, Atom *b, int seqDist) {
  if (seqDist <= 0) return false;
  if (a->res->chain != b->res->chain) return false;
  return std::abs(a->GetSeqNum() - b->GetSeqNum()) < seqDist;
}
void alloc_contacts(std::vector<Contact> &v, PContact &out, int &n) {
  n = (int)v.size();
  out = n ? new Contact[n] : nullptr;          // caller delete[]s (MMDB semantics)
  for (int i = 0; i < n; ++i) out[i] = v[i];
}
// Map a NeighborSearch Mark back to the shim wrapper via the parallel tree.
inline Atom *mark_to_atom(Model *mw, const gemmi::NeighborSearch::Mark *m) {
  return mw->chains[m->chain_idx]->residues[m->residue_idx]->atoms[m->atom_idx];
}

template <class Vec>
void skcombine(Vec &cur, Vec &m, SELECTION_KEY k) {
  std::set<typename Vec::value_type> cs(cur.begin(), cur.end()), ms(m.begin(), m.end());
  if (k == SKEY_NEW) cur = m;
  else if (k == SKEY_OR) { for (auto *x : m) if (!cs.count(x)) cur.push_back(x); }
  else if (k == SKEY_AND) { Vec o; for (auto *x : cur) if (ms.count(x)) o.push_back(x); cur = o; }
  else if (k == SKEY_XOR) { Vec o; for (auto *x : cur) if (!ms.count(x)) o.push_back(x);
                            for (auto *x : m) if (!cs.count(x)) o.push_back(x); cur = o; }
  else if (k == SKEY_CLR) { Vec o; for (auto *x : cur) if (!ms.count(x)) o.push_back(x); cur = o; }
}

} // namespace

// atoms within [d1,d2] of any atom in the given set.
void Manager::SelectNeighbours(int selHnd, SELECTION_TYPE sType, PPAtom atoms,
    int nAtoms, realtype d1, realtype d2, SELECTION_KEY sKey) {
  std::vector<Atom *> mAtoms; std::vector<Residue *> mResidues;
  if (nAtoms > 0) {
    Model *mw = atoms[0]->res->chain->model;
    gemmi::NeighborSearch ns(mw->g(), st.cell, d2);
    ns.populate(true);
    std::set<Atom *> seenA; std::set<Residue *> seenR;
    for (int i = 0; i < nAtoms; ++i)
      for (auto *m : ns.find_atoms(atoms[i]->g().pos, '\0', d1, d2)) {
        if (m->image_idx != 0) continue;
        Atom *b = mark_to_atom(mw, m);
        if (sType == STYPE_ATOM) { if (seenA.insert(b).second) mAtoms.push_back(b); }
        else if (sType == STYPE_RESIDUE) { if (seenR.insert(b->res).second) mResidues.push_back(b->res); }
      }
  }
  Selection &sel = selections[selHnd - 1];
  if (sel.type == STYPE_UNDEFINED) sel.type = sType;
  std::vector<Atom *> oldA = sel.atoms; std::vector<Residue *> oldR = sel.residues;
  if (sType == STYPE_ATOM) skcombine(sel.atoms, mAtoms, sKey);
  else if (sType == STYPE_RESIDUE) skcombine(sel.residues, mResidues, sKey);
  for (Atom *a : oldA) a->_setInSel(selHnd, false);
  for (Atom *a : sel.atoms) a->_setInSel(selHnd, true);
  for (Residue *r : oldR) r->_setInSel(selHnd, false);
  for (Residue *r : sel.residues) r->_setInSel(selHnd, true);
}

void Manager::SeekContacts(PPAtom A1, int n1, PPAtom A2, int n2, realtype d1,
    realtype d2, int seqDist, PContact &contact, int &ncontacts, int /*maxlen*/,
    pmat44 /*TMatrix*/, long group) {
  std::vector<Contact> found;
  if (n1 > 0 && n2 > 0) {
    Model *mw = A1[0]->res->chain->model;   // NeighborSearch is per-model
    gemmi::NeighborSearch ns(mw->g(), st.cell, d2);
    ns.populate(true);
    std::unordered_map<Atom *, int> a2idx;
    a2idx.reserve(n2 * 2);
    for (int j = 0; j < n2; ++j) a2idx.emplace(A2[j], j);

    for (int i = 0; i < n1; ++i) {
      if (A1[i]->res->chain->model != mw) continue;  // single-model contact search
      for (auto *m : ns.find_atoms(A1[i]->g().pos, '\0', d1, d2)) {
        if (m->image_idx != 0) continue;             // exclude symmetry mates
        Atom *b = mark_to_atom(mw, m);
        if (A1[i] == b) continue;
        auto it = a2idx.find(b);
        if (it == a2idx.end()) continue;
        if (seqNeglect(A1[i], b, seqDist)) continue;
        found.push_back({i, it->second, group, A1[i]->g().pos.dist(b->g().pos)});
      }
    }
  }
  alloc_contacts(found, contact, ncontacts);
}

void Manager::SeekContacts(PPAtom A, int n, realtype d1, realtype d2,
    int seqDist, PContact &contact, int &ncontacts, int /*maxlen*/,
    pmat44 /*TMatrix*/, long group) {
  std::vector<Contact> found;
  if (n > 0) {
    Model *mw = A[0]->res->chain->model;
    gemmi::NeighborSearch ns(mw->g(), st.cell, d2);
    ns.populate(true);
    std::unordered_map<Atom *, int> idx;
    idx.reserve(n * 2);
    for (int i = 0; i < n; ++i) idx.emplace(A[i], i);

    for (int i = 0; i < n; ++i) {
      if (A[i]->res->chain->model != mw) continue;
      for (auto *m : ns.find_atoms(A[i]->g().pos, '\0', d1, d2)) {
        if (m->image_idx != 0) continue;
        Atom *b = mark_to_atom(mw, m);
        auto it = idx.find(b);
        if (it == idx.end() || it->second <= i) continue;  // unordered pairs, once
        if (seqNeglect(A[i], b, seqDist)) continue;
        found.push_back({i, it->second, group, A[i]->g().pos.dist(b->g().pos)});
      }
    }
  }
  alloc_contacts(found, contact, ncontacts);
}

} // namespace mmdb
