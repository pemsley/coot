// Test the hardened core: identity cache + localized-patch edits + stability.
#include "backing.hh"
#include <cstdio>

using namespace shim;

static int failures = 0;
#define CHECK(cond, msg)                                                        \
  do {                                                                          \
    if (!(cond)) { std::printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { std::printf("  ok:   %s\n", msg); }                          \
  } while (0)

// identity-encoded coords: x = res*10 + atom
static gemmi::Atom mk(int r, int a) {
  gemmi::Atom at;
  at.name = "A" + std::to_string(r) + "_" + std::to_string(a);
  at.pos = gemmi::Position(r * 10.0 + a, (double)r, (double)a);
  return at;
}

int main() {
  Backing B;
  B.st.models.emplace_back();
  auto &gm = B.st.models.back();
  gm.chains.emplace_back();
  gm.chains.back().name = "A";
  for (int r = 0; r < 3; ++r) {
    gemmi::Residue res; res.name = "GLY"; res.seqid = gemmi::SeqId(r + 1, ' ');
    for (int a = 0; a < 3; ++a) res.atoms.push_back(mk(r, a));
    gm.chains.back().residues.push_back(res);
  }
  B.build_from_gemmi();

  ChainW *chain = B.GetModel(0)->GetChain(0);
  ResidueW *res1 = chain->GetResidue(1);
  AtomW *held = res1->GetAtom(1);            // (res1, atom1) -> x==11

  std::printf("held=%p x=%.1f name=%s\n", (void *)held, held->x(), held->name());

  // --- identity cache: repeated Get* return the SAME pointer ---
  CHECK(chain->GetResidue(1) == res1, "GetResidue(1) twice -> same pointer (identity)");
  CHECK(res1->GetAtom(1) == held, "GetAtom(1) twice -> same pointer (identity)");
  CHECK(held->x() == 11.0, "resolves to correct atom");

  // --- edit 1: AddAtom x5000 (gemmi atoms vector reallocates) ---
  const void *d0 = (void *)res1->g().atoms.data();
  for (int i = 0; i < 5000; ++i) res1->AddAtom(B, mk(1, 100 + i));
  CHECK(d0 != (void *)res1->g().atoms.data(), "gemmi atoms vector reallocated");
  CHECK(held->x() == 11.0, "held correct after AddAtom (append, no shift)");
  CHECK(res1->GetAtom(1) == held, "identity preserved after AddAtom");

  // --- edit 2: InsResidue at front: only THIS chain's residue ri shifts ---
  //     Atom wrappers of other residues must NOT be touched.
  gemmi::Residue newr; newr.name = "ACE"; newr.seqid = gemmi::SeqId(0, ' ');
  newr.atoms.push_back(mk(7, 7));
  chain->InsResidue(B, 0, newr);
  CHECK(res1->ri == 2, "res1 index patched 1 -> 2 (localized to chain)");
  CHECK(held->ai == 1, "held atom index UNCHANGED by residue insert (localized)");
  CHECK(held->x() == 11.0, "held still resolves to same logical atom after mid-insert");
  CHECK(chain->GetResidue(2) == res1, "chain now finds res1 at index 2 (same pointer)");
  CHECK(chain->GetResidue(0)->GetAtom(0)->x() == 77.0, "inserted residue resolves correctly");

  // --- edit 3: DeleteAtom at index 0 of res1: atom indices shift within residue only ---
  res1->DeleteAtom(0);
  CHECK(held->ai == 0, "held atom index patched 1 -> 0 after delete");
  CHECK(held->x() == 11.0, "held still correct after DeleteAtom");
  CHECK(res1->GetAtom(0) == held, "identity preserved after DeleteAtom");

  // --- edit 4: write through reference accessor reaches live gemmi ---
  held->x() = 999.0;
  CHECK(res1->g().atoms[0].pos.x == 999.0, "ref-accessor write reached live gemmi storage");

  std::printf("\n=== %s (%d failures) ===\n",
              failures == 0 ? "CORE PASSED" : "CORE FAILED", failures);
  return failures ? 1 : 0;
}
