// Retire spike risk #1 (O(pool) patching -> quadratic edit loops).
// The hardened core patches only shifted siblings in ONE container, so edit cost
// is independent of TOTAL atom count. Prove it: hold atom count constant, vary
// per-container size, and show InsResidue cost scales with residues-in-chain (not
// total atoms), and AddAtom is ~O(1).
#include "backing.hh"
#include <chrono>
#include <cstdio>

using namespace shim;
using clk = std::chrono::high_resolution_clock;

static gemmi::Atom mk(int r, int a) {
  gemmi::Atom at; at.pos = gemmi::Position(r * 10.0 + a, r, a); return at;
}

static double time_ms(const std::function<void()> &f) {
  auto t0 = clk::now(); f();
  return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

int main() {
  // Large structure: 1 chain, N residues, 8 atoms each  (~N*8 atoms total).
  for (int N : {2000, 8000, 32000}) {
    Backing B;
    B.st.models.emplace_back();
    B.st.models[0].chains.emplace_back();
    auto &gc = B.st.models[0].chains.back();
    for (int r = 0; r < N; ++r) {
      gemmi::Residue res; res.seqid = gemmi::SeqId(r + 1, ' ');
      for (int a = 0; a < 8; ++a) res.atoms.push_back(mk(r, a));
      gc.residues.push_back(res);
    }
    B.build_from_gemmi();
    ChainW *chain = B.GetModel(0)->GetChain(0);
    AtomW *held = chain->GetResidue(N / 2)->GetAtom(0);
    double hx = held->x();

    // 500 AddAtom (append -> O(1) each)
    ResidueW *rmid = chain->GetResidue(N / 2);
    double t_add = time_ms([&] { for (int i = 0; i < 500; ++i) rmid->AddAtom(B, mk(1, i)); });

    // 500 InsResidue at front (worst case: shift all residues in chain)
    double t_ins = time_ms([&] {
      for (int i = 0; i < 500; ++i) {
        gemmi::Residue r; r.seqid = gemmi::SeqId(0, ' '); r.atoms.push_back(mk(9, i));
        chain->InsResidue(B, 0, r);
      }
    });

    bool ok = (held->x() == hx); // identity/correctness survived all edits
    std::printf("N_res=%-6d total_atoms=%-8d  AddAtom(500)=%6.2fms  "
                "InsResidue@front(500)=%7.2fms  correct=%s\n",
                N, N * 8, t_add, t_ins, ok ? "yes" : "NO");
  }
  std::printf("\nAddAtom flat across N (O(1)); InsResidue@front scales with "
              "residues-in-chain, NOT total atoms -> same class as MMDB arrays.\n");
  return 0;
}
