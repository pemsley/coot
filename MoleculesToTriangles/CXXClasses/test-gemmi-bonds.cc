/*
 * MoleculesToTriangles/CXXClasses/test-gemmi-bonds.cc
 *
 * Unit test for gemmi-native bond determination.
 *
 * This file is part of Coot
 */

#include "gemmi-bonds.hh"
#include <gemmi/pdb.hpp>
#include <iostream>
#include <map>
#include <set>
#include <string>

static int g_fail = 0;
static void check(const std::string &name, bool ok, const std::string &detail = "") {
   std::cout << (ok ? "  PASS  " : "  FAIL  ");
   std::cout.width(40); std::cout << std::left << name << " " << detail << "\n";
   if (!ok) g_fail++;
}

int main(int argc, char **argv) {
   std::string path = argc > 1 ? argv[1] : "reference-structures/1e9g.pdb";
   gemmi::Structure st = gemmi::read_pdb_file(path);
   if (st.models.empty()) { std::cout << "no models\n"; return 2; }

   std::vector<coot::m2t::bond_t> bonds = coot::m2t::make_bonds(st);

   // count atoms (and assign-serial side effect already done by make_bonds)
   int n_atoms = 0;
   gemmi::Model &model = st.models[0];
   for (gemmi::Chain &c : model.chains)
      for (gemmi::Residue &r : c.residues)
         n_atoms += (int)r.atoms.size();

   std::cout << "Loaded " << path << " : " << n_atoms << " atoms, "
             << bonds.size() << " bonds\n\n";

   auto adj = coot::m2t::bonds_adjacency(n_atoms, bonds);

   // sanity: protein bond/atom ratio is ~1.0-1.1 (each bond counted once)
   double ratio = n_atoms ? double(bonds.size()) / n_atoms : 0.0;
   check("bond/atom ratio in [0.8, 1.2]", ratio > 0.8 && ratio < 1.2,
         "ratio=" + std::to_string(ratio));

   // find the first standard amino-acid residue and verify backbone connectivity
   auto serial_of = [&](gemmi::Residue &res, const std::string &nm) -> int {
      for (gemmi::Atom &a : res.atoms) if (a.name == nm) return a.serial;
      return -1;
   };
   auto bonded = [&](int s1, int s2) {
      if (s1 < 0 || s2 < 0) return false;
      for (int n : adj[s1]) if (n == s2) return true;
      return false;
   };

   const std::set<std::string> aa = {
      "ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
      "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};

   gemmi::Residue *res_p = nullptr;
   for (gemmi::Chain &c : model.chains) {
      for (gemmi::Residue &r : c.residues)
         if (aa.count(r.name) && serial_of(r, "CB") >= 0) { res_p = &r; break; }
      if (res_p) break;
   }
   check("found a standard residue with CB", res_p != nullptr);
   if (res_p) {
      int N = serial_of(*res_p, "N"), CA = serial_of(*res_p, "CA");
      int C = serial_of(*res_p, "C"), O = serial_of(*res_p, "O"), CB = serial_of(*res_p, "CB");
      std::string where = res_p->name + " " + res_p->seqid.str();
      check("N-CA bonded",  bonded(N, CA),  where);
      check("CA-C bonded",  bonded(CA, C),  where);
      check("C-O bonded",   bonded(C, O),   where);
      check("CA-CB bonded", bonded(CA, CB), where);
      check("N-O NOT bonded (sanity)", !bonded(N, O), where);
   }

   std::cout << "\n" << (g_fail == 0 ? "ALL TESTS PASSED" : "FAILURES: " + std::to_string(g_fail)) << "\n";
   return g_fail == 0 ? 0 : 1;
}
