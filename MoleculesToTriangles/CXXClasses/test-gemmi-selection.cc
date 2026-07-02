/*
 * MoleculesToTriangles/CXXClasses/test-gemmi-selection.cc
 *
 * Unit test for the gemmi-native selection module. Checks that
 * compound_selection_t reproduces expected atom counts (including the
 * boolean & | ! combinations and {} nesting) against a known structure.
 *
 * This file is part of Coot
 */

#include "gemmi-selection.hh"
#include <gemmi/pdb.hpp>
#include <iostream>
#include <string>
#include <set>
#include <functional>

static int g_fail = 0;

static void check(const std::string &name, int got, int expected) {
   bool ok = (got == expected);
   std::cout << (ok ? "  PASS  " : "  FAIL  ");
   std::cout.width(34); std::cout << std::left << name;
   std::cout << " got=" << got << " expected=" << expected << "\n";
   if (!ok) g_fail++;
}

// independent count over the model by a per-atom predicate (chain id, res name, atom name)
static int count_if(gemmi::Model &model,
                    const std::function<bool(const std::string&, const std::string&, const std::string&)> &pred) {
   int n = 0;
   for (gemmi::Chain &c : model.chains)
      for (gemmi::Residue &r : c.residues)
         for (gemmi::Atom &a : r.atoms)
            if (pred(c.name, r.name, a.name)) n++;
   return n;
}

int main(int argc, char **argv) {
   std::string path = argc > 1 ? argv[1] : "reference-structures/1e9g.pdb";
   gemmi::Structure st = gemmi::read_pdb_file(path);
   if (st.models.empty()) { std::cout << "no models in " << path << "\n"; return 2; }
   gemmi::Model &model = st.models[0];

   int total = count_if(model, [](const std::string&, const std::string&, const std::string&){ return true; });
   std::cout << "Loaded " << path << " : " << total << " atoms in model 0\n\n";

   const std::set<std::string> waters = {"WAT","HOH","OH2","H2O"};
   const std::set<std::string> mainch = {"N","CA","C","O","H"};

   int n_water  = count_if(model, [&](const std::string&, const std::string &rn, const std::string&){
      return waters.count(rn) > 0; });
   int n_chainA = count_if(model, [](const std::string &cn, const std::string&, const std::string&){ return cn == "A"; });
   int n_A_water = count_if(model, [&](const std::string &cn, const std::string &rn, const std::string&){
      return cn == "A" && waters.count(rn) > 0; });
   int n_A_or_B  = count_if(model, [](const std::string &cn, const std::string&, const std::string&){
      return cn == "A" || cn == "B"; });
   int n_A_main  = count_if(model, [&](const std::string &cn, const std::string&, const std::string &an){
      return cn == "A" && mainch.count(an) > 0; });

   auto run = [&](const std::string &sel) {
      coot::m2t::compound_selection_t cs(sel);
      return cs.count_matching_atoms(model);
   };

   // basic CID leaves
   check("ALL /*/*/*/*",            run("/*/*/*/*"),            total);
   check("chain A (mmdb CID)",      run("//A/*.*/*:*"),         n_chainA);
   check("WATER subset name",       run("WATER"),               n_water);
   check("ALL subset name",         run("ALL"),                 total);

   // invert
   check("! WATER",                 run("!WATER"),              total - n_water);

   // boolean combinations
   check("//A & WATER (disjoint=0)", run("//A/*.*/*:* & WATER"), n_A_water);
   check("//A & MAIN (nonzero)",    run("//A/*.*/*:* & MAIN"),  n_A_main);
   check("//A | //B",               run("//A/*.*/*:* | //B/*.*/*:*"), n_A_or_B);

   // brace nesting
   check("{ //A } & WATER",         run("{ //A/*.*/*:* } & WATER"), n_A_water);

   std::cout << "\n" << (g_fail == 0 ? "ALL TESTS PASSED" : "FAILURES: " + std::to_string(g_fail))
             << "\n";
   return g_fail == 0 ? 0 : 1;
}
