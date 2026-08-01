/*
 * lidia-core/pdbqt-ligand.cc
 *
 * Flexible-ligand PDBQT writer (RDKit): Gasteiger charges + torsion tree.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <fstream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "utils/coot-utils.hh"
#include "lidia-core/pdbqt-common.hh"
#include "lidia-core/pdbqt-ligand.hh"
#include "lidia-core/rdkit-interface.hh"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>

namespace {

   typedef coot::pdbqt::writable_atom_t wr_atom_t;

   // Is this hydrogen a polar hydrogen (bonded to N, O or S)? Polar Hs are the
   // AutoDock "HD" donor hydrogens and are kept; non-polar Hs are merged away.
   bool is_polar_hydrogen(const RDKit::Atom *at, const RDKit::ROMol &mol) {
      if (at->getAtomicNum() != 1) return false;
      RDKit::ROMol::ADJ_ITER nbr, end;
      boost::tie(nbr, end) = mol.getAtomNeighbors(at);
      for (; nbr != end; ++nbr) {
         int z = mol[*nbr]->getAtomicNum();
         if (z == 7 || z == 8 || z == 16) return true;
      }
      return false;
   }

   // A nitrogen is an AutoDock acceptor (NA) if it carries a lone pair: no attached
   // hydrogens, not positively charged, and not fully substituted. getTotalNumHs(true)
   // counts explicit hydrogen atom neighbours as well as implicit hydrogens.
   bool nitrogen_is_acceptor(const RDKit::Atom *at) {
      if (at->getTotalNumHs(true) > 0) return false;
      if (at->getFormalCharge() > 0) return false;
      if (at->getTotalDegree() >= 4) return false;
      return true;
   }

   // A sulfur is treated as an AutoDock acceptor (SA) when it has a lone pair, i.e.
   // it is a thioether/thiol (degree < 3) rather than a sulfonyl/sulfate sulfur.
   bool sulfur_is_acceptor(const RDKit::Atom *at) {
      return at->getTotalDegree() < 3;
   }

   // Map an RDKit atom to its AutoDock (AD4) atom type.
   std::string autodock_atom_type(const RDKit::Atom *at, const RDKit::ROMol &mol) {
      int z = at->getAtomicNum();
      switch (z) {
      case 1:  return is_polar_hydrogen(at, mol) ? "HD" : "H";
      case 6:  return at->getIsAromatic() ? "A" : "C";
      case 7:  return nitrogen_is_acceptor(at) ? "NA" : "N";
      case 8:  return "OA";
      case 16: return sulfur_is_acceptor(at) ? "SA" : "S";
      case 9:  return "F";
      case 15: return "P";
      case 17: return "Cl";
      case 35: return "Br";
      case 53: return "I";
      case 5:  return "B";
      case 12: return "Mg";
      case 20: return "Ca";
      case 25: return "Mn";
      case 26: return "Fe";
      case 29: return "Cu";
      case 30: return "Zn";
      case 34: return "Se";
      default: return at->getSymbol();
      }
   }

   // Build the list of atoms to write for a single residue: compute Gasteiger
   // charges and AutoDock types via RDKit, keep heavy atoms and polar hydrogens,
   // and merge each non-polar hydrogen's charge into its parent heavy atom. The
   // writable_atom_t::rdkit_idx lets the caller build a torsion tree over the
   // RDKit topology. Returns an empty vector on failure.
   std::vector<wr_atom_t>
   written_atoms_for_residue(mmdb::Residue *residue_p,
                             const coot::dictionary_residue_restraints_t &restraints) {

      std::vector<wr_atom_t> result;
      try {
         RDKit::RWMol mol = coot::rdkit_mol(residue_p, restraints, "", true);
         if (mol.getNumAtoms() == 0)
            return result;
         RDKit::MolOps::sanitizeMol(mol);
         RDKit::computeGasteigerCharges(mol);

         // name -> mmdb atom, so RDKit atoms map back to authoritative coordinates
         std::map<std::string, mmdb::Atom *> name_to_mmdb;
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer())
               name_to_mmdb[std::string(at->name)] = at;
         }

         unsigned int n = mol.getNumAtoms();
         std::vector<double> charge(n, 0.0);
         for (unsigned int i=0; i<n; i++) {
            const RDKit::Atom *a = mol.getAtomWithIdx(i);
            double q = 0.0;
            if (a->hasProp("_GasteigerCharge")) {
               a->getProp("_GasteigerCharge", q);
               if (std::isnan(q) || std::isinf(q)) q = 0.0;
            }
            charge[i] = q;
         }

         // Merge non-polar hydrogens: add their charge to the bonded heavy atom.
         std::vector<bool> keep(n, true);
         for (unsigned int i=0; i<n; i++) {
            const RDKit::Atom *a = mol.getAtomWithIdx(i);
            if (a->getAtomicNum() == 1 && ! is_polar_hydrogen(a, mol)) {
               keep[i] = false;
               RDKit::ROMol::ADJ_ITER nbr, end;
               boost::tie(nbr, end) = mol.getAtomNeighbors(a);
               for (; nbr != end; ++nbr)
                  charge[*nbr] += charge[i];
            }
         }

         for (unsigned int i=0; i<n; i++) {
            if (! keep[i]) continue;
            const RDKit::Atom *a = mol.getAtomWithIdx(i);
            std::string name;
            if (! a->hasProp("name")) continue;
            a->getProp("name", name);
            std::map<std::string, mmdb::Atom *>::const_iterator it = name_to_mmdb.find(name);
            if (it == name_to_mmdb.end()) continue;
            std::string t = autodock_atom_type(a, mol);
            result.push_back(wr_atom_t(i, it->second, charge[i], t));
         }
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: written_atoms_for_residue() " << residue_p->GetResName()
                   << " " << e.what() << std::endl;
         result.clear();
      }
      return result;
   }

   // mmdb-only fallback: no dictionary, so write every atom with element-based
   // typing and zero charge.
   std::vector<wr_atom_t> written_atoms_fallback(mmdb::Residue *residue_p) {
      std::vector<wr_atom_t> result;
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (at->isTer()) continue;
         std::string t = coot::pdbqt::ad_type_from_element(at->element);
         result.push_back(wr_atom_t(-1, at, 0.0, t));
      }
      return result;
   }

   // A rotatable bond of the ligand's written-atom graph: an acyclic single bond
   // between two non-terminal heavy atoms, excluding amide C-N bonds.
   bool bond_is_rotatable(const RDKit::Bond *bond, const RDKit::ROMol &mol,
                          const std::map<int, int> &written_degree) {
      if (bond->getBondType() != RDKit::Bond::SINGLE) return false;
      if (bond->getIsAromatic()) return false;
      if (mol.getRingInfo()->numBondRings(bond->getIdx()) > 0) return false;
      const RDKit::Atom *a = bond->getBeginAtom();
      const RDKit::Atom *b = bond->getEndAtom();
      if (a->getAtomicNum() == 1 || b->getAtomicNum() == 1) return false;
      // both endpoints must be non-terminal once non-polar Hs are merged away
      std::map<int,int>::const_iterator ia = written_degree.find(a->getIdx());
      std::map<int,int>::const_iterator ib = written_degree.find(b->getIdx());
      if (ia == written_degree.end() || ib == written_degree.end()) return false;
      if (ia->second < 2 || ib->second < 2) return false;
      // exclude amide C-N bonds (kept planar by AutoDock convention)
      const RDKit::Atom *c_at = 0; const RDKit::Atom *n_at = 0;
      if (a->getAtomicNum() == 6 && b->getAtomicNum() == 7) { c_at = a; n_at = b; }
      if (b->getAtomicNum() == 6 && a->getAtomicNum() == 7) { c_at = b; n_at = a; }
      if (c_at && n_at) {
         RDKit::ROMol::ADJ_ITER nbr, end;
         boost::tie(nbr, end) = mol.getAtomNeighbors(c_at);
         for (; nbr != end; ++nbr) {
            const RDKit::Bond *cb = mol.getBondBetweenAtoms(c_at->getIdx(), *nbr);
            if (cb && cb->getBondType() == RDKit::Bond::DOUBLE &&
                mol[*nbr]->getAtomicNum() == 8)
               return false; // carbonyl carbon -> amide bond
         }
      }
      return true;
   }

   // union-find over written atoms, joined by every non-rotatable bond
   class union_find {
   public:
      std::map<int,int> parent;
      int find(int x) {
         if (parent.find(x) == parent.end()) parent[x] = x;
         while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
         return x;
      }
      void join(int a, int b) { parent[find(a)] = find(b); }
   };
}

int
coot::pdbqt::write_ligand(mmdb::Residue *residue_p,
                          const coot::dictionary_residue_restraints_t *restraints_p,
                          const std::string &comp_id,
                          const std::string &cid,
                          const std::string &file_name) {

   if (! residue_p) return 0;

   std::ofstream f(file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: write_ligand(): cannot open " << file_name << std::endl;
      return 0;
   }
   f << "REMARK  PDBQT file written by Coot (flexible ligand)\n";
   f << "REMARK  ligand " << comp_id << " " << cid << "\n";

   // Without a dictionary we cannot perceive rotatable bonds: write rigidly.
   std::vector<wr_atom_t> atoms;
   if (restraints_p)
      atoms = written_atoms_for_residue(residue_p, *restraints_p);
   if (atoms.empty()) {
      atoms = written_atoms_fallback(residue_p);
      f << "ROOT\n";
      int serial = 0;
      int n_written = 0;
      for (unsigned int i=0; i<atoms.size(); i++) {
         atoms[i].serial = ++serial;
         f << coot::pdbqt::atom_line(atoms[i], true) << "\n";
         n_written++;
      }
      f << "ENDROOT\nTORSDOF 0\n";
      f.close();
      std::cout << "INFO:: write_ligand(): wrote " << n_written
                << " atoms (rigid, no dictionary) to " << file_name << std::endl;
      return n_written;
   }

   // Rebuild the RDKit mol once so we can reason about its topology.
   RDKit::RWMol mol = coot::rdkit_mol(residue_p, *restraints_p, "", true);
   try { RDKit::MolOps::sanitizeMol(mol); } catch (const std::exception &e) {
      std::cout << "WARNING:: write_ligand(): sanitize " << e.what() << std::endl;
   }

   // written_degree: for each written atom, count of neighbours that are also written
   std::set<int> written_idx;
   std::map<int, wr_atom_t*> idx_to_wr;
   for (unsigned int i=0; i<atoms.size(); i++) {
      written_idx.insert(atoms[i].rdkit_idx);
      idx_to_wr[atoms[i].rdkit_idx] = &atoms[i];
   }
   std::map<int,int> written_degree;
   for (std::set<int>::const_iterator it=written_idx.begin(); it!=written_idx.end(); ++it) {
      int deg = 0;
      const RDKit::Atom *a = mol.getAtomWithIdx(*it);
      RDKit::ROMol::ADJ_ITER nbr, end;
      boost::tie(nbr, end) = mol.getAtomNeighbors(a);
      for (; nbr != end; ++nbr)
         if (written_idx.find(*nbr) != written_idx.end()) deg++;
      written_degree[*it] = deg;
   }

   // Classify each written-written bond as rotatable or rigid; fuse rigid bonds
   // into fragments with union-find.
   union_find uf;
   for (std::set<int>::const_iterator it=written_idx.begin(); it!=written_idx.end(); ++it)
      uf.find(*it); // ensure singletons exist
   std::vector<std::pair<int,int> > rotatable_bonds; // (begin_idx, end_idx), both written
   for (unsigned int ib=0; ib<mol.getNumBonds(); ib++) {
      const RDKit::Bond *bond = mol.getBondWithIdx(ib);
      int ai = bond->getBeginAtomIdx();
      int bi = bond->getEndAtomIdx();
      if (written_idx.find(ai) == written_idx.end()) continue;
      if (written_idx.find(bi) == written_idx.end()) continue;
      if (bond_is_rotatable(bond, mol, written_degree))
         rotatable_bonds.push_back(std::make_pair(ai, bi));
      else
         uf.join(ai, bi);
   }

   // Fragment id per written atom, and atom count per fragment.
   std::map<int,int> frag_of;   // atom idx -> fragment root
   std::map<int,int> frag_size;
   for (std::set<int>::const_iterator it=written_idx.begin(); it!=written_idx.end(); ++it) {
      int root = uf.find(*it);
      frag_of[*it] = root;
      frag_size[root]++;
   }

   // Fragment adjacency across rotatable bonds. Each entry: neighbour fragment and
   // the specific (this-side atom, other-side atom) that define the rotatable bond.
   struct branch_t { int other_frag; int this_atom; int other_atom; };
   std::map<int, std::vector<branch_t> > frag_adj;
   for (unsigned int i=0; i<rotatable_bonds.size(); i++) {
      int a = rotatable_bonds[i].first;
      int b = rotatable_bonds[i].second;
      int fa = frag_of[a];
      int fb = frag_of[b];
      if (fa == fb) continue; // safety
      branch_t ba; ba.other_frag = fb; ba.this_atom = a; ba.other_atom = b;
      branch_t bb; bb.other_frag = fa; bb.this_atom = b; bb.other_atom = a;
      frag_adj[fa].push_back(ba);
      frag_adj[fb].push_back(bb);
   }

   // Root fragment = the largest fragment.
   int root_frag = -1;
   int best = -1;
   for (std::map<int,int>::const_iterator it=frag_size.begin(); it!=frag_size.end(); ++it) {
      if (it->second > best) { best = it->second; root_frag = it->first; }
   }

   // Atoms in a given fragment, with an optional entry atom written first (the
   // atom on this side of the incoming rotatable bond, as PDBQT requires).
   std::map<int, std::vector<int> > frag_atoms;
   for (std::set<int>::const_iterator it=written_idx.begin(); it!=written_idx.end(); ++it)
      frag_atoms[frag_of[*it]].push_back(*it);

   int serial = 0;
   int n_written = 0;
   int n_torsions = 0;
   std::set<int> visited_frag;

   // Recursive tree writer. entry_atom is written first in the fragment (-1 for root).
   std::function<void(int,int)> write_fragment = [&](int frag, int entry_atom) {
      visited_frag.insert(frag);
      std::vector<int> order;
      if (entry_atom >= 0) order.push_back(entry_atom);
      const std::vector<int> &fa = frag_atoms[frag];
      for (unsigned int i=0; i<fa.size(); i++)
         if (fa[i] != entry_atom) order.push_back(fa[i]);
      for (unsigned int i=0; i<order.size(); i++) {
         wr_atom_t *wa = idx_to_wr[order[i]];
         wa->serial = ++serial;
         f << coot::pdbqt::atom_line(*wa, true) << "\n";
         n_written++;
      }
      // Branch out to unvisited neighbouring fragments.
      const std::vector<branch_t> &adj = frag_adj[frag];
      for (unsigned int i=0; i<adj.size(); i++) {
         int nfrag = adj[i].other_frag;
         if (visited_frag.find(nfrag) != visited_frag.end()) continue;
         int this_serial = idx_to_wr[adj[i].this_atom]->serial;
         int child_serial = serial + 1; // the entry atom is written first in the child
         f << "BRANCH " << std::setw(3) << this_serial << " "
           << std::setw(3) << child_serial << "\n";
         n_torsions++;
         write_fragment(nfrag, adj[i].other_atom);
         f << "ENDBRANCH " << std::setw(3) << this_serial << " "
           << std::setw(3) << idx_to_wr[adj[i].other_atom]->serial << "\n";
      }
   };

   f << "ROOT\n";
   // Write only the root fragment's own atoms between ROOT/ENDROOT, then branches.
   {
      const std::vector<int> &fa = frag_atoms[root_frag];
      visited_frag.insert(root_frag);
      for (unsigned int i=0; i<fa.size(); i++) {
         wr_atom_t *wa = idx_to_wr[fa[i]];
         wa->serial = ++serial;
         f << coot::pdbqt::atom_line(*wa, true) << "\n";
         n_written++;
      }
      f << "ENDROOT\n";
      const std::vector<branch_t> &adj = frag_adj[root_frag];
      for (unsigned int i=0; i<adj.size(); i++) {
         int nfrag = adj[i].other_frag;
         if (visited_frag.find(nfrag) != visited_frag.end()) continue;
         int this_serial = idx_to_wr[adj[i].this_atom]->serial;
         int child_serial = serial + 1;
         f << "BRANCH " << std::setw(3) << this_serial << " "
           << std::setw(3) << child_serial << "\n";
         n_torsions++;
         write_fragment(nfrag, adj[i].other_atom);
         f << "ENDBRANCH " << std::setw(3) << this_serial << " "
           << std::setw(3) << idx_to_wr[adj[i].other_atom]->serial << "\n";
      }
   }
   f << "TORSDOF " << n_torsions << "\n";
   f.close();

   std::cout << "INFO:: write_ligand(): wrote " << n_written << " atoms, "
             << n_torsions << " torsions to " << file_name << std::endl;
   return n_written;
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS
