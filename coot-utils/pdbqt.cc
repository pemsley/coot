/*
 * coot-utils/pdbqt.cc
 *
 * Rigid receptor PDBQT writer and shared per-atom formatting/typing (no RDKit).
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>

#include "utils/coot-utils.hh"
#include "coot-utils/pdbqt.hh"
#include "coot-utils/pdbqt-charges.hh"

// ---------------------------------------------------------------------------
//   Native (no-RDKit) typing and charges for the receptor, driven by the
//   CCP4 monomer-library energy type (type_energy from the dictionary) mapped
//   through the hb_type/element information in ener_lib.cif.
// ---------------------------------------------------------------------------

namespace {

   typedef coot::pdbqt::writable_atom_t wr_atom_t;

   // Aromatic-ring carbon energy types (ener_lib.cif name CR..) -> AutoDock "A".
   bool energy_type_is_aromatic_carbon(const std::string &et) {
      return et.size() >= 2 && et[0] == 'C' && et[1] == 'R';
   }

   // AutoDock atom type for a heavy atom from its energy-lib type and element.
   // Falls back to element-only typing when the energy type is unknown/empty.
   // (Nitrogen and cysteine-SG are decided from the hydrogens present by the
   // caller, so this only needs the coarse element split.)
   std::string ad_type_native(const std::string &et, const std::string &element) {
      std::string e = coot::util::remove_leading_spaces(element);
      if (e == "C") return energy_type_is_aromatic_carbon(et) ? "A" : "C";
      if (e == "N") return "N";
      if (e == "O") return "OA";
      if (e == "S") return "SA";
      return coot::pdbqt::ad_type_from_element(element);
   }

   // Partial charge for a heavy atom keyed by CCP4 energy-lib type - a coarse
   // per-type fallback used for residues not covered by the reference table.
   double charge_native(const std::string &et, const std::string &element) {
      static const std::map<std::string, double> table = {
         // carbons
         {"C", 0.50}, {"C1", 0.10}, {"C2", 0.10}, {"CSP", 0.00}, {"CSP1", 0.00},
         {"CH1", 0.10}, {"CH2", 0.00}, {"CH3", -0.05}, {"CT", 0.00},
         {"CR1", 0.00}, {"CR2", 0.00}, {"CR1H", 0.00}, {"CR15", 0.00}, {"CR5", 0.00},
         {"CR56", 0.00}, {"CR55", 0.00}, {"CR16", 0.00}, {"CR6", 0.00}, {"CR66", 0.00},
         // nitrogens
         {"N", -0.30}, {"NH0", -0.30}, {"NH1", -0.35}, {"NH2", -0.40},
         {"NC1", -0.30}, {"NC2", -0.30}, {"NC3", -0.30},
         {"NT", -0.30}, {"NT1", -0.35}, {"NT2", -0.40}, {"NT3", -0.30}, {"NT4", -0.30},
         {"NS", -0.20}, {"NSP", -0.20}, {"NS1", -0.20}, {"NSP1", -0.20},
         {"N20", -0.25}, {"N21", -0.20}, {"N30", -0.25}, {"N31", -0.20}, {"N32", -0.20},
         {"N33", -0.20}, {"NR5", -0.20}, {"NR15", -0.25}, {"NR16", -0.25}, {"NR55", -0.20},
         {"NR56", -0.20}, {"NR6", -0.20}, {"NR66", -0.20}, {"NRD5", -0.25}, {"NRD6", -0.25},
         {"NPA", -0.20}, {"NPB", -0.20},
         // oxygens
         {"O", -0.45}, {"OS", -0.40}, {"OB", -0.45}, {"O2", -0.45},
         {"OC", -0.60}, {"OC2", -0.60}, {"OP", -0.60},
         {"OH1", -0.55}, {"OH2", -0.55}, {"OHA", -0.55}, {"OHB", -0.55}, {"OHC", -0.55},
         // sulfurs
         {"S", -0.15}, {"S1", -0.15}, {"S2", -0.15}, {"S3", -0.10}, {"ST", -0.10}, {"SH1", -0.20},
         // phosphorus
         {"P", 0.50}, {"P1", 0.50}, {"PS", 0.30},
      };
      std::map<std::string, double>::const_iterator it = table.find(et);
      if (it != table.end()) return it->second;
      std::string e = coot::util::remove_leading_spaces(element);
      if (e == "O") return -0.40;
      if (e == "N") return -0.30;
      if (e == "S") return -0.15;
      if (e == "P") return  0.40;
      return 0.00;
   }

   // Polar-hydrogen charge, by the element of the heavy atom it is bonded to.
   double polar_h_charge(const std::string &heavy_element) {
      std::string e = coot::util::remove_leading_spaces(heavy_element);
      if (e == "O") return 0.40;
      if (e == "N") return 0.25;
      if (e == "S") return 0.15;
      return 0.20;
   }
   const double NONPOLAR_H_CHARGE = 0.06; // merged into the parent carbon

   bool element_is_polar_acceptor(const std::string &element) {
      std::string e = coot::util::remove_leading_spaces(element);
      return e == "N" || e == "O" || e == "S";
   }

   bool residue_is_water(const std::string &resname) {
      return resname == "HOH" || resname == "WAT" || resname == "DOD" || resname == "H2O";
   }

   bool residue_is_polymer(mmdb::Residue *res) {
      if (!res) return false;
      return res->isAminoacid() || res->isNucleotide();
   }

   // A cysteine bearing a thiol hydrogen (HG/HG1) is a protonated free thiol.
   bool cysteine_is_thiol(const std::set<std::string> &names) {
      return names.count("HG") || names.count("HG1");
   }

   // Look up a reference (Vina/prepare_receptor) charge for (residue, atom).
   // ok_out is set false on a miss. The backbone amide hydrogen is named "H" in
   // the table, so "HN" is tried as an alias.
   double reference_charge(const std::string &resname, const std::string &atom_name, bool *ok_out) {
      *ok_out = false;
      const std::map<std::string, std::map<std::string, float> > &t = coot::pdbqt_reference_charges();
      std::map<std::string, std::map<std::string, float> >::const_iterator it = t.find(resname);
      if (it == t.end()) return 0.0;
      std::map<std::string, float>::const_iterator jt = it->second.find(atom_name);
      if (jt == it->second.end() && atom_name == "HN") jt = it->second.find("H");
      if (jt == it->second.end()) return 0.0;
      *ok_out = true;
      return jt->second;
   }

   // Type and charge one residue natively (no RDKit): keep heavy atoms and polar
   // hydrogens (as HD); merge each non-polar hydrogen's charge into its parent.
   // Output preserves the residue's atom order. restraints_p may be null.
   std::vector<wr_atom_t>
   native_atoms_for_residue(mmdb::Residue *residue_p,
                            const coot::dictionary_residue_restraints_t *restraints_p) {
      std::vector<wr_atom_t> out;
      mmdb::PPAtom atoms = 0;
      int n_atoms = 0;
      residue_p->GetAtomTable(atoms, n_atoms);

      std::vector<bool> is_h(n_atoms, false), keep(n_atoms, true);
      std::vector<double> charge(n_atoms, 0.0);
      std::vector<std::string> adtype(n_atoms);

      std::set<std::string> names_present;
      for (int i=0; i<n_atoms; i++) {
         if (atoms[i]->isTer()) { keep[i] = false; continue; }
         std::string e = coot::util::remove_leading_spaces(atoms[i]->element);
         if (e == "H" || e == "D") is_h[i] = true;
         names_present.insert(coot::util::remove_whitespace(atoms[i]->name));
      }

      std::string resname(residue_p->GetResName());
      // Standard residues take their charges from the reference (Vina) table, whose
      // heavy-atom charges already include the merged non-polar hydrogens; other
      // residues use the per-energy-type charges and merge non-polar H explicitly.
      // A cysteine with no thiol hydrogen uses the disulfide/deprotonated set (CYX).
      std::string charge_res = resname;
      if (resname == "CYS" && ! cysteine_is_thiol(names_present)) charge_res = "CYX";
      bool use_ref = coot::pdbqt_reference_charges().count(charge_res) > 0;

      // Which heavy atoms bear a hydrogen in the model? Nitrogen typing is decided
      // from the hydrogens actually present, not from the dictionary's (fixed)
      // protonation/tautomer state - the CCP4 HIS monomer, for instance, is the
      // doubly-protonated form, so its dictionary energy types would wrongly make
      // both ring nitrogens donors.
      std::vector<int> parent(n_atoms, -1);
      std::vector<bool> heavy_has_H(n_atoms, false);
      for (int i=0; i<n_atoms; i++) {
         if (!keep[i] || !is_h[i]) continue;
         mmdb::Atom *h = atoms[i];
         int best = -1;
         double best_d2 = 1.0e9;
         for (int j=0; j<n_atoms; j++) {
            if (!keep[j] || is_h[j]) continue;
            double dx = h->x - atoms[j]->x, dy = h->y - atoms[j]->y, dz = h->z - atoms[j]->z;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < best_d2) { best_d2 = d2; best = j; }
         }
         if (best >= 0 && best_d2 <= 1.35*1.35) { parent[i] = best; heavy_has_H[best] = true; }
      }

      // heavy atoms: energy-lib type; reference charge if available, else energy-type
      for (int i=0; i<n_atoms; i++) {
         if (!keep[i] || is_h[i]) continue;
         std::string name(atoms[i]->name);
         std::string name_stripped = coot::util::remove_whitespace(name);
         std::string et = restraints_p ? restraints_p->type_energy(name) : std::string();
         adtype[i] = ad_type_native(et, atoms[i]->element);
         // Nitrogen acceptor/donor from the hydrogens present: the backbone amide N
         // is always a donor (N); a side-chain N carrying a hydrogen is a donor (N);
         // a side-chain N with no hydrogen and a free lone pair (<= 2 heavy
         // neighbours, e.g. a deprotonated histidine ring N) is an acceptor (NA),
         // while 3 heavy neighbours (proline, tertiary amide) stays N.
         std::string el = coot::util::remove_leading_spaces(atoms[i]->element);
         if (el == "N") {
            if (name_stripped == "N" || heavy_has_H[i]) {
               adtype[i] = "N";
            } else {
               int deg = 0;
               for (int j=0; j<n_atoms; j++) {
                  if (j == i || !keep[j] || is_h[j]) continue;
                  double dx = atoms[i]->x - atoms[j]->x, dy = atoms[i]->y - atoms[j]->y,
                         dz = atoms[i]->z - atoms[j]->z;
                  if (dx*dx + dy*dy + dz*dz < 1.9*1.9) deg++;
               }
               adtype[i] = (deg <= 2) ? "NA" : "N";
            }
         }
         // A cysteine SG with no thiol hydrogen - a disulfide, or a deprotonated
         // thiolate - is typed "S"; a protonated free thiol (HG present) is the
         // acceptor "SA". The decision follows the hydrogens present.
         if (resname == "CYS" && name_stripped == "SG" && ! cysteine_is_thiol(names_present))
            adtype[i] = "S";
         bool got_ref = false;
         double rq = reference_charge(charge_res, name_stripped, &got_ref);
         charge[i] = got_ref ? rq : charge_native(et, atoms[i]->element);
      }
      // hydrogens: the precomputed parent atom decides polar (HD) vs merge
      for (int i=0; i<n_atoms; i++) {
         if (!keep[i] || !is_h[i]) continue;
         int best = parent[i];
         bool polar = (best >= 0) ? element_is_polar_acceptor(atoms[best]->element)
                                  : true;   // orphan hydrogen: keep as HD
         if (polar) {
            adtype[i] = "HD";
            std::string hname = coot::util::remove_whitespace(atoms[i]->name);
            bool got_ref = false;
            double rq = reference_charge(charge_res, hname, &got_ref);
            if (got_ref)
               charge[i] = rq;
            else
               charge[i] = (best >= 0) ? polar_h_charge(atoms[best]->element) : 0.20;
         } else {
            keep[i] = false;                      // non-polar hydrogen
            if (! use_ref)                        // reference heavy charges already include it
               charge[best] += NONPOLAR_H_CHARGE;
         }
      }
      for (int i=0; i<n_atoms; i++)
         if (keep[i])
            out.push_back(wr_atom_t(-1, atoms[i], charge[i], adtype[i]));
      return out;
   }
}

// ---------------------------------------------------------------------------
//   Flexible-receptor helpers: which atoms of a flexible residue stay rigid,
//   simple distance bonding, and the dictionary side-chain rotatable bonds.
// ---------------------------------------------------------------------------

namespace {

   bool is_hydrogen_atom(mmdb::Atom *at) {
      std::string e = coot::util::remove_leading_spaces(at->element);
      return e == "H" || e == "D";
   }

   // The main-chain heavy-atom names (CA is deliberately excluded: it is the
   // root of the flexible side chain and so belongs to the flex file).
   bool name_is_rigid_backbone_heavy(const std::string &n) {
      return n == "N" || n == "C" || n == "O" || n == "OXT";
   }

   // For a flexible residue, which atoms stay in the RIGID file: the backbone
   // heavy atoms N/C/O/OXT and any hydrogen bonded to one of them (the amide H,
   // an OXT hydroxyl H). CA, HA and everything side-chain-ward go to the flex
   // file instead.
   bool rigid_keep_for_flex_residue(mmdb::Atom *at, mmdb::Residue *res) {
      std::string name = coot::util::remove_whitespace(at->name);
      if (! is_hydrogen_atom(at))
         return name_is_rigid_backbone_heavy(name);
      mmdb::PPAtom atoms = 0; int n_atoms = 0;
      res->GetAtomTable(atoms, n_atoms);
      int best = -1; double best_d2 = 1.35*1.35;
      for (int j=0; j<n_atoms; j++) {
         if (is_hydrogen_atom(atoms[j])) continue;
         double dx = at->x-atoms[j]->x, dy = at->y-atoms[j]->y, dz = at->z-atoms[j]->z;
         double d2 = dx*dx + dy*dy + dz*dz;
         if (d2 < best_d2) { best_d2 = d2; best = j; }
      }
      if (best < 0) return false;
      std::string hn = coot::util::remove_whitespace(atoms[best]->name);
      return hn == "N" || hn == "O" || hn == "OXT";
   }

   // Distance-based bond test (elements decide the cutoff). Robust to atom-name
   // conventions and matches the hydrogen-parent logic used elsewhere here.
   bool bonded_pdbqt(mmdb::Atom *a, mmdb::Atom *b) {
      double dx = a->x-b->x, dy = a->y-b->y, dz = a->z-b->z;
      double d2 = dx*dx + dy*dy + dz*dz;
      if (is_hydrogen_atom(a) || is_hydrogen_atom(b)) return d2 < 1.35*1.35;
      std::string ea = coot::util::remove_leading_spaces(a->element);
      std::string eb = coot::util::remove_leading_spaces(b->element);
      double lim = (ea == "S" || eb == "S") ? 2.05 : 1.85;
      return d2 < lim*lim;
   }

   typedef std::pair<std::string, std::string> name_pair_t;
   name_pair_t ordered_names(const std::string &a, const std::string &b) {
      return (a < b) ? name_pair_t(a, b) : name_pair_t(b, a);
   }

   // The set of central bonds of the dictionary's non-const, non-ring torsions,
   // keyed by (whitespace-stripped) atom-name pair - the candidate rotatable
   // side-chain bonds. Bonds that are entirely main chain (both central atoms
   // backbone) are dropped, so CA-CB survives but N-CA / CA-C do not. Hydrogen
   // torsions are included, so hydroxyl/thiol/amine rotations are candidates.
   std::set<name_pair_t>
   sidechain_rotatable_bonds(const coot::dictionary_residue_restraints_t &rest) {
      std::set<name_pair_t> bonds;
      std::vector<coot::dict_torsion_restraint_t> tors = rest.get_non_const_torsions(true);
      for (unsigned int i=0; i<tors.size(); i++) {
         std::string b2 = coot::util::remove_whitespace(tors[i].atom_id_2_4c());
         std::string b3 = coot::util::remove_whitespace(tors[i].atom_id_3_4c());
         if (name_is_rigid_backbone_heavy(b2) && name_is_rigid_backbone_heavy(b3)) continue;
         if (b2 == "N" && b3 == "CA") continue;   // phi side of the backbone
         if (b2 == "CA" && b3 == "C")  continue;   // psi side of the backbone
         coot::atom_name_quad quad(tors[i].atom_id_1_4c(), tors[i].atom_id_2_4c(),
                                   tors[i].atom_id_3_4c(), tors[i].atom_id_4_4c());
         if (rest.is_ring_torsion(quad)) continue;
         bonds.insert(ordered_names(b2, b3));
      }
      return bonds;
   }
}

// ---------------------------------------------------------------------------
//                     Rigid receptor / whole-molecule PDBQT
// ---------------------------------------------------------------------------

namespace {

   // Shared rigid-atom writer. Writes every (non-water) residue's native atoms;
   // for a residue whose spec is in flex_specs only the rigid backbone atoms are
   // written (the side chain goes to the flex file). Returns the atom count and
   // advances *serial_p.
   int write_rigid_atoms(std::ostream &f, mmdb::Manager *mol, coot::protein_geometry &geom,
                         int imol_enc, const std::set<coot::residue_spec_t> &flex_specs,
                         int *serial_p) {
      int n_written = 0;
      mmdb::Model *model_p = mol->GetModel(1);
      if (! model_p) return 0;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (! residue_p) continue;
            std::string res_name(residue_p->GetResName());
            if (residue_is_water(res_name)) continue; // receptor prep drops waters
            bool is_flex = flex_specs.count(coot::residue_spec_t(residue_p)) > 0;
            std::pair<bool, coot::dictionary_residue_restraints_t> rp =
               geom.get_monomer_restraints(res_name, imol_enc);
            const coot::dictionary_residue_restraints_t *rest_p = rp.first ? &rp.second : nullptr;
            std::vector<wr_atom_t> atoms = native_atoms_for_residue(residue_p, rest_p);
            bool is_het = ! residue_is_polymer(residue_p);
            for (unsigned int i=0; i<atoms.size(); i++) {
               if (is_flex && ! rigid_keep_for_flex_residue(atoms[i].at, residue_p))
                  continue;   // this atom belongs to the flexible side chain
               atoms[i].serial = ++(*serial_p);
               f << coot::pdbqt::atom_line(atoms[i], is_het) << "\n";
               n_written++;
            }
         }
      }
      return n_written;
   }
}

int
coot::pdbqt::write_receptor(mmdb::Manager *mol, coot::protein_geometry &geom, int imol_enc,
                            const std::string &file_name) {

   if (! mol) return 0;

   std::ofstream f(file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: write_receptor(): cannot open " << file_name << std::endl;
      return 0;
   }

   f << "REMARK  PDBQT file written by Coot (rigid receptor)\n";
   int serial = 0;
   std::set<coot::residue_spec_t> no_flex;
   int n_written = write_rigid_atoms(f, mol, geom, imol_enc, no_flex, &serial);
   f << "TER\n";
   f.close();
   std::cout << "INFO:: write_receptor(): wrote " << n_written
             << " atoms to " << file_name << std::endl;
   return n_written;
}

// ---------------------------------------------------------------------------
//                          Flexible receptor PDBQT
// ---------------------------------------------------------------------------

namespace {

   // Write one flexible side chain (BEGIN_RES ... END_RES) as a CA-rooted torsion
   // tree, using dictionary rotatable bonds. Returns true if a side chain was
   // written. Atom serials are local to the residue (restart at 1), as the
   // AutoDock flex-file convention requires.
   bool write_one_flex_residue(std::ostream &f, mmdb::Residue *residue_p,
                               const coot::dictionary_residue_restraints_t *rest_p) {

      std::string res_name(residue_p->GetResName());

      // Native typing/charges for the whole residue, then keep only the flexible
      // (side-chain + CA) atoms - i.e. those not held rigid in the receptor file.
      std::vector<wr_atom_t> native = native_atoms_for_residue(residue_p, rest_p);
      std::vector<wr_atom_t> fa;     // flexible atoms
      for (unsigned int i=0; i<native.size(); i++)
         if (! rigid_keep_for_flex_residue(native[i].at, residue_p))
            fa.push_back(native[i]);

      int n = static_cast<int>(fa.size());
      if (n == 0) return false;

      // locate CA (the root); record which atoms are heavy and bear a polar H
      int ca = -1;
      std::vector<std::string> name(n);
      std::vector<bool> is_h(n, false);
      for (int i=0; i<n; i++) {
         name[i] = coot::util::remove_whitespace(fa[i].at->name);
         is_h[i] = is_hydrogen_atom(fa[i].at);
         if (name[i] == "CA") ca = i;
      }
      if (ca < 0) {
         std::cout << "WARNING:: write_one_flex_residue(): no CA in "
                   << res_name << " " << residue_p->GetSeqNum() << " - skipped" << std::endl;
         return false;
      }

      // adjacency by distance
      std::vector<std::vector<int> > adj(n);
      for (int i=0; i<n; i++)
         for (int j=i+1; j<n; j++)
            if (bonded_pdbqt(fa[i].at, fa[j].at)) {
               adj[i].push_back(j);
               adj[j].push_back(i);
            }

      // which heavy atoms carry a polar hydrogen (all kept H here are polar)?
      std::vector<bool> heavy_has_polar_h(n, false);
      for (int i=0; i<n; i++)
         if (is_h[i])
            for (unsigned int k=0; k<adj[i].size(); k++)
               heavy_has_polar_h[adj[i][k]] = true;

      std::set<name_pair_t> rot = rest_p ? sidechain_rotatable_bonds(*rest_p)
                                         : std::set<name_pair_t>();

      // Component (BFS) of atom `start`, not crossing the edge (ea,eb).
      auto component_of = [&](int start, int ea, int eb) {
         std::vector<bool> seen(n, false);
         std::vector<int> stack(1, start), out;
         seen[start] = true;
         while (! stack.empty()) {
            int u = stack.back(); stack.pop_back();
            out.push_back(u);
            for (unsigned int k=0; k<adj[u].size(); k++) {
               int v = adj[u][k];
               if ((u==ea && v==eb) || (u==eb && v==ea)) continue;   // skip the cut edge
               if (! seen[v]) { seen[v] = true; stack.push_back(v); }
            }
         }
         return out;
      };

      // A candidate rotatable edge is a REAL cut only if its leaf side (the side
      // away from CA) is worth rotating: >= 2 heavy atoms, or 1 heavy atom that
      // bears a polar hydrogen (hydroxyl/thiol/amine). A lone terminal methyl is
      // therefore not cut - it stays within its parent fragment.
      std::set<std::pair<int,int> > cut_edges;   // (min,max) atom indices
      for (int i=0; i<n; i++) {
         if (is_h[i]) continue;
         for (unsigned int k=0; k<adj[i].size(); k++) {
            int j = adj[i][k];
            if (j < i || is_h[j]) continue;
            if (rot.find(ordered_names(name[i], name[j])) == rot.end()) continue;
            // which endpoint is on the CA side?
            std::vector<int> comp_i = component_of(i, i, j);
            bool i_has_ca = std::find(comp_i.begin(), comp_i.end(), ca) != comp_i.end();
            int leaf = i_has_ca ? j : i;
            std::vector<int> leaf_comp = component_of(leaf, i, j);
            int n_heavy = 0; bool leaf_polar_h = false;
            for (unsigned int t=0; t<leaf_comp.size(); t++) {
               int a = leaf_comp[t];
               if (! is_h[a]) { n_heavy++; if (heavy_has_polar_h[a]) leaf_polar_h = true; }
            }
            if (n_heavy >= 2 || (n_heavy == 1 && leaf_polar_h))
               cut_edges.insert(std::make_pair(std::min(i,j), std::max(i,j)));
         }
      }

      // Fragments = connected components once the cut edges are removed (union-find).
      std::vector<int> uf(n);
      for (int i=0; i<n; i++) uf[i] = i;
      std::function<int(int)> find = [&](int x){ while (uf[x]!=x){ uf[x]=uf[uf[x]]; x=uf[x]; } return x; };
      for (int i=0; i<n; i++)
         for (unsigned int k=0; k<adj[i].size(); k++) {
            int j = adj[i][k];
            std::pair<int,int> e(std::min(i,j), std::max(i,j));
            if (cut_edges.count(e)) continue;
            uf[find(i)] = find(j);
         }

      // Orient each cut: parent endpoint is on the CA side.
      struct cut_t { int parent_atom; int child_atom; int parent_frag; int child_frag; };
      std::vector<cut_t> cuts;
      for (std::set<std::pair<int,int> >::const_iterator it=cut_edges.begin(); it!=cut_edges.end(); ++it) {
         int a = it->first, b = it->second;
         std::vector<int> comp_a = component_of(a, a, b);
         bool a_has_ca = std::find(comp_a.begin(), comp_a.end(), ca) != comp_a.end();
         int pa = a_has_ca ? a : b;
         int ch = a_has_ca ? b : a;
         cut_t c; c.parent_atom = pa; c.child_atom = ch;
         c.parent_frag = find(pa); c.child_frag = find(ch);
         cuts.push_back(c);
      }

      // Emit the tree depth-first from CA's fragment, numbering atoms as written.
      int counter = 0;
      std::vector<int> serial_of(n, 0);
      std::function<void(int, int)> emit_fragment = [&](int frag, int entry_atom) {
         // atoms of this fragment: entry atom first, then the rest in residue order
         std::vector<int> members;
         members.push_back(entry_atom);
         for (int i=0; i<n; i++)
            if (i != entry_atom && find(i) == frag)
               members.push_back(i);
         for (unsigned int m=0; m<members.size(); m++) {
            int idx = members[m];
            serial_of[idx] = ++counter;
            fa[idx].serial = counter;
            f << coot::pdbqt::atom_line(fa[idx], false) << "\n";
         }
         // child branches out of this fragment
         for (unsigned int c=0; c<cuts.size(); c++) {
            if (cuts[c].parent_frag != frag) continue;
            int ps = serial_of[cuts[c].parent_atom];
            f << "BRANCH" << std::setw(4) << ps << std::setw(4) << (counter + 1) << "\n";
            emit_fragment(cuts[c].child_frag, cuts[c].child_atom);
            f << "ENDBRANCH" << std::setw(4) << ps
              << std::setw(4) << serial_of[cuts[c].child_atom] << "\n";
         }
      };

      f << "BEGIN_RES " << res_name << " " << residue_p->GetChainID()
        << std::setw(5) << residue_p->GetSeqNum() << residue_p->GetInsCode() << "\n";
      f << "ROOT\n";
      // The root fragment is CA's fragment, entered at CA.
      // (emit ROOT atoms, then ENDROOT, then its branches - so split the walk.)
      int root_frag = find(ca);
      {
         std::vector<int> members;
         members.push_back(ca);
         for (int i=0; i<n; i++)
            if (i != ca && find(i) == root_frag)
               members.push_back(i);
         for (unsigned int m=0; m<members.size(); m++) {
            int idx = members[m];
            serial_of[idx] = ++counter;
            fa[idx].serial = counter;
            f << coot::pdbqt::atom_line(fa[idx], false) << "\n";
         }
      }
      f << "ENDROOT\n";
      for (unsigned int c=0; c<cuts.size(); c++) {
         if (cuts[c].parent_frag != root_frag) continue;
         int ps = serial_of[cuts[c].parent_atom];
         f << "BRANCH" << std::setw(4) << ps << std::setw(4) << (counter + 1) << "\n";
         emit_fragment(cuts[c].child_frag, cuts[c].child_atom);
         f << "ENDBRANCH" << std::setw(4) << ps
           << std::setw(4) << serial_of[cuts[c].child_atom] << "\n";
      }
      f << "END_RES " << res_name << " " << residue_p->GetChainID()
        << std::setw(5) << residue_p->GetSeqNum() << residue_p->GetInsCode() << "\n";
      return true;
   }
}

int
coot::pdbqt::write_flexible_receptor(mmdb::Manager *mol, coot::protein_geometry &geom, int imol_enc,
                                     const std::vector<coot::residue_spec_t> &flex_residues,
                                     const std::string &rigid_file_name,
                                     const std::string &flex_file_name) {

   if (! mol) return 0;

   std::set<coot::residue_spec_t> flex_specs(flex_residues.begin(), flex_residues.end());

   // 1. rigid file: everything, minus the flexible side chains
   std::ofstream rf(rigid_file_name.c_str());
   if (! rf) {
      std::cout << "WARNING:: write_flexible_receptor(): cannot open " << rigid_file_name << std::endl;
      return 0;
   }
   rf << "REMARK  PDBQT file written by Coot (rigid part of a flexible receptor)\n";
   int serial = 0;
   int n_rigid = write_rigid_atoms(rf, mol, geom, imol_enc, flex_specs, &serial);
   rf << "TER\n";
   rf.close();

   // 2. flex file: one CA-rooted torsion tree per flexible side chain
   std::ofstream ff(flex_file_name.c_str());
   if (! ff) {
      std::cout << "WARNING:: write_flexible_receptor(): cannot open " << flex_file_name << std::endl;
      return 0;
   }
   ff << "REMARK  PDBQT flexible side chains written by Coot\n";
   int n_flex = 0;
   mmdb::Model *model_p = mol->GetModel(1);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (! residue_p) continue;
            if (flex_specs.count(coot::residue_spec_t(residue_p)) == 0) continue;
            std::string res_name(residue_p->GetResName());
            std::pair<bool, coot::dictionary_residue_restraints_t> rp =
               geom.get_monomer_restraints(res_name, imol_enc);
            const coot::dictionary_residue_restraints_t *rest_p = rp.first ? &rp.second : nullptr;
            if (! rest_p)
               std::cout << "WARNING:: write_flexible_receptor(): no dictionary for "
                         << res_name << " - its side chain will have no torsions" << std::endl;
            if (write_one_flex_residue(ff, residue_p, rest_p)) n_flex++;
         }
      }
   }
   ff.close();

   std::cout << "INFO:: write_flexible_receptor(): wrote " << n_rigid << " rigid atoms to "
             << rigid_file_name << " and " << n_flex << " flexible side chain(s) to "
             << flex_file_name << std::endl;
   return n_flex;
}

std::vector<coot::residue_spec_t>
coot::pdbqt::flexible_side_chain_residues_near(mmdb::Manager *mol, double x, double y, double z,
                                               double radius) {

   std::vector<coot::residue_spec_t> v;
   if (! mol) return v;
   double r2 = radius * radius;
   mmdb::Model *model_p = mol->GetModel(1);
   if (! model_p) return v;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int n_res = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<n_res; ires++) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         if (! residue_p) continue;
         if (! residue_is_polymer(residue_p)) continue;
         std::string res_name(residue_p->GetResName());
         // residues with no useful rotatable side chain
         if (res_name == "GLY" || res_name == "ALA" || res_name == "PRO") continue;
         mmdb::PPAtom atoms = 0; int n_atoms = 0;
         residue_p->GetAtomTable(atoms, n_atoms);
         bool within = false;
         for (int i=0; i<n_atoms && !within; i++) {
            double dx = atoms[i]->x - x, dy = atoms[i]->y - y, dz = atoms[i]->z - z;
            if (dx*dx + dy*dy + dz*dz < r2) within = true;
         }
         if (within) v.push_back(coot::residue_spec_t(residue_p));
      }
   }
   return v;
}

// ---------------------------------------------------------------------------
//                     Reading a PDBQT file into an mmdb::Manager
// ---------------------------------------------------------------------------

namespace {

   // Chemical element for an AutoDock atom type (columns 78-79 of a PDBQT file).
   std::string element_for_ad_type(const std::string &t) {
      if (t == "A")  return "C";
      if (t == "C")  return "C";
      if (t == "HD" || t == "HS" || t == "H") return "H";
      if (t == "NA" || t == "NS" || t == "N")  return "N";
      if (t == "OA" || t == "OS" || t == "O")  return "O";
      if (t == "SA" || t == "S")  return "S";
      if (t == "P")  return "P";
      if (t == "F")  return "F";
      if (t == "CL" || t == "Cl") return "Cl";
      if (t == "BR" || t == "Br") return "Br";
      if (t == "I")  return "I";
      if (t == "MG" || t == "Mg") return "Mg";
      if (t == "CA" || t == "Ca") return "Ca";
      if (t == "MN" || t == "Mn") return "Mn";
      if (t == "FE" || t == "Fe") return "Fe";
      if (t == "ZN" || t == "Zn") return "Zn";
      if (t == "CU" || t == "Cu") return "Cu";
      if (t == "SE" || t == "Se") return "Se";
      return t; // best effort
   }

   // Right-justify an element name into the 2-character mmdb convention (" C", "Cl").
   std::string element_2c(const std::string &e) {
      if (e.length() >= 2) return e.substr(0, 2);
      return std::string(" ") + e;
   }

   double to_double(const std::string &s) {
      try { return std::stod(s); } catch (...) { return 0.0; }
   }
}

mmdb::Manager *
coot::pdbqt::read(const std::string &file_name) {

   std::ifstream f(file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: pdbqt::read(): cannot open " << file_name << std::endl;
      return nullptr;
   }

   mmdb::Manager *mol = new mmdb::Manager;
   int h_aff     = mol->RegisterUDReal(mmdb::UDR_MODEL, "vina_affinity");
   int h_rmsd_lb = mol->RegisterUDReal(mmdb::UDR_MODEL, "vina_rmsd_lb");
   int h_rmsd_ub = mol->RegisterUDReal(mmdb::UDR_MODEL, "vina_rmsd_ub");
   int h_inter   = mol->RegisterUDReal(mmdb::UDR_MODEL, "vina_inter");
   int h_intra   = mol->RegisterUDReal(mmdb::UDR_MODEL, "vina_intra");
   int h_unbound = mol->RegisterUDReal(mmdb::UDR_MODEL, "vina_unbound");

   mmdb::Model *model_p = nullptr;   // the current model
   int n_atoms = 0;

   // create model 1 lazily if atoms appear before any MODEL record
   auto ensure_model = [&] () {
      if (! model_p) {
         model_p = new mmdb::Model;
         mol->AddModel(model_p);
      }
   };

   std::string line;
   while (std::getline(f, line)) {
      std::string rec = line.substr(0, 6);
      if (rec == "MODEL " || rec == "MODEL") {
         model_p = new mmdb::Model;
         mol->AddModel(model_p);
         continue;
      }
      if (line.compare(0, 6, "ENDMDL") == 0) { model_p = nullptr; continue; }

      if (line.compare(0, 6, "REMARK") == 0) {
         // REMARK VINA RESULT:   -12.952      0.000      0.000
         std::string::size_type p = line.find("VINA RESULT:");
         if (p != std::string::npos && model_p) {
            std::istringstream iss(line.substr(p + 12));
            double a = 0, b = 0, c = 0;
            iss >> a >> b >> c;
            model_p->PutUDData(h_aff, a);
            model_p->PutUDData(h_rmsd_lb, b);
            model_p->PutUDData(h_rmsd_ub, c);
         } else if (model_p) {
            // energy terms; check the more specific "INTER + INTRA" first
            if (line.find("INTER + INTRA:") != std::string::npos) {
               // (sum, not stored separately)
            } else if ((p = line.find("INTER:")) != std::string::npos) {
               model_p->PutUDData(h_inter, to_double(line.substr(p + 6)));
            } else if ((p = line.find("INTRA:")) != std::string::npos) {
               model_p->PutUDData(h_intra, to_double(line.substr(p + 6)));
            } else if ((p = line.find("UNBOUND:")) != std::string::npos) {
               model_p->PutUDData(h_unbound, to_double(line.substr(p + 8)));
            }
         }
         continue;
      }

      bool is_atom = (line.compare(0, 4, "ATOM") == 0);
      bool is_het  = (line.compare(0, 6, "HETATM") == 0);
      if (! is_atom && ! is_het) continue;   // ROOT/BRANCH/ENDBRANCH/TORSDOF/TER/...
      if (line.length() < 66) continue;

      ensure_model();

      std::string atom_name = line.substr(12, 4);
      std::string alt_loc   = coot::util::remove_whitespace(line.substr(16, 1));
      std::string res_name  = coot::util::remove_whitespace(line.substr(17, 3));
      std::string chain_id  = coot::util::remove_whitespace(line.substr(21, 1));
      int res_seq           = atoi(line.substr(22, 4).c_str());
      std::string ins_code  = coot::util::remove_whitespace(line.substr(26, 1));
      double x = to_double(line.substr(30, 8));
      double y = to_double(line.substr(38, 8));
      double z = to_double(line.substr(46, 8));
      double occ  = (line.length() >= 60) ? to_double(line.substr(54, 6)) : 1.0;
      double bfac = (line.length() >= 66) ? to_double(line.substr(60, 6)) : 0.0;
      std::string ad_type = (line.length() >= 79) ?
                            coot::util::remove_whitespace(line.substr(77, 2)) : std::string();
      std::string element = element_2c(element_for_ad_type(ad_type));
      if (chain_id.empty()) chain_id = "A";

      mmdb::Chain *chain_p = model_p->GetChainCreate(chain_id.c_str(), false);
      mmdb::Residue *res_p = chain_p->GetResidueCreate(res_name.c_str(), res_seq,
                                                       ins_code.c_str(), false);
      mmdb::Atom *at = new mmdb::Atom;
      at->SetCoordinates(x, y, z, occ, bfac);
      at->SetAtomName(atom_name.c_str());
      at->SetElementName(element.c_str());
      strncpy(at->altLoc, alt_loc.c_str(), sizeof(at->altLoc) - 1);
      at->Het = is_het;
      res_p->AddAtom(at);
      n_atoms++;
   }
   f.close();

   if (n_atoms == 0) {
      std::cout << "WARNING:: pdbqt::read(): no atoms read from " << file_name << std::endl;
      delete mol;
      return nullptr;
   }

   mol->FinishStructEdit();
   mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL | mmdb::PDBCLEAN_INDEX);
   std::cout << "INFO:: pdbqt::read(): read " << n_atoms << " atoms in "
             << mol->GetNumberOfModels() << " model(s) from " << file_name << std::endl;
   return mol;
}

std::vector<coot::pdbqt::pose_score_t>
coot::pdbqt::get_scores(mmdb::Manager *mol) {

   std::vector<pose_score_t> scores;
   if (! mol) return scores;

   int h_aff = mol->GetUDDHandle(mmdb::UDR_MODEL, "vina_affinity");
   if (h_aff <= 0) return scores;   // molecule carries no Vina scores

   int h_lb = mol->GetUDDHandle(mmdb::UDR_MODEL, "vina_rmsd_lb");
   int h_ub = mol->GetUDDHandle(mmdb::UDR_MODEL, "vina_rmsd_ub");
   int h_in = mol->GetUDDHandle(mmdb::UDR_MODEL, "vina_inter");
   int h_ia = mol->GetUDDHandle(mmdb::UDR_MODEL, "vina_intra");
   int h_ub2 = mol->GetUDDHandle(mmdb::UDR_MODEL, "vina_unbound");

   auto get = [] (mmdb::Model *m, int handle) -> float {
      mmdb::realtype v = 0.0;
      if (handle > 0 && m->GetUDData(handle, v) == mmdb::UDDATA_Ok) return v;
      return 0.0f;
   };

   for (int i=1; i<=mol->GetNumberOfModels(); i++) {
      mmdb::Model *m = mol->GetModel(i);
      if (! m) continue;
      mmdb::realtype aff = 0.0;
      if (m->GetUDData(h_aff, aff) != mmdb::UDDATA_Ok)
         continue;   // this model has no Vina affinity - skip it
      pose_score_t s;
      s.model_no = i;
      s.affinity = aff;
      s.rmsd_lb  = get(m, h_lb);
      s.rmsd_ub  = get(m, h_ub);
      s.inter    = get(m, h_in);
      s.intra    = get(m, h_ia);
      s.unbound  = get(m, h_ub2);
      scores.push_back(s);
   }
   return scores;
}
