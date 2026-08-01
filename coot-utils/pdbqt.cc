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
#include <iostream>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

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
//                     Rigid receptor / whole-molecule PDBQT
// ---------------------------------------------------------------------------

int
coot::pdbqt::write_receptor(mmdb::Manager *mol, coot::protein_geometry &geom, int imol_enc,
                            const std::string &file_name) {

   int n_written = 0;
   if (! mol) return 0;

   std::ofstream f(file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: write_receptor(): cannot open " << file_name << std::endl;
      return 0;
   }

   f << "REMARK  PDBQT file written by Coot (rigid receptor)\n";

   int serial = 0;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (! residue_p) continue;
            std::string res_name(residue_p->GetResName());
            if (residue_is_water(res_name)) continue; // receptor prep drops waters
            std::pair<bool, coot::dictionary_residue_restraints_t> rp =
               geom.get_monomer_restraints(res_name, imol_enc);
            const coot::dictionary_residue_restraints_t *rest_p = rp.first ? &rp.second : nullptr;
            std::vector<wr_atom_t> atoms = native_atoms_for_residue(residue_p, rest_p);
            bool is_het = ! residue_is_polymer(residue_p);
            for (unsigned int i=0; i<atoms.size(); i++) {
               atoms[i].serial = ++serial;
               f << coot::pdbqt::atom_line(atoms[i], is_het) << "\n";
               n_written++;
            }
         }
      }
   }
   f << "TER\n";
   f.close();
   std::cout << "INFO:: write_receptor(): wrote " << n_written
             << " atoms to " << file_name << std::endl;
   return n_written;
}
