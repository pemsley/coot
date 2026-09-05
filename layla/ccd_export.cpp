/* layla/ccd_export.cpp
 *
 * Copyright 2026 by Global Phasing Ltd.
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

// C++ port of python/convert-pkl-to-mmcif.py: write a CCD-style mmCIF
// for acedrg directly from the RDKit molecule held by the Layla canvas.

#include "ccd_export.hpp"

#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Descriptors/MolDescriptors.h>
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>
#include <rdkit/GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/MolDraw2D/MolDraw2D.h>
#include <rdkit/GraphMol/MolDraw2D/MolDraw2DUtils.h>

#include <gemmi/cifdoc.hpp>
#include <gemmi/to_cif.hpp>

namespace {

   std::string fmt3(double v) {
      char buf[32];
      std::snprintf(buf, sizeof(buf), "%.3f", v);
      return std::string(buf);
   }

   // _chem_comp_bond.value_order - the CCD convention: kekulized orders
   std::string ccd_value_order(RDKit::Bond::BondType bt) {
      switch (bt) {
         case RDKit::Bond::SINGLE:   return "SING";
         case RDKit::Bond::DOUBLE:   return "DOUB";
         case RDKit::Bond::TRIPLE:   return "TRIP";
         case RDKit::Bond::AROMATIC: return "AROM";
         default:                    return "SING";
      }
   }

   // the RDKit python enum names, as pdbeccdutils writes them in the
   // _pdbe_chem_comp_bond_depiction loop
   std::string bond_type_name(RDKit::Bond::BondType bt) {
      switch (bt) {
         case RDKit::Bond::SINGLE:      return "SINGLE";
         case RDKit::Bond::DOUBLE:      return "DOUBLE";
         case RDKit::Bond::TRIPLE:      return "TRIPLE";
         case RDKit::Bond::AROMATIC:    return "AROMATIC";
         case RDKit::Bond::DATIVE:      return "DATIVE";
         case RDKit::Bond::UNSPECIFIED: return "UNSPECIFIED";
         default:                       return "SINGLE";
      }
   }

   std::string bond_dir_name(RDKit::Bond::BondDir bd) {
      switch (bd) {
         case RDKit::Bond::NONE:         return "NONE";
         case RDKit::Bond::BEGINWEDGE:   return "BEGINWEDGE";
         case RDKit::Bond::BEGINDASH:    return "BEGINDASH";
         case RDKit::Bond::ENDDOWNRIGHT: return "ENDDOWNRIGHT";
         case RDKit::Bond::ENDUPRIGHT:   return "ENDUPRIGHT";
         case RDKit::Bond::EITHERDOUBLE: return "EITHERDOUBLE";
         default:                        return "UNKNOWN";
      }
   }

} // namespace

std::string
coot::layla::make_acedrg_input_mmcif(const RDKit::ROMol &mol_2d_in, const std::string &comp_id_in) {

   std::string comp_id = comp_id_in;
   if (comp_id.empty()) comp_id = "LIG";

   if (mol_2d_in.getNumAtoms() == 0)
      throw std::runtime_error("make_acedrg_input_mmcif(): molecule has no atoms");

   // The 2D sketch: Layla's own hand-drawn layout, used verbatim for the
   // depiction loops. If the molecule somehow has no conformer, lay one out.
   RDKit::RWMol mol_2d(mol_2d_in);
   if (mol_2d.getNumConformers() == 0)
      RDDepict::compute2DCoords(mol_2d);

   // The working molecule: add Hs, build a genuine 3D conformer with ETKDG
   // (acedrg regenerates its own 3D anyway, but this makes the *_ideal
   // coordinates sane), light MMFF clean-up.
   RDKit::RWMol mol(mol_2d_in);
   mol.clearConformers();
   RDKit::MolOps::addHs(mol); // add Hs *before* embedding
   RDKit::DGeomHelpers::EmbedParameters params(RDKit::DGeomHelpers::ETKDGv3);
   params.randomSeed = 1;
   if (RDKit::DGeomHelpers::EmbedMolecule(mol, params) != 0) {
      // fallback for awkward molecules: relax the strict experimental-torsion terms
      params.useRandomCoords = true;
      if (RDKit::DGeomHelpers::EmbedMolecule(mol, params) != 0)
         throw std::runtime_error("make_acedrg_input_mmcif(): 3D embedding failed for " + comp_id);
   }
   try {
      RDKit::MMFF::MMFFOptimizeMolecule(mol);
   }
   catch (const std::exception &e) {
      // no MMFF parameters for some atom - the raw ETKDG coordinates will do
   }
   const RDKit::Conformer &conf = mol.getConformer();

   // Match the PDB-CCD convention: aromatic ring bonds carry a kekulized
   // order (SING/DOUB) plus pdbx_aromatic_flag = Y - not a literal "AROM".
   // clearAromaticFlags = false keeps getIsAromatic() true for the flag.
   RDKit::MolOps::Kekulize(mol, false);

   // ---- atom names ------------------------------------------------------
   // Preserving the Coot/Layla atom names is the point of this pipeline, so
   // names carried on the molecule take priority. Atoms without a name (the
   // added Hs) get a generated element+counter name; CCD atom_id must be
   // unique, so generated names are checked against every name in use.
   unsigned int n_atoms = mol.getNumAtoms();
   std::vector<std::string> names(n_atoms);
   std::set<std::string> used;
   for (unsigned int i = 0; i < n_atoms; i++) {
      const RDKit::Atom *at = mol.getAtomWithIdx(i);
      std::string nm;
      if (at->getPropIfPresent<std::string>("name", nm)) {
         // strip whitespace
         auto b = nm.find_first_not_of(" \t");
         auto e = nm.find_last_not_of(" \t");
         nm = (b == std::string::npos) ? std::string() : nm.substr(b, e - b + 1);
         if (! nm.empty()) {
            if (used.count(nm))
               throw std::runtime_error("make_acedrg_input_mmcif(): duplicate atom name '" + nm +
                                        "' on the input molecule; names must be unique");
            names[i] = nm;
            used.insert(nm);
         }
      }
   }
   std::map<std::string, unsigned int> elem_counts;
   for (unsigned int i = 0; i < n_atoms; i++) {
      if (! names[i].empty()) continue;
      const RDKit::Atom *at = mol.getAtomWithIdx(i);
      std::string sym = at->getSymbol();
      for (auto &c : sym) c = toupper(c);
      std::string cand;
      do {
         elem_counts[sym] += 1;
         cand = sym + std::to_string(elem_counts[sym]);
      } while (used.count(cand));
      names[i] = cand;
      used.insert(cand);
   }

   // ---- build the mmCIF -------------------------------------------------
   gemmi::cif::Document doc;
   gemmi::cif::Block &blk = doc.add_new_block(comp_id);
   blk.set_pair("_chem_comp.id", comp_id);
   blk.set_pair("_chem_comp.three_letter_code", comp_id);
   blk.set_pair("_chem_comp.name", gemmi::cif::quote(RDKit::MolToSmiles(mol)));
   blk.set_pair("_chem_comp.type", "non-polymer");
   blk.set_pair("_chem_comp.formula", gemmi::cif::quote(RDKit::Descriptors::calcMolFormula(mol)));

   gemmi::cif::Loop &atom_loop =
      blk.init_loop("_chem_comp_atom.", {"comp_id", "atom_id", "type_symbol", "charge",
                                         "pdbx_model_Cartn_x_ideal",
                                         "pdbx_model_Cartn_y_ideal",
                                         "pdbx_model_Cartn_z_ideal"});
   for (unsigned int i = 0; i < n_atoms; i++) {
      const RDKit::Atom *at = mol.getAtomWithIdx(i);
      const RDGeom::Point3D &p = conf.getAtomPos(i);
      atom_loop.add_row({comp_id, gemmi::cif::quote(names[i]), at->getSymbol(),
                         std::to_string(at->getFormalCharge()),
                         fmt3(p.x), fmt3(p.y), fmt3(p.z)});
   }

   gemmi::cif::Loop &bond_loop =
      blk.init_loop("_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2",
                                         "value_order", "pdbx_aromatic_flag"});
   for (const auto *b : mol.bonds()) {
      bond_loop.add_row({comp_id,
                         gemmi::cif::quote(names[b->getBeginAtomIdx()]),
                         gemmi::cif::quote(names[b->getEndAtomIdx()]),
                         ccd_value_order(b->getBondType()),
                         b->getIsAromatic() ? "Y" : "N"});
   }

   // ---- 2D depiction loops (PDBe extension) -----------------------------
   // Written from Layla's hand-drawn layout, not a recomputed one. mol_2d
   // holds only the drawn heavy atoms (no explicit H), and its atom indices
   // line up with the heavy atoms of the working molecule, so names[] can
   // be reused here.
   const RDKit::Conformer &conf2d = mol_2d.getConformer();
   gemmi::cif::Loop &atom_dep =
      blk.init_loop("_pdbe_chem_comp_atom_depiction.", {"comp_id", "atom_id", "element",
                                                        "model_Cartn_x", "model_Cartn_y",
                                                        "pdbx_ordinal"});
   for (unsigned int i = 0; i < mol_2d.getNumAtoms(); i++) {
      const RDKit::Atom *at = mol_2d.getAtomWithIdx(i);
      const RDGeom::Point3D &p = conf2d.getAtomPos(i);
      atom_dep.add_row({comp_id, gemmi::cif::quote(names[i]), at->getSymbol(),
                        fmt3(p.x), fmt3(p.y), std::to_string(i + 1)});
   }

   // Wedge/hash bond directions for the drawing, derived from the 2D layout
   // + chirality (the same prepareMolForDrawing pdbeccdutils uses). Bonds
   // to H are omitted.
   RDKit::RWMol drawmol(mol_2d);
   try {
      RDKit::MolDraw2DUtils::prepareMolForDrawing(drawmol, true, true, true); // kekulize, addChiralHs, wedgeBonds
   }
   catch (const std::exception &e) {
      drawmol = RDKit::RWMol(mol_2d);
      try {
         RDKit::MolDraw2DUtils::prepareMolForDrawing(drawmol, true, false, false);
      }
      catch (const std::exception &e2) {
         drawmol = RDKit::RWMol(mol_2d);
         RDKit::MolOps::Kekulize(drawmol, true);
      }
   }

   gemmi::cif::Loop &bond_dep =
      blk.init_loop("_pdbe_chem_comp_bond_depiction.", {"comp_id", "atom_id_1", "atom_id_2",
                                                        "value_order", "bond_dir", "pdbx_ordinal"});
   unsigned int dep_ordinal = 0;
   for (const auto *b : drawmol.bonds()) {
      const RDKit::Atom *a1 = b->getBeginAtom();
      const RDKit::Atom *a2 = b->getEndAtom();
      if (a1->getAtomicNum() == 1 || a2->getAtomicNum() == 1) continue;
      if (a1->getIdx() >= mol_2d.getNumAtoms() || a2->getIdx() >= mol_2d.getNumAtoms())
         continue; // atoms added to drawmol by prepareMolForDrawing have no names[] entry
      dep_ordinal++;
      bond_dep.add_row({comp_id,
                        gemmi::cif::quote(names[a1->getIdx()]),
                        gemmi::cif::quote(names[a2->getIdx()]),
                        bond_type_name(b->getBondType()),
                        bond_dir_name(b->getBondDir()),
                        std::to_string(dep_ordinal)});
   }

   std::ostringstream ss;
   gemmi::cif::write_cif_to_stream(ss, doc);
   return ss.str();
}
