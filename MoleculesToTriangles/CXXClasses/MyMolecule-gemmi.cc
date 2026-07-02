/*
 * MoleculesToTriangles/CXXClasses/MyMolecule-gemmi.cc
 *
 * gemmi-native twin of MyMolecule.cpp (mmdb->gemmi migration).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include <string>
#include <set>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "MyMolecule-gemmi.hh"
#include "DiscreteSegment-gemmi.hh"
#include "DishyBase-gemmi.hh"

#include <gemmi/mmread.hpp>   // read_structure_file
#include <gemmi/pdb.hpp>      // read_pdb_string
#include <gemmi/to_pdb.hpp>   // write_pdb
#include <gemmi/resinfo.hpp>  // find_tabulated_residue
#include "gemmi-sse.hh"       // assign_secondary_structure (DSSP)

namespace {
   // amino-acid test matching the old identifySegments set (standard aa plus the
   // modified residues it special-cased).
   bool is_segment_amino_acid(const std::string &resname) {
      if (gemmi::find_tabulated_residue(resname).is_amino_acid()) return true;
      static const std::set<std::string> extra = {"MSE","TPO","THP","SEP","PTR"};
      return extra.count(resname) > 0;
   }
   bool is_nucleic(const std::string &resname) {
      return gemmi::find_tabulated_residue(resname).is_nucleic_acid();
   }
   bool chain_has_nucleic(const gemmi::Chain &chain) {
      for (const gemmi::Residue &r : chain.residues)
         if (is_nucleic(r.name)) return true;
      return false;
   }
}

coot::m2t::MyMolecule::MyMolecule() : doDraw(true) {
}

coot::m2t::MyMolecule::~MyMolecule() {
}

gemmi::Model *coot::m2t::MyMolecule::getModel() {
   if (structure.models.empty()) return nullptr;
   return &structure.models[0];
}

// "feed gemmi now" boundary.
coot::m2t::MyMolecule::MyMolecule(gemmi::Structure structure_in, int secondaryStructureUsageFlag)
   : doDraw(true), structure(std::move(structure_in)) {
   processCoords(secondaryStructureUsageFlag);
}

coot::m2t::MyMolecule::MyMolecule(const std::string &filePath, int secondaryStructureUsageFlag) : doDraw(true) {
   try {
      structure = gemmi::read_structure_file(filePath);
      processCoords(secondaryStructureUsageFlag);
   } catch (const std::exception &e) {
      std::cout << "error could not read file " << filePath << " : " << e.what() << std::endl;
   }
}

int coot::m2t::MyMolecule::loadCoords(const std::string &fileData, int secondaryStructureUsageFlag) {
   try {
      structure = gemmi::read_pdb_string(fileData, "from-string");
   } catch (const std::exception &e) {
      std::cout << "error could not parse PDB string : " << e.what() << std::endl;
      return 1;
   }
   return processCoords(secondaryStructureUsageFlag);
}

int coot::m2t::MyMolecule::processCoords(int secondaryStructureUsageFlag) {

   identifyBonds();

   // TODO GEMMI: per-atom (united-atom) radii were assigned here by
   // CXXUtils::assignUnitedAtomRadius(mmdb) - that is CXXSurface (Phase 5) and
   // is not yet ported. getAtomRadius() will fall back to vdW until then.

   // Secondary structure: gemmi DSSP (dssp.hpp), stamping Residue::flag with
   // 'H'/'E'/'L' (see gemmi-sse). Replaces mmdb CalcSecStructure / residue->SSE.
   // (USE_HEADER_INFO could later read Structure::helices/sheets; for now DSSP
   //  is used whenever SS is wanted.)
   if (secondaryStructureUsageFlag != DONT_USE)
      coot::m2t::assign_secondary_structure(structure);

   return 0;
}

int coot::m2t::MyMolecule::identifyBonds() {
   bonds = coot::m2t::make_bonds(structure); // also assigns atom.serial over model 0
   return 0;
}

int coot::m2t::MyMolecule::identifySegments(std::vector<DiscreteSegment *> &segments,
                                            const coot::m2t::compound_selection_t &selection) {
   gemmi::Model *model = getModel();
   if (!model) return 0;

   for (gemmi::Chain &chain : model->chains) {
      FCXXCoord lastCoord(-1e30, -1e30, -1e30);

      // ----- protein: break the CA trace where consecutive CAs are > 4.1 A apart
      for (gemmi::Residue &residue : chain.residues) {
         if (!is_segment_amino_acid(residue.name)) continue;
         for (gemmi::Atom &atom : residue.atoms) {
            if (atom.name != "CA") continue;
            gemmi::CRA cra{&chain, &residue, &atom};
            if (!selection.matches(cra)) continue;
            if (!(atom.altloc == '\0' || atom.altloc == 'A' || atom.occ > 0.5)) continue;

            FCXXCoord caPos(atom.pos.x, atom.pos.y, atom.pos.z);
            FCXXCoord diff = caPos - lastCoord;
            auto [ax, ay, az] = getRadiiForResidue(chain.name, residue);
            if (diff.get3DLength() > 4.1) {
               DiscreteSegment *segment = new DiscreteSegment();
               segment->addCalpha(cra, ax, ay, az);
               segments.push_back(segment);
            } else {
               segments.back()->addCalpha(cra, ax, ay, az);
            }
            lastCoord = caPos;
         }
      }

      // ----- nucleic acid: break the C3' trace where consecutive C3's are > 8.0 A apart
      if (chain_has_nucleic(chain)) {
         for (gemmi::Residue &residue : chain.residues) {
            if (!is_nucleic(residue.name)) continue;
            for (gemmi::Atom &atom : residue.atoms) {
               if (atom.name != "C3'") continue;
               gemmi::CRA cra{&chain, &residue, &atom};
               if (!selection.matches(cra)) continue;

               FCXXCoord pos(atom.pos.x, atom.pos.y, atom.pos.z);
               FCXXCoord diff = pos - lastCoord;
               auto [ax, ay, az] = getRadiiForResidue(chain.name, residue);
               if (diff.get3DLength() > 8.0) {
                  DiscreteSegment *segment = new DiscreteSegment();
                  segment->addCalpha(cra, ax, ay, az);
                  segments.push_back(segment);
               } else {
                  segments.back()->addCalpha(cra, ax, ay, az);
               }
               lastCoord = pos;
            }
         }
      }
   }
   return 0;
}

int coot::m2t::MyMolecule::identifyDishyBases(std::map<gemmi::Chain *, DishyBaseContainer_t> &dishy_bases_chain_map,
                                              const coot::m2t::compound_selection_t &selection) {
   gemmi::Model *model = getModel();
   if (!model) return 0;

   for (gemmi::Chain &chain : model->chains) {
      if (!chain_has_nucleic(chain)) continue;

      DishyBaseContainer_t dishybases;

      for (gemmi::Residue &residue : chain.residues) {
         if (residue.atoms.empty()) continue;

         gemmi::CRA cra0{&chain, &residue, &residue.atoms[0]};
         if (!selection.matches(cra0)) continue;

         const std::string &res_name = residue.name;

         std::set<std::string> residue_alt_confs_set;
         for (gemmi::Atom &atom : residue.atoms) {
            std::string ac = (atom.altloc == '\0') ? std::string() : std::string(1, atom.altloc);
            residue_alt_confs_set.insert(ac);
         }

         for (std::size_t ialt = 0; ialt < residue_alt_confs_set.size(); ialt++) {

            std::vector<std::string> ref_base_names;
            std::vector<const gemmi::Atom *> base_atoms;
            if (res_name == "C")  ref_base_names = dishybases.cytidine_base_names;
            if (res_name == "DC") ref_base_names = dishybases.cytidine_base_names;
            if (res_name == "T")  ref_base_names = dishybases.thymine_base_names;
            if (res_name == "DT") ref_base_names = dishybases.thymine_base_names;
            if (res_name == "U")  ref_base_names = dishybases.uracil_base_names;
            if (res_name == "DU") ref_base_names = dishybases.uracil_base_names;
            if (res_name == "A")  ref_base_names = dishybases.adenine_base_names;
            if (res_name == "DA") ref_base_names = dishybases.adenine_base_names;
            if (res_name == "G")  ref_base_names = dishybases.guanine_base_names;
            if (res_name == "DG") ref_base_names = dishybases.guanine_base_names;
            float radius = 2.7;
            if (res_name == "A" || res_name == "DA" || res_name == "G" || res_name == "DG")
               radius = 3.4;

            std::vector<gemmi::CRA> ribose_atoms(5, gemmi::CRA{nullptr, nullptr, nullptr});
            for (gemmi::Atom &atom : residue.atoms) {
               std::string atom_name(atom.name);
               std::string atom_alt_conf = (atom.altloc == '\0') ? std::string() : std::string(1, atom.altloc);
               if (atom_alt_conf.empty() ||
                   (residue_alt_confs_set.find(atom_alt_conf) != residue_alt_confs_set.end())) {
                  if (std::find(ref_base_names.begin(), ref_base_names.end(), atom_name) != ref_base_names.end())
                     base_atoms.push_back(&atom);
                  gemmi::CRA cra{&chain, &residue, &atom};
                  if (atom_name == "O4'") ribose_atoms[0] = cra;
                  if (atom_name == "C1'") ribose_atoms[1] = cra;
                  if (atom_name == "C2'") ribose_atoms[2] = cra;
                  if (atom_name == "C3'") ribose_atoms[3] = cra;
                  if (atom_name == "C4'") ribose_atoms[4] = cra;
               }
            }
            bool have_ribose = true;
            for (std::size_t i = 0; i < 5; i++) if (!ribose_atoms[i].atom) have_ribose = false;
            if (!have_ribose) continue;
            if (base_atoms.size() < 4) continue;

            FCXXCoord ribose_centre;
            for (std::size_t i = 0; i < 5; i++) {
               FCXXCoord pos(ribose_atoms[i].atom->pos.x, ribose_atoms[i].atom->pos.y, ribose_atoms[i].atom->pos.z);
               ribose_centre += pos;
            }
            ribose_centre *= 0.2;
            FCXXCoord base_centre;
            for (std::size_t i = 0; i < base_atoms.size(); i++) {
               FCXXCoord pos(base_atoms[i]->pos.x, base_atoms[i]->pos.y, base_atoms[i]->pos.z);
               base_centre += pos;
            }
            base_centre /= float(base_atoms.size());
            std::vector<FCXXCoord> base_atom_positions(base_atoms.size());
            for (unsigned int i = 0; i < base_atoms.size(); i++)
               base_atom_positions[i] = FCXXCoord(base_atoms[i]->pos.x, base_atoms[i]->pos.y, base_atoms[i]->pos.z);
            DishyPlaneLSQ_t lsq(base_atom_positions);
            FCXXCoord base_normal = lsq.normal();
            DishyBase_t db(base_centre, base_normal, radius, ribose_atoms, ribose_centre);
            dishybases.add(db);
         }
      }
      dishy_bases_chain_map[&chain] = dishybases;
   }
   return 0;
}

FCXXCoord coot::m2t::MyMolecule::getCentre() {
   return centreOfSelectionString(std::string("*/*/*/*"));
}

FCXXCoord coot::m2t::MyMolecule::centreOfSelectionString(const std::string &selectionString) {
   coot::m2t::compound_selection_t selection(selectionString);
   gemmi::Model *model = getModel();
   FCXXCoord sum(0., 0., 0., 0.);
   int n = 0;
   if (model) {
      for (gemmi::Chain &chain : model->chains)
         for (gemmi::Residue &residue : chain.residues)
            for (gemmi::Atom &atom : residue.atoms) {
               gemmi::CRA cra{&chain, &residue, &atom};
               if (selection.matches(cra)) {
                  sum += FCXXCoord(atom.pos.x, atom.pos.y, atom.pos.z);
                  n++;
               }
            }
   }
   if (n > 0) sum /= float(n);
   return sum;
}

void coot::m2t::MyMolecule::writePDB(const std::string &filePath) {
   std::ofstream os(filePath.c_str());
   gemmi::write_pdb(structure, os);
}

void coot::m2t::MyMolecule::setResidueRadii(const std::map<std::tuple<std::string, int, std::string>, std::tuple<float, float, float>> &radii) {
   residueRadii = radii;
}

std::tuple<float, float, float>
coot::m2t::MyMolecule::getRadiiForResidue(const std::string &chain_id, const gemmi::Residue &res) const {
   int res_no = res.seqid.num.value;
   char ic = res.seqid.icode;
   std::string ins_code = (ic == ' ' || ic == '\0') ? std::string() : std::string(1, ic);
   auto key = std::make_tuple(chain_id, res_no, ins_code);
   auto it = residueRadii.find(key);
   if (it != residueRadii.end())
      return it->second;
   return std::make_tuple(1.0f, 1.0f, 1.0f);
}

std::ostream& operator<<(std::ostream& o, const coot::m2t::MyMolecule &myMolecule) {
   std::size_t n_atoms = 0;
   const gemmi::Structure &st = myMolecule.getStructure();
   if (!st.models.empty())
      for (const gemmi::Chain &c : st.models[0].chains)
         for (const gemmi::Residue &r : c.residues)
            n_atoms += r.atoms.size();
   o << "Original name:" << myMolecule.getPDBCode() << "\n" << "nAtoms:" << n_atoms;
   return o;
}
