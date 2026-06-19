/*
 * MoleculesToTriangles/CXXClasses/MyMolecule-gemmi.hh
 *
 * gemmi-native twin of MyMolecule.h (mmdb->gemmi migration). Holds a
 * gemmi::Structure + covalent bond list. Lives alongside the original;
 * type is coot::m2t::MyMolecule. "feed gemmi now": construct from a
 * gemmi::Structure.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef MyMolecule_gemmi_hh
#define MyMolecule_gemmi_hh

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <tuple>

#include <gemmi/model.hpp>   // gemmi::Structure, Model, Chain, Residue, Atom, CRA

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include "MoleculesToTriangles/CXXClasses/gemmi-bonds.hh"
#include "MoleculesToTriangles/CXXClasses/gemmi-selection.hh"

namespace coot {
   namespace m2t {

      class DiscreteSegment;
      class DishyBaseContainer_t;

      enum SecondaryStructureUsageType { USE_HEADER_INFO, DONT_USE, CALC_SECONDARY_STRUCTURE };

      class MyMolecule {
      private:
         bool doDraw;
         std::map<std::tuple<std::string, int, std::string>, std::tuple<float, float, float>> residueRadii;
         int processCoords(int secondaryStructureUsageFlag);
      public:
         gemmi::Structure structure;
         std::vector<coot::m2t::bond_t> bonds;   // covalent bonds (atom.serial pairs, model 0)
         std::string PDBCode;

         MyMolecule();
         MyMolecule(const std::string &filePath, int secondaryStructureUsageFlag);
         // "feed gemmi now" entry point - callers pass a gemmi::Structure.
         MyMolecule(gemmi::Structure structure_in, int secondaryStructureUsageFlag);
         ~MyMolecule();

         gemmi::Structure &getStructure() { return structure; }
         const gemmi::Structure &getStructure() const { return structure; }
         gemmi::Model *getModel();                 // model 0, or nullptr

         int loadCoords(const std::string &fileData, int secondaryStructureUsageFlag);
         int identifySegments(std::vector<DiscreteSegment *> &segments,
                              const coot::m2t::compound_selection_t &selection);
         int identifyDishyBases(std::map<gemmi::Chain *, DishyBaseContainer_t> &dishy_bases_chain_map,
                                const coot::m2t::compound_selection_t &selection);
         int identifyBonds();                      // populates `bonds` via gemmi-bonds
         FCXXCoord getCentre();
         FCXXCoord centreOfSelectionString(const std::string &selectionString);

         void setPDBCode(std::string _code) { PDBCode = _code; }
         const std::string getPDBCode() const { return PDBCode; }
         bool getDoDraw() const { return doDraw; }
         void setDoDraw(const bool &yesOrNo) { doDraw = yesOrNo; }
         void writePDB(const std::string &filePath);
         void setResidueRadii(const std::map<std::tuple<std::string, int, std::string>, std::tuple<float, float, float>> &radii);
         std::tuple<float, float, float> getRadiiForResidue(const std::string &chain_id, const gemmi::Residue &res) const;

         static std::shared_ptr<MyMolecule> create(std::string filePathString, int secondaryStructureUsageFlag) {
            return std::make_shared<MyMolecule>(filePathString, secondaryStructureUsageFlag);
         }
         static std::shared_ptr<MyMolecule> createFromString(std::string contentsAsString) {
            auto newMolecule = std::make_shared<MyMolecule>();
            newMolecule->loadCoords(contentsAsString, CALC_SECONDARY_STRUCTURE);
            return newMolecule;
         }
      };
   }
}
std::ostream& operator<<(std::ostream&, const coot::m2t::MyMolecule &);

#endif
