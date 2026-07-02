/*
 * MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance-gemmi.hh
 *
 * gemmi-native twin of MolecularRepresentationInstance.h. Thin wrapper over the
 * (mmdb-free) RepresentationInstance base, holding a coot::m2t::MolecularRepresentation.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef MolecularRepresentationInstance_gemmi_hh
#define MolecularRepresentationInstance_gemmi_hh

#include <memory>
#include <string>

#include "RepresentationInstance.h"
#include "MolecularRepresentation-gemmi.hh"
#include "MyMolecule-gemmi.hh"
#include "ColorScheme-gemmi.hh"
#include "MoleculesToTriangles/CXXClasses/gemmi-selection.hh"

namespace coot {
   namespace m2t {

      class MolecularRepresentationInstance : public RepresentationInstance {
      public:
         MolecularRepresentationInstance() : RepresentationInstance() {
            auto molrep = std::make_shared<MolecularRepresentation>();
            setRepresentation(molrep);
         }
         MolecularRepresentationInstance(std::shared_ptr<MyMolecule> _myMolecule,
                                         std::shared_ptr<ColorScheme> _colorScheme,
                                         std::shared_ptr<compound_selection_t> _compoundSelection,
                                         std::string _renderStyle) : RepresentationInstance() {
            auto molrep = std::make_shared<MolecularRepresentation>(_myMolecule, _colorScheme,
                                                                    _compoundSelection, _renderStyle);
            setRepresentation(molrep);
         }
         MolecularRepresentationInstance(std::shared_ptr<MyMolecule> _myMolecule,
                                         std::shared_ptr<ColorScheme> _colorScheme,
                                         std::string _compoundSelectionString,
                                         std::string _renderStyle) : RepresentationInstance() {
            auto _compoundSelection = std::make_shared<compound_selection_t>(_compoundSelectionString);
            auto molrep = std::make_shared<MolecularRepresentation>(_myMolecule, _colorScheme,
                                                                    _compoundSelection, _renderStyle);
            setRepresentation(molrep);
         }
         static std::shared_ptr<MolecularRepresentationInstance>
         create(std::shared_ptr<MyMolecule> _myMolecule, std::shared_ptr<ColorScheme> _colorScheme,
                std::string _compoundSelectionString, std::string _renderStyle) {
            return std::make_shared<MolecularRepresentationInstance>(_myMolecule, _colorScheme,
                                                                     _compoundSelectionString, _renderStyle);
         }
         std::shared_ptr<MolecularRepresentation> getRepresentation() {
            // use the base's public accessor (representation is private; the base
            // befriends only the original global MolecularRepresentationInstance).
            return std::dynamic_pointer_cast<MolecularRepresentation>(
               RepresentationInstance::getRepresentation());
         }
      };
   }
}

#endif // MolecularRepresentationInstance_gemmi_hh
