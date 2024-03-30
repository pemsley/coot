/*
 * MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.h
 *
 * Copyright 2017 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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

#ifndef MolecularRepresentationInstance_h
#define MolecularRepresentationInstance_h
#include <memory>

#include <stdio.h>
#include "MolecularRepresentation.h"
#include "RepresentationInstance.h"
#include "MyMolecule.h"
#include "ColorScheme.h"
#include "CompoundSelection.h"

class MolecularRepresentationInstance : public RepresentationInstance {
public:
    MolecularRepresentationInstance() : RepresentationInstance() {
        auto molrep = std::shared_ptr<MolecularRepresentation>(new MolecularRepresentation());
        setRepresentation(molrep);
    };
    MolecularRepresentationInstance(std::shared_ptr<MyMolecule> _myMolecule, std::shared_ptr<ColorScheme> _colorScheme, std::shared_ptr<CompoundSelection> _compoundSelection, std::string _renderStyle) : RepresentationInstance(){
        auto molrep = std::shared_ptr<MolecularRepresentation>(new MolecularRepresentation(_myMolecule, _colorScheme, _compoundSelection, _renderStyle));
        setRepresentation(molrep);
    };
    MolecularRepresentationInstance(std::shared_ptr<MyMolecule> _myMolecule, std::shared_ptr<ColorScheme> _colorScheme, std::string _compoundSelectionString, std::string _renderStyle) : RepresentationInstance(){
        auto _compoundSelection=CompoundSelection::create(_compoundSelectionString,"Anon");
        auto molrep = std::shared_ptr<MolecularRepresentation>(new MolecularRepresentation(_myMolecule, _colorScheme, _compoundSelection, _renderStyle));
        setRepresentation(molrep);
    };
    /*
    static std::shared_ptr<MolecularRepresentationInstance> create(std::shared_ptr<MyMolecule> _myMolecule, std::shared_ptr<ColorScheme> _colorScheme, std::shared_ptr<CompoundSelection> _compoundSelection, std::string _renderStyle){
        return std::shared_ptr<MolecularRepresentationInstance>(new MolecularRepresentationInstance(_myMolecule, _colorScheme, _compoundSelection, _renderStyle));
    };
     */
    static std::shared_ptr<MolecularRepresentationInstance> create(std::shared_ptr<MyMolecule> _myMolecule, std::shared_ptr<ColorScheme> _colorScheme, std::string _compoundSelectionString, std::string _renderStyle){
        return std::shared_ptr<MolecularRepresentationInstance>(new MolecularRepresentationInstance(_myMolecule, _colorScheme, _compoundSelectionString, _renderStyle));
    };
    std::shared_ptr<MolecularRepresentation> getRepresentation(){
        return std::dynamic_pointer_cast<MolecularRepresentation>(representation);
    };
};
#endif /* MolecularRepresentationInstance_h */
