//
//  MolecularRepresentationInstance.h
//  MoleculesToTriangles
//
//  Created by Martin Noble on 27/01/2017.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

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
