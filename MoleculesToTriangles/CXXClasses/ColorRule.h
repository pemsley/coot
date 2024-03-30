/*
 *  ColorRule.h
 *  Aesop
 *
 *  Created by Martin Noble on 19/02/2009.
 *  Copyright 2009 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */
#ifndef ColorRule_h
#define ColorRule_h

#include <string>
#include <memory>
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

namespace mmdb {
    class Atom;
}
class CompoundSelection;

class ColorRule  {
protected:
    //Modelled
    std::shared_ptr<CompoundSelection> compoundSelection;
    float rank;
    int type;
    //Local
public:
    ColorRule() : compoundSelection(0), rank(-1) {
    };
    enum ColorRuleType {
        Solid,
        AtomPropertyRamp
    }; 
    float getRank() const {
        return rank;
    };
    void setRank(const float &_rank) {
        rank = _rank;
    };
    std::shared_ptr<CompoundSelection> getCompoundSelection() {
        return compoundSelection;
    };    
    void setCompoundSelection(std::shared_ptr<CompoundSelection> _compoundSelection) {
        compoundSelection = _compoundSelection;
    };
    int getType(){
        return type;
    };
    //virtual ~ColorRule() {};
    virtual FCXXCoord colorForAtom (const mmdb::Atom *atom) = 0;
    static bool compareRank(const std::shared_ptr<ColorRule> rule1, const std::shared_ptr<ColorRule> rule2);
    /*
     virtual void prepareForSelectionInMMDB(int handle, mmdb::Manager *mmdb){
    };
     */
};
#endif
//792916
#include "SolidColorRule.h"
#include "AtomPropertyRampColorRule.h"
