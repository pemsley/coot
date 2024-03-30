/*
 * MoleculesToTriangles/CXXClasses/ColorRule.h
 *
 * Copyright 2009 by Martin Noble, University of Oxford
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
