/*
 * MoleculesToTriangles/CXXClasses/SolidColorRule.h
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
#ifndef SolidColorRule_h
#define SolidColorRule_h

#include "ColorRule.h"
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include <map>
#include <string>
#include <memory>

class SolidColorRule : public ColorRule {
private:
    FCXXCoord color;
public:
    static std::pair<std::string, FCXXCoord  > nameColorPairs[];
    static std::map<std::string, FCXXCoord  > knownColors;
    SolidColorRule () {
        type = Solid;
    };
    FCXXCoord getColor() const{
        return color;
    };
    void setColor(const FCXXCoord _color){
        color = _color;
    }
    virtual FCXXCoord colorForAtom (const mmdb::Atom* atom) {
		return color;
    };
    static std::shared_ptr<SolidColorRule> colorRuleForSelectionStringAndName(std::string selectionString, std::string _name);
    static std::shared_ptr<SolidColorRule> colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection> selection, std::string _name);
    static std::shared_ptr<SolidColorRule> colorRuleForSelectionAndColor(std::shared_ptr<CompoundSelection> selection, FCXXCoord colorCoord);
    static std::shared_ptr<SolidColorRule> colorRuleForSelectionStringAndColor(std::string selectionString, FCXXCoord color);
    std::shared_ptr<FCXXCoord> colorForName(std::string _name);
    static FCXXCoord colorHexToColor(const std::string &colourHex_in);
};

#endif
