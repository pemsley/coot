/*
 *  SolidColorRule.h
 *  Aesop
 *
 *  Created by Martin Noble on 04/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
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
