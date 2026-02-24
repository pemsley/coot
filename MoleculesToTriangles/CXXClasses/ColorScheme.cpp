/*
 * MoleculesToTriangles/CXXClasses/ColorScheme.cpp
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

#include "ColorScheme.h"
#include <algorithm>

std::shared_ptr<ColorScheme> ColorScheme::colorByElementScheme()
{
    std::shared_ptr<ColorScheme> result(new ColorScheme());
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*.*/*:*","OtherElements")),  "#dddddd"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[C].*","/*/*/*/*[C].*")), "#a0d78c"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[N].*","/*/*/*/*[N].*")), "royalblue"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[P].*","/*/*/*/*[P].*")), "orange"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[S].*","/*/*/*/*[S].*")), "yellow"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[O].*","/*/*/*/*[O].*")), "tomato"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[H].*","/*/*/*/*[H].*")), "#dddddd"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[F].*","/*/*/*/*[F].*")), "#aaccaa"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*/*[D].*","/*/*/*/*[D].*")), "pink"));
    return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorBySecondaryScheme()
{
    std::shared_ptr<ColorScheme> result(new ColorScheme());
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "#666666"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None",  "SSE_None"),   "#666666"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"),  "#aa2222"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand","SSE_Strand"), "#aaaa22"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS","NUCLEICACIDS"), "forestgreen"));
    return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorRampChainsScheme(){

   std::cout << "####### M2T colorRampChainsScheme() !!" << std::endl;

   std::string chainIds(" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnpqrstuvwxyz1234567890");
    std::shared_ptr<ColorScheme> result(new ColorScheme());
    
    for (unsigned int i=0; i<chainIds.length(); i++){
        std::string selectionString("/*/");
        selectionString.append(1, chainIds[i]);
        selectionString.append("/*.*/*:*");
        // std::cout << "colorRampChainsScheme() selectionString " << selectionString << std::endl;
        auto colorRule = std::shared_ptr<AtomPropertyRampColorRule>(new AtomPropertyRampColorRule());
        colorRule->setRampType(AtomPropertyRampColorRule::ResidueNumber);
        colorRule->setCompoundSelection(std::shared_ptr<CompoundSelection>(new CompoundSelection(selectionString)));
        // std::cout << "colorRampChainsScheme() selectionString " << colorRule << std::endl;
        result->addRule(colorRule);
    }
    return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorChainsScheme(){

    std::shared_ptr<ColorScheme> result(new ColorScheme());
    std::string colorNames[] = {
                                //    "RED","GREEN","BLUE","CYAN","MAGENTA","YELLOW","WHITE"
                                "Salmon", "Sandy Brown",
                                "Burlywood", "Goldenrod",
                                "tomato", "limegreen", "royalblue", "gold", "aquamarine", "maroon",
                                "lightsalmon", "deeppink", "brown"
    };
    std::string chainIds(" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnpqrstuvwxyz1234567890");
    int nColorNames = sizeof(colorNames)/sizeof(std::string);

    for (unsigned int i=0; i<chainIds.length(); i++){
#if 0
        std::string selectionString("/*/");
        selectionString.append(1, chainIds[i]);
        selectionString.append("/*.*/*:*");
#endif
        std::string selectionString = std::string("//");
        selectionString.append(1, chainIds[i]);

        // std::cout << "### " << selectionString << " " << colorNames[i%nColorNames] << std::endl;
        std::string c =   colorNames[i%nColorNames];
        std::shared_ptr<SolidColorRule> colorRule =
           SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection(selectionString)), c);
        result->addRule(colorRule);
    }
    return result;
}

// static
std::shared_ptr<ColorScheme> ColorScheme::colorChainsSchemeWithColourRules(const std::vector<std::pair<std::string, std::string> > &colour_rules) {

   std::shared_ptr<ColorScheme> result(new ColorScheme());
   std::vector<std::pair<std::string, std::string> >::const_iterator it;
   for (it=colour_rules.begin(); it!=colour_rules.end(); ++it) {
      const std::string &selectionString = it->first;
      const std::string &colorName       = it->second;
      auto colorRule = SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection(selectionString)), colorName);
      result->addRule(colorRule);
   }
   return result;
}


std::shared_ptr<ColorScheme> ColorScheme::colorBFactorScheme(){
    std::shared_ptr<ColorScheme> result(new ColorScheme());
    auto colorRule = std::shared_ptr<AtomPropertyRampColorRule>(new AtomPropertyRampColorRule());
    colorRule->setStartValue(10);
    colorRule->setStartRGB(FCXXCoord (0.,0.,1.));
    colorRule->setEndRGB(FCXXCoord (1.,0.,0.));
    colorRule->setEndValue(60);
    colorRule->setRampType(AtomPropertyRampColorRule::BFactor);
    colorRule->setCompoundSelection(std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*.*/*:*")));
    result->addRule(colorRule);
    return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorSchemeForColorName(const std::string name){     std::shared_ptr<ColorScheme> result(new ColorScheme());
    std::string s(name);
    std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) toupper);
    std::string selectionString("/*/*/*.*/*:*");
    auto colorRule =
    SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection(selectionString)),s);
    result->addRule(colorRule);
    return result;
}

