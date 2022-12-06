/*
 *  ColorScheme.mm
 *  Aesop
 *
 *  Created by Martin Noble on 19/02/2009.
 *  Copyright 2009 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
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
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "white"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None",  "SSE_None"),   "white"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"),  "magenta"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand","SSE_Strand"), "yellow"));
    result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS","NUCLEICACIDS"), "forestgreen"));
    return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorRampChainsScheme(){

   std::string chainIds(" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnpqrstuvwxyz1234567890");
    std::shared_ptr<ColorScheme> result(new ColorScheme());
    
    for (int i=0; i<chainIds.length(); i++){
        std::string selectionString("/*/");
        selectionString.append(1, chainIds[i]);
        selectionString.append("/*.*/*:*");
        auto colorRule = std::shared_ptr<AtomPropertyRampColorRule>(new AtomPropertyRampColorRule());
        colorRule->setRampType(AtomPropertyRampColorRule::ResidueNumber);
        colorRule->setCompoundSelection(std::shared_ptr<CompoundSelection>(new CompoundSelection(selectionString)));
        result->addRule(colorRule);
    }
    return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorChainsScheme(){

   // std::cout << "#######  colorChainsScheme() !" << std::endl;

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
        std::string selectionString("/*/");
        selectionString.append(1, chainIds[i]);
        selectionString.append("/*.*/*:*");
        // std::cout << "### " << selectionString << " " << colorNames[i%nColorNames] << std::endl;
        auto colorRule =
           SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<CompoundSelection>(new CompoundSelection(selectionString)),
                                                        colorNames[i%nColorNames]);
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

