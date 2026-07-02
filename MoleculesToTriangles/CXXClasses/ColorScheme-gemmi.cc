/*
 * MoleculesToTriangles/CXXClasses/ColorScheme-gemmi.cc
 *
 * gemmi-native twin of ColorScheme.cpp - the colour-scheme factories.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "ColorScheme-gemmi.hh"
#include "SolidColorRule-gemmi.hh"
#include "AtomPropertyRampColorRule-gemmi.hh"
#include <algorithm>

using coot::m2t::ColorScheme;
using coot::m2t::SolidColorRule;
using coot::m2t::AtomPropertyRampColorRule;
using coot::m2t::compound_selection_t;

std::shared_ptr<ColorScheme> ColorScheme::colorByElementScheme() {
   std::shared_ptr<ColorScheme> result(new ColorScheme());
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*.*/*:*"), "#dddddd"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[C].*"), "#a0d78c"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[N].*"), "royalblue"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[P].*"), "orange"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[S].*"), "yellow"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[O].*"), "tomato"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[H].*"), "#dddddd"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[F].*"), "#aaccaa"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*/*[D].*"), "pink"));
   return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorBySecondaryScheme() {
   std::shared_ptr<ColorScheme> result(new ColorScheme());
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*.*/*:*"), "#666666"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("SSE_None"),   "#666666"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("SSE_Helix"),  "#aa2222"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("SSE_Strand"), "#aaaa22"));
   result->addRule(SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("NUCLEICACIDS"), "forestgreen"));
   return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorRampChainsScheme() {
   std::string chainIds(" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnpqrstuvwxyz1234567890");
   std::shared_ptr<ColorScheme> result(new ColorScheme());
   for (unsigned int i=0; i<chainIds.length(); i++) {
      std::string selectionString("/*/");
      selectionString.append(1, chainIds[i]);
      selectionString.append("/*.*/*:*");
      auto colorRule = std::make_shared<AtomPropertyRampColorRule>();
      colorRule->setRampType(AtomPropertyRampColorRule::ResidueNumber);
      colorRule->setCompoundSelection(std::make_shared<compound_selection_t>(selectionString));
      result->addRule(colorRule);
   }
   return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorChainsScheme() {
   std::shared_ptr<ColorScheme> result(new ColorScheme());
   std::string colorNames[] = {
      "Salmon", "SandyBrown", "Burlywood", "Goldenrod",
      "tomato", "limegreen", "royalblue", "gold", "aquamarine", "maroon",
      "lightsalmon", "deeppink", "brown"
   };
   std::string chainIds(" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnpqrstuvwxyz1234567890");
   int nColorNames = sizeof(colorNames)/sizeof(std::string);
   for (unsigned int i=0; i<chainIds.length(); i++) {
      std::string selectionString = std::string("//");
      selectionString.append(1, chainIds[i]);
      std::string c = colorNames[i%nColorNames];
      auto colorRule = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>(selectionString), c);
      result->addRule(colorRule);
   }
   return result;
}

std::shared_ptr<ColorScheme>
ColorScheme::colorChainsSchemeWithColourRules(const std::vector<std::pair<std::string, std::string> > &colour_rules) {
   std::shared_ptr<ColorScheme> result(new ColorScheme());
   for (const auto &cr : colour_rules) {
      auto colorRule = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>(cr.first), cr.second);
      result->addRule(colorRule);
   }
   return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorBFactorScheme() {
   std::shared_ptr<ColorScheme> result(new ColorScheme());
   auto colorRule = std::make_shared<AtomPropertyRampColorRule>();
   colorRule->setStartValue(10);
   colorRule->setStartRGB(FCXXCoord(0., 0., 1.));
   colorRule->setEndRGB(FCXXCoord(1., 0., 0.));
   colorRule->setEndValue(60);
   colorRule->setRampType(AtomPropertyRampColorRule::BFactor);
   colorRule->setCompoundSelection(std::make_shared<compound_selection_t>("/*/*/*.*/*:*"));
   result->addRule(colorRule);
   return result;
}

std::shared_ptr<ColorScheme> ColorScheme::colorSchemeForColorName(const std::string name) {
   std::shared_ptr<ColorScheme> result(new ColorScheme());
   std::string s(name);
   std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) toupper);
   auto colorRule = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<compound_selection_t>("/*/*/*.*/*:*"), s);
   result->addRule(colorRule);
   return result;
}
