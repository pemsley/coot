/*
 * MoleculesToTriangles/CXXClasses/SolidColorRule-gemmi.cc
 *
 * gemmi-native twin of SolidColorRule.cpp. The huge colour-name table is NOT
 * duplicated - we reuse the original SolidColorRule's static knownColors map and
 * colorHexToColor (both public, static and mmdb-free).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include <algorithm>
#include "SolidColorRule-gemmi.hh"
#include "SolidColorRule.h"   // original: reuse ::SolidColorRule::knownColors / colorHexToColor

namespace {
   // Resolve a colour name or #hex to an FCXXCoord, reusing the original table.
   // Returns false if the name is unknown (mirrors the original returning a null rule).
   bool resolve_colour(std::string name, FCXXCoord &out) {
      std::transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::toupper);
      auto it = ::SolidColorRule::knownColors.find(name);
      if (it != ::SolidColorRule::knownColors.end()) { out = it->second; return true; }
      if ((name.length() == 7 || name.length() == 9) && name[0] == '#') {
         out = ::SolidColorRule::colorHexToColor(name);
         return true;
      }
      return false;
   }
}

std::shared_ptr<coot::m2t::SolidColorRule>
coot::m2t::SolidColorRule::colorRuleForSelectionAndColor(std::shared_ptr<compound_selection_t> sel, FCXXCoord color) {
   auto result = std::make_shared<SolidColorRule>();
   result->setCompoundSelection(sel);
   result->setColor(color);
   return result;
}

std::shared_ptr<coot::m2t::SolidColorRule>
coot::m2t::SolidColorRule::colorRuleForSelectionAndName(std::shared_ptr<compound_selection_t> sel, std::string _name) {
   FCXXCoord color;
   if (resolve_colour(_name, color))
      return colorRuleForSelectionAndColor(sel, color);
   return nullptr;
}

std::shared_ptr<coot::m2t::SolidColorRule>
coot::m2t::SolidColorRule::colorRuleForSelectionStringAndName(std::string selectionString, std::string _name) {
   auto selection = std::make_shared<compound_selection_t>(selectionString);
   return colorRuleForSelectionAndName(selection, _name);
}

std::shared_ptr<coot::m2t::SolidColorRule>
coot::m2t::SolidColorRule::colorRuleForSelectionStringAndColor(std::string selectionString, FCXXCoord color) {
   auto selection = std::make_shared<compound_selection_t>(selectionString);
   return colorRuleForSelectionAndColor(selection, color);
}
