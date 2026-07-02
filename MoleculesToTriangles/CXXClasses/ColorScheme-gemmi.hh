/*
 * MoleculesToTriangles/CXXClasses/ColorScheme-gemmi.hh
 *
 * gemmi-native twin of ColorScheme.h. colorForAtom takes a gemmi::CRA and tests
 * each rule's compound_selection_t::matches(cra) directly - no mmdb selection
 * handles (prepareForMMDB / freeSelectionHandles are gone). Type is
 * coot::m2t::ColorScheme.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef ColorScheme_gemmi_hh
#define ColorScheme_gemmi_hh

#include <list>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>

#include "ColorRule-gemmi.hh"

namespace coot {
   namespace m2t {

      class ColorScheme {
      private:
         std::string name;
         std::list<std::shared_ptr<ColorRule> > rules;
      public:
         ColorScheme() : name("empty") {}
         ColorScheme(std::shared_ptr<ColorRule> colorRule) { rules.push_back(colorRule); }
         ~ColorScheme() {}

         void addRule(std::shared_ptr<ColorRule> _rule) {
            if (_rule) {
               if (std::find(rules.begin(), rules.end(), _rule) == rules.end()) {
                  if (_rule->getRank() < 0) _rule->setRank(rules.size());
                  rules.push_back(_rule);
                  rules.sort(ColorRule::compareRank);
               }
            } else {
               std::cout << "Error:: ColorScheme::addRule(): null rule" << std::endl;
            }
         }
         void removeRule(std::shared_ptr<ColorRule> _rule) {
            rules.remove(_rule);
            rules.sort(ColorRule::compareRank);
         }
         std::list<std::shared_ptr<ColorRule> > &getColorRules() { return rules; }
         void setName(const std::string &_name) { name = _name; }
         const std::string &getName() const { return name; }

         // last matching rule (in rank order) wins, mirroring the original.
         FCXXCoord colorForAtom(const gemmi::CRA &cra) {
            FCXXCoord result(1., 1., 1., 1.);
            rules.sort(ColorRule::compareRank);
            for (auto &rule : rules) {
               auto cs = rule->getCompoundSelection();
               if (cs && cs->matches(cra))
                  result = rule->colorForAtom(cra);
            }
            return result;
         }

         static std::shared_ptr<ColorScheme> colorByElementScheme();
         static std::shared_ptr<ColorScheme> colorBySecondaryScheme();
         static std::shared_ptr<ColorScheme> colorBFactorScheme();
         static std::shared_ptr<ColorScheme> colorRampChainsScheme();
         static std::shared_ptr<ColorScheme> colorChainsScheme();
         static std::shared_ptr<ColorScheme> colorChainsSchemeWithColourRules(const std::vector<std::pair<std::string, std::string> > &colour_rules);
         static std::shared_ptr<ColorScheme> colorSchemeForColorName(const std::string name);
      };
   }
}

#endif // ColorScheme_gemmi_hh
