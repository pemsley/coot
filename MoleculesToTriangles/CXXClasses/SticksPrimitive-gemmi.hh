/*
 * MoleculesToTriangles/CXXClasses/SticksPrimitive-gemmi.hh
 *
 * gemmi-native twin of SticksPrimitive.h. Renders half-bonds as lines (the base
 * VertexColorNormalPrimitive, DrawAsLines). The original keyed on mmdb::Atom* and
 * evaluated colour internally; this twin takes coord+colour per atom (keyed by a
 * caller-supplied stable id, e.g. gemmi atom.serial) + bond pairs, colour eager.
 * Type is coot::m2t::SticksPrimitive. Also used for the Calpha line trace (the
 * original BondsPrimitive renders via a private-array path that a twin can't reuse).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef SticksPrimitive_gemmi_hh
#define SticksPrimitive_gemmi_hh

#include <map>
#include <vector>
#include <utility>
#include <memory>
#include "VertexColorNormalPrimitive.h"

class Renderer;

namespace coot {
   namespace m2t {

      class SticksPrimitive : public VertexColorNormalPrimitive {
      private:
         struct Node { FCXXCoord coord; FCXXCoord color; };
         std::map<long, Node> atomNodes;            // key = stable atom id (serial)
         std::vector<std::pair<long, long> > bondPairs;
      public:
         SticksPrimitive() {
            _nLines = 0;
            vertexColorArray = 0;
            indexArray = 0;
            emptyArrays();
            drawModeGL = DrawAsLines;
            enableColorGL = true;
         }
         void addAtom(long key, const FCXXCoord &coord, const FCXXCoord &color) {
            atomNodes[key] = Node{coord, color};
         }
         void addBond(long key1, long key2) {
            bondPairs.push_back(std::make_pair(key1, key2));
         }
         void emptyArrays() {
            _nLines = 0;
            delete [] vertexColorArray; vertexColorArray = 0;
            delete [] indexArray; indexArray = 0;
         }
         virtual ~SticksPrimitive() { emptyArrays(); }
         virtual void generateArrays();
         virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);
      };
   }
}

#endif // SticksPrimitive_gemmi_hh
