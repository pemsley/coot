/*
 * MoleculesToTriangles/CXXSurface/CXXCircleNode-gemmi.hh
 *
 * gemmi-native twin of CXXCircleNode.h. Atom token const gemmi::Atom*;
 * CXXCircle used only via pointer (forward-declared).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXCircleNode_gemmi_included
#define CXXCircleNode_gemmi_included
#include <map>
#include <vector>
#include "CXXCoord.h"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXCircle;  // pointer only

      class CXXCircleNode {
      private:
         const CXXCircle *theParent;
         const CXXCircle *theOtherCircle;
         CXXCoord<CXXCoord_ftype> theCoord;
         CXXCoord<CXXCoord_ftype> unitRadius;
         double theAngle;
         int theFlag;
         int thisIsDeleted;
         const gemmi::Atom* atomI;
         const gemmi::Atom* atomJ;
         const gemmi::Atom* atomK;
      public:
         CXXCircleNode();
         CXXCircleNode(const CXXCircle *aParent, const CXXCircle *anOtherCircle, const CXXCoord<CXXCoord_ftype>&crd, int aFlag);
         int setReference(const CXXCoord<CXXCoord_ftype>&referenceVector);

         const CXXCircle *getParent() const { return theParent; }
         CXXCircle *getParent() { return const_cast<CXXCircle *>(theParent); }
         const CXXCircle *getOtherCircle() const { return theOtherCircle; }
         const CXXCoord<CXXCoord_ftype>& getCoord() const { return theCoord; }
         CXXCoord<CXXCoord_ftype> getUnitRadius() const { return unitRadius; }
         const double &getAngle() const { return theAngle; }
         const gemmi::Atom* getAtomI() const { return atomI; }
         const gemmi::Atom* getAtomJ() const { return atomJ; }
         const gemmi::Atom* getAtomK() const { return atomK; }

         const int getFlag() const { return theFlag; }
         void setFlag(int flag) { theFlag = flag; }
         int operator < (const CXXCircleNode &otherOne) const { return (theAngle < (otherOne.getAngle())); }
         int isDeleted() const { return thisIsDeleted; }
         void setDeleted(const int yesOrNo) { thisIsDeleted = yesOrNo; }
         void setAngle(double anAngle) { theAngle = anAngle; }
         void setParent(CXXCircle *parent);
         void setOtherCircle(CXXCircle *parent);
         void setCoord(const CXXCoord<CXXCoord_ftype>&coord);
         static bool shouldDelete(const CXXCircleNode &aNode);

         const CXXCoord_ftype& operator [] (unsigned i) const { return theCoord[i]; }
         CXXCoord_ftype operator [] (int i) const { return theCoord[i]; }

         static int probeContacts(std::vector<CXXCircleNode> &probes, double probeRadius,
                                  std::map<CXXCircleNode *, std::vector<CXXCircleNode *> > &contactMap);
         static bool shouldDeletePointer(CXXCircleNode* &aNodePointer);
         static bool equalsPntr(CXXCircleNode* &node1, CXXCircleNode* &node2);
         static bool equals(CXXCircleNode &node1, CXXCircleNode &node2);
         static void filterContacts(std::map<CXXCircleNode *, std::vector<CXXCircleNode *> > &contactMap);
         static bool angleLessThan(const CXXCircleNode &node1, const CXXCircleNode &node2);
      };
   }
}

#endif
