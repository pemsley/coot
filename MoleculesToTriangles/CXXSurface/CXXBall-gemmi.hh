/*
 * MoleculesToTriangles/CXXSurface/CXXBall-gemmi.hh
 *
 * gemmi-native twin of CXXBall.h. Atom token const gemmi::Atom*; the reentrant
 * probe's selection test uses a std::set<const gemmi::Atom*> (was mmdb selHnd).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXBall_gemmi_included
#define CXXBall_gemmi_included

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "CXXCoord.h"
#include "CXXSphereElement-gemmi.hh"
#include "CXXCircleNode-gemmi.hh"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXSurfaceMaker;

      class CXXBall {
      protected:
         CXXCoord<CXXCoord_ftype> theCoord;
         double theRadius;
      public:
         const CXXCoord_ftype& operator [] (const int i) const { return theCoord[i]; }
         CXXCoord_ftype& operator [] (const int i) { return theCoord[i]; }
         const CXXCoord<CXXCoord_ftype>& getCoord() const { return theCoord; }
         CXXCoord<CXXCoord_ftype>& getCoord() { return theCoord; }
         virtual const double &getRadius() const = 0;
         virtual void initSphereElement(CXXSphereElement &, const double &delta, const CXXSphereElement &unitCellAtOriginForDelta) const = 0;

         static int triangulateBalls(std::vector<const CXXBall*> &ballPntrs,
                                     std::vector<const CXXBall*> &contextBallPntrs,
                                     double delta, CXXSurfaceMaker *aSurface, int insideOrOutside);
         static int ballContacts(std::vector<const CXXBall*> &ballPntrs,
                                 std::vector<const CXXBall*> &contextBallPntrs,
                                 std::map<const CXXBall*, std::vector<const CXXBall*> > &contactMap);
         virtual const gemmi::Atom* getAtomI() const = 0;
      };

      class CXXAtomBall : public CXXBall {
      private:
         const gemmi::Atom* theAtom;
         double theRadius;
         static double mostRecentDelta;
         static CXXSphereElement unitSphereAtOrigin;
      public:
         CXXAtomBall(const gemmi::Atom* theAtom_in, const double &radius_in) : theAtom(theAtom_in), theRadius(radius_in) {
            theCoord = CXXCoord<CXXCoord_ftype>(theAtom->pos.x, theAtom->pos.y, theAtom->pos.z);
         }
         virtual const double &getRadius() const { return theRadius; }
         virtual void initSphereElement(CXXSphereElement &theSphere, const double &delta, const CXXSphereElement &unitCellAtOriginForDelta) const;
         virtual const gemmi::Atom* getAtomI() const { return theAtom; }
      };

      class CXXReentrantProbeBall : public CXXBall {
      private:
         const gemmi::Atom* theAtomI;
         const gemmi::Atom* theAtomJ;
         const gemmi::Atom* theAtomK;
         bool includeAtoms[3];
         double theRadius;
      public:
         CXXReentrantProbeBall() : CXXBall() {}
         CXXReentrantProbeBall(const CXXCircleNode &parentNode, const std::set<const gemmi::Atom*> &selSet, const double &radius_in) :
            theAtomI(parentNode.getAtomI()), theAtomJ(parentNode.getAtomJ()), theAtomK(parentNode.getAtomK()),
            theRadius(radius_in) {
            theCoord = parentNode.getCoord();
            includeAtoms[0] = selSet.count(theAtomI) > 0;
            includeAtoms[1] = selSet.count(theAtomJ) > 0;
            includeAtoms[2] = selSet.count(theAtomK) > 0;
         }
         virtual const double &getRadius() const { return theRadius; }
         virtual void initSphereElement(CXXSphereElement &theSphere, const double &delta, const CXXSphereElement &unitCellAtOriginForDelta) const;
         virtual const gemmi::Atom* getAtomI() const { return theAtomI; }
         const gemmi::Atom* getAtomJ() const { return theAtomJ; }
         const gemmi::Atom* getAtomK() const { return theAtomK; }
         static bool equalsPntr(const CXXBall* &ball1, const CXXBall* &ball2) {
            const CXXReentrantProbeBall &node1(*static_cast<const CXXReentrantProbeBall *>(ball1));
            const CXXReentrantProbeBall &node2(*static_cast<const CXXReentrantProbeBall *>(ball2));
            std::vector<const gemmi::Atom*> ijkCentral(3);
            std::vector<const gemmi::Atom*> ijkOther(3);
            ijkCentral[0] = node1.getAtomI();
            ijkCentral[1] = node1.getAtomJ();
            ijkCentral[2] = node1.getAtomK();
            std::sort(ijkCentral.begin(), ijkCentral.end());
            ijkOther[0] = node2.getAtomI();
            ijkOther[1] = node2.getAtomJ();
            ijkOther[2] = node2.getAtomK();
            std::sort(ijkOther.begin(), ijkOther.end());
            if (ijkCentral[0] != ijkOther[0]) return false;
            if (ijkCentral[1] != ijkOther[1]) return false;
            if (ijkCentral[2] != ijkOther[2]) return false;
            if (!node1.getCoord().isNearly(node2.getCoord(), 0.00001)) return false;
            return true;
         }
      };
   }
}

#endif
