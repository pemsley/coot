/*
 *  CXXBall.h
 *  CXXSurface
 *
 *  Created by Martin Noble on 10/06/2009.
 *  Copyright 2009 LMB, Oxford University. All rights reserved.
 *
 */

#ifndef CXX_mot_CXXBall_included
#define CXX_mot_CXXBall_included

// 20110217 4.3.3 on Ubuntu 9.04 needs this hack too
#ifdef __GNUC__
#ifdef __GNUC_MINOR__
#if (__GNUC__ == 4)
#if    ( (__GNUC_MINOR__ == 1) ||  (__GNUC_MINOR__ == 2) || (__GNUC_MINOR__ == 3))
#define NEED_OPENMP_PRAGMA_HACK
#endif 
#endif 
#endif 
#endif


// #include <string.h>
#include <cstring>
#ifndef  __MMDB_Manager__
#include <mmdb2/mmdb_manager.h>
#endif

#include "CXXCoord.h"
#include "CXXSphereElement.h"
#include "CXXCircleNode.h"

namespace CXX_mot {

class CXXBall {
protected:
	CXXCoord theCoord;
	double theRadius;
public:
	const CXXCoord_ftype& operator [] (const int i) const{
		return theCoord[i];
	};
	CXXCoord_ftype& operator [] (const int i) {
		return theCoord[i];
	};
	const CXXCoord &getCoord () const{
		return theCoord;
	};
	CXXCoord &getCoord () {
		return theCoord;
	};
	virtual const double &getRadius() const = 0;
	virtual void initSphereElement(CXXSphereElement &, const double &delta) const = 0;	
	
	static int triangulateBalls(vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &ballPntrs,
                                vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &contextBallPntrs,
								double delta, CXXSurface *aSurface, int insideOrOutside);
	static int ballContacts(std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &ballPntrs, 
                            std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > &contextBallPntrs, 
							std::map<const CXXBall*, std::vector<const CXXBall*, CXX_old::CXXAlloc<const CXXBall*> > > &contactMap);
	virtual mmdb::PAtom getAtomI() const = 0;
};

class CXXAtomBall: public CXXBall {
private:
	mmdb::PAtom theAtom;
	double theRadius;
	static double mostRecentDelta;
	static CXXSphereElement unitSphereAtOrigin;
public:
	CXXAtomBall(mmdb::PAtom theAtom_in, const double &radius_in) : theAtom (theAtom_in), theRadius(radius_in){
		theCoord=CXXCoord(theAtom->x, theAtom->y, theAtom->z);
	};
	virtual const double &getRadius() const{
		return theRadius;
	};
	virtual void initSphereElement(CXXSphereElement &theSphere, const double &delta) const{
		if (delta != mostRecentDelta){
#pragma omp critical (unitSphereAtOrigin) 
			{
				unitSphereAtOrigin = CXXSphereElement(CXXCoord(0.,0.,0.), 1., delta);
				mostRecentDelta = delta;
			}
		}
		theSphere = unitSphereAtOrigin;
		theSphere.scaleBy(theRadius);		
		theSphere.translateBy (theCoord);
		theSphere.setAtom(theAtom);
	};	
    virtual mmdb::PAtom getAtomI() const{
        return theAtom;
    };
};

class CXXReentrantProbeBall : public CXXBall {
private:
	mmdb::PAtom theAtomI;
	mmdb::PAtom theAtomJ;
	mmdb::PAtom theAtomK;
	bool includeAtoms[3];
	double theRadius;
	static CXX_old::CXXAlloc<CXXReentrantProbeBall> allocator;
public:
	CXXReentrantProbeBall() : CXXBall() {
	};
	CXXReentrantProbeBall(const CXXCircleNode &parentNode, const int &selHnd, const double &radius_in) : 
	theAtomI(parentNode.getAtomI()),theAtomJ(parentNode.getAtomJ()), theAtomK(parentNode.getAtomK()), 
	theRadius(radius_in){
		theCoord = parentNode.getCoord();
		includeAtoms[0] = theAtomI->isInSelection(selHnd);
		includeAtoms[1] = theAtomJ->isInSelection(selHnd);
		includeAtoms[2] = theAtomK->isInSelection(selHnd);
	};
	virtual const double &getRadius() const{
		return theRadius;
	};
	virtual void initSphereElement(CXXSphereElement &theSphere, const double &delta) const {
		theSphere.initWith(theCoord, theAtomI, theAtomJ, theAtomK, 
                           delta, theRadius, includeAtoms);
		
	};
    virtual mmdb::Atom *getAtomI() const{
        return theAtomI;
    };
    mmdb::Atom *getAtomJ() const {
        return theAtomJ;
    };
    mmdb::Atom *getAtomK() const {
        return theAtomK;
    };
    static bool equalsPntr(const CXXBall* &ball1, const CXXBall* &ball2){
        const CXXReentrantProbeBall &node1(*static_cast<const CXXReentrantProbeBall *>(ball1));
        const CXXReentrantProbeBall &node2(*static_cast<const CXXReentrantProbeBall *>(ball2));
        std::vector<mmdb::Atom *>ijkCentral(3);
        std::vector<mmdb::Atom *>ijkOther(3);
        ijkCentral[0] = node1.getAtomI();
        ijkCentral[1] = node1.getAtomJ();
        ijkCentral[2] = node1.getAtomK();
        sort(ijkCentral.begin(), ijkCentral.end());
        ijkOther[0] = node2.getAtomI();
        ijkOther[1] = node2.getAtomJ();
        ijkOther[2] = node2.getAtomK();
        sort(ijkOther.begin(), ijkOther.end());
        if (ijkCentral[0] != ijkOther[0]) return false;
        if (ijkCentral[1] != ijkOther[1]) return false;
        if (ijkCentral[2] != ijkOther[2]) return false;
        if (!node1.getCoord().isNearly(node2.getCoord(), 0.00001)) return false;
        return true;
    };
	/*
	static void *operator new(size_t nObjects) {
		return allocator.allocate(nObjects, 0);
	};
	static void operator delete(void *pntr, size_t objectSize=0) {
		allocator.deallocate(static_cast<CXXReentrantProbeBall *>(pntr), 0);
	};
	 */
};

}
#endif

