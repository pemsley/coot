/*
 *  CXXBall.h
 *  CXXSurface
 *
 *  Created by Martin Noble on 10/06/2009.
 *  Copyright 2009 LMB, Oxford University. All rights reserved.
 *
 */

#ifndef CXXBall_included
#define CXXBall_included

#ifdef __GNUC__
#ifdef __GNUC_MINOR__
#if (__GNUC__ == 4)
#if    ( (__GNUC_MINOR__ == 1) ||  (__GNUC_MINOR__ == 2))
#define const_for_openmp_pragma_arg
#define openmp_pragma_reference_arg
#endif 
#endif 
#endif 
#endif

#ifndef const_for_openmp_pragma_arg
#define const_for_openmp_pragma_arg const
#define openmp_pragma_reference_arg &
#endif

#include <string.h>
#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#endif

#include "CXXCoord.h"
#include "CXXSphereElement.h"
#include "CXXCircleNode.h"

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
	
	static int triangulateBalls(vector<const CXXBall*, CXX::CXXAlloc<const CXXBall*> >  openmp_pragma_reference_arg ballPntrs,
                                vector<const CXXBall*, CXX::CXXAlloc<const CXXBall*> > &contextBallPntrs,
								double delta, CXXSurface *aSurface, int insideOrOutside);
	static int ballContacts(std::vector<const CXXBall*, CXX::CXXAlloc<const CXXBall*> > openmp_pragma_reference_arg ballPntrs, 
                            std::vector<const CXXBall*, CXX::CXXAlloc<const CXXBall*> > &contextBallPntrs, 
							std::map<const CXXBall*, std::vector<const CXXBall*, CXX::CXXAlloc<const CXXBall*> > > openmp_pragma_reference_arg contactMap);
	virtual PCAtom getAtomI() const = 0;
};

class CXXAtomBall: public CXXBall {
private:
	PCAtom theAtom;
	double theRadius;
	static double mostRecentDelta;
	static CXXSphereElement unitSphereAtOrigin;
public:
	CXXAtomBall(PCAtom theAtom_in, const double &radius_in) : theAtom (theAtom_in), theRadius(radius_in){
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
    virtual PCAtom getAtomI() const{
        return theAtom;
    };
};

class CXXReentrantProbeBall : public CXXBall {
private:
	PCAtom theAtomI;
	PCAtom theAtomJ;
	PCAtom theAtomK;
	bool includeAtoms[3];
	double theRadius;
	static CXX::CXXAlloc<CXXReentrantProbeBall> allocator;
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
    virtual CAtom *getAtomI() const{
        return theAtomI;
    };
    CAtom *getAtomJ() const {
        return theAtomJ;
    };
    CAtom *getAtomK() const {
        return theAtomK;
    };
    static bool equalsPntr(const CXXBall* &ball1, const CXXBall* &ball2){
        const CXXReentrantProbeBall &node1(*static_cast<const CXXReentrantProbeBall *>(ball1));
        const CXXReentrantProbeBall &node2(*static_cast<const CXXReentrantProbeBall *>(ball2));
        std::vector<CAtom *>ijkCentral(3);
        std::vector<CAtom *>ijkOther(3);
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

#endif

