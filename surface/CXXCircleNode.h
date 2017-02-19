/*
 *  CXXCircleNode.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXX_mot_CXXCircleNode_included
#define CXX_mot_CXXCircleNode_included
// #include <string.h>
#include <cstring>
#ifndef  __MMDB_Manager__
#include <mmdb2/mmdb_manager.h>
#endif
#include "CXXAlloc.h"

#include "CXXCoord.h"
//#include "CXXAlloc.h"
#include <map>
#include <vector>

namespace CXX_mot {
class CXXCircle;

class CXXCircleNode {
private:
	const CXXCircle *theParent;
	const CXXCircle *theOtherCircle;
	CXXCoord theCoord;
	CXXCoord unitRadius;
	double theAngle;
	int theFlag;
	int thisIsDeleted;
	mmdb::Atom *atomI;
	mmdb::Atom *atomJ;
	mmdb::Atom *atomK;
public:
    CXXCircleNode();    
	CXXCircleNode ( const CXXCircle *aParent, const CXXCircle *anOtherCircle, const CXXCoord &crd, int aFlag);
	int setReference(const CXXCoord &referenceVector);
	
	const CXXCircle *getParent() const{
		return theParent;
	};
	const CXXCircle *getOtherCircle() const{
		return theOtherCircle;
	};
	const CXXCoord &getCoord() const{
		return theCoord;
	}
	CXXCoord getUnitRadius() const{
		return unitRadius;
	};
	const double &getAngle() const {
		return theAngle;
	};
	const mmdb::PAtom getAtomI() const {return atomI;};
	const mmdb::PAtom getAtomJ() const {return atomJ;};
	const mmdb::PAtom getAtomK() const {return atomK;};

	const int getFlag() const {
		return theFlag;
	};
	void setFlag(int flag) {
		theFlag = flag;
	};
	int operator < (const CXXCircleNode &otherOne) const{
		return (theAngle<(otherOne.getAngle()));
	};
	
	int isDeleted() const{
		return thisIsDeleted;
	};
	
	void setDeleted(const int yesOrNo){
		thisIsDeleted = yesOrNo; 
	};

	void setAngle(double anAngle){
		theAngle = anAngle;
	};
	void setParent(CXXCircle *parent);
	void setOtherCircle(CXXCircle *parent);
	void setCoord(const CXXCoord &coord);	
	static bool shouldDelete(const CXXCircleNode &aNode);
	
	const CXXCoord_ftype& operator [] (unsigned i) const{
		return theCoord[i];
	};

	CXXCoord_ftype operator [] (int i) const{
		return theCoord[i];
	};

	static int probeContacts(std::vector<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> > &probes, double probeRadius, 
					  std::map<CXXCircleNode *, std::vector< CXXCircleNode *, CXX_old::CXXAlloc< CXXCircleNode *> > > &contactMap);
	static bool shouldDeletePointer(CXXCircleNode* &aNodePointer);
    static bool equalsPntr(CXXCircleNode* &node1, CXXCircleNode* &node2);
    static bool equals(CXXCircleNode &node1, CXXCircleNode &node2);
	static void filterContacts(std::map<CXXCircleNode *, std::vector< CXXCircleNode *, CXX_old::CXXAlloc< CXXCircleNode *> > > &contactMap);
    static bool angleLessThan(const CXXCircleNode &node1, const CXXCircleNode &node2);
};

}

#endif
