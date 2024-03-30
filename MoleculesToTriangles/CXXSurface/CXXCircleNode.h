/*
 *  CXXCircleNode.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXXCircleNode_included
#define CXXCircleNode_included
#include <string.h>
#ifndef  __MMDB_Manager__
#endif

#include "CXXCoord.h"
#include <map>
#include <vector>
class CXXCircle;
#include "mmdb2/mmdb_manager.h"


class CXXCircleNode {
private:
	const CXXCircle *theParent;
	const CXXCircle *theOtherCircle;
	CXXCoord<CXXCoord_ftype>theCoord;
	CXXCoord<CXXCoord_ftype>unitRadius;
	double theAngle;
	int theFlag;
	int thisIsDeleted;
	mmdb::Atom* atomI;
	mmdb::Atom* atomJ;
	mmdb::Atom* atomK;
public:
    CXXCircleNode();    
	CXXCircleNode ( const CXXCircle *aParent, const CXXCircle *anOtherCircle, const CXXCoord<CXXCoord_ftype>&crd, int aFlag);
	int setReference(const CXXCoord<CXXCoord_ftype>&referenceVector);
	
	const CXXCircle *getParent() const{
		return theParent;
	};
	CXXCircle *getParent() {
        //Dodgy const cast
		return const_cast<CXXCircle *>(theParent);
	};
	const CXXCircle *getOtherCircle() const{
		return theOtherCircle;
	};
	const CXXCoord<CXXCoord_ftype>&getCoord() const{
		return theCoord;
	}
	CXXCoord<CXXCoord_ftype>getUnitRadius() const{
		return unitRadius;
	};
	const double &getAngle() const {
		return theAngle;
	};
	mmdb::Atom* getAtomI() const {return atomI;};
	mmdb::Atom* getAtomJ() const {return atomJ;};
	mmdb::Atom* getAtomK() const {return atomK;};

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
	void setCoord(const CXXCoord<CXXCoord_ftype>&coord);	
	static bool shouldDelete(const CXXCircleNode &aNode);
	
	const CXXCoord_ftype& operator [] (unsigned i) const{
		return theCoord[i];
	};

	CXXCoord_ftype operator [] (int i) const{
		return theCoord[i];
	};

	static int probeContacts(std::vector<CXXCircleNode  > &probes, double probeRadius, 
					  std::map<CXXCircleNode *, std::vector< CXXCircleNode *  > > &contactMap);
	static bool shouldDeletePointer(CXXCircleNode* &aNodePointer);
    static bool equalsPntr(CXXCircleNode* &node1, CXXCircleNode* &node2);
    static bool equals(CXXCircleNode &node1, CXXCircleNode &node2);
	static void filterContacts(std::map<CXXCircleNode *, std::vector< CXXCircleNode *  > > &contactMap);
    static bool angleLessThan(const CXXCircleNode &node1, const CXXCircleNode &node2);
};

#endif


