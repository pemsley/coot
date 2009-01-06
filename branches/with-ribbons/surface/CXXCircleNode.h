/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
/*
 *  CXXCircleNode.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  
 *
 */

#ifndef CXXCircleNode_included
#define CXXCircleNode_included
#include <string.h>
#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#endif

#include "CXXCoord.h"
class CXXCircle;

class CXXCircleNode {
private:
	const CXXCircle *theParent;
	CXXCoord theCoord;
	CXXCoord unitRadius;
	double theAngle;
	int theFlag;
	PCAtom theAtomI;
	PCAtom theAtomJ;
	PCAtom theAtomK;
	int thisIsDeleted;
public:
		CXXCircleNode();
	
	CXXCircleNode ( const CXXCircle *aParent, PCAtom atomK, const CXXCoord &crd, int aFlag);
	
	int setReference(const CXXCoord &referenceVector);
	
	const CXXCircle *getParent() const{
		return theParent;
	};
	const CXXCoord &getCoord() const{
		return theCoord;
	}
	const CXXCoord &getUnitRadius() const{
		return unitRadius;
	};
	const double getAngle() const {
		return theAngle;
	};
	const PCAtom getAtomI() const{
		return theAtomI;
	};
	const PCAtom getAtomJ() const{
		return theAtomJ;
	};
	const PCAtom getAtomK() const{
		return theAtomK;
	};
	void setAtomK(PCAtom anAtom){
		theAtomK = anAtom;
	};
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
	void setParent(CXXCircle *parent){
		theParent = parent;
	};
};
#endif


