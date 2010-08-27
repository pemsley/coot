/*
 *  CXXTorusNode.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXTorusNode.h"


void CXXTorusNode::init(){
	theAtom = 0;
	theta = omega = 0.;
	crd = CXXCoord(0.,0.,0.);
}	

CXXTorusNode::CXXTorusNode()
{
	init();
}
CXXTorusNode::CXXTorusNode(double inTheta, double inOmega){
	init();
	theta = inTheta;
	omega = inOmega;
}
int CXXTorusNode::setTheta(const double inTheta){
	theta = inTheta;
	return 0;
}
int CXXTorusNode::setOmega(const double inOmega){
	omega = inOmega;
	return 0;
}
int CXXTorusNode::setCoord(const CXXCoord &aCoord){
    crd = aCoord;
	return 0;
}
const CXXCoord &CXXTorusNode::coord() const{
	return crd;
}
const double CXXTorusNode::getTheta() const{

	return theta;
}
const double CXXTorusNode::getOmega() const{
	return omega;
}
int CXXTorusNode::setAtom(CAtom *anAtom){
	theAtom = anAtom;
	return 0;
}
CAtom *CXXTorusNode::getAtom() const{
	return theAtom;
}

