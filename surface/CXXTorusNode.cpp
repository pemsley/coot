/*
 *  CXXTorusNode.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXTorusNode.h"


void CXX_mot::CXXTorusNode::init(){
	theAtom = 0;
	theta = omega = 0.;
	crd = CXXCoord(0.,0.,0.);
}	

CXX_mot::CXXTorusNode::CXXTorusNode()
{
	init();
}
CXX_mot::CXXTorusNode::CXXTorusNode(double inTheta, double inOmega){
	init();
	theta = inTheta;
	omega = inOmega;
}
int CXX_mot::CXXTorusNode::setTheta(const double inTheta){
	theta = inTheta;
	return 0;
}
int CXX_mot::CXXTorusNode::setOmega(const double inOmega){
	omega = inOmega;
	return 0;
}
int CXX_mot::CXXTorusNode::setCoord(const CXXCoord &aCoord){
    crd = aCoord;
	return 0;
}
const CXX_mot::CXXCoord &CXX_mot::CXXTorusNode::coord() const{
	return crd;
}
const double CXX_mot::CXXTorusNode::getTheta() const{

	return theta;
}
const double CXX_mot::CXXTorusNode::getOmega() const{
	return omega;
}
int CXX_mot::CXXTorusNode::setAtom(mmdb::Atom *anAtom){
	theAtom = anAtom;
	return 0;
}
mmdb::Atom *CXX_mot::CXXTorusNode::getAtom() const{
	return theAtom;
}

