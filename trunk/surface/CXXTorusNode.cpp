/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
 *  CXXTorusNode.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  
 *
 */

#include "CXXTorusNode.h"


void CXXTorusNode::init(){
	theAtom = 0;
	theta = omega = 0.;
	crd = CXXCoord(0.,0.,0.);
	ind = 0;
}	

CXXTorusNode::CXXTorusNode()
{
	init();
}

CXXTorusNode::~CXXTorusNode()
{
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
const CXXCoord &CXXTorusNode::coord() const{
	return crd;
}
const double CXXTorusNode::getTheta() const{

	return theta;
}
const double CXXTorusNode::getOmega() const{
	return omega;
}
int CXXTorusNode::setIndex(int i){
	ind = i;
	return 0;
}
const int CXXTorusNode::index(void) const{
	return ind;
}
int CXXTorusNode::setCoord(const CXXCoord &c){
	crd = c;
	return 0;
}
int CXXTorusNode::setAtom(CAtom *anAtom){
	theAtom = anAtom;
	return 0;
}
CAtom *CXXTorusNode::getAtom() const{
	return theAtom;
}

