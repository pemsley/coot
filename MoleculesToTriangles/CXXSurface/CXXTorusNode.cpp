/*
 * MoleculesToTriangles/CXXSurface/CXXTorusNode.cpp
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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

#include "CXXTorusNode.h"


void CXXTorusNode::init(){
	theAtom = 0;
	theta = omega = 0.;
	crd = CXXCoord<CXXCoord_ftype>(0.,0.,0.);
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
int CXXTorusNode::setCoord(const CXXCoord<CXXCoord_ftype>&aCoord){
    crd = aCoord;
	return 0;
}
const CXXCoord<CXXCoord_ftype>&CXXTorusNode::coord() const{
	return crd;
}
const double CXXTorusNode::getTheta() const{

	return theta;
}
const double CXXTorusNode::getOmega() const{
	return omega;
}
int CXXTorusNode::setAtom(mmdb::Atom* anAtom){
	theAtom = anAtom;
	return 0;
}
mmdb::Atom* CXXTorusNode::getAtom() const{
	return theAtom;
}

