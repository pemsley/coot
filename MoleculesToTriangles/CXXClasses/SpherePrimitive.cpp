/*
 *  SpherePrimitive.cpp
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 18/07/2008.
 *  Copyright 2008 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */

#include "SpherePrimitive.h"
#include "oglPolyhedron.h"

SpherePrimitive::SpherePrimitive()
{
}
float SpherePrimitive::radius()
{
    return _radius;
}
FCXXCoord &SpherePrimitive::color()
{
    return _color;
}
FCXXCoord &SpherePrimitive::centre()
{
    return _centre;
}

void SpherePrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
	oglPolyhedron *polyhedron = oglPolyhedron::dodecaHedron16();
    polyhedron->setRadius(_radius);
    polyhedron->setCentre(_centre);
    polyhedron->setColor(_color);
	polyhedron->renderWithRenderer(renderer);
}

void SpherePrimitive::generateArrays() {
    
}
