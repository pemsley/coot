/*
 * MoleculesToTriangles/CXXClasses/DisplayPrimitive.cpp
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

#include "DisplayPrimitive.h"
#include "Renderer.h"
#include <iostream>

DisplayPrimitive::~DisplayPrimitive()
{
    liberateAllHandles();
    //std::cout << "In displayPrimitive destructor" << primitiveType;
}

void DisplayPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
    std::cout << "Base class renderWithRenderer" << std::endl;
}

void DisplayPrimitive::liberateAllHandles()
{
    //std::cout << "In liberateAllHandles\n";
    auto renderersIter = renderers.begin();
    for (; renderersIter != renderers.end(); ++renderersIter){
        //std::cout << "Identified a renderer to disarm\n";
        (*renderersIter)->liberateHandlesForDisplayPrimitive(this);
    }
}
