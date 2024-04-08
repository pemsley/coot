/*
 * MoleculesToTriangles/CXXClasses/Light.cpp
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

#include "Light.h"
#include "Renderer.h"

void Light::setupInRendererAsIndexFromViewpointWithScale (std::shared_ptr<Renderer> renderer, int asIndex, FCXXCoord  fromViewpoint, float scale) {
#ifdef DEBUG_MINE
    std::cout << "In light setter upper\n";
#endif
    renderer->setupLightAsIndexFromViewpointWithScale(this, asIndex, fromViewpoint, scale);
}

std::shared_ptr<Light> Light::defaultLight() {
    std::shared_ptr<Light> result(new Light());
    result->setAmbient( FCXXCoord(1.0, 1.0, 1.0, 1.0));
    result->setDiffuse( FCXXCoord(1.0, 1.0, 1.0, 0.1));
    result->setSpecular(FCXXCoord(1.0, 1.0, 1.0, 0.1));
    result->setTranslation(FCXXCoord(40.0,40.0,40.0,1.0));
    result->setRotation(FCXXCoord(0.0,1.0,0.0,0.0));
    result->setIntensity(1.);
    result->setExponent(2);
    result->setLightType(Directional);
    result->setDrawLight(false);
    return result;
}
