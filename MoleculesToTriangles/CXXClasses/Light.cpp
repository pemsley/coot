/*
 *  Light.mm
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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
