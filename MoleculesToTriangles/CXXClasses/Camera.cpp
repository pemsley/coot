/*
 *  Camera.mm
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "Camera.h"
#include "Renderer.h"
#include "SceneSetup.h"

void Camera::renderWithRenderer(std::shared_ptr<Renderer> renderer) {
#ifdef DEBUG_MINE
    std::cout << "In camera render\n";
#endif
    renderer->setupCamera(this);
    sceneSetup->renderWithRendererFromViewpoint(renderer, FCXXCoord (getTranslationX(), getTranslationY(), getTranslationZ()));
};
