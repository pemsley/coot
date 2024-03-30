/*
 *  Renderer.h
 *  Aesop
 *
 *  Created by Martin Noble on 03/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef Renderer_H
#define Renderer_H

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include <map>
#include <vector>

class Light;
class RepresentationInstance;
class Camera;
class SceneSetup;
class DisplayPrimitive;
//class LinesPrimitive;
class BondsPrimitive;
class oglPolyhedron;
class VertexColorNormalPrimitive;

#ifndef RendererHandles_type
#define RendererHandles_type
typedef struct RendererHandles_ {
    unsigned int vertexHandle;
    unsigned int colorHandle;
    unsigned int indexHandle;
    unsigned int arrayObjectHandle;
} RendererHandles;
#endif

class Renderer {
protected:
    virtual void vboRenderVCN(VertexColorNormalPrimitive *prim) = 0;
    virtual void vboRenderVCNFixed(VertexColorNormalPrimitive *prim) = 0;
    virtual void vboRenderVN(VertexColorNormalPrimitive *prim) = 0;
    virtual void vboRenderVC(VertexColorNormalPrimitive *prim) = 0;
    std::map<DisplayPrimitive *, RendererHandles> allocatedHandles;
    Camera *camera;
public:
    virtual void init() = 0;
    virtual void setupCamera(Camera *camera) = 0;
    virtual void setupScene(SceneSetup *sceneSetup) = 0;
    virtual void setupLightAsIndexFromViewpointWithScale(Light *light, int asIndex, FCXXCoord  fromViewpoint, float scale) = 0;
    virtual void setupRepresentationInstance(RepresentationInstance *instance) = 0;
    virtual void restoreModelviewMatrix() = 0;
    //virtual void renderLinesPrimitive(LinesPrimitive *lines) = 0;
    virtual void renderBondsPrimitive(BondsPrimitive *bonds) = 0;
    virtual void renderPolyhedron(oglPolyhedron *polygon) = 0;
    virtual void renderVertexColorNormalPrimitive(VertexColorNormalPrimitive *prim) = 0;
    virtual void renderVertexColorPrimitive(VertexColorNormalPrimitive *prim) = 0;
    virtual FCXXCoord  angstromsForPixels(float x, float y) = 0;
    virtual void liberateHandlesForDisplayPrimitive(DisplayPrimitive *prim) = 0;
    virtual void render(Camera *camera) = 0;
    virtual void clearCameraCanvas(Camera *camera) = 0;
    virtual void drawTestTriangle(const FCXXCoord &) const = 0;
    virtual void resize(int w, int h) = 0;
    //void render();
};

#endif
