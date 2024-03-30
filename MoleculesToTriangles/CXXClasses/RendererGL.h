/*
 *  RendererGL.h
 *  Aesop
 *
 *  Created by Martin Noble on 03/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef RendererGL_H
#define RendererGL_H

#include "Renderer.h"
#include <map>
#include <vector>
#include <memory>

class RendererGL : public Renderer {
protected:
    virtual void vboRenderVCN(VertexColorNormalPrimitive *prim) ;
    void vboRenderVCNFixed(VertexColorNormalPrimitive *prim) ;
    void vboRenderVN(VertexColorNormalPrimitive *prim) ;
    void vboRenderVC(VertexColorNormalPrimitive *prim) ;
public:
    ~RendererGL();
    virtual void init();
    void setupCamera(Camera *camera);
    void setupScene(SceneSetup *sceneSetup);
    void setupLightAsIndexFromViewpointWithScale(Light *light, int asIndex, FCXXCoord  fromViewpoint, float scale);
    void setupRepresentationInstance(RepresentationInstance *instance);
    void restoreModelviewMatrix();
    //virtual void renderLinesPrimitive(LinesPrimitive *lines) ;
    void renderBondsPrimitive(BondsPrimitive *bonds)  ;
    void renderPolyhedron(oglPolyhedron *polygon) ;
    virtual void renderVertexColorNormalPrimitive(VertexColorNormalPrimitive *prim) ;
    void renderVertexColorPrimitive(VertexColorNormalPrimitive *prim) ;
    FCXXCoord  angstromsForPixels(float x, float y);
    void liberateHandlesForDisplayPrimitive(DisplayPrimitive *prim) ;
    void render(Camera *camera);
    void clearCameraCanvas(Camera *camera); 
    
    void myglEnable(int flag);
    void myglDisable(int flag);
    void myglClearDepth(float value);
    void myglGetFloatv(int propertyToGet, float *dest);
    void myglGetIntegerv(int propertyToGet, int *dest);

#pragma mark -
#pragma mark Debug routine
    void drawTestSquare();
    virtual void drawTestTriangle(const FCXXCoord &) const ;
    
    static std::shared_ptr<Renderer> create();
    
    virtual void resize(int w, int h);
};

#endif
