/*
 * MoleculesToTriangles/CXXClasses/RendererGL.h
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
