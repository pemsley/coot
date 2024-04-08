/*
 * MoleculesToTriangles/CXXClasses/RendererGLSL.hpp
 *
 * Copyright 2017 by Martin Noble, University of Oxford
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

#ifndef RendererGLSL_hpp
#define RendererGLSL_hpp

#ifdef WINDOWS_MINGW
typedef unsigned int uint;
#endif

#include <stdio.h>
#include "RendererGL.h"

class RendererGLSL : public RendererGL
{
private:
    int program;
    virtual void vboRenderVCN(VertexColorNormalPrimitive *prim);
    void myglMaterialfv(int faceEnum, int propertyEnum, float *values);
    void myglMaterialf(int faceEnum, int propertyEnum, float value);
    void myglEnableClientState(int clientStateEnum);
    void myglColorPointer(int size, int type, int stride, const void *pointer);
    void myglVertexPointer(int size, int type, int stride, const void *pointer);
    void myglNormalPointer(int type, int stride, const void *pointer);
    std::map<std::string, int>uniforms;
public:
    virtual void renderVertexColorNormalPrimitive(VertexColorNormalPrimitive *prim);
    virtual void init();
    void setProgram(int _program);
    void Use() const;
    void test(){};
    static std::string PointLightFragmentShaderText;
    static std::string PointLightVertexShaderText;
    static std::shared_ptr<Renderer> create();
    std::pair<bool, std::string> readShaderFile(std::string strFilename) const;
    bool loadShaderFile(std::string strFilename, uint iHandle);
    void loadShaders();
};



#endif /* RendererGLSL_hpp */
