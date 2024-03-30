/*
 * MoleculesToTriangles/CXXClasses/CameraPort.h
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

#ifndef CameraPort_h
#define CameraPort_h
#include <memory>

#include <stdio.h>
#include "GLFW/glfw3.h"
#include "Camera.h"
#include "RendererGLSL.hpp"
#include "MolecularRepresentationInstance.h"
#include "mmdb2/mmdb_manager.h"

class ColorScheme;

class CameraPort {
private:
    GLFWwindow* window;
    std::shared_ptr<Camera> camera = 0;
    std::shared_ptr<Renderer> renderer;
    std::shared_ptr<SceneSetup> sceneSetup;
    
    static const std::vector<std::shared_ptr<ColorScheme> >colorSchemes;
    static std::vector<std::shared_ptr<ColorScheme> > createColorSchemes();
    int iCurrentColorScheme;
    
    static const std::vector<std::string>displayStyles;
    static std::vector<std::string> createDisplayStyles();
    int iCurrentDisplayStyle;
    
    double oldX, oldY;
public:
    CameraPort();
    
    void initialiseMoleculesToTriangles();
    
    std::shared_ptr<MolecularRepresentationInstance> addRepresentationInstance(mmdb::Manager *ofMmdb, std::shared_ptr<ColorScheme> colorScheme, std::string selectionString, std::string renderStyle);
    
    void renderSceneSetup();
    
    
    
    void setCamera(std::shared_ptr<Camera> _camera) {camera=_camera;};
    std::shared_ptr<Camera> getCamera() { return camera; };
    void setRenderer(std::shared_ptr<Renderer> _renderer) {renderer=_renderer;};
    GLFWwindow *getWindow(){return window;};
    void runLoop(std::string message);
    static std::shared_ptr<CameraPort> create();
    double getOldX() { return oldX; };
    double getOldY() { return oldY; };
    void setOldX(double newX) { oldX = newX; };
    void setOldY(double newY) { oldY = newY; };
    void zoom(float factor);
    void cycleColorScheme();
    void cycleDisplayStyle();
    
    //Interactivity
    static void cursorMoved(GLFWwindow *win, double newX, double newY);
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
};

#endif /* CameraPort_h */
