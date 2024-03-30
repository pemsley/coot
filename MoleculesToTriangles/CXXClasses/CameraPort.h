//
//  CameraPort.h
//  MoleculesToTriangles
//
//  Created by Martin Noble on 27/01/2017.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

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
